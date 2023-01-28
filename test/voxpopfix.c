#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>

// POSIX (mmap, stat, open,...)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include "world/common.h"
#include "world/harvest.h"
#include "world/cheaptrick.h"
#include "world/d4c.h"
#include "world/synthesis.h"
#include "world/synthesisrealtime.h"

#include <SDL.h>

#include "dr_wav.h"

#define ANALYSIS_MAGIC            (0xa741147a)
#define ANALYSIS_PERIOD_SECONDS   (5e-3)
#define ANALYSIS_F0_FLOOR         (71.0)

#define OTHER_SPECTROGRAM_FFT_LENGTH (1<<12) // problem most noticable here (lower=too low resolution; higher=too blurry)
#define OTHER_SPECTROGRAM_FFT_STORAGE_LENGTH (OTHER_SPECTROGRAM_FFT_LENGTH/2+1)

#define BROKEN_WIDTH_SECONDS_F0               (30e-3)
#define BROKEN_WIDTH_SECONDS_SPECTROGRAM      (16e-3)
#define BROKEN_WIDTH_SECONDS_APERIODICITY     (10e-3)

#define DIVRES (1)

static inline double bessel_I0(double x)
{
	double d = 0.0;
	double ds = 1.0;
	double s = 1.0;
	do {
		d += 2.0;
		ds *= (x*x) / (d*d);
		s += ds;
	} while (ds > s*1e-7);
	return s;
}

static inline double kaiser_bessel_attenuation_to_alpha(double attenuation)
{
	if (attenuation > 50.0f) {
		return 0.1102f * (attenuation - 8.7f);
	} else if (attenuation > 20.0f) {
		return 0.5842f * pow(attenuation - 21.0f, 0.4f) + 0.07886f * (attenuation - 21.0f);
	} else {
		return 0;
	}
}

static inline double kaiser_bessel(double alpha, double x)
{
	return bessel_I0(alpha * sqrt(1.0f - x*x)) / bessel_I0(alpha);
}

static int get_fft_length_from_sample_rate(int sample_rate)
{
	CheapTrickOption option = { 0 };
	InitializeCheapTrickOption(sample_rate, &option);
	option.f0_floor = ANALYSIS_F0_FLOOR;
	return GetFFTSizeForCheapTrick(sample_rate, &option);
}

static int calc_fft_storage_length_from_fft_length(int fft_length)
{
	return fft_length/2+1;
}

struct analysis_header {
	uint32_t magic;
	uint32_t n_channels;
	uint32_t sample_rate;
	uint32_t n_analysis_frames;
};

static int calc_number_of_analysis_frames(int sample_rate, int n_pcm_frames)
{
	return GetSamplesForHarvest(sample_rate, n_pcm_frames, ANALYSIS_PERIOD_SECONDS * 1e3); // ms
}

static size_t calc_analysis_file_size(int n_channels, int n_analysis_frames, int sample_rate)
{
	int fft_storage_length = calc_fft_storage_length_from_fft_length(get_fft_length_from_sample_rate(sample_rate));
	size_t header_size = sizeof(struct analysis_header);
	size_t body_size = n_analysis_frames*OTHER_SPECTROGRAM_FFT_STORAGE_LENGTH*sizeof(double) + n_channels * n_analysis_frames * ((1 + 2*fft_storage_length)*sizeof(double));
	return header_size + body_size;
}

struct analysis {
	int fd;
	void* data;

	int n_channels;
	int sample_rate;
	int n_analysis_frames;

	double f0_min;
	double f0_max;

	double**  other_spectrogram; // [frame_index][fftbin]
	double**  f0;                // [channel][frame_index]
	double*** spectrogram;       // [channel][frame_index][fftbin]
	double*** aperiodicity;      // [channel][frame_index][fftbin]
};

static inline int analysis_get_fft_storage_length(struct analysis* a)
{
	return calc_fft_storage_length_from_fft_length(get_fft_length_from_sample_rate(a->sample_rate));
}

static inline double* analysis_get_other_spectrogram(struct analysis* a, int frame_index)
{
	if (0 <= frame_index && frame_index < a->n_analysis_frames) {
		return a->other_spectrogram[frame_index];
	} else {
		return NULL;
	}
}

static inline double* analysis_get_f0_ptr(struct analysis* a, int channel, int frame_offset)
{
	assert(0 <= channel && channel < a->n_channels);
	if (0 <= frame_offset && frame_offset < a->n_analysis_frames) {
		return &a->f0[channel][frame_offset];
	} else {
		return NULL;
	}
}

static inline double analysis_get_f0(struct analysis* a, int channel, int frame_index)
{
	double* p = analysis_get_f0_ptr(a, channel, frame_index);
	return p ? *p : 0.0;
}

static inline double** analysis_get_spectrogram_ptr(struct analysis* a, int channel, int frame_offset)
{
	if (0 <= frame_offset && frame_offset < a->n_analysis_frames) {
		return &a->spectrogram[channel][frame_offset];
	} else {
		return NULL;
	}
}

static inline double* analysis_get_spectrogram(struct analysis* a, int channel, int frame_index)
{
	double** p = analysis_get_spectrogram_ptr(a, channel, frame_index);
	return p ? *p : NULL;
}

static inline double** analysis_get_aperiodicity_ptr(struct analysis* a, int channel, int frame_offset)
{
	if (0 <= frame_offset && frame_offset < a->n_analysis_frames) {
		return &a->aperiodicity[channel][frame_offset];
	} else {
		return NULL;
	}
}


static inline double* analysis_get_aperiodicity(struct analysis* a, int channel, int frame_index)
{
	double** p = analysis_get_aperiodicity_ptr(a, channel, frame_index);
	return p ? *p : NULL;
}

static int analysis_open(struct analysis* a, const char* path)
{
	memset(a, 0, sizeof *a);

	const int fd = open(path, O_RDWR);
	if (fd == -1 && errno == ENOENT) return 0;
	if (fd == -1) {
		fprintf(stderr, "could not open ''%s'': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	a->fd = fd;

	struct stat st;
	if (fstat(fd, &st) == -1) {
		fprintf(stderr, "could open but not stat ''%s'': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	const size_t sz = st.st_size;

	a->data = mmap(NULL, sz, PROT_READ, MAP_PRIVATE, fd, 0);
	if (a->data == MAP_FAILED) {
		fprintf(stderr, "mmap on ''%s'' failed: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	void* p = a->data;
	struct analysis_header* h = p;
	p += sizeof *h;

	if ((p - a->data) > sz) {
		fprintf(stderr, "file ''%s'' too small\n", path);
		exit(EXIT_FAILURE);
	}

	if (h->magic != ANALYSIS_MAGIC) {
		fprintf(stderr, "invalid magic value in ''%s''\n", path);
		exit(EXIT_FAILURE);
	}

	const int n_channels = h->n_channels;
	const int sample_rate = h->sample_rate;
	const int n_analysis_frames = h->n_analysis_frames;

	printf("\n%s:\n", path);
	printf("   number of channels:        %d\n", n_channels);
	printf("   sample rate:               %d\n", sample_rate);
	printf("   number of analysis frames: %d\n", n_analysis_frames);

	a->n_channels = n_channels;
	a->sample_rate = sample_rate;
	a->n_analysis_frames = n_analysis_frames;

	const size_t expected_sz = calc_analysis_file_size(n_channels, n_analysis_frames, sample_rate);
	if (sz != expected_sz) {
		fprintf(stderr, "''%s'' has wrong file size; expected %zdB, but the file is %zdB\n", path, expected_sz, sz);
		exit(EXIT_FAILURE);
	}

	double* fp = p;

	a->other_spectrogram = calloc(n_analysis_frames, sizeof(*a->other_spectrogram));
	for (int i = 0; i < n_analysis_frames; i++) {
		a->other_spectrogram[i] = fp;
		fp += OTHER_SPECTROGRAM_FFT_STORAGE_LENGTH;
	}

	const int fft_storage_length = calc_fft_storage_length_from_fft_length(get_fft_length_from_sample_rate(sample_rate));

	a->f0 = calloc(n_channels, sizeof(*a->f0));
	a->spectrogram = calloc(n_channels, sizeof(*a->spectrogram));
	a->aperiodicity = calloc(n_channels, sizeof(*a->aperiodicity));

	for (int ch = 0; ch < n_channels; ch++) {
		a->f0[ch] = fp;
		fp += n_analysis_frames;

		a->spectrogram[ch] = calloc(n_analysis_frames, sizeof(**a->spectrogram));
		a->aperiodicity[ch] = calloc(n_analysis_frames, sizeof(**a->aperiodicity));

		for (int i = 0; i < n_analysis_frames; i++) {
			a->spectrogram[ch][i] = fp;
			fp += fft_storage_length;
		}

		for (int i = 0; i < n_analysis_frames; i++) {
			a->aperiodicity[ch][i] = fp;
			fp += fft_storage_length;
		}
	}

	assert(((void*)fp - a->data) == sz);

	{ // calculate f0 range
		int first = 1;
		double f0_min = 0.0;
		double f0_max = 0.0;
		for (int ch = 0; ch < n_channels; ch++) {
			for (int frame_index = 0; frame_index < n_analysis_frames; frame_index++) {
				double f0 = analysis_get_f0(a, ch, frame_index);
				if (f0 > 0.0) {
					if (first) {
						f0_min = f0_max = f0;
						first = 0;
					} else {
						if (f0 < f0_min) f0_min = f0;
						if (f0 > f0_max) f0_max = f0;
					}
				}
			}
		}
		a->f0_min = f0_min;
		a->f0_max = f0_max;
	}

	return 1;
}

struct edit_frame {
	uint32_t is_broken   :1;
	uint32_t is_dry      :1;
};

struct edit {
	int fd;
	size_t sz;
	void* data;
	struct edit_frame* frames;
};

static void edit_open(struct edit* edit, const char* path, struct analysis* analysis)
{
	memset(edit, 0, sizeof *edit);

	int fd = open(path, O_RDWR | O_CREAT, 0644);
	if (fd == -1) {
		fprintf(stderr, "%s: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	edit->fd = fd;

	struct stat st;
	if (fstat(fd, &st) == -1) {
		fprintf(stderr, "could open but not stat ''%s'': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	const int n = analysis->n_analysis_frames;
	const size_t required_sz = n * sizeof(*edit->frames);

	size_t sz = st.st_size;
	if (sz == 0) {
		sz = required_sz;
		if (posix_fallocate(fd, 0, sz) != 0) {
			fprintf(stderr, "%s: posix_fallocate() failed\n", path);
			exit(EXIT_FAILURE);
		}
	} else if (sz != required_sz) {
		fprintf(stderr, "%s: expected %zdB, but file is %zdB\n", path, required_sz, sz);
		exit(EXIT_FAILURE);
	}

	edit->sz = sz;
	edit->data = mmap(NULL, sz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (edit->data == MAP_FAILED) {
		fprintf(stderr, "mmap on ''%s'' failed: %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}

	void* p = edit->data;
	edit->frames = p;
}

static void edit_close(struct edit* edit)
{
	munmap(edit->data, edit->sz);
	close(edit->fd);
}

#define UI_ROWS \
	ROW(BROKEN,             2.0) \
	ROW(DRY,                1.0) \
	ROW(OTHER_SPECTROGRAM, 20.0) \
	ROW(F0,                 6.0) \
	ROW(SPECTROGRAM,        4.0) \
	ROW(APERIODICITY,       4.0)

enum ui_row {
	#define ROW(NAME,WEIGHT) UI_ROW_ ## NAME,
	UI_ROWS
	#undef ROW
	N_UI_ROWS
};

struct ui_nav {
	int64_t x0;
	// frame_index = x0 * slope + constant
	int64_t x2f_slope_48_16_fixed;
	int64_t x2f_constant_48_16_fixed;
};

static inline int ui_nav_x2f(struct ui_nav* nav, int x)
{
	return ((nav->x0 + (uint64_t)x)* nav->x2f_slope_48_16_fixed + nav->x2f_constant_48_16_fixed) >> 16;
}

static inline void ui_nav_apply_pan(struct ui_nav* nav, int dx)
{
	nav->x0 -= dx;
}

static inline double fixed_48_16_to_double(int64_t v)
{
	return (double)v / (double)(1<<16);
}

static inline int64_t double_to_fixed_48_16(double v)
{
	return (int64_t)(v * (double)(1<<16));
}

static inline void ui_nav_apply_d_zoom(struct ui_nav* nav, int mx, double d_zoom) // fun fact: my brain wants to type "z_doom" instead
{
	double slope = fixed_48_16_to_double(nav->x2f_slope_48_16_fixed);
	double constant = fixed_48_16_to_double(nav->x2f_constant_48_16_fixed);

	const double old_slope = slope;
	const double old_constant = constant;
	const double zoom_speed = 1.05;

	slope *= pow(zoom_speed, -d_zoom);

	// must choose new 'constant' so that the frame index under the mouse
	// (mx) remains the same.
	//   frame_index = x * slope + constant
	//   x under mouse is x0+mx
	//   eq: x*old_slope + old_constant = x*slope + constant
	//       constant = (x*old_slope + old_constant) - (x*slope)
	double x = (double)(nav->x0 + mx);
	constant = (x*old_slope + old_constant) - (x*slope);

	nav->x2f_slope_48_16_fixed = double_to_fixed_48_16(slope);
	nav->x2f_constant_48_16_fixed = double_to_fixed_48_16(constant);

	// prevent zero/negative slope
	if (nav->x2f_slope_48_16_fixed <= 1) nav->x2f_slope_48_16_fixed = 1;
}

struct ui {
	int show[N_UI_ROWS];
	int mouseover_frame_index;
	int play_position;
	int painting;
	int paint_frame_min;
	int paint_frame_max;
	int dry_val;
	struct ui_nav nav;
	int last_mx;
	int last_btns;
};

static void ui_init(struct ui* ui)
{
	memset(ui, 0, sizeof *ui);
	for (int i = 0; i < N_UI_ROWS; i++) ui->show[i] = 1;
	ui->nav.x2f_slope_48_16_fixed = 1<<16;
	ui->nav.x2f_constant_48_16_fixed = 0;
}

struct ui_ctx {
	struct analysis* analysis;
	struct edit* edit;
	int y0;
	int w;
	int h;
	uint32_t* pixels;
	struct ui_nav* nav;
	int fft_storage_length;
	int mouseover_frame_index;
	int play_position;
};

static inline int ux_x2f(struct ui_ctx* ctx, int x)
{
	return ui_nav_x2f(ctx->nav, x);
}

static inline int binmap(int y, int h, int n, double scale)
{
	int index = round((1.0 - (double)y / (double)h) * (double)n * scale);
	if (index < 0) index = 0;
	if (index >= n) index = n-1;
	return index;
}


static void ui_ctx_prep(struct ui_ctx* ctx, uint32_t* pixels, int y0, int y1)
{
	ctx->y0 = y0;
	ctx->h = y1-y0;
	ctx->pixels = pixels + y0*ctx->w;
}

static void ux_fill(struct ui_ctx* ctx, uint32_t pixel)
{
	uint32_t* pp = ctx->pixels;
	const int n = ctx->w * ctx->h;
	for (int i = 0; i < n; i++) *(pp++) = pixel;
}

static inline void ux_plot(struct ui_ctx* ctx, int x, int y, uint32_t pixel)
{
	ctx->pixels[x + y*ctx->w] = pixel;
}

static inline uint32_t rgba_add(uint32_t a, uint32_t b)
{
	uint32_t r = 0;
	for (int shift = 0; shift < 32; shift+=8) {
		int va = (a >> shift) & 0xff;
		int vb = (b >> shift) & 0xff;
		int m = va+vb;
		if (m > 0xff) m = 0xff;
		r += m << shift;
	}
	return r;
}

static inline void ux_plot_additive(struct ui_ctx* ctx, int x, int y, uint32_t pixel)
{
	uint32_t* p = &ctx->pixels[x + y*ctx->w];
	*p = rgba_add(*p, pixel);
}

#define MOUSEOVER (1<<0)
#define PLAY_POSITION (1<<1)

static inline int ux_frame_index_flags(struct ui_ctx* ctx, int frame_index)
{
	int flags = 0;
	if (frame_index == ctx->mouseover_frame_index) {
		flags |= MOUSEOVER;
	}
	if (frame_index == ctx->play_position) {
		flags |= PLAY_POSITION;
	}
	return flags;
}

static int get_broken_info(struct analysis* a, struct edit* e, int frame_index, double width_seconds, int* out_first_broken, int* out_last_broken)
{
	const int df = round((width_seconds / ANALYSIS_PERIOD_SECONDS) / 2.0);
	const int f0 = frame_index - df;
	const int f1 = frame_index + df;
	for (int i = f0; i <= f1; i++) {
		if (0 <= i && i < a->n_analysis_frames && e->frames[i].is_broken) {
			if (out_first_broken) {
				int i = frame_index;
				while (get_broken_info(a, e, --i, width_seconds, NULL, NULL)) {}
				i++;
				*out_first_broken = i;
			}
			if (out_last_broken) {
				int i = frame_index;
				while (get_broken_info(a, e, ++i, width_seconds, NULL, NULL)) {}
				i--;
				*out_last_broken = i;
			}
			return 1;
		}
	}
	return 0;
}

static inline int ux_is_broken_in_width(struct ui_ctx* ctx, int frame_index, double width_seconds)
{
	return get_broken_info(ctx->analysis, ctx->edit, frame_index, width_seconds, NULL, NULL);
}

static void ui_render_BROKEN(struct ui_ctx* ctx)
{
	ux_fill(ctx, 0x220000ff);

	struct analysis* a = ctx->analysis;
	struct edit* e = ctx->edit;
	const int n_analysis_frames = a->n_analysis_frames;

	for (int x = 0; x < ctx->w; x++) {
		int frame_index = ux_x2f(ctx, x);
		if (0 <= frame_index && frame_index < n_analysis_frames) {
			const int is_broken = e->frames[frame_index].is_broken;
			const int flags = ux_frame_index_flags(ctx, frame_index);
			uint32_t pixel = 0;
			if (is_broken) pixel |= 0xff0000ff;
			if (flags & MOUSEOVER) pixel |= 0x0000ffff;
			if (flags & PLAY_POSITION) pixel |= 0x222222ff;

			for (int y = 0; y < ctx->h; y++) {
				ux_plot(ctx, x, y, pixel);
			}
		}
	}
}

static void ui_render_DRY(struct ui_ctx* ctx)
{
	ux_fill(ctx, 0x222222ff);

	struct analysis* a = ctx->analysis;
	struct edit* e = ctx->edit;
	const int n_analysis_frames = a->n_analysis_frames;

	for (int x = 0; x < ctx->w; x++) {
		int frame_index = ux_x2f(ctx, x);
		if (0 <= frame_index && frame_index < n_analysis_frames) {
			const int is_dry = e->frames[frame_index].is_dry;
			const int flags = ux_frame_index_flags(ctx, frame_index);
			uint32_t pixel = 0;
			if (is_dry) pixel |= 0xccccccff;
			if (flags & MOUSEOVER) pixel |= 0xffffffff;
			if (flags & PLAY_POSITION) pixel |= 0x222222ff;

			for (int y = 0; y < ctx->h; y++) {
				ux_plot(ctx, x, y, pixel);
			}
		}
	}
}

static inline uint32_t f64rgbenc(double r, double g, double b)
{
	if (r < 0.0) r = 0.0;
	if (r > 1.0) r = 1.0;
	if (g < 0.0) g = 0.0;
	if (g > 1.0) g = 1.0;
	if (b < 0.0) b = 0.0;
	if (b > 1.0) b = 1.0;
	r *= 255.0;
	g *= 255.0;
	b *= 255.0;
	return
		  ((uint32_t)round(r) << 24)
		+ ((uint32_t)round(g) << 16)
		+ ((uint32_t)round(b) << 8)
		+ 0xff; // alpha
}

static void ui_render_OTHER_SPECTROGRAM(struct ui_ctx* ctx)
{
	struct analysis* a = ctx->analysis;
	const int n_analysis_frames = a->n_analysis_frames;
	int last_frame_index = -1;
	const int w = ctx->w;
	const int h = ctx->h;
	uint32_t cached_column[1<<14];
	for (int x = 0; x < w; x++) {
		int frame_index = ux_x2f(ctx, x);
		if (!(0 <= frame_index && frame_index < n_analysis_frames)) continue;

		if (frame_index > last_frame_index) {
			double* xs = analysis_get_other_spectrogram(a, frame_index);
			for (int y = 0; y < h; y++) {
				int bi = binmap(y, ctx->h, OTHER_SPECTROGRAM_FFT_STORAGE_LENGTH, 0.15);
				double v = log(xs[bi] + 1.0) * 0.5;
				double r = v;
				double g = v*v;
				double b = v*v*v;
				cached_column[y] = f64rgbenc(r,g,b);
			}
			last_frame_index = frame_index;
		}

		uint32_t mask = 0;
		const int flags = ux_frame_index_flags(ctx, frame_index);
		if (flags & MOUSEOVER) mask |= 0x0000ff00;
		if (flags & PLAY_POSITION) mask |= 0x222222ff;
		const int is_broken = ux_is_broken_in_width(ctx, frame_index, 0.0);
		if (is_broken) mask |= 0x7f000000;

		for (int y = 0; y < h; y++) ux_plot(ctx, x, y, cached_column[y] | mask);
	}
}

static void ui_render_F0(struct ui_ctx* ctx)
{
	ux_fill(ctx, 0x111100ff);
	struct analysis* a = ctx->analysis;
	const int n_analysis_frames = a->n_analysis_frames;
	for (int x = 0; x < ctx->w; x++) {
		int frame_index = ux_x2f(ctx, x);
		if (0 <= frame_index && frame_index < n_analysis_frames) {
			const int flags = ux_frame_index_flags(ctx, frame_index);
			const int is_broken = ux_is_broken_in_width(ctx, frame_index, BROKEN_WIDTH_SECONDS_F0);

			for (int ch = 0; ch < a->n_channels; ch++) {
				double f0 = analysis_get_f0(a, ch, frame_index);
				if (f0 > 0) {
					double t = (f0 - a->f0_max) / (a->f0_min - a->f0_max);
					double y = t * ctx->h;
					ux_plot_additive(ctx, x, y, 0xaa9944ff);
				} else {
					for (int y = 0; y < ctx->h; y++) {
						ux_plot_additive(ctx, x, y, 0x110022ff);
					}
				}
			}

			uint32_t add = 0;
			if (flags & MOUSEOVER) add |= 0x00007f00;
			if (flags & PLAY_POSITION) add |= 0x222222ff;
			if (is_broken) add |= 0x7f000000;
			if (add) {
				for (int y = 0; y < ctx->h; y++) {
					ux_plot_additive(ctx, x, y, add);
				}
			}
		}
	}
}

static void ui_render_SPECTROGRAM(struct ui_ctx* ctx)
{
	struct analysis* a = ctx->analysis;
	const int n_analysis_frames = a->n_analysis_frames;
	int last_frame_index = -1;
	const int w = ctx->w;
	const int h = ctx->h;
	const int n_channels = a->n_channels;
	uint32_t cached_column[1<<14];
	for (int x = 0; x < w; x++) {
		int frame_index = ux_x2f(ctx, x);
		if (!(0 <= frame_index && frame_index < n_analysis_frames)) continue;

		if (frame_index > last_frame_index) {
			for (int y = 0; y < h; y++) {
				int bi = binmap(y, ctx->h, ctx->fft_storage_length, 0.5);
				double m = 0.0;
				for (int channel = 0; channel < n_channels; channel++) {
					m += analysis_get_spectrogram(a, channel, frame_index)[bi];
				}
				double v = m*3000.0;
				double r = v*v;
				double g = v;
				double b = v*v*v;
				cached_column[y] = f64rgbenc(r,g,b);
			}
			last_frame_index = frame_index;
		}

		uint32_t set = 0;
		const int flags = ux_frame_index_flags(ctx, frame_index);
		if (flags & MOUSEOVER) set |= 0x0000ff00;
		if (flags & PLAY_POSITION) set |= 0x222222ff;

		const int is_broken = ux_is_broken_in_width(ctx, frame_index, BROKEN_WIDTH_SECONDS_SPECTROGRAM);
		if (is_broken) set |= 0x7f000000;

		for (int y = 0; y < h; y++) ux_plot(ctx, x, y, cached_column[y] | set);
	}
}

static void ui_render_APERIODICITY(struct ui_ctx* ctx)
{
	struct analysis* a = ctx->analysis;
	const int n_analysis_frames = a->n_analysis_frames;
	int last_frame_index = -1;
	const int w = ctx->w;
	const int h = ctx->h;
	const int n_channels = a->n_channels;
	uint32_t cached_column[1<<14];
	for (int x = 0; x < w; x++) {
		int frame_index = ux_x2f(ctx, x);
		if (!(0 <= frame_index && frame_index < n_analysis_frames)) continue;

		if (frame_index > last_frame_index) {
			for (int y = 0; y < h; y++) {
				int bi = binmap(y, ctx->h, ctx->fft_storage_length, 1.0);
				double m = 0.0;
				for (int channel = 0; channel < n_channels; channel++) {
					m += analysis_get_aperiodicity(a, channel, frame_index)[bi];
				}
				double v = m;
				double r = v*v;
				double g = v*v*v;
				double b = v;
				cached_column[y] = f64rgbenc(r,g,b);
			}
			last_frame_index = frame_index;
		}

		uint32_t mask = 0xffffffff;
		const int flags = ux_frame_index_flags(ctx, frame_index);
		if (flags & MOUSEOVER) mask &= 0x00ffffff;
		if (flags & PLAY_POSITION) mask &= 0x222222ff;

		const int is_broken = ux_is_broken_in_width(ctx, frame_index, BROKEN_WIDTH_SECONDS_SPECTROGRAM);
		if (is_broken) mask &= 0xffff00ff;

		for (int y = 0; y < h; y++) ux_plot(ctx, x, y, cached_column[y] & mask);
	}
}

static void ui_handle(struct ui* ui, int width, int height, uint32_t* pixels, struct analysis* analysis, struct edit* edit)
{
	double weight_sum = 0.0f;
	{
		int i = 0;
		#define ROW(NAME,WEIGHT) \
			if (ui->show[i]) { weight_sum += WEIGHT; } \
			i++;
		UI_ROWS
		#undef ROW
	}

	if (weight_sum == 0.0) return;

	{
		int i = 0;
		double y0 = 0;
		struct ui_ctx ctx = {
			.analysis = analysis,
			.edit = edit,
			.w = width,
			.nav = &ui->nav,
			.fft_storage_length = analysis_get_fft_storage_length(analysis),
			.mouseover_frame_index = ui->mouseover_frame_index,
			.play_position = ui->play_position,
		};
		#define ROW(NAME,WEIGHT) \
			if (ui->show[i]) { \
				const double dy = (double)height * (WEIGHT / weight_sum); \
				const double y1 = y0 + dy; \
				ui_ctx_prep(&ctx, pixels, floor(y0), floor(y1)); \
				ui_render_##NAME(&ctx); \
				y0 = y1; \
			} \
			i++;

		UI_ROWS
		#undef UI_ROWS
	}
}

static void ui_input(struct ui* ui, int mx, int my, int btns, int mod, float d_zoom, int play_position, struct edit* edit)
{
	const int is_panning = btns & SDL_BUTTON_RMASK;
	if (is_panning) ui_nav_apply_pan(&ui->nav, mx - ui->last_mx);
	if (d_zoom != 0.0) ui_nav_apply_d_zoom(&ui->nav, mx, d_zoom);
	const int mf = ui_nav_x2f(&ui->nav, mx);
	ui->mouseover_frame_index = mf;

	ui->play_position = play_position;

	const int no_mod = (mod == 0);
	const int is_shift = mod & KMOD_SHIFT;
	//const int is_alt = mod & KMOD_ALT;
	const int is_ctrl = mod & KMOD_CTRL;

	const int BROKEN = 1;
	const int ERASE_BROKEN = 2;
	const int DRY_TOGGLE = 3;

	int bdiff = btns ^ ui->last_btns;
	if ((btns & SDL_BUTTON_LMASK) && (bdiff & SDL_BUTTON_LMASK)) {
		if (no_mod) {
			ui->painting = BROKEN;
		} else if (is_shift) {
			ui->painting = ERASE_BROKEN;
		} else if (is_ctrl) {
			ui->painting = DRY_TOGGLE;
			ui->dry_val = !edit->frames[mf].is_dry;
		}
		if (ui->painting) {
			ui->paint_frame_min = ui->paint_frame_max = mf;
		}
	}

	if (ui->painting) {
		if (btns == 0) {
			ui->painting = 0;
		} else {
			if (mf < ui->paint_frame_min) ui->paint_frame_min = mf;
			if (mf > ui->paint_frame_max) ui->paint_frame_max = mf;
			for (int i = ui->paint_frame_min; i <= ui->paint_frame_max; i++) {
				if (ui->painting == BROKEN) {
					edit->frames[i].is_broken = 1;
				} else if (ui->painting == ERASE_BROKEN) {
					edit->frames[i].is_broken = 0;
				} else if (ui->painting == DRY_TOGGLE) {
					edit->frames[i].is_dry = ui->dry_val;
				}
			}
		}
	}

	ui->last_btns = btns;
	ui->last_mx = mx;
}

#define RINGBUF_SZ (1<<6)
#define RINGBUF_CURSOR_MASK (RINGBUF_SZ-1)

int f0_ringbuf_cursor;
double f0_ringbuf_storage[RINGBUF_SZ];


static inline int ringbuf_cursor_advance(int* cursor)
{
	int i = *cursor;
	*cursor = ((*cursor) + 1) & RINGBUF_CURSOR_MASK;
	return i;
}

static double lerp(double a, double b, double u)
{
	return a+(b-a)*u;
}

static double* f0_lerp(double a, double b, double u)
{
	double v = lerp(a,b,u);
	int i = ringbuf_cursor_advance(&f0_ringbuf_cursor);
	f0_ringbuf_storage[i] = v;
	return &f0_ringbuf_storage[i];
}

int fft_ringbuf_cursor;
double* fft_ringbuf_pointers[RINGBUF_SZ];

static double** fft_lerp(int fft_storage_length, double* fft0, double* fft1, double u)
{
	int i = ringbuf_cursor_advance(&fft_ringbuf_cursor);
	double** pp = &fft_ringbuf_pointers[i];
	double* p = *pp;
	for (int i = 0; i < fft_storage_length; i++) {
		p[i] = lerp(fft0[i], fft1[i], u);
	}
	return pp;
}

static double* get_f0_ptr_edited(struct analysis* a, struct edit* e, int channel, int frame_index)
{
	int first_broken, last_broken;
	if (get_broken_info(a, e, frame_index, BROKEN_WIDTH_SECONDS_F0, &first_broken, &last_broken)) {
		const int i0 = first_broken-1;
		const int i1 = last_broken+1;
		double f0_0 = analysis_get_f0(a, channel, i0);
		double f0_1 = analysis_get_f0(a, channel, i1);
		return f0_lerp(f0_0, f0_1, (double)(frame_index-i0) / (double)(i1-i0));
	}
	return analysis_get_f0_ptr(a, channel, frame_index);
}

static double** get_spectrogram_ptr_edited(struct analysis* a, struct edit* e, int channel, int frame_index)
{
	int first_broken, last_broken;
	if (get_broken_info(a, e, frame_index, BROKEN_WIDTH_SECONDS_SPECTROGRAM, &first_broken, &last_broken)) {
		const int i0 = first_broken-1;
		const int i1 = last_broken+1;
		double* spectrogram0 = analysis_get_spectrogram(a, channel, i0);
		double* spectrogram1 = analysis_get_spectrogram(a, channel, i1);
		return fft_lerp(analysis_get_fft_storage_length(a), spectrogram0, spectrogram1, (double)(frame_index-i0) / (double)(i1-i0));
	}
	return analysis_get_spectrogram_ptr(a, channel, frame_index);
}

static double** get_aperiodicity_ptr_edited(struct analysis* a, struct edit* e, int channel, int frame_index)
{
	int first_broken, last_broken;
	if (get_broken_info(a, e, frame_index, BROKEN_WIDTH_SECONDS_APERIODICITY, &first_broken, &last_broken)) {
		const int i0 = first_broken-1;
		const int i1 = last_broken+1;
		double* aperiodicity0 = analysis_get_aperiodicity(a, channel, i0);
		double* aperiodicity1 = analysis_get_aperiodicity(a, channel, i1);
		return fft_lerp(analysis_get_fft_storage_length(a), aperiodicity0, aperiodicity1, (double)(frame_index-i0) / (double)(i1-i0));
	}
	return analysis_get_aperiodicity_ptr(a, channel, frame_index);
}


#define PLAY_EDITED (1)
#define PLAY_DRY   (2)


#define ARRAY_LENGTH(xs) (sizeof(xs) / sizeof(xs[0]))


struct audio {
	SDL_AudioDeviceID device_id;
	struct analysis* analysis;
	struct edit* edit;
	int n_channels;
	int is_playing;
	int position_frame_index;

	WorldSynthesizer synth[2];
} audio;


static void audio_render_next(int n_channels, int n_pcm_frames, float* out)
{
	assert(n_channels <= ARRAY_LENGTH(audio.synth));

	int pos[ARRAY_LENGTH(audio.synth)];

	const int n_analysis_frames = audio.analysis->n_analysis_frames;

	for (int ch = 0; ch < n_channels; ch++) {
		pos[ch] = audio.position_frame_index;
		WorldSynthesizer* synth = &audio.synth[ch];

		for (;;) {
			if (IsLocked(synth)) {
				printf("stopping audio: synth locked?!\n");
				audio.is_playing = 0;
				break;
			}

			assert((synth->synthesized_sample % synth->buffer_size) == 0);

			int synth_success = Synthesis2(synth);
			if (synth_success) {
				float* p = out;
				for (int i = 0; i < n_pcm_frames; i++) {
					*p = synth->buffer[i];
					p += n_channels;
				}
				break;
			}

			int p = pos[ch];
			if (p >= n_analysis_frames) {
				printf("stopping audio: reached end\n");
				audio.is_playing = 0;
				break;
			}

			double* f0;
			double** spectrogram;
			double** aperiodicity;
			if (audio.is_playing == PLAY_EDITED) {
				f0 = get_f0_ptr_edited(audio.analysis, audio.edit, ch, p);
				spectrogram = get_spectrogram_ptr_edited(audio.analysis, audio.edit, ch, p);
				aperiodicity = get_aperiodicity_ptr_edited(audio.analysis, audio.edit, ch, p);
			} else if (audio.is_playing == PLAY_DRY) {
				f0 = analysis_get_f0_ptr(audio.analysis, ch, p);
				spectrogram = analysis_get_spectrogram_ptr(audio.analysis, ch, p);
				aperiodicity = analysis_get_aperiodicity_ptr(audio.analysis, ch, p);
			} else {
				assert(!"playing what?");
			}

			int was_added = AddParameters(
				f0,
				1,
				spectrogram,
				aperiodicity,
				synth);
			if (was_added) {
				pos[ch]++;
			} else {
				printf("stopping audio: ringbuffer full?\n");
				audio.is_playing = 0;
				break;
			}
		}
	}

	for (int ch = 1; ch < n_channels; ch++) {
		if (pos[0] != pos[ch]) {
			printf("desync %d ?!\n", pos[0]-pos[ch]);
			#if 0
			printf("stopping audio: position desync? %d vs %d\n", pos[0], pos[ch]);
			audio.is_playing = 0;
			break;
			#endif
		}
	}

	audio.position_frame_index = pos[0];
}

static void audio_callback(void* usr, uint8_t* out, int n_bytes)
{
	const int n_channels = audio.n_channels;
	const int n = n_bytes / (n_channels * sizeof(float));

	if (audio.is_playing) {
		audio_render_next(n_channels, n, (float*)out);
	}

	if (!audio.is_playing) {
		float* p = (float*)out;
		const int cn = n_channels * n;
		for (int i = 0; i < cn; i++) {
			*(p++) = 0.0f;
		}
	}
}

static void audio_init(struct analysis* analysis, struct edit* edit)
{
	SDL_AudioSpec want = {
		.freq = analysis->sample_rate,
		.format = AUDIO_F32,
		.channels = analysis->n_channels,
		.samples = 1<<10,
		.callback = audio_callback,
	};

	SDL_AudioSpec have = {0};

	audio.device_id = SDL_OpenAudioDevice(NULL, 0, &want, &have, 0);
	assert((audio.device_id > 0) && "could not open audio device");

	assert(have.freq == want.freq);
	assert(have.format == want.format);
	assert(have.channels == want.channels);

	audio.n_channels = analysis->n_channels;
	audio.analysis = analysis;
	audio.edit = edit;

	for (int i = 0; i < ARRAY_LENGTH(audio.synth); i++) {
		InitializeSynthesizer(
			analysis->sample_rate,
			ANALYSIS_PERIOD_SECONDS * 1e3,
			get_fft_length_from_sample_rate(analysis->sample_rate),
			have.samples,
			32, // ringbuffer size
			&audio.synth[i]);
	}

	SDL_PauseAudioDevice(audio.device_id, 0);
}

static inline void audio_lock()
{
	SDL_LockAudioDevice(audio.device_id);
}

static inline void audio_unlock()
{
	SDL_UnlockAudioDevice(audio.device_id);
}

static int audio_get_play_position()
{
	audio_lock();
	int pos = -1;
	if (audio.is_playing) {
		pos = audio.position_frame_index;
	}
	audio_unlock();
	return pos;
}

static void audio_stop()
{
	audio_lock();
	audio.is_playing = 0;
	audio_unlock();
}

static void audio_play(int what, int frame0)
{
	audio_lock();

	audio.is_playing = what;
	audio.position_frame_index = frame0;

	const int n_channels = audio.analysis->n_channels;
	for (int i = 0; i < n_channels; i++) {
		RefreshSynthesizer(&audio.synth[i]);
	}

	audio_unlock();
}

static void lerp_init(int fft_size)
{
	for (int i = 0; i < RINGBUF_SZ; i++) {
		fft_ringbuf_pointers[i] = calloc(fft_size, sizeof **fft_ringbuf_pointers);
	}
}

static void audio_render(const char* path, struct analysis* a, struct edit* e)
{
	audio_lock();

	assert((a->n_channels == 1) && "stereo doesn't work...");

	drwav wav = {0};
	drwav_data_format format = {
		.container = drwav_container_riff,
		.format = DR_WAVE_FORMAT_PCM,
		.channels = a->n_channels,
		.sampleRate = a->sample_rate,
		.bitsPerSample = 16,
	};
	assert(drwav_init_file_write(&wav, path, &format, NULL));

	const int buffer_size = 1<<10;
	float* float_buffer = calloc(buffer_size, sizeof *float_buffer);
	drwav_int16* s16_buffer = calloc(buffer_size, sizeof* s16_buffer);

	const int copy_of_position_frame_index = audio.position_frame_index;
	const int copy_of_is_playing = audio.is_playing;

	audio.position_frame_index = 0;
	audio.is_playing = PLAY_EDITED;

	while (audio.is_playing) {
		audio_render_next(1, buffer_size, float_buffer);
		drwav_f32_to_s16(s16_buffer, float_buffer, buffer_size);
		drwav_uint64 n_written = drwav_write_pcm_frames(&wav, buffer_size, s16_buffer);
	}

	drwav_uninit(&wav);

	audio.is_playing = copy_of_is_playing;
	audio.position_frame_index = copy_of_position_frame_index;

	audio_unlock();
}

int main(int argc, char** argv)
{
	assert(((sizeof(struct analysis_header)&7)==0) && "should prolly be 8-byte aligned, no bro?");

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <path/to/audio.wav>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	const char* path_wav = argv[1];
	char path_analysis[1024];
	char path_edit[1024];
	char path_render[1024];
	snprintf(path_analysis, sizeof path_analysis, "%s.analysis", path_wav);
	snprintf(path_edit, sizeof path_edit, "%s.edit", path_wav);
	snprintf(path_render, sizeof path_render, "%s.render.wav", path_wav);

	printf("\npaths:\n");
	printf("   wav:       %s\n", path_wav);
	printf("   analysis:  %s\n", path_analysis);
	printf("   edit:      %s\n", path_edit);
	printf("   render:    %s\n", path_render);

	uint32_t n_channels;
	uint32_t sample_rate;
	drwav_uint64 n_pcm_frames;
	float* samples = drwav_open_file_and_read_pcm_frames_f32(path_wav, &n_channels, &sample_rate, &n_pcm_frames, NULL);
	if (samples == NULL) {
		fprintf(stderr, "%s: could not open\n", path_wav);
		exit(EXIT_FAILURE);
	}

	printf("\n%s:\n", path_wav);
	printf("   number of channels:  %u\n", n_channels);
	printf("   sample_rate:         %u\n", sample_rate);
	printf("   number of frames:    %llu\n", n_pcm_frames);
	printf("   duration:            %.1fs\n", (double)n_pcm_frames / (double)sample_rate);

	struct analysis analysis;
	if (!analysis_open(&analysis, path_analysis)) {
		const int n_analysis_frames = calc_number_of_analysis_frames(sample_rate, n_pcm_frames);
		const size_t sz = calc_analysis_file_size(n_channels, n_analysis_frames, sample_rate);

		printf("Please wait while generating analysis file (%s; %zdB / %.2fGB) ...\n", path_analysis, sz, (double)sz/(1024.0*1024.0*1024.0));

		int fd = open(path_analysis, O_RDWR | O_CREAT | O_EXCL, 0644);
		if (fd == -1) {
			fprintf(stderr, "Could not create ''%s'': %s\n", path_analysis, strerror(errno));
			exit(EXIT_FAILURE);
		}

		if (posix_fallocate(fd, 0, sz) != 0) {
			fprintf(stderr, "%s: posix_fallocate() failed\n", path_analysis);
			exit(EXIT_FAILURE);
		}

		void* file_data = mmap(NULL, sz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		if (file_data == MAP_FAILED) {
			fprintf(stderr, "mmap on ''%s'' failed: %s\n", path_analysis, strerror(errno));
			exit(EXIT_FAILURE);
		}

		void* p = file_data;

		{
			struct analysis_header h = {
				.magic = ANALYSIS_MAGIC,
				.n_channels = n_channels,
				.sample_rate = sample_rate,
				.n_analysis_frames = n_analysis_frames,
			};
			memcpy(p, &h, sizeof h);
			p += sizeof h;
		}

		double* fp = p;

		double* xs = calloc(n_pcm_frames, sizeof *xs);
		double* time_axis = calloc(n_pcm_frames, sizeof *xs);

		double** fft_dst = calloc(n_analysis_frames, sizeof *fft_dst);

		ForwardRealFFT forward_real_fft = {0};
		InitializeForwardRealFFT(OTHER_SPECTROGRAM_FFT_LENGTH, &forward_real_fft);

		{
			printf(" ... other spectrogram/DFT (mono)\n");

			const double kb_att = 60.0;
			const double kb_alpha = kaiser_bessel_attenuation_to_alpha(kb_att);

			const double ch_scale = 1.0 / (double)n_channels;

			for (int i0 = 0; i0 < n_analysis_frames; i0++) {
				for (int i1 = 0; i1 < OTHER_SPECTROGRAM_FFT_LENGTH; i1++) {
					double t = ((double)i1 * 2.0) / (double)(OTHER_SPECTROGRAM_FFT_LENGTH) - 1.0;
					double w = kaiser_bessel(kb_alpha, t); // could put this expensive func into a LUT...
					int xi = (int)round((double)(i0*sample_rate)*ANALYSIS_PERIOD_SECONDS) + i1 - OTHER_SPECTROGRAM_FFT_LENGTH/2;

					double x = 0.0;
					if (0 <= xi && xi < n_pcm_frames) {
						for (int i2 = 0; i2 < n_channels; i2++) {
							x += samples[xi*n_channels + i2];
						}
					}
					x *= ch_scale;

					forward_real_fft.waveform[i1] = x*w;
				}
				fft_execute(forward_real_fft.forward_fft);

				for (int i1 = 0; i1 < OTHER_SPECTROGRAM_FFT_STORAGE_LENGTH; i1++) {
					double re = forward_real_fft.spectrum[i1][0];
					double im = forward_real_fft.spectrum[i1][1];
					double mag = sqrt(re*re + im*im);
					*(fp++) = mag;
				}
			}
		}

		DestroyForwardRealFFT(&forward_real_fft);

		const int fft_length = get_fft_length_from_sample_rate(sample_rate);
		const int fft_storage_length = calc_fft_storage_length_from_fft_length(fft_length);

		for (int i0 = 0; i0 < n_channels; i0++) {
			printf(" ... channel %d\n", i0);
			for (int i1 = 0; i1 < n_pcm_frames; i1++) {
				xs[i1] = samples[i1*n_channels + i0];
			}

			double* f0 = fp;

			{
				printf(" ... f0/Harvest\n");

				HarvestOption option = { 0 };
				InitializeHarvestOption(&option);
				option.frame_period = ANALYSIS_PERIOD_SECONDS*1e3; // ms
				option.f0_floor = 40.0;

				Harvest(xs, n_pcm_frames, sample_rate, &option, time_axis, f0);
				fp += n_analysis_frames;
			}

			#define PREP_FFT_DST \
				{ \
					for (int ii_ = 0; ii_ < n_analysis_frames; ii_++) { \
						fft_dst[ii_] = fp; \
						fp += fft_storage_length; \
					} \
				}

			{
				printf(" ... spectrogram/CheapTrick\n");
				CheapTrickOption option = { 0 };
				InitializeCheapTrickOption(sample_rate, &option);
				option.f0_floor = ANALYSIS_F0_FLOOR;
				option.fft_size = fft_length;
				PREP_FFT_DST
				CheapTrick(xs, n_pcm_frames, sample_rate, time_axis, f0, n_analysis_frames, &option, fft_dst);
			}

			{
				printf(" ... aperiodicity/D4C\n");
				D4COption option = { 0 };
				InitializeD4COption(&option);
				option.threshold = 0.85;
				PREP_FFT_DST
				D4C(xs, n_pcm_frames, sample_rate, time_axis, f0, n_analysis_frames, fft_length, &option, fft_dst);
			}

			#undef PREP_FFT_DST
		}

		assert(((void*)fp - file_data) == sz);

		printf(" ... done!\n");

		free(fft_dst);
		free(time_axis);
		free(xs);

		munmap(file_data, sz);
		close(fd);

		assert(analysis_open(&analysis, path_analysis));
	}

	struct edit edit;
	edit_open(&edit, path_edit, &analysis);

	lerp_init(analysis_get_fft_storage_length(&analysis));

	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO);

	audio_init(&analysis, &edit);

	SDL_Window* window;
	SDL_Renderer* renderer;
	SDL_CreateWindowAndRenderer(1920, 1080, SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI, &window, &renderer);
	{
		char title[1024];
		snprintf(title, sizeof title, "%s (voxpopfix)", path_wav);
		SDL_SetWindowTitle(window, title);
	}

	SDL_Texture* canvas_texture = NULL;
	int canvas_width = -1;
	int canvas_height = -1;
	uint32_t* canvas_pixels = NULL;
	size_t canvas_sz;

	int exiting = 0;

	struct ui ui;
	ui_init(&ui);

	while (!exiting) {
		float d_zoom = 0.0f;
		int play = 0;

		SDL_Event e;
		while (SDL_PollEvent(&e)) {
			switch (e.type) {

			case SDL_QUIT: {
				exiting = 1;
			} break;

			case SDL_KEYDOWN:
			case SDL_KEYUP: {
				const int is_down = (e.type == SDL_KEYDOWN);
				const int sym = e.key.keysym.sym;
				if (is_down) {
					if (sym == SDLK_ESCAPE) {
						exiting = 1;
					} else if (sym == ' ') {
						play = PLAY_EDITED;
					} else if (sym == 'd') {
						play = PLAY_DRY;
					} else if (sym == 'r') {
						audio_render(path_render, &analysis, &edit);
					}
				}

				if (is_down && '1' <= sym && sym < '1'+N_UI_ROWS) {
					int i = sym-'1';
					ui.show[i] = !ui.show[i];
				}
			} break;

			case SDL_MOUSEWHEEL: {
				d_zoom += e.wheel.preciseY;
			} break;

			}
		}

		const int play_position = audio_get_play_position();
		const int is_playing = (play_position >= 0);

		int mx;
		int my;
		const int btns = SDL_GetMouseState(&mx, &my);
		const int mod_state = SDL_GetModState();
		ui_input(&ui, mx/DIVRES, my/DIVRES, btns, mod_state, d_zoom, play_position, &edit);
		const int mf = ui.mouseover_frame_index;

		if (is_playing && play) {
			audio_stop();
		} else if (play) {
			audio_play(play, mf);
		}

		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
		SDL_RenderClear(renderer);

		int screen_width;
		int screen_height;
		SDL_GL_GetDrawableSize(window, &screen_width, &screen_height);

		if (canvas_texture == NULL || canvas_width != screen_width/DIVRES || canvas_height != screen_height/DIVRES) {
			if (canvas_texture != NULL) SDL_DestroyTexture(canvas_texture);
			if (canvas_pixels != NULL) free(canvas_pixels);
			canvas_width = screen_width/DIVRES;
			canvas_height = screen_height/DIVRES;
			canvas_texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, canvas_width, canvas_height);
			canvas_sz = sizeof(*canvas_pixels) * canvas_width * canvas_height;
			canvas_pixels = malloc(canvas_sz);
		}

		memset(canvas_pixels, 0, canvas_sz);

		ui_handle(&ui, canvas_width, canvas_height, canvas_pixels, &analysis, &edit);

		SDL_UpdateTexture(canvas_texture, NULL, canvas_pixels, canvas_width * sizeof(*canvas_pixels));
		SDL_RenderCopy(renderer, canvas_texture, NULL, NULL);

		SDL_RenderPresent(renderer);
	}

	edit_close(&edit);

	return EXIT_SUCCESS;
}
