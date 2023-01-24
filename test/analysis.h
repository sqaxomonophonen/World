#ifndef ANALYSIS_H

#include <assert.h>
#include <stdio.h>
#include <stdint.h>

#define ANALYSIS_MAGIC (0xa741cafe)

struct analysis {
	int sample_rate;
	int n_frames;
	int fft_size_storage;
	double frame_period_ms;
	double* f0s;
	double** spectrogram;
	double** aperiodicity;
};

static void _r32(FILE* f, void* dst)
{
	assert(fread(dst, 4, 1, f) == 1);
}

static void _r64(FILE* f, void* dst)
{
	assert(fread(dst, 8, 1, f) == 1);
}

static uint32_t _read_u32(FILE* f)
{
	uint32_t v;
	_r32(f, &v);
	return v;
}

static double _read_f64(FILE* f)
{
	double v;
	_r64(f, &v);
	return v;
}

static void _w32(FILE* f, void* p)
{
	assert(fwrite(p, 4, 1, f) == 1);
}

static void _w64(FILE* f, void* p)
{
	assert(fwrite(p, 8, 1, f) == 1);
}

static void _write_u32(FILE* f, uint32_t v)
{
	_w32(f, &v);
}

static void _write_f64(FILE* f, double v)
{
	_w64(f, &v);
}

static double** _alloc_doubles_2d(size_t n0, size_t n1)
{
	double** xs;
	xs = calloc(n0, sizeof *xs);
	for (int i = 0; i < n0; i++) xs[i] = calloc(n1, sizeof **xs);
	return xs;
}


static struct analysis* analysis_read(const char* path)
{
	FILE* f = fopen(path, "rb");

	uint32_t magic = _read_u32(f);
	assert(magic == ANALYSIS_MAGIC);

	struct analysis* a = calloc(1, sizeof *a);

	a->sample_rate = _read_u32(f);
	a->n_frames = _read_u32(f);
	a->fft_size_storage = _read_u32(f);
	a->frame_period_ms = _read_f64(f);

	a->f0s = calloc(a->n_frames, sizeof *a->f0s);
	for (int i = 0; i < a->n_frames; i++) a->f0s[i] = _read_f64(f);

	a->spectrogram = _alloc_doubles_2d(a->n_frames, a->fft_size_storage);
	for (int i = 0; i < a->n_frames; i++) for (int j = 0; j < a->fft_size_storage; j++) a->spectrogram[i][j] = _read_f64(f);

	a->aperiodicity = _alloc_doubles_2d(a->n_frames, a->fft_size_storage);
	for (int i = 0; i < a->n_frames; i++) for (int j = 0; j < a->fft_size_storage; j++) a->aperiodicity[i][j] = _read_f64(f);

	assert(ftell(f) == (4*4 + 8*1 + 8*a->n_frames*(1+2*a->fft_size_storage)));

	fclose(f);

	return a;
}

static void analysis_write(const char* path, int sample_rate, int n_frames, int fft_size_storage, double frame_period_ms, double* f0s, double** spectrogram, double** aperiodicity)
{
	FILE* f = fopen(path, "wb");

	_write_u32(f, ANALYSIS_MAGIC);
	_write_u32(f, sample_rate);
	_write_u32(f, n_frames);
	_write_u32(f, fft_size_storage);
	_write_f64(f, frame_period_ms);

	for (int i = 0; i < n_frames; i++) _write_f64(f, f0s[i]);
	for (int i = 0; i < n_frames; i++) for (int j = 0; j < fft_size_storage; j++) _write_f64(f, spectrogram[i][j]);
	for (int i = 0; i < n_frames; i++) for (int j = 0; j < fft_size_storage; j++) _write_f64(f, aperiodicity[i][j]);

	fclose(f);
}

static void analysis_rewrite(struct analysis* a, const char* path)
{
	analysis_write(path, a->sample_rate, a->n_frames, a->fft_size_storage, a->frame_period_ms, a->f0s, a->spectrogram, a->aperiodicity);
}

static int analysis_get_n_output_samples(struct analysis* a)
{
	return ((a->n_frames - 1) * (a->frame_period_ms/1000.0) * a->sample_rate) + 1;
}

static int analysis_get_true_fft_size(struct analysis* a)
{
	return (a->fft_size_storage-1)*2;
}

static int analysis_get_f0_min_max(struct analysis* a, double* f0_min, double* f0_max)
{
	double vmin, vmax;
	int first = 1;
	for (int i = 0; i < a->n_frames; i++) {
		double f0 = a->f0s[i];
		if (f0 == 0.0) continue;
		assert(f0 > 0.0);
		if (first) {
			vmin = vmax = f0;
			first = 0;
		} else {
			if (f0 < vmin) vmin = f0;
			if (f0 > vmax) vmax = f0;
		}
	}
	if (f0_min) *f0_min = vmin;
	if (f0_max) *f0_max = vmax;
}

#define ANALYSIS_H
#endif
