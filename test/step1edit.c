#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#include <SDL.h>

#include "analysis.h"

static inline double lerp(double a, double b, double t) {
	return a + (b-a)*t;
}

int main(int argc, char** argv)
{
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <in.analysis>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	struct analysis* a = analysis_read(argv[1]);

	double f0_min, f0_max;
	analysis_get_f0_min_max(a, &f0_min, &f0_max);

	printf("f0: [%f;%f]\n", f0_min, f0_max);

	SDL_Window* window;
	SDL_Renderer* renderer;
	SDL_CreateWindowAndRenderer(1920, 1080, SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI, &window, &renderer);
	SDL_SetWindowTitle(window, "step1edit");

	SDL_Texture* canvas_texture = NULL;
	int canvas_width = -1;
	int canvas_height;
	uint32_t* canvas_pixels = NULL;
	size_t canvas_sz;

	int frame = 0;

	int exiting = 0;
	int cursor = 0;
	int selection0 = 0;
	int selection1 = 0;
	int zoom = 2;

	int left = 0;
	int right = 0;
	int vx = 0;
	int speed_mode = 1;

	while (!exiting) {
		SDL_Event e;
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				exiting = 1;
			} else if (e.type == SDL_KEYDOWN) {
				int sym = e.key.keysym.sym;
				if (sym == SDLK_ESCAPE) {
					exiting = 1;
				} else if (sym == SDLK_TAB) {
					speed_mode = !speed_mode;
				} else if (sym == SDLK_LEFT) {
					if (speed_mode) {
						left = 1;
					} else {
						cursor--;
					}
				} else if (sym == SDLK_RIGHT) {
					if (speed_mode) {
						right = 1;
					} else {
						cursor++;
					}
				} else if ('1' <= sym && sym <= '9') {
					zoom = (sym-SDLK_1)+1;
				} else if (sym == '[') {
					selection0 = cursor;
				} else if (sym == ']') {
					selection1 = cursor;
				} else if (sym == SDLK_RETURN) {
					double f0 = a->f0s[selection0];
					double f1 = a->f0s[selection1];
					for (int i = selection0; i <= selection1; i++) {
						double t = (double)(i-selection0) / (double)(selection1-selection0);
						a->f0s[i] = lerp(f0,f1,t);
						for (int what = 0; what < 2; what++) {
							double** data = what == 0 ? a->spectrogram : what == 1 ? a->aperiodicity : 0;
							double* data0 = data[selection0];
							double* data1 = data[selection1];
							double* datai = data[i];
							for (int j = 0; j < a->fft_size_storage; j++) {
								datai[j] = lerp(data0[j], data1[j], t);
							}
						}
					}
				} else if (sym == 's') {
					analysis_rewrite(a, "step1.rewrite.analysis");
				}
			} else if (e.type == SDL_KEYUP) {
				int sym = e.key.keysym.sym;
				if (sym == SDLK_LEFT) {
					left = 0;
				} else if (sym == SDLK_RIGHT) {
					right = 0;
				}
			}
		}

		if (left && !right) {
			vx--;
		} else if (right && !left) {
			vx++;
		} else {
			vx = 0;
		}
		cursor += vx;

		if (cursor < 0) cursor = 0;
		if (cursor >= a->n_frames) cursor = a->n_frames-1;

		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
		SDL_RenderClear(renderer);

		int screen_width;
		int screen_height;
		SDL_GL_GetDrawableSize(window, &screen_width, &screen_height);

		if (canvas_texture == NULL || canvas_width != screen_width || canvas_height != screen_height) {
			if (canvas_texture != NULL) SDL_DestroyTexture(canvas_texture);
			if (canvas_pixels != NULL) free(canvas_pixels);
			canvas_width = screen_width;
			canvas_height = screen_height;
			canvas_texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, canvas_width, canvas_height);
			canvas_sz = sizeof(*canvas_pixels) * canvas_width * canvas_height;
			canvas_pixels = malloc(canvas_sz);
		}

		memset(canvas_pixels, 0, canvas_sz);

		int c3h = canvas_height / 3;

		#define PP(X,Y,R,G,B) \
			if (0 <= (X) && (X) < canvas_width && 0 <= (Y) && (Y) <= canvas_height) { \
				canvas_pixels[(X)+(Y)*canvas_width] = ((R)<<24) + ((G)<<16) + ((B)<<8) + ((0xff)<<0); \
			}

		#define PPF(X,Y,R,G,B) \
			do { \
				int _r = round((R)*255.0); \
				int _g = round((G)*255.0); \
				int _b = round((B)*255.0); \
				if (_r < 0)   _r = 0; \
				if (_r > 255) _r = 255; \
				if (_g < 0)   _g = 0; \
				if (_g > 255) _g = 255; \
				if (_b < 0)   _b = 0; \
				if (_b > 255) _b = 255; \
				PP(X,Y,_r,_g,_b) \
			} while (0);

		{
			int i0 = cursor;
			int x0 = canvas_width/2 - zoom/2;

			for (;;) {
				if (i0 <= 0 || x0 <= 0) break;
				x0 -= zoom;
				i0--;
			}

			{
				int cx = x0;
				int i = i0;
				for (;;) {
					if (cx >= canvas_width || i >= a->n_frames) break;
					assert(0 <= i && i < a->n_frames);
					double f0 = a->f0s[i];

					double pr = 0;
					double pg = 0;
					double pb = 0;
					if (i == cursor) {
						pb += 1.0;
					}
					if (selection0 <= i && i <= selection1) {
						pg += 0.3;
					}

					if (pr>0 || pg>0 || pb>0) {
						for (int y = 0; y < c3h; y++) {
							for (int x = cx; x < (cx+zoom); x++) {
								PPF(x,y,pr,pg,pb)
							}
						}
					}

					if (f0 > 0) {
						const int y = 1*c3h - ((f0-f0_min) / (f0_max-f0_min)) * c3h;
						for (int x = cx; x < (cx+zoom); x++) {
							PP(x,y,255,255,255)
						}
					}

					cx += zoom;
					i++;
				}
			}

			for (int what = 0; what < 2; what++) {
				int cx = x0;
				int i = i0;
				for (;;) {
					if (cx >= canvas_width || i >= a->n_frames) break;
					assert(0 <= i && i < a->n_frames);
					double* data = what == 0 ? a->spectrogram[i] : what == 1 ? a->aperiodicity[i] : NULL;
					assert(data);


					double pr = 0;
					double pg = 0;
					double pb = 0;
					if (i == cursor) {
						pr += -1.0;
						pb += 1.0;
					}
					if (selection0 <= i && i <= selection1) {
						pg += 0.3;
					}

					const int y0 = (what+1)*c3h;
					const int y1 = (what+2)*c3h;
					for (int y = y0; y < y1; y++) {
						double t = (double)(y - y0) / (double)(y1-y0);
						int fi = t * a->fft_size_storage;
						if (fi < 0) fi = 0;
						if (fi >= a->fft_size_storage) fi = a->fft_size_storage-1;
						double v = data[a->fft_size_storage-1-fi];

						if (what == 0) {
							if (v > 0) v = log(v+1.0) * 2000.0;
						}
						for (int x = cx; x < (cx+zoom); x++) {
							double r = v;
							double g = v*v;
							double b = v*v*v;
							r += pr; g += pg; b += pb;
							PPF(x,y,r,g,b)
						}
					}

					cx += zoom;
					i++;
				}
			}
		}

		SDL_UpdateTexture(canvas_texture, NULL, canvas_pixels, canvas_width * sizeof(*canvas_pixels));
		SDL_RenderCopy(renderer, canvas_texture, NULL, NULL);

		SDL_RenderPresent(renderer);

		frame++;
	}

	SDL_DestroyWindow(window);

	return EXIT_SUCCESS;
}
