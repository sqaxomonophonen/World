#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#include "world/synthesis.h"

#include "analysis.h"

int main(int argc, char** argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <in.analysis> <out.wav>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	struct analysis* a = analysis_read(argv[1]);

	int n_samples = analysis_get_n_output_samples(a);

	printf("n samples: %d\n", n_samples);
	double* samples = calloc(n_samples, sizeof *samples);

	printf("synthesising...\n");
	Synthesis(
		a->f0s,
		a->n_frames,
		a->spectrogram,
		a->aperiodicity,
		analysis_get_true_fft_size(a),
		a->frame_period_ms,
		a->sample_rate,
		n_samples,
		samples);

	printf("writing %s ...\n", argv[2]);
	wavwrite(samples, n_samples, a->sample_rate, 16, argv[2]);

	return EXIT_SUCCESS;
}
