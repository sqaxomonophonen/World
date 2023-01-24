#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#include "audioio.h"
#include "world/harvest.h"
#include "world/cheaptrick.h"
#include "world/d4c.h"

#include "analysis.h"

int main(int argc, char** argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <in.wav> <out.analysis>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	const char* in_wav = argv[1];

	int n_samples = GetAudioLength(in_wav);
	printf("n samples: %d\n", n_samples);

	double* xs = calloc(n_samples, sizeof *xs);

	int sample_rate;
	int n_bits;
	wavread(in_wav, &sample_rate, &n_bits, xs);
	printf("sample rate: %d\n", sample_rate);
	printf("n bits: %d\n", n_bits);


	const double frame_period_ms = 5.0;

	int n_frames;
	double* f0s;
	double* time_axis;
	{
		HarvestOption option = { 0 };
		InitializeHarvestOption(&option);
		option.frame_period = frame_period_ms;
		option.f0_floor = 40.0;

		n_frames = GetSamplesForHarvest(sample_rate, n_samples, frame_period_ms);
		printf("n frames: %d\n", n_frames);

		f0s = calloc(n_frames, sizeof *f0s);
		time_axis = calloc(n_frames, sizeof *time_axis);

		printf("Analysing (f0/Harvest) ...\n");
		Harvest(xs, n_samples, sample_rate, &option, time_axis, f0s);
	}

	int fft_size;
	int fft_size_storage;

	double** spectrogram;
	{
		CheapTrickOption option = { 0 };
		InitializeCheapTrickOption(sample_rate, &option);
		option.f0_floor = 71.0;
		fft_size = option.fft_size = GetFFTSizeForCheapTrick(sample_rate, &option);
		fft_size_storage = fft_size/2+1;
		printf("fft size: %d\n", fft_size);

		spectrogram = _alloc_doubles_2d(n_frames, fft_size_storage);

		printf("Analysing (spectral/CheapTrick) ...\n");
		CheapTrick(xs, n_samples, sample_rate, time_axis, f0s, n_frames, &option, spectrogram);
	}

	double** aperiodicity;
	{
		D4COption option = { 0 };
		InitializeD4COption(&option);
		option.threshold = 0.85;
		printf("Analysing (aperiodicity/D4C) ...\n");
		aperiodicity = _alloc_doubles_2d(n_frames, fft_size_storage);
		D4C(xs, n_samples, sample_rate, time_axis, f0s, n_frames, fft_size, &option, aperiodicity);
	}

	analysis_write(
		argv[2],
		sample_rate,
		n_frames,
		fft_size_storage,
		frame_period_ms,
		f0s,
		spectrogram,
		aperiodicity);

	printf("Done!\n");

	return EXIT_SUCCESS;
}
