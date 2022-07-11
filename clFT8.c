#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include<ctype.h>

#include <alsa/asoundlib.h>
#include <pthread.h>
#include <sys/time.h>

#include "ft8/unpack.h"
#include "ft8/ldpc.h"
#include "ft8/decode.h"
#include "ft8/constants.h"
#include "ft8/encode.h"
#include "ft8/crc.h"

#include "fft/kiss_fftr.h"

#include "clFT8.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Global

/* Global FT8 info. */
typedef struct
{
	char Local_CALLSIGN[20];
	char Local_LOCATOR[20];
	float Local_latlon[2];
	bool TX_enable;
	
	char QSO_dist_CALLSIGN[20];
	char QSO_dist_LOCATOR[20];
	char QSO_dist_MESSAGE[20];
	int  QSO_dist_SNR;
	
	char QSO_RESPONSES[6][50];
	int QSO_Index_to_rep;
	
	bool TRX_status;
	pthread_mutex_t TRX_status_lock;
	pthread_cond_t RX_status_cond;
	pthread_cond_t TX_status_cond;
	
	
} FT8info;

FT8info FT8 = {
	.Local_CALLSIGN={'F','4','J','J','J','\0'},
	.Local_LOCATOR={'J','N','3','8','t','g','\0'},
	.TX_enable=1,
	
	.QSO_dist_CALLSIGN="",
	.QSO_dist_LOCATOR="",
	.QSO_dist_MESSAGE="",
	.QSO_dist_SNR=0,
	
	
	.QSO_RESPONSES={'\0','\0','\0','\0','\0','\0'},
	.QSO_Index_to_rep=-1,
	
	.TRX_status = _RX_,
	.TRX_status_lock = PTHREAD_MUTEX_INITIALIZER,
	.RX_status_cond = PTHREAD_COND_INITIALIZER,
	.TX_status_cond = PTHREAD_COND_INITIALIZER
	
	
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Divers
/* Return the environment variable "name" or "def" if it's unset. */
char *getenvDefault(char *name, char *def)
{
	char *val = getenv(name);
	if (val == NULL)
		return def;
	else
		return val;
}

/* Get time now */
double now()
{
	struct timeval tv;
	gettimeofday(&tv, 0);
	return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void waitforstart(){
	double tt = now();
	long long cycle_start = tt - ((long long)tt % 15);
	double WaitTime = (15 - (double)(tt - cycle_start))*1000000;
	usleep(WaitTime);
}

void printDateTime(){
	time_t t = time(NULL);
	struct tm tm = *gmtime(&t);
	int sec = (int)((int)tm.tm_sec / 15)*15;
	printf("%d-%02d-%02d %02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, sec);
}

void latLonForGrid(char * grid, float * latlon) {
	
	latlon[0]=-1;
	latlon[1]=-1;
	
	if(strlen(grid)==4){
		if(isalpha(grid[0]) && isalpha(grid[1]) && isdigit(grid[2]) && isdigit(grid[3])){
			float str_chr_up(char g){
				char strtable[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
				const char *moved_string = strchr(strtable, g);
				if (moved_string) {
					return moved_string - strtable;
				}
				return -1;
			}
			
			float str_num(char g){
				char strtable[] = "0123456789";
				const char *moved_string = strchr(strtable, g);
				if (moved_string) {
					return moved_string - strtable;
				}
				return -1;
			}
			latlon[0] = str_chr_up(grid[1]) * 10;               // 2nd digit: 10deg latitude slot.
			latlon[1] = str_chr_up(grid[0]) * 20;               // 1st digit: 20deg longitude slot.
			latlon[0] += str_num(grid[3]) * 1;                  // 4th digit: 1deg latitude slot.
			latlon[1] += str_num(grid[2]) * 2;                  // 3rd digit: 2deg longitude slot.
			latlon[0] += 0.5 * 1;
			latlon[1] += 0.5 * 2;
			latlon[0] -= 90;                                                  // Locator lat/lon origin shift.
			latlon[1] -= 180;
		}
	}
};


float latLonDist(float * latlonA, float * latlonB){
	return acos(sin(	(latlonA[0] * M_PI / 180.0f)	)*sin(	(latlonA[1] * M_PI / 180.0f)	)+cos(	(latlonA[0] * M_PI / 180.0f)	)*cos(	(latlonA[1] * M_PI / 180.0f)	)*cos(	((latlonB[0]-latlonB[1]) * M_PI / 180.0f)	))*6371. ;
}

void unpackFT8mess(char * message_text, char * unpackeds0, char * unpackeds1, char * unpackeds2){
	char * NullToken = strtok(message_text, " ");
	if(NullToken != 0){
		memcpy(unpackeds0, NullToken, sizeof(NullToken));
		// strcpy(unpackeds0, NullToken);
	}
	
	NullToken = strtok(0, " ");	
	if(NullToken != 0){
		memcpy(unpackeds1, NullToken, sizeof(NullToken));
		// strcpy(unpackeds1, NullToken);
	}
	
	NullToken = strtok(0, " ");		
	if(NullToken != 0){
		memcpy(unpackeds2, NullToken, sizeof(NullToken));
		// strcpy(unpackeds2, NullToken);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Audio

/* Global sound info. */
typedef struct
{
	char *capture_sound_device;
	unsigned int capture_sound_rate;
	snd_pcm_t *capture_handle;
	
	char *playback_sound_device;
	int playback_buffer_frames;
	unsigned int playback_sound_rate;
	unsigned int playback_num_of_samples, playback_buffer_size;
	snd_pcm_t *playback_handle;
	char *playback_buffer;

} soundInfo;


soundInfo sound={
	.capture_sound_device = (char*)"default",
	.capture_sound_rate = 12000,
	
	.playback_sound_device = (char*)"default",
	.playback_buffer_frames = 1024,
	.playback_sound_rate = 12000
};

/* Open and init the default recording device. */
void capture_audioInit(void)
{ 
	int i;
	int err;
	snd_pcm_hw_params_t *capture_hw_params;

	snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

	if ((err = snd_pcm_open (&sound.capture_handle, getenvDefault(capture_sound_device_ENV, sound.capture_sound_device), SND_PCM_STREAM_CAPTURE, 0)) < 0) {
		fprintf (stderr, "cannot open audio device %s (%s)\n", sound.capture_sound_device, snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "audio interface opened\n");
	#endif
	
	if ((err = snd_pcm_hw_params_malloc (&capture_hw_params)) < 0) {
		fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "capture_hw_params allocated\n");
	#endif
				 
	if ((err = snd_pcm_hw_params_any (sound.capture_handle, capture_hw_params)) < 0) {
		fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "capture_hw_params initialized\n");
	#endif
	
	if ((err = snd_pcm_hw_params_set_access (sound.capture_handle, capture_hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
		fprintf (stderr, "cannot set access type (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "capture_hw_params access setted\n");
	#endif
		
	if ((err = snd_pcm_hw_params_set_format (sound.capture_handle, capture_hw_params, format)) < 0) {
		fprintf (stderr, "cannot set sample format (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "capture_hw_params format setted\n");
	#endif
		
	if ((err = snd_pcm_hw_params_set_rate_near (sound.capture_handle, capture_hw_params, &sound.capture_sound_rate, 0)) < 0) {
		fprintf (stderr, "cannot set sample rate (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG	
	fprintf(stdout, "capture_hw_params rate setted\n");
	#endif
	
	if ((err = snd_pcm_hw_params_set_channels (sound.capture_handle, capture_hw_params, 1)) < 0) {
		fprintf (stderr, "cannot set channel count (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "capture_hw_params channels setted\n");
	#endif
		
	if ((err = snd_pcm_hw_params (sound.capture_handle, capture_hw_params)) < 0) {
		fprintf (stderr, "cannot set parameters (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "capture_hw_params setted\n");
	#endif
		
	snd_pcm_hw_params_free (capture_hw_params);

	#if DEBUG
	fprintf(stdout, "capture_hw_params freed\n");
	#endif
		
	if ((err = snd_pcm_prepare (sound.capture_handle)) < 0) {
		fprintf (stderr, "cannot prepare audio interface for use (%s)\n", snd_strerror (err));
		exit (1);
	}

	#if DEBUG
	fprintf(stdout, "audio interface prepared\n\n");
	#endif
	
}

void capture_audioDeInit(void)
{ 
	snd_pcm_close (sound.capture_handle);
}

/* Open and init the default playback device. */
void playback_audioInit(void)
{ 
	int i;
	int err;
	snd_pcm_hw_params_t *playback_hw_params;

	snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

	if ((err = snd_pcm_open (&sound.playback_handle, getenvDefault(playback_sound_device_ENV, sound.playback_sound_device), SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
		fprintf (stderr, "cannot open audio device %s (%s)\n", sound.playback_sound_device, snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "audio interface opened\n");
	#endif
				 
	if ((err = snd_pcm_hw_params_malloc (&playback_hw_params)) < 0) {
		fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params allocated\n");
	#endif
					 
	if ((err = snd_pcm_hw_params_any (sound.playback_handle, playback_hw_params)) < 0) {
		fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params initialized\n");
	#endif
		
	if ((err = snd_pcm_hw_params_set_access (sound.playback_handle, playback_hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
		fprintf (stderr, "cannot set access type (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params access setted\n");
	#endif
		
	if ((err = snd_pcm_hw_params_set_format (sound.playback_handle, playback_hw_params, format)) < 0) {
		fprintf (stderr, "cannot set sample format (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params format setted\n");
	#endif
		
	if ((err = snd_pcm_hw_params_set_rate_near (sound.playback_handle, playback_hw_params, &sound.playback_sound_rate, 0)) < 0) {
		fprintf (stderr, "cannot set sample rate (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG	
	fprintf(stdout, "playback_hw_params rate setted\n");
	#endif
	
	if ((err = snd_pcm_hw_params_set_channels (sound.playback_handle, playback_hw_params, 1)) < 0) {
		fprintf (stderr, "cannot set channel count (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params channels setted\n");
	#endif
		
	if ((err = snd_pcm_hw_params (sound.playback_handle, playback_hw_params)) < 0) {
		fprintf (stderr, "cannot set parameters (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params setted\n");
	#endif
		
	snd_pcm_hw_params_free (playback_hw_params);
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params freed\n");
	#endif
		
	if ((err = snd_pcm_prepare (sound.playback_handle)) < 0) {
		fprintf (stderr, "cannot prepare audio interface for use (%s)\n", snd_strerror (err));
		exit (1);
	}
	
	#if DEBUG
	fprintf(stdout, "audio interface prepared\n\n");
	#endif
}

void playback_audioDeInit(void)
{ 
	free(sound.playback_buffer);
	snd_pcm_close (sound.playback_handle);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FT8_lib encode

#define FT8_SYMBOL_BT 2.0f ///< symbol smoothing filter bandwidth factor (BT)
#define FT4_SYMBOL_BT 1.0f ///< symbol smoothing filter bandwidth factor (BT)

#define GFSK_CONST_K 5.336446f ///< == pi * sqrt(2 / log(2))

/// Computes a GFSK smoothing pulse.
/// The pulse is theoretically infinitely long, however, here it's truncated at 3 times the symbol length.
/// This means the pulse array has to have space for 3*n_spsym elements.
/// @param[in] n_spsym Number of samples per symbol
/// @param[in] b Shape parameter (values defined for FT8/FT4)
/// @param[out] pulse Output array of pulse samples
///
void gfsk_pulse(int n_spsym, float symbol_bt, float* pulse)
{
    for (int i = 0; i < 3 * n_spsym; ++i)
    {
        float t = i / (float)n_spsym - 1.5f;
        float arg1 = GFSK_CONST_K * symbol_bt * (t + 0.5f);
        float arg2 = GFSK_CONST_K * symbol_bt * (t - 0.5f);
        pulse[i] = (erff(arg1) - erff(arg2)) / 2;
    }
}

/// Synthesize waveform data using GFSK phase shaping.
/// The output waveform will contain n_sym symbols.
/// @param[in] symbols Array of symbols (tones) (0-7 for FT8)
/// @param[in] n_sym Number of symbols in the symbol array
/// @param[in] f0 Audio frequency in Hertz for the symbol 0 (base frequency)
/// @param[in] symbol_bt Symbol smoothing filter bandwidth (2 for FT8, 1 for FT4)
/// @param[in] symbol_period Symbol period (duration), seconds
/// @param[in] signal_rate Sample rate of synthesized signal, Hertz
/// @param[out] signal Output array of signal waveform samples (should have space for n_sym*n_spsym samples)
///
void synth_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal)
{
    int n_spsym = (int)(0.5f + signal_rate * symbol_period); // Samples per symbol
    int n_wave = n_sym * n_spsym;                            // Number of output samples
    float hmod = 1.0f;

    printf("n_spsym = %d\n", n_spsym);
    // Compute the smoothed frequency waveform.
    // Length = (nsym+2)*n_spsym samples, first and last symbols extended
    float dphi_peak = 2 * M_PI * hmod / n_spsym;
    float dphi[n_wave + 2 * n_spsym];

    // Shift frequency up by f0
    for (int i = 0; i < n_wave + 2 * n_spsym; ++i)
    {
        dphi[i] = 2 * M_PI * f0 / signal_rate;
    }

    float pulse[3 * n_spsym];
    gfsk_pulse(n_spsym, symbol_bt, pulse);

    for (int i = 0; i < n_sym; ++i)
    {
        int ib = i * n_spsym;
        for (int j = 0; j < 3 * n_spsym; ++j)
        {
            dphi[j + ib] += dphi_peak * symbols[i] * pulse[j];
        }
    }

    // Add dummy symbols at beginning and end with tone values equal to 1st and last symbol, respectively
    for (int j = 0; j < 2 * n_spsym; ++j)
    {
        dphi[j] += dphi_peak * pulse[j + n_spsym] * symbols[0];
        dphi[j + n_sym * n_spsym] += dphi_peak * pulse[j] * symbols[n_sym - 1];
    }

    // Calculate and insert the audio waveform
    float phi = 0;
    for (int k = 0; k < n_wave; ++k)
    { // Don't include dummy symbols
        signal[k] = sinf(phi);
        phi = fmodf(phi + dphi[k + n_spsym], 2 * M_PI);
    }

    // Apply envelope shaping to the first and last symbols
    int n_ramp = n_spsym / 8;
    for (int i = 0; i < n_ramp; ++i)
    {
        float env = (1 - cosf(2 * M_PI * i / (2 * n_ramp))) / 2;
        signal[i] *= env;
        signal[n_wave - 1 - i] *= env;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FT8_lib decode

static float hann_i(int i, int N)
{
	float x = sinf((float)M_PI * i / N);
	return x * x;
}

static float hamming_i(int i, int N)
{
	const float a0 = (float)25 / 46;
	const float a1 = 1 - a0;

	float x1 = cosf(2 * (float)M_PI * i / N);
	return a0 - a1 * x1;
}

static float blackman_i(int i, int N)
{
	const float alpha = 0.16f; // or 2860/18608
	const float a0 = (1 - alpha) / 2;
	const float a1 = 1.0f / 2;
	const float a2 = alpha / 2;

	float x1 = cosf(2 * (float)M_PI * i / N);
	float x2 = 2 * x1 * x1 - 1; // Use double angle formula

	return a0 - a1 * x1 + a2 * x2;
}

void waterfall_init(waterfall_t* me, int max_blocks, int num_bins, int time_osr, int freq_osr)
{
	size_t mag_size = max_blocks * time_osr * freq_osr * num_bins * sizeof(me->mag[0]);
	me->max_blocks = max_blocks;
	me->num_blocks = 0;
	me->num_bins = num_bins;
	me->time_osr = time_osr;
	me->freq_osr = freq_osr;
	me->block_stride = (time_osr * freq_osr * num_bins);
	me->mag = (uint8_t  *)malloc(mag_size);
	#if DEBUG
	printf("Waterfall size = %zu\n", mag_size);
	#endif
}

void waterfall_free(waterfall_t* me)
{
	free(me->mag);
}

/// Configuration options for FT4/FT8 monitor
typedef struct
{
	float f_min;			 ///< Lower frequency bound for analysis
	float f_max;			 ///< Upper frequency bound for analysis
	int sample_rate;		 ///< Sample rate in Hertz
	int time_osr;			///< Number of time subdivisions
	int freq_osr;			///< Number of frequency subdivisions
	ftx_protocol_t protocol; ///< Protocol: FT4 or FT8
} monitor_config_t;

/// FT4/FT8 monitor object that manages DSP processing of incoming audio data
/// and prepares a waterfall object
typedef struct
{
	float symbol_period; ///< FT4/FT8 symbol period in seconds
	int block_size;	  ///< Number of samples per symbol (block)
	int subblock_size;   ///< Analysis shift size (number of samples)
	int nfft;			///< FFT size
	float fft_norm;	  ///< FFT normalization factor
	float* window;	   ///< Window function for STFT analysis (nfft samples)
	float* last_frame;   ///< Current STFT analysis frame (nfft samples)
	waterfall_t wf;	  ///< Waterfall object
	float max_mag;	   ///< Maximum detected magnitude (debug stats)

	// KISS FFT housekeeping variables
	void* fft_work;		///< Work area required by Kiss FFT
	kiss_fftr_cfg fft_cfg; ///< Kiss FFT housekeeping object
} monitor_t;

void monitor_init(monitor_t* me, const monitor_config_t* cfg)
{
	float slot_time = (cfg->protocol == PROTO_FT4) ? FT4_SLOT_TIME : FT8_SLOT_TIME;
	float symbol_period = (cfg->protocol == PROTO_FT4) ? FT4_SYMBOL_PERIOD : FT8_SYMBOL_PERIOD;
	// Compute DSP parameters that depend on the sample rate
	me->block_size = (int)(cfg->sample_rate * symbol_period); // samples corresponding to one FSK symbol
	me->subblock_size = me->block_size / cfg->time_osr;
	me->nfft = me->block_size * cfg->freq_osr;
	me->fft_norm = 2.0f / me->nfft;
	// const int len_window = 1.8f * me->block_size; // hand-picked and optimized

	me->window = (float *)malloc(me->nfft * sizeof(me->window[0]));
	for (int i = 0; i < me->nfft; ++i)
	{
		// window[i] = 1;
		me->window[i] = hann_i(i, me->nfft);
		// me->window[i] = blackman_i(i, me->nfft);
		// me->window[i] = hamming_i(i, me->nfft);
		// me->window[i] = (i < len_window) ? hann_i(i, len_window) : 0;
	}
	me->last_frame = (float *)malloc(me->nfft * sizeof(me->last_frame[0]));

	size_t fft_work_size;
	kiss_fftr_alloc(me->nfft, 0, 0, &fft_work_size);

	#if DEBUG
	printf( "Block size = %d\n", me->block_size);
	printf( "Subblock size = %d\n", me->subblock_size);
	printf( "N_FFT = %d\n", me->nfft);
	printf("FFT work area = %zu\n", fft_work_size);
	#endif

	me->fft_work = malloc(fft_work_size);
	me->fft_cfg = kiss_fftr_alloc(me->nfft, 0, me->fft_work, &fft_work_size);

	const int max_blocks = (int)(slot_time / symbol_period);
	const int num_bins = (int)(cfg->sample_rate * symbol_period / 2);
	waterfall_init(&me->wf, max_blocks, num_bins, cfg->time_osr, cfg->freq_osr);
	me->wf.protocol = cfg->protocol;
	me->symbol_period = symbol_period;

	me->max_mag = -120.0f;
}

void monitor_free(monitor_t* me)
{
	waterfall_free(&me->wf);
	free(me->fft_work);
	free(me->last_frame);
	free(me->window);
}

// Compute FFT magnitudes (log wf) for a frame in the signal and update waterfall data
void monitor_process(monitor_t* me, const float* frame)
{
	// Check if we can still store more waterfall data
	if (me->wf.num_blocks >= me->wf.max_blocks)
		return;

	int offset = me->wf.num_blocks * me->wf.block_stride;
	int frame_pos = 0;

	// Loop over block subdivisions
	for (int time_sub = 0; time_sub < me->wf.time_osr; ++time_sub)
	{
		kiss_fft_scalar timedata[me->nfft];
		kiss_fft_cpx freqdata[me->nfft / 2 + 1];

		// Shift the new data into analysis frame
		for (int pos = 0; pos < me->nfft - me->subblock_size; ++pos)
		{
			me->last_frame[pos] = me->last_frame[pos + me->subblock_size];
		}
		for (int pos = me->nfft - me->subblock_size; pos < me->nfft; ++pos)
		{
			me->last_frame[pos] = frame[frame_pos];
			++frame_pos;
		}

		// Compute windowed analysis frame
		for (int pos = 0; pos < me->nfft; ++pos)
		{
			timedata[pos] = me->fft_norm * me->window[pos] * me->last_frame[pos];
		}

		kiss_fftr(me->fft_cfg, timedata, freqdata);

		// Loop over two possible frequency bin offsets (for averaging)
		for (int freq_sub = 0; freq_sub < me->wf.freq_osr; ++freq_sub)
		{
			for (int bin = 0; bin < me->wf.num_bins; ++bin)
			{
				int src_bin = (bin * me->wf.freq_osr) + freq_sub;
				float mag2 = (freqdata[src_bin].i * freqdata[src_bin].i) + (freqdata[src_bin].r * freqdata[src_bin].r);
				float db = 10.0f * log10f(1E-12f + mag2);
				// Scale decibels to unsigned 8-bit range and clamp the value
				// Range 0-240 covers -120..0 dB in 0.5 dB steps
				int scaled = (int)(2 * db + 240);

				me->wf.mag[offset] = (scaled < 0) ? 0 : ((scaled > 255) ? 255 : scaled);
				++offset;

				if (db > me->max_mag)
					me->max_mag = db;
			}
		}
	}

	++me->wf.num_blocks;
}

void monitor_reset(monitor_t* me)
{
	me->wf.num_blocks = 0;
	me->max_mag = 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Thread FT8
float getFrame(char *buffer, int i)
{
	return ((float)(buffer[(2 * i)] & 0xFF) + ((buffer[(2 * i) + 1] & 0xFF) << 8)) / 32768.0f;
}

void RX_FT8()
{
	
	const int num_samples = 13.6 * sound.capture_sound_rate;
	const int kMax_candidates = 120;
	const int kMin_score = 2; // (10) Minimum sync score threshold for candidates
	const int kMax_decoded_messages = 50;
	const int kLDPC_iterations = 20;
	
	monitor_config_t mon_cfg = {
		.f_min = 100,
		.f_max = 3200,
		.sample_rate = (int)sound.capture_sound_rate,
		.time_osr = 2,
		.freq_osr = 2,
		.protocol = PROTO_FT8
	};
	
	snd_pcm_sframes_t rc;
	
	while(1){
		pthread_mutex_lock(&FT8.TRX_status_lock);
		bool stat = FT8.TRX_status;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if(stat == _RX_){
			// Compute FFT over the whole signal and store it
			monitor_t mon;
			monitor_init(&mon, &mon_cfg);
			
			#if DEBUG
			printf( "Waterfall allocated %d symbols\n", mon.wf.max_blocks);
			#endif

			float signal[mon.block_size];
			char *raw_data = (char *)malloc(mon.block_size*2);
		
			waitforstart();
			snd_pcm_reset(sound.capture_handle);

			for (int frame_pos = 0; frame_pos + mon.block_size <= num_samples; frame_pos += mon.block_size)
			{
				// Process the waveform data frame by frame - you could have a live loop here with data from an audio device

				rc = snd_pcm_readi (sound.capture_handle, raw_data, mon.block_size);
				
				if (rc == -EPIPE)
				{
					/* EPIPE means overrun */
					snd_pcm_recover(sound.capture_handle, rc, 0);
				}
				else if (rc == -EAGAIN)
				{
					/* Not ready yet. Come back again later. */
				}
				else if (rc < 0)
				{
					fprintf(stderr, "error from read: %s\n", snd_strerror(rc));
				}
				else
				{
				
					for (int i = 0; i < mon.block_size; i++)
					{
						signal[i] = getFrame(raw_data,i);
					}
					
					monitor_process(&mon, signal);
				}
				
			}
			
			#if DEBUG
			printf( "Waterfall accumulated %d symbols\n", mon.wf.num_blocks);
			printf( "Max magnitude: %.1f dB\n", mon.max_mag);
			#endif
			
			// Find top candidates by Costas sync score and localize them in time and frequency
			candidate_t candidate_list[kMax_candidates];
			int num_candidates = ft8_find_sync(&mon.wf, kMax_candidates, candidate_list, kMin_score);

			//Creat array for analyse
			char AnalyseArray[num_candidates][3][25];
			int countanalyse=0;

			// Hash table for decoded messages (to check for duplicates)
			int num_decoded = 0;
			message_t decoded[kMax_decoded_messages];
			message_t* decoded_hashtable[kMax_decoded_messages];

			// Initialize hash table pointers
			for (int i = 0; i < kMax_decoded_messages; ++i)
			{
				decoded_hashtable[i] = NULL;
			}

			// Go over candidates and attempt to decode messages
			for (int idx = 0; idx < num_candidates; ++idx)
			{
				
				const candidate_t* cand = &candidate_list[idx];
				if (cand->score < kMin_score)
					continue;

				float freq_hz = (cand->freq_offset + (float)cand->freq_sub / mon.wf.freq_osr) / mon.symbol_period;
				float time_sec = (cand->time_offset + (float)cand->time_sub / mon.wf.time_osr) * mon.symbol_period;

				message_t message;
				decode_status_t status;
				if (!ft8_decode(&mon.wf, cand, &message, kLDPC_iterations, &status))
				{
					// printf("000000 %3d %+4.2f %4.0f ~  ---\n", cand->score, time_sec, freq_hz);
					#if DEBUG
					if (status.ldpc_errors > 0)
					{
						printf( "LDPC decode: %d errors\n", status.ldpc_errors);
					}
					else if (status.crc_calculated != status.crc_extracted)
					{
						printf( "CRC mismatch!\n");
					}
					else if (status.unpack_status != 0)
					{
						printf( "Error while unpacking!\n");
					}
					#endif
					continue;
				}

				#if DEBUG
				printf( "Checking hash table for %4.1fs / %4.1fHz [%d]...\n", time_sec, freq_hz, cand->score);
				#endif
				int idx_hash = message.hash % kMax_decoded_messages;
				bool found_empty_slot = false;
				bool found_duplicate = false;
				do
				{
					if (decoded_hashtable[idx_hash] == NULL)
					{
						#if DEBUG
						printf( "Found an empty slot\n");
						#endif
						found_empty_slot = true;
					}
					else if ((decoded_hashtable[idx_hash]->hash == message.hash) && (0 == strcmp(decoded_hashtable[idx_hash]->text, message.text)))
					{
						#if DEBUG
						printf( "Found a duplicate [%s]\n", message.text);
						#endif
						found_duplicate = true;
					}
					else
					{
						#if DEBUG
						printf( "Hash table clash!\n");
						#endif
						// Move on to check the next entry in hash table
						idx_hash = (idx_hash + 1) % kMax_decoded_messages;
					}
				} while (!found_empty_slot && !found_duplicate);

				if (found_empty_slot)
				{
					// Fill the empty hashtable slot
					memcpy(&decoded[idx_hash], &message, sizeof(message));
					decoded_hashtable[idx_hash] = &decoded[idx_hash];
					++num_decoded;

					
					int snr = 0; // TODO: compute SNR
					
					printDateTime();
					
					if (strncmp(message.text, FT8.Local_CALLSIGN, strlen(FT8.Local_CALLSIGN)) == 0) {
						printf(" %3d %+4.2f %4.0f ~  \e[1;31m%s\e[0m\n", cand->score, time_sec, freq_hz, message.text);
						unpackFT8mess(message.text,AnalyseArray[0][0],AnalyseArray[0][1],AnalyseArray[0][2]);
						countanalyse=-1;
					}
					if ((strncmp(message.text,"CQ",2) == 0) && (countanalyse>-1)) {
						printf(" %3d %+4.2f %4.0f ~  \e[1;34m%s\e[0m\n", cand->score, time_sec, freq_hz, message.text);
						unpackFT8mess(message.text,AnalyseArray[countanalyse][0],AnalyseArray[countanalyse][1],AnalyseArray[countanalyse][2]);
						countanalyse++;
					}else{
						printf(" %3d %+4.2f %4.0f ~  %s\n", cand->score, time_sec, freq_hz, message.text);
					}
				}
			}
			
			int index_from_ope = 0;
			float dist = 0;
			if(countanalyse>0){
				for(int i = 0; i<countanalyse;i++){
					float latlonlocal[2];
					latLonForGrid(AnalyseArray[i][2],latlonlocal);
					#if DEBUG
					printf("%s lat:%f lon:%f\n",AnalyseArray[i][1],latlonlocal[0],latlonlocal[1]);
					#endif
					float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
					if(dist < new_dist){index_from_ope=i;dist=new_dist;}
				}
				printf("*Selected for new QSO: \e[1;31m%s at %f km\e[0m\n", AnalyseArray[index_from_ope][1],dist);
				strcpy(FT8.QSO_dist_CALLSIGN, AnalyseArray[index_from_ope][1]);
				strcpy(FT8.QSO_dist_LOCATOR, AnalyseArray[index_from_ope][2]);
			}
			
			
			if(!num_decoded){printDateTime();printf(" N/A\n");}
			
			#if DEBUG
			printf( "Decoded %d messages\n", num_decoded);
			#endif
			
			monitor_free(&mon);

			if(*FT8.QSO_dist_CALLSIGN != 0) {
				printf( "UnLock TX thread\n");
				pthread_mutex_lock(&FT8.TRX_status_lock);
				FT8.TRX_status = _TX_;
				pthread_cond_signal(&FT8.TX_status_cond);
				pthread_mutex_unlock(&FT8.TRX_status_lock);				
			}
		
		}
		else{
			printf( "Lock RX thread\n");
			pthread_mutex_lock(&FT8.TRX_status_lock);
			pthread_cond_wait(&FT8.RX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
		
	}

}



char* join_strings(char* strings[], char* seperator, int count) {
    char* str = NULL;             /* Pointer to the joined strings  */
    size_t total_length = 0;      /* Total length of joined strings */
    int i = 0;                    /* Loop counter                   */

    /* Find total length of joined strings */
    for (i = 0; i < count; i++) total_length += strlen(strings[i]);
    total_length++;     /* For joined string terminator */
    total_length += strlen(seperator) * (count - 1); // for seperators

    str = (char*) malloc(total_length);  /* Allocate memory for joined strings */
    str[0] = '\0';                      /* Empty string we can append to      */

    /* Append all the strings */
    for (i = 0; i < count; i++) {
        strcat(str, strings[i]);
        if (i < (count - 1)) strcat(str, seperator);
    }

    return str;
}

void gen_FT8_responses()
{
	// strcat(FT8.QSO[0],FT8.Local_CALLSIGN); strcat(FT8.QSO[0],distop); strcat(FT8.QSO[0],FT8.Local_LOCATOR);

	char snr[3],snrr[4]="R";
	sprintf(snr, "%d", FT8.QSO_dist_SNR); 
	strcat(snrr,snr);

	char* strings[3];
	
	strings[0]=FT8.QSO_dist_CALLSIGN;
	strings[1]=FT8.Local_CALLSIGN;
	
	strings[2]=FT8.Local_LOCATOR;
	strcpy(FT8.QSO_RESPONSES[0], join_strings(strings,(char*)" ",3));

	
	strings[2]=snr;
	strcpy(FT8.QSO_RESPONSES[1], join_strings(strings,(char*)" ",3));

	
	strings[2]=snrr;
	strcpy(FT8.QSO_RESPONSES[2], join_strings(strings,(char*)" ",3));

	
	strings[2]=(char*)"RRR";
	strcpy(FT8.QSO_RESPONSES[3], join_strings(strings,(char*)" ",3));

	
	strings[2]=(char*)"73";
	strcpy(FT8.QSO_RESPONSES[4], join_strings(strings,(char*)" ",3));

	
	printf("%s\n",FT8.QSO_RESPONSES[0]);
	printf("%s\n",FT8.QSO_RESPONSES[1]);
	printf("%s\n",FT8.QSO_RESPONSES[2]);
	printf("%s\n",FT8.QSO_RESPONSES[3]);
	printf("%s\n",FT8.QSO_RESPONSES[4]);
	printf("%s\n",FT8.QSO_RESPONSES[5]);

	
}

void Reinit_FT8_QSO()
{
	FT8.QSO_Index_to_rep=-1;
	for(int i = 0; i<5;i++){
		strcpy(FT8.QSO_RESPONSES[i],"");
	}
	strcpy(FT8.QSO_dist_CALLSIGN,"");
	strcpy(FT8.QSO_dist_LOCATOR,"");
	strcpy(FT8.QSO_dist_MESSAGE,"");
	FT8.QSO_dist_SNR=0;
}

void TX_FT8()
{
	while(1){
		pthread_mutex_lock(&FT8.TRX_status_lock);
		bool stat = FT8.TRX_status;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if(stat == _TX_){
			
			if(FT8.QSO_RESPONSES[0][0]=='\0'){gen_FT8_responses();}
			
			// uint8_t packed[FTX_LDPC_K_BYTES];
			// int rc = pack77(message, packed);
			// if (rc < 0)
			// {
				// printf("Cannot parse message!\n");
				// printf("RC = %d\n", rc);
				// return -2;
			// }

			// printf("Packed data: ");
			// for (int j = 0; j < 10; ++j)
			// {
				// printf("%02x ", packed[j]);
			// }
			// printf("\n");

			// int num_tones = FT8_NN;
			// float symbol_period = FT8_SYMBOL_PERIOD;
			// float symbol_bt = FT8_SYMBOL_BT;
			// float slot_time = FT8_SLOT_TIME;

			// // Second, encode the binary message as a sequence of FSK tones
			// uint8_t tones[num_tones]; // Array of 79 tones (symbols)
			// ft8_encode(packed, tones);
			

			// printf("FSK tones: ");
			// for (int j = 0; j < num_tones; ++j)
			// {
				// printf("%d", tones[j]);
			// }
			// printf("\n");

			// // Third, convert the FSK tones into an audio signal
			// int sample_rate = 12000;
			// int num_samples = (int)(0.5f + num_tones * symbol_period * sample_rate); // Number of samples in the data signal
			// int num_silence = (slot_time * sample_rate - num_samples) / 2;           // Silence padding at both ends to make 15 seconds
			// int num_total_samples = num_silence + num_samples + num_silence;         // Number of samples in the padded signal
			// float signal[num_total_samples];
			// for (int i = 0; i < num_silence; i++)
			// {
				// signal[i] = 0;
				// signal[i + num_samples + num_silence] = 0;
			// }

			// // Synthesize waveform data (signal) and save it as WAV file
			// synth_gfsk(tones, num_tones, frequency, symbol_bt, symbol_period, sample_rate, signal + num_silence);
			// save_wav(signal, num_total_samples, sample_rate, wav_path);

			// return 0;	
			
			
			
			//Empty QSO variable
			pthread_mutex_lock(&FT8.TRX_status_lock);
			Reinit_FT8_QSO();
			pthread_mutex_unlock(&FT8.TRX_status_lock);
			
			//unlock thread RX
			printf( "UnLock RX thread\n");
			pthread_mutex_lock(&FT8.TRX_status_lock);
			FT8.TRX_status = _RX_;
			pthread_cond_signal(&FT8.RX_status_cond);
			pthread_mutex_unlock(&FT8.TRX_status_lock);
		}
		else{
			printf( "Lock TX thread\n");
			pthread_mutex_lock(&FT8.TRX_status_lock);
			pthread_cond_wait(&FT8.TX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
	}
}

void * Thread_RX(void *arg) {
	RX_FT8();
	return NULL;
}

void * Thread_TX(void *arg) {
	TX_FT8();
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//main
int main (int argc, char *argv[])
{
	int c;
	while ((c = getopt (argc, argv, "hd:C:L:")) != -1)
		switch (c)
		{
			case 'h':
				printf ("clFT8 -d plughw:CARD=PCH,DEV=0 -c F4JJJ\n");
				exit(0);
				break;
			case 'd':
				sound.capture_sound_device = optarg;
				sound.playback_sound_device = optarg;
				break;
				return 1;
			case 'C':
				strcpy( FT8.Local_CALLSIGN, optarg );
				break;
				return 1;
			case 'L':
				strcpy( FT8.Local_LOCATOR, optarg );
				break;
				return 1;
			default:
				abort ();
		}
	
	latLonForGrid(FT8.Local_LOCATOR,FT8.Local_latlon);
	
	char* strings[3];strings[0]=(char*)"CQ";strings[1]=FT8.Local_CALLSIGN;strings[2]=FT8.Local_LOCATOR;
	strcpy(FT8.QSO_RESPONSES[5], join_strings(strings,(char*)" ",3));	
	
	capture_audioInit();
	playback_audioInit();
	
	pthread_t thread_RX;
	pthread_create(&thread_RX, NULL, Thread_RX, NULL);
	
	pthread_t thread_TX;
	pthread_create(&thread_TX, NULL, Thread_TX, NULL);
	
	while (1)
	{
		sleep(1);
		#if DEBUG
		printf("wait\n");
		#endif
	}
	
	
}
