#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include<ctype.h>
#include <stdbool.h>

#include <alsa/asoundlib.h>
#include <pthread.h>
#include <sys/time.h>

#include "ft8/unpack.h"
#include "ft8/pack.h"
#include "ft8/ldpc.h"
#include "ft8/decode.h"
#include "ft8/constants.h"
#include "ft8/encode.h"
#include "ft8/crc.h"

#include "fft/kiss_fftr.h"

#include "serial/cssl.h"

#include "clFT8.h"

#include "hash/hash.h"




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Global

/* Global FT8 info. */


FT8info FT8 = {
	.Local_CALLSIGN = {'F','4','J','J','J','\0'},
	.Local_LOCATOR = {'J','N','3','8','\0'},
	.TX_enable = 1,
	
	.QSO_dist_CALLSIGN = "",
	.QSO_dist_LOCATOR = "",
	.QSO_dist_MESSAGE = "",
	.QSO_dist_SNR = 0,
	.QSO_dist_FREQUENCY = 0,
	
	.QSO_RESPONSES={'\0','\0','\0','\0','\0','\0'},
	.QSO_Index_to_rep=-1,
	
	.TRX_status = _RX_,
	.TRX_status_lock = PTHREAD_MUTEX_INITIALIZER,
	.RX_status_cond = PTHREAD_COND_INITIALIZER,
	.TX_status_cond = PTHREAD_COND_INITIALIZER,
	
	.Tranceiver_VFOA_Freq = 14074000,
	
	.log_file_name = "QSO.log",
	.infos_to_log[0] = '\0',
	.log_callsigntable_file_name = "Call_Table.log",
	.log_dist_CALLSIGN_for_filter = "",
	.callsigntable_for_filter_index=0,
	.beep_on_log=0,
	
	.filter_on_cq = 0
	
};

HashTable* ht_callsigntable_for_filter;

/* Global sound info. */
soundInfo sound={
	.capture_sound_device = (char*)"default",
	.capture_sound_rate = 12000,
	
	.playback_sound_device = (char*)"default",
	.playback_buffer_frames = 1024,
	.playback_sound_rate = 12000
};

/* Global serial info. */

serial_t serial = {
	.pathname="/dev/ttyACM0",
	.rtscts=0,
	.xonxoff=0,
	.baud=9600,
	.bits=8,
	.parity=0,
	.stopbits=1,
	.finished=0
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

void advance_cursor() {
  static int pos=0;
  char cursor[4]={'/','-','\\','|'};
  printf("%c\b", cursor[pos]);
  fflush(stdout);
  pos = (pos+1) % 4;
}

void waitforstart(){
	double tt = now();
	long long cycle_start = tt - ((long long)tt % 15);
	double WaitTime = (15 - (double)(tt - cycle_start))*1000000;
	usleep(WaitTime);
}

void printDateTime_log(){
	time_t t = time(NULL);	struct tm tm = *gmtime(&t);	int sec = (int)((int)tm.tm_sec / 15)*15;
	printf("%d-%02d-%02d %02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, sec);
}

void printDateTime(){
	time_t t = time(NULL);	struct tm tm = *gmtime(&t);	int sec = (int)((int)tm.tm_sec);
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
	return acos(sin(latlonA[0] * M_PI / 180.0f) * sin(latlonB[0] * M_PI / 180.0f) + cos(latlonA[0] * M_PI / 180.0f) * cos(latlonB[0] * M_PI / 180.0f) * cos((latlonA[1] * M_PI / 180.0f) - (latlonB[1] * M_PI / 180.0f)))*6371.0;
}

int count_occur_str(char * s, char * c){
	int count=0;
	for(int i=0;s[i];i++)  
    {
    	if(s[i]==*c)
    	{
          count++;
		}
 	}
	return count;
}

void unpackFT8mess(char * message_text, char * unpackeds0, char * unpackeds1, char * unpackeds2){
	
	char * NullToken = strtok(message_text, " ");
	if(NullToken != 0){
		memcpy(unpackeds0, NullToken, sizeof(NullToken));
	}else{strcpy(unpackeds0, "");}
	
	NullToken = strtok(0, " ");	
	if(NullToken != 0){
		memcpy(unpackeds1, NullToken, sizeof(NullToken));
	}else{strcpy(unpackeds1, "");}
	
	NullToken = strtok(0, " ");		
	if(NullToken != 0){
		memcpy(unpackeds2, NullToken, sizeof(NullToken));
	}else{strcpy(unpackeds2, "");}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Audio



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

    //printf("n_spsym = %d\n", n_spsym);
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

void unlock_TX_thread(){
	#if DEBUG
	printf( "UnLock TX thread\n");
	#endif
	pthread_mutex_lock(&FT8.TRX_status_lock);
	FT8.TRX_status = _TX_;
	pthread_cond_signal(&FT8.TX_status_cond);
	pthread_mutex_unlock(&FT8.TRX_status_lock);
}

void unlock_RX_thread(){
	#if DEBUG
	printf( "UnLock RX thread\n");
	#endif
	pthread_mutex_lock(&FT8.TRX_status_lock);
	FT8.TRX_status = _RX_;
	pthread_cond_signal(&FT8.RX_status_cond);
	pthread_mutex_unlock(&FT8.TRX_status_lock);
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
			float  AnalyseArrayFreqInfo[num_candidates];
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
					
					printDateTime_log();
					
					if ((strncmp(message.text, FT8.Local_CALLSIGN, strlen(FT8.Local_CALLSIGN)) == 0) && (countanalyse>-1)) {
						
						countanalyse=-1;
						printf(" %3d %+4.2f %4.0f ~  \e[1;31m%s\e[0m\n", cand->score, time_sec, freq_hz, message.text);
						unpackFT8mess(message.text,AnalyseArray[0][0],AnalyseArray[0][1],AnalyseArray[0][2]);
						
						pthread_mutex_lock(&FT8.TRX_status_lock);
								
						strcpy(FT8.QSO_dist_CALLSIGN, AnalyseArray[0][1]);
						strcpy(FT8.QSO_dist_MESSAGE, AnalyseArray[0][2]);
						FT8.QSO_dist_FREQUENCY=freq_hz;
						FT8.QSO_dist_SNR = (cand->score -30)/2;
						FT8.QSO_Index_to_rep = -1;
				
						pthread_mutex_unlock(&FT8.TRX_status_lock);
				
						unlock_TX_thread();
						
					}
					else if ((strncmp(message.text,"CQ",2) == 0) && (countanalyse>-1)) {
						char trashmess[25];
						strcpy(trashmess, message.text);
						unpackFT8mess(trashmess,AnalyseArray[countanalyse][0],AnalyseArray[countanalyse][1],AnalyseArray[countanalyse][2]);
						
						AnalyseArrayFreqInfo[countanalyse]=freq_hz;
						
						if(ht_check(ht_callsigntable_for_filter,AnalyseArray[countanalyse][1]) || strlen(AnalyseArray[countanalyse][2])==0 || count_occur_str(message.text, " ") > 2){
							unpackFT8mess("",AnalyseArray[countanalyse][0],AnalyseArray[countanalyse][1],AnalyseArray[countanalyse][2]);
							printf(" %3d %+4.2f %4.0f ~  \e[1;35m%s\e[0m\n", cand->score, time_sec, freq_hz, message.text);
						}
						else{
							countanalyse++;
							printf(" %3d %+4.2f %4.0f ~  \e[1;34m%s\e[0m\n", cand->score, time_sec, freq_hz, message.text);
							}
							
					
					}else{
						printf(" %3d %+4.2f %4.0f ~  %s\n", cand->score, time_sec, freq_hz, message.text);
					}
				}
			}
			
			
			if(countanalyse>0){
				int index_from_ope = 0;
				float dist = 0;
				
				switch (FT8.filter_on_cq)
				{
				
				case '0':
					index_from_ope = rand() % countanalyse;
					break;
					
				case '1':
					index_from_ope = 0;
					break;				
				
				case '2':
					for(int i = 0; i<countanalyse;i++){
						if(strlen(AnalyseArray[i][2]) == 4){
							float latlonlocal[2];
							latLonForGrid(AnalyseArray[i][2],latlonlocal);
							#if DEBUG
							printf("%s lat:%f lon:%f\n",AnalyseArray[i][1],latlonlocal[0],latlonlocal[1]);
							#endif
							float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
							if(dist < new_dist){index_from_ope=i;dist=new_dist;}
						}
					}				
					break;
					
				case '3':
					for(int i = 0; i<countanalyse;i++){
						if(strlen(AnalyseArray[i][2]) == 4){
							float latlonlocal[2];
							latLonForGrid(AnalyseArray[i][2],latlonlocal);
							#if DEBUG
							printf("%s lat:%f lon:%f\n",AnalyseArray[i][1],latlonlocal[0],latlonlocal[1]);
							#endif
							float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
							if(dist > new_dist){index_from_ope=i;dist=new_dist;}
						}
					}				
					break;
				
				default:
					index_from_ope = rand() % countanalyse;
					break;

				}
				
				pthread_mutex_lock(&FT8.TRX_status_lock);
				
				strcpy(FT8.QSO_dist_CALLSIGN, AnalyseArray[index_from_ope][1]);
				strcpy(FT8.QSO_dist_LOCATOR, AnalyseArray[index_from_ope][2]);
				FT8.QSO_Index_to_rep=0;
				FT8.QSO_dist_FREQUENCY=AnalyseArrayFreqInfo[index_from_ope];
				
				pthread_mutex_unlock(&FT8.TRX_status_lock);
				
				printf("*Selected for new QSO: \e[1;31m%s at %s on %f hz (seq QSO on %d) \e[0m\n", FT8.QSO_dist_CALLSIGN,FT8.QSO_dist_LOCATOR,FT8.QSO_dist_FREQUENCY,FT8.QSO_Index_to_rep);
				
				unlock_TX_thread();
				
			}
			
			
			if(!num_decoded){printDateTime_log();printf(" N/A\n");}
			
			#if DEBUG
			printf( "Decoded %d messages\n", num_decoded);
			#endif
			
			monitor_free(&mon);
		
		}
		else{
			// printf( "Lock RX thread\n");
			pthread_mutex_lock(&FT8.TRX_status_lock);
			pthread_cond_wait(&FT8.RX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
		
	}

}

void gen_FT8_responses()
{
	sprintf(FT8.QSO_RESPONSES[0], "%s %s %s", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN, FT8.Local_LOCATOR); 
	sprintf(FT8.QSO_RESPONSES[1], "%s %s %+1.2d", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN, FT8.QSO_dist_SNR);
	sprintf(FT8.QSO_RESPONSES[2], "%s %s R%+1.2d", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN, FT8.QSO_dist_SNR);
	sprintf(FT8.QSO_RESPONSES[3], "%s %s RRR", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN);
	sprintf(FT8.QSO_RESPONSES[4], "%s %s 73", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN);
	
	// printf("%s\n",FT8.QSO_RESPONSES[0]);
	// printf("%s\n",FT8.QSO_RESPONSES[1]);
	// printf("%s\n",FT8.QSO_RESPONSES[2]);
	// printf("%s\n",FT8.QSO_RESPONSES[3]);
	// printf("%s\n",FT8.QSO_RESPONSES[4]);
	// printf("%s\n",FT8.QSO_RESPONSES[5]);
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

int get_seq_qso_to_rep(char * mess, bool * flaglog)
{
	int rep = -1;*flaglog=false;
	
	if ((mess[0] == '-') && (isdigit(mess[strlen(mess)-1]))) 
	{
		rep = 2;*flaglog=false;
	}
	
	if ((mess[0] == '+') && (isdigit(mess[strlen(mess)-1])))  
	{
		rep = 2;*flaglog=false;
	}

	if ((mess[0] == 'R') && (isdigit(mess[strlen(mess)-1]))) 
	{
		rep = 3;*flaglog=false;
	}

	if ((mess[0] == 'R') && (mess[1] == 'R') && (mess[2] == 'R')) 
	{
		rep = 4;*flaglog=false;
	}
	
	if ((mess[0] == 'R') && (mess[1] == 'R') && (mess[2] == '7') && (mess[3] == '3')) 
	{
		rep = 4;*flaglog=true;
	}

	if ((mess[0] == '7') && (mess[1] == '3')) 
	{
		rep = -1;*flaglog=true;
	}
	
	return rep;
	
}

void log_FT8_QSO()
{
	
	char resultime[40];
	time_t t = time(NULL);	struct tm tm = *gmtime(&t);	int sec = (int)((int)tm.tm_sec);
	sprintf(resultime,"%d-%02d-%02d %02d:%02d:%02d	", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, sec);
	
	printf("\e[1;32mWe can log the qso %s	%s \e[0m\n", resultime, FT8.infos_to_log);
	if(FT8.beep_on_log){putchar('\07');putchar('\a');}
	FILE *fptr;
	fptr = fopen(FT8.log_file_name,"ab");
	if(fptr == NULL)
	{
		printf("Error with QSO log file!");   
		exit(1);             
	}

	fprintf(fptr,resultime);
	fprintf(fptr,FT8.infos_to_log);
	FT8.infos_to_log[0] = '\0';
	fclose(fptr);
	
}

void log_FT8_log_to_callsigntable_ht()
{
	FILE *fptr;
	fptr = fopen(FT8.log_callsigntable_file_name,"ab");
	if(fptr == NULL)
	{
		printf("Error with QSO log file!");   
		exit(1);             
	}

	fprintf(fptr,FT8.log_dist_CALLSIGN_for_filter);
	fprintf(fptr,"\n");
	fclose(fptr);
	ht_insert(ht_callsigntable_for_filter, FT8.log_dist_CALLSIGN_for_filter);
	FT8.log_dist_CALLSIGN_for_filter[0] = '\0';
}

void log_FT8_open_callsigntable_ht()
{
	ht_callsigntable_for_filter = ht_create_table();
	
	FILE *fptr;
	fptr = fopen(FT8.log_callsigntable_file_name,"a+");
	char tmp[20]="";
	while(fgets(tmp, 20, fptr)) 
	{
		tmp[strcspn(tmp, "\n")] = 0;
		tmp[strcspn(tmp, "\r")] = 0;
		ht_insert(ht_callsigntable_for_filter, tmp);
    }
	fclose(fptr);
	printf("Callsign hash table initialised with %d entry.\n",ht_callsigntable_for_filter->count);
	#if DEBUG
	print_table(ht_callsigntable_for_filter);
	#endif
}

void TX_FT8()
{
	while(1){
		pthread_mutex_lock(&FT8.TRX_status_lock);
		bool stat = FT8.TRX_status;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if(stat == _TX_){
			
			bool flaglog = 0;
			
			gen_FT8_responses();
			
			if(FT8.QSO_Index_to_rep == -1){FT8.QSO_Index_to_rep=get_seq_qso_to_rep(FT8.QSO_dist_MESSAGE,&flaglog);}
	
			if(FT8.QSO_Index_to_rep>-1){
				
				printf("Resp to \e[1;31m %s \e[0m with seq \e[1;31m %d \e[0m mess \e[1;31m %s \e[0m\n",FT8.QSO_dist_CALLSIGN, FT8.QSO_Index_to_rep, FT8.QSO_RESPONSES[FT8.QSO_Index_to_rep]);
			
				uint8_t packed[FTX_LDPC_K_BYTES];
				int rc = pack77(FT8.QSO_RESPONSES[FT8.QSO_Index_to_rep], packed);
				if (rc < 0)
				{
					printf("Cannot parse message!\n");
					printf("RC = %d\n", rc);
				}

				#if DEBUG
				printf("Packed data: ");
				for (int j = 0; j < 10; ++j)
				{
					printf("%02x ", packed[j]);
				}
				printf("\n");
				#endif

				int num_tones = FT8_NN;
				float symbol_period = FT8_SYMBOL_PERIOD;
				float symbol_bt = FT8_SYMBOL_BT;
				float slot_time = FT8_SLOT_TIME;

				// Second, encode the binary message as a sequence of FSK tones
				uint8_t tones[num_tones]; // Array of 79 tones (symbols)
				ft8_encode(packed, tones);
				
				#if DEBUG
				printf("FSK tones: ");
				for (int j = 0; j < num_tones; ++j)
				{
					printf("%d", tones[j]);
				}
				printf("\n");
				#endif

				// Third, convert the FSK tones into an audio signal
				int sample_rate = sound.playback_sound_rate;
				int num_samples = (int)(0.5f + num_tones * symbol_period * sample_rate); // Number of samples in the data signal
				int num_silence = (slot_time * sample_rate - num_samples) / 2;
				int num_total_samples = num_silence + num_samples + num_silence;         // Number of samples in the padded signal
				float signal[num_total_samples];
				
				for (int i = 0; i < num_silence; i++)
				{
					signal[i] = 0;
					signal[i + num_samples + num_silence] = 0;
				}
				
				synth_gfsk(tones, num_tones, FT8.QSO_dist_FREQUENCY, symbol_bt, symbol_period, sample_rate, signal);
				
				int16_t * raw_data = (int16_t*)malloc(num_total_samples*sizeof(int16_t)); // num_samples * numChannels * bitsPerSample / 8;

				for (int i = 0; i < num_total_samples; i++)
				{
					float x = signal[i];
					if (x > 1.0)
					x = 1.0;
					else if (x < -1.0)
					x = -1.0;
					int16_t xx = (int16_t)(0.5 + (x * 32767.0));
					raw_data[i] = (int)(0.5 + (x * 32767.0));
				}

				register snd_pcm_uframes_t	count, frames;
				
				waitforstart();

				printDateTime();printf(" start send message\n",num_total_samples);
				
				tranceiver_rtx(_TX_);
								
				int frame_lenght = sound.playback_buffer_frames;
				int frame_lenght_limit = num_total_samples - frame_lenght;
				count = 0;
				do
				{
					frames = snd_pcm_writei(sound.playback_handle, raw_data + count, frame_lenght);

					// If an error, try to recover from it
					if (frames < 0){frames = snd_pcm_recover(sound.playback_handle, frames, 0);}
					if (frames < 0)
					{
						printf("Error playing wave: %s\n", snd_strerror(frames));
						break;
					}
					
					// Update our pointer
					if (frames >= 0)
					{
						count += frames;
						if(count > frame_lenght_limit)
						{
							frame_lenght = num_total_samples - count;
						}
						else{}
					}
					
					

				} while (count < num_total_samples);
				
				if (count == num_total_samples){snd_pcm_drain(sound.playback_handle);}
				
				tranceiver_rtx(_RX_);

				printDateTime();printf(" stop send message\n");

				free(raw_data);
			
			}
			
			unlock_RX_thread();
			
			if(flaglog)
			{
				sprintf(FT8.infos_to_log, "%s	%s	%d	%f\n",FT8.QSO_dist_CALLSIGN,FT8.QSO_dist_LOCATOR,FT8.QSO_dist_SNR,FT8.QSO_dist_FREQUENCY + (float)FT8.Tranceiver_VFOA_Freq);
				sprintf(FT8.log_dist_CALLSIGN_for_filter,"%s",FT8.QSO_dist_CALLSIGN);
				//Empty QSO variable
				pthread_mutex_lock(&FT8.TRX_status_lock);
				Reinit_FT8_QSO();
				pthread_mutex_unlock(&FT8.TRX_status_lock);
			}
			
		}
		else{
			#if DEBUG
			printf( "Lock TX thread\n");
			#endif
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
//serial and tranceiver commands

int serial_init(){
    cssl_start();
    serial.port=cssl_open(serial.pathname,NULL,0,serial.baud,serial.bits,serial.parity,serial.stopbits);
	cssl_setflowcontrol(serial.port,serial.rtscts,serial.xonxoff);
	cssl_settimeout(serial.port,500);
    if (!serial.port) {
	printf("Serial error %s\n",cssl_geterrormsg());
	return -1;
    }else{return 0;}
}

void get_serial_rep(char rep[16])
{
	int i;
	for(i=0;i<15;i++){
		rep[i]=((char)cssl_getchar(serial.port));
		if(rep[i]==';'){break;}
	}
	rep[i+1]='\0';
}

bool tranceiver_set_freq(int freq)
{
	char str[2][16];
	char rep[2][16];

	sprintf(str[0],"FA%11.11d;",freq);
	sprintf(str[1],"FB%11.11d;",freq);
	
	cssl_putstring(serial.port,str[0]);
	cssl_putstring(serial.port,str[1]);
	
	cssl_putstring(serial.port,"FA;");
	get_serial_rep(rep[0]);
	cssl_putstring(serial.port,"FB;");
	get_serial_rep(rep[1]);
	
	return (bool)((strcmp( rep[0], str[0] ) == 0) && (strcmp( rep[1], str[1] ) == 0));
}

void tranceiver_init()
{
	cssl_putstring(serial.port,"FR0;"); //Set receive on VFO_A
	cssl_putstring(serial.port,"FT0;"); //Set Transmit on VFO_A
	cssl_putstring(serial.port,"Q10;"); //Set USB mode
	if(!tranceiver_set_freq(FT8.Tranceiver_VFOA_Freq)){printf("Unable to communicate with tranceiver!");exit(-1);} //Set frequency stored
}

void tranceiver_rtx(bool ptt)
{
	if (ptt == _TX_){
		#if DEBUG
		printf("tranceiver ptt on\n");
		#endif
		cssl_putstring(serial.port,"TX;");
		} //Set receiving
	else if (ptt == _RX_){
		#if DEBUG
		printf("tranceiver ptt off\n");
		#endif
		cssl_putstring(serial.port,"RX;");
		}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//main
int main (int argc, char *argv[])
{
	int c;
	while ((c = getopt (argc, argv, "hbd:C:L:F:x:S:")) != -1)
		switch (c)
		{
			case 'h':
				printf ("clFT8 -d plughw:CARD=PCH,DEV=0 -C F4JJJ -L JN38 -F 14074000 -S /dev/ttyACM0 -x 1 -b\n"
				"-x for set filter\n0 random (default)\n1 best decode score\n2 max distance\n3 min distance\n-b for console beep on log\n"
				"Red color: info, Blue: CQ finded,Magenta: filter by info missing or non standard message or call already made or Empty callsign\n");
				exit(0);
				break;
			case 'b':
				FT8.beep_on_log = 1;
				break;
				return 1;
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
			case 'F':
				FT8.Tranceiver_VFOA_Freq = atoi(optarg);
				break;
				return 1;
			case 'x':
				FT8.filter_on_cq = atoi(optarg);
				break;
				return 1;
			case 'S':
				strcpy( serial.pathname, optarg );
				break;
				return 1;
			default:
				abort ();
		}
	
	latLonForGrid(FT8.Local_LOCATOR,FT8.Local_latlon);
	sprintf(FT8.QSO_RESPONSES[5], "CQ %s %s", FT8.Local_CALLSIGN, FT8.Local_LOCATOR);
	
	log_FT8_open_callsigntable_ht();

	capture_audioInit();
	playback_audioInit();
	
	printf("Starting with this:\n"
		"-set Freq to %d\n"
		"-your callsing is %s\n"
		"-your locator is %s\n"
		"-TRX serial port is %s\n"
		"-Sound device is %s\n"
		"-CQ filter methode %d\n"
		"-Beep on log %d\n",
		FT8.Tranceiver_VFOA_Freq,FT8.Local_CALLSIGN,FT8.Local_LOCATOR,serial.pathname,sound.capture_sound_device,FT8.filter_on_cq,FT8.beep_on_log);
		
	if(FT8.beep_on_log){putchar('\07');putchar('\a');}
		
	int serres = serial_init();
	if(serres==-1){printf("Could not open serial port.");}
	else{tranceiver_init();}
	
	pthread_t thread_RX;
	pthread_create(&thread_RX, NULL, Thread_RX, NULL);
	
	pthread_t thread_TX;
	pthread_create(&thread_TX, NULL, Thread_TX, NULL);
	
	while (1)
	{
		usleep(500000);
		advance_cursor();
		if(FT8.infos_to_log[0] != '\0'){log_FT8_QSO();log_FT8_log_to_callsigntable_ht();}
		#if DEBUG
		printf("wait\n");
		#endif
	}
}
