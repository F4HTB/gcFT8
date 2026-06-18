#ifndef GCFT8_H
#define GCFT8_H

#include <stdbool.h>
#include <pthread.h>
#include <alsa/asoundlib.h>

#include "vendor/cssl/cssl.h"

#define capture_sound_device_ENV (char*)"RTSPECCY_CAPTURE_DEVICE"
#define playback_sound_device_ENV (char*)"RTSPECCY_PLAYBACK_DEVICE"

#define _RX_ 0
#define _TX_ 1

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//gcFT8
typedef struct
{
	char Local_CALLSIGN[20];
	char Local_LOCATOR[10];
	float Local_latlon[2];
	bool TX_enable;
	
	char QSO_dist_CALLSIGN[20];
	char QSO_dist_LOCATOR[10];
	char QSO_dist_MESSAGE[10];
	int  QSO_dist_SNR;
	float  QSO_dist_FREQUENCY;
	
	char QSO_RESPONSES[6][50];
	int QSO_Index_to_rep;
	
	bool TRX_status;
	pthread_mutex_t TRX_status_lock;
	pthread_cond_t RX_status_cond;
	pthread_cond_t TX_status_cond;
	
	int Tranceiver_VFOA_Freq;
	
	char log_file_name[128];
	char infos_to_log [512];
	char log_dist_CALLSIGN_for_filter[20];
	bool beep_on_log;
	
	int filter_on_cq;
	
	
} FT8info;

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

//serial

typedef struct
{
	cssl_t *port;
	char pathname[20];
	int rtscts;
	int xonxoff;
	int baud;
	int bits;
	int parity;
	int stopbits;
	int finished;
} serial_t;

//Functions

void tranceiver_rtx(bool ptt);

#endif
