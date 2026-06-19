#ifndef GCFT8_APP_STATE_H
#define GCFT8_APP_STATE_H

#include <stdbool.h>
#include <pthread.h>
#include <time.h>
#include <alsa/asoundlib.h>

#include "vendor/cssl/cssl.h"

#define capture_sound_device_ENV (char*)"RTSPECCY_CAPTURE_DEVICE"
#define playback_sound_device_ENV (char*)"RTSPECCY_PLAYBACK_DEVICE"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef enum
{
    GCFT8_TRX_RX = 0,
    GCFT8_TRX_TX = 1
} gcft8_trx_state_t;

#define GCFT8_MAX_QSO_SESSIONS 128
#define GCFT8_QSO_SESSION_TTL_SEC 900
#define GCFT8_MAX_SAME_TX_REPEATS 2
#define GCFT8_SNR_INVALID (-2147483647 - 1)

typedef struct
{
	bool in_use;
	bool logged;
	bool log_after_tx_73;
	bool log_pending;
	bool log_in_progress;

	char callsign[20];
	char locator[10];

	float frequency_hz;

	int last_rx_snr;
	int report_sent_snr;

	int next_tx_seq;
	int last_tx_seq;
	int same_tx_repeat_count;

	time_t started_at;
	time_t ended_at;
	time_t last_seen_at;
	time_t last_progress_at;
	time_t last_tx_at;
} gcft8_qso_session_t;

//gcFT8
typedef struct
{
	char Local_CALLSIGN[20];
	char Local_LOCATOR[10];
	float Local_latlon[2];

	gcft8_qso_session_t QSO_sessions[GCFT8_MAX_QSO_SESSIONS];
	
	gcft8_trx_state_t TRX_status;
	pthread_mutex_t TRX_status_lock;
	pthread_cond_t RX_status_cond;
	pthread_cond_t TX_status_cond;
	
	int Tranceiver_VFOA_Freq;
	
	char log_file_name[128];
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

#endif
