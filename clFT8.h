#define capture_sound_device_ENV (char*)"RTSPECCY_CAPTURE_DEVICE"
#define playback_sound_device_ENV (char*)"RTSPECCY_CAPTURE_DEVICE"

#define _RX_ 0
#define _TX_ 1

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//clFT8
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
	float  QSO_dist_FREQUENCY;
	
	char QSO_RESPONSES[6][50];
	int QSO_Index_to_rep;
	
	bool TRX_status;
	pthread_mutex_t TRX_status_lock;
	pthread_cond_t RX_status_cond;
	pthread_cond_t TX_status_cond;
	
	
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

//ft8lib

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