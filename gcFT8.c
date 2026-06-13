#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <alsa/asoundlib.h>
#include <pthread.h>
#include <sys/random.h>
#include <sys/time.h>

#include "ft8/decode.h"
#include "ft8/constants.h"
#include "ft8/encode.h"
#include "ft8/message.h"

#include "fft/kiss_fftr.h"

#include "common/monitor.h"

#include "serial/cssl.h"

#include "gcFT8.h"

#include "hash/hash.h"


#define FT8_TOKEN_TEXT_SIZE 25
#define CALLSIGN_HASH_CACHE_SIZE 256

#define FT8_FILTER_LISTEN_ONLY 0
#define FT8_FILTER_RANDOM_CQ 1
#define FT8_FILTER_BEST_DECODE_SCORE 2
#define FT8_FILTER_MAX_DISTANCE 3
#define FT8_FILTER_MIN_DISTANCE 4
#define FT8_FILTER_MAX_SNR 5
#define FT8_FILTER_MIN_SNR 6

#define FT8_SYMBOL_BT 2.0f ///< symbol smoothing filter bandwidth factor (BT)
#define FT4_SYMBOL_BT 1.0f ///< symbol smoothing filter bandwidth factor (BT)

#define FT8_TX_LEAD_SILENCE_SEC 0.5f
#define FT8_TX_TAIL_SILENCE_SEC 0.1f
#define FT4_TX_LEAD_SILENCE_SEC 0.3f
#define FT4_TX_TAIL_SILENCE_SEC 0.1f

typedef enum
{
	GCFT8_MODE_FT8,
	GCFT8_MODE_FT4
} gcft8_mode_t;

typedef struct
{
	const char* band;
	int ft8_frequency_hz;
	int ft4_frequency_hz;
} gcft8_band_frequency_t;

typedef struct
{
	gcft8_mode_t mode;
	const char* name;
	ftx_protocol_t protocol;
	int num_tones;
	float symbol_period;
	float symbol_bt;
	float slot_time;
	float tx_lead_silence;
	float tx_tail_silence;
	float rx_capture_time;
} gcft8_mode_config_t;

static const gcft8_band_frequency_t gcft8_band_frequencies[] = {
	{ "80",  3573000,  3575000 },
	{ "60",  5357000,  5357000 },
	{ "40",  7074000,  7047500 },
	{ "30", 10136000, 10140000 },
	{ "20", 14074000, 14080000 },
	{ "17", 18100000, 18104000 },
	{ "15", 21074000, 21140000 },
	{ "12", 24915000, 24919000 },
	{ "11", 27245000,        0 },
	{ "10", 28074000, 28180000 }
};

static const gcft8_mode_config_t gcft8_mode_configs[] = {
	{
		.mode = GCFT8_MODE_FT8,
		.name = "ft8",
		.protocol = FTX_PROTOCOL_FT8,
		.num_tones = FT8_NN,
		.symbol_period = FT8_SYMBOL_PERIOD,
		.symbol_bt = FT8_SYMBOL_BT,
		.slot_time = FT8_SLOT_TIME,
		.tx_lead_silence = FT8_TX_LEAD_SILENCE_SEC,
		.tx_tail_silence = FT8_TX_TAIL_SILENCE_SEC,
		.rx_capture_time = 13.6f
	},
	{
		.mode = GCFT8_MODE_FT4,
		.name = "ft4",
		.protocol = FTX_PROTOCOL_FT4,
		.num_tones = FT4_NN,
		.symbol_period = FT4_SYMBOL_PERIOD,
		.symbol_bt = FT4_SYMBOL_BT,
		.slot_time = FT4_SLOT_TIME,
		.tx_lead_silence = FT4_TX_LEAD_SILENCE_SEC,
		.tx_tail_silence = FT4_TX_TAIL_SILENCE_SEC,
		.rx_capture_time = FT4_SLOT_TIME - 0.4f
	}
};

typedef struct
{
	char callsign[12];
	uint32_t hash;
} callsign_hash_entry_t;

static callsign_hash_entry_t callsign_hash_cache[CALLSIGN_HASH_CACHE_SIZE];
static int callsign_hash_cache_size;
static volatile sig_atomic_t shutdown_requested;
static gcft8_mode_t gcft8_mode = GCFT8_MODE_FT8;
static bool gcft8_filter_has_band;
static char gcft8_filter_band[8];
static int gcft8_filter_frequency_mhz;

static bool gcft8_shutdown_requested(void)
{
	return shutdown_requested != 0;
}

static void gcft8_signal_handler(int signal_number)
{
	(void)signal_number;
	shutdown_requested = 1;
}

static void copy_text(char* dst, size_t dst_size, const char* src)
{
	if (dst_size == 0)
		return;

	if (src == NULL)
		src = "";

	size_t len = strlen(src);
	if (len >= dst_size)
		len = dst_size - 1;

	memcpy(dst, src, len);
	dst[len] = '\0';
}

static int grid_letter_index(char c)
{
	const char strtable[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	const char* found = strchr(strtable, toupper((unsigned char)c));

	if (found == NULL)
		return -1;

	return (int)(found - strtable);
}

static int grid_digit_index(char c)
{
	if (!isdigit((unsigned char)c))
		return -1;

	return c - '0';
}

static const gcft8_mode_config_t* gcft8_mode_config(gcft8_mode_t mode)
{
	for (size_t idx = 0; idx < sizeof(gcft8_mode_configs) / sizeof(gcft8_mode_configs[0]); ++idx)
	{
		if (gcft8_mode_configs[idx].mode == mode)
			return &gcft8_mode_configs[idx];
	}

	return &gcft8_mode_configs[0];
}

static const gcft8_mode_config_t* gcft8_current_mode_config(void)
{
	return gcft8_mode_config(gcft8_mode);
}

static bool gcft8_parse_mode(const char* value, gcft8_mode_t* mode)
{
	char normalized_mode[4];
	size_t len;

	if ((value == NULL) || (mode == NULL))
		return false;

	len = strlen(value);
	if ((len == 0) || (len >= sizeof(normalized_mode)))
		return false;

	for (size_t idx = 0; idx < len; ++idx)
		normalized_mode[idx] = (char)tolower((unsigned char)value[idx]);
	normalized_mode[len] = '\0';

	for (size_t idx = 0; idx < sizeof(gcft8_mode_configs) / sizeof(gcft8_mode_configs[0]); ++idx)
	{
		if (strcmp(normalized_mode, gcft8_mode_configs[idx].name) == 0)
		{
			*mode = gcft8_mode_configs[idx].mode;
			return true;
		}
	}

	return false;
}

static bool gcft8_frequency_for_band(const char* band, gcft8_mode_t mode, int* frequency_hz)
{
	char normalized_band[4];
	size_t len;

	if ((band == NULL) || (frequency_hz == NULL))
		return false;

	len = strlen(band);
	if ((len > 0) && ((band[len - 1] == 'm') || (band[len - 1] == 'M')))
		--len;

	if ((len == 0) || (len >= sizeof(normalized_band)))
		return false;

	for (size_t idx = 0; idx < len; ++idx)
	{
		if (!isdigit((unsigned char)band[idx]))
			return false;

		normalized_band[idx] = band[idx];
	}
	normalized_band[len] = '\0';

	for (size_t idx = 0; idx < sizeof(gcft8_band_frequencies) / sizeof(gcft8_band_frequencies[0]); ++idx)
	{
		if (strcmp(normalized_band, gcft8_band_frequencies[idx].band) == 0)
		{
			int frequency = (mode == GCFT8_MODE_FT4) ? gcft8_band_frequencies[idx].ft4_frequency_hz : gcft8_band_frequencies[idx].ft8_frequency_hz;
			if (frequency <= 0)
				return false;

			*frequency_hz = frequency;
			return true;
		}
	}

	return false;
}

static const char* gcft8_adif_mode_name(void)
{
	return (gcft8_mode == GCFT8_MODE_FT4) ? "FT4" : "FT8";
}

static bool gcft8_normalize_band_adif(const char* band, char* out, size_t out_size)
{
	char normalized_band[4];
	size_t len;

	if ((band == NULL) || (out == NULL) || (out_size == 0))
		return false;

	len = strlen(band);
	if ((len > 0) && ((band[len - 1] == 'm') || (band[len - 1] == 'M')))
		--len;

	if ((len == 0) || (len >= sizeof(normalized_band)))
		return false;

	for (size_t idx = 0; idx < len; ++idx)
	{
		if (!isdigit((unsigned char)band[idx]))
			return false;

		normalized_band[idx] = band[idx];
	}
	normalized_band[len] = '\0';

	if ((len + 2) > out_size)
		return false;

	snprintf(out, out_size, "%sM", normalized_band);
	return true;
}

static void gcft8_set_filter_band_context(const char* band)
{
	if (gcft8_normalize_band_adif(band, gcft8_filter_band, sizeof(gcft8_filter_band)))
	{
		gcft8_filter_has_band = true;
	}
}

static void gcft8_set_filter_frequency_context(int frequency_hz)
{
	gcft8_filter_has_band = false;
	gcft8_filter_band[0] = '\0';
	gcft8_filter_frequency_mhz = frequency_hz / 1000000;
}

static void gcft8_update_adif_log_filename(const char* local_callsign, char* log_file_name, size_t log_file_name_size)
{
	char safe_callsign[20];
	size_t out_idx = 0;

	if ((local_callsign == NULL) || (log_file_name == NULL) || (log_file_name_size == 0))
		return;

	for (size_t idx = 0; (local_callsign[idx] != '\0') && (out_idx + 1 < sizeof(safe_callsign)); ++idx)
	{
		unsigned char ch = (unsigned char)local_callsign[idx];
		if (isalnum(ch))
			safe_callsign[out_idx++] = (char)toupper(ch);
		else if ((ch == '-') || (ch == '_'))
			safe_callsign[out_idx++] = (char)ch;
		else
			safe_callsign[out_idx++] = '_';
	}

	if (out_idx == 0)
		copy_text(safe_callsign, sizeof(safe_callsign), "UNKNOWN");
	else
		safe_callsign[out_idx] = '\0';

	snprintf(log_file_name, log_file_name_size, "QSO_%s.adif", safe_callsign);
}

static bool ft8_parse_filter_mode(const char* value, int* filter_mode)
{
	int result = 0;

	if ((value == NULL) || (value[0] == '\0') || (filter_mode == NULL))
		return false;

	for (size_t idx = 0; value[idx] != '\0'; ++idx)
	{
		if (!isdigit((unsigned char)value[idx]))
			return false;

		result = (result * 10) + (value[idx] - '0');
	}

	if ((result < FT8_FILTER_LISTEN_ONLY) || (result > FT8_FILTER_MIN_SNR))
		return false;

	*filter_mode = result;
	return true;
}

static bool ft8_getrandom_u32(uint32_t* value)
{
	uint8_t* bytes = (uint8_t*)value;
	size_t remaining = sizeof(*value);

	if (value == NULL)
		return false;

	while (remaining > 0)
	{
		ssize_t rc = getrandom(bytes, remaining, 0);
		if (rc < 0)
		{
			if (errno == EINTR)
				continue;

			return false;
		}

		if (rc == 0)
			return false;

		bytes += rc;
		remaining -= (size_t)rc;
	}

	return true;
}

static uint32_t ft8_fallback_random_u32(void)
{
	static uint32_t state;

	if (state == 0)
	{
		state = (uint32_t)time(NULL) ^ (uint32_t)getpid() ^ (uint32_t)(uintptr_t)&state;
		if (state == 0)
			state = 0x6D2B79F5u;
	}

	state ^= state << 13;
	state ^= state >> 17;
	state ^= state << 5;
	return state;
}

static uint32_t ft8_random_u32(void)
{
	uint32_t value;

	if (ft8_getrandom_u32(&value))
		return value;

	return ft8_fallback_random_u32();
}

static int ft8_random_index(int count)
{
	uint32_t limit;
	uint32_t threshold;

	if (count <= 1)
		return 0;

	limit = (uint32_t)count;
	threshold = (uint32_t)(-limit % limit);

	for (;;)
	{
		uint32_t value = ft8_random_u32();
		if (value >= threshold)
			return (int)(value % limit);
	}
}

static void callsign_hash_cache_init(void)
{
	callsign_hash_cache_size = 0;
	memset(callsign_hash_cache, 0, sizeof(callsign_hash_cache));
}

static void callsign_hash_cache_insert(const char* callsign, uint32_t stored_hash)
{
	uint16_t hash10 = ((stored_hash & 0x3FFFFFu) >> 12) & 0x3FFu;
	int idx_hash = (hash10 * 23) % CALLSIGN_HASH_CACHE_SIZE;
	int first_idx = idx_hash;

	for (int probes = 0; probes < CALLSIGN_HASH_CACHE_SIZE; ++probes)
	{
		if (callsign_hash_cache[idx_hash].callsign[0] == '\0')
		{
			copy_text(callsign_hash_cache[idx_hash].callsign, sizeof(callsign_hash_cache[idx_hash].callsign), callsign);
			callsign_hash_cache[idx_hash].hash = stored_hash;
			++callsign_hash_cache_size;
			return;
		}

		if (((callsign_hash_cache[idx_hash].hash & 0x3FFFFFu) == (stored_hash & 0x3FFFFFu)) && (strcmp(callsign_hash_cache[idx_hash].callsign, callsign) == 0))
		{
			callsign_hash_cache[idx_hash].hash = stored_hash;
			return;
		}

		idx_hash = (idx_hash + 1) % CALLSIGN_HASH_CACHE_SIZE;
	}

	copy_text(callsign_hash_cache[first_idx].callsign, sizeof(callsign_hash_cache[first_idx].callsign), callsign);
	callsign_hash_cache[first_idx].hash = stored_hash;
}

static void callsign_hash_cache_cleanup(uint8_t max_age)
{
	callsign_hash_entry_t old_cache[CALLSIGN_HASH_CACHE_SIZE];

	memcpy(old_cache, callsign_hash_cache, sizeof(old_cache));
	callsign_hash_cache_init();

	for (int idx = 0; idx < CALLSIGN_HASH_CACHE_SIZE; ++idx)
	{
		if (old_cache[idx].callsign[0] != '\0')
		{
			uint8_t age = (uint8_t)(old_cache[idx].hash >> 24);
			if (age <= max_age)
			{
				uint32_t stored_hash = (((uint32_t)age + 1u) << 24) | (old_cache[idx].hash & 0x3FFFFFu);
				callsign_hash_cache_insert(old_cache[idx].callsign, stored_hash);
			}
		}
	}
}

static void callsign_hash_cache_save(const char* callsign, uint32_t hash)
{
	callsign_hash_cache_insert(callsign, hash & 0x3FFFFFu);
}

static bool callsign_hash_cache_lookup(ftx_callsign_hash_type_t hash_type, uint32_t hash, char* callsign)
{
	uint8_t hash_shift = (hash_type == FTX_CALLSIGN_HASH_10_BITS) ? 12 : (hash_type == FTX_CALLSIGN_HASH_12_BITS ? 10 : 0);
	uint16_t hash10 = (hash >> (12 - hash_shift)) & 0x3FFu;
	int idx_hash = (hash10 * 23) % CALLSIGN_HASH_CACHE_SIZE;

	for (int probes = 0; probes < CALLSIGN_HASH_CACHE_SIZE; ++probes)
	{
		if (callsign_hash_cache[idx_hash].callsign[0] == '\0')
			break;

		if (((callsign_hash_cache[idx_hash].hash & 0x3FFFFFu) >> hash_shift) == hash)
		{
			copy_text(callsign, 12, callsign_hash_cache[idx_hash].callsign);
			return true;
		}

		idx_hash = (idx_hash + 1) % CALLSIGN_HASH_CACHE_SIZE;
	}

	callsign[0] = '\0';
	return false;
}

static ftx_callsign_hash_interface_t callsign_hash_if = {
	.lookup_hash = callsign_hash_cache_lookup,
	.save_hash = callsign_hash_cache_save
};

static float db_to_power(float db)
{
	return powf(10.0f, db / 10.0f);
}

static int round_to_int(float value)
{
	return (int)((value >= 0.0f) ? value + 0.5f : value - 0.5f);
}

static int gcft8_estimate_snr(const ftx_waterfall_t* wf, const ftx_candidate_t* cand, const ftx_message_t* message)
{
	uint8_t tones[FT4_NN];
	int num_symbols;
	int num_tones;
	float signal_power = 0.0f;
	float noise_power = 0.0f;
	int signal_count = 0;
	int noise_count = 0;

	if (wf->protocol == FTX_PROTOCOL_FT4)
	{
		num_symbols = FT4_NN;
		num_tones = 4;
		ft4_encode(message->payload, tones);
	}
	else
	{
		num_symbols = FT8_NN;
		num_tones = 8;
		ft8_encode(message->payload, tones);
	}

	if (cand->freq_offset < 0 || cand->freq_offset + num_tones > wf->num_bins)
		return -99;

	for (int symbol = 0; symbol < num_symbols; ++symbol)
	{
		int block = cand->time_offset + symbol;
		if (block < 0 || block >= wf->num_blocks)
			continue;

		int offset = block;
		offset = (offset * wf->time_osr) + cand->time_sub;
		offset = (offset * wf->freq_osr) + cand->freq_sub;
		offset = (offset * wf->num_bins) + cand->freq_offset;

		const WF_ELEM_T* bins = wf->mag + offset;
		int active_tone = tones[symbol];
		if (active_tone < 0 || active_tone >= num_tones)
			continue;

		signal_power += db_to_power(WF_ELEM_MAG(bins[active_tone]));
		++signal_count;

		for (int tone = 0; tone < num_tones; ++tone)
		{
			if (tone == active_tone)
				continue;

			noise_power += db_to_power(WF_ELEM_MAG(bins[tone]));
			++noise_count;
		}
	}

	if (signal_count == 0 || noise_count == 0 || signal_power <= 0.0f || noise_power <= 0.0f)
		return -99;

	float signal_db = 10.0f * log10f(signal_power / signal_count);
	float noise_db = 10.0f * log10f(noise_power / noise_count);
	return round_to_int(signal_db - noise_db - 26.0f);
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Global

/* Global FT8 info. */


FT8info FT8 = {
	.Local_CALLSIGN = {'F','4','J','J','J',0},
	.Local_LOCATOR = {'J','N','3','8',0},
	.TX_enable = 1,
	
	.QSO_dist_CALLSIGN = "",
	.QSO_dist_LOCATOR = "",
	.QSO_dist_MESSAGE = "",
	.QSO_dist_SNR = 0,
	.QSO_dist_FREQUENCY = 0,
	
	.QSO_RESPONSES = {{0}},
	.QSO_Index_to_rep=-1,
	
	.TRX_status = _RX_,
	.TRX_status_lock = PTHREAD_MUTEX_INITIALIZER,
	.RX_status_cond = PTHREAD_COND_INITIALIZER,
	.TX_status_cond = PTHREAD_COND_INITIALIZER,
	
	.Tranceiver_VFOA_Freq = 14074000,
	
	.log_file_name = "QSO_F4JJJ.adif",
	.infos_to_log = {0},
	.log_dist_CALLSIGN_for_filter = "",
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

static bool ft8_listen_only_filter_active(void)
{
	return FT8.filter_on_cq == FT8_FILTER_LISTEN_ONLY;
}

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

static void clear_status_line(void)
{
	if (isatty(STDOUT_FILENO))
	{
		printf("\r%24s\r", "");
		fflush(stdout);
	}
}

void advance_cursor(float slot_time) {
  static int pos=0;
  char cursor[4]={'/','-','\\','|'};
  double slot = slot_time;
  double within_slot;

  if (!isatty(STDOUT_FILENO))
    return;

  if (slot <= 0.0)
    slot = 15.0;

  within_slot = fmod(now(), slot);
  if (within_slot < 0.0)
    within_slot += slot;

  if (fabs(slot - floor(slot + 0.5)) < 0.001)
  {
    int slot_seconds = (int)(slot + 0.5);
    int current_second = (int)floor(within_slot) + 1;

    if (current_second < 1)
      current_second = 1;
    if (current_second > slot_seconds)
      current_second = slot_seconds;

    printf("\r%c %02d/%02d", cursor[pos], current_second, slot_seconds);
  }
  else
  {
    printf("\r%c %04.1f/%04.1f", cursor[pos], within_slot, slot);
  }

  fflush(stdout);
  pos = (pos+1) % 4;
}

bool wait_for_slot_start(float slot_time){
	double slot = slot_time;
	double target_time;

	if (slot <= 0.0)
		return !gcft8_shutdown_requested();

	target_time = (floor(now() / slot) + 1.0) * slot;

	while (!gcft8_shutdown_requested())
	{
		double wait_time = target_time - now();
		useconds_t sleep_time;

		if (wait_time <= 0.0)
			return true;

		sleep_time = (useconds_t)(wait_time * 1000000.0);
		if (sleep_time > 100000)
			sleep_time = 100000;

		usleep(sleep_time);
	}

	return false;
}

void printDateTime_log(){
	clear_status_line();
	const gcft8_mode_config_t* mode_cfg = gcft8_current_mode_config();
	long long slot_ms = (long long)(mode_cfg->slot_time * 1000.0f + 0.5f);
	long long now_ms = (long long)(now() * 1000.0);
	long long slot_start_ms = (now_ms / slot_ms) * slot_ms;
	time_t t = (time_t)(slot_start_ms / 1000);
	int ms = (int)(slot_start_ms % 1000);
	struct tm tm = *gmtime(&t);
	if (ms == 0)
		printf("%d-%02d-%02d %02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	else
		printf("%d-%02d-%02d %02d:%02d:%02d.%03d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, ms);
}

void printDateTime(){
	clear_status_line();
	time_t t = time(NULL);	struct tm tm = *gmtime(&t);	int sec = (int)((int)tm.tm_sec);
	printf("%d-%02d-%02d %02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, sec);
}

void printDateTime_ms(){
	clear_status_line();
	struct timeval tv;
	gettimeofday(&tv, 0);
	time_t t = tv.tv_sec;
	struct tm tm = *gmtime(&t);
	printf("%d-%02d-%02d %02d:%02d:%02d.%03ld", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, (long)(tv.tv_usec / 1000));
}

void latLonForGrid(char * grid, float * latlon) {
	int lat_field;
	int lon_field;
	int lat_square;
	int lon_square;
	
	latlon[0]=-1;
	latlon[1]=-1;
	
	if(strlen(grid)==4){
		if(isalpha((unsigned char)grid[0]) && isalpha((unsigned char)grid[1]) && isdigit((unsigned char)grid[2]) && isdigit((unsigned char)grid[3])){
			lat_field = grid_letter_index(grid[1]);
			lon_field = grid_letter_index(grid[0]);
			lat_square = grid_digit_index(grid[3]);
			lon_square = grid_digit_index(grid[2]);
			if (lat_field < 0 || lon_field < 0 || lat_square < 0 || lon_square < 0)
				return;

			latlon[0] = lat_field * 10;               // 2nd digit: 10deg latitude slot.
			latlon[1] = lon_field * 20;               // 1st digit: 20deg longitude slot.
			latlon[0] += lat_square * 1;              // 4th digit: 1deg latitude slot.
			latlon[1] += lon_square * 2;              // 3rd digit: 2deg longitude slot.
			latlon[0] += 0.5 * 1;
			latlon[1] += 0.5 * 2;
			latlon[0] -= 90;                                                  // Locator lat/lon origin shift.
			latlon[1] -= 180;
			}
	}
};


float latLonDist(float * latlonA, float * latlonB){
	double value = sin(latlonA[0] * M_PI / 180.0f) * sin(latlonB[0] * M_PI / 180.0f) + cos(latlonA[0] * M_PI / 180.0f) * cos(latlonB[0] * M_PI / 180.0f) * cos((latlonA[1] * M_PI / 180.0f) - (latlonB[1] * M_PI / 180.0f));
	if (value > 1.0)
		value = 1.0;
	else if (value < -1.0)
		value = -1.0;

	return (float)(acos(value)*6371.0);
}

int count_occur_str(char * s, const char * c){
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

void unpackFT8mess(const char * message_text, char * unpackeds0, char * unpackeds1, char * unpackeds2){
	char message_copy[FTX_MAX_MESSAGE_LENGTH];
	char * NullToken;
	
	copy_text(message_copy, sizeof(message_copy), message_text);

	NullToken = strtok(message_copy, " ");
	if(NullToken != 0){
		copy_text(unpackeds0, FT8_TOKEN_TEXT_SIZE, NullToken);
	}else{unpackeds0[0] = '\0';}
	
	NullToken = strtok(NULL, " ");
	if(NullToken != 0){
		copy_text(unpackeds1, FT8_TOKEN_TEXT_SIZE, NullToken);
	}else{unpackeds1[0] = '\0';}
	
	NullToken = strtok(NULL, " ");
	if(NullToken != 0){
		copy_text(unpackeds2, FT8_TOKEN_TEXT_SIZE, NullToken);
	}else{unpackeds2[0] = '\0';}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Audio



/* Open and init the default recording device. */
void capture_audioInit(void)
{ 
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
	if (sound.capture_handle != NULL)
	{
		snd_pcm_drop(sound.capture_handle);
		snd_pcm_close(sound.capture_handle);
		sound.capture_handle = NULL;
	}
}

/* Open and init the default playback device. */
void playback_audioInit(void)
{ 
	int err;
	snd_pcm_hw_params_t *playback_hw_params;
	snd_pcm_sw_params_t *playback_sw_params;
	snd_pcm_uframes_t period_frames = (snd_pcm_uframes_t)sound.playback_buffer_frames;
	snd_pcm_uframes_t buffer_frames = period_frames * 8;

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

	if ((err = snd_pcm_hw_params_set_period_size_near(sound.playback_handle, playback_hw_params, &period_frames, 0)) < 0) {
		fprintf (stderr, "cannot set playback period size (%s)\n", snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_set_buffer_size_near(sound.playback_handle, playback_hw_params, &buffer_frames)) < 0) {
		fprintf (stderr, "cannot set playback buffer size (%s)\n", snd_strerror (err));
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

	if ((err = snd_pcm_hw_params_get_period_size(playback_hw_params, &period_frames, 0)) < 0) {
		fprintf (stderr, "cannot get playback period size (%s)\n", snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_hw_params_get_buffer_size(playback_hw_params, &buffer_frames)) < 0) {
		fprintf (stderr, "cannot get playback buffer size (%s)\n", snd_strerror (err));
		exit (1);
	}

	sound.playback_buffer_frames = (int)period_frames;
		
	snd_pcm_hw_params_free (playback_hw_params);
	
	#if DEBUG
	fprintf(stdout, "playback_hw_params freed\n");
	#endif

	if ((err = snd_pcm_sw_params_malloc (&playback_sw_params)) < 0) {
		fprintf (stderr, "cannot allocate playback software parameter structure (%s)\n", snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_sw_params_current(sound.playback_handle, playback_sw_params)) < 0) {
		fprintf (stderr, "cannot initialize playback software parameters (%s)\n", snd_strerror (err));
		exit (1);
	}

	snd_pcm_uframes_t start_threshold = period_frames;
	if ((err = snd_pcm_sw_params_set_start_threshold(sound.playback_handle, playback_sw_params, start_threshold)) < 0) {
		fprintf (stderr, "cannot set playback start threshold (%s)\n", snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_sw_params_set_avail_min(sound.playback_handle, playback_sw_params, period_frames)) < 0) {
		fprintf (stderr, "cannot set playback avail min (%s)\n", snd_strerror (err));
		exit (1);
	}

	if ((err = snd_pcm_sw_params(sound.playback_handle, playback_sw_params)) < 0) {
		fprintf (stderr, "cannot set playback software parameters (%s)\n", snd_strerror (err));
		exit (1);
	}

	snd_pcm_sw_params_free(playback_sw_params);
		
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
	if (sound.playback_handle != NULL)
	{
		snd_pcm_drop(sound.playback_handle);
		snd_pcm_close(sound.playback_handle);
		sound.playback_handle = NULL;
	}
}

static bool playback_prepare_for_tx(void)
{
	int err = snd_pcm_prepare(sound.playback_handle);
	if (err < 0)
	{
		fprintf(stderr, "cannot prepare playback interface for TX (%s)\n", snd_strerror(err));
		return false;
	}

	return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FT8_lib encode

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
//Thread FT8
float getFrame(char *buffer, int i)
{
	uint16_t lo = (uint8_t)buffer[2 * i];
	uint16_t hi = (uint8_t)buffer[(2 * i) + 1];
	int16_t sample = (int16_t)(lo | (hi << 8));
	return sample / 32768.0f;
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
	const gcft8_mode_config_t* mode_cfg = gcft8_current_mode_config();
	const int num_samples = (int)(mode_cfg->rx_capture_time * sound.capture_sound_rate);
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
		.protocol = mode_cfg->protocol
	};
	
	snd_pcm_sframes_t rc;
	
	while(!gcft8_shutdown_requested()){
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
			char *raw_data = (char *)malloc((size_t)mon.block_size * 2);
			if (raw_data == NULL)
			{
				monitor_free(&mon);
				fprintf(stderr, "Out of memory while allocating RX buffer\n");
				exit(1);
			}
		
			if (!wait_for_slot_start(mode_cfg->slot_time))
			{
				free(raw_data);
				monitor_free(&mon);
				break;
			}
			
			snd_pcm_reset(sound.capture_handle);
			
			for (int frame_pos = 0; !gcft8_shutdown_requested() && frame_pos + mon.block_size <= num_samples; frame_pos += mon.block_size)
			{
				// Process the waveform data frame by frame - you could have a live loop here with data from an audio device

				int frames_read = 0;
				while (!gcft8_shutdown_requested() && frames_read < mon.block_size)
				{
					rc = snd_pcm_readi(sound.capture_handle, raw_data + (frames_read * 2), (snd_pcm_uframes_t)(mon.block_size - frames_read));
					if (rc < 0)
					{
						rc = snd_pcm_recover(sound.capture_handle, (int)rc, 1);
						if (rc < 0)
							break;
						continue;
					}
					if (rc == 0)
						continue;

					frames_read += (int)rc;
				}

				if (frames_read == mon.block_size)
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
			ftx_candidate_t candidate_list[kMax_candidates];
			int num_candidates = ftx_find_candidates(&mon.wf, kMax_candidates, candidate_list, kMin_score);

			//Creat array for analyse
			char AnalyseArray[kMax_candidates][3][FT8_TOKEN_TEXT_SIZE];
			float  AnalyseArrayFreqInfo[kMax_candidates];
			int  AnalyseArraySNRInfo[kMax_candidates];
			
			int countanalyse=0;

			// Hash table for decoded messages (to check for duplicates)
			int num_decoded = 0;
			ftx_message_t decoded[kMax_decoded_messages];
			ftx_message_t* decoded_hashtable[kMax_decoded_messages];

			// Initialize hash table pointers
			for (int i = 0; i < kMax_decoded_messages; ++i)
			{
				decoded_hashtable[i] = NULL;
			}

			// Go over candidates and attempt to decode messages
			for (int idx = 0; !gcft8_shutdown_requested() && idx < num_candidates; ++idx)
			{
				
				const ftx_candidate_t* cand = &candidate_list[idx];
				if (cand->score < kMin_score)
					continue;

				float freq_hz = (mon.min_bin + cand->freq_offset + (float)cand->freq_sub / mon.wf.freq_osr) / mon.symbol_period;
				float time_sec = (cand->time_offset + (float)cand->time_sub / mon.wf.time_osr) * mon.symbol_period;

				ftx_message_t message;
				ftx_decode_status_t status;
				if (!ftx_decode_candidate(&mon.wf, cand, kLDPC_iterations, &message, &status))
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
					#endif
					continue;
				}

				char message_text[FTX_MAX_MESSAGE_LENGTH];
				ftx_message_offsets_t message_offsets;
				ftx_message_rc_t message_status = ftx_message_decode(&message, &callsign_hash_if, message_text, &message_offsets);
				if (message_status != FTX_MESSAGE_RC_OK)
				{
					#if DEBUG
					printf("Error while unpacking message: %d\n", (int)message_status);
					#endif
					continue;
				}

				int snr = gcft8_estimate_snr(&mon.wf, cand, &message);

				#if DEBUG
				printf( "Checking hash table for %4.1fs / %4.1fHz [%d]...\n", time_sec, freq_hz, cand->score);
				#endif
				int idx_hash = message.hash % kMax_decoded_messages;
				bool found_empty_slot = false;
				bool found_duplicate = false;
				int probes = 0;
				do
				{
					if (decoded_hashtable[idx_hash] == NULL)
					{
						#if DEBUG
						printf( "Found an empty slot\n");
						#endif
						found_empty_slot = true;
					}
					else if ((decoded_hashtable[idx_hash]->hash == message.hash) && (0 == memcmp(decoded_hashtable[idx_hash]->payload, message.payload, sizeof(message.payload))))
					{
						#if DEBUG
						printf( "Found a duplicate [%s]\n", message_text);
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
					++probes;
				} while (!found_empty_slot && !found_duplicate && probes < kMax_decoded_messages);

				if (!found_empty_slot && !found_duplicate)
					continue;

				if (found_empty_slot)
				{
					// Fill the empty hashtable slot
					memcpy(&decoded[idx_hash], &message, sizeof(message));
					decoded_hashtable[idx_hash] = &decoded[idx_hash];
					++num_decoded;
					
					printDateTime_log();
					
					if ((strncmp(message_text, FT8.Local_CALLSIGN, strlen(FT8.Local_CALLSIGN)) == 0) && (countanalyse>-1)) {
						
						countanalyse=-1;
						printf(" %d %3d %+4.2f %4.0f ~  \033[1;31m%s\033[0m\n", snr, cand->score, time_sec, freq_hz, message_text);
						if (!ft8_listen_only_filter_active())
						{
							unpackFT8mess(message_text,AnalyseArray[0][0],AnalyseArray[0][1],AnalyseArray[0][2]);

							pthread_mutex_lock(&FT8.TRX_status_lock);

							copy_text(FT8.QSO_dist_CALLSIGN, sizeof(FT8.QSO_dist_CALLSIGN), AnalyseArray[0][1]);
							copy_text(FT8.QSO_dist_MESSAGE, sizeof(FT8.QSO_dist_MESSAGE), AnalyseArray[0][2]);
							FT8.QSO_dist_FREQUENCY=freq_hz;
							FT8.QSO_dist_SNR = snr;
							FT8.QSO_Index_to_rep = -1;

							pthread_mutex_unlock(&FT8.TRX_status_lock);

							unlock_TX_thread();
						}
						
					}
					else if ((strncmp(message_text,"CQ",2) == 0) && (countanalyse>-1)) {
						unpackFT8mess(message_text,AnalyseArray[countanalyse][0],AnalyseArray[countanalyse][1],AnalyseArray[countanalyse][2]);
						
						AnalyseArrayFreqInfo[countanalyse]=freq_hz;
						AnalyseArraySNRInfo[countanalyse]=snr;
						
						if(ht_check(ht_callsigntable_for_filter,AnalyseArray[countanalyse][1]) || strlen(AnalyseArray[countanalyse][2])==0 || count_occur_str(message_text, " ") > 2){
							unpackFT8mess("",AnalyseArray[countanalyse][0],AnalyseArray[countanalyse][1],AnalyseArray[countanalyse][2]);
							printf(" %d %3d %+4.2f %4.0f ~  \033[1;35m%s\033[0m\n", snr, cand->score, time_sec, freq_hz, message_text);
						}
						else{
							countanalyse++;
							printf(" %d %3d %+4.2f %4.0f ~  \033[1;34m%s\033[0m\n", snr, cand->score, time_sec, freq_hz, message_text);
							}
							
					
					}else{
						printf(" %d %3d %+4.2f %4.0f ~  %s\n", snr, cand->score, time_sec, freq_hz, message_text);
					}
				}
			}
						
			if((countanalyse>0) && !ft8_listen_only_filter_active() && !gcft8_shutdown_requested()){
				int index_from_ope = 0;
				float dist = 0;
				int actusnr;
				
				switch (FT8.filter_on_cq)
				{
				
				case FT8_FILTER_RANDOM_CQ:
					index_from_ope = ft8_random_index(countanalyse);
					break;
					
				case FT8_FILTER_BEST_DECODE_SCORE:
					index_from_ope = 0;
					break;				
				
				case FT8_FILTER_MAX_DISTANCE:
					dist = 0;
					for(int i = 0; i<countanalyse;i++){
						if(strlen(AnalyseArray[i][2]) == 4){
							float latlonlocal[2];
							latLonForGrid(AnalyseArray[i][2],latlonlocal);
							float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
							#if DEBUG
							printf("%s lat:%f lon:%f dist:%f\n",AnalyseArray[i][1],latlonlocal[0],latlonlocal[1],new_dist);
							#endif
							if(dist < new_dist){index_from_ope=i;dist=new_dist;}
						}
					}				
					break;
					
				case FT8_FILTER_MIN_DISTANCE:
					dist = 6372;
					for(int i = 0; i<countanalyse;i++){
						if(strlen(AnalyseArray[i][2]) == 4){
							float latlonlocal[2];
							latLonForGrid(AnalyseArray[i][2],latlonlocal);
							float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
							#if DEBUG
							printf("%s lat:%f lon:%f dist:%f\n",AnalyseArray[i][1],latlonlocal[0],latlonlocal[1],new_dist);
							#endif
							if(dist > new_dist){index_from_ope=i;dist=new_dist;}
						}
					}				
					break;

				case FT8_FILTER_MAX_SNR:
					actusnr = -100;
					for(int i = 0; i<countanalyse;i++){
						if(AnalyseArraySNRInfo[i] > actusnr){index_from_ope=i;actusnr=AnalyseArraySNRInfo[i];}
					}				
					break;

				case FT8_FILTER_MIN_SNR:
					actusnr = 100;
					for(int i = 0; i<countanalyse;i++){
						if(AnalyseArraySNRInfo[i] < actusnr){index_from_ope=i;actusnr=AnalyseArraySNRInfo[i];}
					}				
					break;
				
								
				default:
					index_from_ope = ft8_random_index(countanalyse);
					break;

				}
				
				pthread_mutex_lock(&FT8.TRX_status_lock);
				
				copy_text(FT8.QSO_dist_CALLSIGN, sizeof(FT8.QSO_dist_CALLSIGN), AnalyseArray[index_from_ope][1]);
				copy_text(FT8.QSO_dist_LOCATOR, sizeof(FT8.QSO_dist_LOCATOR), AnalyseArray[index_from_ope][2]);
				FT8.QSO_Index_to_rep=0;
				FT8.QSO_dist_FREQUENCY=AnalyseArrayFreqInfo[index_from_ope];
				
				pthread_mutex_unlock(&FT8.TRX_status_lock);
				
				clear_status_line();
				printf("*Selected for new QSO: \033[1;31m%s at %s on %f hz (seq QSO on %d) \033[0m\n", FT8.QSO_dist_CALLSIGN,FT8.QSO_dist_LOCATOR,FT8.QSO_dist_FREQUENCY,FT8.QSO_Index_to_rep);
				
				unlock_TX_thread();
				
			}
			
			
			if(!num_decoded){printDateTime_log();printf(" N/A\n");}
			
			#if DEBUG
			printf( "Decoded %d messages\n", num_decoded);
			#endif
			callsign_hash_cache_cleanup(10);
			
			free(raw_data);
			monitor_free(&mon);
		
		}
		else{
			// printf( "Lock RX thread\n");
			pthread_mutex_lock(&FT8.TRX_status_lock);
			while (!gcft8_shutdown_requested() && FT8.TRX_status != _RX_)
				pthread_cond_wait(&FT8.RX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
		
	}

}

void gen_FT8_responses()
{
	snprintf(FT8.QSO_RESPONSES[0], sizeof(FT8.QSO_RESPONSES[0]), "%s %s %s", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN, FT8.Local_LOCATOR);
	snprintf(FT8.QSO_RESPONSES[1], sizeof(FT8.QSO_RESPONSES[1]), "%s %s %+d", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN, FT8.QSO_dist_SNR);
	snprintf(FT8.QSO_RESPONSES[2], sizeof(FT8.QSO_RESPONSES[2]), "%s %s R%+d", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN, FT8.QSO_dist_SNR);
	snprintf(FT8.QSO_RESPONSES[3], sizeof(FT8.QSO_RESPONSES[3]), "%s %s RRR", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN);
	snprintf(FT8.QSO_RESPONSES[4], sizeof(FT8.QSO_RESPONSES[4]), "%s %s 73", FT8.QSO_dist_CALLSIGN, FT8.Local_CALLSIGN);
}

void Reinit_FT8_QSO()
{
	FT8.QSO_Index_to_rep=-1;
	for(int i = 0; i<5;i++){
		FT8.QSO_RESPONSES[i][0] = '\0';
	}

	FT8.QSO_dist_CALLSIGN[0]=0;
	FT8.QSO_dist_LOCATOR[0]=0;
	FT8.QSO_dist_MESSAGE[0]=0;
	
	FT8.QSO_dist_SNR=0;
}

int get_seq_qso_to_rep(const char * mess, bool * flaglog)
{
	int rep = -1;
	size_t len;
	*flaglog=false;

	if (mess == NULL || mess[0] == '\0')
		return -1;

	len = strlen(mess);
	
	if ((mess[0] == '-') && (isdigit((unsigned char)mess[len-1])))
	{
		rep = 2;*flaglog=false;
	}
	
	if ((mess[0] == '+') && (isdigit((unsigned char)mess[len-1])))
	{
		rep = 2;*flaglog=false;
	}

	if ((mess[0] == 'R') && (isdigit((unsigned char)mess[len-1])))
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

static bool gcft8_append_text(char* dst, size_t dst_size, size_t* offset, const char* text)
{
	size_t len;

	if ((dst == NULL) || (offset == NULL) || (text == NULL) || (*offset >= dst_size))
		return false;

	len = strlen(text);
	if (len >= dst_size - *offset)
		return false;

	memcpy(dst + *offset, text, len + 1);
	*offset += len;
	return true;
}

static bool gcft8_append_adif_field(char* dst, size_t dst_size, size_t* offset, const char* name, const char* value)
{
	char field_header[64];
	int written;

	if ((value == NULL) || (value[0] == '\0'))
		return true;

	written = snprintf(field_header, sizeof(field_header), "<%s:%zu>", name, strlen(value));
	if ((written < 0) || ((size_t)written >= sizeof(field_header)))
		return false;

	return gcft8_append_text(dst, dst_size, offset, field_header) &&
		gcft8_append_text(dst, dst_size, offset, value) &&
		gcft8_append_text(dst, dst_size, offset, " ");
}

static void gcft8_format_adif_datetime(time_t t, char date[9], char time_on[7])
{
	struct tm tm = *gmtime(&t);
	strftime(date, 9, "%Y%m%d", &tm);
	strftime(time_on, 7, "%H%M%S", &tm);
}

static bool gcft8_build_adif_qso_record(char* dst, size_t dst_size, const char* callsign, const char* locator, int rst_sent, double rf_frequency_hz, time_t timestamp)
{
	size_t offset = 0;
	char qso_date[9];
	char time_on[7];
	char freq_mhz[16];
	char rst_sent_text[8];

	if ((dst == NULL) || (dst_size == 0) || (callsign == NULL) || (callsign[0] == '\0'))
		return false;

	dst[0] = '\0';
	gcft8_format_adif_datetime(timestamp, qso_date, time_on);
	snprintf(freq_mhz, sizeof(freq_mhz), "%.6f", rf_frequency_hz / 1000000.0);
	snprintf(rst_sent_text, sizeof(rst_sent_text), "%+03d", rst_sent);

	if (!gcft8_append_adif_field(dst, dst_size, &offset, "CALL", callsign))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "MODE", gcft8_adif_mode_name()))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "QSO_DATE", qso_date))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "TIME_ON", time_on))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "FREQ", freq_mhz))
		return false;
	if (gcft8_filter_has_band && !gcft8_append_adif_field(dst, dst_size, &offset, "BAND", gcft8_filter_band))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "GRIDSQUARE", locator))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "RST_SENT", rst_sent_text))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "STATION_CALLSIGN", FT8.Local_CALLSIGN))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "MY_GRIDSQUARE", FT8.Local_LOCATOR))
		return false;

	return gcft8_append_text(dst, dst_size, &offset, "<EOR>\n");
}

static bool gcft8_write_adif_header(FILE* fptr)
{
	char qso_date[9];
	char time_on[7];
	char created_timestamp[16];

	if (fptr == NULL)
		return false;

	gcft8_format_adif_datetime(time(NULL), qso_date, time_on);
	snprintf(created_timestamp, sizeof(created_timestamp), "%s %s", qso_date, time_on);

	fprintf(fptr, "Generated by gcFT8\n");
	fprintf(fptr, "<ADIF_VER:5>3.1.7 <PROGRAMID:5>gcFT8 <CREATED_TIMESTAMP:15>%s <EOH>\n", created_timestamp);
	return true;
}

static void gcft8_ensure_adif_log_file(void)
{
	FILE* fptr = fopen(FT8.log_file_name, "rb");
	long size = 0;

	if (fptr != NULL)
	{
		if ((fseek(fptr, 0, SEEK_END) == 0))
			size = ftell(fptr);
		fclose(fptr);
	}

	if ((fptr == NULL) || (size == 0))
	{
		fptr = fopen(FT8.log_file_name, "ab");
		if (fptr == NULL)
		{
			printf("Error with ADIF log file!\n");
			exit(1);
		}

		gcft8_write_adif_header(fptr);
		fclose(fptr);
	}
}

void log_adif_qso(void)
{
	char adif_record[sizeof(FT8.infos_to_log)];
	FILE *fptr;

	pthread_mutex_lock(&FT8.TRX_status_lock);
	copy_text(adif_record, sizeof(adif_record), FT8.infos_to_log);
	FT8.infos_to_log[0] = 0;
	pthread_mutex_unlock(&FT8.TRX_status_lock);

	if (adif_record[0] == '\0')
		return;

	gcft8_ensure_adif_log_file();
	fptr = fopen(FT8.log_file_name,"ab");
	if(fptr == NULL)
	{
		printf("Error with ADIF log file!\n");
		exit(1);
	}

	fprintf(fptr, "%s", adif_record);
	fclose(fptr);
	clear_status_line();
	printf("\033[1;32mLogged QSO to %s: %s\033[0m", FT8.log_file_name, adif_record);
	if(FT8.beep_on_log){putchar('\07');putchar('\a');}
}

static void gcft8_uppercase_text(char* text);

void log_qso_to_filter_table(void)
{
	char log_dist_CALLSIGN_for_filter[sizeof(FT8.log_dist_CALLSIGN_for_filter)];

	pthread_mutex_lock(&FT8.TRX_status_lock);
	copy_text(log_dist_CALLSIGN_for_filter, sizeof(log_dist_CALLSIGN_for_filter), FT8.log_dist_CALLSIGN_for_filter);
	FT8.log_dist_CALLSIGN_for_filter[0] = 0;
	pthread_mutex_unlock(&FT8.TRX_status_lock);

	if (log_dist_CALLSIGN_for_filter[0] == '\0')
		return;

	gcft8_uppercase_text(log_dist_CALLSIGN_for_filter);
	ht_insert(ht_callsigntable_for_filter, log_dist_CALLSIGN_for_filter);
}

typedef struct
{
	char call[32];
	char mode[8];
	char band[8];
	char freq[24];
} gcft8_adif_record_t;

static bool gcft8_ascii_equal_ignore_case(const char* a, const char* b)
{
	if ((a == NULL) || (b == NULL))
		return false;

	while ((*a != '\0') && (*b != '\0'))
	{
		if (toupper((unsigned char)*a) != toupper((unsigned char)*b))
			return false;
		++a;
		++b;
	}

	return (*a == '\0') && (*b == '\0');
}

static void gcft8_uppercase_text(char* text)
{
	if (text == NULL)
		return;

	for (size_t idx = 0; text[idx] != '\0'; ++idx)
		text[idx] = (char)toupper((unsigned char)text[idx]);
}

static void gcft8_adif_copy_value(char* dst, size_t dst_size, const char* src, size_t src_len)
{
	if ((dst == NULL) || (dst_size == 0))
		return;

	if (src_len >= dst_size)
		src_len = dst_size - 1;

	memcpy(dst, src, src_len);
	dst[src_len] = '\0';
}

static bool gcft8_adif_record_matches_filter(const gcft8_adif_record_t* record)
{
	if ((record == NULL) || (record->call[0] == '\0') || (record->mode[0] == '\0'))
		return false;

	if (!gcft8_ascii_equal_ignore_case(record->mode, gcft8_adif_mode_name()))
		return false;

	if (gcft8_filter_has_band)
	{
		if (record->band[0] != '\0')
			return gcft8_ascii_equal_ignore_case(record->band, gcft8_filter_band);

		if (record->freq[0] != '\0')
		{
			double freq_mhz = strtod(record->freq, NULL);
			return ((int)floor(freq_mhz)) == gcft8_filter_frequency_mhz;
		}

		return false;
	}

	if (record->freq[0] != '\0')
	{
		double freq_mhz = strtod(record->freq, NULL);
		return ((int)floor(freq_mhz)) == gcft8_filter_frequency_mhz;
	}

	return false;
}

static void gcft8_process_adif_record(const gcft8_adif_record_t* record)
{
	char callsign[sizeof(((gcft8_adif_record_t*)0)->call)];

	if (!gcft8_adif_record_matches_filter(record))
		return;

	copy_text(callsign, sizeof(callsign), record->call);
	gcft8_uppercase_text(callsign);
	ht_insert(ht_callsigntable_for_filter, callsign);
}

void load_qso_filter_from_adif(void)
{
	FILE *fptr;
	long file_size;
	char* content;
	size_t bytes_read;
	gcft8_adif_record_t record;

	ht_callsigntable_for_filter = ht_create_table();
	gcft8_ensure_adif_log_file();

	fptr = fopen(FT8.log_file_name,"rb");
	if(fptr == NULL)
	{
		printf("Error with ADIF log file!\n");
		exit(1);
	}

	if (fseek(fptr, 0, SEEK_END) != 0)
	{
		fclose(fptr);
		printf("Error reading ADIF log file!\n");
		exit(1);
	}

	file_size = ftell(fptr);
	if (file_size < 0)
	{
		fclose(fptr);
		printf("Error reading ADIF log file!\n");
		exit(1);
	}
	rewind(fptr);

	content = (char*)malloc((size_t)file_size + 1u);
	if (content == NULL)
	{
		fclose(fptr);
		fprintf(stderr, "Out of memory while reading ADIF log\n");
		exit(1);
	}

	bytes_read = fread(content, 1, (size_t)file_size, fptr);
	fclose(fptr);
	content[bytes_read] = '\0';

	memset(&record, 0, sizeof(record));
	for (size_t idx = 0; idx < bytes_read; )
	{
		if (content[idx] != '<')
		{
			++idx;
			continue;
		}

		size_t tag_end = idx + 1;
		while ((tag_end < bytes_read) && (content[tag_end] != '>'))
			++tag_end;

		if (tag_end >= bytes_read)
			break;

		char field_name[32];
		size_t descriptor_start = idx + 1;
		size_t descriptor_len = tag_end - descriptor_start;
		size_t field_len = 0;
		size_t name_len = 0;

		while ((name_len < descriptor_len) && (content[descriptor_start + name_len] != ':') && (name_len + 1 < sizeof(field_name)))
		{
			field_name[name_len] = (char)toupper((unsigned char)content[descriptor_start + name_len]);
			++name_len;
		}
		field_name[name_len] = '\0';

		if (strcmp(field_name, "EOR") == 0)
		{
			gcft8_process_adif_record(&record);
			memset(&record, 0, sizeof(record));
			idx = tag_end + 1;
			continue;
		}

		if (strcmp(field_name, "EOH") == 0)
		{
			idx = tag_end + 1;
			continue;
		}

		if ((name_len < descriptor_len) && (content[descriptor_start + name_len] == ':'))
		{
			field_len = (size_t)atoi(content + descriptor_start + name_len + 1);
		}

		size_t value_start = tag_end + 1;
		if ((field_len == 0) || (value_start + field_len > bytes_read))
		{
			idx = tag_end + 1;
			continue;
		}

		if (strcmp(field_name, "CALL") == 0)
			gcft8_adif_copy_value(record.call, sizeof(record.call), content + value_start, field_len);
		else if (strcmp(field_name, "MODE") == 0)
			gcft8_adif_copy_value(record.mode, sizeof(record.mode), content + value_start, field_len);
		else if (strcmp(field_name, "BAND") == 0)
			gcft8_adif_copy_value(record.band, sizeof(record.band), content + value_start, field_len);
		else if (strcmp(field_name, "FREQ") == 0)
			gcft8_adif_copy_value(record.freq, sizeof(record.freq), content + value_start, field_len);

		idx = value_start + field_len;
	}

	free(content);
	printf("QSO filter table initialised with %d entry from %s.\n", ht_callsigntable_for_filter->count, FT8.log_file_name);
	#if DEBUG
	print_table(ht_callsigntable_for_filter);
	#endif
}

void TX_FT8()
{
	while(!gcft8_shutdown_requested()){
		pthread_mutex_lock(&FT8.TRX_status_lock);
		bool stat = FT8.TRX_status;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if((stat == _TX_) && !gcft8_shutdown_requested()){
			
			bool flaglog = 0;
			char tx_message[sizeof(FT8.QSO_RESPONSES[0])];
			char tx_callsign[sizeof(FT8.QSO_dist_CALLSIGN)];
			float tx_frequency = 0.0f;
			int tx_index = -1;
			
			pthread_mutex_lock(&FT8.TRX_status_lock);
			gen_FT8_responses();
			
			if(FT8.QSO_Index_to_rep == -1){FT8.QSO_Index_to_rep=get_seq_qso_to_rep(FT8.QSO_dist_MESSAGE,&flaglog);}
			tx_index = FT8.QSO_Index_to_rep;
			if((tx_index>-1) && !gcft8_shutdown_requested()){
				copy_text(tx_message, sizeof(tx_message), FT8.QSO_RESPONSES[tx_index]);
				copy_text(tx_callsign, sizeof(tx_callsign), FT8.QSO_dist_CALLSIGN);
				tx_frequency = FT8.QSO_dist_FREQUENCY;
			}
			pthread_mutex_unlock(&FT8.TRX_status_lock);
	
			if((tx_index>-1) && !gcft8_shutdown_requested()){
				
				const gcft8_mode_config_t* mode_cfg = gcft8_current_mode_config();
				clear_status_line();
				printf("Resp to \033[1;31m %s \033[0m with seq \033[1;31m %d \033[0m mess \033[1;31m %s \033[0m\n", tx_callsign, tx_index, tx_message);
			
				ftx_message_t tx_msg;
				ftx_message_init(&tx_msg);
				ftx_message_rc_t rc = ftx_message_encode(&tx_msg, &callsign_hash_if, tx_message);
				if (rc != FTX_MESSAGE_RC_OK)
				{
					printf("Cannot parse message!\n");
					printf("RC = %d\n", (int)rc);
					unlock_RX_thread();
					continue;
				}

				#if DEBUG
				printf("Packed data: ");
				for (int j = 0; j < 10; ++j)
				{
					printf("%02x ", tx_msg.payload[j]);
				}
				printf("\n");
				#endif

				int num_tones = mode_cfg->num_tones;
				float symbol_period = mode_cfg->symbol_period;
				float symbol_bt = mode_cfg->symbol_bt;
				float slot_time = mode_cfg->slot_time;

				// Second, encode the binary message as a sequence of FSK tones
				uint8_t tones[num_tones];
				if (mode_cfg->protocol == FTX_PROTOCOL_FT4)
					ft4_encode(tx_msg.payload, tones);
				else
					ft8_encode(tx_msg.payload, tones);
				
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
				int num_lead_silence = (int)(0.5f + mode_cfg->tx_lead_silence * sample_rate);
				int num_tail_silence = (int)(0.5f + mode_cfg->tx_tail_silence * sample_rate);
				int num_total_samples = num_lead_silence + num_samples + num_tail_silence;
				float signal[num_total_samples];
				
				for (int i = 0; i < num_total_samples; i++)
				{
					signal[i] = 0.0f;
				}
				
				synth_gfsk(tones, num_tones, tx_frequency, symbol_bt, symbol_period, sample_rate, signal + num_lead_silence);
				
				int16_t * raw_data = (int16_t*)malloc(num_total_samples*sizeof(int16_t)); // num_samples * numChannels * bitsPerSample / 8;
				if (raw_data == NULL)
				{
					fprintf(stderr, "Out of memory while allocating TX buffer\n");
					exit(1);
				}

				for (int i = 0; i < num_total_samples; i++)
				{
					float x = signal[i];
					if (x > 1.0)
					x = 1.0;
					else if (x < -1.0)
					x = -1.0;
					raw_data[i] = (int16_t)(x * 32767.0f);
				}

				snd_pcm_sframes_t count;

				if (!playback_prepare_for_tx())
				{
					free(raw_data);
					unlock_RX_thread();
					continue;
				}
				
				if (!wait_for_slot_start(slot_time))
				{
					free(raw_data);
					tranceiver_rtx(_RX_);
					break;
				}

				printDateTime_ms();printf(" start send message\n");
				
				tranceiver_rtx(_TX_);
								
				count = 0;
				do
				{
					snd_pcm_uframes_t frame_length = (snd_pcm_uframes_t)(num_total_samples - count);
					snd_pcm_sframes_t frames;

					if (frame_length > (snd_pcm_uframes_t)sound.playback_buffer_frames)
						frame_length = (snd_pcm_uframes_t)sound.playback_buffer_frames;

					frames = snd_pcm_writei(sound.playback_handle, raw_data + count, frame_length);

					// If an error, try to recover from it
					if (frames < 0){frames = snd_pcm_recover(sound.playback_handle, (int)frames, 0);}
					if (frames < 0)
					{
						printf("Error playing wave: %s\n", snd_strerror((int)frames));
						break;
					}
					if (frames == 0)
						continue;
					
					// Update our pointer
					count += frames;
					
					

				} while (!gcft8_shutdown_requested() && count < num_total_samples);
				
				if ((count == num_total_samples) && !gcft8_shutdown_requested())
				{
					snd_pcm_drain(sound.playback_handle);
				}
				else
				{
					snd_pcm_drop(sound.playback_handle);
				}
				
				tranceiver_rtx(_RX_);

				printDateTime_ms();printf(" stop send message\n");

				free(raw_data);
			
			}
			
			if(flaglog)
			{
				pthread_mutex_lock(&FT8.TRX_status_lock);
				if(FT8.QSO_dist_CALLSIGN[0]!=0 && FT8.QSO_dist_LOCATOR[0]!=0 && FT8.QSO_dist_MESSAGE[0]!=0){
					if (!gcft8_build_adif_qso_record(FT8.infos_to_log, sizeof(FT8.infos_to_log), FT8.QSO_dist_CALLSIGN, FT8.QSO_dist_LOCATOR, FT8.QSO_dist_SNR, (double)FT8.Tranceiver_VFOA_Freq + (double)FT8.QSO_dist_FREQUENCY, time(NULL)))
						FT8.infos_to_log[0] = '\0';
					snprintf(FT8.log_dist_CALLSIGN_for_filter, sizeof(FT8.log_dist_CALLSIGN_for_filter), "%s",FT8.QSO_dist_CALLSIGN);
					//Empty QSO variable
					Reinit_FT8_QSO();
				}
				pthread_mutex_unlock(&FT8.TRX_status_lock);
			}

			unlock_RX_thread();
			
		}
		else{
			#if DEBUG
			printf( "Lock TX thread\n");
			#endif
			pthread_mutex_lock(&FT8.TRX_status_lock);
			while (!gcft8_shutdown_requested() && FT8.TRX_status != _TX_)
				pthread_cond_wait(&FT8.TX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
	}

	tranceiver_rtx(_RX_);
}

void * Thread_RX(void *arg) {
	(void)arg;
	RX_FT8();
	return NULL;
}

void * Thread_TX(void *arg) {
	(void)arg;
	TX_FT8();
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//serial and tranceiver commands

int serial_init(){
    cssl_start();
    serial.port=cssl_open(serial.pathname,NULL,0,serial.baud,serial.bits,serial.parity,serial.stopbits);
	if (!serial.port) {
		printf("Serial error %s\n",cssl_geterrormsg());
		return -1;
	}

	cssl_setflowcontrol(serial.port,serial.rtscts,serial.xonxoff);
	cssl_settimeout(serial.port,500);
	return 0;
}

void serial_deinit(void)
{
	if (serial.port != NULL)
	{
		cssl_putstring(serial.port,"RX;");
		cssl_close(serial.port);
		serial.port = NULL;
	}
	cssl_stop();
}

void get_serial_rep(char rep[16])
{
	int i;
	for(i=0;i<15;i++){
		rep[i]=((char)cssl_getchar(serial.port));
		if(rep[i]==';'){break;}
	}
	rep[(i < 15) ? i + 1 : 15]=0;
}

bool tranceiver_set_freq(int freq)
{
	char str[2][16];
	char rep[2][16];

	snprintf(str[0], sizeof(str[0]), "FA%11.11d;",freq);
	snprintf(str[1], sizeof(str[1]), "FB%11.11d;",freq);
	
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
	if (serial.port == NULL)
		return;

	cssl_putstring(serial.port,"FR0;"); //Set receive on VFO_A
	cssl_putstring(serial.port,"FT0;"); //Set Transmit on VFO_A
	cssl_putstring(serial.port,"Q10;"); //Set USB mode
	if(!tranceiver_set_freq(FT8.Tranceiver_VFOA_Freq)){printf("Unable to communicate with tranceiver!");exit(-1);} //Set frequency stored
}

void tranceiver_rtx(bool ptt)
{
	if (serial.port == NULL)
		return;

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

static void install_signal_handlers(void)
{
	struct sigaction action;

	memset(&action, 0, sizeof(action));
	action.sa_handler = gcft8_signal_handler;
	sigemptyset(&action.sa_mask);
	sigaction(SIGINT, &action, NULL);
	sigaction(SIGTERM, &action, NULL);
}

static void wake_workers_for_shutdown(void)
{
	pthread_mutex_lock(&FT8.TRX_status_lock);
	FT8.TRX_status = _RX_;
	pthread_cond_broadcast(&FT8.RX_status_cond);
	pthread_cond_broadcast(&FT8.TX_status_cond);
	pthread_mutex_unlock(&FT8.TRX_status_lock);

	if (sound.capture_handle != NULL)
		snd_pcm_drop(sound.capture_handle);
	if (sound.playback_handle != NULL)
		snd_pcm_drop(sound.playback_handle);
}

static void gcft8_cleanup(void)
{
	clear_status_line();
	tranceiver_rtx(_RX_);
	serial_deinit();
	playback_audioDeInit();
	capture_audioDeInit();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//main

static void print_usage(const char* program_name, FILE* stream)
{
	if ((program_name == NULL) || (program_name[0] == '\0'))
		program_name = "gcFT8";

	fprintf(stream,
		"Usage:\n"
		"  %s [options]\n"
		"\n"
		"Example:\n"
		"  %s --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --beep\n"
		"\n"
		"Options:\n"
		"  --help                     Show this help and exit\n"
		"  --mode <ft8|ft4>           Digital mode (default: ft8)\n"
		"  --sound-device <device>    ALSA capture and playback device, prefer plughw (default: default)\n"
		"  --callsign <callsign>      Your callsign (default: F4JJJ)\n"
		"  --locator <locator>        Your Maidenhead locator (default: JN38)\n"
		"  --frequency <hz>           TRX frequency in Hz, exclusive with --band\n"
		"  --band <band>              Mode-specific band frequency, exclusive with --frequency; suffix m is allowed\n"
		"  --serial-device <device>   Transceiver serial device (default: /dev/ttyACM0)\n"
		"  --filter <mode>            Operating/filter mode (default: 0, listen only)\n"
		"  --beep                     Enable console beep when a QSO is logged\n"
		"\n"
		"Filters:\n"
		"  0  Listen only, no automatic TX\n"
		"  1  Random CQ selection\n"
		"  2  Best decode score\n"
		"  3  Maximum distance\n"
		"  4  Minimum distance\n"
		"  5  Maximum SNR\n"
		"  6  Minimum SNR\n"
		"\n"
		"Bands:\n"
		"  Band       FT8 Hz    FT4 Hz\n",
		program_name,
		program_name);

	for (size_t idx = 0; idx < sizeof(gcft8_band_frequencies) / sizeof(gcft8_band_frequencies[0]); ++idx)
	{
		if (gcft8_band_frequencies[idx].ft4_frequency_hz > 0)
			fprintf(stream, "  %2sm  %8d  %8d\n", gcft8_band_frequencies[idx].band, gcft8_band_frequencies[idx].ft8_frequency_hz, gcft8_band_frequencies[idx].ft4_frequency_hz);
		else
			fprintf(stream, "  %2sm  %8d       n/a\n", gcft8_band_frequencies[idx].band, gcft8_band_frequencies[idx].ft8_frequency_hz);
	}

	fprintf(stream,
		"\n"
		"Display colors:\n"
		"  Red      Local station related message\n"
		"  Blue     CQ candidate\n"
		"  Magenta  Filtered CQ, missing info, non-standard message, already worked callsign or empty callsign\n");
}

enum
{
	CLI_OPTION_HELP = 1000,
	CLI_OPTION_BEEP,
	CLI_OPTION_SOUND_DEVICE,
	CLI_OPTION_CALLSIGN,
	CLI_OPTION_LOCATOR,
	CLI_OPTION_MODE,
	CLI_OPTION_FREQUENCY,
	CLI_OPTION_BAND,
	CLI_OPTION_FILTER,
	CLI_OPTION_SERIAL_DEVICE
};

int main (int argc, char *argv[])
{
	int c;
	bool frequency_option_used = false;
	bool band_option_used = false;
	char selected_band[16] = "";
	static const struct option long_options[] = {
		{ "help", no_argument, NULL, CLI_OPTION_HELP },
		{ "beep", no_argument, NULL, CLI_OPTION_BEEP },
		{ "sound-device", required_argument, NULL, CLI_OPTION_SOUND_DEVICE },
		{ "callsign", required_argument, NULL, CLI_OPTION_CALLSIGN },
		{ "locator", required_argument, NULL, CLI_OPTION_LOCATOR },
		{ "mode", required_argument, NULL, CLI_OPTION_MODE },
		{ "frequency", required_argument, NULL, CLI_OPTION_FREQUENCY },
		{ "band", required_argument, NULL, CLI_OPTION_BAND },
		{ "filter", required_argument, NULL, CLI_OPTION_FILTER },
		{ "serial-device", required_argument, NULL, CLI_OPTION_SERIAL_DEVICE },
		{ NULL, 0, NULL, 0 }
	};

	while ((c = getopt_long(argc, argv, "", long_options, NULL)) != -1)
		switch (c)
		{
			case CLI_OPTION_HELP:
				print_usage(argv[0], stdout);
				exit(0);
				break;
			case CLI_OPTION_BEEP:
				FT8.beep_on_log = 1;
				break;
				return 1;
			case CLI_OPTION_SOUND_DEVICE:
				sound.capture_sound_device = optarg;
				sound.playback_sound_device = optarg;
				break;
				return 1;
			case CLI_OPTION_CALLSIGN:
				copy_text(FT8.Local_CALLSIGN, sizeof(FT8.Local_CALLSIGN), optarg);
				break;
				return 1;
			case CLI_OPTION_LOCATOR:
				copy_text(FT8.Local_LOCATOR, sizeof(FT8.Local_LOCATOR), optarg);
				break;
				return 1;
			case CLI_OPTION_MODE:
				if (!gcft8_parse_mode(optarg, &gcft8_mode))
				{
					fprintf(stderr, "Invalid mode '%s'. Allowed modes: ft8, ft4.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
				return 1;
			case CLI_OPTION_FREQUENCY:
				frequency_option_used = true;
				FT8.Tranceiver_VFOA_Freq = atoi(optarg);
				break;
				return 1;
			case CLI_OPTION_BAND:
				band_option_used = true;
				copy_text(selected_band, sizeof(selected_band), optarg);
				break;
				return 1;
			case CLI_OPTION_FILTER:
				if (!ft8_parse_filter_mode(optarg, &FT8.filter_on_cq))
				{
					fprintf(stderr, "Invalid filter '%s'. Allowed filters: 0, 1, 2, 3, 4, 5, 6.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
				return 1;
			case CLI_OPTION_SERIAL_DEVICE:
				copy_text(serial.pathname, sizeof(serial.pathname), optarg);
				break;
				return 1;
			default:
				print_usage(argv[0], stderr);
				exit(1);
		}

	if (optind < argc)
	{
		fprintf(stderr, "Unexpected argument '%s'.\n", argv[optind]);
		print_usage(argv[0], stderr);
		exit(1);
	}

	if (frequency_option_used && band_option_used)
	{
		fprintf(stderr, "Use either --frequency or --band, not both.\n");
		print_usage(argv[0], stderr);
		exit(1);
	}

	if (band_option_used)
	{
		if (!gcft8_frequency_for_band(selected_band, gcft8_mode, &FT8.Tranceiver_VFOA_Freq))
		{
			fprintf(stderr, "Invalid %s band '%s'. Use --help to list supported bands.\n", gcft8_current_mode_config()->name, selected_band);
			print_usage(argv[0], stderr);
			exit(1);
		}
		gcft8_set_filter_band_context(selected_band);
		gcft8_filter_frequency_mhz = FT8.Tranceiver_VFOA_Freq / 1000000;
	}
	else if (frequency_option_used)
	{
		gcft8_set_filter_frequency_context(FT8.Tranceiver_VFOA_Freq);
	}
	else if (!frequency_option_used)
	{
		(void)gcft8_frequency_for_band("20", gcft8_mode, &FT8.Tranceiver_VFOA_Freq);
		gcft8_set_filter_band_context("20");
		gcft8_filter_frequency_mhz = FT8.Tranceiver_VFOA_Freq / 1000000;
	}

	install_signal_handlers();
	
	latLonForGrid(FT8.Local_LOCATOR,FT8.Local_latlon);
	snprintf(FT8.QSO_RESPONSES[5], sizeof(FT8.QSO_RESPONSES[5]), "CQ %s %s", FT8.Local_CALLSIGN, FT8.Local_LOCATOR);
	callsign_hash_cache_init();
	gcft8_update_adif_log_filename(FT8.Local_CALLSIGN, FT8.log_file_name, sizeof(FT8.log_file_name));
	load_qso_filter_from_adif();
	
	capture_audioInit();
	playback_audioInit();
	
	const gcft8_mode_config_t* startup_mode_cfg = gcft8_current_mode_config();
	printf("Starting with this:\n"
		"-mode is %s\n"
		"-set Freq to %d\n"
		"-your callsing is %s\n"
		"-your locator is %s\n"
		"-TRX serial port is %s\n"
		"-Sound device is %s\n"
		"-ADIF log file is %s\n"
		"-CQ filter methode %d\n"
		"-Beep on log %d\n",
		startup_mode_cfg->name,FT8.Tranceiver_VFOA_Freq,FT8.Local_CALLSIGN,FT8.Local_LOCATOR,serial.pathname,sound.capture_sound_device,FT8.log_file_name,FT8.filter_on_cq,FT8.beep_on_log);
		
	if(FT8.beep_on_log){putchar('\07');putchar('\a');}
		
	int serres = serial_init();
	if(serres==-1){printf("Could not open serial port.");}
	else{tranceiver_init();}
	
	pthread_t thread_RX;
	pthread_create(&thread_RX, NULL, Thread_RX, NULL);
	
	pthread_t thread_TX;
	pthread_create(&thread_TX, NULL, Thread_TX, NULL);
	
	while (!gcft8_shutdown_requested())
	{
		usleep(500000);
		advance_cursor(startup_mode_cfg->slot_time);
		pthread_mutex_lock(&FT8.TRX_status_lock);
		bool log_pending = FT8.infos_to_log[0] != 0;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if(log_pending){log_adif_qso();log_qso_to_filter_table();}
		#if DEBUG
		printf("wait\n");
		#endif
	}

	shutdown_requested = 1;
	wake_workers_for_shutdown();
	pthread_join(thread_RX, NULL);
	pthread_join(thread_TX, NULL);
	gcft8_cleanup();
	printf("Stopped cleanly.\n");
	return 0;
}
