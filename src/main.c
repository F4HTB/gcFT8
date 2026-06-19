#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>
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

#include "protocol/ftx/decode.h"
#include "protocol/ftx/constants.h"
#include "protocol/ftx/encode.h"
#include "protocol/ftx/message.h"

#include "dsp/monitor.h"
#include "dsp/gfsk.h"

#include "protocol/ft2/ft2_waveform.h"

#include "vendor/cssl/cssl.h"

#include "app/state.h"

#include "util/hash_table.h"


#define FT8_TOKEN_TEXT_SIZE 25
#define CALLSIGN_HASH_CACHE_SIZE 256
#define GCFT8_MAX_ONLY_PREFIXES 64
#define GCFT8_PREFIX_TEXT_SIZE 16
#define GCFT8_MAX_SP_TAGS 64
#define GCFT8_SP_TAG_TEXT_SIZE 8
#define GCFT8_MAX_LOCATOR_ZONES 64
#define GCFT8_TX_MESSAGE_TEXT_SIZE 50
#define GCFT8_ADIF_RECORD_TEXT_SIZE 512

#define FT8_FILTER_LISTEN_ONLY 0
#define FT8_FILTER_RANDOM_CQ 1
#define FT8_FILTER_BEST_DECODE_SCORE 2
#define FT8_FILTER_MAX_DISTANCE 3
#define FT8_FILTER_MIN_DISTANCE 4
#define FT8_FILTER_MAX_SNR 5
#define FT8_FILTER_MIN_SNR 6

#define FT8_SYMBOL_BT 2.0f ///< symbol smoothing filter bandwidth factor (BT)
#define FT4_SYMBOL_BT 1.0f ///< symbol smoothing filter bandwidth factor (BT)
#define FT2_SYMBOL_BT 1.0f ///< symbol smoothing filter bandwidth factor (BT)

#define FT8_TX_LEAD_SILENCE_SEC 0.5f
#define FT8_TX_TAIL_SILENCE_SEC 0.1f
#define FT4_TX_LEAD_SILENCE_SEC 0.3f
#define FT4_TX_TAIL_SILENCE_SEC 0.1f
#define FT2_TX_LEAD_SILENCE_SEC 0.1f
#define FT2_TX_TAIL_SILENCE_SEC 0.1f
#define FT2_RX_CAPTURE_MARGIN_SEC 0.1f
#define FT2_TX_LATE_GRACE_SEC 0.5f

typedef enum
{
	GCFT8_MODE_FT8,
	GCFT8_MODE_FT4,
	GCFT8_MODE_FT2
} gcft8_mode_t;

typedef struct
{
	const char* band;
	int ft8_frequency_hz;
	int ft4_frequency_hz;
	int ft2_frequency_hz;
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
	float tx_late_grace;
} gcft8_mode_config_t;

typedef struct
{
	int lon_min;
	int lon_max;
	int lat_min;
	int lat_max;
} gcft8_locator_zone_t;

typedef struct
{
	bool is_cq;
	bool is_for_local;
	bool has_locator;
	char cq_tag[FT8_TOKEN_TEXT_SIZE];
	char to_call[FT8_TOKEN_TEXT_SIZE];
	char from_call[FT8_TOKEN_TEXT_SIZE];
	char locator[FT8_TOKEN_TEXT_SIZE];
	char message[FT8_TOKEN_TEXT_SIZE];
} gcft8_decoded_message_view_t;

typedef enum
{
	GCFT8_REJECT_NONE,
	GCFT8_REJECT_MISSING_CALLSIGN,
	GCFT8_REJECT_MISSING_LOCATOR,
	GCFT8_REJECT_ALREADY_WORKED,
	GCFT8_REJECT_SNR,
	GCFT8_REJECT_PREFIX,
	GCFT8_REJECT_LOCATOR_ZONE,
	GCFT8_REJECT_SP_TAG
} gcft8_reject_reason_t;

typedef enum
{
	GCFT8_DISPLAY_NORMAL,
	GCFT8_DISPLAY_FOR_LOCAL,
	GCFT8_DISPLAY_CQ_CANDIDATE,
	GCFT8_DISPLAY_ALREADY_WORKED,
	GCFT8_DISPLAY_FILTERED
} gcft8_display_class_t;

typedef struct
{
	ftx_decoded_message_t decoded;
	gcft8_decoded_message_view_t view;
	float frequency_hz;
	float time_sec;
	int snr;
	int score;
	gcft8_reject_reason_t reject_reason;
	gcft8_display_class_t display_class;
} gcft8_rx_candidate_t;

typedef enum
{
	GCFT8_QSO_ACTION_NONE,
	GCFT8_QSO_ACTION_TX,
	GCFT8_QSO_ACTION_LOG_NOW,
	GCFT8_QSO_ACTION_TX_AND_LOG_AFTER
} gcft8_qso_action_type_t;

typedef struct
{
	gcft8_qso_action_type_t type;
	int tx_seq;
} gcft8_qso_action_t;

typedef struct
{
	char adif_record[GCFT8_ADIF_RECORD_TEXT_SIZE];
	char callsign[20];
	int session_index;
} gcft8_pending_log_t;

typedef struct
{
	float* signal;
	size_t signal_capacity;
	int16_t* pcm;
	size_t pcm_capacity;
	gfsk_scratch_t gfsk_scratch;
} gcft8_tx_context_t;

typedef struct
{
	monitor_t monitor;
	monitor_config_t config;
	bool monitor_initialized;
	float* signal;
	size_t signal_capacity;
	char* raw_data;
	size_t raw_data_capacity;
	ftx_candidate_t* candidate_list;
	size_t candidate_capacity;
	gcft8_rx_candidate_t* cq_candidates;
	size_t cq_candidate_capacity;
	ftx_message_t* decoded;
	size_t decoded_capacity;
	ftx_message_t** decoded_hashtable;
	size_t decoded_hashtable_capacity;
} gcft8_rx_context_t;

static const gcft8_band_frequency_t gcft8_band_frequencies[] = {
	{ "80",  3573000,  3575000,  3578000 },
	{ "60",  5357000,  5357000,  5360000 },
	{ "40",  7074000,  7047500,  7062000 },
	{ "30", 10136000, 10140000, 10144000 },
	{ "20", 14074000, 14080000, 14084000 },
	{ "17", 18100000, 18104000, 18108000 },
	{ "15", 21074000, 21140000, 21144000 },
	{ "12", 24915000, 24919000, 24923000 },
	{ "11", 27245000,        0,        0 },
	{ "10", 28074000, 28180000, 28184000 }
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
		.rx_capture_time = 13.6f,
		.tx_late_grace = 0.0f
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
		.rx_capture_time = FT4_SLOT_TIME - 0.4f,
		.tx_late_grace = 0.0f
	},
	{
		.mode = GCFT8_MODE_FT2,
		.name = "ft2",
		.protocol = FTX_PROTOCOL_FT2,
		.num_tones = FT2_NN,
		.symbol_period = FT2_SYMBOL_PERIOD,
		.symbol_bt = FT2_SYMBOL_BT,
		.slot_time = FT2_SLOT_TIME,
		.tx_lead_silence = FT2_TX_LEAD_SILENCE_SEC,
		.tx_tail_silence = FT2_TX_TAIL_SILENCE_SEC,
		.rx_capture_time = FT2_SLOT_TIME - FT2_RX_CAPTURE_MARGIN_SEC,
		.tx_late_grace = FT2_TX_LATE_GRACE_SEC
	}
};

typedef struct
{
	char callsign[12];
	uint32_t hash;
} callsign_hash_entry_t;

static callsign_hash_entry_t callsign_hash_cache[CALLSIGN_HASH_CACHE_SIZE];
static int callsign_hash_cache_size;
static pthread_mutex_t callsign_hash_cache_lock = PTHREAD_MUTEX_INITIALIZER;
static volatile sig_atomic_t shutdown_requested;
static gcft8_mode_t gcft8_mode = GCFT8_MODE_FT8;
static bool gcft8_filter_has_band;
static char gcft8_filter_band[8];
static int gcft8_filter_frequency_mhz;
static bool gcft8_snr_min_enabled;
static int gcft8_snr_min;
static char gcft8_only_prefixes[GCFT8_MAX_ONLY_PREFIXES][GCFT8_PREFIX_TEXT_SIZE];
static size_t gcft8_only_prefix_count;
static char gcft8_only_sp_tags[GCFT8_MAX_SP_TAGS][GCFT8_SP_TAG_TEXT_SIZE];
static size_t gcft8_only_sp_tag_count;
static gcft8_locator_zone_t gcft8_locator_zones[GCFT8_MAX_LOCATOR_ZONES];
static size_t gcft8_locator_zone_count;
static HashTable* ht_callsigntable_for_filter;
static gcft8_tx_context_t gcft8_tx_context;
static gcft8_rx_context_t gcft8_rx_context;

static bool gcft8_snr_filter_rejects(int snr);
static bool gcft8_prefix_filter_rejects(const char* callsign);
static bool gcft8_locator_zone_filter_rejects(const char* locator);
static void printDateTime_log(void);
static void tranceiver_rtx(gcft8_trx_state_t ptt);
static void gcft8_qso_sessions_init(void);
static void gcft8_qso_prune_expired(time_t now);
static bool gcft8_qso_update_from_direct_candidate(const gcft8_rx_candidate_t* candidate, time_t now, bool* tx_needed);
static int gcft8_qso_update_from_cq_candidate(const gcft8_rx_candidate_t* candidate, time_t now);
static int gcft8_qso_select_tx_session(void);
static bool gcft8_qso_build_tx_message(const gcft8_qso_session_t* session, char* dst, size_t dst_size);
static void gcft8_qso_mark_tx_sent(int session_idx, int sent_seq, bool tx_ok, time_t now);
static bool gcft8_qso_take_pending_log(gcft8_pending_log_t* pending_log);
static void gcft8_qso_finish_pending_log(const gcft8_pending_log_t* pending_log, bool success);
static bool log_adif_qso(const char* adif_record);
static void log_qso_to_filter_table(const char* callsign);
static void gcft8_flush_pending_logs(void);

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

static bool gcft8_starts_with(const char* text, const char* prefix)
{
	if ((text == NULL) || (prefix == NULL))
		return false;

	return strncmp(text, prefix, strlen(prefix)) == 0;
}

static void gcft8_decoded_message_view_init(gcft8_decoded_message_view_t* view)
{
	memset(view, 0, sizeof(*view));
}

static bool gcft8_is_signed_number_text(const char* text)
{
	if ((text == NULL) || ((text[0] != '-') && (text[0] != '+')) || !isdigit((unsigned char)text[1]))
		return false;

	for (size_t idx = 2; text[idx] != '\0'; ++idx)
	{
		if (!isdigit((unsigned char)text[idx]))
			return false;
	}

	return true;
}

static bool gcft8_field_is_locator(const ftx_decoded_field_t* field)
{
	if ((field == NULL) || (field->type != FTX_FIELD_GRID) || (strlen(field->text) != 4))
		return false;

	return true;
}

static void gcft8_analyze_decoded_message(const ftx_decoded_message_t* decoded, const char* local_callsign, gcft8_decoded_message_view_t* view)
{
	const ftx_decoded_field_t* field0;
	const ftx_decoded_field_t* field1;
	const ftx_decoded_field_t* field2;

	if (view == NULL)
		return;

	gcft8_decoded_message_view_init(view);
	if ((decoded == NULL) || (decoded->field_count == 0))
		return;

	field0 = &decoded->fields[0];
	field1 = (decoded->field_count > 1) ? &decoded->fields[1] : NULL;
	field2 = (decoded->field_count > 2) ? &decoded->fields[2] : NULL;

	if (((field0->type == FTX_FIELD_TOKEN) || (field0->type == FTX_FIELD_TOKEN_WITH_ARG)) &&
		((strcmp(field0->text, "CQ") == 0) || gcft8_starts_with(field0->text, "CQ ")))
	{
		view->is_cq = true;
		copy_text(view->to_call, sizeof(view->to_call), field0->text);
		if (gcft8_starts_with(field0->text, "CQ "))
			copy_text(view->cq_tag, sizeof(view->cq_tag), field0->text + 3);

		if ((field1 != NULL) && (field1->type == FTX_FIELD_CALL))
			copy_text(view->from_call, sizeof(view->from_call), field1->text);

		if (gcft8_field_is_locator(field2))
		{
			copy_text(view->locator, sizeof(view->locator), field2->text);
			view->has_locator = true;
		}
		return;
	}

	if ((decoded->field_count >= 2) && (field0->type == FTX_FIELD_CALL) && (field1 != NULL) && (field1->type == FTX_FIELD_CALL))
	{
		copy_text(view->to_call, sizeof(view->to_call), field0->text);
		copy_text(view->from_call, sizeof(view->from_call), field1->text);

		if (field2 != NULL)
		{
			copy_text(view->message, sizeof(view->message), field2->text);
			if (gcft8_field_is_locator(field2))
			{
				copy_text(view->locator, sizeof(view->locator), field2->text);
				view->has_locator = true;
			}
		}

		view->is_for_local = (local_callsign != NULL) && (strcmp(view->to_call, local_callsign) == 0);
	}
}

static void gcft8_rx_candidate_init(gcft8_rx_candidate_t* candidate)
{
	if (candidate == NULL)
		return;

	memset(candidate, 0, sizeof(*candidate));
	candidate->snr = -99;
	candidate->reject_reason = GCFT8_REJECT_NONE;
	candidate->display_class = GCFT8_DISPLAY_NORMAL;
}

static bool gcft8_find_decoded_slot(ftx_message_t* const decoded_hashtable[], int table_size, const ftx_message_t* message, int* slot_idx, bool* found_duplicate)
{
	int idx_hash;
	int probes = 0;

	if ((decoded_hashtable == NULL) || (table_size <= 0) || (message == NULL) || (slot_idx == NULL) || (found_duplicate == NULL))
		return false;

	idx_hash = message->hash % table_size;
	*found_duplicate = false;

	do
	{
		if (decoded_hashtable[idx_hash] == NULL)
		{
			*slot_idx = idx_hash;
			return true;
		}

		if ((decoded_hashtable[idx_hash]->hash == message->hash) && (memcmp(decoded_hashtable[idx_hash]->payload, message->payload, sizeof(message->payload)) == 0))
		{
			*slot_idx = idx_hash;
			*found_duplicate = true;
			return true;
		}

		idx_hash = (idx_hash + 1) % table_size;
		++probes;
	} while (probes < table_size);

	return false;
}

static bool gcft8_sp_tag_filter_rejects(const char* tag)
{
	if (gcft8_only_sp_tag_count == 0)
		return false;

	if ((tag == NULL) || (tag[0] == '\0'))
		return true;

	for (size_t idx = 0; idx < gcft8_only_sp_tag_count; ++idx)
	{
		if (strcmp(tag, gcft8_only_sp_tags[idx]) == 0)
			return false;
	}

	return true;
}

static gcft8_reject_reason_t gcft8_cq_reject_reason(const gcft8_rx_candidate_t* candidate)
{
	if (candidate == NULL)
		return GCFT8_REJECT_MISSING_CALLSIGN;

	if (candidate->view.from_call[0] == '\0')
		return GCFT8_REJECT_MISSING_CALLSIGN;

	if (ht_check(ht_callsigntable_for_filter, candidate->view.from_call))
		return GCFT8_REJECT_ALREADY_WORKED;

	if (!candidate->view.has_locator)
		return GCFT8_REJECT_MISSING_LOCATOR;

	if (gcft8_snr_filter_rejects(candidate->snr))
		return GCFT8_REJECT_SNR;

	if (gcft8_prefix_filter_rejects(candidate->view.from_call))
		return GCFT8_REJECT_PREFIX;

	if (gcft8_locator_zone_filter_rejects(candidate->view.locator))
		return GCFT8_REJECT_LOCATOR_ZONE;

	if (gcft8_sp_tag_filter_rejects(candidate->view.cq_tag))
		return GCFT8_REJECT_SP_TAG;

	return GCFT8_REJECT_NONE;
}

static gcft8_display_class_t gcft8_display_class_for_candidate(const gcft8_rx_candidate_t* candidate)
{
	if (candidate == NULL)
		return GCFT8_DISPLAY_NORMAL;

	if (candidate->view.is_for_local)
		return GCFT8_DISPLAY_FOR_LOCAL;

	if (candidate->view.is_cq)
	{
		if (candidate->reject_reason == GCFT8_REJECT_ALREADY_WORKED)
			return GCFT8_DISPLAY_ALREADY_WORKED;

		if (candidate->reject_reason != GCFT8_REJECT_NONE)
			return GCFT8_DISPLAY_FILTERED;

		return GCFT8_DISPLAY_CQ_CANDIDATE;
	}

	return GCFT8_DISPLAY_NORMAL;
}

static void gcft8_display_rx_candidate(const gcft8_rx_candidate_t* candidate)
{
	const char* color_start = "";
	const char* color_end = "";

	if (candidate == NULL)
		return;

	switch (candidate->display_class)
	{
	case GCFT8_DISPLAY_FOR_LOCAL:
		color_start = "\033[1;31m";
		color_end = "\033[0m";
		break;
	case GCFT8_DISPLAY_CQ_CANDIDATE:
		color_start = "\033[1;34m";
		color_end = "\033[0m";
		break;
	case GCFT8_DISPLAY_ALREADY_WORKED:
		color_start = "\033[1;33m";
		color_end = "\033[0m";
		break;
	case GCFT8_DISPLAY_FILTERED:
		color_start = "\033[1;35m";
		color_end = "\033[0m";
		break;
	case GCFT8_DISPLAY_NORMAL:
	default:
		break;
	}

	printDateTime_log();
	printf(" %d %3d %+4.2f %4.0f ~  %s%s%s\n", candidate->snr, candidate->score, candidate->time_sec, candidate->frequency_hz, color_start, candidate->decoded.text, color_end);
}

static bool gcft8_ensure_float_buffer(float** buffer, size_t* capacity, size_t required)
{
	float* resized;

	if ((buffer == NULL) || (capacity == NULL))
		return false;

	if (required <= *capacity)
		return true;

	resized = (float*)realloc(*buffer, required * sizeof((*buffer)[0]));
	if (resized == NULL)
		return false;

	*buffer = resized;
	*capacity = required;
	return true;
}

static bool gcft8_ensure_i16_buffer(int16_t** buffer, size_t* capacity, size_t required)
{
	int16_t* resized;

	if ((buffer == NULL) || (capacity == NULL))
		return false;

	if (required <= *capacity)
		return true;

	resized = (int16_t*)realloc(*buffer, required * sizeof((*buffer)[0]));
	if (resized == NULL)
		return false;

	*buffer = resized;
	*capacity = required;
	return true;
}

static bool gcft8_ensure_buffer(void** buffer, size_t* capacity, size_t required, size_t element_size)
{
	void* resized;

	if ((buffer == NULL) || (capacity == NULL) || (element_size == 0))
		return false;

	if (required <= *capacity)
		return true;

	if (required > (SIZE_MAX / element_size))
		return false;

	resized = realloc(*buffer, required * element_size);
	if (resized == NULL)
		return false;

	*buffer = resized;
	*capacity = required;
	return true;
}

static void gcft8_tx_context_init(gcft8_tx_context_t* ctx)
{
	if (ctx == NULL)
		return;

	memset(ctx, 0, sizeof(*ctx));
	gfsk_scratch_init(&ctx->gfsk_scratch);
}

static bool gcft8_tx_context_ensure(gcft8_tx_context_t* ctx, size_t sample_count)
{
	if (ctx == NULL)
		return false;

	return gcft8_ensure_float_buffer(&ctx->signal, &ctx->signal_capacity, sample_count) &&
		gcft8_ensure_i16_buffer(&ctx->pcm, &ctx->pcm_capacity, sample_count);
}

static void gcft8_tx_context_free(gcft8_tx_context_t* ctx)
{
	if (ctx == NULL)
		return;

	free(ctx->signal);
	free(ctx->pcm);
	gfsk_scratch_free(&ctx->gfsk_scratch);
	gcft8_tx_context_init(ctx);
}

static void gcft8_rx_context_init(gcft8_rx_context_t* ctx)
{
	if (ctx == NULL)
		return;

	memset(ctx, 0, sizeof(*ctx));
}

static bool gcft8_monitor_config_equal(const monitor_config_t* a, const monitor_config_t* b)
{
	return (a != NULL) && (b != NULL) &&
		(a->f_min == b->f_min) &&
		(a->f_max == b->f_max) &&
		(a->sample_rate == b->sample_rate) &&
		(a->time_osr == b->time_osr) &&
		(a->freq_osr == b->freq_osr) &&
		(a->protocol == b->protocol);
}

static void gcft8_rx_context_free(gcft8_rx_context_t* ctx)
{
	if (ctx == NULL)
		return;

	if (ctx->monitor_initialized)
		monitor_free(&ctx->monitor);
	free(ctx->signal);
	free(ctx->raw_data);
	free(ctx->candidate_list);
	free(ctx->cq_candidates);
	free(ctx->decoded);
	free(ctx->decoded_hashtable);
	gcft8_rx_context_init(ctx);
}

static bool gcft8_rx_context_prepare(gcft8_rx_context_t* ctx, const monitor_config_t* config, int max_candidates, int max_decoded_messages)
{
	if ((ctx == NULL) || (config == NULL) || (max_candidates <= 0) || (max_decoded_messages <= 0))
		return false;

	if (ctx->monitor_initialized && !gcft8_monitor_config_equal(&ctx->config, config))
	{
		monitor_free(&ctx->monitor);
		ctx->monitor_initialized = false;
	}

	if (!ctx->monitor_initialized)
	{
		monitor_init(&ctx->monitor, config);
		ctx->config = *config;
		ctx->monitor_initialized = true;
	}

	return gcft8_ensure_float_buffer(&ctx->signal, &ctx->signal_capacity, (size_t)ctx->monitor.block_size) &&
		gcft8_ensure_buffer((void**)&ctx->raw_data, &ctx->raw_data_capacity, (size_t)ctx->monitor.block_size * 2u, sizeof(ctx->raw_data[0])) &&
		gcft8_ensure_buffer((void**)&ctx->candidate_list, &ctx->candidate_capacity, (size_t)max_candidates, sizeof(ctx->candidate_list[0])) &&
		gcft8_ensure_buffer((void**)&ctx->cq_candidates, &ctx->cq_candidate_capacity, (size_t)max_candidates, sizeof(ctx->cq_candidates[0])) &&
		gcft8_ensure_buffer((void**)&ctx->decoded, &ctx->decoded_capacity, (size_t)max_decoded_messages, sizeof(ctx->decoded[0])) &&
		gcft8_ensure_buffer((void**)&ctx->decoded_hashtable, &ctx->decoded_hashtable_capacity, (size_t)max_decoded_messages, sizeof(ctx->decoded_hashtable[0]));
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
			int frequency;

			switch (mode)
			{
				case GCFT8_MODE_FT2:
					frequency = gcft8_band_frequencies[idx].ft2_frequency_hz;
					break;
				case GCFT8_MODE_FT4:
					frequency = gcft8_band_frequencies[idx].ft4_frequency_hz;
					break;
				case GCFT8_MODE_FT8:
				default:
					frequency = gcft8_band_frequencies[idx].ft8_frequency_hz;
					break;
			}

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
	switch (gcft8_mode)
	{
		case GCFT8_MODE_FT2:
			return "FT2";
		case GCFT8_MODE_FT4:
			return "FT4";
		case GCFT8_MODE_FT8:
		default:
			return "FT8";
	}
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

static const char* ft8_filter_mode_name(int filter_mode)
{
	switch (filter_mode)
	{
	case FT8_FILTER_LISTEN_ONLY:
		return "Listen only, no automatic TX";
	case FT8_FILTER_RANDOM_CQ:
		return "Random CQ selection";
	case FT8_FILTER_BEST_DECODE_SCORE:
		return "Best decode score";
	case FT8_FILTER_MAX_DISTANCE:
		return "Maximum distance";
	case FT8_FILTER_MIN_DISTANCE:
		return "Minimum distance";
	case FT8_FILTER_MAX_SNR:
		return "Maximum SNR";
	case FT8_FILTER_MIN_SNR:
		return "Minimum SNR";
	default:
		return "Unknown";
	}
}

static bool gcft8_parse_snr_min(const char* value, int* snr_min)
{
	char* endptr = NULL;
	long result;

	if ((value == NULL) || (value[0] == '\0') || (snr_min == NULL))
		return false;

	errno = 0;
	result = strtol(value, &endptr, 10);
	if ((errno == ERANGE) || (endptr == value) || (endptr == NULL) || (*endptr != '\0'))
		return false;

	if ((result < INT_MIN) || (result > INT_MAX))
		return false;

	*snr_min = (int)result;
	return true;
}

static bool gcft8_snr_filter_rejects(int snr)
{
	return gcft8_snr_min_enabled && (snr < gcft8_snr_min);
}

static bool gcft8_parse_prefix_token(const char* start, size_t len, char* dst, size_t dst_size)
{
	size_t out_len = 0;

	while ((len > 0) && isspace((unsigned char)*start))
	{
		++start;
		--len;
	}

	while ((len > 0) && isspace((unsigned char)start[len - 1]))
		--len;

	if ((len == 0) || (len >= dst_size))
		return false;

	for (size_t idx = 0; idx < len; ++idx)
	{
		unsigned char ch = (unsigned char)start[idx];

		if (!isalnum(ch))
			return false;

		dst[out_len++] = (char)toupper(ch);
	}

	dst[out_len] = '\0';
	return true;
}

static bool gcft8_parse_only_prefixes(const char* value)
{
	char parsed[GCFT8_MAX_ONLY_PREFIXES][GCFT8_PREFIX_TEXT_SIZE];
	size_t parsed_count = 0;
	size_t token_start = 0;
	size_t idx = 0;

	if ((value == NULL) || (value[0] == '\0'))
		return false;

	memset(parsed, 0, sizeof(parsed));

	for (;;)
	{
		if ((value[idx] == ',') || (value[idx] == '\0'))
		{
			if (parsed_count >= GCFT8_MAX_ONLY_PREFIXES)
				return false;

			if (!gcft8_parse_prefix_token(value + token_start, idx - token_start, parsed[parsed_count], sizeof(parsed[parsed_count])))
				return false;

			++parsed_count;

			if (value[idx] == '\0')
				break;

			token_start = idx + 1;
		}

		++idx;
	}

	memcpy(gcft8_only_prefixes, parsed, sizeof(parsed));
	gcft8_only_prefix_count = parsed_count;
	return true;
}

static bool gcft8_parse_sp_tag_token(const char* start, size_t len, char* dst, size_t dst_size)
{
	size_t out_len = 0;
	int digits = 0;
	int letters = 0;

	while ((len > 0) && isspace((unsigned char)*start))
	{
		++start;
		--len;
	}

	while ((len > 0) && isspace((unsigned char)start[len - 1]))
		--len;

	if ((len == 0) || (len >= dst_size) || (len > 4))
		return false;

	for (size_t idx = 0; idx < len; ++idx)
	{
		unsigned char ch = (unsigned char)start[idx];

		if (isalpha(ch))
		{
			++letters;
			dst[out_len++] = (char)toupper(ch);
		}
		else if (isdigit(ch))
		{
			++digits;
			dst[out_len++] = (char)ch;
		}
		else
		{
			return false;
		}
	}

	if (!(((letters >= 1) && (letters <= 4) && (digits == 0)) || ((digits == 3) && (letters == 0))))
		return false;

	dst[out_len] = '\0';
	return true;
}

static bool gcft8_parse_only_sp_tags(const char* value)
{
	char parsed[GCFT8_MAX_SP_TAGS][GCFT8_SP_TAG_TEXT_SIZE];
	size_t parsed_count = 0;
	size_t token_start = 0;
	size_t idx = 0;

	if ((value == NULL) || (value[0] == '\0'))
		return false;

	memset(parsed, 0, sizeof(parsed));

	for (;;)
	{
		if ((value[idx] == ',') || (value[idx] == '\0'))
		{
			if (parsed_count >= GCFT8_MAX_SP_TAGS)
				return false;

			if (!gcft8_parse_sp_tag_token(value + token_start, idx - token_start, parsed[parsed_count], sizeof(parsed[parsed_count])))
				return false;

			++parsed_count;

			if (value[idx] == '\0')
				break;

			token_start = idx + 1;
		}

		++idx;
	}

	memset(gcft8_only_sp_tags, 0, sizeof(gcft8_only_sp_tags));
	memcpy(gcft8_only_sp_tags, parsed, parsed_count * sizeof(parsed[0]));
	gcft8_only_sp_tag_count = parsed_count;
	return true;
}

static bool gcft8_is_simple_portable_suffix(const char* suffix)
{
	return (strcmp(suffix, "P") == 0) ||
	       (strcmp(suffix, "M") == 0) ||
	       (strcmp(suffix, "MM") == 0) ||
	       (strcmp(suffix, "AM") == 0) ||
	       (strcmp(suffix, "QRP") == 0);
}

static void gcft8_clean_callsign_for_prefix(const char* callsign, char* dst, size_t dst_size)
{
	size_t out_len = 0;

	if (dst_size == 0)
		return;

	if (callsign != NULL)
	{
		for (size_t idx = 0; callsign[idx] != '\0' && out_len + 1 < dst_size; ++idx)
			dst[out_len++] = callsign[idx];
	}

	dst[out_len] = '\0';

	for (;;)
	{
		char* slash = strrchr(dst, '/');

		if ((slash == NULL) || !gcft8_is_simple_portable_suffix(slash + 1))
			break;

		*slash = '\0';
	}
}

static bool gcft8_prefix_filter_rejects(const char* callsign)
{
	char cleaned_callsign[FT8_TOKEN_TEXT_SIZE];

	if (gcft8_only_prefix_count == 0)
		return false;

	gcft8_clean_callsign_for_prefix(callsign, cleaned_callsign, sizeof(cleaned_callsign));
	if (cleaned_callsign[0] == '\0')
		return true;

	for (size_t idx = 0; idx < gcft8_only_prefix_count; ++idx)
	{
		size_t prefix_len = strlen(gcft8_only_prefixes[idx]);

		if (strncmp(cleaned_callsign, gcft8_only_prefixes[idx], prefix_len) == 0)
			return false;
	}

	return true;
}

static bool gcft8_parse_locator_zone_endpoint(const char* start, size_t len, int* lon, int* lat)
{
	char lon_char;
	char lat_char;

	while ((len > 0) && isspace((unsigned char)*start))
	{
		++start;
		--len;
	}

	while ((len > 0) && isspace((unsigned char)start[len - 1]))
		--len;

	if ((len != 2) || (lon == NULL) || (lat == NULL))
		return false;

	lon_char = (char)toupper((unsigned char)start[0]);
	lat_char = (char)toupper((unsigned char)start[1]);

	if ((lon_char < 'A') || (lon_char > 'R') || (lat_char < 'A') || (lat_char > 'R'))
		return false;

	*lon = lon_char - 'A';
	*lat = lat_char - 'A';
	return true;
}

static bool gcft8_parse_locator_zone_token(const char* start, size_t len, gcft8_locator_zone_t* zone)
{
	bool found_colon = false;
	size_t colon_pos = 0;
	int lon_a;
	int lat_a;
	int lon_b;
	int lat_b;

	if (zone == NULL)
		return false;

	for (size_t idx = 0; idx < len; ++idx)
	{
		if (start[idx] == ':')
		{
			if (found_colon)
				return false;

			found_colon = true;
			colon_pos = idx;
		}
	}

	if (!found_colon)
		return false;

	if (!gcft8_parse_locator_zone_endpoint(start, colon_pos, &lon_a, &lat_a))
		return false;

	if (!gcft8_parse_locator_zone_endpoint(start + colon_pos + 1, len - colon_pos - 1, &lon_b, &lat_b))
		return false;

	zone->lon_min = (lon_a < lon_b) ? lon_a : lon_b;
	zone->lon_max = (lon_a > lon_b) ? lon_a : lon_b;
	zone->lat_min = (lat_a < lat_b) ? lat_a : lat_b;
	zone->lat_max = (lat_a > lat_b) ? lat_a : lat_b;
	return true;
}

static bool gcft8_parse_locator_zones(const char* value)
{
	gcft8_locator_zone_t parsed[GCFT8_MAX_LOCATOR_ZONES];
	size_t parsed_count = 0;
	size_t token_start = 0;
	size_t idx = 0;

	if ((value == NULL) || (value[0] == '\0'))
		return false;

	memset(parsed, 0, sizeof(parsed));

	for (;;)
	{
		if ((value[idx] == ',') || (value[idx] == '\0'))
		{
			if (parsed_count >= GCFT8_MAX_LOCATOR_ZONES)
				return false;

			if (!gcft8_parse_locator_zone_token(value + token_start, idx - token_start, &parsed[parsed_count]))
				return false;

			++parsed_count;

			if (value[idx] == '\0')
				break;

			token_start = idx + 1;
		}

		++idx;
	}

	memset(gcft8_locator_zones, 0, sizeof(gcft8_locator_zones));
	memcpy(gcft8_locator_zones, parsed, parsed_count * sizeof(parsed[0]));
	gcft8_locator_zone_count = parsed_count;
	return true;
}

static bool gcft8_locator_zone_filter_rejects(const char* locator)
{
	int lon;
	int lat;

	if (gcft8_locator_zone_count == 0)
		return false;

	if ((locator == NULL) || (locator[0] < 'A') || (locator[0] > 'R') || (locator[1] < 'A') || (locator[1] > 'R'))
		return true;

	lon = locator[0] - 'A';
	lat = locator[1] - 'A';

	for (size_t idx = 0; idx < gcft8_locator_zone_count; ++idx)
	{
		const gcft8_locator_zone_t* zone = &gcft8_locator_zones[idx];

		if ((lon >= zone->lon_min) && (lon <= zone->lon_max) && (lat >= zone->lat_min) && (lat <= zone->lat_max))
			return false;
	}

	return true;
}

static void gcft8_format_only_prefixes(char* dst, size_t dst_size)
{
	size_t offset = 0;

	if (dst_size == 0)
		return;

	if (gcft8_only_prefix_count == 0)
	{
		copy_text(dst, dst_size, "off");
		return;
	}

	dst[0] = '\0';
	for (size_t idx = 0; idx < gcft8_only_prefix_count; ++idx)
	{
		int written = snprintf(dst + offset, dst_size - offset, "%s%s", (idx == 0) ? "" : ",", gcft8_only_prefixes[idx]);

		if (written < 0)
			break;

		if ((size_t)written >= dst_size - offset)
		{
			dst[dst_size - 1] = '\0';
			break;
		}

		offset += (size_t)written;
	}
}

static void gcft8_format_only_sp_tags(char* dst, size_t dst_size)
{
	size_t offset = 0;

	if (dst_size == 0)
		return;

	if (gcft8_only_sp_tag_count == 0)
	{
		copy_text(dst, dst_size, "off");
		return;
	}

	dst[0] = '\0';
	for (size_t idx = 0; idx < gcft8_only_sp_tag_count; ++idx)
	{
		if (idx > 0)
		{
			if ((offset + 1) >= dst_size)
				break;

			dst[offset++] = ',';
			dst[offset] = '\0';
		}

		for (size_t src_idx = 0; gcft8_only_sp_tags[idx][src_idx] != '\0'; ++src_idx)
		{
			if ((offset + 1) >= dst_size)
				break;

			dst[offset++] = gcft8_only_sp_tags[idx][src_idx];
			dst[offset] = '\0';
		}
	}
}

static void gcft8_format_locator_zones(char* dst, size_t dst_size)
{
	size_t offset = 0;

	if (dst_size == 0)
		return;

	if (gcft8_locator_zone_count == 0)
	{
		copy_text(dst, dst_size, "off");
		return;
	}

	dst[0] = '\0';
	for (size_t idx = 0; idx < gcft8_locator_zone_count; ++idx)
	{
		const gcft8_locator_zone_t* zone = &gcft8_locator_zones[idx];
		int written = snprintf(dst + offset, dst_size - offset, "%s%c%c:%c%c",
			(idx == 0) ? "" : ",",
			(char)('A' + zone->lon_min),
			(char)('A' + zone->lat_max),
			(char)('A' + zone->lon_max),
			(char)('A' + zone->lat_min));

		if (written < 0)
			break;

		if ((size_t)written >= dst_size - offset)
		{
			dst[dst_size - 1] = '\0';
			break;
		}

		offset += (size_t)written;
	}
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

static void callsign_hash_cache_init_unlocked(void)
{
	callsign_hash_cache_size = 0;
	memset(callsign_hash_cache, 0, sizeof(callsign_hash_cache));
}

static void callsign_hash_cache_init(void)
{
	pthread_mutex_lock(&callsign_hash_cache_lock);
	callsign_hash_cache_init_unlocked();
	pthread_mutex_unlock(&callsign_hash_cache_lock);
}

static void callsign_hash_cache_insert_unlocked(const char* callsign, uint32_t stored_hash)
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

	pthread_mutex_lock(&callsign_hash_cache_lock);
	memcpy(old_cache, callsign_hash_cache, sizeof(old_cache));
	callsign_hash_cache_init_unlocked();

	for (int idx = 0; idx < CALLSIGN_HASH_CACHE_SIZE; ++idx)
	{
		if (old_cache[idx].callsign[0] != '\0')
		{
			uint8_t age = (uint8_t)(old_cache[idx].hash >> 24);
			if (age <= max_age)
			{
				uint32_t stored_hash = (((uint32_t)age + 1u) << 24) | (old_cache[idx].hash & 0x3FFFFFu);
				callsign_hash_cache_insert_unlocked(old_cache[idx].callsign, stored_hash);
			}
		}
	}
	pthread_mutex_unlock(&callsign_hash_cache_lock);
}

static void callsign_hash_cache_save(const char* callsign, uint32_t hash)
{
	pthread_mutex_lock(&callsign_hash_cache_lock);
	callsign_hash_cache_insert_unlocked(callsign, hash & 0x3FFFFFu);
	pthread_mutex_unlock(&callsign_hash_cache_lock);
}

static bool callsign_hash_cache_lookup(ftx_callsign_hash_type_t hash_type, uint32_t hash, char* callsign)
{
	uint8_t hash_shift = (hash_type == FTX_CALLSIGN_HASH_10_BITS) ? 12 : (hash_type == FTX_CALLSIGN_HASH_12_BITS ? 10 : 0);
	uint16_t hash10 = (hash >> (12 - hash_shift)) & 0x3FFu;
	int idx_hash = (hash10 * 23) % CALLSIGN_HASH_CACHE_SIZE;
	bool found = false;

	pthread_mutex_lock(&callsign_hash_cache_lock);

	for (int probes = 0; probes < CALLSIGN_HASH_CACHE_SIZE; ++probes)
	{
		if (callsign_hash_cache[idx_hash].callsign[0] == '\0')
			break;

		if (((callsign_hash_cache[idx_hash].hash & 0x3FFFFFu) >> hash_shift) == hash)
		{
			copy_text(callsign, 12, callsign_hash_cache[idx_hash].callsign);
			found = true;
			break;
		}

		idx_hash = (idx_hash + 1) % CALLSIGN_HASH_CACHE_SIZE;
	}

	if (!found)
		callsign[0] = '\0';

	pthread_mutex_unlock(&callsign_hash_cache_lock);
	return found;
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
	uint8_t tones[FT2_NN];
	int num_symbols;
	int num_tones;
	float signal_power = 0.0f;
	float noise_power = 0.0f;
	int signal_count = 0;
	int noise_count = 0;

	if (wf->protocol == FTX_PROTOCOL_FT2)
	{
		num_symbols = FT2_NN;
		num_tones = 4;
		ft2_encode(message->payload, tones);
	}
	else if (wf->protocol == FTX_PROTOCOL_FT4)
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
	
	.TRX_status = GCFT8_TRX_RX,
	.TRX_status_lock = PTHREAD_MUTEX_INITIALIZER,
	.RX_status_cond = PTHREAD_COND_INITIALIZER,
	.TX_status_cond = PTHREAD_COND_INITIALIZER,
	
	.Tranceiver_VFOA_Freq = 14074000,
	
	.log_file_name = "QSO_F4JJJ.adif",
	.beep_on_log=0,
	
	.filter_on_cq = 0
	
};

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
static double now()
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

static void advance_cursor(float slot_time) {
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

static bool wait_for_slot_start(float slot_time){
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

static bool wait_for_tx_slot_start(float slot_time, float late_grace)
{
	double slot = slot_time;
	double grace = late_grace;
	double current_time;
	double slot_start;
	double since_slot_start;

	if (slot <= 0.0)
		return !gcft8_shutdown_requested();

	current_time = now();
	slot_start = floor(current_time / slot) * slot;
	since_slot_start = current_time - slot_start;

	if ((grace > 0.0) && (since_slot_start >= 0.0) && (since_slot_start <= grace))
		return !gcft8_shutdown_requested();

	return wait_for_slot_start(slot_time);
}

static void printDateTime_log(){
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

static void print_slot_separator(const gcft8_mode_config_t* mode_cfg, const char* direction)
{
	long long slot_ms;
	long long now_ms;
	long long slot_start_ms;
	time_t t;
	int ms;
	struct tm tm;

	if (mode_cfg == NULL)
		mode_cfg = gcft8_current_mode_config();
	if (direction == NULL)
		direction = "";

	clear_status_line();
	slot_ms = (long long)(mode_cfg->slot_time * 1000.0f + 0.5f);
	now_ms = (long long)(now() * 1000.0);
	if (slot_ms <= 0)
		slot_start_ms = now_ms;
	else
		slot_start_ms = (now_ms / slot_ms) * slot_ms;

	t = (time_t)(slot_start_ms / 1000);
	ms = (int)(slot_start_ms % 1000);
	tm = *gmtime(&t);

	printf("\033[1;36m-- %s %s ", mode_cfg->name, direction);
	if (ms == 0)
		printf("%d-%02d-%02d %02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	else
		printf("%d-%02d-%02d %02d:%02d:%02d.%03d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, ms);
	printf(" --\033[0m\n");
}

static void printDateTime_ms(){
	clear_status_line();
	struct timeval tv;
	gettimeofday(&tv, 0);
	time_t t = tv.tv_sec;
	struct tm tm = *gmtime(&t);
	printf("%d-%02d-%02d %02d:%02d:%02d.%03ld", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, (long)(tv.tv_usec / 1000));
}

static void latLonForGrid(char * grid, float * latlon) {
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


static float latLonDist(float * latlonA, float * latlonB){
	double value = sin(latlonA[0] * M_PI / 180.0f) * sin(latlonB[0] * M_PI / 180.0f) + cos(latlonA[0] * M_PI / 180.0f) * cos(latlonB[0] * M_PI / 180.0f) * cos((latlonA[1] * M_PI / 180.0f) - (latlonB[1] * M_PI / 180.0f));
	if (value > 1.0)
		value = 1.0;
	else if (value < -1.0)
		value = -1.0;

	return (float)(acos(value)*6371.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Audio



/* Open and init the default recording device. */
static void capture_audioInit(void)
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

static void capture_audioDeInit(void)
{ 
	if (sound.capture_handle != NULL)
	{
		snd_pcm_drop(sound.capture_handle);
		snd_pcm_close(sound.capture_handle);
		sound.capture_handle = NULL;
	}
}

/* Open and init the default playback device. */
static void playback_audioInit(void)
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

static void playback_audioDeInit(void)
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
//Thread FT8
static float getFrame(char *buffer, int i)
{
	uint16_t lo = (uint8_t)buffer[2 * i];
	uint16_t hi = (uint8_t)buffer[(2 * i) + 1];
	int16_t sample = (int16_t)(lo | (hi << 8));
	return sample / 32768.0f;
}

static void unlock_TX_thread(){
	#if DEBUG
	printf( "UnLock TX thread\n");
	#endif
	pthread_mutex_lock(&FT8.TRX_status_lock);
	FT8.TRX_status = GCFT8_TRX_TX;
	pthread_cond_signal(&FT8.TX_status_cond);
	pthread_mutex_unlock(&FT8.TRX_status_lock);
}

static void unlock_RX_thread(){
	#if DEBUG
	printf( "UnLock RX thread\n");
	#endif
	pthread_mutex_lock(&FT8.TRX_status_lock);
	FT8.TRX_status = GCFT8_TRX_RX;
	pthread_cond_signal(&FT8.RX_status_cond);
	pthread_mutex_unlock(&FT8.TRX_status_lock);
}



static void RX_FT8()
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
	
	if (!gcft8_rx_context_prepare(&gcft8_rx_context, &mon_cfg, kMax_candidates, kMax_decoded_messages))
	{
		fprintf(stderr, "Out of memory while preparing RX context\n");
		exit(1);
	}

	snd_pcm_sframes_t rc;
	
	while(!gcft8_shutdown_requested()){
		pthread_mutex_lock(&FT8.TRX_status_lock);
		gcft8_trx_state_t stat = FT8.TRX_status;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if(stat == GCFT8_TRX_RX){
			monitor_t* mon = &gcft8_rx_context.monitor;
			float* signal = gcft8_rx_context.signal;
			char* raw_data = gcft8_rx_context.raw_data;
			ftx_candidate_t* candidate_list = gcft8_rx_context.candidate_list;
			gcft8_rx_candidate_t* cq_candidates = gcft8_rx_context.cq_candidates;
			ftx_message_t* decoded = gcft8_rx_context.decoded;
			ftx_message_t** decoded_hashtable = gcft8_rx_context.decoded_hashtable;

			monitor_reset(mon);
			
			#if DEBUG
			printf( "Waterfall allocated %d symbols\n", mon->wf.max_blocks);
			#endif

			if (!wait_for_slot_start(mode_cfg->slot_time))
			{
				break;
			}
			
			snd_pcm_reset(sound.capture_handle);
			
			for (int frame_pos = 0; !gcft8_shutdown_requested() && frame_pos + mon->block_size <= num_samples; frame_pos += mon->block_size)
			{
				// Process the waveform data frame by frame - you could have a live loop here with data from an audio device

				int frames_read = 0;
				while (!gcft8_shutdown_requested() && frames_read < mon->block_size)
				{
					rc = snd_pcm_readi(sound.capture_handle, raw_data + (frames_read * 2), (snd_pcm_uframes_t)(mon->block_size - frames_read));
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

				if (frames_read == mon->block_size)
				{
					for (int i = 0; i < mon->block_size; i++)
					{
						signal[i] = getFrame(raw_data,i);
					}
					
					monitor_process(mon, signal);
				}
				
			}
			
			#if DEBUG
			printf( "Waterfall accumulated %d symbols\n", mon->wf.num_blocks);
			printf( "Max magnitude: %.1f dB\n", mon->max_mag);
			#endif
			
			print_slot_separator(mode_cfg, "RX");

			// Find top candidates by Costas sync score and localize them in time and frequency
			int num_candidates = ftx_find_candidates(&mon->wf, kMax_candidates, candidate_list, kMin_score);

			int cq_candidate_count = 0;
			bool rx_requested_tx = false;
			bool should_unlock_tx = false;
			time_t rx_now = time(NULL);

			pthread_mutex_lock(&FT8.TRX_status_lock);
			gcft8_qso_prune_expired(rx_now);
			pthread_mutex_unlock(&FT8.TRX_status_lock);

			// Hash table for decoded messages (to check for duplicates)
			int num_decoded = 0;
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

				float freq_hz = (mon->min_bin + cand->freq_offset + (float)cand->freq_sub / mon->wf.freq_osr) / mon->symbol_period;
				float time_sec = (cand->time_offset + (float)cand->time_sub / mon->wf.time_osr) * mon->symbol_period;

				ftx_message_t message;
				ftx_decode_status_t status;
				if (!ftx_decode_candidate(&mon->wf, cand, kLDPC_iterations, &message, &status))
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

				#if DEBUG
				printf( "Checking hash table for %4.1fs / %4.1fHz [%d]...\n", time_sec, freq_hz, cand->score);
				#endif
				int idx_hash = 0;
				bool found_duplicate = false;
				if (!gcft8_find_decoded_slot(decoded_hashtable, kMax_decoded_messages, &message, &idx_hash, &found_duplicate))
				{
						#if DEBUG
						printf( "Decoded message table is full\n");
						#endif
					continue;
				}

				if (found_duplicate)
				{
					#if DEBUG
					printf( "Found a duplicate payload/hash\n");
					#endif
					continue;
				}

				gcft8_rx_candidate_t rx_candidate;
				gcft8_rx_candidate_init(&rx_candidate);
				ftx_message_rc_t message_status = ftx_message_decode(&message, &callsign_hash_if, &rx_candidate.decoded);
				if (message_status != FTX_MESSAGE_RC_OK)
				{
					#if DEBUG
					printf("Error while unpacking message: %d\n", (int)message_status);
					#endif
					continue;
				}

				memcpy(&decoded[idx_hash], &message, sizeof(message));
				decoded_hashtable[idx_hash] = &decoded[idx_hash];
				++num_decoded;

				rx_candidate.frequency_hz = freq_hz;
				rx_candidate.time_sec = time_sec;
				rx_candidate.score = cand->score;
				rx_candidate.snr = gcft8_estimate_snr(&mon->wf, cand, &message);
				gcft8_analyze_decoded_message(&rx_candidate.decoded, FT8.Local_CALLSIGN, &rx_candidate.view);

				if (rx_candidate.view.is_cq)
					rx_candidate.reject_reason = gcft8_cq_reject_reason(&rx_candidate);
				rx_candidate.display_class = gcft8_display_class_for_candidate(&rx_candidate);
				gcft8_display_rx_candidate(&rx_candidate);

				if (rx_candidate.view.is_for_local)
				{
					if (!ft8_listen_only_filter_active())
					{
						bool tx_needed = false;
						pthread_mutex_lock(&FT8.TRX_status_lock);
						if (gcft8_qso_update_from_direct_candidate(&rx_candidate, rx_now, &tx_needed))
						{
							if (tx_needed)
								rx_requested_tx = true;
						}
						pthread_mutex_unlock(&FT8.TRX_status_lock);
					}
				}
				else if (rx_candidate.view.is_cq && (rx_candidate.reject_reason == GCFT8_REJECT_NONE) && (cq_candidate_count < kMax_candidates))
				{
					cq_candidates[cq_candidate_count++] = rx_candidate;
				}
			}

			if((cq_candidate_count>0) && !rx_requested_tx && !ft8_listen_only_filter_active() && !gcft8_shutdown_requested()){
				int index_from_ope = 0;
				float dist = 0;
				int actusnr;
				
				switch (FT8.filter_on_cq)
				{
				
				case FT8_FILTER_RANDOM_CQ:
					index_from_ope = ft8_random_index(cq_candidate_count);
					break;
					
				case FT8_FILTER_BEST_DECODE_SCORE:
					index_from_ope = 0;
					break;				
				
				case FT8_FILTER_MAX_DISTANCE:
					dist = 0;
					for(int i = 0; i<cq_candidate_count;i++){
						if(strlen(cq_candidates[i].view.locator) == 4){
							float latlonlocal[2];
							latLonForGrid(cq_candidates[i].view.locator,latlonlocal);
							float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
							#if DEBUG
							printf("%s lat:%f lon:%f dist:%f\n",cq_candidates[i].view.from_call,latlonlocal[0],latlonlocal[1],new_dist);
							#endif
							if(dist < new_dist){index_from_ope=i;dist=new_dist;}
						}
					}				
					break;
					
				case FT8_FILTER_MIN_DISTANCE:
					dist = 6372;
					for(int i = 0; i<cq_candidate_count;i++){
						if(strlen(cq_candidates[i].view.locator) == 4){
							float latlonlocal[2];
							latLonForGrid(cq_candidates[i].view.locator,latlonlocal);
							float new_dist = latLonDist(latlonlocal, FT8.Local_latlon);
							#if DEBUG
							printf("%s lat:%f lon:%f dist:%f\n",cq_candidates[i].view.from_call,latlonlocal[0],latlonlocal[1],new_dist);
							#endif
							if(dist > new_dist){index_from_ope=i;dist=new_dist;}
						}
					}				
					break;

				case FT8_FILTER_MAX_SNR:
					actusnr = -100;
					for(int i = 0; i<cq_candidate_count;i++){
						if(cq_candidates[i].snr > actusnr){index_from_ope=i;actusnr=cq_candidates[i].snr;}
					}				
					break;

				case FT8_FILTER_MIN_SNR:
					actusnr = 100;
					for(int i = 0; i<cq_candidate_count;i++){
						if(cq_candidates[i].snr < actusnr){index_from_ope=i;actusnr=cq_candidates[i].snr;}
					}				
					break;
				
								
				default:
					index_from_ope = ft8_random_index(cq_candidate_count);
					break;

				}
				
				pthread_mutex_lock(&FT8.TRX_status_lock);
				int selected_session = gcft8_qso_update_from_cq_candidate(&cq_candidates[index_from_ope], rx_now);
				pthread_mutex_unlock(&FT8.TRX_status_lock);
				
				if (selected_session >= 0)
				{
					clear_status_line();
					printf("*Selected for new QSO: \033[1;31m%s at %s on %f hz (seq QSO on %d) \033[0m\n",
						cq_candidates[index_from_ope].view.from_call,
						cq_candidates[index_from_ope].view.locator,
						cq_candidates[index_from_ope].frequency_hz,
						0);
				}
				
			}

			if (!ft8_listen_only_filter_active() && !gcft8_shutdown_requested())
			{
				pthread_mutex_lock(&FT8.TRX_status_lock);
				int selected_session = gcft8_qso_select_tx_session();
				if (selected_session >= 0)
					should_unlock_tx = true;
				pthread_mutex_unlock(&FT8.TRX_status_lock);
			}


			gcft8_flush_pending_logs();

			if (should_unlock_tx)
				unlock_TX_thread();
			
			
			if(!num_decoded){printDateTime_log();printf(" N/A\n");}
			
			#if DEBUG
			printf( "Decoded %d messages\n", num_decoded);
			#endif
			callsign_hash_cache_cleanup(10);
		
		}
		else{
			// printf( "Lock RX thread\n");
			pthread_mutex_lock(&FT8.TRX_status_lock);
			while (!gcft8_shutdown_requested() && FT8.TRX_status != GCFT8_TRX_RX)
				pthread_cond_wait(&FT8.RX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
		
	}

}

static void gcft8_qso_session_clear(gcft8_qso_session_t* session)
{
	if (session == NULL)
		return;

	memset(session, 0, sizeof(*session));
	session->last_rx_snr = GCFT8_SNR_INVALID;
	session->report_sent_snr = GCFT8_SNR_INVALID;
	session->next_tx_seq = -1;
	session->last_tx_seq = -1;
}

static void gcft8_qso_session_begin(gcft8_qso_session_t* session, const char* callsign, time_t now)
{
	if ((session == NULL) || (callsign == NULL))
		return;

	gcft8_qso_session_clear(session);
	session->in_use = true;
	copy_text(session->callsign, sizeof(session->callsign), callsign);
	session->started_at = now;
	session->last_seen_at = now;
	session->last_progress_at = now;
}

static void gcft8_qso_sessions_init(void)
{
	for (int idx = 0; idx < GCFT8_MAX_QSO_SESSIONS; ++idx)
		gcft8_qso_session_clear(&FT8.QSO_sessions[idx]);
}

static bool gcft8_qso_session_expired(const gcft8_qso_session_t* session, time_t now)
{
	if ((session == NULL) || !session->in_use || session->log_pending)
		return false;

	return difftime(now, session->last_seen_at) > GCFT8_QSO_SESSION_TTL_SEC;
}

static int gcft8_qso_find_session(const char* callsign)
{
	if ((callsign == NULL) || (callsign[0] == '\0'))
		return -1;

	for (int idx = 0; idx < GCFT8_MAX_QSO_SESSIONS; ++idx)
	{
		const gcft8_qso_session_t* session = &FT8.QSO_sessions[idx];
		if (session->in_use && (strcmp(session->callsign, callsign) == 0))
			return idx;
	}

	return -1;
}

static int gcft8_qso_allocate_session(time_t now)
{
	int expired_idx = -1;
	int logged_oldest_idx = -1;
	time_t logged_oldest_time = now;

	for (int idx = 0; idx < GCFT8_MAX_QSO_SESSIONS; ++idx)
	{
		gcft8_qso_session_t* session = &FT8.QSO_sessions[idx];

		if (!session->in_use)
			return idx;

		if ((expired_idx < 0) && gcft8_qso_session_expired(session, now))
			expired_idx = idx;

		if (session->logged && !session->log_pending && ((logged_oldest_idx < 0) || (session->last_seen_at < logged_oldest_time)))
		{
			logged_oldest_idx = idx;
			logged_oldest_time = session->last_seen_at;
		}
	}

	if (expired_idx >= 0)
		return expired_idx;

	return logged_oldest_idx;
}

static int gcft8_qso_find_or_create_session(const char* callsign, time_t now)
{
	int idx = gcft8_qso_find_session(callsign);
	gcft8_qso_session_t* session;

	if (idx >= 0)
		return idx;

	idx = gcft8_qso_allocate_session(now);
	if (idx < 0)
		return -1;

	session = &FT8.QSO_sessions[idx];
	gcft8_qso_session_begin(session, callsign, now);
	return idx;
}

static void gcft8_qso_update_locator(gcft8_qso_session_t* session, const char* locator)
{
	if ((session == NULL) || (locator == NULL) || (strlen(locator) != 4))
		return;

	copy_text(session->locator, sizeof(session->locator), locator);
}

static gcft8_qso_action_t gcft8_qso_message_to_action(const gcft8_decoded_message_view_t* view)
{
	gcft8_qso_action_t action;
	action.type = GCFT8_QSO_ACTION_NONE;
	action.tx_seq = -1;

	if ((view == NULL) || !view->is_for_local)
		return action;

	if (strcmp(view->message, "RRR") == 0)
	{
		action.type = GCFT8_QSO_ACTION_TX_AND_LOG_AFTER;
		action.tx_seq = 4;
		return action;
	}

	if (strcmp(view->message, "RR73") == 0)
	{
		action.type = GCFT8_QSO_ACTION_TX_AND_LOG_AFTER;
		action.tx_seq = 4;
		return action;
	}

	if (strcmp(view->message, "73") == 0)
	{
		action.type = GCFT8_QSO_ACTION_LOG_NOW;
		return action;
	}

	if ((view->message[0] == 'R') && gcft8_is_signed_number_text(view->message + 1))
	{
		action.type = GCFT8_QSO_ACTION_TX;
		action.tx_seq = 3;
		return action;
	}

	if (gcft8_is_signed_number_text(view->message))
	{
		action.type = GCFT8_QSO_ACTION_TX;
		action.tx_seq = 2;
		return action;
	}

	if (view->has_locator)
	{
		action.type = GCFT8_QSO_ACTION_TX;
		action.tx_seq = 1;
		return action;
	}

	return action;
}

static bool gcft8_qso_action_requires_existing_session(gcft8_qso_action_t action)
{
	return (action.type == GCFT8_QSO_ACTION_LOG_NOW) || (action.type == GCFT8_QSO_ACTION_TX_AND_LOG_AFTER);
}

static bool gcft8_qso_update_from_direct_candidate(const gcft8_rx_candidate_t* candidate, time_t now, bool* tx_needed)
{
	gcft8_qso_action_t action;
	int idx;
	gcft8_qso_session_t* session;

	if (tx_needed != NULL)
		*tx_needed = false;

	if ((candidate == NULL) || !candidate->view.is_for_local || (candidate->view.from_call[0] == '\0'))
		return false;

	action = gcft8_qso_message_to_action(&candidate->view);
	if (action.type == GCFT8_QSO_ACTION_NONE)
		return false;

	if (gcft8_qso_action_requires_existing_session(action))
		idx = gcft8_qso_find_session(candidate->view.from_call);
	else
		idx = gcft8_qso_find_or_create_session(candidate->view.from_call, now);
	if (idx < 0)
		return false;

	session = &FT8.QSO_sessions[idx];
	session->last_seen_at = now;
	session->last_rx_snr = candidate->snr;
	session->frequency_hz = candidate->frequency_hz;
	if (candidate->view.has_locator && (session->locator[0] == '\0'))
		gcft8_qso_update_locator(session, candidate->view.locator);

	if (action.type == GCFT8_QSO_ACTION_LOG_NOW)
	{
		session->next_tx_seq = -1;
		session->ended_at = now;
		if (!session->logged)
			session->log_pending = true;
		session->last_progress_at = now;
		return true;
	}

	if (session->next_tx_seq != action.tx_seq)
		session->last_progress_at = now;
	if ((action.tx_seq >= 0) && (session->last_tx_seq == action.tx_seq) && (session->same_tx_repeat_count >= GCFT8_MAX_SAME_TX_REPEATS))
	{
		session->next_tx_seq = -1;
		session->log_after_tx_73 = false;
		return true;
	}
	session->next_tx_seq = action.tx_seq;
	session->log_after_tx_73 = action.type == GCFT8_QSO_ACTION_TX_AND_LOG_AFTER;
	if (tx_needed != NULL)
		*tx_needed = action.tx_seq >= 0;
	return true;
}

static int gcft8_qso_update_from_cq_candidate(const gcft8_rx_candidate_t* candidate, time_t now)
{
	int idx;
	gcft8_qso_session_t* session;

	if ((candidate == NULL) || !candidate->view.is_cq || (candidate->view.from_call[0] == '\0'))
		return -1;

	idx = gcft8_qso_find_session(candidate->view.from_call);
	if (idx >= 0)
	{
		session = &FT8.QSO_sessions[idx];
		if (session->logged || session->log_pending)
			return -1;
		session->started_at = now;
		session->ended_at = 0;
		session->report_sent_snr = GCFT8_SNR_INVALID;
		session->log_after_tx_73 = false;
	}
	else
	{
		idx = gcft8_qso_allocate_session(now);
		if (idx < 0)
			return -1;
		session = &FT8.QSO_sessions[idx];
		gcft8_qso_session_begin(session, candidate->view.from_call, now);
	}

	session->last_seen_at = now;
	session->last_rx_snr = candidate->snr;
	session->frequency_hz = candidate->frequency_hz;
	gcft8_qso_update_locator(session, candidate->view.locator);
	session->last_progress_at = now;

	if ((session->last_tx_seq == 0) && (session->same_tx_repeat_count >= GCFT8_MAX_SAME_TX_REPEATS))
	{
		session->next_tx_seq = -1;
		return -1;
	}
	session->next_tx_seq = 0;
	session->log_after_tx_73 = false;
	return idx;
}

static int gcft8_qso_tx_priority(const gcft8_qso_session_t* session)
{
	int score;

	if ((session == NULL) || !session->in_use || (session->next_tx_seq < 0))
		return -1000000;

	score = session->next_tx_seq * 100;
	if ((session->next_tx_seq == session->last_tx_seq) && (session->same_tx_repeat_count >= GCFT8_MAX_SAME_TX_REPEATS))
		score -= 1000;
	return score;
}

static int gcft8_qso_select_tx_session(void)
{
	int best_idx = -1;
	int best_score = -1000000;

	for (int idx = 0; idx < GCFT8_MAX_QSO_SESSIONS; ++idx)
	{
		const gcft8_qso_session_t* session = &FT8.QSO_sessions[idx];
		int score = gcft8_qso_tx_priority(session);

		if (score <= -1000000)
			continue;

		if ((best_idx < 0) || (score > best_score) ||
			((score == best_score) && (session->last_tx_at < FT8.QSO_sessions[best_idx].last_tx_at)) ||
			((score == best_score) && (session->last_tx_at == FT8.QSO_sessions[best_idx].last_tx_at) && (session->last_progress_at < FT8.QSO_sessions[best_idx].last_progress_at)))
		{
			best_idx = idx;
			best_score = score;
		}
	}

	return best_idx;
}

static bool gcft8_qso_build_tx_message(const gcft8_qso_session_t* session, char* dst, size_t dst_size)
{
	int snr;
	int written;

	if ((session == NULL) || !session->in_use || (dst == NULL) || (dst_size == 0) || (session->callsign[0] == '\0'))
		return false;

	snr = session->last_rx_snr == GCFT8_SNR_INVALID ? 0 : session->last_rx_snr;
	switch (session->next_tx_seq)
	{
	case 0:
		written = snprintf(dst, dst_size, "%s %s %s", session->callsign, FT8.Local_CALLSIGN, FT8.Local_LOCATOR);
		break;
	case 1:
		written = snprintf(dst, dst_size, "%s %s %+03d", session->callsign, FT8.Local_CALLSIGN, snr);
		break;
	case 2:
		written = snprintf(dst, dst_size, "%s %s R%+03d", session->callsign, FT8.Local_CALLSIGN, snr);
		break;
	case 3:
		written = snprintf(dst, dst_size, "%s %s RRR", session->callsign, FT8.Local_CALLSIGN);
		break;
	case 4:
		written = snprintf(dst, dst_size, "%s %s 73", session->callsign, FT8.Local_CALLSIGN);
		break;
	default:
		dst[0] = '\0';
		return false;
	}

	if ((written < 0) || ((size_t)written >= dst_size))
	{
		dst[0] = '\0';
		return false;
	}

	return true;
}

static void gcft8_qso_mark_tx_sent(int session_idx, int sent_seq, bool tx_ok, time_t now)
{
	gcft8_qso_session_t* session;

	if (!tx_ok || (session_idx < 0) || (session_idx >= GCFT8_MAX_QSO_SESSIONS) || !FT8.QSO_sessions[session_idx].in_use)
		return;

	session = &FT8.QSO_sessions[session_idx];
	session->last_tx_at = now;
	if (sent_seq == session->last_tx_seq)
		++session->same_tx_repeat_count;
	else
	{
		session->last_tx_seq = sent_seq;
		session->same_tx_repeat_count = 1;
	}

	if ((sent_seq == 1) || (sent_seq == 2))
		session->report_sent_snr = session->last_rx_snr;

	if ((sent_seq == 4) && session->log_after_tx_73)
	{
		session->ended_at = now;
		if (!session->logged)
			session->log_pending = true;
		session->log_after_tx_73 = false;
	}

	if (session->next_tx_seq == sent_seq)
		session->next_tx_seq = -1;
}

static void gcft8_qso_prune_expired(time_t now)
{
	for (int idx = 0; idx < GCFT8_MAX_QSO_SESSIONS; ++idx)
	{
		if (gcft8_qso_session_expired(&FT8.QSO_sessions[idx], now))
			gcft8_qso_session_clear(&FT8.QSO_sessions[idx]);
	}
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

static bool gcft8_build_adif_qso_record(char* dst, size_t dst_size, const char* callsign, const char* locator, int rst_sent, double rf_frequency_hz, time_t started_at, time_t ended_at)
{
	size_t offset = 0;
	char qso_date[9];
	char time_on[7];
	char qso_date_off[9];
	char time_off[7];
	char freq_mhz[16];
	char rst_sent_text[16];

	if ((dst == NULL) || (dst_size == 0) || (callsign == NULL) || (callsign[0] == '\0'))
		return false;
	if (ended_at == 0)
		ended_at = started_at;

	dst[0] = '\0';
	gcft8_format_adif_datetime(started_at, qso_date, time_on);
	gcft8_format_adif_datetime(ended_at, qso_date_off, time_off);
	snprintf(freq_mhz, sizeof(freq_mhz), "%.6f", rf_frequency_hz / 1000000.0);
	rst_sent_text[0] = '\0';
	if (rst_sent != GCFT8_SNR_INVALID)
		snprintf(rst_sent_text, sizeof(rst_sent_text), "%+03d", rst_sent);

	if (!gcft8_append_adif_field(dst, dst_size, &offset, "CALL", callsign))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "MODE", gcft8_adif_mode_name()))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "QSO_DATE", qso_date))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "TIME_ON", time_on))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "QSO_DATE_OFF", qso_date_off))
		return false;
	if (!gcft8_append_adif_field(dst, dst_size, &offset, "TIME_OFF", time_off))
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

static bool gcft8_build_adif_qso_record_from_session(const gcft8_qso_session_t* session, char* dst, size_t dst_size)
{
	if ((session == NULL) || !session->in_use || (session->callsign[0] == '\0'))
		return false;

	return gcft8_build_adif_qso_record(dst, dst_size, session->callsign, session->locator, session->report_sent_snr,
		(double)FT8.Tranceiver_VFOA_Freq + (double)session->frequency_hz, session->started_at, session->ended_at);
}

static bool gcft8_qso_take_pending_log(gcft8_pending_log_t* pending_log)
{
	if (pending_log == NULL)
		return false;

	pending_log->adif_record[0] = '\0';
	pending_log->callsign[0] = '\0';
	pending_log->session_index = -1;

	for (int idx = 0; idx < GCFT8_MAX_QSO_SESSIONS; ++idx)
	{
		gcft8_qso_session_t* session = &FT8.QSO_sessions[idx];
		if (!session->in_use || !session->log_pending || session->log_in_progress)
			continue;

		if (gcft8_build_adif_qso_record_from_session(session, pending_log->adif_record, sizeof(pending_log->adif_record)))
		{
			copy_text(pending_log->callsign, sizeof(pending_log->callsign), session->callsign);
			pending_log->session_index = idx;
			session->log_in_progress = true;
			return true;
		}
	}

	return false;
}

static void gcft8_qso_finish_pending_log(const gcft8_pending_log_t* pending_log, bool success)
{
	gcft8_qso_session_t* session;

	if ((pending_log == NULL) || (pending_log->session_index < 0) || (pending_log->session_index >= GCFT8_MAX_QSO_SESSIONS))
		return;

	pthread_mutex_lock(&FT8.TRX_status_lock);
	session = &FT8.QSO_sessions[pending_log->session_index];
	if (session->in_use && session->log_in_progress && (strcmp(session->callsign, pending_log->callsign) == 0))
	{
		if (success)
		{
			session->logged = true;
			session->log_pending = false;
		}
		session->log_in_progress = false;
	}
	pthread_mutex_unlock(&FT8.TRX_status_lock);
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

static bool log_adif_qso(const char* adif_record)
{
	FILE *fptr;
	bool write_ok;

	if ((adif_record == NULL) || (adif_record[0] == '\0'))
		return false;

	gcft8_ensure_adif_log_file();
	fptr = fopen(FT8.log_file_name,"ab");
	if(fptr == NULL)
	{
		printf("Error with ADIF log file!\n");
		exit(1);
	}

	write_ok = fprintf(fptr, "%s", adif_record) >= 0;
	if (fclose(fptr) != 0)
		write_ok = false;
	if (!write_ok)
	{
		fprintf(stderr, "Error writing ADIF log file %s\n", FT8.log_file_name);
		return false;
	}

	clear_status_line();
	printf("\033[1;32mLogged QSO to %s: %s\033[0m", FT8.log_file_name, adif_record);
	if(FT8.beep_on_log){putchar('\07');putchar('\a');}
	return true;
}

static void gcft8_uppercase_text(char* text);

static void log_qso_to_filter_table(const char* callsign)
{
	char callsign_for_filter[20];

	if ((callsign == NULL) || (callsign[0] == '\0'))
		return;

	copy_text(callsign_for_filter, sizeof(callsign_for_filter), callsign);
	gcft8_uppercase_text(callsign_for_filter);
	ht_insert(ht_callsigntable_for_filter, callsign_for_filter);
}

static void gcft8_flush_pending_logs(void)
{
	for (;;)
	{
		gcft8_pending_log_t pending_log;
		bool has_log;

		pthread_mutex_lock(&FT8.TRX_status_lock);
		has_log = gcft8_qso_take_pending_log(&pending_log);
		pthread_mutex_unlock(&FT8.TRX_status_lock);

		if (!has_log)
			break;

		bool logged = log_adif_qso(pending_log.adif_record);
		if (logged)
			log_qso_to_filter_table(pending_log.callsign);
		gcft8_qso_finish_pending_log(&pending_log, logged);
		if (!logged)
			break;
	}
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

static void load_qso_filter_from_adif(void)
{
	FILE *fptr;
	long file_size;
	char* content;
	size_t bytes_read;
	gcft8_adif_record_t record;

	ht_callsigntable_for_filter = ht_create_table();
	if (ht_callsigntable_for_filter == NULL)
	{
		fprintf(stderr, "Out of memory while creating QSO filter table\n");
		exit(1);
	}
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
	printf("QSO filter table initialised with %zu entry from %s.\n", ht_callsigntable_for_filter->count, FT8.log_file_name);
	#if DEBUG
	print_table(ht_callsigntable_for_filter);
	#endif
}

static void TX_FT8()
{
	while(!gcft8_shutdown_requested()){
		pthread_mutex_lock(&FT8.TRX_status_lock);
		gcft8_trx_state_t stat = FT8.TRX_status;
		pthread_mutex_unlock(&FT8.TRX_status_lock);
		if((stat == GCFT8_TRX_TX) && !gcft8_shutdown_requested()){
			
			bool tx_completed = false;
			char tx_message[GCFT8_TX_MESSAGE_TEXT_SIZE];
			char tx_callsign[sizeof(FT8.QSO_sessions[0].callsign)];
			float tx_frequency = 0.0f;
			int tx_index = -1;
			int tx_session_index = -1;
			
			tx_message[0] = '\0';
			tx_callsign[0] = '\0';

			pthread_mutex_lock(&FT8.TRX_status_lock);
			tx_session_index = gcft8_qso_select_tx_session();
			if (tx_session_index >= 0)
			{
				gcft8_qso_session_t* session = &FT8.QSO_sessions[tx_session_index];
				tx_index = session->next_tx_seq;
				if ((tx_index > -1) && gcft8_qso_build_tx_message(session, tx_message, sizeof(tx_message)) && !gcft8_shutdown_requested())
				{
					copy_text(tx_callsign, sizeof(tx_callsign), session->callsign);
					tx_frequency = session->frequency_hz;
				}
				else
				{
					tx_index = -1;
				}
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
				uint8_t tones[FT2_NN];
				if (mode_cfg->protocol == FTX_PROTOCOL_FT2)
					ft2_encode(tx_msg.payload, tones);
				else if (mode_cfg->protocol == FTX_PROTOCOL_FT4)
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
				float* signal;
				int16_t* raw_data;
				bool waveform_ok;

				if (!gcft8_tx_context_ensure(&gcft8_tx_context, (size_t)num_total_samples))
				{
					fprintf(stderr, "Out of memory while allocating TX buffers\n");
					exit(1);
				}

				signal = gcft8_tx_context.signal;
				raw_data = gcft8_tx_context.pcm;
				
				for (int i = 0; i < num_total_samples; i++)
				{
					signal[i] = 0.0f;
				}
				
				if (mode_cfg->protocol == FTX_PROTOCOL_FT2)
					waveform_ok = synth_ft2_gfsk(tones, num_tones, tx_frequency, symbol_bt, symbol_period, sample_rate, signal + num_lead_silence, &gcft8_tx_context.gfsk_scratch);
				else
					waveform_ok = synth_gfsk(tones, num_tones, tx_frequency, symbol_bt, symbol_period, sample_rate, signal + num_lead_silence, &gcft8_tx_context.gfsk_scratch);
				if (!waveform_ok)
				{
					fprintf(stderr, "Out of memory while synthesizing TX waveform\n");
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
					unlock_RX_thread();
					continue;
				}
				
				if (!wait_for_tx_slot_start(slot_time, mode_cfg->tx_late_grace))
				{
					tranceiver_rtx(GCFT8_TRX_RX);
					break;
				}

				print_slot_separator(mode_cfg, "TX");
				printDateTime_ms();printf(" start send message\n");
				
				tranceiver_rtx(GCFT8_TRX_TX);
								
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
				tx_completed = (count == num_total_samples) && !gcft8_shutdown_requested();
				
				tranceiver_rtx(GCFT8_TRX_RX);

				printDateTime_ms();printf(" stop send message\n");

			}

			pthread_mutex_lock(&FT8.TRX_status_lock);
			gcft8_qso_mark_tx_sent(tx_session_index, tx_index, tx_completed, time(NULL));
			pthread_mutex_unlock(&FT8.TRX_status_lock);

			gcft8_flush_pending_logs();

			unlock_RX_thread();
			
		}
		else{
			#if DEBUG
			printf( "Lock TX thread\n");
			#endif
			pthread_mutex_lock(&FT8.TRX_status_lock);
			while (!gcft8_shutdown_requested() && FT8.TRX_status != GCFT8_TRX_TX)
				pthread_cond_wait(&FT8.TX_status_cond, &FT8.TRX_status_lock);
			pthread_mutex_unlock(&FT8.TRX_status_lock);		
		}
	}

	tranceiver_rtx(GCFT8_TRX_RX);
}

static void * Thread_RX(void *arg) {
	(void)arg;
	RX_FT8();
	return NULL;
}

static void * Thread_TX(void *arg) {
	(void)arg;
	TX_FT8();
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//serial and tranceiver commands

static int serial_init(){
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

static void serial_deinit(void)
{
	if (serial.port != NULL)
	{
		cssl_putstring(serial.port,"RX;");
		cssl_close(serial.port);
		serial.port = NULL;
	}
	cssl_stop();
}

static void get_serial_rep(char rep[16])
{
	int i;
	for(i=0;i<15;i++){
		rep[i]=((char)cssl_getchar(serial.port));
		if(rep[i]==';'){break;}
	}
	rep[(i < 15) ? i + 1 : 15]=0;
}

static bool tranceiver_set_freq(int freq)
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

static void tranceiver_init()
{
	if (serial.port == NULL)
		return;

	cssl_putstring(serial.port,"FR0;"); //Set receive on VFO_A
	cssl_putstring(serial.port,"FT0;"); //Set Transmit on VFO_A
	cssl_putstring(serial.port,"Q10;"); //Set USB mode
	if(!tranceiver_set_freq(FT8.Tranceiver_VFOA_Freq)){printf("Unable to communicate with tranceiver!");exit(-1);} //Set frequency stored
}

static void tranceiver_rtx(gcft8_trx_state_t ptt)
{
	if (serial.port == NULL)
		return;

	if (ptt == GCFT8_TRX_TX){
		#if DEBUG
		printf("tranceiver ptt on\n");
		#endif
		cssl_putstring(serial.port,"TX;");
		} //Set receiving
	else if (ptt == GCFT8_TRX_RX){
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
	FT8.TRX_status = GCFT8_TRX_RX;
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
	tranceiver_rtx(GCFT8_TRX_RX);
	serial_deinit();
	playback_audioDeInit();
	capture_audioDeInit();
	gcft8_tx_context_free(&gcft8_tx_context);
	gcft8_rx_context_free(&gcft8_rx_context);
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
		"  --mode <ft8|ft4|ft2>       Digital mode (default: ft8)\n"
		"  --sound-device <device>    ALSA capture and playback device, prefer plughw (default: default)\n"
		"  --callsign <callsign>      Your callsign (default: F4JJJ)\n"
		"  --locator <locator>        Your Maidenhead locator (default: JN38)\n"
		"  --frequency <hz>           TRX frequency in Hz, exclusive with --band\n"
		"  --band <band>              Mode-specific band frequency, exclusive with --frequency; suffix m is allowed\n"
		"  --serial-device <device>   Transceiver serial device (default: /dev/ttyACM0)\n"
		"  --filter <mode>            Operating/filter mode (default: 0, listen only)\n"
		"  --snr-min <snr>            Reject CQ candidates below this SNR, for example -18\n"
		"  --only-prefix <list>       Only auto-select CQ callsigns with these prefixes, for example JA,VK,ZL\n"
		"  --only-sp-tag <list>       Only auto-select special CQ tags, for example POTA,SOTA,DX\n"
		"  --only-locator-zone <list> Only auto-select CQ locators in two-letter zones, for example BP:FL\n"
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
		"  Band       FT8 Hz    FT4 Hz    FT2 Hz\n",
		program_name,
		program_name);

	for (size_t idx = 0; idx < sizeof(gcft8_band_frequencies) / sizeof(gcft8_band_frequencies[0]); ++idx)
	{
		char ft4_text[16];
		char ft2_text[16];

		if (gcft8_band_frequencies[idx].ft4_frequency_hz > 0)
			snprintf(ft4_text, sizeof(ft4_text), "%d", gcft8_band_frequencies[idx].ft4_frequency_hz);
		else
			copy_text(ft4_text, sizeof(ft4_text), "n/a");

		if (gcft8_band_frequencies[idx].ft2_frequency_hz > 0)
			snprintf(ft2_text, sizeof(ft2_text), "%d", gcft8_band_frequencies[idx].ft2_frequency_hz);
		else
			copy_text(ft2_text, sizeof(ft2_text), "n/a");

		fprintf(stream, "  %2sm  %8d  %8s  %8s\n", gcft8_band_frequencies[idx].band, gcft8_band_frequencies[idx].ft8_frequency_hz, ft4_text, ft2_text);
	}

	fprintf(stream,
		"\n"
		"Display colors:\n"
		"  Cyan     RX/TX slot separator\n"
		"  Red      Local station related message\n"
		"  Blue     CQ candidate\n"
		"  Yellow   Already worked callsign\n"
		"  Magenta  Filtered CQ, missing locator/callsign, or CQ rejected by optional filters\n");
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
	CLI_OPTION_SNR_MIN,
	CLI_OPTION_ONLY_PREFIX,
	CLI_OPTION_ONLY_SP_TAG,
	CLI_OPTION_ONLY_LOCATOR_ZONE,
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
		{ "snr-min", required_argument, NULL, CLI_OPTION_SNR_MIN },
		{ "only-prefix", required_argument, NULL, CLI_OPTION_ONLY_PREFIX },
		{ "only-sp-tag", required_argument, NULL, CLI_OPTION_ONLY_SP_TAG },
		{ "only-locator-zone", required_argument, NULL, CLI_OPTION_ONLY_LOCATOR_ZONE },
		{ "serial-device", required_argument, NULL, CLI_OPTION_SERIAL_DEVICE },
		{ NULL, 0, NULL, 0 }
	};

	while ((c = getopt_long(argc, argv, "", long_options, NULL)) != -1)
		switch (c)
		{
			case CLI_OPTION_HELP:
				print_usage(argv[0], stdout);
				exit(0);
			case CLI_OPTION_BEEP:
				FT8.beep_on_log = 1;
				break;
			case CLI_OPTION_SOUND_DEVICE:
				sound.capture_sound_device = optarg;
				sound.playback_sound_device = optarg;
				break;
			case CLI_OPTION_CALLSIGN:
				copy_text(FT8.Local_CALLSIGN, sizeof(FT8.Local_CALLSIGN), optarg);
				break;
			case CLI_OPTION_LOCATOR:
				copy_text(FT8.Local_LOCATOR, sizeof(FT8.Local_LOCATOR), optarg);
				break;
			case CLI_OPTION_MODE:
				if (!gcft8_parse_mode(optarg, &gcft8_mode))
				{
					fprintf(stderr, "Invalid mode '%s'. Allowed modes: ft8, ft4, ft2.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
			case CLI_OPTION_FREQUENCY:
				frequency_option_used = true;
				FT8.Tranceiver_VFOA_Freq = atoi(optarg);
				break;
			case CLI_OPTION_BAND:
				band_option_used = true;
				copy_text(selected_band, sizeof(selected_band), optarg);
				break;
			case CLI_OPTION_FILTER:
				if (!ft8_parse_filter_mode(optarg, &FT8.filter_on_cq))
				{
					fprintf(stderr, "Invalid filter '%s'. Allowed filters: 0, 1, 2, 3, 4, 5, 6.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
			case CLI_OPTION_SNR_MIN:
				if (!gcft8_parse_snr_min(optarg, &gcft8_snr_min))
				{
					fprintf(stderr, "Invalid SNR minimum '%s'. Use an integer value, for example -18.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				gcft8_snr_min_enabled = true;
				break;
			case CLI_OPTION_ONLY_PREFIX:
				if (!gcft8_parse_only_prefixes(optarg))
				{
					fprintf(stderr, "Invalid prefix list '%s'. Use comma-separated alphanumeric prefixes, for example JA,VK,ZL.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
			case CLI_OPTION_ONLY_SP_TAG:
				if (!gcft8_parse_only_sp_tags(optarg))
				{
					fprintf(stderr, "Invalid special CQ tag list '%s'. Use comma-separated 1-4 letter tags or 3-digit tags, for example POTA,SOTA,DX.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
			case CLI_OPTION_ONLY_LOCATOR_ZONE:
				if (!gcft8_parse_locator_zones(optarg))
				{
					fprintf(stderr, "Invalid locator zone list '%s'. Use comma-separated LL:LL ranges with A-R letters, for example BP:FL,IO:KM.\n", optarg);
					print_usage(argv[0], stderr);
					exit(1);
				}
				break;
			case CLI_OPTION_SERIAL_DEVICE:
				copy_text(serial.pathname, sizeof(serial.pathname), optarg);
				break;
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
	gcft8_tx_context_init(&gcft8_tx_context);
	gcft8_rx_context_init(&gcft8_rx_context);
	gcft8_qso_sessions_init();
	
	latLonForGrid(FT8.Local_LOCATOR,FT8.Local_latlon);
	callsign_hash_cache_init();
	gcft8_update_adif_log_filename(FT8.Local_CALLSIGN, FT8.log_file_name, sizeof(FT8.log_file_name));
	load_qso_filter_from_adif();
	
	capture_audioInit();
	playback_audioInit();
	
	const gcft8_mode_config_t* startup_mode_cfg = gcft8_current_mode_config();
	char startup_snr_min_text[32];
	char startup_prefix_filter[256];
	char startup_sp_tag_filter[256];
	char startup_locator_zone_filter[256];
	if (gcft8_snr_min_enabled)
		snprintf(startup_snr_min_text, sizeof(startup_snr_min_text), "%d", gcft8_snr_min);
	else
		copy_text(startup_snr_min_text, sizeof(startup_snr_min_text), "off");
	gcft8_format_only_prefixes(startup_prefix_filter, sizeof(startup_prefix_filter));
	gcft8_format_only_sp_tags(startup_sp_tag_filter, sizeof(startup_sp_tag_filter));
	gcft8_format_locator_zones(startup_locator_zone_filter, sizeof(startup_locator_zone_filter));
	printf("Starting with this:\n"
		"-mode is %s\n"
		"-set Freq to %d\n"
		"-your callsign is %s\n"
		"-your locator is %s\n"
		"-TRX serial port is %s\n"
		"-Sound device is %s\n"
		"-ADIF log file is %s\n"
		"-CQ filter method is %d (%s)\n"
		"-SNR minimum filter is %s\n"
		"-Prefix filter is %s\n"
		"-Special CQ tag filter is %s\n"
		"-Locator zone filter is %s\n"
		"-Beep on log %d\n",
		startup_mode_cfg->name,
		FT8.Tranceiver_VFOA_Freq,
		FT8.Local_CALLSIGN,
		FT8.Local_LOCATOR,
		serial.pathname,
		sound.capture_sound_device,
		FT8.log_file_name,
		FT8.filter_on_cq,
		ft8_filter_mode_name(FT8.filter_on_cq),
		startup_snr_min_text,
		startup_prefix_filter,
		startup_sp_tag_filter,
		startup_locator_zone_filter,
		FT8.beep_on_log);
		
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
