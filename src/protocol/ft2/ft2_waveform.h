#ifndef GCFT8_PROTOCOL_FT2_WAVEFORM_H
#define GCFT8_PROTOCOL_FT2_WAVEFORM_H

#include <stdbool.h>
#include <stdint.h>

#include "dsp/gfsk.h"

#ifdef __cplusplus
extern "C"
{
#endif

bool synth_ft2_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal, gfsk_scratch_t* scratch);

#ifdef __cplusplus
}
#endif

#endif
