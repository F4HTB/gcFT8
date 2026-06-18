#ifndef GCFT8_PROTOCOL_FT2_WAVEFORM_H
#define GCFT8_PROTOCOL_FT2_WAVEFORM_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

void synth_ft2_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal);

#ifdef __cplusplus
}
#endif

#endif
