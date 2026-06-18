#ifndef GCFT8_DSP_GFSK_H
#define GCFT8_DSP_GFSK_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

void gfsk_pulse(int n_spsym, float symbol_bt, float* pulse);
void synth_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal);

#ifdef __cplusplus
}
#endif

#endif
