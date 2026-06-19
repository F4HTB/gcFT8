#ifndef GCFT8_DSP_GFSK_H
#define GCFT8_DSP_GFSK_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
    float* dphi;
    size_t dphi_capacity;
    float* pulse;
    size_t pulse_capacity;
} gfsk_scratch_t;

void gfsk_scratch_init(gfsk_scratch_t* scratch);
void gfsk_scratch_free(gfsk_scratch_t* scratch);
bool gfsk_scratch_ensure(gfsk_scratch_t* scratch, size_t dphi_capacity, size_t pulse_capacity);
void gfsk_pulse(int n_spsym, float symbol_bt, float* pulse);
bool synth_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal, gfsk_scratch_t* scratch);

#ifdef __cplusplus
}
#endif

#endif
