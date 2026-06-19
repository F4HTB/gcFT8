#include "dsp/gfsk.h"

#include "util/math_constants.h"

#include <math.h>
#include <stdlib.h>

#define GFSK_CONST_K 5.336446f ///< == pi * sqrt(2 / log(2))

void gfsk_scratch_init(gfsk_scratch_t* scratch)
{
    if (scratch == NULL)
        return;

    scratch->dphi = NULL;
    scratch->dphi_capacity = 0;
    scratch->pulse = NULL;
    scratch->pulse_capacity = 0;
}

void gfsk_scratch_free(gfsk_scratch_t* scratch)
{
    if (scratch == NULL)
        return;

    free(scratch->dphi);
    free(scratch->pulse);
    gfsk_scratch_init(scratch);
}

bool gfsk_scratch_ensure(gfsk_scratch_t* scratch, size_t dphi_capacity, size_t pulse_capacity)
{
    if (scratch == NULL)
        return false;

    if (dphi_capacity > scratch->dphi_capacity)
    {
        float* dphi = (float*)realloc(scratch->dphi, dphi_capacity * sizeof(scratch->dphi[0]));
        if (dphi == NULL)
            return false;

        scratch->dphi = dphi;
        scratch->dphi_capacity = dphi_capacity;
    }

    if (pulse_capacity > scratch->pulse_capacity)
    {
        float* pulse = (float*)realloc(scratch->pulse, pulse_capacity * sizeof(scratch->pulse[0]));
        if (pulse == NULL)
            return false;

        scratch->pulse = pulse;
        scratch->pulse_capacity = pulse_capacity;
    }

    return true;
}

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

bool synth_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal, gfsk_scratch_t* scratch)
{
    int n_spsym;
    int n_wave;
    size_t dphi_len;
    size_t pulse_len;
    float* dphi;
    float* pulse;
    float hmod = 1.0f;
    float dphi_peak;

    if ((symbols == NULL) || (signal == NULL) || (scratch == NULL) || (n_sym <= 0) || (symbol_period <= 0.0f) || (signal_rate <= 0))
        return false;

    n_spsym = (int)(0.5f + signal_rate * symbol_period);
    if (n_spsym <= 0)
        return false;

    n_wave = n_sym * n_spsym;
    dphi_len = (size_t)n_wave + (2u * (size_t)n_spsym);
    pulse_len = 3u * (size_t)n_spsym;
    if (!gfsk_scratch_ensure(scratch, dphi_len, pulse_len))
        return false;

    dphi = scratch->dphi;
    pulse = scratch->pulse;
    dphi_peak = 2 * M_PI * hmod / n_spsym;

    for (size_t i = 0; i < dphi_len; ++i)
    {
        dphi[i] = 2 * M_PI * f0 / signal_rate;
    }

    gfsk_pulse(n_spsym, symbol_bt, pulse);

    for (int i = 0; i < n_sym; ++i)
    {
        int ib = i * n_spsym;
        for (int j = 0; j < 3 * n_spsym; ++j)
        {
            dphi[j + ib] += dphi_peak * symbols[i] * pulse[j];
        }
    }

    for (int j = 0; j < 2 * n_spsym; ++j)
    {
        dphi[j] += dphi_peak * pulse[j + n_spsym] * symbols[0];
        dphi[j + n_sym * n_spsym] += dphi_peak * pulse[j] * symbols[n_sym - 1];
    }

    float phi = 0;
    for (int k = 0; k < n_wave; ++k)
    {
        signal[k] = sinf(phi);
        phi = fmodf(phi + dphi[k + n_spsym], 2 * M_PI);
    }

    int n_ramp = n_spsym / 8;
    for (int i = 0; i < n_ramp; ++i)
    {
        float env = (1 - cosf(2 * M_PI * i / (2 * n_ramp))) / 2;
        signal[i] *= env;
        signal[n_wave - 1 - i] *= env;
    }

    return true;
}
