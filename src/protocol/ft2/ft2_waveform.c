#include "protocol/ft2/ft2_waveform.h"

#include "util/math_constants.h"

#include <math.h>

bool synth_ft2_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal, gfsk_scratch_t* scratch)
{
    int n_spsym;
    int n_wave;
    int n_payload_sym;
    size_t dphi_len;
    size_t pulse_len;
    float* dphi;
    float* pulse;
    float hmod = 1.0f;
    float dphi_peak;

    if (n_sym < 3)
        return synth_gfsk(symbols, n_sym, f0, symbol_bt, symbol_period, signal_rate, signal, scratch);

    if ((symbols == NULL) || (signal == NULL) || (scratch == NULL) || (symbol_period <= 0.0f) || (signal_rate <= 0))
        return false;

    n_spsym = (int)(0.5f + signal_rate * symbol_period);
    if (n_spsym <= 0)
        return false;

    n_wave = n_sym * n_spsym;
    n_payload_sym = n_sym - 2;
    dphi_len = (size_t)n_wave;
    pulse_len = 3u * (size_t)n_spsym;
    if (!gfsk_scratch_ensure(scratch, dphi_len, pulse_len))
        return false;

    dphi = scratch->dphi;
    pulse = scratch->pulse;
    dphi_peak = 2 * M_PI * hmod / n_spsym;

    for (int i = 0; i < n_wave; ++i)
    {
        dphi[i] = 2 * M_PI * f0 / signal_rate;
    }

    gfsk_pulse(n_spsym, symbol_bt, pulse);

    // Decodium FT2 encodes 103 sync/data symbols and emits a 105-symbol waveform span.
    for (int i = 0; i < n_payload_sym; ++i)
    {
        int ib = i * n_spsym;
        for (int j = 0; j < 3 * n_spsym && ib + j < n_wave; ++j)
        {
            dphi[ib + j] += dphi_peak * symbols[i + 1] * pulse[j];
        }
    }

    float phi = 0;
    for (int k = 0; k < n_wave; ++k)
    {
        signal[k] = sinf(phi);
        phi = fmodf(phi + dphi[k], 2 * M_PI);
    }

    for (int i = 0; i < n_spsym; ++i)
    {
        float env_up = (1.0f - cosf(2 * M_PI * i / (2 * n_spsym))) / 2.0f;
        float env_down = (1.0f + cosf(2 * M_PI * i / (2 * n_spsym))) / 2.0f;
        signal[i] *= env_up;
        signal[((n_sym - 1) * n_spsym) + i] *= env_down;
    }

    return true;
}
