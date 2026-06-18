#include "protocol/ft2/ft2_waveform.h"

#include "dsp/gfsk.h"
#include "util/math_constants.h"

#include <math.h>

void synth_ft2_gfsk(const uint8_t* symbols, int n_sym, float f0, float symbol_bt, float symbol_period, int signal_rate, float* signal)
{
    if (n_sym < 3)
    {
        synth_gfsk(symbols, n_sym, f0, symbol_bt, symbol_period, signal_rate, signal);
        return;
    }

    int n_spsym = (int)(0.5f + signal_rate * symbol_period);
    int n_wave = n_sym * n_spsym;
    int n_payload_sym = n_sym - 2;
    float hmod = 1.0f;
    float dphi_peak = 2 * M_PI * hmod / n_spsym;
    float dphi[n_wave];
    float pulse[3 * n_spsym];

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
}
