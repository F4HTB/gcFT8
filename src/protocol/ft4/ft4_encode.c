#include "protocol/ft4/ft4_encode.h"

#include "protocol/ftx/encode_4tone.h"

void ft4_encode(const uint8_t* payload, uint8_t* tones)
{
    ftx_encode_4tone_105(payload, tones);
}
