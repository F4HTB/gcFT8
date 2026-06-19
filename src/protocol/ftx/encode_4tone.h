#ifndef GCFT8_PROTOCOL_FTX_ENCODE_4TONE_H
#define GCFT8_PROTOCOL_FTX_ENCODE_4TONE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

void ftx_encode_4tone_105(const uint8_t* payload, uint8_t* tones);

#ifdef __cplusplus
}
#endif

#endif
