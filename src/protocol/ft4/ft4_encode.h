#ifndef GCFT8_PROTOCOL_FT4_ENCODE_H
#define GCFT8_PROTOCOL_FT4_ENCODE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

void ft4_encode(const uint8_t* payload, uint8_t* tones);

#ifdef __cplusplus
}
#endif

#endif
