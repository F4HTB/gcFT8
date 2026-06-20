#ifndef GCFTX_PROTOCOL_FT2_ENCODE_H
#define GCFTX_PROTOCOL_FT2_ENCODE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

void ft2_encode(const uint8_t* payload, uint8_t* tones);

#ifdef __cplusplus
}
#endif

#endif
