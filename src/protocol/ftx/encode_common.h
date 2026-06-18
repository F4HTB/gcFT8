#ifndef GCFT8_PROTOCOL_FTX_ENCODE_COMMON_H
#define GCFT8_PROTOCOL_FTX_ENCODE_COMMON_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

void ftx_encode174(const uint8_t* message, uint8_t* codeword);

#ifdef __cplusplus
}
#endif

#endif
