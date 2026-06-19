#include "protocol/ft2/ft2_encode.h"

#include "protocol/ftx/constants.h"
#include "protocol/ftx/encode_4tone.h"

void ft2_encode(const uint8_t* payload, uint8_t* tones)
{
#if FT2_NN != FT4_NN
#error "FT2 encoder expects the shared 105-symbol 4-tone layout"
#endif
    ftx_encode_4tone_105(payload, tones);
}
