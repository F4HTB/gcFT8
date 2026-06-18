#include "protocol/ft2/ft2_encode.h"

#include "protocol/ft4/ft4_encode.h"
#include "protocol/ftx/constants.h"

void ft2_encode(const uint8_t* payload, uint8_t* tones)
{
#if FT2_NN != FT4_NN
#error "FT2 encoder reuses the FT4 105-symbol layout"
#endif
    ft4_encode(payload, tones);
}
