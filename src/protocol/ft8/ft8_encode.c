#include "protocol/ft8/ft8_encode.h"

#include "protocol/ftx/constants.h"
#include "protocol/ftx/crc.h"
#include "protocol/ftx/encode_common.h"

void ft8_encode(const uint8_t* payload, uint8_t* tones)
{
    uint8_t a91[FTX_LDPC_K_BYTES];

    ftx_add_crc(payload, a91);

    uint8_t codeword[FTX_LDPC_N_BYTES];
    ftx_encode174(a91, codeword);

    uint8_t mask = 0x80u;
    int i_byte = 0;
    for (int i_tone = 0; i_tone < FT8_NN; ++i_tone)
    {
        if ((i_tone >= 0) && (i_tone < 7))
        {
            tones[i_tone] = kFT8_Costas_pattern[i_tone];
        }
        else if ((i_tone >= 36) && (i_tone < 43))
        {
            tones[i_tone] = kFT8_Costas_pattern[i_tone - 36];
        }
        else if ((i_tone >= 72) && (i_tone < 79))
        {
            tones[i_tone] = kFT8_Costas_pattern[i_tone - 72];
        }
        else
        {
            uint8_t bits3 = 0;

            if (codeword[i_byte] & mask)
                bits3 |= 4;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }
            if (codeword[i_byte] & mask)
                bits3 |= 2;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }
            if (codeword[i_byte] & mask)
                bits3 |= 1;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }

            tones[i_tone] = kFT8_Gray_map[bits3];
        }
    }
}
