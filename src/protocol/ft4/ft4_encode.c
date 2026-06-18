#include "protocol/ft4/ft4_encode.h"

#include "protocol/ftx/constants.h"
#include "protocol/ftx/crc.h"
#include "protocol/ftx/encode_common.h"

void ft4_encode(const uint8_t* payload, uint8_t* tones)
{
    uint8_t a91[FTX_LDPC_K_BYTES];
    uint8_t payload_xor[10];

    for (int i = 0; i < 10; ++i)
    {
        payload_xor[i] = payload[i] ^ kFT4_XOR_sequence[i];
    }

    ftx_add_crc(payload_xor, a91);

    uint8_t codeword[FTX_LDPC_N_BYTES];
    ftx_encode174(a91, codeword);

    uint8_t mask = 0x80u;
    int i_byte = 0;
    for (int i_tone = 0; i_tone < FT4_NN; ++i_tone)
    {
        if ((i_tone == 0) || (i_tone == 104))
        {
            tones[i_tone] = 0;
        }
        else if ((i_tone >= 1) && (i_tone < 5))
        {
            tones[i_tone] = kFT4_Costas_pattern[0][i_tone - 1];
        }
        else if ((i_tone >= 34) && (i_tone < 38))
        {
            tones[i_tone] = kFT4_Costas_pattern[1][i_tone - 34];
        }
        else if ((i_tone >= 67) && (i_tone < 71))
        {
            tones[i_tone] = kFT4_Costas_pattern[2][i_tone - 67];
        }
        else if ((i_tone >= 100) && (i_tone < 104))
        {
            tones[i_tone] = kFT4_Costas_pattern[3][i_tone - 100];
        }
        else
        {
            uint8_t bits2 = 0;

            if (codeword[i_byte] & mask)
                bits2 |= 2;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }
            if (codeword[i_byte] & mask)
                bits2 |= 1;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }
            tones[i_tone] = kFT4_Gray_map[bits2];
        }
    }
}
