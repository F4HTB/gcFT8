#include "encode_common.h"
#include "constants.h"

// Returns 1 if an odd number of bits are set in x, zero otherwise
static uint8_t parity8(uint8_t x)
{
    x ^= x >> 4;  // a b c d ae bf cg dh
    x ^= x >> 2;  // a b ac bd cae dbf aecg bfdh
    x ^= x >> 1;  // a ab bac acbd bdcae caedbf aecgbfdh
    return x % 2; // modulo 2
}

// Encode via LDPC a 91-bit message and return a 174-bit codeword.
// The generator matrix has dimensions (87,87).
// The code is a (174,91) regular LDPC code with column weight 3.
// Arguments:
// [IN] message   - array of 91 bits stored as 12 bytes (MSB first)
// [OUT] codeword - array of 174 bits stored as 22 bytes (MSB first)
void ftx_encode174(const uint8_t* message, uint8_t* codeword)
{
    // This implementation accesses the generator bits straight from the packed binary representation in kFTX_LDPC_generator

    // Fill the codeword with message and zeros, as we will only update binary ones later
    for (int j = 0; j < FTX_LDPC_N_BYTES; ++j)
    {
        codeword[j] = (j < FTX_LDPC_K_BYTES) ? message[j] : 0;
    }

    // Compute the byte index and bit mask for the first checksum bit
    uint8_t col_mask = (0x80u >> (FTX_LDPC_K % 8u)); // bitmask of current byte
    uint8_t col_idx = FTX_LDPC_K_BYTES - 1;          // index into byte array

    // Compute the LDPC checksum bits and store them in codeword
    for (int i = 0; i < FTX_LDPC_M; ++i)
    {
        // Fast implementation of bitwise multiplication and parity checking
        // Normally nsum would contain the result of dot product between message and kFTX_LDPC_generator[i],
        // but we only compute the sum modulo 2.
        uint8_t nsum = 0;
        for (int j = 0; j < FTX_LDPC_K_BYTES; ++j)
        {
            uint8_t bits = message[j] & kFTX_LDPC_generator[i][j]; // bitwise AND (bitwise multiplication)
            nsum ^= parity8(bits);                                 // bitwise XOR (addition modulo 2)
        }

        // Set the current checksum bit in codeword if nsum is odd
        if (nsum % 2)
        {
            codeword[col_idx] |= col_mask;
        }

        // Update the byte index and bit mask for the next checksum bit
        col_mask >>= 1;
        if (col_mask == 0)
        {
            col_mask = 0x80u;
            ++col_idx;
        }
    }
}
