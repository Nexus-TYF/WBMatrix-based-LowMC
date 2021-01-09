#include "WBMatrix/WBMatrix.h"
#ifdef __GNUC__
#include <x86intrin.h>
#endif
#ifdef _MSC_VER
#include <intrin.h>
#endif
#pragma intrinsic(__rdtsc)

//Repeat test times and calculate on average for accuracy
#define TEST 100

//CPU cycles set start;
uint64_t start_rdtsc()
{
	return __rdtsc();
}

//CPU cycles set end;
uint64_t end_rdtsc()
{
    return __rdtsc();
}

#define blocksize 256
#define keysize 80
#define rounds 12
#define numofboxes 49

M256 LinMatrices[rounds];
M256 invLinMatrices[rounds];
V256 roundconstants[rounds];
M256 KeyMatrices[rounds + 1];
V256 roundkeys[rounds + 1];
V256 key;

const uint8_t Sbox[] =
        {0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02};
const uint8_t invSbox[] =
        {0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04};

void keyschedule () {
    for (unsigned r = 0; r <= rounds; ++r) {
        MatMulVecM256(KeyMatrices[r], key, &roundkeys[r]);
    }
    return;
}
void set_key (V256 k) {
    key.V[0] = k.V[0];
    key.V[1] = k.V[1];
    key.V[2] = k.V[2];
    key.V[3] = k.V[3];
    keyschedule();
}
void instantiate_LowMC () {
    // Create LinMatrices and invLinMatrices
    for (unsigned r = 0; r < rounds; ++r) {
        // Create matrix
        genMatpairM256(&LinMatrices[r], &invLinMatrices[r]);

        // Create roundconstants
        randV256(&roundconstants[r]);
    }

    // Create KeyMatrices
    for (unsigned r = 0; r <= rounds; ++r) {
        // Create matrix
        do {
            randM256(&KeyMatrices[r]);
        // Repeat if matrix is not of maximal rank
        } while ( isinvertM256(KeyMatrices[r]) );
    }
    // key schedule
    V256 k;
    k.V[0] = 0x0123456789abcdef;
    set_key(k);
    
    return;
}
V256 Substitution (V256 message) {
    V256 temp;
    initV256(&temp);
    //Get the identity part of the message
    if(numofboxes <= 21)
    {
        temp.V[0] = (message.V[0] >> 3 * numofboxes);
        temp.V[1] ^= message.V[1];
        temp.V[2] ^= message.V[2];
        temp.V[3] ^= message.V[3];
        //Get the rest through the Sboxes
        for (unsigned i = 1; i <= numofboxes; ++i) {
            temp.V[0] <<= 3;
            temp.V[0] ^= Sbox[ ((message.V[0] >> 3 * (numofboxes - i)) & 0x7) ];
        }
    }
    else if (numofboxes <= 42)
    {
        temp.V[1] ^= (message.V[1] >> (3 * (numofboxes - 22) + 2));
        temp.V[2] ^= message.V[2];
        temp.V[3] ^= message.V[3];
        //Get the rest through the Sboxes
        temp.V[0] ^= Sbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ];
        for (unsigned i = 1; i <= 21; ++i) {
            temp.V[0] <<= 3;
            temp.V[0] ^= Sbox[ ((message.V[0] >> 3 * (21 - i)) & 0x7) ];
        }
        for (unsigned i = 1; i <= (numofboxes - 22); ++i) {
            temp.V[1] <<= 3;
            temp.V[1] ^= Sbox[ ((message.V[1] >> 3 * (numofboxes - 22 - i) + 2) & 0x7) ];
        }
        temp.V[1] <<= 2;
        temp.V[1] ^= ((Sbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ] & 0x6) >> 1);
    }
    else if (numofboxes <= 64)
    {
        temp.V[3] ^= message.V[3];
        //Get the rest through the Sboxes
        temp.V[0] ^= Sbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ];
        for (unsigned i = 1; i <= 21; ++i) {
            temp.V[0] <<= 3;
            temp.V[0] ^= Sbox[ ((message.V[0] >> 3 * (21 - i)) & 0x7) ];
        }
        temp.V[1] ^= Sbox[ (message.V[2] & 0x1 << 2) ^ (message.V[1] >> 62 & 0x3) ];
        for (unsigned i = 1; i <= 20; ++i) {
            temp.V[1] <<= 3;
            temp.V[1] ^= Sbox[ ((message.V[1] >> 3 * (numofboxes - 22 - i) + 2) & 0x7) ];
        }
        temp.V[1] <<= 2;
        temp.V[1] ^= ((Sbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ] & 0x6) >> 1);
        for (unsigned i = 1; i <= (numofboxes - 43); ++i) {
            temp.V[2] <<= 3;
            temp.V[2] ^= Sbox[ ((message.V[2] >> 3 * (numofboxes - 43 - i) + 1) & 0x7) ];
        }
        temp.V[2] <<= 1;
        temp.V[2] ^= ((Sbox[ (message.V[2] & 0x1 << 2) ^ (message.V[1] >> 62 & 0x3) ] & 0x4) >> 2);
    }
    
    return temp;
}
V256 invSubstitution (V256 message) {
    V256 temp;
    initV256(&temp);
    //Get the identity part of the message
    if(numofboxes <= 21)
    {
        temp.V[0] = (message.V[0] >> 3 * numofboxes);
        temp.V[1] ^= message.V[1];
        temp.V[2] ^= message.V[2];
        temp.V[3] ^= message.V[3];
        //Get the rest through the Sboxes
        for (unsigned i = 1; i <= numofboxes; ++i) {
            temp.V[0] <<= 3;
            temp.V[0] ^= invSbox[ ((message.V[0] >> 3 * (numofboxes - i)) & 0x7) ];
        }
    }
    else if (numofboxes <= 42)
    {
        temp.V[1] ^= (message.V[1] >> (3 * (numofboxes - 22) + 2));
        temp.V[2] ^= message.V[2];
        temp.V[3] ^= message.V[3];
        //Get the rest through the Sboxes
        temp.V[0] ^= invSbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ];
        for (unsigned i = 1; i <= 21; ++i) {
            temp.V[0] <<= 3;
            temp.V[0] ^= invSbox[ ((message.V[0] >> 3 * (21 - i)) & 0x7) ];
        }
        for (unsigned i = 1; i <= (numofboxes - 22); ++i) {
            temp.V[1] <<= 3;
            temp.V[1] ^= invSbox[ ((message.V[1] >> 3 * (numofboxes - 22 - i) + 2) & 0x7) ];
        }
        temp.V[1] <<= 2;
        temp.V[1] ^= ((invSbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ] & 0x6) >> 1);
    }
    else if (numofboxes <= 64)
    {
        temp.V[3] ^= message.V[3];
        //Get the rest through the Sboxes
        temp.V[0] ^= invSbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ];
        for (unsigned i = 1; i <= 21; ++i) {
            temp.V[0] <<= 3;
            temp.V[0] ^= invSbox[ ((message.V[0] >> 3 * (21 - i)) & 0x7) ];
        }
        temp.V[1] ^= invSbox[ (message.V[2] & 0x1 << 2) ^ (message.V[1] >> 62 & 0x3) ];
        for (unsigned i = 1; i <= 20; ++i) {
            temp.V[1] <<= 3;
            temp.V[1] ^= invSbox[ ((message.V[1] >> 3 * (numofboxes - 22 - i) + 2) & 0x7) ];
        }
        temp.V[1] <<= 2;
        temp.V[1] ^= ((invSbox[ (message.V[1] & 0x3 << 1) ^ (message.V[0] >> 63 & 0x1) ] & 0x6) >> 1);
        for (unsigned i = 1; i <= (numofboxes - 43); ++i) {
            temp.V[2] <<= 3;
            temp.V[2] ^= invSbox[ ((message.V[2] >> 3 * (numofboxes - 43 - i) + 1) & 0x7) ];
        }
        temp.V[2] <<= 1;
        temp.V[2] ^= ((invSbox[ (message.V[2] & 0x1 << 2) ^ (message.V[1] >> 62 & 0x3) ] & 0x4) >> 2);
    }
    return temp;
}
V256 encrypt (V256 message) 
{
    V256 c;
    VecAddVecV256(message, roundkeys[0], &c);
    for (unsigned r = 1; r <= rounds; ++r) {
        c = Substitution(c);
        MatMulVecM256(LinMatrices[r-1], c, &c);
        VecAddVecV256(c, roundconstants[r-1], &c);
        VecAddVecV256(c, roundkeys[r], &c);
    }
    return c;
}
V256 decrypt (V256 message) {
    V256 c = message;
    for (unsigned r = rounds; r > 0; --r) {
        VecAddVecV256(c, roundkeys[r], &c);
        VecAddVecV256(c, roundconstants[r-1], &c);
        MatMulVecM256(invLinMatrices[r-1], c, &c);
        c = invSubstitution(c);
    }
    VecAddVecV256(c, roundkeys[0], &c);
    return c;
}

int main()
{
    uint64_t begin;
    uint64_t end;
    uint64_t ans = 0;
    int i;

    printf("\nLowMC WBMatrix Version\n");
    begin = start_rdtsc();
    for (i = 0; i < TEST; i++)
    {
        instantiate_LowMC();
    }
    end = end_rdtsc();
    ans = (end - begin);
    printf("The initialization of LowMC 256-80 cost %llu CPU cycles\n", (ans) / TEST);

    V256 m;
    m.V[1] = 0xFFD5;
    begin = start_rdtsc();
    for (i = 0; i < TEST; i++)
    {
        m = encrypt( m );
    }
    end = end_rdtsc();
    ans = (end - begin);
    printf("The Encryption of LowMC 256-80 cost %llu CPU cycles\n", (ans) / TEST);

    begin = start_rdtsc();
    for (i = 0; i < TEST; i++)
    {
        m = decrypt( m );
    }
    end = end_rdtsc();
    ans = (end - begin);
    printf("The Decryption of LowMC 256-80 cost %llu CPU cycles\n", (ans) / TEST);

    return 0;
}