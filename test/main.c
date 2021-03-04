#include "lowmc.h"
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

int main()
{
    uint64_t begin;
    uint64_t end;
    uint64_t ans = 0;
    int i;

    printf("LowMC WBMatrix Version\n");
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