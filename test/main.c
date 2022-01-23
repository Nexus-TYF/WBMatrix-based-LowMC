#include "lowmc.h"

int main()
{
    int i;
    V256 ptx, ctx, k;
    initV256(&ptx);

    k.V[0] = 0x0123456789abcdef;
    instantiate_LowMC(k);

    printf("plaintext:\n");
    printV256(ptx);
    ctx = encrypt(ptx);
    printf("ciphertext:\n");
    printV256(ctx);

    ptx = decrypt(ctx);
    printf("decrypted text:\n");
    printV256(ptx);

    return 0;
}