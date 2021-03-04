#include "WBMatrix/WBMatrix.h"

#ifndef LOWMC_H
#define LOWMC_H

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

void keyschedule();
void set_key(V256 k);
void instantiate_LowMC();
V256 Substitution(V256 message);
V256 invSubstitution(V256 message);
V256 encrypt(V256 message);
V256 decrypt(V256 message);
#endif