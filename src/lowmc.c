#include "lowmc.h"

const uint8_t Sbox[] =
        {0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02};
const uint8_t invSbox[] =
        {0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04};
        
void keyschedule() 
{
    for (unsigned r = 0; r <= rounds; ++r) 
    {
        MatMulVecM256(KeyMatrices[r], key, &roundkeys[r]);
    }
    return;
}
void set_key(V256 k) 
{
    key.V[0] = k.V[0];
    key.V[1] = k.V[1];
    key.V[2] = k.V[2];
    key.V[3] = k.V[3];
    keyschedule();
}
void instantiate_LowMC(V256 k)
{
    // Create LinMatrices and invLinMatrices
    for (unsigned r = 0; r < rounds; ++r) 
    {
        // Create matrix
        genMatpairM256(&LinMatrices[r], &invLinMatrices[r]);

        // Create roundconstants
        randV256(&roundconstants[r]);
    }

    // Create KeyMatrices
    for (unsigned r = 0; r <= rounds; ++r) 
    {
        // Create matrix
        do {
            randM256(&KeyMatrices[r]);
        // Repeat if matrix is not of maximal rank
        } while(isinvertM256(KeyMatrices[r]));
    }
    // key schedule
    set_key(k);
    
    return;
}
V256 Substitution(V256 message) 
{
    V256 temp;
    initV256(&temp);
    //Get the identity part of the message
    if(numofboxes <= 21)
    {
        temp.V[0] = (message.V[0] >> 3 * numofboxes);
        temp.V[1] = message.V[1];
        temp.V[2] = message.V[2];
        temp.V[3] = message.V[3];
        //Get the rest through the Sboxes
        for (unsigned i = 1; i <= numofboxes; ++i) 
        {
            temp.V[0] <<= 3;
            temp.V[0] ^= Sbox[ ((message.V[0] >> 3 * (numofboxes - i)) & 0x7) ];
        }
    }
    else temp = message;
    return temp;
}
V256 invSubstitution(V256 message)
{
    V256 temp;
    initV256(&temp);
    //Get the identity part of the message
    if(numofboxes <= 21)
    {
        temp.V[0] = (message.V[0] >> 3 * numofboxes);
        temp.V[1] = message.V[1];
        temp.V[2] = message.V[2];
        temp.V[3] = message.V[3];
        //Get the rest through the Sboxes
        for (unsigned i = 1; i <= numofboxes; ++i) 
        {
            temp.V[0] <<= 3;
            temp.V[0] ^= invSbox[ ((message.V[0] >> 3 * (numofboxes - i)) & 0x7) ];
        }
    }
    else temp = message;
    return temp;
}
V256 encrypt(V256 message) 
{
    V256 c;
    VecAddVecV256(message, roundkeys[0], &c);
    for (unsigned r = 1; r <= rounds; ++r) 
    {
        c = Substitution(c);
        MatMulVecM256(LinMatrices[r-1], c, &c);
        VecAddVecV256(c, roundconstants[r-1], &c);
        VecAddVecV256(c, roundkeys[r], &c);
    }
    return c;
}
V256 decrypt(V256 message)
{
    V256 c = message;
    for (unsigned r = rounds; r > 0; --r) 
    {
        VecAddVecV256(c, roundkeys[r], &c);
        VecAddVecV256(c, roundconstants[r-1], &c);
        MatMulVecM256(invLinMatrices[r-1], c, &c);
        c = invSubstitution(c);
    }
    VecAddVecV256(c, roundkeys[0], &c);
    return c;
}