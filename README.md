# WBMatrix based LowMC

An implementation of [LowMC block cipher](https://eprint.iacr.org/2016/687).

Supports for the block size $n=256$, the key
size $k=256$, and the number of Sboxes $m\in[1,21]$. The number of rounds $r$ needed to reach the security can be referenced to the source code of LowMC. $m$ and $r$ can be modified in [lowmc.h](https://github.com/Nexus-TYF/WBMatrix-based-LowMC/blob/main/include/lowmc.h).

The exampled parameter sets for LowMC instantiations are listed in the following table. $d$ is the allowed data complexity of attacks.

|     | $n$  | $m$  | $k$  | $d$  | $r$  |
|  ----  | ----  | ----  | ----  | ----  | ----  |
| LowMC v2  | 256 | 12  | 256  | 128  | 35  |
| LowMC v3  | 256 | 12  | 256  | 128  | 38  |

## Build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## Run

```
$ ./LMC
```

## Included libraries
1. [WBMatrix](https://github.com/Nexus-TYF/WBMatrix)<br>
2. [LowMC](https://github.com/LowMC/lowmc)<br>
