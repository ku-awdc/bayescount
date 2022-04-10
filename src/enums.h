#ifndef ENUMS_H
#define ENUMS_H

enum class dists { nbinom, mvnbinom, lnormpois, mvlnormpois, poisson };
enum class containers { stdarray, stdvector, rcppvector };
enum class ktypes { fix, mm, ql, ml };
enum class optswitch { never, sometimes, always};


// TODO: remove:
enum class analyses { bnb_estcv, all_estcv, bnb_fixcv, all_fixcv };


#endif // ENUMS_H
