/////////////////////////////////////
// This C code is a companion to the paper
//
// Reference: Owen Kaser and Daniel Lemire, Strongly universal string hashing is fast,
// Computer Journal
// http://arxiv.org/abs/1202.4961
//
// It shows that we can compute strongly universal hash functions very quickly.
/////////////////////////////////////

//
// this code will hash strings of 32-bit characters. To use on
// strings of 8-bit characters, you may need some adequate padding.
//
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <cassert>
#include <unistd.h>
#ifdef __AVX__
#define __PCLMUL__ 1
#endif

#include <iostream>
#include <iomanip>

typedef unsigned long long ticks;

// Taken from stackoverflow (see
// http://stackoverflow.com/questions/3830883/cpu-cycle-count-based-profiling-in-c-c-linux-x86-64)
// Can give nonsensical results on multi-core AMD processors.
ticks rdtsc() {
  unsigned int lo, hi;
  asm volatile(
      "cpuid \n" /* serializing */
      "rdtsc"
      : "=a"(lo), "=d"(hi) /* outputs */
      : "a"(0)             /* inputs */
      : "%ebx", "%ecx");
  /* clobbers*/
  return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}

ticks startRDTSC(void) { return rdtsc(); }

ticks stopRDTSCP(void) { return rdtsc(); }
// start and stop are as recommended by
// Gabriele Paoloni, How to Benchmark Code Execution Times on Intel IA-32 and IA-64
// Instruction Set Architectures
// September 2010
// http://edc.intel.com/Link.aspx?id=3954

/*static __inline__ ticks fancystartRDTSC(void) {
        unsigned cycles_low, cycles_high;
        asm volatile ("CPUID\n\t"
                        "RDTSC\n\t"
                        "mov %%edx, %0\n\t"
                        "mov %%eax, %1\n\t": "=r" (cycles_high), "=r" (cycles_low)::
                        "%rax", "%rbx", "%rcx", "%rdx");
        return ((ticks) cycles_high << 32) | cycles_low;
}

static __inline__ ticks fancystopRDTSCP(void) {
        unsigned cycles_low, cycles_high;
/// This should work fine on most machines, if the RDTSCP thing
/// fails for you, use the  rdtsc() call instead.
        asm volatile("RDTSCP\n\t"
                        "mov %%edx, %0\n\t"
                        "mov %%eax, %1\n\t"
                        "CPUID\n\t": "=r" (cycles_high), "=r" (cycles_low):: "%rax",
                        "%rbx", "%rcx", "%rdx");
        return ((ticks) cycles_high << 32) | cycles_low;
}*/

extern "C" {
#include "hashfunctions32bits.h"
#include "hashfunctions64bits.h"
}

#include <string>

template <typename T>
struct Named {
  T val;
  ::std::string name;
};

template <typename T>
Named<T> named(const T &x, const ::std::string &name) {
  return Named<T>{x, name};
}

#define NAME(x) named(x, #x)

#ifdef __PCLMUL__

extern "C" {
#include "clmulhashfunctions32bits.h"
#include "clmulhashfunctions64bits.h"
#include "clmulpoly64bits.h"
#include "ghash.h"
#include "clmulhierarchical64bits.h"
}

#include "treehash/treehash04.hpp"

const Named<hashFunction64> funcArr64[] = {
  NAME(&hashCity), NAME(&hashVHASH64),
  NAME(&CLHASH), NAME(&hashGaloisFieldfast64_precomp_unroll),
  NAME(&hashGaloisFieldfast64halfunrolled_precomp), NAME(&hashSipHash),
  NAME(&GHASH64bit),
  NAME(&TreeHash<8>::hash<2>)
};

const Named<hashFunction> funcArr[] = {NAME(&hashGaloisFieldMultilinear),
    NAME(&hashGaloisFieldMultilinearHalfMultiplications), NAME(&hashMultilinear),
    NAME(&hashMultilinear2by2), NAME(&hashMultilinearhalf), NAME(&hashMultilineardouble),
    NAME(&hashNH), NAME(&hashRabinKarp), NAME(&hashFNV1), NAME(&hashFNV1a),
    NAME(&hashSAX), NAME(&pyramidal_Multilinear)};

#else

const Named<hashFunction> funcArr[] = {NAME(&hashMultilinear), NAME(&hashMultilinear2by2),
    NAME(&hashMultilinearhalf), NAME(&hashMultilineardouble), NAME(&hashNH),
    NAME(&hashRabinKarp), NAME(&hashFNV1), NAME(&hashFNV1a), NAME(&hashSAX),
    NAME(&pyramidal_Multilinear)};

const Named<hashFunction64> funcArr64[] = {
    NAME(&hashCity), NAME(&hashMMH_NonPyramidal), NAME(&hashNH64)};

#endif

const int HowManyFunctions64 = sizeof(funcArr64) / sizeof(funcArr64[0]);
const int HowManyFunctions = sizeof(funcArr) / sizeof(funcArr[0]);

void force_computation(uint32_t forcedValue) {
  // make sure forcedValue has to be computed, but avoid output (unless unlucky)
  if (forcedValue % 277387 == 17) printf("wow, what a coincidence! (in benchmark.c)");
  // printf("# ignore this #%d\n", forcedValue);
}

void printusage(char *command) { printf(" Usage: %s -b (32|64)", command); }

int main(int argc, char **arg) {
  int N = 1 << 6;  // should be divisible by two!
  int SHORTTRIALS = 1 << 10;
  int HowManyRepeats = 1;
  int bit = 64;
  int i, k, j;
  int elapsed;
  Named<hashFunction> thisfunc;
  ticks bef, aft;
  struct timeval start, finish;
  uint64_t *randbuffer;
  int err = posix_memalign(reinterpret_cast<void **>(&randbuffer), 16,
      sizeof(uint64_t) * (N + 3));  //, ] __attribute__ ((aligned (16)));
  if (err) exit(err);
  uint32_t sumToFoolCompiler = 0;
  uint32_t *intstring;  //[N] __attribute__ ((aligned (16))); // // could force 16-byte
  // alignment with  __attribute__ ((aligned (16)));
  err = posix_memalign(reinterpret_cast<void **>(&intstring), 16, sizeof(uint32_t) * N);
  if (err) exit(err);
  int c;
  while ((c = getopt(argc, arg, "hb:")) != -1) {
    switch (c) {
      case 'h':
        printusage(arg[0]);
        return 0;
      case 'b':
        bit = atoi(optarg);
        if ((bit != 32) && (bit != 64)) {
          printusage(arg[0]);
          return -1;
        }
        break;
      default:
        abort();
    }
  }
  for (i = 0; i < N + 3; ++i) {
    randbuffer[i] = rand() | ((uint64_t)(rand()) << 32);
  }
  for (i = 0; i < N; ++i) {
    intstring[i] = rand();
  }
  // printf(
  //     "For documentation, see Strongly universal string hashing is fast at "
  //     "http://arxiv.org/abs/1202.4961 \n");
  // printf(
  //     "Reporting the number of cycles per byte and the billions of bytes processed per "
  //     "second.\n");
  size_t max_function_name_length = 0;
  for (i = 0; i < HowManyFunctions64; ++i) {
    max_function_name_length =
      ::std::max(funcArr64[i].name.size(), max_function_name_length);
  }
  for (i = 0; i < HowManyFunctions; ++i) {
    max_function_name_length =
      ::std::max(funcArr[i].name.size(), max_function_name_length);
  }
  for (k = 0; k < HowManyRepeats; ++k) {
    if (bit == 64) {
      // printf("test #%d (64-bit hash values) ", k + 1);
      // printf("(%d bytes) \n", N * 4);
      // ::std::cout << "Bytes per cycle" << std::endl;
      Named<hashFunction64> thisfunc64;
      for (i = 0; i < HowManyFunctions64; ++i) {
        sumToFoolCompiler = 0;
        thisfunc64 = funcArr64[i];
        ::std::cout << ::std::setw(1 + max_function_name_length) << ::std::setfill(' ')
                    << ::std::left << thisfunc64.name;
        fflush(stdout);
        sumToFoolCompiler += thisfunc64.val(&randbuffer[0], (uint64_t *)&intstring[0],
            N / 2);  // we do not count the first run
        gettimeofday(&start, 0);

        assert(N / 2 * 2 == N);
        ticks best = -1;
        for (j = 0; j < SHORTTRIALS; ++j) {
          bef = startRDTSC();
          sumToFoolCompiler +=
              thisfunc64.val(&randbuffer[0], (uint64_t *)&intstring[0], N / 2);
          aft = stopRDTSCP();
          const ticks diff = aft - bef;
          if (diff < best) best = diff;
        }
        gettimeofday(&finish, 0);
        elapsed =
            (1000000 * (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec));
        /*
          printf(
            "CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
            (4.0 * SHORTTRIALS * N) / (1000. * elapsed));
        */
        // printf("Bytes per CPU cycle: %f\n",
        //        (8.0*N/2)/(best * 1.0));
        printf("%f\n", (8.0 * N / 2) / (best * 1.0));
        force_computation(sumToFoolCompiler);
      }
    } else {
      printf("test #%d (32-bit hash values)\n", k + 1);
      size_t max_function_name_length = 0;
      for (i = 0; i < HowManyFunctions64; ++i) {
        max_function_name_length =
            ::std::max(funcArr[i].name.size(), max_function_name_length);
      }

      for (i = 0; i < HowManyFunctions; ++i) {
        sumToFoolCompiler = 0;
        thisfunc = funcArr[i];
        ::std::cout << ::std::setw(1 + max_function_name_length) << ::std::setfill(' ')
                    << ::std::left << thisfunc.name;
        fflush(stdout);
        sumToFoolCompiler += thisfunc.val(&randbuffer[0], &intstring[0],
            N);  // we do not count the first pass
        gettimeofday(&start, 0);
        bef = startRDTSC();
        for (j = 0; j < SHORTTRIALS; ++j)
          sumToFoolCompiler += thisfunc.val(&randbuffer[0], &intstring[0], N);
        aft = stopRDTSCP();
        gettimeofday(&finish, 0);
        elapsed =
            (1000000 * (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec));
        printf("CPU cycle/byte = %f \t billions of bytes per second =  %f    \n",
            (aft - bef) * 1.0 / (4.0 * SHORTTRIALS * N),
            (4.0 * SHORTTRIALS * N) / (1000. * elapsed));
        force_computation(sumToFoolCompiler);
      }
    }
    // printf("\n");
  }
  force_computation(
      sumToFoolCompiler);  // printf("# ignore this #%d\n", sumToFoolCompiler);
}
