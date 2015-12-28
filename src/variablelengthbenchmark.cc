/////////////////////////////////////
// This C code is a companion to the paper
//
/////////////////////////////////////

//
// this code will hash strings of 64-bit characters. To use on
// strings of 8-bit characters, you may need some adequate padding.
//
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <cassert>

#ifdef __AVX__
#define __PCLMUL__ 1
#endif

typedef unsigned long long ticks;

// Taken from stackoverflow (see http://stackoverflow.com/questions/3830883/cpu-cycle-count-based-profiling-in-c-c-linux-x86-64)
// Can give nonsensical results on multi-core AMD processors.
ticks rdtsc() {
    unsigned int lo, hi;
    asm volatile (
        "cpuid \n" /* serializing */
        "rdtsc"
        : "=a"(lo), "=d"(hi) /* outputs */
        : "a"(0) /* inputs */
        : "%ebx", "%ecx");
    /* clobbers*/
    return ((unsigned long long) lo) | (((unsigned long long) hi) << 32);
}

ticks startRDTSC(void) {
    return rdtsc();
}

ticks stopRDTSCP(void) {
    return rdtsc();
}
// start and stop are as recommended by
// Gabriele Paoloni, How to Benchmark Code Execution Times on Intel IA-32 and IA-64 Instruction Set Architectures
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

#include "clmulhashfunctions32bits.h"
#include "clmulhashfunctions64bits.h"
#include "clmulpoly64bits.h"
#include "clmulhierarchical64bits.h"
#include "ghash.h"
  //#include "bigendianuniversal.h"
}

#include "treehash/treehash04.hpp"

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

const Named<hashFunction64> funcArr64[] =
  {
    NAME(&CLHASH),
    NAME(&TreeHash<8>::hash<2>),
    NAME(&TreeHash<11>::hash<2>),
  };

#include <iostream>

using namespace std;

template <hashFunction64 F>
struct H {
  const string name;
};

template <hashFunction64 F>
inline void one(const uint64_t *randbuffer, uint32_t *sumToFoolCompiler, size_t n,
    uint64_t *intstring, const H<F> &) {
  static const size_t SHORTTRIALS = 1 << 16;
  *sumToFoolCompiler += F(randbuffer, intstring, n);  // we do not count the first one
  ticks lowest = ~(ticks)0;
  for (size_t j = 0; j < SHORTTRIALS; ++j) {
    const ticks bef = startRDTSC();
    *sumToFoolCompiler += F(randbuffer, intstring, n);
    const ticks aft = stopRDTSCP();
    const ticks diff = aft - bef;
    lowest = (lowest < diff) ? lowest : diff;
  }
  printf(" %f ", (8.0 * n) / (lowest * 1.0));
  fflush(stdout);
}

inline void many(int, const uint64_t *, uint32_t *, size_t, uint64_t *intstring) {
  printf("\n");
  fflush(stdout);
  free(intstring);
}

template <typename F, typename... FS>
inline void many(int selector, const uint64_t *randbuffer, uint32_t *sumToFoolCompiler,
    size_t n, uint64_t *intstring, const F &f, const FS &... fs) {
  if (selector & 1) {
    one(randbuffer, sumToFoolCompiler, n, intstring, f);
  }
  selector >>= 1;
  many(selector, randbuffer, sumToFoolCompiler, n, intstring, fs...);
}

template <typename F>
inline void printHashName(const F& f) {
  cout << f.name << ' ';
}


inline void printHashNames(int) {
  cout << endl;
}

template <typename F, typename... FS>
inline void printHashNames(int selector, const F& f, const FS &... fs) {
  if (selector & 1) printHashName(f);
  selector >>= 1;
  printHashNames(selector, fs...);
}

template <typename... FS>
inline void help(bool print, int selector, uint64_t *randbuffer, uint32_t *sumToFoolCompiler,
    size_t n, const FS &... fs) {
  if (print) {
    cout << 0 << ' ';
    printHashNames(selector, fs...);
  }
  uint64_t *intstring;
  int err =
      posix_memalign(reinterpret_cast<void **>(&intstring), 16, sizeof(uint64_t) * n);
  if (err) exit(err);
  for (size_t i = 0; i < n; ++i) {
    intstring[i] = rand() | ((uint64_t)(rand()) << 32);
  }
  cout << n << ' ';
  fflush(stdout);
  many(selector, randbuffer, sumToFoolCompiler, n, intstring, fs...);
}

/*&hashVHASH64
                                                 ,&hashCity
                                                 , &hashSipHash
                                                 ,&GHASH64bit

                                                 ,&treeCL9
                                                 ,&treeCTZ8
                                                 ,&treeCTZquad13
                                                 ,&treeCTZoct13*/

//  &CLHASH,
  // &TreeHash04<1,  Wrap<3> >::hash,
  //&TreeHash04<4,  Wrap<3> >::hash,
  //&TreeHash04<8,  Wrap<3> >::hash,
  //&TreeHash04<16,  Wrap<3> >::hash,
  //&TreeHash04<32,  Wrap<3> >::hash,
  //&TreeHash04<64,  Wrap<3> >::hash,
  //&TreeHash04<2,  Wrap<1> >::hash,
  //&TreeHash04<3,  Wrap<1> >::hash,
  //&TreeHash04<5,  Wrap<3> >::hash,
  //&TreeHash04<7,  Wrap<7> >::hash,
    //&TreeHash04<11,  Wrap<2> >::hash,
    // &TreeHash04<11,  Wrap<2> >::hash,
    // &TreeHash04<11,  Wrap<3> >::hash,
    // &TreeHash04<11,  Wrap<4> >::hash,
    // &TreeHash04<11,  Wrap<5> >::hash,
    // &TreeHash04<11,  Wrap<6> >::hash,
    // &TreeHash04<11,  Wrap<7> >::hash,
    // &TreeHash04<11,  Wrap<8> >::hash,
    // &TreeHash04<11,  Wrap<9> >::hash,
    // &TreeHash04<7,  Wrap<3> >::hash,
    // &TreeHash04<8,  Wrap<3> >::hash,
    // &TreeHash04<9,  Wrap<3> >::hash,
    // &TreeHash04<10,  Wrap<3> >::hash,
    // &TreeHash04<11,  Wrap<3> >::hash,
   // &TreeHash04<10, Wrap<3> >::hash,
  // &TreeHash<11>::hash<2>,
  //  &TreeHash<2>::hash<3>,
  //  &TreeHash<8>::hash<2>,
  // &TreeHash<11>::hash<1>,
  // &TreeHash<2>::hash<1>,
  //&TreeHash<11>::hash<0>,
  //&TreeHash<2>::hash<0>,
   // &TreeHash04<4,  Wrap<3> >::hash,
   // &TreeHash04<7,  Wrap<3> >::hash,
   // &TreeHash04<11, Wrap<3> >::hash,
  //&treeCTZdoct12,
  //&piecewiseTreeHash,
  //&piecewiseNearHash,

  //&treeMalloc, &treeAlloca
  //&hashAlloc
//};

// const char* functionnames64[HowManyFunctions64] = { /*"64-bit VHASH        "
//                                                     , "Google's City       ", "SipHash             "
//                                                     ,"GHASH          "
//                                                     ,"treeCL9"
//                                                     ,"treeCTZ8"
//                                                     ,"treeCTZquad13"
//                                                     ,"treeCTZoct13"*/
//   //"hashAlloc",
//   //"treeCTZdoct12",
//   //  "piecewiseTreeHash",
//   //"piecewiseNearHash",
//     "CLHASH",
//   //"treeMalloc", "treeAlloca"
//   // "&TreeHash04<2,  Wrap<3> >::hash",
//   // "&TreeHash04<4,  Wrap<3> >::hash",
//   // "&TreeHash04<8,  Wrap<3> >::hash",
//   // "&TreeHash04<16,  Wrap<3> >::hash",
//   // "&TreeHash04<32,  Wrap<3> >::hash",
//   // "&TreeHash04<64,  Wrap<3> >::hash",
//   //"&TreeHash04<2,  Wrap<1> >::hash",

//   //"&TreeHash04<5, Wrap<3> >::hash",
//   //"&TreeHash04<, Wrap<7> >::hash",
//   "&TreeHash<11>::hash<2>",
//        "&TreeHash<2>::hash<3>",
//     "&TreeHash<8>::hash<2>",
//   //   "&TreeHash04<2>::hash<1>",
//     //"&TreeHash04<11>::hash<0>",
//     //"&TreeHash04<2>::hash<0>",

//   // "&TreeHash04<11,  Wrap<2> >::hash",
//   // "&TreeHash04<11,  Wrap<3> >::hash",
//   // "&TreeHash04<11,  Wrap<4> >::hash",
//   // "&TreeHash04<11,  Wrap<5> >::hash",
//   // "&TreeHash04<11,  Wrap<6> >::hash",
//   // "&TreeHash04<11,  Wrap<7> >::hash",
//   // "&TreeHash04<11,  Wrap<8> >::hash",
//   // "&TreeHash04<11,  Wrap<9> >::hash",
//   // "&TreeHash04<7, Wrap<3> >::hash",
//   // "&TreeHash04<8, Wrap<3> >::hash",
//   // "&TreeHash04<9,  Wrap<3> >::hash",
//   // "&TreeHash04<10,  Wrap<3> >::hash",
//   // "&TreeHash04<11, Wrap<3> >::hash",
// };

const int HowManyFunctions64 = sizeof(funcArr64)/sizeof(funcArr64[0]);

int main(int c, char ** arg) {
    (void) (c);
    (void) (arg);
    const int N = 2048; // should be divisible by two!
    int which_algos = 0xffffffff;
    assert(HowManyFunctions64 <= 32);
    if (c > 1)
        which_algos = atoi(arg[1]);  // bitmask
    int lengthStart = 1, lengthEnd = N, lengthIncrement = 1; // inclusive
    if (c > 2)
        lengthStart = atoi(arg[2]);
    if (c > 3)
        lengthEnd = atoi(arg[3]);
    if (c > 4)
        lengthIncrement = atoi(arg[4]);
    assert((lengthEnd & 1) == 0);

    //int i, j;
    //int length;
    //int SHORTTRIALS;
    //struct timeval start, finish;
    uint64_t randbuffer[150] __attribute__ ((aligned (16)));// 150 should be plenty
    uint32_t sumToFoolCompiler = 0;
    //uint64_t * intstring = (uint64_t *)malloc(lengthEnd * sizeof(uint64_t)); // // could force 16-byte alignment with  __attribute__ ((aligned (16)));

    for (int i = 0; i < 150; ++i) {
        randbuffer[i] = rand() | ((uint64_t)(rand()) << 32);
    }
    //printf("#Reporting the number of bytes per cycle.\n");
    //printf("#First number is input length in  8-byte words.\n");
    // printf("0          ");
    // for (int i = 0; i < HowManyFunctions64; ++i) {
    //     if (which_algos & (0x1 << i))
    //       printf("\"%s\" ", funcArr64[i].name.c_str());
    // }
    printf("\n");
    fflush(stdout);
    //uint64_t * intstring = NULL;
    for (int length = lengthStart; length <= lengthEnd; length += (c > 4) ? lengthIncrement : (1 + (length >> 7))) {
      //uint64_t * randbuffer, uint32_t * sumToFoolCompiler,
      //size_t n, const FS&& ... fs
#define HF(f) H< f > { #f }
      help(length == lengthStart, which_algos, randbuffer, &sumToFoolCompiler, length, 
           HF(CLHASH),
           HF(TreeHash<11>::hash<3>),
           HF(TreeHash<18>::hash<2>),
           HF(TreeHash<18>::hash<3>)
           // H<TreeHash<11>::hash<2> >(),
           // H<TreeHash<2>::hash<3> >()
           );
#undef HF
      /*
      if (intstring) free(intstring);
      int err = posix_memalign(reinterpret_cast<void **>(&intstring), 16,
                               sizeof(uint64_t) * length);  //, ] __attribute__ ((aligned (16)));
      if (err) exit(err);
      for (i = 0; i < length; ++i) {
        intstring[i] = rand() | ((uint64_t)(rand()) << 32);
      }

      //SHORTTRIALS = 80000000 / length;
      //SHORTTRIALS = (SHORTTRIALS < 1024) ? 1024 : SHORTTRIALS;
        //SHORTTRIALS = (SHORTTRIALS > (1 << 16)) ? (1 << 16) : SHORTTRIALS;
      SHORTTRIALS = 1 << 10;
        printf("%8d \t\t", length);
        fflush(stdout);
        Named<hashFunction64> thisfunc64;
        for (i = 0; i < HowManyFunctions64; ++i) {
          //if (!(which_algos & (0x1 << i)))
          //    continue;  // skip unselected algos
            thisfunc64 = funcArr64[i];
            sumToFoolCompiler += thisfunc64.val(randbuffer, intstring, length); // we do not count the first one
            gettimeofday(&start, 0);
            ticks lowest = ~(ticks)0;
            for (j = 0; j < SHORTTRIALS; ++j) {
                const ticks bef = startRDTSC();
                sumToFoolCompiler += thisfunc64.val(randbuffer, intstring, length);
                const ticks aft = stopRDTSCP();
                const ticks diff = aft-bef;
                lowest = (lowest < diff) ? lowest : diff;
            }
            gettimeofday(&finish, 0);
            printf(" %f ", (8.0 * length)/ (lowest * 1.0));
            fflush(stdout);
        }
        printf("\n");
      */
    }
    //free(intstring);
      
    printf("# ignore this #%d\n", sumToFoolCompiler);

}

