//#include "treehash03.hpp"

//#include "/home/jbapple/progs/iaca/iaca-lin64/include/iacaMarks.h"

#include <x86intrin.h>

// Two different representations of 128 bits:
typedef uint64_t ui128[2];
typedef __m128i up128;

// Multiply two unsigned 64-bit ints, producing a 128-bit int but returning only the
// high-order bits.
static inline uint64_t mulHi(uint64_t x, uint64_t y) {
  uint64_t lo, hi;
  __asm__("mulq %3" : "=a,a"(lo), "=d,d"(hi) : "%0,0"(x), "r,m"(y));
  return hi;
}

static inline uint64_t mulLo(uint64_t x, uint64_t y) {
  uint64_t lo, hi;
  __asm__("mulq %3" : "=a,a"(lo), "=d,d"(hi) : "%0,0"(x), "r,m"(y));
  return lo;
}

// Universal hashing of x and y
static inline uint64_t deltaHash(const ui128 h, const uint64_t x, const uint64_t y) {
  return x + mulLo(y, h[1]) + mulHi(y, h[0]);
}

// hash x and y into out.
static inline void hashfar(
    uint64_t *out, const ui128 h, const uint64_t x, const uint64_t y) {
  *out = deltaHash(h, x, y);
}

// Just something needed for class template specialization.
template <size_t N>
struct Wrap {
  static const size_t val = N;
};

// Hash
template <typename T, typename WN>
struct CLTriangleHash {
  typedef typename T::Atom Atom;
  static const size_t LTS = WN::val;
  static const size_t TS = 1 << LTS;
  static const size_t HTS = 1 << (LTS - 1);

  static inline void far(Atom out, const up128 *r128, const Atom in[TS]) {

    Atom workspace[HTS - 1];

    T::atom_hashfar(out, *r128, in[0], in[HTS]);
    for (size_t i = 1; i < HTS; ++i) {
      T::atom_hashfar(workspace[i - 1], *r128, in[i], in[i + HTS]);
    }
    CLTriangleHash<T, Wrap<LTS - 1> >::near(out, r128 + 1, workspace);

    // T::atom_hashfar(out, *r128, in[0], in[1]);
    // for (size_t i = 1; i + 2 < TS; i += 2) {
    //   T::atom_hashfar(workspace[i/2], *r128, in[i + 1], in[i + 2]);
    // }
    // CLTriangleHash<T, Wrap<LTS - 1> >::near(out, r128 + 1, workspace);

  }

  static inline void near(Atom out, const up128 *r128, Atom workspace[TS - 1]) {

    T::atom_hashfar(out, *r128, out, workspace[TS - 2]);
    for (size_t i = 0; i < HTS - 1; ++i) {
      T::atom_hashfar(workspace[i], *r128, workspace[i], workspace[i + HTS - 1]);
    }
    CLTriangleHash<T, Wrap<LTS - 1> >::near(out, r128 + 1, workspace);

    // T::atom_hashfar(out, *r128, out, workspace[0]);
    // for (size_t i = 1; i + 2 < TS; i += 2) {
    //   T::atom_hashfar(workspace[i/2], *r128, workspace[i], workspace[i + 1]);
    // }
    // CLTriangleHash<T, Wrap<LTS - 1> >::near(out, r128 + 1, workspace);

  }
};

template <typename T, typename WN>
struct CLTreeTriangleHash {
  typedef typename T::Atom Atom;
  static const size_t LTS = WN::val;
  static const size_t TS = 1 << LTS;

  static inline up128 atomhash(const Atom *data, const size_t length, const up128 *rest,
      const size_t word_length, Atom workspace[65], const up128 *r128, size_t i) {
    Atom *scratch = &workspace[0];
    Atom *accum = &workspace[1];

    while (i + TS <= length) {

      // if (i & TS) {
      //   CLTriangleHash<T, Wrap<LTS> >::far(*scratch, r128, &data[i]);
      //   i += TS;
      //   T::atom_cascade(i >> (LTS - 1), *scratch, &accum[LTS - 1], &r128[LTS - 1]);
      // } else {
      //   CLTriangleHash<T, Wrap<LTS> >::far(accum[LTS - 1], r128, &data[i]);
      //   i += TS;
      // }

      if (!(i & TS)) {
        CLTriangleHash<T, Wrap<LTS> >::far(accum[LTS - 1], r128, &data[i]);
        i += TS;
      } else {
        CLTriangleHash<T, Wrap<LTS> >::far(*scratch, r128, &data[i]);
        i += TS;
        T::atom_cascade(i >> (LTS - 1), *scratch, &accum[LTS - 1], &r128[LTS - 1]);
      }

    }

    return CLTreeTriangleHash<T, Wrap<LTS - 1> >::atomhash(
        data, length, rest, word_length, workspace, r128, i);
  }
};

template <typename T>
struct CLTriangleHash<T, Wrap<1> > {
  typedef typename T::Atom Atom;
  static const size_t LTS = 1;
  static const size_t TS = 2;
  static const size_t HTS = 1;

  static inline void far(Atom out, const up128 *r128, const Atom in[2]) {
    T::atom_hashfar(out, *r128, in[0], in[1]);
  }

  static inline void near(Atom out, const up128 *r128, Atom workspace[1]) {
    T::atom_hashnear(out, *r128, workspace[0]);
  }
};

template <typename T>
struct CLTreeTriangleHash<T, Wrap<0> > {
  typedef typename T::Atom Atom;
  static const size_t LTS = 0;
  static const size_t TS = 1;

  static inline up128 atomhash(const Atom *data, const size_t length, const up128 *rest,
      const size_t word_length, Atom workspace[65], const up128 *r128, size_t) {
    int place = __builtin_ctzll(length);
    const Atom *latest = place ? &workspace[place] : &data[length - 1];
    for (; (length >> place) > 2; ++place) {
      if ((length >> place) & 2) {
        T::atom_hashnear(workspace[place + 1], r128[place + 1], *latest);
        latest = &workspace[place + 1];
      }
    }
    return T::atom_hashdown(&r128[place + 1], *latest, rest,
        word_length - T::ATOM_SIZE * length, word_length);
  }
};

static inline up128 clhashfar(const up128 r128, const up128 x, const up128 y) {
  up128 result = _mm_xor_si128(x, r128);
  result = _mm_clmulepi64_si128(result, result, 1);
  result = _mm_xor_si128(result, y);
  return result;
}

static inline void clhashnear(const up128 r128, up128 * x, const up128 y) {
  *x = _mm_xor_si128(*x, r128);
  *x = _mm_clmulepi64_si128(*x, *x, 1);
  *x = _mm_xor_si128(*x, y);
}

template <size_t ATOM_SIZE_PARAM>
struct TreeHash {
  static const size_t ATOM_SIZE = ATOM_SIZE_PARAM;

  typedef up128 Atom[ATOM_SIZE];

  static inline void atom_hashfar(
      Atom out, const up128 r128, const Atom x, const Atom y) {
    for (size_t i = 0; i < ATOM_SIZE; ++i) {
      out[i] = clhashfar(r128, x[i], y[i]);
    }
  }

  static inline void atom_hashnear(Atom out, const up128 r128, const Atom y) {
    for (size_t i = 0; i < ATOM_SIZE; ++i) {
      clhashnear(r128, &out[i], y[i]);
    }
  }

  static inline void atom_cascade(
      size_t i, Atom scratch, Atom accum[64], const up128 *__restrict__ r128) {
    int j = 0;
    for (; j + 1 < __builtin_ctzll(i / 2); ++j) {
      atom_hashnear(scratch, r128[j + 1], accum[j]);
    }
    atom_hashfar(accum[j + 1], r128[j + 1], scratch, accum[j]);
  }

  static inline up128 atom_hashdown(const up128 *r128, const Atom hashed,
      const up128 *rest, const size_t rest_length, const size_t) {
    Atom tmp;
    for (size_t i = 0; i < ATOM_SIZE; ++i) {
      if (i < rest_length) {
        tmp[i] = clhashfar(r128[0], hashed[i], rest[i]);
      } else {
        tmp[i] = hashed[i];
      }
    }
    for (size_t current = ATOM_SIZE; current > 1; current = current / 2 + (current & 1)) {
      ++r128;
      for (size_t i = current & 1; i < current / 2; ++i) {
        clhashnear(*r128, &tmp[i], tmp[i + current / 2]);
      }
    }
    return tmp[0];
  }

  template <size_t LOG_TRIANGLE_SIZE>
  static inline up128 atomhash(const uint64_t *r64, const Atom *data, const size_t length,
      const up128 *rest, const size_t word_length) {
    Atom workspace[65];
    const up128 *r128 = reinterpret_cast<const up128 *>(r64);
    size_t i = 0;

    return CLTreeTriangleHash<TreeHash, Wrap<LOG_TRIANGLE_SIZE> >::atomhash(
        data, length, rest, word_length, workspace, r128, i);
  }

  template <size_t LOG_TRIANGLE_SIZE>
  static uint64_t hash(const void *rvoid, const uint64_t *data, const size_t length) {
    const up128 *data128 = reinterpret_cast<const up128 *>(data);
    const size_t length128 = length / 2;
    const uint64_t *r64 = reinterpret_cast<const uint64_t *>(rvoid);
    up128 result = atomhash<LOG_TRIANGLE_SIZE>(r64 + 4,
        reinterpret_cast<const Atom *>(data128), length128 / ATOM_SIZE,
        data128 + ATOM_SIZE * (length128 / ATOM_SIZE), length128);
    int64_t d64[2] = {_mm_extract_epi64(result, 0), _mm_extract_epi64(result, 1)};
    uint64_t *u64 = reinterpret_cast<uint64_t *>(d64);
    u64[0] = deltaHash(r64, u64[0], u64[1]);
    const uint64_t rem = (length & 1) ? deltaHash(r64, data[length - 1], length) : length;
    return u64[0] * r64[2] + rem * r64[3] + mulHi(r64[2], rem);
  }
};
