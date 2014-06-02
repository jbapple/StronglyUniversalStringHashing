#ifndef _CLMUL_H_
#define _CLMUL_H_


#ifdef __PCLMUL__

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <x86intrin.h>
//#include <wmmintrin.h>


/////////////////////////////////////////////////////////////////
// working from 
// "Modular Reduction in GF(2n) without Pre-computational Phase"
// by M. Knezevic, K. Sakiyama, J. Fan, and I. Verbauwhede (2008)
// algo. 4
//
// A is input, M is irred. (n=32)
//
//(((A div x^n) * M ) div x^n) * M) mod x^n
//+(A mod x^n)
//////////////////////////////////////////////////
uint32_t barrettWithoutPrecomputation32( __m128i A) {
    ///http://www.jjj.de/mathdata/minweight-primpoly.txt
    const uint64_t irredpoly = 1UL+(1UL<<2)+(1UL<<6)+(1UL<<7)+(1UL<<32);
    // it is important, for the algo. we have chosen that 7 is smaller 
    // equal than 16=32/2
    const int n = 32;// degree of the polynomial
    const __m128i C = _mm_set_epi64x(0,irredpoly);// C is the irreducible poly.
    /////////////////
    /// This algo. requires two multiplications (_mm_clmulepi64_si128)
    /// They are probably the bottleneck.
    /// Note: Barrett's original algorithm also requires two multiplications.
    ////////////////
    const __m128i Q1 = _mm_srli_epi64 (A, n);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_epi64 (Q2, n);
    // commenting out the long way derived from the paper (following two lines are enough)
    //__m128i R1 = _mm_and_si128 (maskm128,A);
    //__m128i R2 = _mm_and_si128 (maskm128,_mm_clmulepi64_si128( Q3, C, 0x00));
    //__m128i final  = _mm_xor_si128 (R1, R2);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);        
    return _mm_cvtsi128_si32(final); 
}
uint32_t barrettWithoutPrecomputation64( __m128i A) {
    ///http://www.jjj.de/mathdata/minweight-primpoly.txt
    // it is important, for the algo. we have chosen that 4 is smaller 
    // equal than 32=64/2
    const int n = 64;// degree of the polynomial
    const __m128i C = _mm_set_epi64x(1U,(1U<<4)+(1U<<3)+(1U<<1)+(1U<<0));// C is the irreducible poly. (64,4,3,1,0)
    /////////////////
    /// This algo. requires two multiplications (_mm_clmulepi64_si128)
    /// They are probably the bottleneck.
    /// Note: Barrett's original algorithm also required two multiplications.
    ////////////////
    assert(n/8==8);
    const __m128i Q1 = _mm_srli_si128 (A, 8);
    const __m128i Q2 = _mm_clmulepi64_si128( Q1, C, 0x00);// A div x^n
    const __m128i Q3 = _mm_srli_si128 (Q2, 8);
    // commenting out the long way derived from the paper (following two lines are enough)
    //__m128i R1 = _mm_and_si128 (maskm128,A);
    //__m128i R2 = _mm_and_si128 (maskm128,_mm_clmulepi64_si128( Q3, C, 0x00));
    //__m128i final  = _mm_xor_si128 (R1, R2);
    const __m128i Q4 = _mm_clmulepi64_si128( Q3, C, 0x00);
    const __m128i final  = _mm_xor_si128 (A, Q4);        
    return _mm_cvtsi128_si64(final); 
}


uint32_t hashGaloisFieldMultilinear(const uint64_t *  randomsource, const uint32_t *  string, const size_t length) {
    const uint32_t * const endstring = string + length;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32++));
    for(; string!= endstring; ++randomsource32,++string ) {
        __m128i temp = _mm_set_epi64x(*randomsource32,*string);
        __m128i clprod  = _mm_clmulepi64_si128( temp, temp, 0x10);
        acc = _mm_xor_si128 (clprod,acc);   
    }
    return barrettWithoutPrecomputation32(acc);
}



uint32_t hashGaloisFieldMultilinearHalfMultiplications(const uint64_t*  randomsource, const uint32_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint32_t * const endstring = string + length;
    const uint32_t *  randomsource32 = ( const uint32_t * )randomsource;
    __m128i acc = _mm_set_epi64x(0,*(randomsource32));
    randomsource32 += 1;
     for(; string!= endstring; randomsource32+=2,string+=2 ) {
        __m128i temp1 = _mm_set_epi64x(*randomsource32,*(randomsource32+1));
        __m128i temp2 = _mm_set_epi64x(*string,*(string+1));
        __m128i twosums = _mm_xor_si128(temp1,temp2); 
        __m128i clprod  = _mm_clmulepi64_si128( twosums, twosums, 0x10);
        acc = _mm_xor_si128 (clprod,acc);   
    }
    return barrettWithoutPrecomputation32(acc);
}


// a 64-bit version
uint64_t hashGaloisFieldfast64(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length/2*2;
     __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i clprod1  = _mm_clmulepi64_si128( temp1, temp2, 0x00);
        const __m128i clprod2  = _mm_clmulepi64_si128( temp1, temp2, 0x11);
        acc = _mm_xor_si128 (clprod1,acc);   
        acc = _mm_xor_si128 (clprod2,acc);   
    }
    return barrettWithoutPrecomputation64(acc);
}


// a 64-bit version with half the number of multiplications
uint64_t hashGaloisFieldfast64half(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 2 * 2 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length*2/2;
     __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string!= endstring; randomsource+=2,string+=2 ) {
        const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);   
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);   
    }
    return barrettWithoutPrecomputation64(acc);
}


// a 64-bit version with half the number of multiplications
uint64_t hashGaloisFieldfast64halfunrolled(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
    assert(length / 4 * 4 == length); // if not, we need special handling (omitted)
    const uint64_t * const endstring = string + length*4/4;
     __m128i acc = _mm_set_epi64x(0,*(randomsource));
    randomsource += 1;
    for(; string!= endstring; randomsource+=4,string+=4 ) {
        {const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);   
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);}   
        {const __m128i temp1 = _mm_lddqu_si128((__m128i * )randomsource);
        const __m128i temp2 = _mm_lddqu_si128((__m128i *) string);
        const __m128i add1 =  _mm_xor_si128 (temp1,temp2);   
        const __m128i clprod1  = _mm_clmulepi64_si128( add1, add1, 0x10);
        acc = _mm_xor_si128 (clprod1,acc);}
    }
    return barrettWithoutPrecomputation64(acc);
}


// like MHH
uint64_t referenceproduct(const uint64_t*  randomsource, const uint64_t *  string, const size_t length) {
	
	uint64_t low = 0;
	uint64_t high = 0;
	for(size_t i = 0; i<length; ++i) {
    __asm__ ("mulq %[v]\n"
             "addq %%rax,  %[rl]\n"
        "adcq %%rdx,  %[rh]\n"
             :  [rh] "+r" (high), [rl] "+r" (low)  : [u] "a" (randomsource[i]), [v] "r" (string[i])  :"rdx","cc");
	}
	return low+high;// should be modulo
}


#endif

#endif
