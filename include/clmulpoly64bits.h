#ifndef CLMULPOLY64BITS_H_
#define CLMULPOLY64BITS_H_

#include "clmul.h"

// simple 64-bit polynomial hashing, uses only one key
// not expected to be fast!
uint64_t hashGaloisFieldPoly64(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i key = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	for (; string < endstring; ++string) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, key, 0x00);
		acc = barrettWithoutPrecomputation64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

uint64_t precomphashGaloisFieldPoly64(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i key = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	for (; string < endstring; ++string) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, key, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

uint64_t fasthashGaloisFieldPoly64_2_noprecomp(const void* rs,
		const uint64_t * string, const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = barrettWithoutPrecomputation64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = barrettWithoutPrecomputation64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = barrettWithoutPrecomputation64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

uint64_t fasthashGaloisFieldPoly64_2(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = precompReduction64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = precompReduction64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

// fast 64-bit polynomial hashing, uses only one key
// expected to be fast!
//TODO: can use more keys for increased universality
uint64_t fasthashGaloisFieldPoly64_4(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = precompReduction64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		// powers of the keys are packed into two registers
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		if (string + 3 < endstring) {
			__m128i tkey3 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
			__m128i tkey4 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
			__m128i key2 = _mm_xor_si128(
					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
					_mm_slli_si128(tkey4, 8));
			for (; string + 3 < endstring; string += 4) {
				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
				acc = precompReduction64_si128(
						_mm_xor_si128(acc,
								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
										_mm_xor_si128(clprod2, clprod3))));
			}
		}
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = precompReduction64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

// fast 64-bit polynomial hashing, uses only one key
// expected to be fast!
uint64_t fasthashGaloisFieldPoly64_8(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = precompReduction64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		// powers of the keys are packed into two registers
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		if (string + 3 < endstring) {
			__m128i tkey3 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
			__m128i tkey4 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
			__m128i key2 = _mm_xor_si128(
					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
					_mm_slli_si128(tkey4, 8));
			if (string + 7 < endstring) {
				__m128i tkey5 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
				__m128i tkey6 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
				__m128i tkey7 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
				__m128i tkey8 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
				__m128i key3 = _mm_xor_si128(
						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey6, 8));
				__m128i key4 = _mm_xor_si128(
						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey8, 8));
				for (; string + 7 < endstring; string += 8) {
					__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4)); //a5 a6
					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6)); //a7 a8
					const __m128i x1 = _mm_srli_si128(temp4, 8); //a8
					const __m128i clprod1 = _mm_clmulepi64_si128(temp4, key,
							0x00); //a7*k
					const __m128i clprod2 = _mm_clmulepi64_si128(temp3, key,
							0x11); //a6*k^2
					const __m128i clprod3 = _mm_clmulepi64_si128(temp3, key2,
							0x00); //a5*k^3
					const __m128i clprod4 = _mm_clmulepi64_si128(temp2, key2,
							0x11); //a4*k^4
					const __m128i clprod5 = _mm_clmulepi64_si128(temp2, key3,
							0x00); //a3*k^5
					const __m128i clprod6 = _mm_clmulepi64_si128(temp, key3,
							0x11); //a2*k^6
					const __m128i clprod7 = _mm_clmulepi64_si128(temp, key4,
							0x00); //a1*k^7
					acc = _mm_clmulepi64_si128(acc, key4, 0x10); //k^8

					const __m128i t1 = _mm_xor_si128(x1, clprod1);
					const __m128i t2 = _mm_xor_si128(clprod2, clprod3);
					const __m128i t3 = _mm_xor_si128(clprod4, clprod5);
					const __m128i t4 = _mm_xor_si128(clprod6, clprod7);

					const __m128i z1 = _mm_xor_si128(t1, t2);
					const __m128i z2 = _mm_xor_si128(t3, t4);

					acc = precompReduction64_si128(
							_mm_xor_si128(acc, _mm_xor_si128(z1, z2)));
				}

			}
			for (; string + 3 < endstring; string += 4) {
				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
				acc = precompReduction64_si128(
						_mm_xor_si128(acc,
								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
										_mm_xor_si128(clprod2, clprod3))));
			}
		}
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = precompReduction64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

// experimental
uint64_t fasthashGaloisFieldPoly64_16(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = precompReduction64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		// powers of the keys are packed into two registers
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		if (string + 3 < endstring) {
			__m128i tkey3 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
			__m128i tkey4 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
			__m128i key2 = _mm_xor_si128(
					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
					_mm_slli_si128(tkey4, 8));
			if (string + 15 < endstring) {
				__m128i tkey5 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
				__m128i tkey6 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
				__m128i tkey7 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
				__m128i tkey8 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
				__m128i tkey9 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey5, tkey4, 0x00));
				__m128i tkey10 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey5, tkey5, 0x00));
				__m128i tkey11 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey6, tkey5, 0x00));
				__m128i tkey12 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey6, tkey6, 0x00));
				__m128i tkey13 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey6, tkey7, 0x00));
				__m128i tkey14 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey7, tkey7, 0x00));
				__m128i tkey15 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey8, tkey7, 0x00));
				__m128i tkey16 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey8, tkey8, 0x00));
				__m128i key3 = _mm_xor_si128(
						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey6, 8));
				__m128i key4 = _mm_xor_si128(
						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey8, 8));
				__m128i key5 = _mm_xor_si128(
						_mm_and_si128(tkey9, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey10, 8));
				__m128i key6 = _mm_xor_si128(
						_mm_and_si128(tkey11, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey12, 8));
				__m128i key7 = _mm_xor_si128(
						_mm_and_si128(tkey13, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey14, 8));
				__m128i key8 = _mm_xor_si128(
						_mm_and_si128(tkey15, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey16, 8));

				for (; string + 15 < endstring; string += 16) {
					__m128i temp = _mm_lddqu_si128((__m128i *) string);
					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2));
					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4));
					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6));
					__m128i temp5 = _mm_lddqu_si128((__m128i *) (string + 8));
					__m128i temp6 = _mm_lddqu_si128((__m128i *) (string + 10));
					__m128i temp7 = _mm_lddqu_si128((__m128i *) (string + 12));
					__m128i temp8 = _mm_lddqu_si128((__m128i *) (string + 14));
					 __m128i x1 = _mm_srli_si128(temp8, 8);
					 __m128i clprod1 = _mm_clmulepi64_si128(temp8, key,
							0x00);
					 __m128i clprod2 = _mm_clmulepi64_si128(temp7, key,
							0x11);
					 __m128i clprod3 = _mm_clmulepi64_si128(temp7, key2,
							0x00);
					 __m128i clprod4 = _mm_clmulepi64_si128(temp6, key2,
							0x11);
					 __m128i clprod5 = _mm_clmulepi64_si128(temp6, key3,
							0x00);
					 __m128i clprod6 = _mm_clmulepi64_si128(temp5, key3,
							0x11);
					 __m128i clprod7 = _mm_clmulepi64_si128(temp5, key4,
							0x00);
					 __m128i clprod8 = _mm_clmulepi64_si128(temp4, key4,
							0x11);
					 __m128i clprod9 = _mm_clmulepi64_si128(temp4, key5,
							0x00);
					 __m128i clprod10 = _mm_clmulepi64_si128(temp3, key5,
							0x11);
					 __m128i clprod11 = _mm_clmulepi64_si128(temp3, key6,
							0x00);
					 __m128i clprod12 = _mm_clmulepi64_si128(temp2, key6,
							0x11);
					 __m128i clprod13 = _mm_clmulepi64_si128(temp2, key7,
							0x00);
					 __m128i clprod14 = _mm_clmulepi64_si128(temp, key7,
							0x11);
					 __m128i clprod15 = _mm_clmulepi64_si128(temp, key8,
							0x00);
					acc = _mm_clmulepi64_si128(acc, key8, 0x10);
					 __m128i t1 = _mm_xor_si128(x1, clprod1);
					 __m128i t2 = _mm_xor_si128(clprod2, clprod3);
					 __m128i t3 = _mm_xor_si128(clprod4, clprod5);
					 __m128i t4 = _mm_xor_si128(clprod6, clprod7);
					 __m128i t5 = _mm_xor_si128(clprod8, clprod9);
					 __m128i t6 = _mm_xor_si128(clprod10, clprod11);
					 __m128i t7 = _mm_xor_si128(clprod12, clprod13);
					 __m128i t8 = _mm_xor_si128(clprod14, clprod15);

					 __m128i z1 = _mm_xor_si128(t1, t2);
					 __m128i z2 = _mm_xor_si128(t3, t4);
					 __m128i z3 = _mm_xor_si128(t5, t6);
					 __m128i z4 = _mm_xor_si128(t7, t8);
					 __m128i Z1 = _mm_xor_si128(z1, z2);
					 __m128i Z2 = _mm_xor_si128(z3, z4);
					acc = precompReduction64_si128(
							_mm_xor_si128(acc, _mm_xor_si128(Z1, Z2)));

				}

			}
			for (; string + 3 < endstring; string += 4) {
				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
				acc = precompReduction64_si128(
						_mm_xor_si128(acc,
								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
										_mm_xor_si128(clprod2, clprod3))));
			}
		}
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = precompReduction64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

// experimental, excepted to be faster than fasthashGaloisFieldPoly64_8
// but is not
uint64_t halfhashGaloisFieldPoly64_8(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = precompReduction64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		// powers of the keys are packed into two registers
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		if (string + 3 < endstring) {
			__m128i tkey3 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
			__m128i tkey4 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
			__m128i key2 = _mm_xor_si128(
					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
					_mm_slli_si128(tkey4, 8));
			if (string + 7 < endstring) {
				__m128i tkey5 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
				__m128i tkey6 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
				__m128i tkey7 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
				__m128i tkey8 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
				__m128i key3 = _mm_xor_si128(
						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey6, 8));
				__m128i key4 = _mm_xor_si128(
						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey8, 8));

				/**
				 * For the half multiplication thing, we want to have the keys organized
				 * as follows:
				 * (0,k),(k^2,k^3)...
				 */
				// TODO: these keys are probably wrong
				__m128i hkey1 = _mm_slli_si128(key,8);
				__m128i hkey2 = _mm_xor_si128(_mm_srli_si128(key,8), _mm_slli_si128(key2,8));
				__m128i hkey3 = _mm_xor_si128(_mm_srli_si128(key2,8), _mm_slli_si128(key3,8));
				__m128i hkey4 = _mm_xor_si128(_mm_srli_si128(key3,8), _mm_slli_si128(key4,8));
				for (; string + 7 < endstring; string += 8) {
					__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4)); //a5 a6
					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6)); //a7 a8
					__m128i t1 = _mm_xor_si128(hkey1,temp);
					__m128i t2 = _mm_xor_si128(hkey2,temp2);
					__m128i t3 = _mm_xor_si128(hkey3,temp3);
					__m128i t4 = _mm_xor_si128(hkey4,temp4);
					__m128i clprod1 = _mm_clmulepi64_si128(t1,t1,0x10);
					__m128i clprod2 = _mm_clmulepi64_si128(t2,t2,0x10);
					__m128i clprod3 = _mm_clmulepi64_si128(t3,t3,0x10);
					__m128i clprod4 = _mm_clmulepi64_si128(t4,t4,0x10);
					__m128i tacc = _mm_clmulepi64_si128(acc, key4, 0x10); //k^8
					const __m128i b1 = _mm_xor_si128(clprod2, clprod1);
					const __m128i b2 = _mm_xor_si128(clprod3, clprod4);
					const __m128i z1 = _mm_xor_si128(b1, b2);
					acc = precompReduction64_si128
							(
							_mm_xor_si128(tacc, z1));
				}

			}
			for (; string + 3 < endstring; string += 4) {
				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
				acc = precompReduction64_si128
						(
						_mm_xor_si128(acc,
								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
										_mm_xor_si128(clprod2, clprod3))));
			}
		}
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = precompReduction64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}


uint64_t halfhashGaloisFieldPoly64_16(const void* rs, const uint64_t * string,
		const size_t length) {
	const uint64_t * randomsource = (const uint64_t *) rs;
	assert(*randomsource != 0); //otherwise silly
	const uint64_t * const endstring = string + length;
	__m128i tkey1 = _mm_set_epi64x(0, *(randomsource));
	__m128i acc = _mm_set_epi64x(0, *string);
	++string;
	if (string + 1 < endstring) {
		// we start by precomputing the powers of the key
		__m128i tkey2 = precompReduction64_si128(
				_mm_clmulepi64_si128(tkey1, tkey1, 0x00));
		// powers of the keys are packed into two registers
		__m128i key = _mm_xor_si128(_mm_and_si128(tkey1, _mm_set_epi64x(0, -1)),
				_mm_slli_si128(tkey2, 8));
		if (string + 3 < endstring) {
			__m128i tkey3 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey1, 0x00));
			__m128i tkey4 = precompReduction64_si128(
					_mm_clmulepi64_si128(tkey2, tkey2, 0x00));
			__m128i key2 = _mm_xor_si128(
					_mm_and_si128(tkey3, _mm_set_epi64x(0, -1)),
					_mm_slli_si128(tkey4, 8));
			if (string + 15 < endstring) {
				__m128i tkey5 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey2, tkey3, 0x00));
				__m128i tkey6 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey3, 0x00));
				__m128i tkey7 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey3, tkey4, 0x00));
				__m128i tkey8 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey4, tkey4, 0x00));
				__m128i tkey9 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey5, tkey4, 0x00));
				__m128i tkey10 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey5, tkey5, 0x00));
				__m128i tkey11 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey6, tkey5, 0x00));
				__m128i tkey12 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey6, tkey6, 0x00));
				__m128i tkey13 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey6, tkey7, 0x00));
				__m128i tkey14 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey7, tkey7, 0x00));
				__m128i tkey15 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey8, tkey7, 0x00));
				__m128i tkey16 = precompReduction64_si128(
						_mm_clmulepi64_si128(tkey8, tkey8, 0x00));
				__m128i key3 = _mm_xor_si128(
						_mm_and_si128(tkey5, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey6, 8));
				__m128i key4 = _mm_xor_si128(
						_mm_and_si128(tkey7, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey8, 8));
				__m128i key5 = _mm_xor_si128(
						_mm_and_si128(tkey9, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey10, 8));
				__m128i key6 = _mm_xor_si128(
						_mm_and_si128(tkey11, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey12, 8));
				__m128i key7 = _mm_xor_si128(
						_mm_and_si128(tkey13, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey14, 8));
				__m128i key8 = _mm_xor_si128(
						_mm_and_si128(tkey15, _mm_set_epi64x(0, -1)),
						_mm_slli_si128(tkey16, 8));

				/**
				 * For the half multiplication thing, we want to have the keys organized
				 * as follows:
				 * (0,k),(k^2,k^3)...
				 */
				// TODO: these keys are probably wrong
				__m128i hkey1 = _mm_slli_si128(key,8);
				__m128i hkey2 = _mm_xor_si128(_mm_srli_si128(key,8), _mm_slli_si128(key2,8));
				__m128i hkey3 = _mm_xor_si128(_mm_srli_si128(key2,8), _mm_slli_si128(key3,8));
				__m128i hkey4 = _mm_xor_si128(_mm_srli_si128(key3,8), _mm_slli_si128(key4,8));
				__m128i hkey5 = _mm_xor_si128(_mm_srli_si128(key4,8), _mm_slli_si128(key5,8));
				__m128i hkey6 = _mm_xor_si128(_mm_srli_si128(key5,8), _mm_slli_si128(key6,8));
				__m128i hkey7 = _mm_xor_si128(_mm_srli_si128(key6,8), _mm_slli_si128(key7,8));
				__m128i hkey8 = _mm_xor_si128(_mm_srli_si128(key7,8), _mm_slli_si128(key8,8));

				for (; string + 15 < endstring; string += 16) {
					__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
					__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
					__m128i temp3 = _mm_lddqu_si128((__m128i *) (string + 4)); //a5 a6
					__m128i temp4 = _mm_lddqu_si128((__m128i *) (string + 6)); //a7 a8
					__m128i temp5 = _mm_lddqu_si128((__m128i *) (string + 8)); //a7 a8
					__m128i temp6 = _mm_lddqu_si128((__m128i *) (string + 10)); //a7 a8
					__m128i temp7 = _mm_lddqu_si128((__m128i *) (string + 12)); //a7 a8
					__m128i temp8 = _mm_lddqu_si128((__m128i *) (string + 14)); //a7 a8
					__m128i t1 = _mm_xor_si128(hkey1,temp);
					__m128i t2 = _mm_xor_si128(hkey2,temp2);
					__m128i t3 = _mm_xor_si128(hkey3,temp3);
					__m128i t4 = _mm_xor_si128(hkey4,temp4);
					__m128i t5 = _mm_xor_si128(hkey5,temp);
					__m128i t6 = _mm_xor_si128(hkey6,temp2);
					__m128i t7 = _mm_xor_si128(hkey7,temp3);
					__m128i t8 = _mm_xor_si128(hkey8,temp4);
					__m128i clprod1 = _mm_clmulepi64_si128(t1,t1,0x10);
					__m128i clprod2 = _mm_clmulepi64_si128(t2,t2,0x10);
					__m128i clprod3 = _mm_clmulepi64_si128(t3,t3,0x10);
					__m128i clprod4 = _mm_clmulepi64_si128(t4,t4,0x10);
					__m128i clprod5 = _mm_clmulepi64_si128(t5,t5,0x10);
					__m128i clprod6 = _mm_clmulepi64_si128(t6,t6,0x10);
					__m128i clprod7 = _mm_clmulepi64_si128(t7,t7,0x10);
					__m128i clprod8 = _mm_clmulepi64_si128(t8,t8,0x10);
					__m128i tacc = _mm_clmulepi64_si128(acc, key8, 0x10);
					const __m128i b1 = _mm_xor_si128(clprod2, clprod1);
					const __m128i b2 = _mm_xor_si128(clprod3, clprod4);
					const __m128i b3 = _mm_xor_si128(clprod5, clprod6);
					const __m128i b4 = _mm_xor_si128(clprod7, clprod8);
					const __m128i z1 = _mm_xor_si128(b1, b2);
					const __m128i z2 = _mm_xor_si128(b3, b4);
					const __m128i Z1 = _mm_xor_si128(z1, z2);
					acc = precompReduction64_si128
							(
							_mm_xor_si128(tacc, Z1));
				}

			}
			for (; string + 3 < endstring; string += 4) {
				__m128i temp = _mm_lddqu_si128((__m128i *) string); //a1 a2
				__m128i temp2 = _mm_lddqu_si128((__m128i *) (string + 2)); //a3 a4
				const __m128i x1 = _mm_srli_si128(temp2, 8); //a4
				const __m128i clprod1 = _mm_clmulepi64_si128(temp, key2, 0x00); //a1*k^3
				const __m128i clprod2 = _mm_clmulepi64_si128(temp, key, 0x11); //a2*k^2
				const __m128i clprod3 = _mm_clmulepi64_si128(temp2, key, 0x00); //a3*k
				acc = _mm_clmulepi64_si128(acc, key2, 0x10); //k^4
				acc = precompReduction64_si128
						(
						_mm_xor_si128(acc,
								_mm_xor_si128(_mm_xor_si128(x1, clprod1),
										_mm_xor_si128(clprod2, clprod3))));
			}
		}
		for (; string + 1 < endstring; string += 2) {
			__m128i temp = _mm_lddqu_si128((__m128i *) string);
			acc = _mm_clmulepi64_si128(acc, key, 0x10);
			const __m128i clprod1 = _mm_clmulepi64_si128(temp, key, 0x00);
			acc = _mm_xor_si128(clprod1, acc);
			acc = _mm_xor_si128(acc, _mm_srli_si128(temp, 8));
			acc = precompReduction64_si128(acc);
		}
	}
	if (string < endstring) {
		const __m128i temp = _mm_set_epi64x(0, *string);
		const __m128i multi = _mm_clmulepi64_si128(acc, tkey1, 0x00);
		acc = precompReduction64_si128(multi);
		acc = _mm_xor_si128(acc, temp);
	}
	return _mm_cvtsi128_si64(acc);
}

#endif /* CLMULPOLY64BITS_H_ */