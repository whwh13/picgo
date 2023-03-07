/*
AUTHOR: G.A. Papoian, Date: Nov 28, 2018

A workaround around the absence of _mm256_permutevar8x32_epi32 in AVX.
*/

#ifndef DIST_AVX_AUX
#define DIST_AVX_AUX

#ifdef __AVX__

#include <iostream>
#include <random>
#include <array>
#include <cmath>
#include <bitset>

//#include "dist_moduleV2/dist_simd.h"
#include "dist_simd.h"


namespace medyan::dist {

    extern uint _lut_avx      [256*8];
    extern uint _lut_avx_step2[5*8]; 
    extern uint _lut_avx_blend[5*8];
		
	// In AVX, _mm256_permutevar8x32_epi32 does not exist (it was introduced in AVX2)
	// We have to implement it using more elementary instructions.
	// _mm256_permutevar8x32_epi32 is more effecient, when available. 
	// Maybe the code below could be further optimized.
	inline void swizzlea_8x32_avx_impl(uvec8_i &ind, int icond){
		__m256i vmask = _mm256_load_si256((__m256i*)(&_lut_avx[8*icond]));
		
		__m256i ind_p = _mm256_castps_si256(_mm256_permutevar_ps(*(__m256*)(&(ind.mVec)), vmask));
		__m128i ind_p_l128 = _mm256_extractf128_si256(ind_p,1);		
		__m256i ind_p_l256 = _mm256_castps_si256(_mm256_broadcast_ps((__m128*)(&ind_p_l128)));
		
		uint c_r = _mm_popcnt_u32(icond & 0b0000'1111);
		__m256i vmask_step2 = _mm256_load_si256((__m256i*)(&_lut_avx_step2[8*c_r]));
		
		__m256i ind_p_l256_p = _mm256_castps_si256(_mm256_permutevar_ps(*(__m256*)(&ind_p_l256), vmask_step2));
		
		uint c_l = _mm_popcnt_u32(icond & 0b1111'0000);
		
		__m256i vmask_blend = _mm256_load_si256((__m256i*)(&_lut_avx_blend[8*c_r]));
							
		ind.mVec = _mm256_castps_si256(_mm256_blendv_ps(*(__m256*)(&(ind_p_l256_p)), *(__m256*)(&(ind_p)), *(__m256*)(&(vmask_blend)) ) );
	}
	
} // end-of-namespace medyan::dist

#endif // __AVX__

#endif // DIST_AVX_AUX
