/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

SIMD specific typedefs and many functions for bitwise or elementwise printing of SIMD variables.
*/

#ifndef DIST_SIMD_UTILS
#define DIST_SIMD_UTILS

#include "dist_common.h"
/*#include "dist_moduleV2/dist_coords.h"
#include "dist_moduleV2/dist_serial.h"*/
#include "dist_coords.h"
#include "dist_serial.h"

namespace medyan::dist {


	template <typename Algo> 
	constexpr inline uint get_simd_size(tag_simd<Algo,float> tag){
		return 8;
	}

	inline void print_bits_simd_epi16(__m256i c)
	{
		int16_t *x = (int16_t*)(&c);
	    for(int i = 15; i>=0 ; --i){
			std::bitset<16> y(x[i]);
			std::cout << y << " ";
		}
		cout << endl;
	}

	inline void print_bits_simd_epi32(__m256i c)
	{
		int *x = (int*)(&c);
		for(int i = 7; i>=0; --i) {
			std::bitset<32> y(x[i]);
			std::cout << y << " ";
		}
		cout << endl;
	}

	inline void print_m128i_as_int32(__m128i m){
		int32_t *x = (int32_t *)(&m);
		cout << x[3] << " " << x[2] << " " << x[1] << " " << x[0] << endl;
	}
	
	inline void print_m256i_as_int32(__m256i m){
		int32_t *x = (int32_t *)(&m);
		cout << x[7] << " " << x[6] << " " << x[5] << " " << x[4] << " "
	         << x[3] << " " << x[2] << " " << x[1] << " " << x[0] << endl;
	}

	inline void print_m256i_as_int16(__m256i m){
		int16_t *x = (int16_t *)(&m);
		cout << x[15] << " " << x[14] << " " << x[13] << " " << x[12] << " "
			 << x[11] << " " << x[10] << " " << x[9] << " " << x[8] << " "
		 	 << x[7] << " " << x[6] << " " << x[5] << " " << x[4] << " "
	         << x[3] << " " << x[2] << " " << x[1] << " " << x[0] << endl;
	}

	inline void print_m256_as_float(__m256 m){
		float *x = (float *)(&m);
		cout << x[7] << " " << x[6] << " " << x[5] << " " << x[4] << " "
			 << x[3] << " " << x[2] << " " << x[1] << " " << x[0] << endl;
	}

} // end-of-namespace dist

#endif // DIST_SIMD_UTILS

