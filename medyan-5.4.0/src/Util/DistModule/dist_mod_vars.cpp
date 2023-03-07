/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Initializing various module variables.
*/

//#include "dist_moduleV2/dist_out.h"
//#include "dist_moduleV2/dist_simd.h"
//#include "dist_moduleV2/dist_simd_utils.h"

#ifdef SIMDBINDINGSEARCH
#include "dist_out.h"
#include "dist_simd.h"
#include "dist_simd_utils.h"

namespace medyan::dist {
	
	std::once_flag avx_initialized, avx2_initialized;

	alignas(32) uint _lut_avx[256*8]; 
	alignas(32) uint _lut_avx_step2[5*8]; 
	alignas(32) uint _lut_avx_blend[5*8];
	alignas(32) uint _lut_avx2[256*8]; 
	alignas(32) uint _lut_128[16*4]; 

	__m256i zero_m256i; // all bits initialized to 0

	void init_lut_avx2()
	{
		// cout << "init__lut()..." << endl;

		for(uint mask=0; mask<256; ++mask) {
			uint cnt = 0;
			for(uint i=0;i<8;++i){
				if(mask & (1<<i)){
					_lut_avx2[mask*8+cnt++]=i;
				}
			}
			// cout << "mask: "; print_bits(mask);
			// cout << view(_lut_avx2,range(mask*8,(mask+1)*8)) << endl;;
		}
	}
	
	void init_lut_128()
	{
		// cout << "init_lut_128()..." << endl;
		// lut = lut -1;
		for(uint mask=0; mask<16; ++mask) {
			uint cnt = 0;
			for(uint i=0;i<4;++i){
				if(mask & (1<<i)){
					_lut_128[mask*4+cnt++]=i;
				}
			}
			// cout << "mask: "; print_bits(mask,4);
			// cout << xt::view(_lut_128,xt::range(mask*4,(mask+1)*4)) << endl << endl;
		}
	}

	void init_lut_avx()
	{
		// cout << "init__lut()..." << endl;

		int maskr = 0b000'1111;
		int maskl = 0b1111'0000;
		
		for(uint mask=0; mask<256; ++mask) {

			int condr = mask & maskr;
			int condl = (mask & maskl) >> 4;

			std::copy(&_lut_128[4*condr],&_lut_128[4*(condr+1)],&_lut_avx[8*mask]);
			std::copy(&_lut_128[4*condl],&_lut_128[4*(condl+1)],&_lut_avx[8*mask+4]);
		}
	}
	
	void init_lut_avx_step2()
	{		
		std::array<int,5> mask_compact;

		mask_compact[0] = 0b00'00'00'00'11'10'01'00;		
		mask_compact[1] = 0b00'00'00'11'10'01'00'00;		
		mask_compact[2] = 0b00'00'11'10'01'00'00'00;		
		mask_compact[3] = 0b00'11'10'01'00'00'00'00;		
		mask_compact[4] = 0b11'10'01'00'00'00'00'00;
		
		for(uint i=0; i<5; ++i){
			for(uint m=0; m<8; ++m){
				_lut_avx_step2[8*i+m] = (mask_compact[i] & (0b11 << 2*m)) >> 2*m; 
			}				
		}
	}
	
	void init_lut_avx_blendv()
	{	
		for(uint i=0; i<5; ++i)
			for(uint m=0; m<8; ++m)
				_lut_avx_blend[8*i+m] = (m>=i)? 0 : ~0;
	}
	
	void init_avx_module(uint N)
	{	
		init_lut_128();
		init_lut_avx();
		init_lut_avx_step2();
		init_lut_avx_blendv();
		
		zero_m256i = _mm256_setzero_si256();
		// one_m256i = ~zero_m256i;
	}
	
	void init_avx2_module(uint N)
	{	
		init_lut_avx2();
		zero_m256i = _mm256_setzero_si256();	
	}

} // end-of-namespace dist

#endif