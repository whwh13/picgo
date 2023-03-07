/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Many common typedefs and some miscellaneous functions for printing bits in variables.
*/

#ifndef DIST_COMMON
#define DIST_COMMON

#include <bitset>
#include <vector>

#include "umesimd/UMESimd.h"

typedef unsigned int uint;

namespace medyan::dist {

	// Algo examples below; T can curently be float for avx and avx2 and int16_t for avx2
	template <typename Algo, typename T>  
	struct tag_simd{};

	struct simd_no{}; // Serial
	struct simd_avx{}; // AVX or AVX2
	struct simd_avx_par{}; // AVX/2 multithreaded
	struct cuda{}; //Algo

	// UME::SIMD typedefs
	
	using uvec8_f = UME::SIMD::SIMD8_32f;
	using uvec8_i = UME::SIMD::SIMD8_32i;
	using uvec8_ui = UME::SIMD::SIMD8_32u;
	
	using umask8 = UME::SIMD::SIMDMask8;

	using std::cout;
	using std::endl;

	// Requires C++14 - for debugging
	#define show_type_name(_t) \
	    std::system(("echo " + std::string(typeid(_t).name()) + " | c++filt -t").c_str())

	// Mainly for debugging. Some custom functions to print the contents of vectors/matricies.
	template <typename X>
	void print_vec(const std::vector<X> &v){
		uint N = v.size();
		for(uint i=0; i<N; ++i)
			std::cout << v[i] << " ";
		std::cout << std::endl;
	}

	template <typename X>
	void print_mat(const std::vector<std::vector<X>> &m){
		uint N = m.shape()[0];
		for(uint i=0; i<N; ++i){
			for(uint j=0; j<N; ++j){
				std::cout << m[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

	// May be helpful in debugging when manipulating bits
	template <typename int_type>					 
	void print_bits(int_type c, int nlastbits=8){
	     std::bitset<8*sizeof(c)> x(c);
	     for(int i = nlastbits-1; i>=0 ; --i){
			std::cout << std::noboolalpha << x[i];
			if(i%8==0)
				cout << " ";
		}
		// cout << endl;
	}

} // end-of-namespace dist

#endif // DIST_COMMON
