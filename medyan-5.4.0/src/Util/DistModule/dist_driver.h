/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

This is the main header that should be included in the client code.
*/

#ifndef DIST_DRIVER
#define DIST_DRIVER

#include <random>

/*#include "dist_moduleV2/dist_common.h"
#include "dist_moduleV2/dist_out.h"
#include "dist_moduleV2/dist_coords.h"
#include "dist_moduleV2/dist_serial.h"
#include "dist_moduleV2/dist_simd.h"
#include "dist_moduleV2/dist_avx.h"
#include "dist_moduleV2/dist_avx_par.h"
#include "dist_moduleV2/dist_cuda.h"*/

#include "dist_common.h"
#include "dist_out.h"
#include "dist_coords.h"
#include "dist_serial.h"
#include "dist_simd.h"
#include "dist_avx.h"
#include "dist_avx_par.h"
#include "dist_cuda.h"

namespace medyan::dist {

	// This is a simpler, default interface to the library, which tries to pick 
	// the best algorithm given the available SIMD instruction set.
	template <uint D, bool SELF>
	void find_distances(dOut<D,SELF> &out, Coords &c)
	{
		find_distances(out, c, default_simd_algo());
	}
	
	// This is a simpler, default interface to the library, which tries to pick
	// the best algorithm given the available SIMD instruction set.
	template <uint D, bool SELF>
	void find_distances(dOut<D,SELF> &out, Coords &c1, Coords &c2)
	{
		find_distances(out, c1, c2, default_simd_algo());
	}

	
	// // these are needed to launch std::thread
	// inline void find_distances_1(dOut<1> &out, Coords &c){
	// 	find_distances<1>(out,c);
	// }
	//
	// // these are needed to launch std::thread
	// inline void find_distances_2(dOut<2> &out, Coords &c){
	// 	find_distances<2>(out,c);
	// }
	//
	// // these are needed to launch std::thread
	// inline void find_distances_cmp2_1(dOut<1> &out, Coords &c1, Coords &c2){
	// 	find_distances<1>(out,c1, c2);
	// }
	//
	// // these are needed to launch std::thread
	// inline void find_distances_cmp2_2(dOut<2> &out, Coords &c1, Coords &c2){
	// 	find_distances<2>(out,c1, c2);
	// }


	// General functions that can be called from main to test various functionalities
	void benchmark_dist_module(uint N1 = 6000, uint N2 = 5600);
	void examples_dist_module();
	void dist_example_avx_veccoord();

	// After find_distances(...) returns, this function can be used enquire aboutthe calculated contacts
	template <uint D>
	void report_contact_stats(const dOut<D> &out)
	{
		std::random_device rd;
		std::mt19937 mt(rd());
	
		for(uint d=0; d<D; ++d){
			std::uniform_int_distribution<uint> dist_i(0, 100);
			uint ii = dist_i(mt) & ~1; //make ii even
	
		    cout << endl << "Num contacts found=" << out.counter[d] << endl;
		    cout         << " ii=" << ii << ": " << out.dout[2*d][ii] << " " << out.dout[2*d+1][ii] << endl;
		}
	}

} // end-of-namespace dist

#endif // DIST_DRIVER
