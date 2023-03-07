/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Simple serial implementation for finding contacts.
*/

#ifndef DIST_CUDA
#define DIST_CUDA

/*#include "dist_moduleV2/dist_common.h"
#include "dist_moduleV2/dist_out.h"
#include "dist_moduleV2/dist_coords.h"*/

#include "dist_common.h"
#include "dist_out.h"
#include "dist_coords.h"

namespace medyan::dist {	
			
	template <uint D, bool SELF>
	void find_distances_cuda(dOut<D,SELF> &out, Coords &c1, Coords &c2);
	
	
	
	template <uint D, bool SELF>
	inline void find_distances(dOut<D,SELF> &out, Coords &c, tag_simd<cuda, float> tag){
						
		std:: cout << "find_distances(dOut<D> &out, Coords &c1, Coords &c2, tag_simd<cuda,float> tag)" << std::endl;
		find_distances_cuda<D,true>(out, c, c);
	}

	
	template <uint D, bool SELF>
	inline void find_distances(dOut<D,SELF> &out, Coords &c1, Coords &c2, tag_simd<cuda, float> tag){
						
		std:: cout << "find_distances(dOut<D> &out, Coords &c1, Coords &c2, tag_simd<cuda,float> tag)" << std::endl;
		find_distances_cuda(out, c1, c2);
	}

} // end-of-namespace dist

#endif // DIST_CUDA
