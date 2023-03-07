/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Main public SIMD functions are: find_distances(dOut<D> &out, Coords &c) and 
find_distances(dOut<D> &out, Coords &c, tag_simd<Algo,T> tag).

The first should be considered as the default interface unless hand tuning is required 
by picking another SIMD algorithm.

Most other functions/variables below should be considered private to the module.

*/

#ifndef DIST_SIMD
#define DIST_SIMD

#include <thread>
#include <mutex>

//#include "dist_moduleV2/dist_simd_utils.h"

#include "dist_simd_utils.h"

namespace medyan::dist {

	extern std::once_flag avx_initialized, avx2_initialized; // thread-safe initialization of module variables

	void init_avx_module(uint N);

	void init_avx2_module(uint N);

	inline auto default_simd_algo(){
		return tag_simd<simd_avx,float>();
	}
	
	// Do cases not covered by SIMD serially.
	template <uint D, bool SELF, typename Algo, typename T>
	void dist_simd_corner_scalar(dOut<D,SELF> &out, Coords &c, tag_simd<Algo,T> tag){
	
		const uint simd_size = get_simd_size(tag);
	
		const uint N = c.size();
		const uint vec_size = N - N % simd_size;
		uint iend = vec_size - simd_size;
		uint jend = vec_size - simd_size;

		for(uint i=0; i<vec_size; ++i){
			uint irem = i%simd_size;
			for(uint j=jend+irem; j<N; ++j){
				dist_scalar_ij(out, c, c, i, j);
			}
		}
		
		for(uint i=vec_size; i<N; ++i){
			for(uint j=i+1; j<N; ++j){
				dist_scalar_ij(out, c, c, i, j);
			}
		}
	
	}
	
	// Do cases not covered by SIMD serially.
	template <uint D, bool SELF, typename Algo, typename T>
	void dist_simd_corner_scalar(dOut<D,SELF> &out, Coords &c1, Coords &c2, tag_simd<Algo,T> tag){

		const uint simd_size = get_simd_size(tag);

		const uint N1 = c1.size();
		const uint N2 = c2.size();
		const uint vec_size1 = N1 - N1 % simd_size;
		const uint vec_size2 = N2 - N2 % simd_size;
		uint iend = vec_size1 - simd_size;
		uint jend = vec_size2 - simd_size;

		for(uint i=0; i<vec_size1; ++i){
			uint irem = i%simd_size;
			for(uint j=jend+irem; j<N2; ++j){
				dist_scalar_ij(out, c1, c2, i, j);
			}
		}
		
		for(uint i=0; i<vec_size1; ++i){
			uint irem = i%simd_size;
			for(uint j=0; j<irem; ++j){
				dist_scalar_ij(out, c1, c2, i, j);
			}
		}
		
		for(uint i=vec_size1; i<N1; ++i){
			for(uint j=0; j<vec_size2; ++j){
				dist_scalar_ij(out, c1, c2, i, j);
			}
		}
		
		for(uint i=vec_size1; i<N1; ++i){
			for(uint j=vec_size2; j<N2; ++j){
				dist_scalar_ij(out, c1, c2, i, j);
			}
		}	
	}


	// This is the outer SIMD algorithm which iterates over i,j batches, and 
	// delegates the specific calculation based on the algo/type dual TAG.
	template <uint D, bool SELF, typename Algo, typename T>
	void find_distances(dOut<D,SELF> &out, Coords &c, tag_simd<Algo,T> tag)
	{	
		// show_type_name(tag); // uncomment on the left to check which algorithm got actually used

		const uint N = c.size();
	
		std::call_once(avx_initialized, init_avx_module, N);
		std::call_once(avx2_initialized, init_avx2_module, N);
		
		// init_special(out, c, tag);
	
		const uint simd_size = get_simd_size(tag);
	
		const uint vec_size = N - N % simd_size;
	
		uint iend = vec_size - simd_size;
		uint jend = vec_size - simd_size;

		// cout << "vec_size=" << vec_size << " simd_size=" << simd_size
		// 	 << " iend, jend: " << iend << endl;

		for(uint i = 0; i < vec_size; i += simd_size){
			uint jstart = i+1;
			// c.load_xyz_i(i,tag);
			uvec8_f c_vxi(&c.x[i]);
			uvec8_f c_vyi(&c.y[i]);
			uvec8_f c_vzi(&c.z[i]);
			uvec8_i c_finfo(&c.filinfo[i]);
			
			for(uint j = jstart; j < jend; ++j){
				dist_simd_ij(out, c.indices, c_vxi, c_vyi, c_vzi, c_finfo, c, i, j, tag);
			}
		}
	
		dist_simd_corner_scalar(out, c, tag);
	}

	// This is the outer SIMD algorithm which iterates over i,j batches, and
	// delegates the specific calculation based on the algo/type dual TAG.
	template <uint D, bool SELF, typename Algo, typename T>
	void find_distances(dOut<D,SELF> &out, Coords &c1, Coords &c2, tag_simd<Algo,T> tag)
	{
		// show_type_name(tag); // uncomment on the left to check which algorithm got actually used

		uint N1 = c1.size();
		uint N2 = c2.size();
		
		int N = std::max({N1,N2});

		std::call_once(avx_initialized, init_avx_module, N);
		std::call_once(avx2_initialized, init_avx2_module, N);

		// init_special(out, c1, c2, tag);

		const uint simd_size = get_simd_size(tag);

		const uint vec_size1 = N1 - N1 % simd_size;
		const uint vec_size2 = N2 - N2 % simd_size;

		uint iend = vec_size1 - simd_size;
		uint jend = vec_size2 - simd_size;

		// cout << "vec_size=" << vec_size << " simd_size=" << simd_size
		// 	 << " iend, jend: " << iend << endl;

		for(uint i = 0; i < vec_size1; i += simd_size){
			// c1.load_xyz_i(i,tag);
			
			uvec8_f c1_vxi(&c1.x[i]);
			uvec8_f c1_vyi(&c1.y[i]);
			uvec8_f c1_vzi(&c1.z[i]);
			uvec8_i c1_finfo(&c1.filinfo[i]);
			
			for(uint j = 0; j < jend; ++j){
				dist_simd_ij(out, c1.indices, c1_vxi, c1_vyi, c1_vzi, c1_finfo, c2, i, j,
						tag);
			}
		}

		dist_simd_corner_scalar(out, c1, c2, tag);
	}
	

} // end-of-namespace dist


#endif // DIST_SIMD



