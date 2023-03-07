/*
AUTHOR: G.A. Papoian, Date: Jan 5, 2019

A multithreaded versio of the code for both AVX and AVX2 calculations.
*/

#ifndef DIST_AVX_PAR
#define DIST_AVX_PAR

#ifdef __AVX__

#include <iostream>
#include <random>
#include <array>
#include <cmath>
#include <bitset>
#include <queue> 
#include <atomic>
#include <stdexcept>
#include <thread>


#include "dist_simd.h"
#include "dist_avx_aux.h"
#include "Util/Io/Log.hpp"

// # define PBLOCKDIMX 4
// # define PBLOCKDIMY 4
// # define PSTRIDE 4

# define PBLOCKDIMX 1
# define PBLOCKDIMY 64
# define PSTRIDE 128

// template <typename X>
// void print_array(const xt::xarray<X> &m){
// 	uint N = m.shape()[0];
// 	for(uint i=0; i<N; ++i){
// 		for(uint j=0; j<N; ++j){
// 			std::cout << m(i,j) << " ";
// 		}
// 		std::cout << std::endl;
// 	}
// }

namespace medyan::dist {

	// extern xt::xarray<int> mat;
	
	struct PGrid {
		uint x;
		uint y;
	};
	
	struct PBlock {
		uint x;
		uint y;
	};
	
	extern PGrid griddim;
	extern PBlock blockdim;
		
	extern std::array<std::atomic<uint>, 2> ncontacts_global_h;
	// extern std::array<uint, 4> ncontacts_global_h;
	
	extern std::queue<std::pair <uint, uint>> block_queue;
	
	template <uint D, bool SELF>
	void kernel_dist_simd_avx_par(dOut<D,SELF> &out, Coords &c1, Coords &c2);
	
	
	// Do cases not covered by SIMD serially.
	template <uint D, bool SELF, typename Algo, typename T>
	void dist_avx_par_corner_scalar(dOut<D,SELF> &out, Coords &c1, Coords &c2, tag_simd<Algo,T> tag){

		const uint simd_size = get_simd_size(tag);

		const uint N1 = c1.size();
		const uint N2 = c2.size();
		uint jstart = N2 - N2 % simd_size;

		// cout << "dist_avx_par_corner_scalar(...): N1=" << N1 << " N2=" << N2 << " jstart=" << jstart << endl;

		for(uint i=0; i<N1; ++i){
			for(uint j=jstart; j<N2; ++j){
				dist_scalar_ij(out, c1, c2, i, j);
			}
		}
	}
	
	
    template <uint D, bool SELF>
    inline void find_distances(dOut<D,SELF> &out, Coords &c1, Coords &c2, tag_simd<simd_avx_par, float> tag){

		LOG(ERROR) << "Parallel implementation causes invalid memory access and should not be used. Read https://github.com/medyan-dev/medyan/issues/67.";
		throw std::runtime_error("Parallel implementation is wrong and not yet fixed.");
		
		uint N1 = c1.size();
		uint N2 = c2.size();
		
//		cout << "find_distances(simd_avx_par): N1=" << N1 << " N2=" << N2 << endl;
		
		int N = std::max({N1,N2});
			
		std::call_once(avx_initialized, init_avx_module, N);
		std::call_once(avx2_initialized, init_avx2_module, N);
		
		const uint simd_size = get_simd_size(tag_simd<simd_avx_par,float>());
		
		blockdim.x = PBLOCKDIMX; 
		blockdim.y = PBLOCKDIMY;
		griddim.x = std::ceil(static_cast<float>((N1 + blockdim.x - 1)) / blockdim.x / PSTRIDE / simd_size);
		griddim.y = std::ceil(static_cast<float>((N2 + blockdim.y - 1)) / blockdim.y);
				
		constexpr uint simd_and_stride_length_x = PSTRIDE*simd_size;
		
		uint nblocks=0;
		for(uint by=0; by<griddim.y; ++by){
			uint py = by * PBLOCKDIMY; 
			for(uint bx=0; bx<griddim.x; ++bx){
			    uint px = (bx+1) * PBLOCKDIMX * simd_and_stride_length_x;
				if(!SELF || px >= py){
					// cout  << "by=" << by << " py=" << py << " bx=" << bx << " px=" << px << endl;
					block_queue.emplace(bx,by);
					++nblocks;
				}
			}			
		}
		
//		cout << "PBlockdim:" << blockdim.x << " " << blockdim.y << endl;
//		cout << "PGriddim:" << griddim.x << " " << griddim.y << endl;
//		cout << "nblocks:" << nblocks << endl;
		
		
		// mat = xt::zeros<int>({N1, N2});
		// xt::xarray<int> mat1 = xt::ones<int>({N1, N2});
		
		for(uint i=0; i<2; ++i)
			ncontacts_global_h[i]=0;

		unsigned int nthreads = std::thread::hardware_concurrency();

		nthreads = 1;
//		cout<<"#threads "<<nthreads<<endl;
//		cout << "find_distances(...tag_simd<simd_avx_par, float>...) nthreads=" << nthreads << endl;
		
		std::vector<std::thread> threads_avx;
		for(uint i=0; i<nthreads; ++i){
			threads_avx.push_back(std::thread(kernel_dist_simd_avx_par<D,SELF>, std::ref(out), std::ref(c1), std::ref(c2)));
		}

		for(auto &t : threads_avx)
			t.join();
	  
		// for(uint i=0; i<nblocks; ++i)
		// 	kernel_dist_simd_avx_par(out, c1, c2);
		
		// std::thread t1(kernel_dist_simd_avx_par, out, c1, c2, tag)
		
		std::copy(ncontacts_global_h.begin(), ncontacts_global_h.begin()+D, out.counter.begin());
		
		dist_avx_par_corner_scalar(out, c1, c2, tag);
		
		// print_array(mat);
		
		// for(uint d=0; d<D; ++d)
		// 	cout << "kernel_dist_simd_avx_par d=" << d << " ncontacts=" << out.counter[d] << endl;
		
		// cout << "Contacts found: " << out.counter[0] << ", capacity=" << out.capacity() << endl;
		// for(uint i=0; i<out.counter[0]; ++i)
		// 	cout << "find_distances_simd_avx_par: i=" << out.dout[0][i] << " j=" << out.dout[1][i] << endl;
		// cout << endl;
		
		// auto comp = mat == mat1;
		
		// cout << comp;
	}
	
    template <uint D>
    inline void find_distances(dOut<D,true> &out, Coords &c, tag_simd<simd_avx_par, float> tag){
		find_distances(out, c, c, tag);
	}
	
} // end-of-namespace dist

#endif // __AVX__

#endif // DIST_AVX_PAR
