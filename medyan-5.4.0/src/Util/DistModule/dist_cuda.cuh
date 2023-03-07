#ifndef DIST_CUDA_IMPL
#define DIST_CUDA_IMPL

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>

/*#include "dist_moduleV2/dist_out.h"
#include "dist_moduleV2/dist_coords.h"

#include "dist_moduleV2/dist_cuda_simple.cuh"*/

#include "dist_out.h"
#include "dist_coords.h"

#include "dist_cuda_simple.cuh"

#define gpuAssert( condition ) { if( (condition) != 0 ) { fprintf( stderr, "\n FAILURE %s in %s, line %d\n", cudaGetErrorString(condition), __FILE__, __LINE__ ); exit( 1 ); } }

namespace medyan::dist {

	typedef unsigned int uint;

	using namespace std;
	
	uint ncontacts_h[4]; // assuming max of 4 comparisons	
	
	template <uint D, bool SELF>
	void find_distances_cuda(dOut<D,SELF> &out, Coords &c1, Coords &c2)
	{
		
		cout << "find_distances_cuda(...): Self=" << boolalpha << SELF << endl;
		
		uint N1 = c1.size();
		uint N2 = c2.size();
		uint Ni = out.capacity();

		uint zero[4] = {0,0,0,0};

		auto tic = std::chrono::steady_clock::now();
		thrust::device_vector<float> x1_d(c1.x);
		float *x1 = thrust::raw_pointer_cast(x1_d.data());
		thrust::device_vector<float> y1_d(c1.y);
		float *y1 = thrust::raw_pointer_cast(y1_d.data());
		thrust::device_vector<float> z1_d(c1.z);
		float *z1 = thrust::raw_pointer_cast(z1_d.data());
		thrust::device_vector<int> indices1_d(c1.indices);
		int *indices1 = thrust::raw_pointer_cast(indices1_d.data());


		thrust::device_vector<float> x2_d(c2.x);
		float *x2 = thrust::raw_pointer_cast(x2_d.data());
		thrust::device_vector<float> y2_d(c2.y);
		float *y2 = thrust::raw_pointer_cast(y2_d.data());
		thrust::device_vector<float> z2_d(c2.z);
		float *z2 = thrust::raw_pointer_cast(z2_d.data());
		thrust::device_vector<int> indices2_d(c2.indices);
		int *indices2 = thrust::raw_pointer_cast(indices2_d.data());		

		dim3 block(BLOCKX, BLOCKY);
		dim3 grid((N1 + block.x - 1) / block.x, (N2 + block.y - 1) / block.y);

		cout << "Block:" << block.x << " " << block.y << endl;
		cout << "Grid:" << grid.x << " " << grid.y << endl;
		
		gpuAssert( cudaMemcpyToSymbol(ncontacts_global, zero, D*sizeof(uint)) );
				
		thrust::device_vector<int*> indices_out_d(2*D);
		thrust::device_vector<int> *i_or_j_vecs[2*D];
		for(uint d=0; d<2*D; ++d){
			i_or_j_vecs[d] = new thrust::device_vector<int>(Ni);
			indices_out_d[d] = thrust::raw_pointer_cast(i_or_j_vecs[d]->data());
		}

		thrust::device_vector<float> dlim(4);
		thrust::copy(&out.dt[0],&out.dt[0]+2*D, dlim.begin());

   		auto tac = std::chrono::steady_clock::now();
   		auto us = std::chrono::duration_cast<std::chrono::microseconds>(tac-tic).count();
   		std::cout << "Time(CUDA: Before the kernel launch): " << us << " us" << std::endl << std::endl;

	    tic = std::chrono::steady_clock::now();

		find_distances_kernel_simple<D,SELF><<<grid, block>>>(N1, N2, Ni, thrust::raw_pointer_cast(dlim.data()),
		                                              x1, y1, z1, x2, y2, z2,
		                                              indices1, indices2,
							                          thrust::raw_pointer_cast(indices_out_d.data()));
		gpuAssert( cudaDeviceSynchronize() );
		
   		tac = std::chrono::steady_clock::now();
   		us = std::chrono::duration_cast<std::chrono::microseconds>(tac-tic).count();
   		std::cout << "Time(CUDA: Kernel time): " << us << " us" << std::endl << std::endl;
		
	    tic = std::chrono::steady_clock::now();
		
		gpuAssert( cudaMemcpyFromSymbol(ncontacts_h, ncontacts_global, 2*sizeof(uint)) );

		// cout << "d=" << 0 << " ncontacts_h[0]=" << ncontacts_h[0] << endl;
		// cout << "d=" << 1 << " ncontacts_h[1]=" << ncontacts_h[1] << endl;

		std::copy(&ncontacts_h[0],&ncontacts_h[0]+D,&out.counter[0]);

		for(uint d=0; d<2*D; ++d){
			thrust::copy(i_or_j_vecs[d]->begin(), i_or_j_vecs[d]->begin()+ncontacts_h[d/2], out.dout[d].data());
		    // gpuAssert( cudaMemcpy(out.dout[d].data(), indices_out_d[d], ncontacts_h[d/2] * sizeof(int), cudaMemcpyDeviceToHost) );
		}

		// cout << "Contacts found: " << ncontacts_h[0] << ", capacity=" << out.capacity() << endl;
		// for(uint i=0; i<ncontacts_h[0]; ++i)
		// 	cout << "find_distances_cuda: i=" << out.dout[0][i] << " j=" << out.dout[1][i] << endl;
		// cout << endl;

		
		for(uint d=0; d<2*D; ++d){
			delete i_or_j_vecs[d];
		}
		
   		tac = std::chrono::steady_clock::now();
   		us = std::chrono::duration_cast<std::chrono::microseconds>(tac-tic).count();
   		std::cout << "Time(CUDA: After the kernel launch): " << us << " us" << std::endl << std::endl;
	}
	
} // end-of-namespace dist

#endif // DIST_CUDA_IMPL
