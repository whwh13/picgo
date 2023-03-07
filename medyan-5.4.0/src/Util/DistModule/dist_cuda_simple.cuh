#include "dist_moduleV2/dist_cuda_common.cuh"

namespace medyan::dist {

	// extern __device__ uint ncontacts_global[2];
	
	template <uint D, bool Self>
    __global__ void find_distances_kernel_simple(uint N1, uint N2, uint Ni, float *dlim,
                               float *x1, float *y1, float *z1, float *x2, float *y2, float *z2,
                               int *indices1, int *indices2,
							   int **i_or_j_vecs)
	{	
		// printf("blockIdx.x=%d blockIdx.y=%d threadIdx.x=%d threadIdx.y=%d cond=%d\n",
		//         blockIdx.x,   blockIdx.y,   threadIdx.x,   threadIdx.y, (Self && blockIdx.x< blockDim.y) );
		// primitive cell
		if(Self && blockIdx.x < blockIdx.y)
			return;
		
	    uint px = blockIdx.x * blockDim.x + threadIdx.x; 
		uint py = blockIdx.y * blockDim.y + threadIdx.y; 
		
		bool self_ij_toexclude = false;
		
		if(Self && px <= py)
			self_ij_toexclude = true;
						
		float dx = x1[py]-x2[px];
		float dy = y1[py]-y2[px];
		float dz = z1[py]-z2[px];

		float dist = dx*dx + dy*dy + dz*dz;
				
		bool cond[D];
		uint ncontacts_thread[D];
		
		// printf("dlim: %f %f %f %f\n", dlim[0], dlim[1], dlim[2], dlim[3]);

		for(uint d=0; d<D; ++d){
			cond[d] = (dist > dlim[2*d] && dist < dlim[2*d+1]) && (!self_ij_toexclude);
			if(cond[d]){
				ncontacts_thread[d] =  atomicAdd(&ncontacts_global[d], 1);
				i_or_j_vecs[2*d][ncontacts_thread[d]] = indices2[py];
				i_or_j_vecs[2*d+1][ncontacts_thread[d]] = indices1[px];
				// printf("blockIdx.x=%d blockIdx.y=%d px=%d py=%d d=%d cond[d]=%d ncontacts_thread[d]=%d\n",
				//         blockIdx.x,   blockIdx.y,   px,   py,   d,   cond[d],   ncontacts_thread[d]);
			}
		}
	}
	

	
} // end-of-namespace dist