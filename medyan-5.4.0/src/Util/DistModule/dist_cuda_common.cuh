#ifndef DIST_CUDA_COMMON
#define DIST_CUDA_COMMON

namespace medyan::dist {

#define BLOCKX 16
#define BLOCKY 16
	
	__device__ uint ncontacts_global[4];
	
}

#endif // DIST_CUDA_COMMON
