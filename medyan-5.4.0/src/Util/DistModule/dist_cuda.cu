#ifdef CUDAACCL

//#include "dist_moduleV2/dist_cuda.cuh"
#include "dist_cuda.cuh"

namespace medyan::dist {

	template void find_distances_cuda(dOut<1,true> &out, Coords &c1, Coords &c2);
	template void find_distances_cuda(dOut<1,false> &out, Coords &c1, Coords &c2);
	template void find_distances_cuda(dOut<2,true> &out, Coords &c1, Coords &c2);
	template void find_distances_cuda(dOut<2,false> &out, Coords &c1, Coords &c2);
	
} // end of namespace dist

#endif // CUDAACCL
