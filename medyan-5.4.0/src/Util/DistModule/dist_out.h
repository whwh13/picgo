/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Both input distance thresholds and the output arrays that hold found contact indices are here. 
*/

#ifndef DIST_OUT
#define DIST_OUT

#include <initializer_list>
#include <algorithm>
#include <array>
#include <vector>

#include "dist_common.h"
#include "dist_coords.h"
#include "SysParams.h"

/*#include "dist_moduleV2/dist_common.h"
#include "dist_moduleV2/dist_coords.h"*/

//#ifdef __CUDACC__
//#include <thrust/host_vector.h>
//#include <thrust/system/cuda/experimental/pinned_allocator.h>
//#endif


namespace medyan::dist {
	
	// This class holds both found contacts data and also the input specifying the distance thresholds.
	// Template parameter D is the number of threshold - i.e. D=2 if we need to apply 2 distance predicates.
	template <uint D, bool SELF=true>
	struct dOut {
		// variables
        alignas(32) std::array<float,2*D> dt; // distance thresholds' vector - lower/upper
        alignas(32) std::array<uvec8_f, 2*D> v_dt; // UME::SIMD versions of dt: i.e. dt elements broadcasted to __m256 registers and wrapped
		alignas(32) std::array<uvec8_i, 1> maxneighbors;
		int N1, N2;

		std::array<uint,D> counter; // number of found contacts
		
	  //#ifdef __CUDACC__
		//	std::array<thrust::host_vector<int,thrust::system::cuda::experimental::pinned_allocator<int>>,2*D> dout;
	  //#else
		std::array<std::vector<int>,2*D> dout; // output vectors for each threshold window
	  //#endif
	
		void reset_counters()
		{
			for(uint i=0; i<D; ++i)
				counter[i] =  0;
		}
		
		void resize_output_indices(){
			uint max_pairs;
			
			if(SELF==true){
				max_pairs = N1*(N1-1)/2;
			}
			else{
				max_pairs = N1*N2;				
			}
			
			for(auto &xv: dout)
				xv.resize(max_pairs);
		}
		
		void init_dout(uint n1, int n2, std::initializer_list<float> vals) {
			N1 = n1;
			N2 = n2;		
			reset_counters();
		
			std::copy(vals.begin(),vals.end(),dt.begin());


            maxneighbors[0] = ChemParams::minCylinderDistanceSameFilament;


			for(uint i=0; i<2*D; ++i){
				v_dt[i] = dt[i];
			}
			
			resize_output_indices();
		}
		
		void init_dout(uint n, std::initializer_list<float> vals) {
			init_dout(n, n, vals);
		}
	
		dOut()
		{}

		dOut(uint n, std::initializer_list<float> vals) {
			init_dout(n, vals);
		}
		
		dOut(uint n1, uint n2, std::initializer_list<float> vals) {
			init_dout(n1, n2, vals);
		}
	
		uint size(uint index) const {return counter[index];}
		uint dim() const {return D;}
		int capacity() const {return dout[0].size();}

	};

} // end-of-namespace dist

#endif // DIST_OUT
