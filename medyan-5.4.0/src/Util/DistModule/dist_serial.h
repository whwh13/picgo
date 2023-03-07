/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Simple serial implementation for finding contacts.
*/

#ifndef DIST_SERIAL
#define DIST_SERIAL

/*#include "dist_moduleV2/dist_common.h"
#include "dist_moduleV2/dist_out.h"
#include "dist_moduleV2/dist_coords.h"*/

#include "dist_common.h"
#include "dist_out.h"
#include "dist_coords.h"

namespace medyan::dist {

	// Scalary implementation that is used in part by the SIMD code too (for corner cases)
	template <uint D, bool SELF>
	inline void dist_scalar_ij(dOut<D,SELF> &out, const Coords &c1, const Coords &c2, uint i, uint j){
	
		float dx = c2.x[j]-c1.x[i];
		float dy = c2.y[j]-c1.y[i];
		float dz = c2.z[j]-c1.z[i];
		int diff = c1.filinfo[i] - c2.filinfo[j];
		diff = abs(diff);
//		diff = 3;

		float dist_sq = dx*dx + dy*dy + dz*dz;

		for(uint d=0; d<D; ++d){
			float d_l = out.dt[2*d];
			float d_h = out.dt[2*d+1];
			uint &counter(out.counter[d]);
            if(dist_sq > d_l && dist_sq<d_h && diff > ChemParams::minCylinderDistanceSameFilament){
//			if(dist_sq > d_l && dist_sq<d_h){
				(out.dout[2*d])[counter] = c1.indices[i];
				(out.dout[2*d+1])[counter] = c2.indices[j];
				++counter;
			}
		}	
	}

	// The main scalar outer double loop.
	template <uint D, bool SELF>
	inline void find_distances(dOut<D,SELF> &out, Coords &c, tag_simd<simd_no,float> tag){
		uint N = c.size();
		for(uint i=0; i<N; ++i){
			for(uint j=(i+1); j<N; ++j){
				dist_scalar_ij(out, c, c, i, j);			
				// float d_l = out.dt[0];
				// float d_h = out.dt[1];
				//
				// float dx = c.x[j]-c.x[i];
				// float dy = c.y[j]-c.y[i];
				// float dz = c.z[j]-c.z[i];
				//
				// float dist_sq = dx*dx + dy*dy + dz*dz;
				//
				// if(dist_sq > d_l && dist_sq<d_h){
				// 	cout << "dist_scalar_ij: i=" << i << " j=" << j << " d=" << dist_sq << endl;
				// }
			}
		}
	}
	
	// The main scalar outer double loop.
	template <uint D, bool SELF>
	inline void find_distances(dOut<D,SELF> &out, Coords &c1, Coords &c2, tag_simd<simd_no,float> tag){
		uint N1 = c1.size();
		uint N2 = c2.size();
		for(uint i=0; i<N1; ++i){
			for(uint j=0; j<N2; ++j){
				dist_scalar_ij(out, c1, c2, i, j);
			}
		}
	}

} // end-of-namespace dist

#endif // DIST_SERIAL
