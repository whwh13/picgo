/*
AUTHOR: G.A. Papoian, Date: Jan 6, 2019

A multithreaded versio of the code for both AVX and AVX2 calculations.
*/

#include <queue>
#include <mutex>


//#include "dist_moduleV2/dist_avx_par.h"

#ifdef SIMDBINDINGSEARCH

#include "dist_avx_par.h"

namespace medyan::dist {
	using namespace std;

	// xt::xarray<int> mat;

	PGrid griddim;
	PBlock blockdim;

	queue<pair <uint, uint>> block_queue;

	std::mutex mutex_block;
	// std::atomic_flag lock_block = ATOMIC_FLAG_INIT;
	std::array<std::atomic<uint>, 2> ncontacts_global_h;
	// std::array<uint, 4> ncontacts_global_h;

	uvec8_i zero_to_seven{0,1,2,3,4,5,6,7};
	extern uint _lut_avx2[256*8];

	// uvec8_i zero_to_eight{0,1,2,3,4,5,6,7};

	void probe_queue(){
		uint i=0;
		while (!block_queue.empty())
		{
			auto this_block = block_queue.front();
			std::cout << "block # " << i << " block = " << this_block.first << " " << this_block.second << endl;
			block_queue.pop();
			++i;
		}
	}

	template <uint D, bool SELF>
	inline void dist_simd_ij_warp(PBlock block, uint tx, uint ty, uint N1, uint N2,
	                              dOut<D,SELF> &out, uvec8_i &c1_i_ind, uvec8_i &c1_finfo, uvec8_f &c1_vxi,
	                              uvec8_f &c1_vyi, uvec8_f &c1_vzi, Coords &c2){

		// mat(py,px) += 1;

		constexpr uint simd_size = get_simd_size(tag_simd<simd_avx_par,float>());
		constexpr uint simd_and_stride_length_x = PSTRIDE*simd_size;

		// cout << "simd_and_stride_length_x=" << simd_and_stride_length_x << endl;

		uint bid = block.x + block.y * griddim.x; // block ID (global)

		uint px = (block.x * PBLOCKDIMX + tx) * simd_and_stride_length_x;
		uint py = block.y * PBLOCKDIMY + ty;

		constexpr uint size_shared_arr = simd_and_stride_length_x; // divide by 4 or 8 or something

		int shared_i_out[D*size_shared_arr]; // fix it!
		int shared_j_out[D*size_shared_arr];

		// for(uint d=0; d<D; ++d){
		// 	fill(&shared_i_out[d][0], &shared_i_out[d][size_shared_arr], 0);
		// 	fill(&shared_j_out[d][0], &shared_j_out[d][size_shared_arr], 0);
		// }


		uint ncontacts_shared[2] = {0,0};

		for(uint stride=0; stride<PSTRIDE; ++stride){

			uint sx = px + stride*simd_size;
			uint sy = py;

			if(sy>=N1 || (sx+8)>N2 || (SELF && (sx+8)<sy))
				continue;

			uvec8_f c2_vxj(&c2.x[sx]);
			uvec8_f c2_vyj(&c2.y[sx]);
			uvec8_f c2_vzj(&c2.z[sx]);
			uvec8_i c2_finfo(&c2.filinfo[sx]);

			//uvec8_i c1_i_ind;

			// cout << "dist_simd_ij(...)" << endl;
			// for(uint i=0; i<8; ++i)
			// 	cout << "px=" << px << " py=" << py << " c1_vxi=" << c1_vxi[i] << " c2_vxj=" << c2_vxj[i] << endl;
			// cout << endl;

			uvec8_f vdx = c2_vxj-c1_vxi;
			uvec8_f vdy = c2_vyj-c1_vyi;
			uvec8_f vdz = c2_vzj-c1_vzi;
			uvec8_f vsum = vdx*vdx + vdy*vdy + vdz*vdz;

			//compute finfo
			uvec8_i diff = c2_finfo - c1_finfo;
			diff = UME::SIMD::FUNCTIONS::abs(diff);
			uvec8_i &reffinfo = out.maxneighbors[0];

			#ifdef __AVX2__
			__m256i vcond2 = (_mm256_cmpgt_epi32(diff.mVec,reffinfo.mVec));
			#else

			/* AVX does not support 256bit packed integer comparisons. So cast the data to
		 * float and compare*/

		/*vcond2 =  compare( convert_to_float(diff.mVec), convert_to_float (reffinfo.mVec)
		 * ) */
		__m256 vcond2 = _mm256_cmp_ps(
				 		_mm256_cvtepi32_ps(diff.mVec),
				 		_mm256_cvtepi32_ps(reffinfo.mVec),
				 		_CMP_GT_OQ);
			#endif



			// cout << "dist_simd_ij(...)" << endl;
			// for(uint i=0; i<8; ++i)
			// 	cout << "px=" << px << " py=" << py << " c1_vxi=" << c1_vxi[i] << " c2_vxj=" << c2_vxj[i] << " vsum=" << vsum[i] << endl;
			// cout << endl;

			for(uint d=0; d<D; ++d){

				// cout << "\n\n\ndist_simd_ij_warp bid=" << bid << " block.x=" << block.x << " block.y=" << block.y
				// 	                                   << " tx=" << tx << " ty=" << ty
				// 									   << " px=" << px << " py=" << py
				// 									   << " sx=" << sx << " sy=" << sy << " d=" << d << endl;

				uvec8_f &vdl = out.v_dt[2*d];
				uvec8_f &vdh = out.v_dt[2*d+1];

				// The commented line below is preferable, but DTII gcc 6.1 and UME::SIMD interact badly, slowing down
				// the program by 4 fold (because UMD::SIMD starts using scalar emulation, for some reason).
				// On my Macbook Pro it is fine, even with gcc (8.x)

				// auto vcond = (vsum > vdl) && (vsum < vdh);
//				__m256i vcond = _mm256_castps_si256(_mm256_and_ps(_mm256_cmp_ps(vsum.mVec,vdl.mVec,_CMP_GT_OQ),
//				                                _mm256_cmp_ps(vdh.mVec,vsum.mVec,_CMP_GT_OQ)));

				__m256i vcond;

#ifdef __AVX2__
				__m256i vcond3 = _mm256_castps_si256(_mm256_and_ps(
						/*Cndn2*/_mm256_cmp_ps(vsum.mVec,vdl.mVec,_CMP_GT_OQ),
						/*Cndn3*/_mm256_cmp_ps(vdh.mVec,vsum.mVec,_CMP_GT_OQ)));

				vcond = _mm256_and_si256(vcond2, vcond3);
#else
				/* AVX does not support bitwise AND of 256bit integers. So do not cast
				 * vcond3 as an integer*/
			__m256 vcond3 = _mm256_and_ps(
					/*Cndn2*/_mm256_cmp_ps(vsum.mVec,vdl.mVec,_CMP_GT_OQ),
					/*Cndn3*/_mm256_cmp_ps(vdh.mVec,vsum.mVec,_CMP_GT_OQ));

			vcond = _mm256_castps_si256(_mm256_and_ps(vcond2,vcond3));

#endif


				int icond = _mm256_movemask_ps(_mm256_castsi256_ps(vcond));

				if(!icond)
					continue;

				if(SELF && sx <= sy){
					// cout << "Crossed a diagonal: sy=" << sy << " sx=" << sx << endl;
					// print_bits(icond);
					uvec8_i sindex_x(sx);
					uvec8_i sindex_y(sy);
					sindex_x = sindex_x + zero_to_seven;
					// print_m256i_as_int32(sindex_x.mVec);
					auto mask_within_bounds = sindex_x > sindex_y;
					int ic2 = _mm256_movemask_ps(_mm256_castsi256_ps(mask_within_bounds.mMask));
					// print_bits(ic2);
					vcond = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(vcond), _mm256_castsi256_ps(mask_within_bounds.mMask)));
					int ic3 = _mm256_movemask_ps(_mm256_castsi256_ps(vcond));
					// print_bits(ic3);
					// cout << endl << endl;
					icond = ic3;
				}


				uint ncontacts_warp = 0;
				uvec8_i c2_j_ind(&c2.indices[sx]);

				// cout << "c2_j_ind: ";
				// for(uint k=0; k<8; ++k){
				// 	cout << c2_j_ind[k] << " ";
				// }
				// cout << endl;
				// cout << "c1_i_ind: ";
				// for(uint k=0; k<8; ++k){
				// 	cout << c1_i_ind[k] << " ";
				// }
				// cout << endl;
#ifdef __AVX2__
				UME::SIMD::SIMDSwizzle<8> vmask(&_lut_avx2[8*icond]);
				c2_j_ind.swizzlea(vmask);
#else
				swizzlea_8x32_avx_impl(c2_j_ind, icond);
#endif
				// for(uint q=0; q<8; ++q)
				// 	c1_i_ind[q]=999;
				c1_i_ind.store(&(shared_i_out[size_shared_arr*d+ncontacts_shared[d]]));
				c2_j_ind.store(&(shared_j_out[size_shared_arr*d+ncontacts_shared[d]]));

				ncontacts_warp = _mm_popcnt_u32(icond);
				ncontacts_shared[d] += ncontacts_warp;

				// uint ncontacts_global_h_current = ncontacts_global_h[d];
				// ncontacts_global_h[d]+=ncontacts_warp;

				// cout << "stride=" << stride << " px=" << px << " sx=" << sx << " sy=" << sy << " d=" << d << " "; print_bits(icond);
				// cout << "ncontacts_warp=" << ncontacts_warp << " ncontacts_shared[d]=" << ncontacts_shared[d] << endl;
				// cout << "c2_j_ind: ";
				// for(uint k=0; k<8; ++k){
				// 	cout << c2_j_ind[k] << " ";
				// }
				// cout << endl;
				// cout << "c1_i_ind: ";
				// for(uint k=0; k<8; ++k){
				// 	cout << c1_i_ind[k] << " ";
				// }
				// cout << endl;
				// for(uint k=0; k<size_shared_arr; ++k){
				// 	cout << shared_j_out[d][k] << " ";
				// }
				// cout << endl;
				// for(uint k=0; k<size_shared_arr; ++k){
				// 	cout << shared_i_out[d][k] << " ";
				// }
				// cout << endl;
			}
		}

		for(uint d=0; d<D; ++d){

			if(ncontacts_shared[d]==0)
				continue;

			uint ncontacts_global_h_current = std::atomic_fetch_add(&ncontacts_global_h[d], ncontacts_shared[d]);

			// memcpy (&(out.dout[2*d])[ncontacts_global_h_current], &(shared_i_out[d][0]), ncontacts_shared[d]*sizeof(int) );
			// memcpy (&(out.dout[2*d+1])[ncontacts_global_h_current], &(shared_j_out[d][0]), ncontacts_shared[d]*sizeof(int) );

			std::copy(&(shared_i_out[size_shared_arr*d]), &shared_i_out[size_shared_arr*d+ncontacts_shared[d]], &(out.dout[2*d])[ncontacts_global_h_current]);
			std::copy(&(shared_j_out[size_shared_arr*d]), &shared_j_out[size_shared_arr*d+ncontacts_shared[d]], &(out.dout[2*d+1])[ncontacts_global_h_current]);
		}

		// uint ncontacts_global_h_current = std::atomic_fetch_add(&ncontacts_global_h[d], ncontacts_warp);
		// c1_i_ind.store(&);
		// c2_j_ind.store(&(out.dout[2*d+1])[ncontacts_global_h_current]);	

	}

	template <uint D, bool SELF>
	void kernel_dist_simd_avx_par(dOut<D,SELF> &out, Coords &c1, Coords &c2){

		while(true){
			mutex_block.lock();
			if(block_queue.empty()){
				mutex_block.unlock();
				return;
			}
			auto nextblock = block_queue.front();
			block_queue.pop();
			mutex_block.unlock();


			PBlock block{nextblock.first,nextblock.second};

			uint N1 = c1.size();
			uint N2 = c2.size();

			for(uint tx=0; tx<PBLOCKDIMX; ++tx){
				for(uint ty=0; ty<PBLOCKDIMY; ++ty){
					// mat(px,py) =1;

					uint py = block.y * PBLOCKDIMY + ty;

					uvec8_f c1_vxi(c1.x[py]);
					uvec8_f c1_vyi(c1.y[py]);
					uvec8_f c1_vzi(c1.z[py]);
					uvec8_i c1_i_ind(c1.indices[py]);
					uvec8_i c1_finfo(c1.filinfo[py]);

					dist_simd_ij_warp(block, tx, ty, N1, N2, out, c1_i_ind, c1_finfo,
					                  c1_vxi, c1_vyi, c1_vzi, c2);
					// cout << "blockx=" << blockx << " blocky=" << blocky << " x=" << x << " px=" << px << " y=" << y << " py=" << py << endl;
				}
			}
			// cout << endl;
		}

	}

	template void kernel_dist_simd_avx_par(dOut<1U,true> &out, Coords &c1, Coords &c2);
	template void kernel_dist_simd_avx_par(dOut<1U,false> &out, Coords &c1, Coords &c2);
	template void kernel_dist_simd_avx_par(dOut<2U,true> &out, Coords &c1, Coords &c2);
	template void kernel_dist_simd_avx_par(dOut<2U,false> &out, Coords &c1, Coords &c2);
}

#endif