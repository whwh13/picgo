/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

A simple testing suite. Mainly checks the SIMD version results again the serial ones.
*/

#include <iostream>

#include <catch2/catch.hpp>

#include "Util/DistModule/dist_coords.h"
#include "Util/DistModule/dist_driver.h"
#include "Util/DistModule/dist_out.h"

TEST_CASE("Dist module", "[Dist]") {
    using namespace std;
    using namespace medyan;
	using namespace medyan::dist;

	const auto pack_two_ints = [](int32_t high, int32_t low) -> int32_t {
		return (high << 16) | low; // pack lower 16 bits of each integer into a new integer; low, high < 65536 (i.e N < 65536)
	};

	const auto int_low = [](int32_t i) -> int32_t {
		int mask_l = (1<<16)-1;
		return i & mask_l;
	};

	const auto print_int_low = [&](int32_t i) {
		cout << int_low(i) << " ";
	};

    const auto int_high = [](int32_t i) {
        int32_t mask_h = ((1<<16)-1) << 16;
        return ((i & mask_h) >> 16);
    };

	const auto print_int_high = [&](int32_t i) {
		cout << int_high(i) << " ";
	};

	const auto convert_dout_unique = [&](auto &out, unsigned d=0) {

		const auto ncontacts = out.counter[d];
	
		vector<int> vout(ncontacts);
        for(int k=0; k < ncontacts; ++k){
			vout[k] = pack_two_ints(out.dout[2*d][k], out.dout[2*d+1][k]);
		}
	
		// int mask_l = (1<<16)-1;
		// int mask_h = ((1<<16)-1) << 16;
		//
		// print_bits(mask_l,32);
		// print_bits(mask_h,32);
	
		// for(uint k=1000; k<1010; ++k){
		// 	cout << "i, j = [" << out.dout[0][k] << ":" << out.dout[1][k] << "], vout=" << vout[k] << endl;
		// 	cout << "i from vout: " << (vout[k] & mask_l) << endl;
		// 	cout << "j from vout: " << ((vout[k] & mask_h) >> 16) << endl;
		// }
	
		sort(vout.begin(), vout.end());
		return vout;
	};

    // vout1 and vout2 must be sorted lists.
    const auto test_compare_douts = [](const vector<int> &vout1, const vector<int> &vout2) {

		size_t n_max = max(vout1.size(), vout2.size());
		size_t res_size = vout1.size() + vout2.size();
		vector<int> res(res_size);
		auto newend = std::set_symmetric_difference(vout1.begin(),vout1.end(),
									   vout2.begin(),vout2.end(),
									   res.begin());
	    size_t n_different = std::distance(res.begin(), newend);

		double fraction_failure = (double)n_different/n_max;
		return 1-fraction_failure;
	};

    const auto test_algo1 = [&](Coords &coords, auto &tag, string algo_name) {
		INFO("test_algo1(...), with " << algo_name);

		const double min_fract_succ = 0.99;
		float d_l{5.0f}, d_h{12.0f};
	
		const auto N = coords.size();

		dOut<1> out_serial(N, {d_l,d_h});
		find_distances(out_serial,coords,tag_simd<simd_no,float>());	
		vector<int> vout_serial = convert_dout_unique(out_serial);
		
		dOut<1> out(N, {d_l,d_h});
		find_distances(out,coords,tag);	
		vector<int> vout = convert_dout_unique(out);
		const auto fraction_success = test_compare_douts(vout,vout_serial);
		bool success = fraction_success > min_fract_succ;
	
        INFO("serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0]);
        if(!success) {
			FAIL_CHECK("test_algo(...): " << algo_name << " did not pass. fraction_success=" << fraction_success);
			// for(uint k=0; k<out_serial.counter[0]; ++k){
			for(uint k=0; k<20; ++k) {
				// if(int_high(vout_serial[k])!=0 && int_high(vout[k])!=0)
				// 	continue;
				cout << "k=" << k << " serial i and j: ";
				print_int_high(vout_serial[k]); print_int_low(vout_serial[k]);
				cout << " and " << algo_name << " i and j: ";
				print_int_high(vout[k]); print_int_low(vout[k]);
				cout << endl;
				// cout << "vout_serial=" << vout_serial[k] << endl;
				// cout << "vout=" << vout[k] << endl;
			}
		}
    };
	
    const auto test_algo2 = [&](Coords &coords, auto &tag, string algo_name) {
		INFO("test_algo2(...), with " << algo_name);

        const double min_fract_succ = 0.99;
	
		const auto N = coords.size();

		dOut<2> out_serial(N, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out_serial,coords,tag_simd<simd_no,float>());
		dOut<2> out(N, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out,coords,tag);	
	
		bool success = true;	
		for(uint d=0; d<2; ++d){
			vector<int> vout_serial = convert_dout_unique(out_serial,d);
	
			vector<int> vout = convert_dout_unique(out,d);
			const auto fraction_success = test_compare_douts(vout,vout_serial);
			success = fraction_success > min_fract_succ;

            INFO("serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0]);
            if(!success) {
				FAIL_CHECK("test_algo2(...): " << algo_name << " d=" << d << " did not pass. fraction_success=" << fraction_success);
				for(uint k=0; k<20; ++k){
					cout << "serial i and j: ";
					print_int_low(vout_serial[k]); print_int_high(vout_serial[k]);
					cout << " and " << algo_name << " i and j: ";
					print_int_low(vout[k]); print_int_high(vout[k]);
					cout << endl;
					// cout << "vout_serial=" << vout_serial[k] << endl;
					// cout << "vout=" << vout[k] << endl;
				}
			}
		}
	};

	const auto test_algo1_betwn_comps = [&](Coords &c1, Coords &c2, auto &tag, string algo_name) {
		INFO("test_algo1_betwn_comps(...), with " << algo_name);
		const double min_fract_succ = 0.999;
		float d_l{5.0f}, d_h{12.0f};
	
		const auto N1 = c1.size();
		const auto N2 = c2.size();

		// uint N = std::max({N1,N2});
	
		dOut<1,false> out_serial(N1, N2, {d_l,d_h});
		find_distances(out_serial,c1,c2,tag_simd<simd_no,float>());	
		vector<int> vout_serial = convert_dout_unique(out_serial);
		
		dOut<1,false> out(N1, N2, {d_l,d_h});
		find_distances(out,c1,c2,tag);	
		vector<int> vout = convert_dout_unique(out);
		const auto fraction_success = test_compare_douts(vout,vout_serial);
		bool success = fraction_success > min_fract_succ;
	
        INFO("serial ncontacts=" << out_serial.counter[0] << " " << algo_name << " ncontacts=" << out.counter[0]);
		if(!success) {
			FAIL_CHECK("test_algo1_betwn_comps(...): " << algo_name << " did not pass. fraction_success=" << fraction_success);
			for(uint k=0; k<20; ++k){
				cout << "serial i and j: ";
				print_int_low(vout_serial[k]); print_int_high(vout_serial[k]);
				cout << " and " << algo_name << " i and j: ";
				print_int_low(vout[k]); print_int_high(vout[k]);
				cout << endl;
				// cout << "vout_serial=" << vout_serial[k] << endl;
				// cout << "vout=" << vout[k] << endl;
			}
		}
	};
	
    const auto test_algo2_betwn_comps = [&](Coords &c1, Coords &c2, auto &tag, string algo_name) {
        INFO("\ntest_algo2_betwn_comps(...), with " << algo_name);
		
		const double min_fract_succ = 0.999;
	
		const auto N1 = c1.size();
		const auto N2 = c2.size();		
	
		dOut<2,false> out_serial(N1, N2, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out_serial,c1,c2,tag_simd<simd_no,float>());
		dOut<2,false> out(N1, N2, {5.0f, 12.0f, 6.0f, 13.0f});
		find_distances(out,c1,c2,tag);	
	
		bool success = true;	
		for(uint d=0; d<2; ++d){
			vector<int> vout_serial = convert_dout_unique(out_serial,d);
	
			vector<int> vout = convert_dout_unique(out,d);
			const auto fraction_success = test_compare_douts(vout,vout_serial);
			success = fraction_success > min_fract_succ;
	
            INFO("serial ncontacts=" << out_serial.counter[d] << " " << algo_name << " ncontacts=" << out.counter[d]);
			if(!success) {
				FAIL_CHECK("test_algo2_betwn_comps(...): " << algo_name << " d=" << d << " did not pass. fraction_success=" << fraction_success);
				for(uint k=0; k<20; ++k){
					cout << "serial i and j: ";
					print_int_low(vout_serial[k]); print_int_high(vout_serial[k]);
					cout << " and " << algo_name << " i and j: ";
					print_int_low(vout[k]); print_int_high(vout[k]);
					cout << endl;
					// cout << "vout_serial=" << vout_serial[k] << endl;
					// cout << "vout=" << vout[k] << endl;
				}
			}
		}
	};



    const int n1 = 600;
    const int n2 = 560;

    Coords c1(n1), c2(n2);

    tag_simd<simd_no,  float>           t_serial;
    tag_simd<simd_avx, float>           t_avx;		
    tag_simd<dist::simd_avx_par,float>  t_avx_par;
#ifdef __CUDACC__
    tag_simd<cuda,     float>           t_cuda;
#endif

    {
        INFO("single comparisions");
        test_algo1(c1, t_avx, "AVX");
        // test_algo1(c1, t_avx_par, "AVX-PARALLEL");
#ifdef __CUDACC__
        test_algo1(c1, t_cuda, "CUDA");
#endif
    }

    {
        INFO("two comparisions");
        test_algo2(c1, t_avx, "AVX");
        // test_algo2(c1, t_avx_par, "AVX-PARALLEL");
#ifdef __CUDACC__
        test_algo2(c1, t_cuda, "CUDA");
#endif
    }

    {
        INFO("two-compartment functions:");
#ifdef __CUDACC__
        test_algo1_betwn_comps(c1, c2, t_cuda, "CUDA");
#endif
        test_algo1_betwn_comps(c1, c2, t_avx, "AVX");
        // test_algo1_betwn_comps(c1, c2, t_avx_par, "AVX-PARALLEL");
    }

    {
        INFO("two-compartment functions with two comparisions:");
        test_algo2_betwn_comps(c1, c2, t_avx, "AVX");
        // test_algo2_betwn_comps(c1, c2, t_avx_par, "AVX-PARALLEL");
#ifdef __CUDACC__
        test_algo2_betwn_comps(c1, c2, t_cuda, "CUDA");
#endif
    }

}
