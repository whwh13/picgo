#ifdef SIMDBINDINGSEARCH
#include <random>
#include <iostream>
#include <chrono>

/*#include "dist_moduleV2/dist_coords.h"
#include "dist_moduleV2/dist_driver.h"
#include "dist_moduleV2/dist_out.h"*/

#include "dist_coords.h"
#include "dist_driver.h"
#include "dist_out.h"

namespace medyan::dist {

	using namespace std;

	void benchmark_dist_module(uint N1, uint N2)
	{
		// const uint N = 6000; // number of sites/atoms
		Coords c1(N1), c2(N2);


		// cout << out_serial.dt << endl;
		// cout << out_serial.dout << endl;

		auto tic = std::chrono::steady_clock::now();
		dOut<1> out_serial1(N1, {5.0f,12.0f});
		find_distances(out_serial1,c1,tag_simd<simd_no,float>());
		report_contact_stats(out_serial1);
		auto tac = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SERIAL1): " << ms << " ms" << std::endl << std::endl;
	
		dOut<2> out_serial2(N1, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_serial2,c1,tag_simd<simd_no,float>());
		report_contact_stats(out_serial2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SERIAL2): " << ms << " ms" << std::endl << std::endl;

		dOut<1> out_simd_avx_float_1(N1, {5.0f,12.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_simd_avx_float_1,c1,tag_simd<simd_avx,float>());
		report_contact_stats(out_simd_avx_float_1);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SIMD-AVX-F1): " << ms << " ms" << std::endl << std::endl;
		
		dOut<1> out_simd_avx_par_float_1(N1, {5.0f,12.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_simd_avx_par_float_1,c1,tag_simd<simd_avx_par,float>());
		report_contact_stats(out_simd_avx_par_float_1);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SIMD-AVX-PAR-F1): " << ms << " ms" << std::endl << std::endl;

		dOut<2> out_simd_avx_float_2(N1, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_simd_avx_float_2,c1,tag_simd<simd_avx,float>());
		report_contact_stats(out_simd_avx_float_2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SIMD-AVX-F2): " << ms << " ms" << std::endl << std::endl;
		
		dOut<2> out_simd_avx_par_float_2(N1, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_simd_avx_par_float_2,c1,tag_simd<simd_avx_par,float>());
		report_contact_stats(out_simd_avx_par_float_2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SIMD-AVX-PAR-F2): " << ms << " ms" << std::endl << std::endl;

#ifdef __CUDACC__
		dOut<1> out_cuda_float_1(N1, {5.0f,12.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_cuda_float_1,c1,tag_simd<cuda,float>());
		report_contact_stats(out_cuda_float_1);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(CUDA-F1): " << ms << " ms" << std::endl << std::endl;

		dOut<2> out_cuda_float_2(N1, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_cuda_float_2,c1,tag_simd<cuda,float>());
		report_contact_stats(out_cuda_float_2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(CUDA-F2): " << ms << " ms" << std::endl << std::endl;
#endif

		// Assuming that only 1/2 of each compartment will be used;
		Coords c3(N1/2), c4(N2/2);
		
		dOut<2> out_2cmp_simd_avx_float_2(N1, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_2cmp_simd_avx_float_2,c3,c4,tag_simd<simd_avx,float>());
		report_contact_stats(out_2cmp_simd_avx_float_2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SIMD-CMP2-AVX-F2): " << ms << " ms" << std::endl << std::endl;
		
		dOut<2> out_2cmp_simd_avx_par_float_2(N1, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_2cmp_simd_avx_par_float_2,c3,c4,tag_simd<simd_avx_par,float>());
		report_contact_stats(out_2cmp_simd_avx_par_float_2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SIMD-CMP2-AVX-PAR-F2): " << ms << " ms" << std::endl << std::endl;
		
#ifdef __CUDACC__
		dOut<1> out_2cmp_cuda_float(N1/2, {5.0f, 12.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_2cmp_cuda_float,c3,c4,tag_simd<cuda,float>());
		report_contact_stats(out_2cmp_cuda_float);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(CUDA-CMP2-F1): " << ms << " ms" << std::endl << std::endl;

		dOut<2> out_2cmp_cuda_float_2(N1/2, {5.0f, 12.0f, 6.0f, 13.0f});
		tic = std::chrono::steady_clock::now();
		find_distances(out_2cmp_cuda_float_2,c3,c4,tag_simd<cuda,float>());
		report_contact_stats(out_2cmp_cuda_float_2);
		tac = std::chrono::steady_clock::now();
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(CUDA-CMP2-F2): " << ms << " ms" << std::endl << std::endl;
#endif
	}

} // end-of-namespace dist
#endif
