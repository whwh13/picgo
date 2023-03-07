/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

Example functions for various scenarios how to run find_distances(...)
*/


#ifdef SIMDBINDINGSEARCH

#include <random>
#include <iostream>
#include <iterator>
#include <chrono>
#include <vector>
#include <algorithm>
#include <thread>


/*#include "dist_moduleV2/dist_coords.h"
#include "dist_moduleV2/dist_driver.h"
#include "dist_moduleV2/dist_out.h"*/

#include "dist_coords.h"
#include "dist_driver.h"
#include "dist_out.h"

namespace medyan::dist {

	using namespace std;

	void dist_example_simplest()
	{
		const uint N = 6000;
	
		std::random_device rd;
		std::mt19937 mt(rd());
	    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
		std::uniform_int_distribution<int> fil_f(0, 100);
		
		vector<double> x(N), y(N), z(N);
		vector<int> finfo(N);
			
		for(uint i=0; i<N; ++i){
			x[i] = dist_d(mt);
			y[i] = dist_d(mt);
			z[i] = dist_d(mt);
			finfo[i] = fil_f(mt);
		}
	
		vector<int> indices(N);
		uint q=0;
		generate(indices.begin(),indices.end(),[&q](){return q++;});
		std::shuffle(indices.begin(),indices.end(), mt);
	
		Coords coords(x,y,z,indices,finfo);
		dOut<2U> out(N, {5.0f, 12.0f, 6.0f, 13.0f});

		// Note: If the third argument is not provided, then the best SIMD (but non-parallel algo will be chosen)
		find_distances(out,coords); 
		
		report_contact_stats(out);
	}
	
	void dist_example_between_compartments()
	{
		uint N1 = 6000;
		uint N2 = 5600;
	
		std::random_device rd;
		std::mt19937 mt(rd());
	    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
		std::uniform_int_distribution<int> fil_f(0, 100);
		
		vector<double> x1(N1), y1(N1), z1(N1);
		vector<double> x2(N2), y2(N2), z2(N2);
		vector<int> finfo1(N1), finfo2(N2);
			
		generate(x1.begin(),x1.end(),[&](){return dist_d(mt);});
		generate(y1.begin(),y1.end(),[&](){return dist_d(mt);});
		generate(z1.begin(),z1.end(),[&](){return dist_d(mt);});
		
		generate(x2.begin(),x2.end(),[&](){return dist_d(mt);});
		generate(y2.begin(),y2.end(),[&](){return dist_d(mt);});
		generate(z2.begin(),z2.end(),[&](){return dist_d(mt);});

		generate(finfo1.begin(),finfo1.end(),[&](){return fil_f(mt);});

		generate(finfo2.begin(),finfo2.end(),[&](){return fil_f(mt);});
		
	
		vector<int> indices1(N1);
		uint q1=0;
		generate(indices1.begin(),indices1.end(),[&q1](){return q1++;});
		std::shuffle(indices1.begin(),indices1.end(), mt);
		
		vector<int> indices2(N2);
		uint q2=0;
		generate(indices2.begin(),indices2.end(),[&q2](){return q2++;});
		std::shuffle(indices2.begin(),indices2.end(), mt);
		
		// We are done with generating fake coordinates in two distinct compartments
		//The next few lines below do the distance calculations
	
		Coords c1(x1,y1,z1,indices1,finfo1);
		Coords c2(x2,y2,z2,indices2,finfo2);
		uint N = std::max({N1,N2});
		dOut<2U> out(N, {5.0f, 12.0f, 6.0f, 13.0f});

		// Note: If the third argument is not provided, then the best SIMD (but non-parallel algo will be chosen)
		find_distances(out,c1,c2); 
		
		report_contact_stats(out);
	}

	void dist_example_nonparallel()
	{
		const uint N = 6000;
	
		std::random_device rd;
		std::mt19937 mt(rd());
	    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
		std::uniform_int_distribution<int> fil_f(0, 100);
		
		vector<double> x(N), y(N), z(N);
		vector<int> finfo(N);
			
		for(uint i=0; i<N; ++i){
			x[i] = dist_d(mt);
			y[i] = dist_d(mt);
			z[i] = dist_d(mt);
			finfo[i] = fil_f(mt);
		}
	
		vector<int> indices(N);
		uint q=0;
		generate(indices.begin(),indices.end(),[&q](){return q++;});
		std::shuffle(indices.begin(),indices.end(), mt);
	
		Coords coords(x,y,z,indices,finfo);
		dOut<2U> out(N, {5.0f, 12.0f, 6.0f, 13.0f});
		auto tic = std::chrono::steady_clock::now();
		find_distances(out,coords,tag_simd<simd_avx,float>());
		auto tac = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(SERIAL-AVX-EXAMPLE): " << ms << " ms" << std::endl << std::endl;
		report_contact_stats(out);
	}

	void dist_example_avx_veccoord()
	{
		
		cout << "dist_example_avx_veccoord()" << endl;
		
		const uint ncomps = 16; // number of compartments
		const uint N = 6000;
		
		std::random_device rd;
		std::mt19937 mt(rd());
	    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
		std::uniform_int_distribution<int> fil_f(0, 100);
		
		vector<float> x1, y1, z1;
	    vector<int> idx1, finfo;
	    dist::dOut<2U> bspairsoutX;
	    bspairsoutX.init_dout(N,{5.0f, 12.0f, 6.0f, 13.0f});
	    dist::tag_simd<dist::simd_avx_par,  float>   t_avx_par;		
	
		vector<dist::Coords> veccoord(ncomps);
		
		for(uint c=0; c<ncomps; ++c){
			x1.resize(N); 
			y1.resize(N); 
			z1.resize(N); 
			idx1.resize(N);
			finfo.resize(N);

			
			for(uint i=0; i<N; ++i){
				x1[i] = dist_d(mt); 
				y1[i] = dist_d(mt); 
				z1[i] = dist_d(mt); 
				idx1[i] = i;
				finfo[i] = fil_f(mt);
			}
			
			veccoord[c].init_coords(x1, y1, z1, idx1,finfo);
		}
		
		cout << "comp 0: findnig interactions within a compartment" << endl;
        dist::find_distances(bspairsoutX, veccoord[0], t_avx_par);
		report_contact_stats(bspairsoutX);
		bspairsoutX.reset_counters();
		
		for(uint c=1; c<ncomps; ++c){

			cout << "comp " << c << " and " << c-1 << " finding interactions between the compartments" << endl;
            dist::find_distances(bspairsoutX, veccoord[c], veccoord[c-1], t_avx_par);
			report_contact_stats(bspairsoutX);
			bspairsoutX.reset_counters();
			
			cout << "comp " << c << " finding interactions within a compartment" << endl;
	        dist::find_distances(bspairsoutX, veccoord[0], t_avx_par);
			bspairsoutX.reset_counters();

		}
	}

	// void dist_example_parallel()
	// {
	// 	const uint ncomps = 128; // number of compartments
	// 	const uint nthreads = 16;
	// 	const uint N = 6000;
	//
	// 	std::random_device rd;
	// 	std::mt19937 mt(rd());
	//     std::uniform_real_distribution<float> dist_d(1.0, 10.0);
	//
	// 	array<vector<double>,ncomps> X;
	// 	array<vector<double>,ncomps> Y;
	// 	array<vector<double>,ncomps> Z;
	// 	array<vector<int>,ncomps> INDICES;
	//
	//
	// 	array<Coords,nthreads> COORDS; // note nthreads << ncomps - we will be reusing COORDS and OUTS
	// 	array<Coords,nthreads> COORDS2; // note nthreads << ncomps - we will be reusing COORDS and OUTS
	// 	array<dOut<2U>,nthreads> OUTS;
	//
	// 	auto tic = std::chrono::steady_clock::now();
	// 	auto tac = std::chrono::steady_clock::now();
	// 	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
	//
	//
	// 	for(uint c=0; c<ncomps; ++c){
	// 		X[c].resize(N);
	// 		Y[c].resize(N);
	// 		Z[c].resize(N);
	// 		for(uint i=0; i<N; ++i){
	// 			X[c][i] = dist_d(mt);
	// 			Y[c][i] = dist_d(mt);
	// 			Z[c][i] = dist_d(mt);
	// 		}
	//
	// 		INDICES[c].resize(N);
	// 		uint q=0;
	// 		generate(INDICES[c].begin(),INDICES[c].end(),[&q](){return q++;});
	// 		random_shuffle(INDICES[c].begin(),INDICES[c].end());
	// 	}
	//
	//
	// 	// We are done with generating fake coordinates in 128 compartments
	// 	//The next few lines below do the distance calculations WITHIN compartments
	//
	// 	tic = std::chrono::steady_clock::now();
	// 	// uint chunks = ncomps/nthreads; // in a production code handle corner cases - when not divisible
	// 	for(uint c=0; c<(ncomps-nthreads); c+=nthreads){
	//
	// 		std::vector<std::thread> threads;
	//
	// 		for(uint t=0; t<nthreads; ++t){
	// 			uint q = c+t;
	// 			if(c==0)
	// 				OUTS[t].init_dout(N, {5.0f, 12.0f, 6.0f, 13.0f}); // this is an expensive call: 30 ms per iteration
	// 			else
	// 				OUTS[t].reset_counters();
	//
	// 			COORDS[t].init_coords(X[q],Y[q],Z[q],INDICES[q]);
	// 			threads.push_back(thread(find_distances_2, ref(OUTS[t]), ref(COORDS[t])));
	// 		}
	// 		for(auto &t: threads)
	// 			t.join();
	//
	// 		for(auto &o: OUTS)
	// 			report_contact_stats(o);
	// 	}
	//
	// 	tac = std::chrono::steady_clock::now();
	// 	ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
	// 	std::cout << "Time(PARALLEL-AVX2-EXAMPLE): " << (double)ms/ncomps << " ms" << std::endl << std::endl;
	//
	//
	// 	//The next few lines below do the distance calculations BETWEEN neighboring compartments
	//
	// 	tic = std::chrono::steady_clock::now();
	// 	// uint chunks = ncomps/nthreads; // in a production code handle corner cases - when not divisible
	// 	for(uint c=0; c<(ncomps-nthreads)-1; c+=nthreads){
	//
	// 		std::vector<std::thread> threads;
	//
	// 		for(uint t=0; t<nthreads; ++t){
	// 			uint q = c+t;
	// 			if(c==0)
	// 				OUTS[t].init_dout(N, {5.0f, 12.0f, 6.0f, 13.0f}); // this is an expensive call: 30 ms per iteration
	// 			else
	// 				OUTS[t].reset_counters();
	//
	// 			COORDS[t].init_coords(X[q],Y[q],Z[q],INDICES[q]);
	// 			COORDS2[t].init_coords(X[q+1],Y[q+1],Z[q+1],INDICES[q+1]);
	// 			threads.push_back(thread(find_distances_cmp2_2, ref(OUTS[t]), ref(COORDS[t]), ref(COORDS2[t])));
	// 		}
	// 		for(auto &t: threads)
	// 			t.join();
	//
	// 		for(auto &o: OUTS)
	// 			report_contact_stats(o);
	// 	}
	//
	// 	tac = std::chrono::steady_clock::now();
	// 	ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
	// 	std::cout << "Time(PARALLEL-CMP2-AVX2-EXAMPLE): " << (double)ms/ncomps << " ms" << std::endl << std::endl;
	// }

#ifdef __CUDACC__	
	void dist_example_between_compartments_cuda()
	{
		uint N1 = 6000;
		uint N2 = 5600;
	
		std::random_device rd;
		std::mt19937 mt(rd());
	    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
		
		vector<double> x1(N1), y1(N1), z1(N1);
		vector<double> x2(N2), y2(N2), z2(N2);
			
		generate(x1.begin(),x1.end(),[&](){return dist_d(mt);});
		generate(y1.begin(),y1.end(),[&](){return dist_d(mt);});
		generate(z1.begin(),z1.end(),[&](){return dist_d(mt);});
		
		generate(x2.begin(),x2.end(),[&](){return dist_d(mt);});
		generate(y2.begin(),y2.end(),[&](){return dist_d(mt);});
		generate(z2.begin(),z2.end(),[&](){return dist_d(mt);});
		
	
		vector<int> indices1(N1);
		uint q1=0;
		generate(indices1.begin(),indices1.end(),[&q1](){return q1++;});
		std::shuffle(indices1.begin(),indices1.end(), mt);
		
		vector<int> indices2(N2);
		uint q2=0;
		generate(indices2.begin(),indices2.end(),[&q2](){return q2++;});
		std::shuffle(indices2.begin(),indices2.end(), mt);
		
		// We are done with generating fake coordinates in two distinct compartments
		//The next few lines below do the distance calculations
	
		Coords c1(x1,y1,z1,indices1);
		Coords c2(x2,y2,z2,indices2);
		uint N = std::max({N1,N2});
		
		dOut<1U> out(N, {5.0f, 12.0f, 6.0f});
		auto tic = std::chrono::steady_clock::now();
		find_distances(out,c1,c2,tag_simd<cuda,float>());
		auto tac = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tac-tic).count();
		std::cout << "Time(CUDA-EXAMPLE): " << ms << " ms" << std::endl << std::endl;
		report_contact_stats(out);
	}
#endif
	
	void examples_dist_module()
	{	
		cout << "dist_example_between_compartments():" << endl;
		dist_example_between_compartments();

		cout << "dist_example_simplest():" << endl;
		dist_example_simplest();

		cout << "dist_example_nonparallel():" << endl;
		dist_example_nonparallel();

		// cout << "\ndist_example_parallel():" << endl;
		// dist_example_parallel();

#ifdef __CUDACC__
		cout << "\ndist_example_between_compartments_cuda():" << endl;
		dist_example_between_compartments_cuda();
#endif

		cout << "\ndist_example_avx_veccoord():" << endl;
		dist_example_avx_veccoord();
		
	}
} // end-of-namespace dist

#endif