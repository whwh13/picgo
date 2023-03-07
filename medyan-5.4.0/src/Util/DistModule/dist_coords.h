/*
AUTHOR: G.A. Papoian, Date: Nov 22, 2018

A Coords structure that holds coordinates, bead indices and a number of auxiliary variables.
*/

#ifndef DIST_COORDS
#define DIST_COORDS

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <initializer_list>

//#include "dist_moduleV2/dist_common.h"
#include "dist_common.h"
#include "Rand.h"

namespace medyan::dist {

	// This structure accepts coordinates and indicies (or could mock generate them for testing purposes)
	struct Coords {
		// --- Variables ------------------	
		// 
		// public

		std::vector<int> indices;
		std::vector<float> x, y, z;
		std::vector<int> filinfo;
	
		// --- Methods ------------------
	
		template <typename T, typename U>
		void init_coords(const std::vector<T> &xx, const std::vector<T> &yy, const
		std::vector<T> &zz, const std::vector<U> &indx, const std::vector<U> &finfo)
		{
			uint N = xx.size();
			this->resize(N);

			copy(xx.begin(),xx.end(),x.begin());
			copy(yy.begin(),yy.end(),y.begin());
			copy(zz.begin(),zz.end(),z.begin());
			copy(indx.begin(),indx.end(),indices.begin());
			copy(finfo.begin(),finfo.end(),filinfo.begin());
		}
		
		void _init_coords_mock(uint N)
		{
			this->resize(N);	
			create_mock_values();
			uint i=0;
			std::generate(indices.begin(),indices.end(),[&i](){return i++;});
		}
					
		Coords() = default;

		Coords(uint N)
		{
			_init_coords_mock(N);
		}

		template <typename T, typename U>
		Coords(const std::vector<T> &xx, const std::vector<T> &yy, const std::vector<T>
		        &zz, const std::vector<U> &indx, const std::vector<U> &finfo)
		{
			init_coords(xx,yy,zz,indx,finfo);
		}
	
		void resize(uint N){
			x.resize(N);
			y.resize(N);
			z.resize(N);
			indices.resize(N);
			filinfo.resize(N);
		}
			
		uint size() const {return x.size();}
	
		void create_mock_values() {	
		    std::uniform_real_distribution<float> dist_d(1.0, 10.0);
			std::uniform_int_distribution<int> fil_f(0, 100);
		
			uint N = x.size();
	
			for(uint i=0; i<N; ++i){
				x[i] = dist_d(Rand::eng);
				y[i] = dist_d(Rand::eng);
				z[i] = dist_d(Rand::eng);
				filinfo[i] = fil_f(Rand::eng);
			}
		}	
	};

} // end-of-namespace dist

#endif // DIST_COORDS
	
