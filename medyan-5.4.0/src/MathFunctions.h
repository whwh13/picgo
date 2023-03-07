//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_MathFunctions_h
#define MEDYAN_MathFunctions_h

#include <cmath>
#include <vector>
#include <bitset>
#include "common.h"
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "Util/Math/Vec.hpp"

/// @namespace mathfunc is used for the mathematics module for the entire codebase
/// mathfunc includes functions to calculate distances, products, and midpoints

namespace medyan {
namespace mathfunc {

/// Vector and array converter. Need to ensure the vector has size of _Size
// No need for move semantics because normally we use this for copying integers or doubles
template< size_t dim, typename Float >
inline auto vector2Vec(const vector<Float>& v) {
    // Assert v.size() == Size
    medyan::Vec<dim, Float> res;
    for(size_t idx = 0; idx < dim; ++idx){
        res[idx] = v[idx];
    }
    return res;
}
template< typename VecType, size_t = VecType::vec_size >
inline auto vec2Vector(const VecType& a) {
    return vector<typename VecType::float_type>(a.begin(), a.end());
}

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

#else
    static __inline__ __device__ floatingpoint atomicAdd(floatingpoint *address, floatingpoint val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    if (val==0.0)
      return __longlong_as_floatingpoint(old);
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __floatingpoint_as_longlong(val +__longlong_as_floatingpoint(assumed)));
    } while (assumed != old);
    return __longlong_as_floatingpoint(old);
  }
    static __inline__ __device__ floatingpoint atomicExch(floatingpoint *address, floatingpoint val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    if (val==0.0)
      return __longlong_as_floatingpoint(old);
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __floatingpoint_as_longlong(__longlong_as_floatingpoint(0.0)));
    } while (assumed != old);
    return __longlong_as_floatingpoint(old);
  }
#endif
#ifdef CUDAACCLareParallel
     __host__ inline int nextPowerOf2( int n)
    {
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;
        return n;
    }
    __host__ inline vector<int> getaddred2bnt(int nint){
        vector<int> bnt;
        int THREADSPERBLOCK = 0;
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        THREADSPERBLOCK = prop.maxThreadsPerBlock;
        auto maxthreads = 8 * THREADSPERBLOCK;
        int M = max(nextPowerOf2(nint),64*4);
        int blocks, threads;
        if(M > THREADSPERBLOCK){
            if(M > maxthreads) {
                blocks = 8;
                blocks = M /(4 * THREADSPERBLOCK) +1;
                threads = THREADSPERBLOCK;
                std::cout<<"MaxF Number of elements is greater than number of "
                        "threads"<<endl;
                cout<<"Choosing "<<blocks<<" blocks and "<<THREADSPERBLOCK<<" threads per"
                        " block"<<endl;
            }
            else if(M > THREADSPERBLOCK){
                blocks = M /(4 * THREADSPERBLOCK) +1;
                threads = THREADSPERBLOCK;
            }
        }
        else
        { blocks = 1; threads = M/4;}
        //0 M
        //1 THREADSPERBLOCK
        //2 blocks
        //3 threads
        bnt.clear();
        bnt.push_back(M);
        bnt.push_back(THREADSPERBLOCK);
        bnt.push_back(blocks);
        bnt.push_back(threads);
        return bnt;
    }
    __global__ void resetintvariableCUDA(int *variable);
    __global__ void resetfloatingpointvariableCUDA(floatingpoint *variable);
//     __global__ void addvector(floatingpoint *U, int *params, floatingpoint *U_sum);
//    __global__ void addvector(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot, int* culpritID, char* culpritFF,
//                              char* culpritinteraction, char *FF, char *interaction);
    __global__ void addvectorred(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot);
    __global__ void addvectorred2(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot);
    __global__ void addvectorred3(floatingpoint *U, int *params, floatingpoint *U_sum);
    #endif
    /// Normalize a vector
    inline void normalize(vector<floatingpoint> &v) {

        floatingpoint norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }

    /// Return normalized vector not in place
    inline vector<floatingpoint> normalizeVector(const vector<floatingpoint> &v) {

        floatingpoint x2 = v[0] * v[0];
        floatingpoint y2 = v[1] * v[1];
        floatingpoint z2 = v[2] * v[2];

        floatingpoint norm = sqrt( x2+y2+z2);

        vector<floatingpoint> v1;

        v1.push_back(v[0] / norm);
        v1.push_back(v[1] / norm);
        v1.push_back(v[2] / norm);

        return v1;
    }

    /// Return normalized vector
    /// ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
	#endif

    inline uint32_t nextPowerOf2( uint32_t n)
    {
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;
        return n;
    }

    inline void normalizeVector(floatingpoint *v) {

        floatingpoint x2 = v[0] * v[0];
        floatingpoint y2 = v[1] * v[1];
        floatingpoint z2 = v[2] * v[2];

        floatingpoint norm = sqrt( x2+y2+z2);
        *(v) = *(v) / norm;
        *(v + 1) = *(v + 1) / norm;
        *(v + 2) = *(v + 2) / norm;
    }

    /// Get the magnitude of a vector
    inline floatingpoint magnitude(const vector<floatingpoint> &v) {

        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    inline floatingpoint sqmagnitude(const vector<floatingpoint> &v) {

        return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    inline floatingpoint sqmagnitude(const floatingpoint* v) {

        floatingpoint x2 = v[0] * v[0];
        floatingpoint y2 = v[1] * v[1];
        floatingpoint z2 = v[2] * v[2];

        return (x2+y2+z2);
    }

    inline doubleprecision sqmagnitude(const floatingpoint* c1, const floatingpoint* c2) {

        doubleprecision x = c1[0] - c2[0];
        doubleprecision y = c1[1] - c2[1];
        doubleprecision z = c1[2] - c2[2];

        return (x*x+y*y+z*z);
    }
    ///ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint magnitude(floatingpoint const *v) {

        return sqrt((*(v)) * (*(v)) + (*(v + 1)) * (*(v + 1)) + (*(v + 2)) * (*(v + 2)));
    }

    #ifdef CUDAACCL
    __host__ __device__
	#endif
    inline floatingpoint getdistancefromplane(const floatingpoint *coord, const floatingpoint * plane,int id){
        return (plane[0] * coord[id] + plane[1] * coord[id + 1] + plane[2] * coord[id + 2] + plane[3]) /
               sqrt(pow(plane[0], 2) + pow(plane[1], 2) + pow(plane[2], 2));
    }

#ifdef CUDAACCL
    __host__ __device__
#endif
    inline floatingpoint getdistancefromplane(const floatingpoint *coord, const floatingpoint * plane){
        int id = 0;
        return (plane[0] * coord[id] + plane[1] * coord[id + 1] + plane[2] * coord[id + 2] + plane[3]) /
               sqrt(pow(plane[0], 2) + pow(plane[1], 2) + pow(plane[2], 2));
    }

    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint getstretcheddistancefromplane(floatingpoint *coord,
    		floatingpoint *force, floatingpoint * plane, floatingpoint z, int id){
        floatingpoint movedPoint[3] = {coord[id] + z*force[id],
                                     coord[id + 1] + z*force[id + 1],
                                     coord[id + 2] + z*force[id + 2]};
        return (plane[0] * movedPoint[0] + plane[1] * movedPoint[1] + plane[2] * movedPoint[2] + plane[3]) /
               sqrt(pow(plane[0], 2) + pow(plane[1], 2) + pow(plane[2], 2));
    }
    //@{
    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3)
    inline floatingpoint twoPointDistance(const vector<floatingpoint> &v1, const vector<floatingpoint> &v2) {

        return sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                    (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    inline floatingpoint twoPointDistancesquared(const vector<floatingpoint> &v1, const
    vector<floatingpoint> &v2) {

	    floatingpoint d0 = v2[0] - v1[0];
	    floatingpoint d1 = v2[1] - v1[1];
	    floatingpoint d2 = v2[2] - v1[2];

        return d0*d0 + d1*d1 + d2*d2;
    }

    inline floatingpoint twoPointDistance(floatingpoint const *v1, const vector<floatingpoint> &v2) {

        return sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                    (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    inline floatingpoint twoPointDistance(const vector<floatingpoint> &v1, floatingpoint const *v2) {

        return sqrt((*(v2) - v1[0]) * (*(v2) - v1[0]) +
                    (*(v2 + 1) - v1[1]) * (*(v2 + 1) - v1[1]) +
                    (*(v2 + 2) - v1[2]) * (*(v2 + 2) - v1[2]));
    }
    //@}

    /// Compute distance between two points with coordinates: (x1,y1,z1) and (x2,y2,z3)
    /// ARRAY VERSION
    #ifdef CUDAACCL
     __host__ __device__
     #endif
    inline floatingpoint twoPointDistance(floatingpoint const *v1, floatingpoint const *v2) {

        return sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    inline floatingpoint twoPointDistancesquared(floatingpoint const *v1, floatingpoint const *v2) {
        floatingpoint dx = (v2[0] - v1[0]);
        floatingpoint dy = (v2[1] - v1[1]);
        floatingpoint dz = (v2[2] - v1[2]);
        return(dx*dx+dy*dy+dz*dz);

//        return ((v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]) +
//                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }

    inline floatingpoint twoPointDistancesquared(floatingpoint const *v1, floatingpoint const *v2, int id1,
                                          int id2) {
        floatingpoint dx = (v2[id1] - v1[id2]);
        floatingpoint dy = (v2[id1 + 1] - v1[id2 + 1]);
        floatingpoint dz = (v2[id1 + 2] - v1[id2 + 2]);
        return(dx*dx+dy*dy+dz*dz);

//        return ((v2[0] - v1[0]) * (v2[0] - v1[0]) + (v2[1] - v1[1]) * (v2[1] - v1[1]) +
//                    (v2[2] - v1[2]) * (v2[2] - v1[2]));
    }
    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint twoPointDistance(floatingpoint const *v1, floatingpoint const *v2, int const id) {

        return sqrt((v2[id] - v1[id]) * (v2[id] - v1[id]) + (v2[id + 1] - v1[id + 1]) * (v2[id + 1] - v1[id + 1])
                    + (v2[id + 2] - v1[id + 2]) * (v2[id + 2] - v1[id + 2]));
    }
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint twoPointDistancemixedID(floatingpoint const *v1, floatingpoint const *v2, int const id1, int const id2) {

        return sqrt((v2[id2] - v1[id1]) * (v2[id2] - v1[id1]) + (v2[id2 + 1] - v1[id1 + 1]) * (v2[id2 + 1] - v1[id1 +
                1]) + (v2[id2 + 2] - v1[id1 + 2]) * (v2[id2 + 2] - v1[id1 + 2]));
    }
    /// Compute distance between two points with coordinates
    /// (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z)
    inline floatingpoint twoPointDistanceStretched(const vector<floatingpoint> &v1,
                                            const vector<floatingpoint> &p1,
                                            const vector<floatingpoint> &v2,
                                            const vector<floatingpoint> &p2, floatingpoint d) {

        return sqrt(((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) +
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) +
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])));
    }

    /// Compute distance between two points with coordinates
    /// (x1 -d*p1x,y1-d*p1y,z1-d*p1z) and (x2-d*p2x,y2-d*p2y,z2-d*p2z)
    /// CUDA & ARRAY VERSION
#ifdef CUDAACCL
    __host__ __device__
#endif
    inline floatingpoint twoPointDistanceStretched(floatingpoint const *v1,
                                                   floatingpoint const *p1,
                                                    floatingpoint const *v2,
                                                   floatingpoint const *p2,
                                                   floatingpoint d) {
        return sqrt(((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) +
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) +
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])));
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint twoPointDistanceStretched(floatingpoint const *v1,
                                                   floatingpoint const *p1,
                                                    floatingpoint const *v2,
                                                   floatingpoint const *p2, floatingpoint d, int const id) {
        return sqrt(((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                    ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) +
                    ((v2[id + 1] + d * p2[id + 1]) - (v1[id + 1] + d * p1[id + 1])) *
                    ((v2[id + 1] + d * p2[id + 1]) - (v1[id + 1] + d * p1[id + 1])) +
                    ((v2[id + 2] + d * p2[id + 2]) - (v1[id + 2] + d * p1[id + 2])) *
                    ((v2[id + 2] + d * p2[id + 2]) - (v1[id + 2] + d * p1[id + 2])));
    }
    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint twoPointDistanceStretchedmixedID(floatingpoint const *v1,
                                                          floatingpoint const *p1,
                                            floatingpoint const *v2,
                                                          floatingpoint const *p2, floatingpoint d, int const id1, const int id2) {
        return sqrt(((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                    ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) +
                    ((v2[id2 + 1] + d * p2[id2 + 1]) - (v1[id1 + 1] + d * p1[id1 + 1])) *
                    ((v2[id2 + 1] + d * p2[id2 + 1]) - (v1[id1 + 1] + d * p1[id1 + 1])) +
                    ((v2[id2 + 2] + d * p2[id2 + 2]) - (v1[id1 + 2] + d * p1[id1 + 2])) *
                    ((v2[id2 + 2] + d * p2[id2 + 2]) - (v1[id1 + 2] + d * p1[id1 + 2])));
    }
    //@{
    /// Calculates a normal to a line starting at (x1,y1,z1) and ending at (x2,y2,z2)
    inline vector<floatingpoint> twoPointDirection(const vector<floatingpoint> &v1,
                                            const vector<floatingpoint> &v2) {
        vector<floatingpoint> tau(3, 0);
        floatingpoint invD = 1 / twoPointDistance(v1, v2);
        tau[0] = invD * (v2[0] - v1[0]);
        tau[1] = invD * (v2[1] - v1[1]);
        tau[2] = invD * (v2[2] - v1[2]);
        return tau;
    }

    inline vector<floatingpoint> twoPointDirection(floatingpoint const *v1,
                                            const vector<floatingpoint> &v2) {
        vector<floatingpoint> tau(3, 0);
        floatingpoint invD = 1 / twoPointDistance(v1, v2);
        tau[0] = invD * (v2[0] - v1[0]);
        tau[1] = invD * (v2[1] - v1[1]);
        tau[2] = invD * (v2[2] - v1[2]);
        return tau;
    }

    inline void twoPointDirection(floatingpoint *tau,
                                  floatingpoint const *v1,
                                  floatingpoint const *v2) {

        floatingpoint invD = 1 / twoPointDistance(v1, v2);
        tau[0] = invD * (*(v2) - *(v1));
        tau[1] = invD * (*(v2 + 1) - *(v1 + 1));
        tau[2] = invD * (*(v2 + 2) - *(v1 + 2));
    }
    //@}

    /// Scalar product of two vectors v1(x,y,z) and v2(x,y,z)
    inline floatingpoint dotProduct(const vector<floatingpoint> &v1, const vector<floatingpoint> &v2) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    /// Scalar product of two vectors v1(x,y,z) and v2(x,y,z)
    /// ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
	#endif
	template <class dataType = floatingpoint>
    inline dataType dotProduct(dataType const *v1, dataType const *v2) {
        return (*(v1)) * (*(v2)) + (*(v1 + 1)) * (*(v2 + 1)) + (*(v1 + 2)) * (*(v2 + 2));

//        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }

    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    template<class dataType = floatingpoint>
    inline dataType scalarProduct(const vector<dataType> &v1, const vector<dataType>
            &v2, const vector<dataType> &v3, const vector<dataType> &v4) {

        return ((v2[0] - v1[0]) * (v4[0] - v3[0]) +
                (v2[1] - v1[1]) * (v4[1] - v3[1]) +
                (v2[2] - v1[2]) * (v4[2] - v3[2]));
    }

    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3)
    /// ARRAY VERSION & CUDA version
#ifdef CUDAACCL
    __host__ __device__
#endif
    template<class dataType = floatingpoint>
    inline dataType scalarProduct(dataType const *v1, dataType const *v2,
                                  dataType const *v3, dataType const *v4) {
//        return((*(v2)-*(v1))*(*(v4)-*(v3))
//               +(*(v2+1)-*(v1+1))*(*(v4+1)-*(v3+1))
//               +(*(v2+2)-*(v1+2))*(*(v4+2)-*(v3+2)));
        return ((v2[0] - v1[0]) * (v4[0] - v3[0]) +
                (v2[1] - v1[1]) * (v4[1] - v3[1]) +
                (v2[2] - v1[2]) * (v4[2] - v3[2]));
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline floatingpoint scalarProduct(floatingpoint const *v1, floatingpoint const *v2,
                                floatingpoint const *v3, floatingpoint const *v4, int const id) {

        return ((v2[id] - v1[id]) * (v4[id] - v3[id]) +
                (v2[id + 1] - v1[id + 1]) * (v4[id + 1] - v3[id + 1]) +
                (v2[id + 2] - v1[id + 2]) * (v4[id + 2] - v3[id + 2]));
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
	#endif

	inline doubleprecision scalarProduct(doubleprecision const *v1, floatingpoint const
	*v2){
        return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    }


    inline floatingpoint scalarProductmixedID(floatingpoint const *v1, floatingpoint const *v2,
                                floatingpoint const *v3, floatingpoint const *v4, int
                                const id1, int const id2, int const id3, int
                                       const id4) {
        return ((v2[id2] - v1[id1]) * (v4[id4] - v3[id3]) +
                (v2[id2 + 1] - v1[id1 + 1]) * (v4[id4 + 1] - v3[id3 + 1]) +
                (v2[id2 + 2] - v1[id1 + 2]) * (v4[id4 + 2] - v3[id3 + 2]));
    }


    /// Scalar product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3) but with x+d*p coordinates
    ///
    //CUDA version & ARRAY VERSION
#ifdef CUDAACCL
    __host__ __device__
#endif
    inline floatingpoint scalarProductStretched(floatingpoint const *v1,
                                                floatingpoint const *p1,
                                         floatingpoint const *v2,
                                                floatingpoint const *p2,
                                         floatingpoint const *v3,
                                                floatingpoint const *p3,
                                         floatingpoint const *v4,
                                                floatingpoint const *p4,
                                         floatingpoint d) {
        floatingpoint xx = ((v2[0] + d * p2[0]) - (v1[0] + d * p1[0])) *
                    ((v4[0] + d * p4[0]) - (v3[0] + d * p3[0]));
        floatingpoint yy = ((v2[1] + d * p2[1]) - (v1[1] + d * p1[1])) *
                    ((v4[1] + d * p4[1]) - (v3[1] + d * p3[1]));
        floatingpoint zz = ((v2[2] + d * p2[2]) - (v1[2] + d * p1[2])) *
                    ((v4[2] + d * p4[2]) - (v3[2] + d * p3[2]));
        return xx + yy + zz;
//        floatingpoint xx = ((*(v2) + d*(*(p2))) - (*(v1) + d*(*(p1)))) *
//        ((*(v4)) + d*(*(p4))) - ((*(v3)) + d*(*(p3)));
//        floatingpoint yy = ((*(v2+1) + d*(*(p2+1))) - (*(v1+1) + d*(*(p1+1))))*
//        ((*(v4+1) + d*(*(p4+1))) - (*(v3+1) + d*(*(p3+1))));
//        floatingpoint zz = ((*(v2+2) + d*(*(p2+2))) - (*(v1+2) + d*(*(p1+2))))*
//        ((*(v4+2) + d*(*(p4+2))) - (*(v3+2) + d*(*(p3+2))));
//        return xx + yy + zz;

    }

    //CUDA version
    #ifdef CUDAACCL
     __host__ __device__
     #endif
    inline floatingpoint scalarProductStretched(floatingpoint const *v1,
                                         floatingpoint const *p1,
                                                floatingpoint const *v2,
                                         floatingpoint const *p2,
                                                floatingpoint const *v3,
                                         floatingpoint const *p3,
                                                floatingpoint const *v4,
                                         floatingpoint const *p4,
                                         floatingpoint d,int const id ) {
        floatingpoint xx = ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                    ((v4[id] + d * p4[id]) - (v3[id] + d * p3[id]));
        floatingpoint yy = ((v2[id + 1] + d * p2[id + 1]) - (v1[id + 1] + d * p1[id + 1])) *
                    ((v4[id + 1] + d * p4[id + 1]) - (v3[id + 1] + d * p3[id + 1]));
        floatingpoint zz = ((v2[id + 2] + d * p2[id + 2]) - (v1[id + 2] + d * p1[id + 2])) *
                    ((v4[id + 2] + d * p4[id + 2]) - (v3[id + 2] + d * p3[id + 2]));
        return xx + yy + zz;
    }

    #ifdef CUDAACCL
    __host__ __device__
     #endif
    inline floatingpoint scalarProductStretchedmixedIDv2(floatingpoint const *v1,

                                         floatingpoint const *v2,

                                         floatingpoint const *v3,

                                         floatingpoint const *v4,

                                         floatingpoint d, int const id, floatingpoint const *p, int const id1, int const id2, int
                                                  const id3, int const id4) {
        floatingpoint xx = ((v2[id] + d * p[id2]) - (v1[id] + d * p[id1])) *
                    ((v4[id] + d * p[id4]) - (v3[id] + d * p[id3]));
        floatingpoint yy = ((v2[id + 1] + d * p[id2 + 1]) - (v1[id + 1] + d * p[id1 + 1])) *
                    ((v4[id + 1] + d * p[id4 + 1]) - (v3[id + 1] + d * p[id3 + 1]));
        floatingpoint zz = ((v2[id + 2] + d * p[id2 + 2]) - (v1[id + 2] + d * p[id1 + 2])) *
                    ((v4[id + 2] + d * p[id4 + 2]) - (v3[id + 2] + d * p[id3 + 2]));
        return xx + yy + zz;
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
     #endif
    inline floatingpoint scalarProductStretchedmixedID(floatingpoint const *v1,
                                                       floatingpoint const *p1,
                                         floatingpoint const *v2,
                                                       floatingpoint const *p2,
                                         floatingpoint const *v3,
                                                       floatingpoint const *p3,
                                         floatingpoint const *v4,
                                                       floatingpoint const *p4,
                                         floatingpoint d,
                                         int const id1,
                                         int const id2,
                                         int const id3,
                                         int const id4) {
        floatingpoint xx = ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                    ((v4[id4] + d * p4[id4]) - (v3[id3] + d * p3[id3]));
        floatingpoint yy = ((v2[id2 + 1] + d * p2[id2 + 1]) - (v1[id1 + 1] + d * p1[id1 + 1])) *
                    ((v4[id4 + 1] + d * p4[id4 + 1]) - (v3[id3 + 1] + d * p3[id3 + 1]));
        floatingpoint zz = ((v2[id2 + 2] + d * p2[id2 + 2]) - (v1[id1 + 2] + d * p1[id1 + 2])) *
                    ((v4[id4 + 2] + d * p4[id4 + 2]) - (v3[id3 + 2] + d * p3[id3 + 2]));
        return xx + yy + zz;
    }

    inline floatingpoint scalarprojection(vector<floatingpoint> a, vector<floatingpoint> b){
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    }

    inline floatingpoint scalarprojection(floatingpoint* a, floatingpoint* b){
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    }
    inline floatingpoint maxdistbetweencylinders(const vector<floatingpoint> &v1,
                                   const vector<floatingpoint> &v2,
                                   const vector<floatingpoint> &v3,
                                   const vector<floatingpoint> &v4) {
        floatingpoint maxdist = 0.0;
        floatingpoint d1 = twoPointDistancesquared(v1, v3);maxdist = max(maxdist,d1);
        floatingpoint d2 = twoPointDistancesquared(v1, v4);maxdist = max(maxdist,d2);
        floatingpoint d3 = twoPointDistancesquared(v2, v3);maxdist = max(maxdist,d3);
        floatingpoint d4 = twoPointDistancesquared(v2, v4);maxdist = max(maxdist,d4);
        return maxdist;
    }

    inline floatingpoint maxdistbetweencylinders(const floatingpoint* v1,
                                          const floatingpoint* v2,
                                          const floatingpoint* v3,
                                          const floatingpoint* v4) {
        floatingpoint maxdist = 0.0;
        floatingpoint d1 = twoPointDistancesquared(v1, v3);maxdist = max(maxdist,d1);
        floatingpoint d2 = twoPointDistancesquared(v1, v4);maxdist = max(maxdist,d2);
        floatingpoint d3 = twoPointDistancesquared(v2, v3);maxdist = max(maxdist,d3);
        floatingpoint d4 = twoPointDistancesquared(v2, v4);maxdist = max(maxdist,d4);
        return maxdist;
    }

    inline floatingpoint maxdistsqbetweenbindingsites(const floatingpoint* v1,
                                          const floatingpoint* v2) {
        floatingpoint maxdist = 0.0;
        floatingpoint d1 = twoPointDistancesquared(v1, v2,0,0);maxdist = max(maxdist,d1);
        floatingpoint d2 = twoPointDistancesquared(v1, v2,3,3);maxdist = max(maxdist,d2);
        floatingpoint d3 = twoPointDistancesquared(v2, v2,6,6);maxdist = max(maxdist,d3);
        floatingpoint d4 = twoPointDistancesquared(v2, v2,9,9);maxdist = max(maxdist,d4);
        return maxdist;
    }

    /// Scalar product of two vectors with coordinates: v1[z,y,z] + d*p1[x,y,z] and
    /// v2[x,y,z] + d*p2[x,y,z]
    inline floatingpoint dotProductStretched(const vector<floatingpoint> &v1,
                                      const vector<floatingpoint> &p1,
                                      const vector<floatingpoint> &v2,
                                      const vector<floatingpoint> &p2,
                                      floatingpoint d) {

        floatingpoint xx = (v1[0] + d * p1[0]) * (v2[0] + d * p2[0]);
        floatingpoint yy = (v1[1] + d * p1[1]) * (v2[1] + d * p2[1]);
        floatingpoint zz = (v1[2] + d * p1[2]) * (v2[2] + d * p2[2]);
        return xx + yy + zz;

    }


    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3). Returns a 3d vector.
    inline vector<floatingpoint> vectorProduct(const vector<floatingpoint> &v1,
                                        const vector<floatingpoint> &v2,
                                        const vector<floatingpoint> &v3,
                                        const vector<floatingpoint> &v4) {
        vector<floatingpoint> v;

        floatingpoint vx = (v2[1] - v1[1]) * (v4[2] - v3[2]) - (v2[2] - v1[2]) * (v4[1] - v3[1]);
        floatingpoint vy = (v2[2] - v1[2]) * (v4[0] - v3[0]) - (v2[0] - v1[0]) * (v4[2] - v3[2]);
        floatingpoint vz = (v2[0] - v1[0]) * (v4[1] - v3[1]) - (v2[1] - v1[1]) * (v4[0] - v3[0]);

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;
    };

    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3). Returns a 3d vector.
    /// ARRAY VERSION
    template <class dataType = floatingpoint>
    inline void vectorProduct(dataType *v,
                              dataType const *v1,
                              dataType const *v2,
                              dataType const *v3,
                              dataType const *v4) {
	    dataType v211 = (v2[1] - v1[1]);
	    dataType v430 = (v4[0] - v3[0]);
	    dataType v212 = (v2[2] - v1[2]);
	    dataType v432 = (v4[2] - v3[2]);
	    dataType v431 = (v4[1] - v3[1]);
	    dataType v210 = (v2[0] - v1[0]);
	    dataType vx =  v211 * v432 - v212 * v431;
	    dataType vy = v212 * v430 - v210 * v432;
	    dataType vz = v210 * v431 - v211 * v430;

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };

    /// CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProduct(floatingpoint *v,
                              floatingpoint const *v1,
                              floatingpoint const *v2,
                              floatingpoint const *v3,
                              floatingpoint const *v4, int const id) {
        floatingpoint vx =
                (v2[id+1] - v1[id+1]) * (v4[id+2] - v3[id+2]) - (v2[id+2] - v1[id+2]) * (v4[id+1] - v3[id+1]);
        floatingpoint vy = (v2[id+2] - v1[id+2]) * (v4[id] - v3[id]) - (v2[id] - v1[id]) * (v4[id+2] - v3[id+2]);
        floatingpoint vz = (v2[id] - v1[id]) * (v4[id+1] - v3[id+1]) - (v2[id+1] - v1[id+1]) * (v4[id] - v3[id]);

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    /// CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProductmixedID(floatingpoint *v,
                              floatingpoint const *v1,
                              floatingpoint const *v2,
                              floatingpoint const *v3,
                              floatingpoint const *v4, int const id1, int const id2,int const id3, int const id4) {
        floatingpoint vx =
                (v2[id2+1] - v1[id1+1]) * (v4[id4+2] - v3[id3+2]) - (v2[id2+2] - v1[id1+2]) * (v4[id4+1] - v3[id3+1]);
        floatingpoint vy = (v2[id2+2] - v1[id1+2]) * (v4[id4] - v3[id3]) - (v2[id2] - v1[id1]) * (v4[id4+2] - v3[id3+2]);
        floatingpoint vz = (v2[id2] - v1[id1]) * (v4[id4+1] - v3[id3+1]) - (v2[id2+1] - v1[id1+1]) * (v4[id4] - v3[id3]);

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };

    /// Vector product of two vectors with coordinates: (x2-x1,y2-y1,z2-z1) and
    /// (x4-x3,y4-y3,z4-z3), but with v -> v+d*p. Returns a 3d vector.
    /// ARRAY VERSION
    inline vector<floatingpoint> vectorProductStretched (const vector<floatingpoint>& v1,
                                                  const vector<floatingpoint>& p1,
                                                  const vector<floatingpoint>& v2,
                                                  const vector<floatingpoint>& p2,
                                                  const vector<floatingpoint>& v3,
                                                  const vector<floatingpoint>& p3,
                                                  const vector<floatingpoint>& v4,
                                                  const vector<floatingpoint>& p4,
                                                  floatingpoint d){
        vector<floatingpoint> v;

        floatingpoint vx =
                ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]))
                - ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]));

        floatingpoint vy =
                ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]))
                - ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]));

        floatingpoint vz =
                ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]))
                - ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]));

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;


    };
    //ARRAY VERSION
    //BranchingDihedralCosineV1
    inline void vectorProductStretched(floatingpoint *v,
                                       floatingpoint const *v1,
                                       floatingpoint const *p1,
                                       floatingpoint const *v2,
                                       floatingpoint const *p2,
                                       floatingpoint const *v3,
                                       floatingpoint const *p3,
                                       floatingpoint const *v4,
                                       floatingpoint const *p4,
                                       floatingpoint d) {
        floatingpoint vx =
                ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]))
                - ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]));

        floatingpoint vy =
                ((v2[2]+d*p2[2])-(v1[2]+d*p1[2]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]))
                - ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[2]+d*p4[2])-(v3[2]+d*p3[2]));

        floatingpoint vz =
                ((v2[0]+d*p2[0])-(v1[0]+d*p1[0]))*((v4[1]+d*p4[1])-(v3[1]+d*p3[1]))
                - ((v2[1]+d*p2[1])-(v1[1]+d*p1[1]))*((v4[0]+d*p4[0])-(v3[0]+d*p3[0]));
        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    }

    //
    ///CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProductStretched(floatingpoint *v,
                                       floatingpoint const *v1,
                                       floatingpoint const *p1,
                                       floatingpoint const *v2,
                                       floatingpoint const *p2,
                                       floatingpoint const *v3,
                                       floatingpoint const *p3,
                                       floatingpoint const *v4,
                                       floatingpoint const *p4,
                                       floatingpoint d, int const id) {
        floatingpoint vx =
                ((v2[id+1] + d * p2[id+1]) - (v1[id+1] + d * p1[id+1])) *
                ((v4[id+2] + d * p4[id+2]) - (v3[id+2] + d * (p3[id+2])))
                - ((v2[id+2] + d * (p2[id+2])) - (v1[id+2] + d * (p1[id+2]))) *
                  ((v4[id+1] + d * (p4[id+1])) - (v3[id+1] + d * (p3[id+1])));

        floatingpoint vy =
                ((v2[id+2] + d * p2[id+2]) - (v1[id+2] + d * p1[id+2])) *
                (v4[id] + d * p4[id] - (v3[id] + d * p3[id]))
                - ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                  ((v4[id+2] + d * (p4[id+2])) - (v3[id+2] + d * p3[id+2]));

        floatingpoint vz =
                ((v2[id] + d * p2[id]) - (v1[id] + d * p1[id])) *
                ((v4[id+1] + d * p4[id+1]) - (v3[id+1] + d * p3[id+1]))
                - ((v2[id+1] + d * p2[id+1]) - (v1[id+1] + d * p1[id+1])) *
                  ((v4[id] + d * p4[id]) - (v3[id] + d * p3[id]));

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    ///CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void vectorProductStretchedmixedID(floatingpoint *v,
                                       floatingpoint const *v1,
                                              floatingpoint const *p1,
                                       floatingpoint const *v2,
                                              floatingpoint const *p2,
                                       floatingpoint const *v3,
                                              floatingpoint const *p3,
                                       floatingpoint const *v4,
                                              floatingpoint const *p4,
                                       floatingpoint d, int const id1, int const id2, int const id3, int const id4) {
        floatingpoint vx =
                ((v2[id2+1] + d * p2[id2+1]) - (v1[id1+1] + d * p1[id1+1])) *
                ((v4[id4+2] + d * p4[id4+2]) - (v3[id3+2] + d * (p3[id3+2])))
                - ((v2[id2+2] + d * (p2[id2+2])) - (v1[id1+2] + d * (p1[id1+2]))) *
                  ((v4[id4+1] + d * (p4[id4+1])) - (v3[id3+1] + d * (p3[id3+1])));

        floatingpoint vy =
                ((v2[id2+2] + d * p2[id2+2]) - (v1[id1+2] + d * p1[id1+2])) *
                (v4[id4] + d * p4[id4] - (v3[id3] + d * p3[id3]))
                - ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                  ((v4[id4+2] + d * (p4[id4+2])) - (v3[id3+2] + d * p3[id3+2]));

        floatingpoint vz =
                ((v2[id2] + d * p2[id2]) - (v1[id1] + d * p1[id1])) *
                ((v4[id4+1] + d * p4[id4+1]) - (v3[id3+1] + d * p3[id3+1]))
                - ((v2[id2+1] + d * p2[id2+1]) - (v1[id1+1] + d * p1[id1+1])) *
                  ((v4[id4] + d * p4[id4]) - (v3[id3] + d * p3[id3]));

        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    };
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z]. Returns a 3d vector.

    inline vector<floatingpoint> crossProduct(const vector<floatingpoint> &v1,
                                       const vector<floatingpoint> &v2) {
        vector<floatingpoint> v;

        floatingpoint vx = v1[1] * v2[2] - v1[2] * v2[1];
        floatingpoint vy = v1[2] * v2[0] - v1[0] * v2[2];
        floatingpoint vz = v1[0] * v2[1] - v1[1] * v2[0];

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;
    };
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z].
    /// ARRAY VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif

	template <class dataType>
	inline void crossProduct(dataType *cp,
	                         floatingpoint const *v1,
	                         floatingpoint const *v2) {
		cp[0] = v1[1] * v2[2] - v1[2] * v2[1];
		cp[1] = v1[2] * v2[0] - v1[0] * v2[2];
		cp[2] = v1[0] * v2[1] - v1[1] * v2[0];
	};
    /// Vector product of two vectors v1[x,y,z] and v2[x,y,z]. Returns a 3d vector.
    inline vector<floatingpoint> crossProductStretched(const vector<floatingpoint> &v1,
                                                const vector<floatingpoint> &p1,
                                                const vector<floatingpoint> &v2,
                                                const vector<floatingpoint> &p2,
                                                floatingpoint d) {
        vector<floatingpoint> v;

        floatingpoint vx = (v1[1] + d * p1[1]) * (v2[2] + d * p2[2]) - (v1[2] + d * p1[2]) * (v2[1] + d * p2[1]);
        floatingpoint vy = (v1[2] + d * p1[2]) * (v2[0] + d * p2[0]) - (v1[0] + d * p1[0]) * (v2[2] + d * p2[2]);
        floatingpoint vz = (v1[0] + d * p1[0]) * (v2[1] + d * p2[1]) - (v1[1] + d * p1[1]) * (v2[0] + d * p2[0]);

        v.push_back(vx);
        v.push_back(vy);
        v.push_back(vz);

        return v;
    };

    /// Projection of a new point based on a given direction and starting point
    inline vector<floatingpoint> nextPointProjection(const vector<floatingpoint> &coordinate,
                                              floatingpoint d, const vector<floatingpoint> &tau) {
        vector<floatingpoint> v;
        v.push_back(coordinate[0] + d * tau[0]);
        v.push_back(coordinate[1] + d * tau[1]);
        v.push_back(coordinate[2] + d * tau[2]);
        return v;
    }

    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v1| = alpha.
    inline vector<floatingpoint> midPointCoordinate(const vector<floatingpoint> &v1,
                                             const vector<floatingpoint> &v2, floatingpoint alpha) {
        vector<floatingpoint> v;
        v.push_back(v1[0] * (1.0 - alpha) + alpha * v2[0]);
        v.push_back(v1[1] * (1.0 - alpha) + alpha * v2[1]);
        v.push_back(v1[2] * (1.0 - alpha) + alpha * v2[2]);
        return v;
    }
    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v1| = alpha. ARRAY VERSION
    //CUDA & ARRAY version
#ifdef CUDAACCL
    __host__ __device__
#endif
    template <class dataType = floatingpoint>
    inline void midPointCoordinate(dataType *v, dataType const *v1, dataType const *v2, dataType alpha) {

        dataType beta = 1 - alpha;
        v[0] = (v1[0] * beta + alpha * v2[0]);
        v[1] = (v1[1] * beta + alpha * v2[1]);
        v[2] = (v1[2] * beta + alpha * v2[2]);
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void midPointCoordinate(floatingpoint *v, floatingpoint const *v1, floatingpoint const *v2, floatingpoint alpha, int id) {
        v[0] = v1[id] * (1.0 - alpha) + alpha * v2[id];
        v[1] = v1[id + 1] * (1.0 - alpha) + alpha * v2[id + 1];
        v[2] = v1[id + 2] * (1.0 - alpha) + alpha * v2[id + 2];
    }

    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v| = alpha, but with x-d*p coordinates
    inline vector<floatingpoint> midPointCoordinateStretched(const vector<floatingpoint> &v1,
                                                      const vector<floatingpoint> &p1,
                                                      const vector<floatingpoint> &v2,
                                                      const vector<floatingpoint> &p2,
                                                      floatingpoint alpha, floatingpoint d) {

        vector<floatingpoint> v;
        v.push_back((v1[0] + d * p1[0]) * (1.0 - alpha) + alpha * (v2[0] + d * p2[0]));
        v.push_back((v1[1] + d * p1[1]) * (1.0 - alpha) + alpha * (v2[1] + d * p2[1]));
        v.push_back((v1[2] + d * p1[2]) * (1.0 - alpha) + alpha * (v2[2] + d * p2[2]));
        return v;
    }

    /// Returns coordinates of a point v located on a line between v1 and v2.
    /// |v-v1|/|v2-v| = alpha, but with x-d*p coordinates
    /// ARRAY & CUDA VERSION
#ifdef CUDAACCL
    __host__ __device__
#endif
    inline void midPointCoordinateStretched(floatingpoint *v,
                                            floatingpoint const *v1,
                                            floatingpoint const *p1,
                                            floatingpoint const *v2,
                                            floatingpoint const *p2,
                                            floatingpoint alpha, floatingpoint d) {

        v[0] = (v1[0] + d * p1[0]) * (1.0 - alpha) + alpha * (v2[0] + d * p2[0]);
        v[1] = ((v1[1] + d * p1[1]) * (1.0 - alpha) + alpha * (v2[1] + d * p2[1]));
        v[2] = ((v1[2] + d * p1[2]) * (1.0 - alpha) + alpha * (v2[2] + d * p2[2]));
    }

    //CUDA version
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void midPointCoordinateStretched(floatingpoint *v,
                                            floatingpoint *v1,
                                            floatingpoint *p1,
                                            floatingpoint *v2,
                                            floatingpoint *p2,
                                            floatingpoint alpha, floatingpoint d, int id) {

        v[0] = (v1[id] + d * p1[id]) * (1.0 - alpha) + alpha * (v2[id] + d * p2[id]);
        v[1] = ((v1[id + 1] + d * p1[id + 1]) * (1.0 - alpha) + alpha * (v2[id + 1] + d * p2[id + 1]));
        v[2] = ((v1[id + 2] + d * p1[id + 2]) * (1.0 - alpha) + alpha * (v2[id + 2] + d * p2[id + 2]));
//        printf("%f \n",(v1[id] + d * p1[id]) * (1.0 - alpha) + alpha * (v2[id] + d * p2[id]));
    }

    /// Returns true if two vectors (p1->p2 and p3->p4) are parallel
    inline bool areParallel(const vector<floatingpoint> &p1, const vector<floatingpoint> &p2,
                            const vector<floatingpoint> &p3, const vector<floatingpoint> &p4) {

        auto v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        auto v2 = {p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]};

        return areEqual(magnitude(crossProduct(v1, v2)), 0.0);
    }

    /// ARRAY VERSION
    inline bool areParallel(floatingpoint const *p1, floatingpoint const *p2,
                            floatingpoint const *p3, floatingpoint const *p4) {

        floatingpoint *v1 = new floatingpoint[3];
        floatingpoint *v2 = new floatingpoint[3];
        floatingpoint *cp = new floatingpoint[3];


        v1[0] = p2[0] - p1[0];
        v1[1] = p2[1] - p1[1];
        v1[2] = p2[2] - p1[2];

        v2[0] = p4[0] - p3[0];
        v2[1] = p4[1] - p3[1];
        v2[2] = p4[2] - p3[2];

        crossProduct<floatingpoint>(cp, v1, v2);

        auto retVal = areEqual(magnitude(cp), 0.0);
        delete [] v1;
        delete [] v2;
        delete [] cp;

        return retVal;
    }

    /// CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline bool areParallel(floatingpoint const *p1, floatingpoint const *p2,
                            floatingpoint const *p3, floatingpoint const *p4, int id) {

        floatingpoint v1[3];
        floatingpoint v2[3];
        floatingpoint cp[3];


        v1[0] = p2[id] - p1[id];
        v1[1] = p2[id + 1] - p1[id + 1];
        v1[2] = p2[id + 2] - p1[id + 2];

        v2[0] = p4[id] - p3[id];
        v2[1] = p4[id + 1] - p3[id + 1];
        v2[2] = p4[id + 2] - p3[id + 2];

        crossProduct<floatingpoint>(cp, v1, v2);
//        printf("cp %d %f %f %f\n", id, cp[0], cp[1], cp[2]);
        auto retVal = areEqual(magnitude(cp), 0.0);
//        printf("aE %d %d \n",id, retVal);
//        delete [] v1;
//        delete [] v2;
//        delete [] cp;

        return retVal;
    }

    /// Returns true if two vectors (p1->p2 and p3->p4) are in the same plane
    /// ARRAY VERSION
    inline bool areInPlane(floatingpoint const *p1, floatingpoint const *p2,
                           floatingpoint const *p3, floatingpoint const *p4) {

        floatingpoint *v1 = new floatingpoint[3];
        floatingpoint *v2 = new floatingpoint[3];
        floatingpoint *v3 = new floatingpoint[3];
        floatingpoint *cp = new floatingpoint[3];

        *(v1) = *(p2) - *(p1);
        *(v1 + 1) = *(p2 + 1) - *(p1 + 1);
        *(v1 + 2) = *(p2 + 2) - *(p1 + 2);

        *(v2) = *(p3) - *(p1);
        *(v2 + 1) = *(p3 + 1) - *(p1 + 1);
        *(v2 + 2) = *(p3 + 2) - *(p1 + 2);

        *(v3) = *(p4) - *(p1);
        *(v3 + 1) = *(p4 + 1) - *(p1 + 1);
        *(v3 + 2) = *(p4 + 2) - *(p1 + 2);

        crossProduct<floatingpoint>(cp, v1, v2);
        //Normalize before checking the angle.
	    normalizeVector(cp);
	    normalizeVector(v3);
	    //Check for equality with a lowered threshold
        auto retVal = areEqualLT(dotProduct(v3, cp), (floatingpoint)0.0);
        delete [] v1;
        delete [] v2;
        delete [] v3;
        delete [] cp;

        return retVal;
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline bool areInPlane(floatingpoint const *p1, floatingpoint const *p2,
                           floatingpoint const *p3, floatingpoint const *p4, int const id) {

        floatingpoint v1[3];
        floatingpoint v2[3];
        floatingpoint v3[3];
        floatingpoint cp[3];

        v1[0] = p2[id] - p1[id];
        v1[1] = p2[id + 1] - p1[id + 1];
        v1[2] = p2[id + 2] - p1[id + 2];

        v2[0] = p3[id] - p1[id];
        v2[1] = p3[id + 1] - p1[id + 1];
        v2[2] = p3[id + 2] - p1[id + 2];

        v3[0] = p4[id] - p1[id];
        v3[1] = p4[id + 1] - p1[id + 1];
        v3[2] = p4[id + 2] - p1[id + 2];

        crossProduct<floatingpoint>(cp, v1, v2);
	    //Normalize before checking the angle.
	    normalizeVector(cp);
	    normalizeVector(v3);
        auto retVal = areEqualLT(dotProduct(v3, cp), (floatingpoint)0.0);
        return retVal;
    }

    /// Function to move bead out of plane by specified amount
    inline vector<floatingpoint> movePointOutOfPlane(const vector<floatingpoint> &p1,
                                              const vector<floatingpoint> &p2,
                                              const vector<floatingpoint> &p3,
                                              const vector<floatingpoint> &p4,
                                              int i, floatingpoint d) {
        vector<floatingpoint> norm;
        vector<floatingpoint> v1;
        vector<floatingpoint> v2;

        //get plane
        v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        v2 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};

        norm = normalizeVector(crossProduct(v1, v2));

        //move bead 1
        if (i == 1) {
            vector<floatingpoint> newP1;
            newP1.push_back(p1[0] + norm[0] * d);
            newP1.push_back(p1[1] + norm[1] * d);
            newP1.push_back(p1[2] + norm[2] * d);
            return newP1;
        }

            //move bead 2
        else if (i == 2) {
            vector<floatingpoint> newP2;
            newP2.push_back(p2[0] + norm[0] * d);
            newP2.push_back(p2[1] + norm[1] * d);
            newP2.push_back(p2[2] + norm[2] * d);
            return newP2;
        }

            //move bead 3
        else if (i == 3) {
            vector<floatingpoint> newP3;
            newP3.push_back(p3[0] + norm[0] * d);
            newP3.push_back(p3[1] + norm[1] * d);
            newP3.push_back(p3[2] + norm[2] * d);
            return newP3;
        }

            //move bead 4
        else {
            vector<floatingpoint> newP4;
            newP4.push_back(p4[0] + norm[0] * d);
            newP4.push_back(p4[1] + norm[1] * d);
            newP4.push_back(p4[2] + norm[2] * d);
            return newP4;
        }
    }

    ///IN-PLACE ARRAY VERSION
    inline void movePointOutOfPlane(floatingpoint *p1,
                                    floatingpoint *p2,
                                    floatingpoint *p3,
                                    floatingpoint *p4,floatingpoint *newp,
                                    int i, floatingpoint d) {

        floatingpoint *norm = new floatingpoint[3];
        floatingpoint *v1 = new floatingpoint[3];
        floatingpoint *v2 = new floatingpoint[3];

        //get plane
        v1[0] = p2[0] - p1[0];
        v1[1] = p2[1] - p1[1];
        v1[2] = p2[2] - p1[2];

        v2[0] = p3[0] - p1[0];
        v2[1] = p3[1] - p1[1];
        v2[2] = p3[2] - p1[2];

        crossProduct<floatingpoint>(norm, v1, v2);
        normalizeVector(norm);

        //move bead 1
        if (i == 1) {

            newp[0] = (p1[0] + norm[0] * d);
            newp[1] = (p1[1] + norm[1] * d);
            newp[2] = (p1[2] + norm[2] * d);
        }

            //move bead 2
        else if (i == 2) {
            newp[0] = (p2[0] + norm[0] * d);
            newp[1] = (p2[1] + norm[1] * d);
            newp[2] = (p2[2] + norm[2] * d);
        }

            //move bead 3
        else if (i == 3) {
            newp[0] = (p3[0] + norm[0] * d);
            newp[1] = (p3[1] + norm[1] * d);
            newp[2] = (p3[2] + norm[2] * d);
        }

            //move bead 4
        else {
            newp[0] = (p4[0] + norm[0] * d);
            newp[1] = (p4[1] + norm[1] * d);
            newp[2] = (p4[2] + norm[2] * d);
        }
        delete [] norm;
        delete [] v1;
        delete [] v2;
    }

    ///CUDA VERSION
    #ifdef CUDAACCL
    __host__ __device__
    #endif
    inline void movePointOutOfPlane(floatingpoint *p1,
                                    floatingpoint *p2,
                                    floatingpoint *p3,
                                    floatingpoint *p4,
                                    int i, floatingpoint d, int id) {

        floatingpoint norm[3];
        floatingpoint v1[3];
        floatingpoint v2[3];

        //get plane
        v1[0] = p2[id] - p1[id];
        v1[1] = p2[id + 1] - p1[id + 1];
        v1[2] = p2[id + 2] - p1[id + 2];

        v2[0] = p3[id] - p1[id];
        v2[1] = p3[id + 1] - p1[id + 1];
        v2[2] = p3[id + 2] - p1[id + 2];

        crossProduct<floatingpoint>(norm, v1, v2);
        normalizeVector(norm);

        //move bead 1
        if (i == 1) {

            p1[id] = (p1[id] + norm[0] * d);
            p1[id + 1] = (p1[id + 1] + norm[1] * d);
            p1[id + 2] = (p1[id + 2] + norm[2] * d);
        }

            //move bead 2
        else if (i == 2) {
            p2[id] = (p2[id] + norm[id] * d);
            p2[id + 1] = (p2[id + 1] + norm[id + 1] * d);
            p2[id + 2] = (p2[id + 2] + norm[id + 2] * d);
        }

            //move bead 3
        else if (i == 3) {
            p3[id] = (p3[0] + norm[0] * d);
            p3[id + 1] = (p3[id + 1] + norm[1] * d);
            p3[id + 2] = (p3[id + 2] + norm[2] * d);
        }

            //move bead 4
        else {
            p4[id] = (p4[id] + norm[0] * d);
            p4[id + 1] = (p4[id + 1] + norm[1] * d);
            p4[id + 2] = (p4[id + 2] + norm[2] * d);
        }
//        delete [] norm;
//        delete [] v1;
//        delete [] v2;
    }

    /// Function to create a initial branching point and direction, given an
    /// initial normal vector and point.
    /// @param l - the distance of the branch from the original point
    /// @param m - the size of the branch projection
    /// @param theta - the angle of branching
    /// @return a vector describing the initial branching direction and point
    tuple<vector<floatingpoint>, vector<floatingpoint>> branchProjection(                                                                      const vector<floatingpoint> &n,
                                                           const vector<floatingpoint> &p,
                                                           floatingpoint l, floatingpoint m, floatingpoint theta);


    /// Returns true if two vectors (p1->p2 and p3->p4) are in the same plane
    inline bool areInPlane(const vector<floatingpoint> &p1, const vector<floatingpoint> &p2,
                           const vector<floatingpoint> &p3, const vector<floatingpoint> &p4) {

        auto v1 = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
        auto v2 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};
        auto v3 = {p4[0] - p1[0], p4[1] - p1[1], p4[2] - p1[2]};

        //TODO normalize cp and v3
        auto cp = crossProduct(v1, v2);
        return areEqual(dotProduct(v3, cp), 0.0);
    }

    inline size_t blockToSmemZero(int blockSize){return 0.0 * sizeof(floatingpoint);}
    inline size_t blockToSmem(int blockSize){return 12 * blockSize * sizeof(floatingpoint);}
    inline size_t blockToSmemez(int blockSize){return 24 * blockSize * sizeof(floatingpoint);}
    inline size_t blockToSmemF(int blockSize){return 6 * blockSize * sizeof(floatingpoint);}
    inline size_t blockToSmemFB(int blockSize){return 9 * blockSize * sizeof(floatingpoint);}
    inline size_t blockToSmemFB2(int blockSize){return 18 * blockSize * sizeof(floatingpoint);}
    inline size_t blockToSmemFB3(int blockSize){return 3 * blockSize * sizeof(floatingpoint);}


//    inline bool checkNaN_INF(floatingpoint *x, int startpoint, int endpoint){
//    	for(int i = startpoint; i <= endpoint; i++){
//    		if(isnan(x[i])||isinf(x[i]))
//    			return true;
//    	}
//    	return false;
//    }

    template <class dataType>
    inline bool checkNaN_INF(dataType *x, int startpoint, int endpoint){
        for(int i = startpoint; i <= endpoint; i++){
            if(fabs(x[i]) == numeric_limits<floatingpoint>::infinity()||isnan(fabs(x[i]))
            ||isinf(fabs(x[i]))||fabs(x[i])>1e15)
                return true;
        }
        return false;
    }

    template <class dataType>
    inline void printvariablebinary(dataType *x, int startpoint, int endpoint){
    	string str;
    	for(int i = startpoint; i <= endpoint; i++){
    		str.clear();
		    union { float f; uint32_t i; } u;
		    u.f = x[i];
		    for (int i = 0; i < 32; i++)
		    {
			    if (u.i % 2)  str.push_back('1');
			    else str.push_back('0');
			    u.i >>= 1;
		    }

		    // Reverse the string since now it's backwards
		    string temp(str.rbegin(), str.rend());
		    str = temp;
		    cout<<str<<" ";
    	}
    	cout<<endl;
    }


    /// Function to move bead out of plane by specified amount
    vector<floatingpoint> movePointOutOfPlane(const vector<floatingpoint> &p1,
                                       const vector<floatingpoint> &p2,
                                       const vector<floatingpoint> &p3,
                                       const vector<floatingpoint> &p4,
                                       int i, floatingpoint d);


    float delGGenChem(float delGZero, vector<species_copy_t> reacN, vector<int> reacNu, vector<species_copy_t> prodN, vector<int> prodNu);

    float delGGenChemI(float delGZero, vector<species_copy_t> reacN, vector<int> reacNu, vector<species_copy_t> prodN, vector<int> prodNu);

    float delGDifChem(species_copy_t reacN ,species_copy_t prodN);


    float delGPolyChem(float delGzero, species_copy_t reacN, string whichWay);

    float delGMyoChem(float nh, float rn);



    // Auxiliary type to dispatch distance squared for FENE function.
    template< typename Float > struct FeneDistSq { Float value = 0; };
    template< typename Float > FeneDistSq(Float) -> FeneDistSq< Float >;

    // Energy of Finite Extensible Nonlinear Elastic (FENE) potential.
    //
    // FENE potential:
    //     E = -(1/2) * k * rmax^2 * ln(1 - (r - r0)^2 / rmax^2).
    template<
        typename Float1,
        typename Float2,
        typename Float3,
        std::enable_if_t<
            std::is_floating_point_v<Float1>
            && std::is_floating_point_v<Float2>
            && std::is_floating_point_v<Float3>
        >* = nullptr
    >
    inline auto fene(FeneDistSq<Float1> dist2, Float2 k, Float3 rmax) {
        using ResType = std::common_type_t< Float1, Float2, Float3 >;
        const auto rmax2 = rmax * rmax;
        const auto inLog = std::max<ResType>(0, 1 - dist2.value / rmax2);
        return -k * rmax2 * std::log(inLog) / 2;
    }

    // Part of derivative of Finite Extensible Nonlinear Elastic (FENE) potential, wrt dist (or, equivalently, r - r0).
    // The full derivative is the result of this function, multiplied by dist.
    template<
        typename Float1,
        typename Float2,
        typename Float3,
        std::enable_if_t<
            std::is_floating_point_v<Float1>
            && std::is_floating_point_v<Float2>
            && std::is_floating_point_v<Float3>
        >* = nullptr
    >
    inline auto dFeneCoeff(FeneDistSq<Float1> dist2, Float2 k, Float3 rmax) {
        using ResType = std::common_type_t< Float1, Float2, Float3 >;
        const auto rmax2 = rmax * rmax;
        const auto inLog = std::max<ResType>(0, 1 - dist2.value / rmax2);
        return k / inLog;
    }

} // namespace mathfunc


// Ceiling integer division implementation of https://stackoverflow.com/a/22417111/7120360.
template< typename Int1, typename Int2, std::enable_if_t<std::is_integral_v<Int1> && std::is_signed_v<Int1> && std::is_integral_v<Int2> && std::is_signed_v<Int2>>* = nullptr >
constexpr auto ceildiv(Int1 a, Int2 b) {
    return a / b + ((a % b != 0) ? !((a > 0) ^ (b > 0)) : 0);
}

} // namespace medyan

#endif
