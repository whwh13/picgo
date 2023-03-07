
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
//-----------------------------------------------------------------

#include "Compartment.h"

#include <catch2/catch.hpp>

#include "Controller/GController.h"
#include "MathFunctions.h"
#include "Visitor.h"

#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Util/Math/CuboidSlicing.hpp"
#include "Filament.h"
#include "Cylinder.h"
#include <stdint.h>

namespace medyan {
using namespace mathfunc;

namespace {

// Get packed integer containing filament ID and cylinder index (position from minus end, starting 0).
//
// The result is needed for SIMD neighbor list search to filter out filament binding on nearby cylinders on the same filament.
//
// Parameters
//   - filId: Filament ID.
//   - cylIndex: Cylinder's position from minus end of the filament, starting zero.
//   - allowSameFilamentBinding: If set to false, all cylinder indices will be treated as zero.
inline std::uint32_t packFilamentIdAndCylinderIndexForSimd(std::uint32_t filId, std::uint32_t cylIndex, bool allowSameFilamentBinding) {
    // Precondition check.
    //---------------------------------
    constexpr int numBitsCylIndex = 7;                     // Lower bits
    constexpr int numBitsFilId = 32 - numBitsCylIndex;     // Higher bits
    constexpr std::uint32_t maxFilId = (1 << numBitsFilId) - 1;
    constexpr std::uint32_t maxCylIndex = (1 << numBitsCylIndex) - 1 - ChemParams::minCylinderDistanceSameFilament;

    if(allowSameFilamentBinding && cylIndex > maxCylIndex) {
        throw std::runtime_error("Cylinder index is greater than maximum allowed.");
    }
    if(filId > maxFilId) {
        throw std::runtime_error("Filament ID is greater than maximum allowed.");
    }

    // Actual computation.
    //---------------------------------
    return (filId << numBitsCylIndex) | (allowSameFilamentBinding ? cylIndex : 0);
}

TEST_CASE("Filament ID and cylinder index packing for SIMD", "[SimdUtility]") {
    // The result depends on some internal implementation of the function. If the constants used by the function changes, the expected results should change as well.
    CHECK(packFilamentIdAndCylinderIndexForSimd(0, 125, true) == 125);
    CHECK_THROWS(packFilamentIdAndCylinderIndexForSimd(0, 126, true));
    CHECK_THROWS(packFilamentIdAndCylinderIndexForSimd(1 << 26, 0, true));

    CHECK(packFilamentIdAndCylinderIndexForSimd(1, 20, true) == 148);
    CHECK(packFilamentIdAndCylinderIndexForSimd(1, 20, false) == 128);
}

} // namespace


#ifdef SIMDBINDINGSEARCH

void Compartment::SIMDcoordinates_section(){

    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    int N = getCylinders().size() * SysParams::Chemistry().maxbindingsitespercylinder;
    if(N) {
//        Cyldcindexvec.resize(_cylinders.size());
        CylcIDvec.resize(getCylinders().size());

        short _filamentType = 0;
        uint32_t _fID = 0;
        uint32_t _fpos = 0;
        bool checkftype = false;
        if (SysParams::Chemistry().numFilaments > 1)
            checkftype = true;
        unsigned int i = 0;

        bscoords_section.resize(SysParams::Chemistry().numFilaments * 27);
        //Loop through all filament Types
        for(short filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
          //Loop through cylinders in the compartment
            for (auto cyl:getCylinders()) {
              //clear partitioned coordiante vectors
                for(short i =0; i < 27; i++) {
                    partitionedcoordx[i].clear();
                    partitionedcoordy[i].clear();
                    partitionedcoordz[i].clear();
                    cindex_bs_section[i].clear();
                    finfo_bs_section[i].clear();
                }

                int cindex = cyl->getStableIndex();

                _filamentType = Cylinder::getDbDataConst()[cindex].type;
                _fID = Cylinder::getDbDataConst()[cindex].filamentId;
                _fpos = Cylinder::getDbDataConst()[cindex].positionOnFilament -
                        Cylinder::getDbDataConst()[cindex].filamentFirstEntry;

                // Packed integer containing filament ID and cylinder index.
                const auto cylfinfo = packFilamentIdAndCylinderIndexForSimd(_fID, _fpos, SysParams::Chemistry().allowPairBindingOnSameFilament);

                //Only consider cylinders that are filType
                if (checkftype && _filamentType != filType) continue;

                auto x1 = cyl->getFirstBead()->vcoordinate();
                auto x2 = cyl->getSecondBead()->vcoordinate();
//            uint32_t shiftedindex = (i << SysParams::Chemistry().shiftbybits);
                uint32_t shiftedindex = (cyl->getStableIndex() << SysParams::Chemistry().shiftbybits);

//                Cyldcindexvec[i] = cyl->_dcIndex;
                CylcIDvec[i] = cyl->getId();
                uint32_t j = 0;
                if(cyl->isMinusEnd() == false) {
                    float cylsizesquared = SysParams::Geometry().cylinderSize[_filamentType]
                                           *
                                           SysParams::Geometry().cylinderSize[_filamentType];
                    float maxmp = sqrt(twoPointDistancesquared(x1, x2) / cylsizesquared);

                    for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                         it !=
                         SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {

                        auto mp = (float) *it /
                                  SysParams::Geometry().cylinderNumMon[_filamentType];

                        if (mp <= maxmp) {
                            auto coord = midPointCoordinate(x1, x2, cyl->adjustedrelativeposition(mp));
//                            auto coord = midPointCoordinate(x1, x2, mp);
                            //last 4 bits are binding site while first 12 bits are cylinder index.
                            uint32_t index = shiftedindex | j;
                            int pindices[3];
                            getpartition3Dindex(pindices, coord);
                            addcoordtopartitons(pindices, coord, index, cylfinfo);
                        }
                        j++;
                    }
                }
                else{
                    /* If it is the minus end Cylinder, add the binding sites that are
                     * species Filament*/
                    for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                         it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                         auto sf = Cylinder::getDbDataConst()[cindex]
                                .chemCylinder->getCMonomer(*it)->activeSpeciesFilament();
                         if(sf !=-1){
                             auto mp = (float) *it /
                                       SysParams::Geometry().cylinderNumMon[_filamentType];
                             auto coord = midPointCoordinate(x1, x2, cyl->adjustedrelativeposition(mp));
//                             auto coord = midPointCoordinate(x1, x2, mp);
                             //last 4 bits are binding site while first 12 bits are cylinder index.
                             uint32_t index = shiftedindex | j;
                             int pindices[3];
                             getpartition3Dindex(pindices, coord);
                             addcoordtopartitons(pindices, coord, index, cylfinfo);
                         }
                         j++;
                    }
                }
            }

//    assert(k<65536);
            //Create input vector for SIMD calculations
//        cout<<bscoords.size()<<" "<<partitionedcoordx[0].size()<<endl;

            for (short i = 0; i < 27; i++) {
//            cout<<partitionedcoordx[i].size()<<" ";
                bscoords_section[filType*27 + i].init_coords(partitionedcoordx[i],
                        partitionedcoordy[i], partitionedcoordz[i], cindex_bs_section[i], finfo_bs_section[i]);
            }
//        cout<<endl;
        }
    }
}

void Compartment::SIMDcoordinates4linkersearch_section(bool isvectorizedgather){

	if(areEqual(medyan::HybridBindingSearchManager::largestlinkerdistance, 0.0)) return;

    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 1;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    floatingpoint _rMax = medyan::HybridBindingSearchManager::largestlinkerdistance;
    floatingpoint searchdist = SysParams::Geometry().largestCylinderSize/2 + _rMax;
    bool rMaxvsCmpSize = (SysParams::Geometry().largestCylinderSize/2 + _rMax) <
                            SysParams::Geometry().largestCompartmentSide/2;

    floatingpoint coord_bounds[6] = {
            centerCoord[0] - SysParams::Geometry().compartmentSizeX/2 + searchdist,
            centerCoord[1] - SysParams::Geometry().compartmentSizeY/2 + searchdist,
            centerCoord[2] - SysParams::Geometry().compartmentSizeZ/2 + searchdist,
            centerCoord[0] + SysParams::Geometry().compartmentSizeX/2 - searchdist,
            centerCoord[1] + SysParams::Geometry().compartmentSizeY/2 - searchdist,
            centerCoord[2] + SysParams::Geometry().compartmentSizeZ/2 - searchdist};

    int N = getCylinders().size() * maxnbs;
    //Paritioned coordinates are created for each unique filamentType
	  bscoords_section_linker.resize(SysParams::Chemistry().numFilaments * 27);

    for (short filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
        if(N) {
            for(short i =0; i < 27; i++) {
                partitionedcoordx[i].clear();
                partitionedcoordy[i].clear();
                partitionedcoordz[i].clear();
                cindex_bs_section[i].clear();
                finfo_bs_section[i].clear();
            }

//            Cyldcindexvec.resize(_cylinders.size());
            short _filamentType = 0;
            uint32_t _fID = 0;
            uint32_t _fpos = 0;


            bool checkftype = false;
            if (SysParams::Chemistry().numFilaments > 1)
                        checkftype = true;
//            unsigned int i = 0;
            for (auto cyl:getCylinders()) {
                uint32_t cindex = cyl->getStableIndex();

                _filamentType = Cylinder::getDbData()[cindex].type;
                _filamentType = Cylinder::getDbDataConst()[cindex].type;
                _fID = Cylinder::getDbDataConst()[cindex].filamentId;
                _fpos = Cylinder::getDbDataConst()[cindex].positionOnFilament-
                        Cylinder::getDbDataConst()[cindex].filamentFirstEntry;

                // Packed integer containing filament ID and cylinder index.
                const auto cylfinfo = packFilamentIdAndCylinderIndexForSimd(_fID, _fpos, SysParams::Chemistry().allowPairBindingOnSameFilament);

                //Consider only cylinders of filamentType fType
                if (checkftype && _filamentType != filType) continue;

                auto x1 = cyl->getFirstBead()->vcoordinate();
                auto x2 = cyl->getSecondBead()->vcoordinate();
//            uint16_t shiftedindex = (i << SysParams::Chemistry().shiftbybits);
                uint32_t shiftedindex = (cindex << SysParams::Chemistry().shiftbybits);
//                Cyldcindexvec[i] = cindex;
//                i++;
                uint32_t j = 0;
                int dBI = SysParams::Chemistry().linkerbindingskip-1;
                if(cyl->isMinusEnd() == false) {
                    for (int itI = 0; itI < SysParams::Chemistry().bindingSites[_filamentType].size(); itI += dBI) {

                        auto it = SysParams::Chemistry().bindingSites[_filamentType].begin() + itI;

                        bool state = false;
                        if (isvectorizedgather)
                            state = checkoccupancy(boundstate, bstatepos,
                                                   maxnbs * cindex + j);
                        else
                            state = checkoccupancy(cyl, *it, _filamentType, bstatepos);
                        if (state) {
                            auto mp = (float) *it /
                                      SysParams::Geometry().cylinderNumMon[_filamentType];
                            float cylsizesquared =
                                    SysParams::Geometry().cylinderSize[_filamentType]
                                    * SysParams::Geometry().cylinderSize[_filamentType];
                            float maxmp = sqrt(
                                    twoPointDistancesquared(x1, x2) / cylsizesquared);
                            if (mp <= maxmp) {
                                auto coord = midPointCoordinate(x1, x2, cyl->adjustedrelativeposition(mp));
//                                auto coord = midPointCoordinate(x1, x2, mp);
                                //last 4 bits are binding site while first 12 bits are cylinder index.
                                uint32_t index = shiftedindex | j;
                                int pindices[3];
                                if (rMaxvsCmpSize) {
                                    getpartitionindex<true>(pindices, coord, coord_bounds);
                                    addcoordtorMaxbasedpartitons<true>(pindices, coord,
                                                                       index,
                                                                       cylfinfo);
                                } else {
                                    getpartitionindex<false>(pindices, coord, coord_bounds);
                                    addcoordtorMaxbasedpartitons<false>(pindices, coord,
                                                                        index,
                                                                        cylfinfo);
                                }
                            }
                        }
                        j+=dBI;
                    }
                }
                else{
                    /* If it is the minus end Cylinder, add the binding sites that are
                     * species Filament*/
                    for (int itI = 0; itI < SysParams::Chemistry().bindingSites[_filamentType].size(); itI += dBI) {

                        auto it = SysParams::Chemistry().bindingSites[_filamentType].begin() + itI;

                        auto sf = Cylinder::getDbDataConst()[cindex]
                                .chemCylinder->getCMonomer(*it)->activeSpeciesFilament();
                        if(sf !=-1){
                            bool state = false;
                            if (isvectorizedgather)
                                state = checkoccupancy(boundstate, bstatepos,
                                                       maxnbs * cindex + j);
                            else
                                state = checkoccupancy(cyl, *it, _filamentType, bstatepos);
                            if (state) {
                                auto mp = (float) *it /
                                          SysParams::Geometry().cylinderNumMon[_filamentType];
                                auto coord = midPointCoordinate(x1, x2, cyl->adjustedrelativeposition(mp));
//                                auto coord = midPointCoordinate(x1, x2, mp);
                                //last 4 bits are binding site while first 12 bits are cylinder index.
                                uint32_t index = shiftedindex | j;
                                int pindices[3];
                                if (rMaxvsCmpSize) {
                                    getpartitionindex<true>(pindices, coord, coord_bounds);
                                    addcoordtorMaxbasedpartitons<true>(pindices, coord,
                                                                       index,
                                                                       cylfinfo);
                                } else {
                                    getpartitionindex<false>(pindices, coord, coord_bounds);
                                    addcoordtorMaxbasedpartitons<false>(pindices, coord,
                                                                        index,
                                                                        cylfinfo);
                                }
                            }
                        }
                        j+=dBI;
                    }
                }
            }
//        cout<<endl;
//    assert(k<65536);
            //Create input vector for SIMD calculations
//        cout<<"Linker coord size ";
            for (short i = 0; i < 27; i++) {
                bscoords_section_linker[filType * 27 + i].init_coords(partitionedcoordx[i],
                                                                      partitionedcoordy[i],
                                                                      partitionedcoordz[i],
                                                                      cindex_bs_section[i],
                                                                      finfo_bs_section[i]);
            }
        }//if(N)
        else{
            for (short i = 0; i < 27; i++) {
                    bscoords_section_linker[filType * 27 + i].resize(0);
            }
        }
    }
}

void Compartment::SIMDcoordinates4motorsearch_section(bool isvectorizedgather){

	if(areEqual(medyan::HybridBindingSearchManager::largestmotordistance,0.0)) return;

    for(short i =0; i < 27; i++) {
        partitionedcoordx[i].clear();
        partitionedcoordy[i].clear();
        partitionedcoordz[i].clear();
        cindex_bs_section[i].clear();
        finfo_bs_section[i].clear();
    }

    //setting size to the number of maximum binding sites per cylinder * number of
    // cylinders in compartment.
    short bstatepos = 2;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    short maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;

    floatingpoint _rMax = medyan::HybridBindingSearchManager::largestmotordistance;
    floatingpoint searchdist = SysParams::Geometry().largestCylinderSize/2 + _rMax;
    bool rMaxvsCmpSize = (searchdist) <
                         SysParams::Geometry().largestCompartmentSide/2;

//    cout<<"rMaxvsCmpSize "<<rMaxvsCmpSize<<endl;
    floatingpoint coord_bounds[6] = {
            centerCoord[0] - SysParams::Geometry().compartmentSizeX/2 + searchdist,
            centerCoord[1] - SysParams::Geometry().compartmentSizeY/2 + searchdist,
            centerCoord[2] - SysParams::Geometry().compartmentSizeZ/2 + searchdist,
            centerCoord[0] + SysParams::Geometry().compartmentSizeX/2 - searchdist,
            centerCoord[1] + SysParams::Geometry().compartmentSizeY/2 - searchdist,
            centerCoord[2] + SysParams::Geometry().compartmentSizeZ/2 - searchdist};

    int N = getCylinders().size() * maxnbs;
	bscoords_section_motor.resize(SysParams::Chemistry().numFilaments * 27);
    for (short filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
        if (N) {
            for(short i =0; i < 27; i++) {
                partitionedcoordx[i].clear();
                partitionedcoordy[i].clear();
                partitionedcoordz[i].clear();
                cindex_bs_section[i].clear();
                finfo_bs_section[i].clear();
            }

//            Cyldcindexvec.resize(_cylinders.size());

            short _filamentType = 0;
            uint32_t _fID = 0;
            uint32_t _fpos = 0;

            bool checkftype = false;
            if (SysParams::Chemistry().numFilaments > 1)
                checkftype = true;
            unsigned int i = 0;

            for (auto cyl:getCylinders()) {
                uint32_t cindex = cyl->getStableIndex();

                _filamentType = Cylinder::getDbDataConst()[cindex].type;
                _fID = Cylinder::getDbDataConst()[cindex].filamentId;
                _fpos = Cylinder::getDbDataConst()[cindex].positionOnFilament-
                        Cylinder::getDbDataConst()[cindex].filamentFirstEntry;

                // Packed integer containing filament ID and cylinder index.
                const auto cylfinfo = packFilamentIdAndCylinderIndexForSimd(_fID, _fpos, SysParams::Chemistry().allowPairBindingOnSameFilament);

                if (checkftype && _filamentType != filType) continue;

                auto x1 = cyl->getFirstBead()->vcoordinate();
                auto x2 = cyl->getSecondBead()->vcoordinate();
//                uint16_t shiftedindex = (i << 4);
                uint32_t shiftedindex = (cindex << SysParams::Chemistry().shiftbybits);
//                Cyldcindexvec[i] = cindex;
                i++;
                uint32_t j = 0;
                if(cyl->isMinusEnd() == false) {
                    for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                         it !=
                         SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                        bool state = false;
                        if (isvectorizedgather)
                            state = checkoccupancy(boundstate, bstatepos,
                                                   maxnbs * cindex + j);
                        else
                            state = checkoccupancy(cyl, *it, _filamentType, bstatepos);
                        if (state) {
                            auto mp = (float) *it /
                                      SysParams::Geometry().cylinderNumMon[_filamentType];
                            auto coord = midPointCoordinate(x1, x2, cyl->adjustedrelativeposition(mp));
//                            auto coord = midPointCoordinate(x1, x2, mp);
                            //last 4 bits are binding site while first 12 bits are cylinder index.
                            uint32_t index = shiftedindex | j;
                            //split and crosscheck
//						cout<<index<<" ";
                            int pindices[3];
                            if (rMaxvsCmpSize) {
                                getpartitionindex<true>(pindices, coord, coord_bounds);
                                addcoordtorMaxbasedpartitons<true>(pindices, coord, index,
                                                                   cylfinfo);
                            } else {
                                getpartitionindex<false>(pindices, coord, coord_bounds);
                                addcoordtorMaxbasedpartitons<false>(pindices, coord, index,
                                                                    cylfinfo);
                            }
/*                        getpartition3Dindex(pindices, coord, coord_bounds);
                        addcoordtopartitons_smallrmax(pindices, coord, index);
                        getpartition3Dindex(pindices, coord);
                        addcoordtopartitons(pindices, coord, index);*/
                        }
                        j++;
                    }
                }
                else{
                    /* If it is the minus end Cylinder, add the binding sites that are
                     * species Filament*/
                    for (auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
                         it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
                        auto sf = Cylinder::getDbDataConst()[cindex]
                                .chemCylinder->getCMonomer(*it)->activeSpeciesFilament();
                        if(sf !=-1){
                            bool state = false;
                            if (isvectorizedgather)
                                state = checkoccupancy(boundstate, bstatepos,
                                                       maxnbs * cindex + j);
                            else
                                state = checkoccupancy(cyl, *it, _filamentType, bstatepos);
                            if (state) {
                                auto mp = (float) *it /
                                          SysParams::Geometry().cylinderNumMon[_filamentType];
                                auto coord = midPointCoordinate(x1, x2, cyl->adjustedrelativeposition(mp));
//                                auto coord = midPointCoordinate(x1, x2, mp);
                                //last 4 bits are binding site while first 12 bits are cylinder index.
                                uint32_t index = shiftedindex | j;
                                int pindices[3];
                                if (rMaxvsCmpSize) {
                                    getpartitionindex<true>(pindices, coord, coord_bounds);
                                    addcoordtorMaxbasedpartitons<true>(pindices, coord,
                                                                       index,
                                                                       cylfinfo);
                                } else {
                                    getpartitionindex<false>(pindices, coord, coord_bounds);
                                    addcoordtorMaxbasedpartitons<false>(pindices, coord,
                                                                        index,
                                                                        cylfinfo);
                                }
                            }
                        }
                        j++;
                    }
                }
            }
//            cout<<endl;
//    assert(k<65536);
            //Create input vector for SIMD calculations
//        cout<<"Motor coords size ";
            for (short i = 0; i < 27; i++) {
//            cout<<partitionedcoordx[i].size()<<" ";
                bscoords_section_motor[filType * 27 + i].init_coords(partitionedcoordx[i],
                                                                     partitionedcoordy[i],
                                                                     partitionedcoordz[i],
                                                                     cindex_bs_section[i],
                                                                     finfo_bs_section[i]);
            }

        } //if(N)
        else {
            for (short i = 0; i < 27; i++)
                bscoords_section_motor[filType * 27 + i].resize(0);
        }
    }
}

void Compartment::getpartition3Dindex(int (&indices)[3], vector<floatingpoint> coord){
    short i = 0;

    for(auto x:coord){
        indices[i] =(x > centerCoord[i]);
        i++;
    }
}

//if rMax+Cylsize/2+delta is less than CmpSize/2
template<>
void Compartment::getpartitionindex<true>(int (&indices)[3], vector<floatingpoint> coord,
                       floatingpoint (&cmpcornercoords)[6]){
    short i = 0;

    for(auto x:coord){
        if(x <= cmpcornercoords[i])
            indices[i] = 0;
        else if(x>=cmpcornercoords[i+3])
            indices[i] = 2;
        else
            indices[i] = 1;

        i++;
    }
}

//if rMax+Cylsize/2+delta is greater than CmpSize/2
template<>
void Compartment::getpartitionindex<false>(int (&indices)[3], vector<floatingpoint> coord,
                                          floatingpoint (&cmpcornercoords)[6]){
    short i = 0;

    for(auto x:coord){
        if(x >= cmpcornercoords[i])
            indices[i] = 2;
        else if(x<=cmpcornercoords[i+3])
            indices[i] = 0;
        else
            indices[i] = 1;

        i++;
    }
}

void Compartment::addcoord(vector<floatingpoint> coord, uint32_t index, uint32_t
cylfinfo, short i){
    partitionedcoordx[i].push_back(coord[0]);
    partitionedcoordy[i].push_back(coord[1]);
    partitionedcoordz[i].push_back(coord[2]);
    cindex_bs_section[i].push_back(index);
    finfo_bs_section[i].push_back(cylfinfo);
}

void Compartment::deallocateSIMDcoordinates(){
	for(short i =0; i < 27; i++) {
		partitionedcoordx[i].clear();
		partitionedcoordy[i].clear();
		partitionedcoordz[i].clear();
		cindex_bs_section[i].clear();
        finfo_bs_section[i].clear();
		if(bscoords_section_linker.size() > 0) bscoords_section_linker[i].resize(0);
		if(bscoords_section_motor.size() > 0) bscoords_section_motor[i].resize(0);
	}

}

bool Compartment::checkoccupancy(vector<vector<bool>>& boundstate, short bstatepos, int pos){
    return boundstate[bstatepos][pos];
}

bool Compartment::checkoccupancy(Cylinder* cyl, short it, short _filamentType,
                                 short bstatepos){
    bool status;
    if(bstatepos == 1)
        status = cyl->getCCylinder()->getCMonomer(it)->speciesBound(
            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN() == 1.0;
    else if(bstatepos == 2)
        status = cyl->getCCylinder()->getCMonomer(it)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN() == 1.0;
    return status;
}

void Compartment::addcoordtopartitons(int (&pindices)[3], vector<floatingpoint> coord, uint32_t
index, uint32_t cylfinfo){
    addcoord(coord, index, cylfinfo, 0);

    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 26);
    }
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 24);
    }
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
    }
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 21);
    }
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 14);
        //Vertex
        addcoord(coord, index, cylfinfo, 22);
    }
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);

    }
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 8);
        //Vertex
        addcoord(coord, index, cylfinfo, 23);
    }
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 13);
        //Vertex
        addcoord(coord, index, cylfinfo, 25);
    }
}

void Compartment::addcoordtopartitons_smallrmax(int (&pindices)[3], vector<floatingpoint> coord,
                                                uint16_t index, uint32_t cylfinfo){
    addcoord(coord, index, cylfinfo, 0);
    //111
    if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1) {
        return;
    }
    //000
    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 26);
    }
        //001
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
    }
        //002
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 24);
    }
        //010
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
    }
        //011
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        //Vertex
    }
        //012
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
    }
        //020
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
    }
        //021
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        //Vertex
    }
        //022
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 21);
    }
        //100
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 14);
        //Vertex
    }
        //101
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        //Vertex
    }
        //102
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 12);
        //Vertex
    }
        //110
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        //Vertex
    }
        //112
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        //Vertex
    }
        //120
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 11);
        //Vertex
    }
        //121
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        //Vertex
    }
        //122
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 13);
        //Vertex
    }
        //200
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 14);
        //Vertex
        addcoord(coord, index, cylfinfo, 22);
    }
        //201
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        //Vertex
    }
        //202
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);

    }
        //210
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        //Vertex

    }
        //211
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        //Vertex
    }
        //212
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
    }
        //220
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 8);
        //Vertex
        addcoord(coord, index, cylfinfo, 23);
    }
        //221
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
    }
        //222
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 13);
        //Vertex
        addcoord(coord, index, cylfinfo, 25);
    }
}

//there are 27 partitions made based on rMax. These partitions are composites of
// underlying 27 unique voxels that are non-verlapping.

//if rMax+Cylsize/2+delta is lesser than CmpSize/2
template<>
void Compartment::addcoordtorMaxbasedpartitons<true>(int (&pindices)[3], vector<floatingpoint>
coord, uint32_t index, uint32_t cylfinfo){
    addcoord(coord, index, cylfinfo, 0);
    //111
    if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1) {
        return;
    }
    //000
    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 26);
    }
        //001
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
    }
        //002
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 24);
    }
        //010
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
    }
        //011
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        //Vertex
    }
        //012
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
    }
        //020
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
    }
        //021
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        //Vertex
    }
        //022
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 21);
    }
        //100
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 14);
        //Vertex
    }
        //101
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        //Vertex
    }
        //102
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 12);
        //Vertex
    }
        //110
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        //Vertex
    }
        //112
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        //Vertex
    }
        //120
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 11);
        //Vertex
    }
        //121
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        //Vertex
    }
        //122
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 13);
        //Vertex
    }
        //200
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 14);
        //Vertex
        addcoord(coord, index, cylfinfo, 22);
    }
        //201
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        //Vertex
    }
        //202
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);

    }
        //210
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        //Vertex

    }
        //211
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        //Vertex
    }
        //212
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
    }
        //220
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 8);
        //Vertex
        addcoord(coord, index, cylfinfo, 23);
    }
        //221
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
    }
        //222
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 13);
        //Vertex
        addcoord(coord, index, cylfinfo, 25);
    }
}

//if rMax+Cylsize/2+delta is greater than CmpSize/2
template<>
void Compartment::addcoordtorMaxbasedpartitons<false>(int (&pindices)[3], vector<floatingpoint>
coord, uint32_t index, uint32_t cylfinfo){

    addcoord(coord, index, cylfinfo, 0);
    //111
    if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 1) {
        for(int part = 1; part < 27; part++){
            addcoord(coord, index, cylfinfo, part);
        }
        return;
    }
    //000
    if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 26);
    }
        //001
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 24);
        addcoord(coord, index, cylfinfo, 26);
    }
        //002
    else if(pindices[0] ==0 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 18);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 24);
    }
        //010
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
        addcoord(coord, index, cylfinfo, 26);
    }
        //011
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
        addcoord(coord, index, cylfinfo, 21);
        addcoord(coord, index, cylfinfo, 24);
        addcoord(coord, index, cylfinfo, 26);
    }
        //012
    else if(pindices[0] ==0 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 24);
        addcoord(coord, index, cylfinfo, 21);
    }
        //020
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
    }
        //021
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 5);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 16);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
        addcoord(coord, index, cylfinfo, 21);
    }
        //022
    else if(pindices[0] ==0 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 7);
        //Vertex
        addcoord(coord, index, cylfinfo, 21);
    }
        //100
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 22);
        addcoord(coord, index, cylfinfo, 26);
    }
        //101
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);
        addcoord(coord, index, cylfinfo, 22);
        addcoord(coord, index, cylfinfo, 24);
        addcoord(coord, index, cylfinfo, 26);
    }
        //102
    else if(pindices[0] ==1 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);
        addcoord(coord, index, cylfinfo, 24);
    }
        //110
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
        addcoord(coord, index, cylfinfo, 22);
        addcoord(coord, index, cylfinfo, 23);
        addcoord(coord, index, cylfinfo, 26);
    }
        //112
    else if(pindices[0] ==1 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 18);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);
        addcoord(coord, index, cylfinfo, 21);
        addcoord(coord, index, cylfinfo, 24);
        addcoord(coord, index, cylfinfo, 25);
    }
        //120
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
        addcoord(coord, index, cylfinfo, 23);
    }
        //121
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 16);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 19);
        addcoord(coord, index, cylfinfo, 21);
        addcoord(coord, index, cylfinfo, 23);
        addcoord(coord, index, cylfinfo, 25);
    }
        //122
    else if(pindices[0] ==1 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 5);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 7);
        addcoord(coord, index, cylfinfo, 9);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 21);
        addcoord(coord, index, cylfinfo, 25);
    }
        //200
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 14);
        //Vertex
        addcoord(coord, index, cylfinfo, 22);
    }
        //201
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);
        addcoord(coord, index, cylfinfo, 22);
    }
        //202
    else if(pindices[0] ==2 && pindices[1] == 0 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 15);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);

    }
        //210
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 22);
        addcoord(coord, index, cylfinfo, 23);
    }
        //211
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 14);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);
        addcoord(coord, index, cylfinfo, 22);
        addcoord(coord, index, cylfinfo, 23);
        addcoord(coord, index, cylfinfo, 25);
    }
        //212
    else if(pindices[0] ==2 && pindices[1] == 1 && pindices[2] == 2){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 4);
        addcoord(coord, index, cylfinfo, 6);
        //Edge
        addcoord(coord, index, cylfinfo, 10);
        addcoord(coord, index, cylfinfo, 12);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 20);
        addcoord(coord, index, cylfinfo, 25);
    }
        //220
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 0){
        //Plane
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 2);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 8);
        //Vertex
        addcoord(coord, index, cylfinfo, 23);
    }
        //221
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 1){
        //Plane
        addcoord(coord, index, cylfinfo, 1);
        addcoord(coord, index, cylfinfo, 2);
        addcoord(coord, index, cylfinfo, 6);
        addcoord(coord, index, cylfinfo, 3);
        //Edge
        addcoord(coord, index, cylfinfo, 8);
        addcoord(coord, index, cylfinfo, 11);
        addcoord(coord, index, cylfinfo, 13);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 17);
        //Vertex
        addcoord(coord, index, cylfinfo, 23);
        addcoord(coord, index, cylfinfo, 25);
    }
        //222
    else if(pindices[0] ==2 && pindices[1] == 2 && pindices[2] == 2){
        //Plane
        addcoord(coord, index,cylfinfo,  6);
        addcoord(coord, index, cylfinfo, 3);
        addcoord(coord, index, cylfinfo, 1);
        //Edge
        addcoord(coord, index, cylfinfo, 17);
        addcoord(coord, index, cylfinfo, 15);
        addcoord(coord, index, cylfinfo, 13);
        //Vertex
        addcoord(coord, index, cylfinfo, 25);
    }
}

#endif

Compartment& Compartment::operator=(const Compartment &other) {

    _internal_reactions.clear();
    _diffusion_reactions.clear();
    _species.clear();
    other.cloneSpecies(this);
    other.cloneReactions(this);
    diffusionCoefficients_ = other.diffusionCoefficients_;

    return *this;
}

bool Compartment::apply_impl(SpeciesVisitor &v) {
    for(auto &s : _species.species()) {
        v.visit(s.get());
    }
    return true;
}

bool Compartment::apply_impl(ReactionVisitor &v) {
    for(auto &r : _internal_reactions.reactions()) {
        v.visit(r.get());
    }
    return true;
}

void Compartment::computeNonSlicedVolumeArea() {
    auto sizex = SysParams::Geometry().compartmentSizeX;
    auto sizey = SysParams::Geometry().compartmentSizeY;
    auto sizez = SysParams::Geometry().compartmentSizeZ;

    _partialArea = {{sizey * sizez, sizey * sizez, sizex * sizez, sizex * sizez, sizex * sizey, sizex * sizey}};

    _volumeFrac = 1.0;
}

void Compartment::computeSlicedVolumeArea(SubSystem& sys, SliceMethod sliceMethod) {

    switch(sliceMethod) {
    case SliceMethod::membrane:
        {
            // The calculation requires the
            //  - The position calculation of triangles
            //  - The area calculation of triangles
            //  - The unit normal vector of triangles
            const size_t numTriangle = getTriangles().size();
            if(numTriangle) {
                double sumArea = 0.0;
                Vec< 3, floatingpoint > sumNormal {};
                Vec< 3, floatingpoint > sumPos {};
                for(auto tiSys : getTriangles()) {
                    auto& t = sys.triangles[tiSys];
                    const auto& mesh = t.getParent(sys).getMesh();
                    const auto ti = medyan::Membrane::MeshType::triangleIndex( t.getTopoIndex() );
                    const auto area = mesh.attribute(ti).gTriangle.area;
                    const auto& unitNormal = mesh.attribute(ti).gTriangle.unitNormal;
                    sumNormal += unitNormal * area;
                    sumPos += t.coordinate * area;
                    sumArea += area;
                }
                normalize(sumNormal);
                sumPos /= sumArea;

                auto res = PlaneCuboidSlicer() (
                    sumPos, sumNormal,
                    {
                        centerCoord[0] - SysParams::Geometry().compartmentSizeX * (floatingpoint)0.5,
                        centerCoord[1] - SysParams::Geometry().compartmentSizeY * (floatingpoint)0.5,
                        centerCoord[2] - SysParams::Geometry().compartmentSizeZ * (floatingpoint)0.5
                    },
                    {{
                        SysParams::Geometry().compartmentSizeX,
                        SysParams::Geometry().compartmentSizeY,
                        SysParams::Geometry().compartmentSizeZ
                    }}
                );

                _volumeFrac = res.volumeIn / sys.getCompartmentGrid()->compartmentVolume;
                _partialArea = res.areaIn;
            }
        }
        break;

    case SliceMethod::cylinderBoundary:
        {

            // The calculation requires the
            //  - The position calculation of triangles
            //  - The area calculation of triangles
            //  - The unit normal vector of triangles
            // ASSUMPTIONS:
            //  - This compartment is a CUBE
            //get compartment sizes in X,Y and the radius of cylinder
            auto sizex = SysParams::Geometry().compartmentSizeX;
            auto sizey = SysParams::Geometry().compartmentSizeY;
            auto sizez = SysParams::Geometry().compartmentSizeZ;
            auto r = SysParams::Boundaries().diameter / 2; //radius

            //get geometry center of the compartment
            auto x = centerCoord[0];
            auto y = centerCoord[1];

            auto leftx = x - sizex / 2;
            auto rightx = x + sizex / 2;
            auto lowy = y - sizey / 2;
            auto upy = y + sizey / 2;

            float pleft, pright, plow, pup, lleft, lright, llow, lup;
            vector<float> edge; //left0,right1,low2,up3

            // find intersections for 4 edges of a compartment
            // if no intersection for an edge, write -1
            pleft = r - sqrt(r * r - (leftx - r) * (leftx - r));

            if(pleft > upy || pleft < lowy) {

                pleft = r + sqrt(r * r - (leftx - r) * (leftx - r));

                if(pleft > upy || pleft < lowy) edge.push_back(-1);
                else edge.push_back(pleft);
            }
            else edge.push_back(pleft);


            pright = r - sqrt(r * r - (rightx - r) * (rightx - r));

            if(pright > upy || pright < lowy){

                pright = r + sqrt(r * r - (rightx - r) * (rightx - r));

                if(pright > upy || pright < lowy) edge.push_back(-1);
                else edge.push_back(pright);
            }
            else edge.push_back(pright);


            plow = r - sqrt(r * r - (lowy - r) * (lowy - r));

            if(plow > rightx || plow < leftx){

                plow = r + sqrt(r * r - (lowy - r) * (lowy - r));

                if(plow > rightx || plow < leftx) edge.push_back(-1);
                else edge.push_back(plow);
            }
            else edge.push_back(plow);

            pup = r - sqrt(r * r - (upy - r) * (upy - r));

            if(pup > rightx || pup < leftx){

                pup = r + sqrt(r * r - (upy - r) * (upy - r));

                if(pup > rightx || pup < leftx) edge.push_back(-1);
                else edge.push_back(pup);
            }
            else edge.push_back(pup);

            vector<int> internum;

            for(int i=0; i < 4; i++){
                //if intersections are opposite
                if(edge[i] > -0.5) internum.push_back(i);
            }

            //3 intersections -> two of them are the same vertex, remove the extra one
            if(internum.size() < 4 && internum.size() > 2){

                if(internum[0] == 0 && internum[1] == 1) internum.erase(internum.end()-1);
                else internum.erase(internum.begin());

            }


            float distsq = (x - r) * (x - r) + (y - r) * (y - r);
            float totalVol = sizex * sizey;
            float fraction_i;

            //internum = 0 if no intersection
            if(internum.size() < 1){
                //outside network
                if(distsq > r * r){
                    lleft = 0;
                    lright = 0;
                    llow = 0;
                    lup = 0;
                    _volumeFrac = 0.0;
                }
                    //inside network
                else{
                    lleft = sizey;
                    lright = sizey;
                    llow = sizex;
                    lup = sizex;
                    _volumeFrac = 1.0;
                }
            }
                //internum = 4 if both intersection are vertices
            else if (internum.size() > 3){
                lleft = sizey;
                lright = sizey;
                llow = sizex;
                lup = sizex;
                _volumeFrac = 0.5;
            }
            else {
                //1. intersect left0 and lower2 planes
                if(internum[0] == 0 && internum[1] == 2){
                    fraction_i = 0.5 * (edge[0] - lowy) * (edge[2] - leftx) / totalVol;

                    if(distsq < r * r) {
                        _volumeFrac = 1 - fraction_i;
                        lleft = upy - edge[0];
                        lright = sizey;
                        llow = rightx - edge[2];
                        lup = sizex;
                    }
                    else{
                        _volumeFrac = fraction_i;
                        lleft = edge[0] - lowy;
                        lright = sizey;
                        llow = edge[2] - leftx;
                        lup = sizex;
                    }
                }
                    //2. intersect left0 and right1 planes
                else if(internum[0] == 0 && internum[1] == 1){
                    fraction_i = 0.5 * (edge[0] - lowy + edge[1] - lowy) * sizex / totalVol;

                    if(y > r) {
                        _volumeFrac = fraction_i;
                        lleft = edge[0] - lowy;
                        lright = edge[1] - lowy;
                        llow = sizex;
                        lup = sizex;
                    }
                    else{
                        _volumeFrac = 1 - fraction_i;
                        lleft = upy - edge[0];
                        lright = upy - edge[1];
                        llow = sizex;
                        lup = sizex;
                    }
                }
                    //3. intersect left0 and upper3 planes
                else if(internum[0] == 0 && internum[1] == 3){
                    fraction_i = 0.5 * (upy - edge[0]) * (edge[3] - leftx) / totalVol;

                    if(distsq < r * r) {
                        _volumeFrac = 1 - fraction_i;
                        lleft = edge[0] - lowy;
                        lright = sizey;
                        llow = sizex;
                        lup = rightx - edge[3];
                    }
                    else{
                        _volumeFrac = fraction_i;
                        lleft = upy - edge[0];
                        lright = sizey;
                        llow = sizex;
                        lup = edge[3] - leftx;
                    }
                }
                    //4. intersect lower2 and right1 planes
                else if(internum[0] == 1 && internum[1] == 2){
                    fraction_i = 0.5 * (edge[1] - lowy) * (rightx - edge[2]) / totalVol;

                    if(distsq < r * r) {
                        _volumeFrac = 1 - fraction_i;
                        lleft = sizey;
                        lright = upy - edge[1];
                        llow = edge[2] - leftx;
                        lup = sizex;
                    }
                    else{
                        _volumeFrac = fraction_i;
                        lleft = sizey;
                        lright = edge[1] - lowy;
                        llow = rightx - edge[2];
                        lup = sizex;

                    }
                }
                    //5. intersect right1 and up3 planes
                else if(internum[0] == 1 && internum[1] == 3){
                    fraction_i = 0.5 * (upy - edge[1]) * (rightx - edge[3]) / totalVol;

                    if(distsq < r * r) {
                        _volumeFrac = 1 - fraction_i;
                        lleft = sizey;
                        lright = edge[1] - lowy;
                        llow = sizex;
                        lup = edge[3] - leftx;
                    }
                    else{
                        _volumeFrac = fraction_i;
                        lleft = sizey;
                        lright = upy - edge[1];
                        llow = sizex;
                        lup = rightx - edge[3];

                    }
                }
                    //6. intersect lower2 and up3 planes (internum[0] == 2 && internum[1] == 3)
                else {
                    fraction_i = 0.5 * (edge[2] - leftx + edge[3] - leftx) * sizey / totalVol;

                    if(x > r) {
                        _volumeFrac = fraction_i;
                        lleft = sizey;
                        lright = sizey;
                        llow = edge[2] - leftx;
                        lup = edge[3] - leftx;
                    }
                    else{
                        _volumeFrac = 1 - fraction_i;
                        lleft = sizey;
                        lright = sizey;
                        llow = rightx - edge[2];
                        lup = rightx - edge[3];
                    }
                }
            }


            _partialArea = {{lleft * sizez, lright * sizez, llow * sizez, lup * sizez, _volumeFrac * sizex * sizey, _volumeFrac * sizex * sizey}};
        }
        break;
    }
}





bool operator==(const Compartment& a, const Compartment& b) {
    if(a.numberOfSpecies()!=b.numberOfSpecies() ||
       a.numberOfInternalReactions()!=b.numberOfInternalReactions())
        return false;

    if(typeid(a)!=typeid(b))
        return false;

    bool spec_bool = false;
    auto sit_pair = mismatch(a._species.species().begin(),
                             a._species.species().end(),
                             b._species.species().begin(),
            [](const unique_ptr<Species> &A, const unique_ptr<Species> &B)
            {return (*A)==(*B); });
    if(sit_pair.first==a._species.species().end())
        spec_bool=true;


    bool reac_bool = false;
    auto rit_pair = mismatch(a._internal_reactions.reactions().begin(),
                             a._internal_reactions.reactions().end(),
                             b._internal_reactions.reactions().begin(),
            [](const unique_ptr<ReactionBase> &A, const unique_ptr<ReactionBase> &B)
            {return (*A)==(*B);});
    if(rit_pair.first==a._internal_reactions.reactions().end())
        reac_bool=true;

    return spec_bool && reac_bool;
}

} // namespace medyan
