
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
#include "FilamentStretchingandBending.h"

#include "FilamentStretchingHarmonicandBendingHarmonic.h"
#include "FilamentStretchingHarmonicandBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "SysParams.h"
#include "SubSystem.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "cross_check.h"

namespace medyan {
using namespace mathfunc;

template <class FBendingInteractionType>
void FilamentStretchingandBending<FBendingInteractionType>::precomputevars(FP *coord) {

    int cylcount = 0;
    int hingecount = 0;
    floatingpoint *coord1, *coord2, *coord3;
    bool fetchcoordsstatus = true;
    if(Bead::numElements() >= 1) {
	    for (auto b = 0; b < Bead::numElements()-1; b++) {
		    if (fetchcoordsstatus) {
				coord1 = &coord[beadSetall[b]];
				coord2 = &coord[beadSetall[b + 1]];
		    } else {
			    coord1 = coord2;
			    coord2 = coord3;
		    }
		    if (beadpaircyllengthstatus[b]) {
			    cyllengthset[cylcount] = twoPointDistance(coord1, coord2);
			    cylcount++;
		    }
		    fetchcoordsstatus = true;
		    if (beadtriplet_hingestatus[b]) {
			    coord3 = &coord[beadSetall[b + 2]];
			    cylbenddotproduct[hingecount] = scalarProduct(coord1, coord2, coord2,
			                                                  coord3);
			    hingecount++;
			    //If there are connections between b+1, b+2, and b+3.
			    if (beadtriplet_hingestatus[b + 1]) {
				    fetchcoordsstatus = false;
			    }
		    }
	    }
    }
}

template <class FBendingInteractionType>
void FilamentStretchingandBending<FBendingInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {

	// Count number of interactions that involve both a stretching and bending calculation
	_numhybridInteractions = 0;
	for(auto f : Filament::getFilaments())
		if(f->getCylinderVector().size() > 1) _numhybridInteractions += f->getCylinderVector().size() - 1;

	totalenergy.assign(3*numthreads, 0);
	for(int i = 0; i < 3*numthreads; i++)
	    totalenergy[i] = (floatingpoint)0.0;
	// The first cylinder in each filament is not considered in the hybrid stretching
	// bending paradigm. Need a function to calculate energies separately for those
	// cylinders.
	_strnumInteractions = Filament::getFilaments().size();
    //These two numbers should match
/*	cout<<"Str NumInt Method  1 "<<Cylinder::getCylinders().size()<<" Method 2 "
		<<_strnumInteractions+_numhybridInteractions<<endl;*/
	beadSetcylsansbending.assign(nstr * _strnumInteractions, 0);
	kstrsansbending.assign(_strnumInteractions, 0);
	eqlsansbending.assign(_strnumInteractions, 0);

	kstr.assign(_numhybridInteractions, 0);
	eql.assign(_numhybridInteractions, 0);
	beadSet.assign(n * _numhybridInteractions, 0);
	kbend.assign(_numhybridInteractions, 0);
	eqt.assign(_numhybridInteractions, 0);

	//Testing Cylset stores IDs
	cylSet.assign(2*_numhybridInteractions, 0);
	cylSetcylsansbending.assign(_strnumInteractions, 0);

	//Datasets to help calculate precomputed vars.
    auto Totalcyl = Cylinder::getCylinders().size();
    auto Totalbeads = Bead::getBeads().size();
    beadSetall.assign(Totalbeads, 0);
    beadpaircyllengthstatus.assign(Totalbeads, 0);
    beadtriplet_hingestatus.assign(Totalbeads, 0);
	cyllengthset.assign(Totalcyl, 0);
    cylbenddotproduct.assign(_numhybridInteractions, 0);
    int beadcount = 0;
    int hingecount = 0;

	int i = 0;

	int istr = 0;

	int cylcount = 0;

	for (auto f: Filament::getFilaments()) {

		auto cyl = *f->getCylinderVector().begin();
		beadSetcylsansbending[nstr * istr] = cyl->getFirstBead()->getIndex() * 3 + si.bead;
		beadSetcylsansbending[nstr * istr + 1] = cyl->getSecondBead()->getIndex() * 3 + si.bead;
		kstrsansbending[istr] = cyl->getMCylinder()->getStretchingConst();
		eqlsansbending[istr] = cyl->getMCylinder()->getEqLength();
		cylSetcylsansbending[istr] = cylcount;
		beadSetall[beadcount] = beadSetcylsansbending[nstr * istr];
//        cout<<beadcount<<endl;
		beadcount++;
        beadpaircyllengthstatus[beadcount-1] = 1;
//        cout<<beadcount-1<<endl;
		beadSetall[beadcount] = beadSetcylsansbending[nstr * istr + 1];
		beadcount++;
		cylcount++;
		istr++;

		if (f->getCylinderVector().size() > 1){

			for (auto it = f->getCylinderVector().begin()+1;
			     it != f->getCylinderVector().end(); it++){

				auto it2 = it - 1;
				beadSet[n * i] = (*it2)->getFirstBead()->getIndex() * 3 + si.bead;
				beadSet[n * i + 1] = (*it)->getFirstBead()->getIndex() * 3 + si.bead;
				beadSet[n * i + 2] = (*it)->getSecondBead()->getIndex() * 3 + si.bead;

                cylSet[ncylperint * i] = cylcount - 1;
                cylSet[ncylperint * i + 1] = cylcount;
				cylcount++;

				beadSetall[beadcount] = beadSet[n * i + 2];
				beadpaircyllengthstatus[beadcount-1] = 1;
//                cout<<beadcount-1<<endl;
				beadcount++;

                beadtriplet_hingestatus[hingecount] = 1;
//                cout<<hingecount<<endl;
                hingecount++;

				kbend[i] = (*it)->getMCylinder()->getBendingConst();
				eqt[i]  = (*it)->getMCylinder()->getEqTheta();

				kstr[i] = (*it)->getMCylinder()->getStretchingConst();
				eql[i]  = (*it)->getMCylinder()->getEqLength();

				i++;
			}
		}
		beadpaircyllengthstatus[beadcount-1] = 0;
        beadtriplet_hingestatus[hingecount] = 0;
        beadtriplet_hingestatus[hingecount+1] = 0;
        hingecount = hingecount + 2;
//        cout<<hingecount<<endl;
	}

}


//Needs to have a return value.
template <class FStretchingandBendingInteractionType>
FP FilamentStretchingandBending<FStretchingandBendingInteractionType>::
        computeEnergy(FP *coord){


	const int startID = 0;
	int threadID = 0;
	precomputevars(coord);
    _FFType.energy(coord, _numhybridInteractions, cylSet.data(), cyllengthset.data(), cylbenddotproduct.data(),
            kstr.data(), kbend.data(), eql.data(), eqt.data(), totalenergy.data(), startID, _numhybridInteractions, threadID);
    _FFType.energy(coord, cylSetcylsansbending.data(), cyllengthset.data(), kstrsansbending.data(),
           eqlsansbending.data(), totalenergy.data(), startID, _strnumInteractions, threadID);

	floatingpoint U = (floatingpoint) 0.0;
	for(int t = 0 ; t < numthreads; t++){
		for(int j = 0; j <3;j++) {
			if(totalenergy[3*t+j] >= (floatingpoint) 0.0)
				U +=  totalenergy[3*t+j];
			else
				return (floatingpoint) -1.0;
		}
	}

	return U;

}

template <class FStretchingandBendingInteractionType>
void FilamentStretchingandBending<FStretchingandBendingInteractionType>::computeForces
(FP *coord, FP *f) {
	 precomputevars(coord);
     const int startID = 0;
     int threadID = 0;
    _FFType.forces(coord, f, _numhybridInteractions, beadSet.data(), cylSet.data(), cyllengthset.data(), cylbenddotproduct.data(),
         kstr.data(), kbend.data(), eql.data(), eqt.data());
    _FFType.forces(coord, f,  beadSet.data(), cylSetcylsansbending.data(), cyllengthset.data(), beadSetcylsansbending.data(),
         kstrsansbending.data(), eqlsansbending.data(), startID, _strnumInteractions, threadID);

#ifdef DETAILEDOUTPUT
	floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
//        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
//                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<getName()<<" "<<maxF<<endl;
#endif

}

// Explicit instantiation.
template class FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>;
template class FilamentStretchingandBending<FilamentStretchingHarmonicandBendingCosine>;

} // namespace medyan
