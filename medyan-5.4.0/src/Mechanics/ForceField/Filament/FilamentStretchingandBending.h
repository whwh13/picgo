
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

#ifndef MEDYAN_FilamentStretchingandBending_h
#define MEDYAN_FilamentStretchingandBending_h

#include "common.h"
#include "Mechanics/ForceField/ForceField.h"

namespace medyan {
//FORWARD DECLARATIONS
class Filament;

/// Represents a Filament bending interaction
template <class FStretchingandBendingInteractionType>
class FilamentStretchingandBending : public ForceField {

private:
	FStretchingandBendingInteractionType _FFType;
	const int numthreads = 1;

	// Cache of vectorized data
	Size _numhybridInteractions;
	Size _strnumInteractions;
	std::vector<int> beadSet;
	///Array describing the constants in calculation
	std::vector<FP> kbend;
	std::vector<FP> eqt;

	std::vector<int> beadSetcylsansbending;
	std::vector<FP> kstrsansbending;
	std::vector<FP> eqlsansbending;
	std::vector<FP> kstr;
	std::vector<FP> eql;
	std::vector<FP> totalenergy;//Will have 3 entries. The total of stretching, bending
	// and sum of stretching + bending
	std::vector<int> cylSet;//Informs which element in cyllength set to consider.
	std::vector<int> cylSetcylsansbending;//Informs which element in cyllength set to consider.

	std::vector<int> beadSetall;//element i helps determine position in coord array to look at.
	std::vector<char> beadpaircyllengthstatus;//element i informs if beads i and i+1 form a bond.
	std::vector<char> beadtriplet_hingestatus;//element i informs if bead i, i+1, and i+2 form a hinge.
	std::vector<FP> cyllengthset;
	std::vector<FP> cylbenddotproduct;

	//precomputes necessary cylinder length and dot products
	void precomputevars(FP *coord);

public:

	///Array describing indexed set of interactions
	///For filaments, this is a 3-bead potential
	const static int n = 3;
	const static int nstr = 2;
	const static int ncylperint = 2;

	virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;

	virtual FP computeEnergy(FP *coord) override;
	virtual void computeForces(FP *coord, FP *f) override;

	virtual std::string getName() override {return "FilamentStretchingAndBending";}
};

} // namespace medyan

#endif
