
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

#ifndef MEDYAN_ForceFieldManager_h
#define MEDYAN_ForceFieldManager_h

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <unordered_map>

#include "common.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Types.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"


namespace medyan {
typedef Eigen::Triplet<double> Triplet;

// Forward declarations
class Cylinder;


/// A class to store and iterate over all [ForceFields](@ref ForceField).
/*!
 *  The ForceFieldManager is used to store all [ForceFields](@ref ForceField)
 *  initialized by the system, as well as iterate over these potentials and calculate
 *  total forces and energies. This class contains functions for the said calculations.
 */
class ForceFieldManager {

public:
    std::vector<std::unique_ptr<ForceField>> forceFields; ///< All forcefields in the system
    SubSystem* ps = nullptr;

    // Requirement for geometry computation.
    SurfaceGeometrySettings surfaceGeometrySettings {};


    ForceFieldManager(SubSystem* ps) : ps(ps) {}

    // Accessors to force fields.
    auto& getForceFields() const { return forceFields; }

    /// Vectorize all interactions involved in calculation
    void vectorizeAllForceFields(const FFCoordinateStartingIndex&, const SimulConfig&);

    // Auxiliary function to update dependent coordinates.
    void computeDependentCoordinates(floatingpoint* coord) const {
        for(auto& pff : forceFields) {
            pff->computeDependentCoordinates(coord);
        }
    }
    // Auxiliary function to propagate forces on dependent coordinates to independent coordinates.
    void propagateDependentForces(const floatingpoint* coord, floatingpoint* force) const {
        for(auto& pff : forceFields) {
            pff->propagateDependentForces(coord, force);
        }
    }
    // Auxiliary function to push forward an independent tangent vector to all coordinate space.
    void pushForwardIndependentTangentVector(const floatingpoint* coord, floatingpoint* vec) const {
        for(auto& pff : forceFields) {
            pff->pushForwardIndependentTangentVector(coord, vec);
        }
    }
    // (Temporary) auxiliary function before calling energy computation for force fields.
    // When membrane geometry use dependent coordinates, this function can be removed.
    void beforeComputeEnergy(floatingpoint* coord) const {
        computeDependentCoordinates(coord);
        for(auto& m : ps->membranes) {
            updateGeometryValue(m.getMesh(), coord, surfaceGeometrySettings);
        }
    }
    // (Temporary) auxiliary function before calling force computation for force fields.
    // When membrane geometry use dependent coordinates, this function can be removed.
    void beforeComputeForce(floatingpoint* coord) const {
        computeDependentCoordinates(coord);
        for(auto& m : ps->membranes) {
            updateGeometryValueWithDerivative(m.getMesh(), coord, surfaceGeometrySettings);
        }
    }

    /// Compute the energy using all available force fields
    /// @return Returns infinity if there was a problem with a ForceField
    /// energy calculation, such that beads will not be moved to this
    /// problematic configuration.
    floatingpoint computeEnergy(floatingpoint *coord, bool verbose = false) const;

    // Get all individual energy names. Must correspond to the result of computeIndividualEnergies function.
    std::vector<std::string> getIndividualEnergyNames() const {
        std::vector<std::string> names;
        for(auto& pff : forceFields) {
            names.push_back(pff->getName());
        }
        return names;
    }
    // Compute energies from each force field.
    EnergyReport computeIndividualEnergies(floatingpoint* coord) const;
    // Compute energies from each force field, but the energies have unit of kT.
    EnergyReport computeEnergyHRMD(floatingpoint *coord) const;
    
    
    /// Compute the forces of all force fields 
    // numVar is the number of all variables, including independent and dependent variables.
    void computeForces(floatingpoint *coord, floatingpoint* force, int numVar);
    
    // Compute the Hessian matrix if the feature is enabled.
    // Requires valid vectorization.
    void computeHessian(const std::vector<floatingpoint>& coord, int total_DOF, float delta);
    
    void setCurrBeadMap(const FFCoordinateStartingIndex& si);
    
    // compute the displacement projections along the eigenvectors.
    // Warning: this function only works if all bead coordinates are independent coordinates. Otherwise, out-of-bound access may result.
    void computeProjections(const FFCoordinateStartingIndex&, const std::vector<floatingpoint>& currCoords);
    
    void clearHessian(int a){
        if(a == 0){
            hessianVector.clear();
        }else if(a==1){
            evaluesVector.clear();
            IPRIVector.clear();
            IPRIIVector.clear();
            
        }else{
            projectionsVector.clear();
            tauVector.clear();
        };
    }
    
    vector<floatingpoint> HRMDenergies;
    
    vector<vector<vector<floatingpoint>>> hessianVector;
    
    vector<Eigen::VectorXcd> evaluesVector;
    vector<Eigen::VectorXcd> IPRIVector;
    vector<Eigen::VectorXcd> IPRIIVector;
    Eigen::VectorXcd evalues;
    Eigen::MatrixXcd evectors;
    vector<floatingpoint> tauVector;
    vector<Eigen::VectorXcd> projectionsVector;
    
    int hessCounter;

    // Map bead pointer to coordinate index in the vectorized data.
    std::unordered_map<Bead*, int> prevBeadMap;
    std::unordered_map<Bead*, int> currBeadMap;
    // Previous coordinates during last eigenvector projection.
    std::vector< floatingpoint > prevCoords;



#ifdef CUDAACCL
        cudaStream_t  streamF = NULL;
    /// CUDA Copy forces from f to fprev
    void CUDAcopyForces(cudaStream_t  stream, floatingpoint *f, floatingpoint *fprev);
#endif

    /// Compute the load forces on the beads. This does not update the force (xyz) vector
    /// contained by Bead, but updates the loadForce vector which contains precalculated
    /// load values based on the bead's directionality of growth in a filament.
    void computeLoadForces() const;

    // Compute all auxiliary force-related parameters for the system.
    // For example, some forces on certain elements can be used for changing chemical rates.
    void computeAuxParams(SubSystem& sys) const {
        computeLoadForces();

        for(auto& pff : forceFields) {
            pff->computeAuxParams(sys);
        }
    }

    // Compute the load forces on the bead for a specific cylinder.
    void computeLoadForce(Cylinder* c, ForceField::LoadForceEnd end) const;
#ifdef CROSSCHECK
    /// Reset the forces of all objects
    void resetForces();
#endif
#ifdef CUDAACCL
    vector<int> blocksnthreads;
    int *gpu_nint;
    //@{
    vector<int> bntaddvec2;
    int *gpu_params;
    vector<int> params;
    //@}

#endif
    // Compute all auxiliary force-related parameters for the system.
    // Notes:
    // - Requires valid vectorization.
    // - Can only be called once per minimization step.
    void assignallforcemags();

private:
    chrono::high_resolution_clock::time_point tbegin, tend;
};

} // namespace medyan

#endif
