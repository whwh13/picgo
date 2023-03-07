
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

#ifndef MEDYAN_MController_h
#define MEDYAN_MController_h

#include <memory> // unique_ptr
#include <optional> // optional
#include <type_traits>

#include "common.h"
#include "Mechanics/ForceField/ForceFieldDiagnose.hpp"
#include "Mechanics/ForceField/ForceFieldManager.h"
#include "Mechanics/Minimizer/CGMethod.hpp"
#include "Rand.h"
#include "Structure/DofSerializer.hpp"
#include "Structure/SubSystem.h"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {

template< typename InterruptFunc >
struct CGInterruptSettings {
    int forceCallLimitPerInterrupt = 0;
    int interruptCallLimit = 0; // 0 means no interruption, in which case force call limit is ignored.
    InterruptFunc func {};
};
template< typename InterruptFunc >
CGInterruptSettings(int, int, InterruptFunc) -> CGInterruptSettings<InterruptFunc>;

// Default interruption.
constexpr CGInterruptSettings defaultCGInterruptSettings { 0, 0, [] {} };
using DefaultCGInterruptFunc = decltype(defaultCGInterruptSettings.func);

/// Used to initialize, control, and run the mechanical components of a simulation

/*!
 *  MechanicalController is a class used by the SubSystem to initialize [ForceFields]
 *  (@ref ForceField), given an initial selection of which force fields should be 
 *  included. It can compute forces and energies over all force fields using the
 *  ForceFieldManager, as well as run energy [Minimization] (@ref Minimizer) algorithms.
 */
class MController {
    
public:
    ForceFieldManager ffm;  ///< Container and methods for all force
                                   ///< fields in system

    // Parameters used in conjugate gradient minimization.
    std::optional<ConjugateGradientParams> cgParams;
    medyan::CGMethod cgMinimizer; // Conjugate gradient minimizer.

private:
    SubSystem* _subSystem; ///< A pointer to the subsystem

    /// Initialize the force-fields used in the simualtion
    void initializeFF (const SimulConfig&);
    
    /// Initialize the minimization algorithms used in the simulation
    void initializeMinAlgorithm (const MechParams::MechanicsAlgorithm& minimizer);

public:
    /// Constructor which sets a subsystem pointer
    MController(SubSystem* s) : ffm(s) {
        _subSystem = s;

        if(_subSystem->getCylinderLoadForceFunc()) {
            log::warn("The cylinder load force function has already been set.");
        }
        else {
            _subSystem->setCylinderLoadForceFunc([this](Cylinder* c, ForceFieldTypes::LoadForceEnd end) {
                ffm.computeLoadForce(c, end);
            });
        }
    }
    ~MController() {
        if(_subSystem->getCylinderLoadForceFunc()) {
            _subSystem->setCylinderLoadForceFunc(nullptr);
        }
        else {
            log::warn("The cylinder load force function has already been deleted.");
        }
    }
    
    /// Initialze the force fields and minimizers used
    void initialize(const SimulConfig& conf) {
        initializeFF(conf);
        initializeMinAlgorithm(conf.mechParams.mechanicsAlgorithm);

        // Initialize protein curvature mismatch.
        //------------------------------
        _subSystem->initializeProteinCurvatureMismatch(conf);
    }

    /// Run a minimization on the system using the chosen algorithm
    // Returns energy minimization report.
    template< typename InterruptFunc = DefaultCGInterruptFunc >
    auto run(const SimulConfig& conf, const CGInterruptSettings<InterruptFunc>& cgInterruptSettings = defaultCGInterruptSettings) {

        MinimizationResult res;

        //-------------------------------------------------
        // Energy minimization
        //-------------------------------------------------

        // Before minimization, serialize all the system data and prepare force fields.
        std::vector<floatingpoint> coord;
        std::vector<floatingpoint> force;
        FFCoordinateStartingIndex si {};

        const auto callMinimization = [&, this](const ConjugateGradientParams& cgParams) {
            si = serializeDof(*_subSystem, coord);
            ffm.vectorizeAllForceFields(si, conf);

            // Minimize mechanical energy
            const auto energyFunc = [this](floatingpoint* coord) {
                return ffm.computeEnergy(coord);
            };
            const auto forceFunc = [this](floatingpoint* coord, floatingpoint* force, int numVal) {
                ffm.computeForces(coord, force, numVal);
            };
            const auto energyFuncIndividual = [this](floatingpoint* coord) {
                return ffm.computeIndividualEnergies(coord);
            };
            const auto funcRecoveryOnError = [&, this](std::vector<floatingpoint>& coord, const std::vector<floatingpoint>& searchDir) {
                moveAlongDirectionLimitingBeadMove_(coord, searchDir, cgParams.maxDist, si);
            };
            const auto res = cgMinimizer.minimize(cgParams, coord, si.ndof, energyFunc, forceFunc, energyFuncIndividual, funcRecoveryOnError, &force);

            // Deal with minimization errors.
            if(!res.success()) {
                if(res.errorFlags & MinimizationResult::errorLineSearch) {
                    log::error("Line search error occurred. Running diagnosis...");
                    // Run line search diagnosis.
                    diagnoseForLineSearch(
                        ffm,
                        coord,
                        cgMinimizer.searchDir,
                        si.ndof,
                        cgMinimizer.LAMBDATOL
                    );
                    // Throw error because we cannot proceed.
                    throw std::runtime_error("Error in line search.");
                }
            }

            // Copy the coordinate and force data back to the system
            deserializeDof(*_subSystem, coord, force);

            return res;
        };

        if(cgInterruptSettings.interruptCallLimit == 0) {
            // Do not interrupt. Call minimization normally.
            res = callMinimization(this->cgParams.value());
        }
        else {
            ConjugateGradientParams cgParamsEach = this->cgParams.value();
            cgParamsEach.forceCallLimit = cgInterruptSettings.forceCallLimitPerInterrupt;
            cgParamsEach.reportOnCallLimit = false;

            bool converged = false;

            for(int iter = 0; iter < cgInterruptSettings.interruptCallLimit; ++iter) {
                // Call minimization.
                const auto eachRes = callMinimization(cgParamsEach);

                // Append or overwrite to res.
                res.updateWithMinimizationResult(eachRes);

                if(res.converged()) {
                    converged = true;
                    break;
                } else {
                    // Run interrupt if not converged.
                    cgInterruptSettings.func();
                }
            }

            // Last iteration.
            if(!converged) {
                cgParamsEach.reportOnCallLimit = true;

                // Call minimization.
                const auto eachRes = callMinimization(cgParamsEach);

                // Append or overwrite to res.
                res.updateWithMinimizationResult(eachRes);
            }
        }

        //-------------------------------------------------
        // Computations after energy minimization.
        //-------------------------------------------------
        afterMinimization_(conf, si, coord);


        // Return minimization report
        return res;
    }

    // Will vectorize system and compute energy and force once, but no energy minimization is performed.
    void updateMechanics(const SimulConfig& conf) {
        std::vector<floatingpoint> coord;
        std::vector<floatingpoint> force;
        const auto si = serializeDof(*_subSystem, coord);
        force.assign(coord.size(), 0);
        ffm.vectorizeAllForceFields(si, conf);

        // Compute energy and force once.
        ffm.computeEnergy(coord.data());
        ffm.computeForces(coord.data(), force.data(), force.size());

        // Deserialization.
        deserializeDof(*_subSystem, coord, force);

        // After minimization.
        afterMinimization_(conf, si, coord);
    }

    
    ForceFieldManager* getForceFieldManager(){
        return &ffm;
    }


private:
    // Auxiliary function to do computations after minimization.
    void afterMinimization_(const SimulConfig& conf, const FFCoordinateStartingIndex& si, const std::vector<floatingpoint>& coord) {
        // compute the Hessian matrix at this point if the feature is enabled.
        if(conf.mechParams.hessTracking) {
            if(ffm.hessCounter % conf.mechParams.hessSkip == 0 || tau() > conf.chemParams.chemistryAlgorithm.runTime){
                
                if(tau() > 0.0){
                    ffm.computeProjections(si, coord);
                };
                ffm.computeHessian(coord, si.ndof, conf.mechParams.hessDelta);
            }
            ffm.hessCounter += 1;
        }

        // Compute auxiliary parameters that require valid vectorization.
        ffm.assignallforcemags();

        // Invalidate membrane index cache.
        for(auto& m : _subSystem->membranes) {
            invalidateIndexCacheForFF(m.getMesh());
        }

        //------------------------------
        // From here, the vectorization is no longer valid.
        //------------------------------

        // Update auxiliary parameters, such as load forces.
        // Some parameters might be used by the chemical rate changer.
        ffm.computeAuxParams(*_subSystem);
    }

    // Auxiliary function to attempt to move beads to get system out of stuck.
    // Note:
    //   - Must be used during energy minimization, where a valid vectorization is available.
    void moveAlongDirectionLimitingBeadMove_(
        std::vector<floatingpoint>&       coord,
        const std::vector<floatingpoint>& moveDir,
        floatingpoint maxDistBead,
        const FFCoordinateStartingIndex& si
    ) {
        using namespace medyan;

        // Create a temporary vector for pushing forward move direction.
        auto moveDirFull = moveDir;
        moveDirFull.resize(coord.size());

        // Make sure that target coordinates have up-to-date dependent coordinates.
        ffm.computeDependentCoordinates(coord.data());
        const auto debugCoordPrev = coord;

        // Push forward the move direction to all coordinates.
        ffm.pushForwardIndependentTangentVector(coord.data(), moveDirFull.data());

        // Determine step size.
        auto step = std::numeric_limits<floatingpoint>::infinity();
        for(Bead* pb : Bead::getBeads()) {
            auto moveDirBead = makeRefVec<3>(moveDirFull.data() + medyan::findBeadCoordIndex(*pb, si));
            for(auto& eachMoveDir : moveDirBead) {
                step = std::min(step, maxDistBead / std::abs(eachMoveDir));
            }
        }
        for(auto& m : _subSystem->membranes) {
            const auto& mesh = m.getMesh();
            for(auto& v : mesh.getVertices()) {
                auto moveDirVertex = makeRefVec<3>(moveDirFull.data() + v.attr.cachedCoordIndex);
                for(auto& eachMoveDir : moveDirVertex) {
                    step = std::min(step, maxDistBead / std::abs(eachMoveDir));
                }
            }
        }

        // If step size is still infinity, something is wrong. Maybe the number of beads/vertices is 0.
        if(step == std::numeric_limits<floatingpoint>::infinity()) {
            log::error("During an attempt to recover by moving along search direction, step size is infinity.");
            log::info("num_beads={}, num_vertices={}", Bead::getBeads().size(), _subSystem->vertices.size());
            throw std::runtime_error("Invalid step size.");
        }

        // Move the system according to the move direction.
        log::debug("Try to move system with step size {}", step);
        for(int i = 0; i < si.ndof; ++i) {
            coord[i] += step * moveDir[i];
        }

        // Update dependent coordinates again.
        ffm.computeDependentCoordinates(coord.data());

        // Compare move result.
        {
            floatingpoint maxActualMove = 0;
            int maxi = -1;
            for (int i = 0; i < coord.size(); ++i) {
                auto diff = std::abs(debugCoordPrev[i] - coord[i]);
                if (diff > maxActualMove) {
                    maxActualMove = diff;
                    maxi = i;
                }
            }
            log::debug("Max move magnitude {} at index {}", maxActualMove, maxi);
        }
    }

};

} // namespace medyan

#endif
