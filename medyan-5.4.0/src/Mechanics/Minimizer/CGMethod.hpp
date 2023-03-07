
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

#ifndef MEDYAN_Mechanics_Minimizer_CGMethod_hpp
#define MEDYAN_Mechanics_Minimizer_CGMethod_hpp

#include <algorithm>
#include <cmath>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <utility> // forward, move
#include <vector>

#include "common.h"
#include "CUDAcommon.h"
#include "Mechanics/ForceField/ForceFieldManager.h"
#include "Mechanics/Minimizer/LineSearch.hpp"
#include "Mechanics/Minimizer/MinimizationTypes.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// For performing a conjugate gradient minimization method


enum class ConjugateGradientDescentSearch {
    steepest,
    fletcherRieves,
    polakRibiere,
};

struct ConjugateGradientParams {
    // Descent direction search method
    ConjugateGradientDescentSearch descentSearchMethod = ConjugateGradientDescentSearch::polakRibiere;

    floatingpoint maxDist;   ///< Max distance used to move
    floatingpoint lambdaMax; ///< Maximum lambda that can be returned
    floatingpoint lambdaRunningAverageProbability = 0.0; //running average probability.
    std::string lineSearchAlgorithm = "BACKTRACKING";

    // Convergence criteria.
    floatingpoint gradTol      = 1e-1; // Gradient tolerance infinity norm.
    floatingpoint energyRelTol = 1e-9; // Energy change relative tolerance.

    int           energyCallLimit = 0;     // 0 means no limit.
    int           forceCallLimit  = 10000; // 0 means no limit.
    bool          reportOnCallLimit = true;

    // If true, safe mode kicks in to try to restore line search errors.
    bool tryToRecoverInLineSearchError = true;
};


/*!
 *  CGMethod is an abstract class that contains methods for conjugate gradient
 *  minimization. It has functions for various line search techniques, including golden
 *  section, backtracking, and quadratic line search. This base class also contains
 *  parameters for tolerance criterion and other helpful parameters.
 */
class CGMethod {

public:


    ///< Data vectors for calculation
    //-------------------------------------------------------------------------
    // The vectorized data is separated into independent variable section and dependent variable section. The independent variables occupy indices { 0, 1, ..., numDof-1 }, while the dependent variables occupy the rest of the indices.
    // Before each energy/force computation, the dependent coordinates should be updated from independent variables.
    // After each force computation, the forces accumulated on dependent variables should be passed onto the independent variables using the chain rule.

    // The number of independent variables (ie degree of freedom). It is not necessarily the size of the vectors, because dependent variables may be stored as well.
    int                          numDof = 0;
    std::vector< floatingpoint > coordLineSearch; // Temporary coords used during line search
    std::vector< floatingpoint > coordMinE; // Temporary coords used to record position for minimumE

    std::vector< floatingpoint > force; // Negative gradient
    std::vector< floatingpoint > forcePrev; // Previous force
    std::vector< floatingpoint > searchDir; // The current search direction in conjugate gradient method
    std::vector< floatingpoint > forceTol; // The force tolerance (in each dimension); must be positive

    chrono::high_resolution_clock::time_point tbegin, tend;


    //@{
    /// Parameter used in backtracking line search
    const floatingpoint LAMBDAREDUCE = 0.5;     ///< Lambda reduction parameter for backtracking
    floatingpoint LAMBDATOL = 1e-8;       ///< Lambda tolerance parameter

    const floatingpoint SAFELAMBDAREDUCE = 0.9;  ///< Lambda reduction parameter for conservative backtracking

    floatingpoint BACKTRACKSLOPE = 0.4;   ///< Backtracking slope

    const floatingpoint QUADTOL = 0.1; // Used in Quadratic line search
    const floatingpoint LAMBDAQUADTOL = 1e-3; //if values change less than 0.1% between
    // successive runs, consider it converged.

    //additional parameters to help store additional parameters
    floatingpoint minimumE = (floatingpoint) 1e10;
    floatingpoint TotalEnergy = (floatingpoint)0.0;
    floatingpoint maxForcebackup = (floatingpoint)0.0;
    //@}


    //@{

#ifdef CUDAACCL
    int *gpu_mutexlock;
    vector<int> blocksnthreads;
    vector<int> bntaddvector;
    vector<int> bnt;
    int *gpu_nint;
    floatingpoint *gpu_g, *gpu_maxF;
    floatingpoint *gSum;
    floatingpoint *gSum2;
    floatingpoint *gpu_fmax;
    floatingpoint *g_currentenergy;
    floatingpoint *gpu_params = NULL;
    floatingpoint *gpu_FDotF;//curGrad
    floatingpoint *gpu_FADotFA;//newGrad
    floatingpoint *gpu_FADotFAP;//prevGrad
    floatingpoint *gpu_FDotFA;
    floatingpoint *gpu_initlambdalocal;
    bool *gpu_convergencecheck;
    bool *convergencecheck;
    floatingpoint gpuFDotF(floatingpoint *f1, floatingpoint *f2);
    void CUDAresetlambda(cudaStream_t stream);
    void CUDAinitializeLambda(cudaStream_t stream1, bool *check_in, bool *check_out, bool
            *Polaksafestate, int *gpu_state);
    void CUDAfindLambda(cudaStream_t  stream1, cudaStream_t stream2, cudaEvent_t event, bool *checkin, bool
            *checkout, bool *gpu_safestate, int *gpu_state);
    void CUDAprepforbacktracking(cudaStream_t stream, bool *check_in, bool *check_out);
    void CUDAprepforsafebacktracking(cudaStream_t stream, bool *check_in, bool *check_out);
    void CUDAallFDotF(cudaStream_t stream);
    void CUDAallFADotFA(cudaStream_t stream);
    void CUDAallFADotFAP(cudaStream_t stream);
    void CUDAallFDotFA(cudaStream_t stream);
    void CUDAshiftGradient(cudaStream_t stream, bool *Mcheckin);
    void CUDAshiftGradientifSafe(cudaStream_t stream, bool *Mcheckin, bool *Scheckin);
//    void CUDAgetPolakvars(bool calc_safestate,cudaStream_t streamcalc, floatingpoint* gpu_GRADTOL, bool *gminstatein,
//                                    bool *gminstateout, bool *gsafestateout, volatile bool *cminstate);
    void CUDAgetPolakvars(cudaStream_t streamcalc, floatingpoint* gpu_GRADTOL, bool *gminstatein,
    bool *gminstateout, volatile bool *cminstate);
    void CUDAgetPolakvars2(cudaStream_t streamcalc, bool *gsafestateout);
    void CUDAinitializePolak(cudaStream_t stream, bool *minstatein, bool *minstateout, bool *safestatein, bool
    *safestateout);
    void CUDAmoveBeads(cudaStream_t stream, bool *gpu_checkin );
//    void getmaxFCUDA(floatingpoint *gpu_forceAux, int *gpu_nint, floatingpoint *gpu_fmax);
    //PING PONG for backtracking (both normal and safe)
//    struct backtrackingvars {
//        floatingpoint currentEnergy;
//        floatingpoint energyLambda;
//        floatingpoint lambda;
//    };
    bool *g_stop1, *g_stop2, *g_s1, *g_s2, *g_ss;
//    backtrackingvars *bvar, *gpu_bvar1, *gpu_bvar2, *g_b1, *g_b2, *g_bs;
    cudaStream_t s1 = NULL, s2 = NULL, s3 = NULL, *sp1, *sp2, *sps, stream_bt = NULL;
    cudaStream_t stream_startmin = NULL;
    cudaEvent_t e1 = NULL, e2 = NULL, *ep1, *ep2, *eps;
    cudaEvent_t  e = NULL;
    int *gpu_state;
    // @PING PONG ENDS

#endif
    volatile bool *cconvergencecheck;
    bool *h_stop;

    void calcCoordLineSearch(const std::vector<floatingpoint>& coord, floatingpoint d);

    double searchDirDotSearchDir() const;
    double forceDotForce() const;
    double forceDotForcePrev() const;
    double searchDirDotForce() const;

    // Check whether the force is no larger than tolerance in every dimension
    bool forceBelowTolerance() const {
        for(std::size_t i = 0; i < numDof; ++i) {
            if(std::abs(force[i]) > forceTol[i]) return false;
        }
        return true;
    }
    // Check whether the force is no larger than a relatex tolerance in every dimension
    bool forceBelowRelaxedTolerance(floatingpoint factor) const {
        for(std::size_t i = 0; i < numDof; ++i) {
            if(std::abs(force[i]) > factor * forceTol[i]) return false;
        }
        return true;
    }

    /// Get the max force in the system
    floatingpoint maxF() const {
        floatingpoint magMax = 0;
        for(Index i = 0; i < numDof; ++i) {
            magMax = std::max(magMax, std::abs(force[i]));
        }
        if(!std::isfinite(magMax)) {
            log::error("maxF is infinity. Check parameters. Exiting.");
            throw std::runtime_error("maxF error");
        }
        return magMax;
    }

    /// Transfers data to lightweight arrays for min
    void startMinimization();
    /// Transfers updated coordinates and force to bead members
    void endMinimization();

    /// Move beads in search direction by d
    void moveAlongSearchDir(std::vector<floatingpoint>& coord, floatingpoint d);

    /// shift the gradient by d
    void shiftSearchDir(double d);


#ifdef CUDAACCL
    //@{
    floatingpoint backtrackingLineSearchCUDA(ForceFieldManager& FFM, floatingpoint MAXDIST,
                                  floatingpoint LAMBDAMAX, bool *gpu_safestate);
    //@}
#endif // CUDAACCL

    //@{
    /// Linear search methods
    /// A simple backtracking search method that computes an optimal energy change and compares the backtracked energy to it.
    // Note:
    //   - EnergyFunc is a function of signature (floatingpoint* coord) -> floatingpoint.
    template< typename EnergyFunc >
    LineSearchResult backtrackingLineSearch(
        std::vector<floatingpoint>&          coord,
        EnergyFunc&&                         energyFunc,     // Energy computation function.
        floatingpoint                        MAXDIST,
        floatingpoint                        maxForce,
        floatingpoint                        LAMBDAMAX,
        medyan::LambdaRunningAverageManager* ptrLambdaAvg,   // Lambda running average manager. Can be nullptr if running average is not used.
        floatingpoint                        armijoRuleFactor // Within [0, 1). A factor of zero means any configuration with a lower energy is sufficient for this line search.
    ) {
        using namespace std;

        LineSearchResult res;

        //return zero if no forces
        if(maxForce == 0.0) {
            // Zero lambda, zero energy change.
            return res;
        }

        //calculate first lambda
        floatingpoint lambda = min(LAMBDAMAX, MAXDIST / maxForce);
        if(ptrLambdaAvg && ptrLambdaAvg->runningaveragestatus) {
            lambda = min(lambda, ptrLambdaAvg->suggestLambdaMax());
        }

        //@} Lambda phase 1
#ifdef DETAILEDOUTPUT_LAMBDA
        std::cout<<"SL lambdamax "<<LAMBDAMAX<<" serial_lambda "<<lambda<<" fmax "<<maxForce<<" state "<<sconvergencecheck<<endl;
#endif
        auto tbegin = chrono::high_resolution_clock::now();
        const floatingpoint currentEnergy = energyFunc(coord.data());
        ++res.numEnergyCall;
        CUDAcommon::tmin.computeenerycallszero++;
        auto tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
        CUDAcommon::tmin.computeenergy+= elapsed_energy.count();
        CUDAcommon::tmin.computeenergyzero+= elapsed_energy.count();

        #ifdef TRACKDIDNOTMINIMIZE
        SysParams::Mininimization().TotalE.push_back(currentEnergy);
        #endif

#ifdef DETAILEDOUTPUT_ENERGY
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        floatingpoint cuda_energy[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(floatingpoint),
                                        cudaMemcpyDeviceToHost));
        std::cout<<"Total Energy CE pN.nm CUDA "<<cuda_energy[0]<<" SERL "<<currentEnergy<<endl;
        std::cout<<endl;
#endif

        const auto expectedSlope = searchDirDotForce();
        const auto reducedSlope = armijoRuleFactor * expectedSlope;

        bool sconvergencecheck = false;
        floatingpoint energyChange = (floatingpoint)0.0;
        floatingpoint energyLambda = (floatingpoint)0.0;
        int iter = 0;
        while(!(sconvergencecheck)) {
            iter++;
            //TODO let each forcefield calculate energy IFF conv state = false. That will help
            // them avoid unnecessary iterations.
            //let each forcefield also add energies to two different energy variables.

            auto tbegin = chrono::high_resolution_clock::now();
            calcCoordLineSearch(coord, lambda);
            energyLambda = energyFunc(coordLineSearch.data());
            ++res.numEnergyCall;
            CUDAcommon::tmin.computeenerycallsnonzero++;
            auto tend = chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
            CUDAcommon::tmin.computeenergy+= elapsed_energy.count();
            CUDAcommon::tmin.computeenergynonzero+= elapsed_energy.count();

#ifdef DETAILEDOUTPUT_ENERGY
            CUDAcommon::handleerror(cudaDeviceSynchronize());
            floatingpoint cuda_energy[1];
            CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(floatingpoint),
                                            cudaMemcpyDeviceToHost));
            std::cout<<"Total Energy EL pN.nm CUDA "<<cuda_energy[0]<<" SERL "
                    ""<<energyLambda<<endl;
            std::cout<<endl;
#endif

            //@{ Lambda phase 2
            const floatingpoint idealEnergyChange = -reducedSlope * lambda;
            if(lambda<0 || idealEnergyChange>0) {
                sconvergencecheck = true;
                lambda=0.0;
            }
            energyChange = energyLambda - currentEnergy;
#ifdef DETAILEDOUTPUT_LAMBDA
            std::cout<<"BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA"
                    " "<<searchDirDotForce()<<endl;
            std::cout<<"SL energyChange "<<energyChange<<" idealEnergyChange "
                    ""<<idealEnergyChange<<endl;
#endif
            //return if ok
            //Armijo conditon
            if(energyChange <= idealEnergyChange) {
                sconvergencecheck = true;
            }
            else {
                //reduce lambda
                lambda *= LAMBDAREDUCE;
            }

            if(lambda <= 0.0 || lambda <= LAMBDATOL) {
                sconvergencecheck = true;

                // Set lambda to zero, indicating line search failure.
                lambda = 0.0;
                energyLambda = currentEnergy;
                energyChange = 0;
            }
#ifdef DETAILEDOUTPUT_LAMBDA
            std::cout<<"SL2 BACKTRACKSLOPE "<<BACKTRACKSLOPE<<" lambda "<<lambda<<" allFDotFA "
                                                                                <<searchDirDotForce()<<endl;
            std::cout<<"SL2 energyChange "<<energyChange<<" idealEnergyChange "
                    ""<<idealEnergyChange
                    <<" lambda "<<lambda<<" state "<<sconvergencecheck<<endl;
#endif
        }

        //running average
        if(ptrLambdaAvg) {
            ptrLambdaAvg->record(lambda);
        }

        TotalEnergy = energyLambda;

        res.lambda = lambda;
        res.newEnergy = energyLambda;
        return res;
    }
    


    ///Quadratic line search introduced from LAMMPS based on Dennis and Schnabel
    // Note:
    //   - EnergyFunc is a function of signature (floatingpoint* coord) -> floatingpoint.
    //   - ForceFunc is a function of signature (floatingpoint* coord, floatingpoint* force, int numVar) -> void.
    template< typename EnergyFunc, typename ForceFunc >
    medyan::LineSearchResult quadraticLineSearch(
        std::vector<floatingpoint>& coord,
        EnergyFunc&&                energyFunc,
        ForceFunc&&                 forceFunc,
        floatingpoint               MAXDIST,
        floatingpoint               maxForce,
        floatingpoint               LAMBDAMAX
    ) {
        using namespace std;

        LineSearchResult res;

        //@{ Lambda phase 1
        floatingpoint lambda;
        bool sconvergencecheck = false;
        //return zero if no forces
        if(maxForce == 0.0) {
            return res;
        }

        //calculate initial guess for lambda
        floatingpoint lambdacap = min(LAMBDAMAX, MAXDIST / maxForce);
        lambda = lambdacap;

        //@} Lambda phase 1
        tbegin = chrono::high_resolution_clock::now();
        floatingpoint currentEnergy = energyFunc(coord.data());
        ++res.numEnergyCall;
        floatingpoint Energyi_1 = currentEnergy;
        CUDAcommon::tmin.computeenerycallszero++;
        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
        CUDAcommon::tmin.computeenergy+= elapsed_energy.count();
        CUDAcommon::tmin.computeenergyzero+= elapsed_energy.count();

        #ifdef TRACKDIDNOTMINIMIZE
        SysParams::Mininimization().TotalE.push_back(currentEnergy);
        #endif

        double FDotFA = searchDirDotForce();
        double FDotFAprev = FDotFA;
        double FDotFAnext;
        double Del_FDotFA;
        FDotFAprev = FDotFA;
        floatingpoint lambdaprev = 0.0;
        int iter = 0;
        floatingpoint energyChange = (floatingpoint)0.0;
        floatingpoint energyLambda = (floatingpoint)0.0;

        while(!(sconvergencecheck)) {
    //		cout<<"starting with lambda "<<lambda<<endl;
            iter++;
        // @{ Lambda phase 2
            moveAlongSearchDir(coord, lambda);
            energyLambda = energyFunc(coordLineSearch.data());
            ++res.numEnergyCall;
            CUDAcommon::tmin.computeenerycallsnonzero++;
        //Step1: Calculate Forces & Dot products
            forceFunc.computeForces(coordLineSearch.data(), force.data(), force.size());
            ++res.numForceCall;
            FDotFAnext = searchDirDotForce();
            Del_FDotFA = FDotFAnext - FDotFAprev;
            floatingpoint idealEnergyChange = -BACKTRACKSLOPE * lambda * FDotFA;
        //Step3: Calculate Lambdaquad
            floatingpoint relerr = fabs(1.0-(0.5*(lambda-lambdaprev)*(FDotFAnext+FDotFAprev)
                    +energyLambda) /Energyi_1);

            floatingpoint lambdaquad = lambda - FDotFAprev*(lambda-lambdaprev)/Del_FDotFA;

        //Step 4. Test if the energy given by alphaquad is within the acceptable range
        //prescribed by Armijo condition on alpha
            if(relerr <= QUADTOL && lambdaquad > 0 && lambdaquad < lambdacap &&
            lambdaquad > LAMBDATOL) {
                tbegin = chrono::high_resolution_clock::now();
                moveAlongSearchDir(coord, lambdaquad);
                floatingpoint energyLambdaquad = energyFunc(coordLineSearch.data());
                ++res.numEnergyCall;
                CUDAcommon::tmin.computeenerycallsnonzero++;
                tend = chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
                CUDAcommon::tmin.computeenergy += elapsed_energy.count();
                CUDAcommon::tmin.computeenergynonzero += elapsed_energy.count();

                energyChange = energyLambdaquad - currentEnergy;
//				cout << "Ideal " << idealEnergyChange << " quad energy " << energyChange <<endl;
                //if it satisfies, set lambda to lambdaquad and return.
                if (energyChange <= idealEnergyChange) {
                    sconvergencecheck = true;
                    lambda = lambdaquad;
                    moveAlongSearchDir(coord, lambda);
                    forceFunc(coordLineSearch.data(), force.data(), force.size());
                    ++res.numForceCall;
                }
                else
                    Energyi_1 = energyLambdaquad;
            }
            else{
                //check if normal backtracking works
                energyChange = energyLambda - currentEnergy;
//				cout << "Ideal " << idealEnergyChange << " bt energy " << energyChange << endl;
                //if it satisfies, set lambda to lambdaquad and return.
                if (energyChange <= idealEnergyChange) {
                    sconvergencecheck = true;
                }
                else
                    Energyi_1 = energyLambda;
            }

            //if not try again
            if(!sconvergencecheck) {
                FDotFAprev = FDotFAnext;
                lambdaprev = lambda;
                //reduce lambda
                lambda *= LAMBDAREDUCE;

            }
            //If backtracked all the way, exit
            if(lambda <= 0.0 || lambda <= LAMBDATOL) {
                sconvergencecheck = true;
                lambda = 0.0;
                moveAlongSearchDir(coord, lambda);
                forceFunc(coordLineSearch.data(), force.data(), force.size());
                ++res.numForceCall;
            }
            //@{ Lambda phase 2
        }

        if(lambda > 0)
            TotalEnergy = energyLambda;
        else
            TotalEnergy = currentEnergy;

        res.lambda = lambda;
        res.newEnergy = energyLambda;
        return res;
    }

    // Note:
    //   - EnergyFunc is a function of signature (floatingpoint* coord) -> floatingpoint.
    template< typename EnergyFunc >
	medyan::LineSearchResult quadraticLineSearchV2(
        std::vector<floatingpoint>& coord,
        EnergyFunc&&                energyFunc,
        floatingpoint               MAXDIST,
        floatingpoint               maxForce,
        floatingpoint               LAMBDAMAX
    ) {
        using namespace std;
        using namespace medyan;

        LineSearchResult res;

        //@{ Lambda phase 1
        floatingpoint lambda;
        bool sconvergencecheck = false;
        //return zero if no forces
        if(maxForce == 0.0) {
            return res;
        }

        //calculate first lambda
        lambda = min(LAMBDAMAX, MAXDIST / maxForce);

        //@} Lambda phase 1
        tbegin = chrono::high_resolution_clock::now();
        floatingpoint currentEnergy = energyFunc(coord.data());
        ++res.numEnergyCall;
        CUDAcommon::tmin.computeenerycallszero++;
        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
        CUDAcommon::tmin.computeenergy+= elapsed_energy.count();
        CUDAcommon::tmin.computeenergyzero+= elapsed_energy.count();

        #ifdef TRACKDIDNOTMINIMIZE
        SysParams::Mininimization().TotalE.push_back(currentEnergy);
        #endif

        floatingpoint energyChange = (floatingpoint)0.0;
        floatingpoint energyLambda = (floatingpoint)0.0;
        floatingpoint prevEnergy = (floatingpoint) 0.0;
        floatingpoint prevLambda = (floatingpoint) 0.0;
        vector<floatingpoint> lambdavec;
        vector<floatingpoint> energyvec;

        int iter = 0;
        while(!(sconvergencecheck)) {

            iter++;

            tbegin = chrono::high_resolution_clock::now();
            moveAlongSearchDir(coord, lambda);
            energyLambda = energyFunc(coordLineSearch.data());
            ++res.numEnergyCall;
            CUDAcommon::tmin.computeenerycallsnonzero++;
            tend = chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
            CUDAcommon::tmin.computeenergy+= elapsed_energy.count();
            CUDAcommon::tmin.computeenergynonzero+= elapsed_energy.count();

            //@{ Lambda phase 2
            floatingpoint idealEnergyChange = -BACKTRACKSLOPE * lambda * searchDirDotForce();

            if(idealEnergyChange>0) {
                cout<<"Ideal Energy Change is positive. Exiting."<<endl;
                exit(EXIT_FAILURE);
            }

            energyChange = energyLambda - currentEnergy;

            //return if ok
            //Armijo conditon
            if(energyChange <= idealEnergyChange)
                sconvergencecheck = true;
            else {
                prevLambda = lambda;
                prevEnergy = energyLambda;
                //reduce lambda
                lambda *= LAMBDAREDUCE;
            }

            if(lambda <= 0.0 || lambda <= LAMBDATOL) {
                sconvergencecheck = true;
                lambda = 0.0;
                lambdavec.push_back(lambda);
                energyvec.push_back(currentEnergy);
            }
            else if(sconvergencecheck){
                lambdavec.push_back(lambda);
                energyvec.push_back(energyLambda);
            }
            else{
                lambdavec.push_back(prevLambda);
                energyvec.push_back(prevEnergy);
            }
        }

        // Try quadratic optimization

        if(lambda > 0.0 && iter > 1 && energyLambda < currentEnergy){
            bool quadstatus = false;
            int i = 0;
            for(i=lambdavec.size()-1; i >=1; i--){
                if(energyLambda < energyvec[i]){
                    quadstatus = true;
                }
            }
            if(quadstatus) {
                vector<floatingpoint> quadlambdavec;
                vector<floatingpoint> quadenergyvec;
                quadlambdavec.push_back(floatingpoint(0.0));
                quadlambdavec.push_back(lambda);
                quadlambdavec.push_back(lambdavec[i]);
                quadenergyvec.push_back(currentEnergy);
                quadenergyvec.push_back(energyLambda);
                quadenergyvec.push_back(energyvec[i]);
                floatingpoint lambdaquad = quadraticoptimization(
                    coord, energyFunc, res.numEnergyCall,
                    quadlambdavec, quadenergyvec
                );
                if(lambdaquad <= 0.0 || lambdaquad <= LAMBDATOL)
                    lambda = 0.0;
                else
                    lambda = lambdaquad;
            }
        }


        if(lambda > 0)
            TotalEnergy = energyLambda;
        else
            TotalEnergy = currentEnergy;

        res.lambda = lambda;
        res.newEnergy = energyLambda;
        return res;
    }

    // Note:
    //   - EnergyFunc is a function of signature (floatingpoint* coord) -> floatingpoint.
    template< typename EnergyFunc >
	floatingpoint quadraticoptimization(
        std::vector<floatingpoint>& coord,
        EnergyFunc&&                energyFunc,
        int&                        numEnergyCall,
        const vector<floatingpoint>& lambdavec,
        const vector<floatingpoint>& energyvec
    ) {
        using namespace std;

        floatingpoint x0, x1, x2, y0, y1, y2;
        x0 = lambdavec[0];
        x1 = lambdavec[1];
        x2 = lambdavec[2];
        y0 = energyvec[0];
        y1 = energyvec[1];
        y2 = energyvec[2];
        double btlambda = lambdavec[1];
        double lambdaquad = lambdavec[1];
        double energyQuad = 0.0;
        //Mq(l) = a*l*l + b*l +c;
        //Refer quadratic optimization doc.
        int iter = 0;
    //	cout<<"btlambda "<<btlambda<<endl;
        while(true) {
    /*		cout<<"Trial t"<<iter<<"Lambda = ["<<x0<<" "<<x1<<" "<<x2<<"];t"<<iter<<"Energies = ["<<
            y0<<" "<<y1<<" "<<y2<<"];"<<endl;*/
            iter++;
            double d = (x0-x1)*(x0-x2)*(x1-x2);
            double a = y0 * (x1 - x2) - y1 * (x0 - x2) + y2 * (x0 - x1);
            double b = -(y0 * (x1 + x2) * (x1 - x2) - y1 * (x0 + x2) * (x0 - x2) +
                        y2 * (x0 + x1) * (x0 - x1));
    /*		double c = y0 * x1 * x2 * (x1 - x2) - y1 * x0 * x2 * (x0 - x2) +
                    y2 * x0 * x1 * (x0 - x1);

            cout<<"Mq"<<iter-1<<" = "<<a/d<<"*x.*x + "<<b/d<<"*x + "<<c/d<<";"<<endl;
            cout<<"-b/2a = "<<-b/(2*a)<<";"<<endl;*/

            lambdaquad = -b/(2*a);
    //		cout<<"["<<abs(lambdaquad-x1)/x1<<" "<<LAMBDAQUADTOL<<"];"<<endl;


            //Check if lambda has converged
            if (abs(lambdaquad-x1)/x1 < LAMBDAQUADTOL) {
                return lambdaquad;
            }
            else if(a/d < 0){
                cout<<"WARNING! Quadratic model does not have a minima. Returning "
            "BACKTRACKING Lambda "<<endl;
                return btlambda;
            }
            //else compute Energy at lambdaquad
            else {
                moveAlongSearchDir(coord, lambdaquad);
                energyQuad = energyFunc(coordLineSearch.data());
                ++numEnergyCall;
                CUDAcommon::tmin.computeenerycallsnonzero++;
    //			cout<<"["<<abs(energyQuad-y1)/y1<<" "<<LAMBDAQUADTOL<<"];"<<endl;
                if(abs(energyQuad-y1)/y1 < LAMBDAQUADTOL){
                    return lambdaquad;
                }

                //Determine points to use in next iteration.
                //determine if lambdaquad is to the left or right of x1
                if(lambdaquad<x1){
                    if(energyQuad<y0 && energyQuad<y1){
                        x1 = lambdaquad;
                        y1 = energyQuad;
                    }
                    else{
    //					if(energyQuad<y0 && energyQuad>y1){
                        x0 = lambdaquad;
                        y0 = energyQuad;
                    }
    /*				else{
                        cout<<"WARNING! Unreasonable lambdaquad during quadratic optimization"
                                ". Returning BACKTRACKING lambda"<<endl;
                        return btlambda;
                    }*/
                }
                else if(x1 < lambdaquad){
                    if(energyQuad<y1 && energyQuad<y2){
                        x1 = lambdaquad;
                        y1 = energyQuad;
                    }
                    else {
    //					if(energyQuad>y1 && energyQuad<y2){
                        x2 = lambdaquad;
                        y2 = energyQuad;
                    }
    /*				else{
                        cout<<"WARNING! Unreasonable lambdaquad during quadratic optimization"
                            ". Returning BACKTRACKING lambda"<<endl;
                        return btlambda;
                    }*/
                }
                else{
                    cout<<"WARNING! Lambda out of range. Returning quadlambda from "
                            "previous iteration"<<endl;
                    return x1;
				}
            }
        }

    }

    void setLAMBDATOL(int maxF_order){

        int orderdimension = 3; ///1000s of nm
        //Float gives you a minimum of 9 sig figs. If you operate in a 10^3nm system, the
        // decimal part can go upto 10^-6. In our system, Lambda*F_i should not be
        // greater than 10^-6. We would like to be cautious and ensure that all numbers
        // have  to the order of 10^-3. Hence, O(Lambda*F_i) >= 10^-(9-3-3) = 10^-3
        int LAMBDATOLorder = -(9-orderdimension -3) - maxF_order;
        LAMBDATOL = 1;
        if(LAMBDATOLorder > 0){
            for(int i =0; i < LAMBDATOLorder; i ++)
                LAMBDATOL *= 10;
        }
        else{
            for(int i =0; i > LAMBDATOLorder; i --)
                LAMBDATOL *= 0.1;
        }
		//Since, wthe force threshold are in 10^0 of pN at the lowest range, our lambda
		// should be in the order of 10^-5.
        LAMBDATOL = max<floatingpoint>(1e-8, LAMBDATOL);
        LAMBDATOL = min<floatingpoint>(1e-5, LAMBDATOL);

//        cout<<"maxF order "<<maxF_order<<" lambdatol "<<LAMBDATOL<<endl;
    }

    void setBACKTRACKSLOPE(floatingpoint _btslope){BACKTRACKSLOPE = _btslope;}

    void copycoordsifminimumE(floatingpoint maxForce, const std::vector<floatingpoint>& coord){

    	if(TotalEnergy <= minimumE){
    		//update minimum energy
    		minimumE = TotalEnergy;
    		maxForcebackup = maxForce;
    		//take backup of coordinates.
            coordMinE = coord;
    	}
    }

    void copybackupcoordinates(std::vector<floatingpoint>& coord){

    	if(coordMinE.size()) {
		    cout<<"Copying coordinates with the lowest energy during minimization "<<endl;
		    cout<<"Energy = "<<minimumE<<" pN.nm"<<endl;
		    cout<<"MaxForce = "<<maxForcebackup<<" pN "<<endl;
            coord = std::move(coordMinE);
            coordMinE.clear();
	    }
    }

    //@}

public:

    /// Minimize the system
    // Note:
    //   - EnergyFunc:           floatingpoint* coord     ->   floatingpoint
    //   - ForceFunc:            floatingpoint* coord, floatingpoint* force, int numVar   ->   void
    //   - EnergyFuncIndividual: floatingpoint*           ->   EnergyReport
    //   - FuncRecoveryOnError:  vector<floatingpoint>& coord, const vector<floatingpoint>& searchDir  ->   void
    template< ConjugateGradientDescentSearch method, typename EnergyFunc, typename ForceFunc, typename EnergyFuncIndividual, typename FuncRecoveryOnError >
    MinimizationResult minimize(
        const ConjugateGradientParams& cgParams,
        std::vector< floatingpoint >&  coord,
        int                            numDof,
        EnergyFunc&&                   energyFunc,
        ForceFunc&&                    forceFunc,
        EnergyFuncIndividual&&         energyFuncIndividual,
        FuncRecoveryOnError&&          funcRecoveryOnError,
        std::vector< floatingpoint >*  ptrOutForce
    ) {
        using namespace std;

        #ifdef TRACKDIDNOTMINIMIZE
        SysParams::Mininimization().beta.clear();
        SysParams::Mininimization().Lambda.clear();
        SysParams::Mininimization().Energyvec.clear();
        SysParams::Mininimization().TotalE.clear();
        SysParams::Mininimization().maxF.clear();
        SysParams::Mininimization().safeModeORnot.clear();
        SysParams::Mininimization().tempEnergyvec.clear();
        SysParams::Mininimization().gradientvec.clear();
        #endif

        MinimizationResult result;
        result.energyCallLimit = cgParams.energyCallLimit;
        result.gradTol         = cgParams.gradTol;
        result.energyRelTol    = cgParams.energyRelTol;

        {
            int beadMaxStep = 3 * Bead::numBeads();
            result.forceCallLimit = std::max<int>(cgParams.forceCallLimit, beadMaxStep);
        }

        // Temporary variables for debug purposes.
        //---------------------------------
        std::vector<floatingpoint> coordBackup;
        std::vector<floatingpoint> forceBackup;
        floatingpoint prevlambda = 0;
        floatingpoint prevbeta = 0;
        chrono::high_resolution_clock::time_point tbegin, tend;

        //@@@{ STEP 1: Start minimization
        tbegin = chrono::high_resolution_clock::now();
        startMinimization();


#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        // Caches all the data to the CGMethod interval data vectors.
        //---------------------------------
        this->numDof = numDof;
        forcePrev.assign(numDof, 0);
        searchDir.assign(numDof, 0);
        forceTol.resize(numDof, cgParams.gradTol); // Future: can set different default grad tol

        const auto nvar = coord.size();
        coordLineSearch.assign(nvar, 0);
        force.assign(nvar, 0);


        // Initial energy and force computation.
        //------------------------------
        result.energiesBefore = energyFuncIndividual(coord.data());
        ++result.numEnergyCall;
        forceFunc(coord.data(), force.data(), force.size());
        ++result.numForceCall;

        result.gradInfNorm = maxF();

        auto prevEnergy = result.energiesBefore.total;


#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
#ifdef OPTIMOUT
        floatingpoint lambdatime = 0.0;
        int safestatuscount = 0;
#endif

#ifdef CUDAACCL
        cross_checkclass::Aux=false;
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        #ifdef OPTIMOUT
        cout<<"Energy before minimization"<<endl;
        energyFunc(coord.data());
        CUDAcommon::tmin.computeenerycallszero++;
        cout<<endl;
        #endif
        #ifdef ADDITIONALINFO
        energyFunc(coord.data(), true);
        if(SysParams::RUNSTATE == false){
        int counter = 0;
        auto individualenergiesvec = MotorGhostInteractions::individualenergies;
        auto tpdistvec = MotorGhostInteractions::tpdistvec;
        auto eqlvec = MotorGhostInteractions::eqlvec;
        auto kstrvec = MotorGhostInteractions::kstrvec;
        for(auto l :MotorGhost::getMotorGhosts()) {
            Cylinder *cyl1 = l->getFirstCylinder();
            Cylinder *cyl2 = l->getSecondCylinder();
            float pos1 = l->getFirstPosition() * SysParams::Geometry()
                    .cylinderNumMon[cyl1->getType()];
            float pos2 = l->getSecondPosition() * SysParams::Geometry()
                    .cylinderNumMon[cyl2->getType()];
            cout << l->getId() << " " << l->getType() << " " << cyl1->getStableIndex() << " "
                << cyl2->getStableIndex() << " " << pos1 << " " << pos2 << " "
                << l->getMMotorGhost()->getEqLength() << " " << l->getCMotorGhost()
                        ->getDiffusingSpecies()->getName() << " " << l->getNumHeads() << " "
                << l->getnumBoundHeads() << " " << kstrvec[counter] << " " <<
                eqlvec[counter] << " " << tpdistvec[counter] << " "
                                                                ""
                << individualenergiesvec[counter] << endl;

            counter++;
        }
        }
        #endif

        // State variables during minimization.
        //-----------------------------
        tbegin = chrono::high_resolution_clock::now();
        // Safe mode can be enabled to relax the line search ending criteria.
        bool safeMode = false;
        int numIter = 0;
        // force-dot-force for last iteration.
        floatingpoint lastFDotF = 0;


        // Auxiliary function to reset conjugate gradient states.
        const auto resetCGParams = [&, this] {
            shiftSearchDir(0);
            if constexpr(method == ConjugateGradientDescentSearch::fletcherRieves || method == ConjugateGradientDescentSearch::polakRibiere) {
                lastFDotF = forceDotForce();
            }
            if constexpr(method == ConjugateGradientDescentSearch::polakRibiere) {
                forcePrev = force;
            }
        };
        resetCGParams();


#ifdef CUDAACCL
        volatile bool *Mc_isminimizationstate;
        volatile bool *Mc_issafestate;
#endif


#ifdef  CUDAACCL
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));


        if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&stream1));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&stream2));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || stream1 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&stream3));
        FFM.CUDAcopyForces(stream1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);//pass a
        // stream
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        FFM.CUDAcopyForces(stream2, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_force);//pass a
        // stream
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        floatingpoint *gpu_GRADTOL;
        floatingpoint gradtol[1];
        gradtol[0]= cgParams.gradTol;
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_GRADTOL, sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemcpy(gpu_GRADTOL, gradtol, sizeof(floatingpoint), cudaMemcpyHostToDevice));
        CGMethod::CUDAallFDotF(stream3);//curGrad //pass a stream

        //synchronize streams
        CUDAcommon::handleerror(cudaStreamSynchronize(stream1));
        CUDAcommon::handleerror(cudaStreamSynchronize(stream2));
        CUDAcommon::handleerror(cudaStreamSynchronize(stream3));
        if(!(CUDAcommon::getCUDAvars().conservestreams)) {
            CUDAcommon::handleerror(cudaStreamDestroy(stream1));
            CUDAcommon::handleerror(cudaStreamDestroy(stream2));
            CUDAcommon::handleerror(cudaStreamDestroy(stream3));
        }
    //PING PONG
        bool  *Mmh_stop, *Mmg_stop1, *Mmg_stop2, *Mmg_s1, *Mmg_s2, *Mmg_ss;//minimization state
        bool  *Msh_stop, *Msg_stop1, *Msg_stop2, *Msg_s1, *Msg_s2, *Msg_ss;//safe state

        //PING PONG
        //minimization state
        cudaMalloc(&Mmg_stop1, sizeof(bool));
        cudaMalloc(&Mmg_stop2, sizeof(bool));
        cudaHostAlloc(&Mmh_stop, sizeof(bool), cudaHostAllocDefault);
        //safe state
        cudaMalloc(&Msg_stop1, sizeof(bool));
        cudaMalloc(&Msg_stop2, sizeof(bool));
        cudaHostAlloc(&Msh_stop, sizeof(bool), cudaHostAllocDefault);
        //@

        if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms1 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&Ms1));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms2 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&Ms2));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms3 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&Ms3));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || Ms4 == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&Ms4));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || Me1 == NULL)
            CUDAcommon::handleerror(cudaEventCreate(&Me1));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || Me2 == NULL)
            CUDAcommon::handleerror(cudaEventCreate(&Me2));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || stream_shiftsafe == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&stream_shiftsafe));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || stream_dotcopy == NULL)
            CUDAcommon::handleerror(cudaStreamCreate(&stream_dotcopy));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || event_safe == NULL)
            CUDAcommon::handleerror(cudaEventCreate(&event_safe));
        if(!(CUDAcommon::getCUDAvars().conservestreams) || event_dot == NULL)
            CUDAcommon::handleerror(cudaEventCreate(&event_dot));

        Mmh_stop[0] = true; //Minimizationstate //Yes = Minimize. No = Don't minimize.
        Msh_stop[0] = false; //safe state
        Mc_isminimizationstate = Mmh_stop;//points to address of Mmh_stop
        Mc_issafestate = Msh_stop;//points to address of Msh_stop

        Msp1 = &Ms1;
        Msp2 = &Ms2;
        Mep1 = &Me1;
        Mep2 = &Me2;
        Mmg_s1 = Mmg_stop1;
        Mmg_s2 = Mmg_stop2;
        Msg_s1 = Msg_stop1;
        Msg_s2 = Msg_stop2;
    //set Mmg_stop1, Mmg_stop2 to true and Msg_stop1, Msg_stop2 to false.

    // stick to single stream going forward.
    //@CUDA Get minimizaton state{
        //calculate MAXF and allFDotFA
        CGMethod::CUDAinitializePolak(*Msp1, Mmg_s1, Mmg_s2, Msg_s1, Msg_s2);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
        CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGPolakRibiereMethod.cu");
#ifdef SERIAL_CUDACROSSCHECK
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        std::cout<<"FMAX SERL "<<maxF()<<endl;
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
        //Copy to host
        CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
        cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
    //@}
#endif

#ifdef DETAILEDOUTPUT
        std::cout<<"printing beads & forces"<<endl;
        long i = 0;
        long index = 0;
        for(auto b:Bead::getBeads()){
            std::cout<<b->getId()<<" "<< b->coord <<" "
                    "" << b->force <<endl;
        }
        std::cout<<"printed beads & forces"<<endl;
#endif


        //CUDA based minimization begins
#ifdef CUDAACCL
        while (/* Iteration criterion */  numIter < N &&
            /* Gradient tolerance  */  (Mc_isminimizationstate[0])) {

#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
    //PING PONG SWAP
    //        CUDAcommon::handleerror(cudaStreamWaitEvent(*Msp2, *Mep1, 0));

            CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));

            Msps = Msp1;
            Msp1 = Msp2;
            Msp2 = Msps;
            Meps = Mep1;
            Mep1 = Mep2;
            Mep2 = Meps;
            Mmg_ss = Mmg_s1;
            Mmg_s1 = Mmg_s2;
            Mmg_s2 = Mmg_ss;
            Msg_ss = Msg_s1;
            Msg_s1 = Msg_s2;
            Msg_s2 = Msg_ss;
    //PING ENDS
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            numIter++;
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
            std::cout<<"SL safestate "<<safeMode<<endl;
#endif
            if(Mc_issafestate[0]) {
                safeMode = false;
            }
            //find lambda by line search, move beads
            lambda = backtrackingLineSearchCUDA(FFM, cgParams.maxDist, cgParams.lambdaMax, Msg_s1);
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif

            CUDAmoveBeads(*Msp1, Mmg_s1);
    //        CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));//This seems unnecessary.
            //wait for movebeads to finish before calculating forces1

#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            //Synchronize to ensure forces are calculated on moved beads.
            CUDAcommon::handleerror(cudaStreamSynchronize(*Msp1));
    //        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));
            //@CUDA Get minimizaton state
            // @{
    //        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
            CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));
            //@}

#if defined(CROSSCHECK) || defined(CUDAACCL)
            cross_checkclass::Aux=true;
#endif
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif

            cvars = CUDAcommon::getCUDAvars();
            cvars.streamvec.clear();
            CUDAcommon::cudavars = cvars;
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_dotcopy));

            //compute new forces
            forceFunc(coord.data(), force.data(), force.size());//split and synchronize
#ifdef DETAILEDOUTPUT
            // wARNING This output is no longer safe because it assumes bead
            // coordinates start with index 0
            std::cout<<"MB printing beads & forces L "<<lambda<<endl;
            long i = 0;
            long index = 0;
            for(auto b:Bead::getBeads()){
                index = 3 * b->getIndex();

                std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
                        ""<<coord[index + 2]<<" "
                        ""<<forceAux[index]<<" "
                        ""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
                        b->getIndex()<<endl;
            }
            std::cout<<"MB printed beads & forces"<<endl;
#endif

#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            //wait for forces to be calculated
            for(auto strm:CUDAcommon::getCUDAvars().streamvec)
                CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
            //compute direction CUDA
            CGMethod::CUDAallFADotFA(stream_dotcopy); //newGrad
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            CGMethod::CUDAallFADotFAP(stream_dotcopy); //prevGrad
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            CUDAcommon::handleerror(cudaEventRecord(event_dot,stream_dotcopy));
            CUDAcommon::handleerror(cudaStreamWaitEvent(stream_shiftsafe, event_dot,0));
    //Copy forces
            FFM.CUDAcopyForces(stream_dotcopy, CUDAcommon::getCUDAvars().gpu_forceAuxP,CUDAcommon::getCUDAvars().gpu_forceAux);
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            //Polak-Ribieri update beta & shift gradient
            CUDAshiftGradient(stream_shiftsafe, Mmg_s1);
#ifdef SERIAL_CUDACROSSCHECK
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_shiftsafe));
#endif
    /*        CUDAcommon::handleerror(cudaStreamSynchronize(*Msp2));
            //@CUDA Get minimizaton state{
    //        CGMethod::CUDAgetPolakvars(true, *Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Msg_s2, Mc_isminimizationstate);
            CGMethod::CUDAgetPolakvars(*Msp1, gpu_GRADTOL, Mmg_s1, Mmg_s2, Mc_isminimizationstate);
    #ifdef ALLSYNC
            cudaDeviceSynchronize();
    #endif
            CUDAcommon::handleerror(cudaEventRecord(*Mep1, *Msp1));*/
            CGMethod::CUDAgetPolakvars2(stream_shiftsafe, Msg_s2);
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            CUDAshiftGradientifSafe(stream_shiftsafe, Mmg_s1, Msg_s1);
#ifdef ALLSYNC
            cudaDeviceSynchronize();
#endif
            if(Mc_isminimizationstate[0]  == true){
                //Copy to host
                CUDAcommon::handleerror(cudaStreamWaitEvent(Ms3, *Mep1, 0));
    //            CUDAcommon::handleerror(cudaStreamWaitEvent(Ms4, event_safe, 0));//event_safe is not attached to any event
#ifdef ALLSYNC
                cudaDeviceSynchronize();
#endif
                cudaMemcpyAsync(Mmh_stop, Mmg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms3);
#ifdef ALLSYNC
                cudaDeviceSynchronize();
#endif
                //TODO remove later... June 7, 2018. Removed.
    //            std::cout<<"safe state copy"<<endl;
    //            cudaMemcpyAsync(Msh_stop, Msg_s2, sizeof(bool), cudaMemcpyDeviceToHost, Ms4);

            }
        }

        std::cout<<"CUDA Total number of iterations "<<numIter<<endl;
#endif //CUDAACCL
//CUDA based minimization ends

        #ifdef TRACKDIDNOTMINIMIZE
        SysParams::Mininimization().maxF.push_back(maxForce);
        #endif

        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_other(tend - tbegin);
        CUDAcommon::tmin.tother+= elapsed_other.count();
        //@@@} STEP 4 OTHER

        while(
            !result.callLimitExceeded() &&
            !result.converged()
        ) {

            //@@@{ STEP 5 OTHER
            tbegin = chrono::high_resolution_clock::now();
            if (std::is_same<floatingpoint, float>::value) {
                //set the floor of lambda (lowest lambda allowed based on maxf
                int maxForder = static_cast<int>(floor(log10(result.gradInfNorm)));
                if (maxForder < 0) maxForder--;
                CGMethod::setLAMBDATOL(maxForder);
            }
#ifdef CROSSCHECK_CYLINDER
            CGMethod::_crosscheckdumpMechFile<<"Lambda order set"<<endl;
#endif

            tend = chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_other2(tend - tbegin);
            CUDAcommon::tmin.tother += elapsed_other2.count();
            //@@@} OTHER

            // Conjugate gradient related state variables.
            //-------------------------
            double beta = 0;
            double newGrad = 0, prevGrad = 0;

            numIter++;
#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_LAMBDA)
            std::cout<<"SL safestate "<<safeMode<<endl;
#endif
#ifdef TRACKDIDNOTMINIMIZE
            SysParams::Mininimization().safeModeORnot.push_back(safeMode);
#endif

            //@@@{ STEP 6 FIND LAMBDA
            tbegin = chrono::high_resolution_clock::now();
            medyan::LambdaRunningAverageManager lambdaAvg;
            lambdaAvg.lambdaRunningAverageProbability = cgParams.lambdaRunningAverageProbability;
            medyan::LineSearchResult lsr;
            if (cgParams.lineSearchAlgorithm == "BACKTRACKING") {
                lsr = backtrackingLineSearch(coord, energyFunc, cgParams.maxDist, result.gradInfNorm, cgParams.lambdaMax, &lambdaAvg, safeMode ? 0 : BACKTRACKSLOPE);

            } else if (cgParams.lineSearchAlgorithm == "QUADRATIC") {
                lsr = safeMode
                    ? backtrackingLineSearch(coord, energyFunc, cgParams.maxDist, result.gradInfNorm, cgParams.lambdaMax, &lambdaAvg, 0)
                    : quadraticLineSearchV2(coord, energyFunc, cgParams.maxDist, result.gradInfNorm, cgParams.lambdaMax);

            } else {
                lsr = safeMode
                    ? backtrackingLineSearch(coord, energyFunc, cgParams.maxDist, result.gradInfNorm, cgParams.lambdaMax, &lambdaAvg, 0)
                    : quadraticLineSearchV2(coord, energyFunc, cgParams.maxDist, result.gradInfNorm, cgParams.lambdaMax);

            }
            const auto lambda = lsr.lambda;
            result.numEnergyCall += lsr.numEnergyCall;
            result.numForceCall += lsr.numForceCall;

            if(lambda == 0) {
                // Line search error.
                if(cgParams.tryToRecoverInLineSearchError) {
                    if(safeMode) {
                        // Even in safe mode, line search failed to find a positive step size.
                        log::debug("During energy minimization, line search failed to find a lower energy state.");

                        // Try to kick the system out of current configuration.
                        funcRecoveryOnError(coord, searchDir);

                        // Reset force and conjugate gradient parameters.
                        forceFunc(coord.data(), force.data(), force.size());
                        ++result.numForceCall;
                        result.gradInfNorm = maxF();
                        resetCGParams();

                        safeMode = false;
                    }
                    else {
                        // Reset conjugate gradient parameters and enable safe mode.
                        resetCGParams();
                        safeMode = true;
                    }

                }
                else {
                    // Abort minimization.
                    result.errorFlags |= MinimizationResult::errorLineSearch;
                    return result;
                }
            }
            else {
                // Valid positive lambda.

                #ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"Lambda found"<<endl;
                #endif

                tend = chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_lambda(tend - tbegin);
                CUDAcommon::tmin.findlambda += elapsed_lambda.count();
                #ifdef OPTIMOUT
                lambdatime += elapsed_lambda.count();
                #endif

#ifdef TRACKDIDNOTMINIMIZE
                SysParams::Mininimization().Lambda.push_back(lambda);
#endif


#ifdef SERIAL_CUDACROSSCHECK
                CUDAcommon::handleerror(cudaDeviceSynchronize());
                floatingpoint cuda_lambda[1];
                CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(floatingpoint),
                                                    cudaMemcpyDeviceToHost));
                std::cout<<"Lambda CUDA "<<cuda_lambda[0]<<" SERL "<<lambda<<endl;
#endif
                //@@@{ STEP7 OTHER
                {
#if defined(TRACKDIDNOTMINIMIZE) || defined(EVSALPHA)
                    //Backup coordinate
                    coordBackup = coord;
                    forceBackup = force;
                    calculateEvsalpha(energyFunc, lambda, cgParams.lambdaMax, searchDirDotForce(), prevlambda, coordBackup, forceBackup);
                    prevlambda = lambda;
                    cout<<endl;

#endif
                    tbegin = chrono::high_resolution_clock::now();
                    moveAlongSearchDir(coord, lambda);
                    tend = chrono::high_resolution_clock::now();
                    chrono::duration<floatingpoint> elapsed_other3(tend - tbegin);
                    CUDAcommon::tmin.tother += elapsed_other3.count();
                }
                #ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"Beads moved"<<endl;
                #endif
                //@@@} OTHER
#if defined(CROSSCHECK) || defined(CUDAACCL)
                cross_checkclass::Aux=true;
#endif
                ///@@@{ STEP 8 compute new forces
                tbegin = chrono::high_resolution_clock::now();
                forceFunc(coord.data(), force.data(), force.size());
                ++result.numForceCall;
                tend = chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_force(tend - tbegin);
                CUDAcommon::tmin.computeforces += elapsed_force.count();
                #ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"Compute force set"<<endl;
                #endif

                // Update convergence criteria.
                result.gradInfNorm = maxF();
                result.energyRelChange = std::abs(lsr.newEnergy - prevEnergy) / std::abs(prevEnergy);
                prevEnergy = lsr.newEnergy;


#ifdef DETAILEDOUTPUT
                // wARNING This output is no longer safe because it assumes bead
                // coordinates start with index 0
                std::cout<<"MB printing beads & forces L "<<lambda<<endl;
                long i = 0;
                long index = 0;
                for(auto b:Bead::getBeads()){
                    index = 3 * b->getIndex();

                    std::cout<<b->getId()<<" "<<coord[index]<<" "<<coord[index + 1]<<" "
                            ""<<coord[index + 2]<<" "
                            ""<<forceAux[index]<<" "
                            ""<<forceAux[index + 1]<<" "<<forceAux[index + 2]<<" "<<3 *
                            b->getIndex()<<endl;
                }
                std::cout<<"MB printed beads & forces"<<endl;
#endif


                if(method == ConjugateGradientDescentSearch::fletcherRieves) {
                    newGrad = forceDotForce();
                }
                else if constexpr(method == ConjugateGradientDescentSearch::polakRibiere) {
                    newGrad = forceDotForce();
                    prevGrad = forceDotForcePrev();
                }


                // Find beta (for updating search direction).
                if constexpr(method == ConjugateGradientDescentSearch::fletcherRieves) {
                    //Fletcher-Rieves update
                    beta = newGrad / lastFDotF;
                }
                else if constexpr(method == ConjugateGradientDescentSearch::polakRibiere) {
                    //Polak-Ribiere update
                    //Max(0,betaPR) allows us to reset the direction under non-ideal circumstances.
                    //The direction is reset of steepest descent direction (-gk).
#ifdef FLOAT_PRECISION
                    double betaPR = max<double>((double) 0.0, (newGrad - prevGrad) / lastFDotF);
                    double betaFR = max<double>((double) 0.0, newGrad / lastFDotF);
                    //Efficient hybrid Conjugate gradient techniques, Eq 21
                    prevbeta = beta;
                    if (betaPR == 0.0)
                        beta = betaFR;
                    else if (betaPR < 1.25 * betaFR)
                        beta = betaPR;
                    else
                        beta = betaFR;
#else
                    beta = max<double>(0.0, (newGrad - prevGrad) / lastFDotF);
#endif
                }
                // otherwise, beta is 0 by default.

#ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"Beta set"<<endl;
#endif
                //Global convergence properties of conjugate gradient methods for optimization Eq
                // 3.8
                //A SURVEY OF NONLINEAR CONJUGATE GRADIENT METHODS Section 9.
                //Allows for negative beta values.
                /*double betaPR = (newGrad - prevGrad) / curGrad;
                double betaFR = newGrad/ curGrad;
                beta = max<double>(-betaFR, min<double>(betaPR, betaFR));
                cout<<"betaPR "<<betaPR<<" betaFR "<<betaFR<<" beta "<<beta<<endl;*/

                // Shift gradient.
                shiftSearchDir(beta);

                // Reset search direction if conjugacy is lost.
                if(searchDirDotForce() <= 0) {
                    resetCGParams();
                }


                #ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"shiftGradient"<<endl;
                #endif

                tend = chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_other4(tend - tbegin);
                CUDAcommon::tmin.tother += elapsed_other4.count();
                //@@@} OTHER

#ifdef TRACKDIDNOTMINIMIZE
                SysParams::Mininimization().beta.push_back(beta);
                SysParams::Mininimization().maxF.push_back(maxForce);
#endif

#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_BETA)
                std::cout<<"Shift Gradient "<<beta<<endl;
                CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
                std::cout<<"Beta serial "<<beta<<endl;
                std::cout<<"newGrad "<<newGrad<<" prevGrad "<<prevGrad<<" curGrad "<<curGrad<<endl;
#endif

                //@@@{ STEP 10 vectorized copy
                if constexpr(method == ConjugateGradientDescentSearch::polakRibiere) {
                    tbegin = chrono::high_resolution_clock::now();
                    forcePrev = force;
                    tend = chrono::high_resolution_clock::now();
                    chrono::duration<floatingpoint> elapsed_copy2(tend - tbegin);
                    CUDAcommon::tmin.copyforces += elapsed_copy2.count();
                }
                //@@@}

                // For any iteration "k",
                // "searchDir" is the conjugate gradient direction (dk)
                // "force" is the force/steepest descent direction (-gk = -grad E)
                // <-gk+1,dk+1> < 0 => gk and dk are at acute angles with one another
                /*Note: Gradient and conjugate direction should be at obtuse angles for effective
                    * descent*/
                //curGrad = newGrad => gk+1.gk+1 and gk.gk are equal. Gradient has not
                // changed in magnitude.
#ifdef TRACKDIDNOTMINIMIZE
                vector<floatingpoint>gradlocal;
                gradlocal.push_back(searchDirDotForce());
                gradlocal.push_back(lastFDotF);
                gradlocal.push_back(newGrad);
                gradlocal.push_back(prevGrad);

                SysParams::Mininimization().gradientvec.push_back(gradlocal);
#endif

                #ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"safeMode set"<<endl;
                #endif
                //Create back up coordinates to go to in case Energy minimization fails at an
                // undeisrable state.
                if (forceBelowRelaxedTolerance(10) && numIter > result.forceCallLimit / 2) {
                    copycoordsifminimumE(result.gradInfNorm, coord);
                }
                #ifdef CROSSCHECK_CYLINDER
                CGMethod::_crosscheckdumpMechFile<<"copycoordsinminimumE"<<endl;
                #endif

                if constexpr(method == ConjugateGradientDescentSearch::fletcherRieves || method == ConjugateGradientDescentSearch::polakRibiere) {
                    lastFDotF = newGrad;
                }

            }
        } // End minimization


#ifdef OPTIMOUT
        std::cout << "SERL Total number of iterations " <<lineSearchAlgorithm<<" "<<
        numIter << endl;

#endif
        if (result.energyRelChangeConverged()) {
            log::info("WARNING: Minimization exited when Energy Tolerance was reached at {} steps.", numIter);
            log::info("Maximum force in system = {}.", result.gradInfNorm);
        }

        if (cgParams.reportOnCallLimit && result.callLimitExceeded()) {

            log::warn("Did not minimize in limited steps.");
            log::info("{} energy calls, {} force calls, max force = {}.",
                result.numEnergyCall, result.numForceCall, result.gradInfNorm);

#ifdef CUDAACCL
            auto cvars = CUDAcommon::getCUDAvars();
            cvars.streamvec.clear();
            CUDAcommon::cudavars = cvars;
#endif
            log::info("System energy... {}", energyFuncIndividual(coord.data()));
            CUDAcommon::tmin.computeenerycallszero++;
#ifdef CUDAACCL
            for(auto strm:CUDAcommon::getCUDAvars().streamvec)
                CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif
            cout << endl;
            //Copy back coordinates that correspond to minimum energy
            copybackupcoordinates(coord);

        }


        #ifdef TRACKDIDNOTMINIMIZE
        if(numIter) {
            auto tempparams = SysParams::Mininimization();

            cout << "Obegin maxForce Lambda Beta SafeModestatus FDotFA curGrad NewGrad prevGrad TotalE Evec ()" << endl;
            for (auto i = 0; i < tempparams.maxF.size()-1; i++) {
                cout << tempparams.maxF[i] << " " << tempparams.Lambda[i] << " " << tempparams
                        .beta[i] << " " << tempparams.safeModeORnot[i] <<" ";
                for(auto j:tempparams.gradientvec[i])
                    cout<< j <<" ";
                cout << tempparams.TotalE[i] << " ";
                for (auto j:tempparams.Energyvec[i]) {
                    cout << j << " ";
                }
                cout << endl;
            }
            cout<<"End maxF "<<tempparams.maxF[tempparams.maxF.size()-1]<<endl;
            cout << "Oend ------------------" << endl;
        }
        #endif

        #ifdef TRACKDIDNOTMINIMIZE
        SysParams::Mininimization().beta.clear();
        SysParams::Mininimization().Lambda.clear();
        SysParams::Mininimization().Energyvec.clear();
        SysParams::Mininimization().TotalE.clear();
        SysParams::Mininimization().maxF.clear();
        SysParams::Mininimization().safeModeORnot.clear();
        SysParams::Mininimization().tempEnergyvec.clear();
        SysParams::Mininimization().gradientvec.clear();
        energyFunc(coord.data());

        #endif

#if defined(CROSSCHECK) || defined(CUDAACCL)
        cross_checkclass::Aux=false;
#endif
#ifdef CUDAACCL
        cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
        CUDAcommon::cudavars = cvars;
#endif
        result.energiesAfter = energyFuncIndividual(coord.data());
        ++result.numEnergyCall;
#ifdef OPTIMOUT
        cout<<"Energy after minimization"<<endl;
        energyFunc(coord.data());
        CUDAcommon::tmin.computeenerycallszero++;
        cout<<endl;
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif
#ifdef CUDAACCL
        for(auto strm:CUDAcommon::getCUDAvars().streamvec)
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm));
#endif

#ifdef CUDAACCL
        FFM.CUDAcopyForces(*Msp1, CUDAcommon::getCUDAvars().gpu_forceAux,CUDAcommon::getCUDAvars().gpu_force);
        //copy back forces and calculate load forces in CPU.
#endif
#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

#ifdef CUDAACCL
        CUDAcommon::handleerror(cudaFreeHost(Mmh_stop));
        CUDAcommon::handleerror(cudaFree(Mmg_stop1));
        CUDAcommon::handleerror(cudaFree(Mmg_stop2));
        CUDAcommon::handleerror(cudaFree(Msg_stop1));
        CUDAcommon::handleerror(cudaFree(Msg_stop2));
        CUDAcommon::handleerror(cudaFree(gpu_GRADTOL));
        CUDAcommon::handleerror(cudaFreeHost(Msh_stop));
        CUDAcommon::handleerror(cudaStreamSynchronize(Ms1));
        CUDAcommon::handleerror(cudaStreamSynchronize(Ms2));
        CUDAcommon::handleerror(cudaStreamSynchronize(Ms3));
        CUDAcommon::handleerror(cudaStreamSynchronize(Ms4));
        if(!(CUDAcommon::getCUDAvars().conservestreams)) {
            CUDAcommon::handleerror(cudaStreamDestroy(Ms1));
            CUDAcommon::handleerror(cudaStreamDestroy(Ms2));
            CUDAcommon::handleerror(cudaStreamDestroy(Ms3));
            CUDAcommon::handleerror(cudaStreamDestroy(Ms4));
            CUDAcommon::handleerror(cudaEventDestroy(Me1));
            CUDAcommon::handleerror(cudaEventDestroy(Me2));
            CUDAcommon::handleerror(cudaEventDestroy(event_safe));
            CUDAcommon::handleerror(cudaEventDestroy(event_dot));
            CUDAcommon::handleerror(cudaStreamDestroy(stream_dotcopy));
            CUDAcommon::handleerror(cudaStreamDestroy(stream_shiftsafe));
        }

#endif
        #ifdef OPTIMOUT
            cout<<"Safestatuscount "<<safestatuscount<<endl;
        #endif

        //@ STEP 11 END MINIMIZATION
        tbegin = chrono::high_resolution_clock::now();
        endMinimization();


        #ifdef OPTIMOUT
        std::cout<<"End Minimization************"<<endl;
        cout<<"Time taken for lambda "<<lambdatime<<endl;
        std::cout << "----------------------------------------" << endl;
        #endif


        // Move the forces to the output if required.
        if(ptrOutForce) {
            *ptrOutForce = std::move(force);
        }

        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_end(tend - tbegin);
        CUDAcommon::tmin.endminimization+= elapsed_end.count();
        //@} END MINIMIZTION

#ifdef DETAILEDOUTPUT
        std::cout<<"printing beads & forces"<<endl;
        for(auto b:Bead::getBeads()){
            std::cout<<b->getId()<<" "<< b->coord <<" "
                    ""<<b->force <<endl;
        }
        std::cout<<"printed beads & forces"<<endl;
#endif


        return result;

    } // minimize(...)

    // Select the correct version to run.
    template< typename... Args >
    auto minimize(const ConjugateGradientParams& cgParams, Args&&... args) {
        switch(cgParams.descentSearchMethod) {
            case ConjugateGradientDescentSearch::steepest:
                return minimize<ConjugateGradientDescentSearch::steepest      >(cgParams, std::forward<Args>(args)...);
            case ConjugateGradientDescentSearch::fletcherRieves:
                return minimize<ConjugateGradientDescentSearch::fletcherRieves>(cgParams, std::forward<Args>(args)...);
            case ConjugateGradientDescentSearch::polakRibiere:
                return minimize<ConjugateGradientDescentSearch::polakRibiere  >(cgParams, std::forward<Args>(args)...);
            default:
                log::error("Unrecognized descent search method: {}", medyan::underlying(cgParams.descentSearchMethod));
                throw std::runtime_error("Unrecognized descent search method");
        }
    }

    template< typename EnergyFunc >
    void calculateEvsalpha(
        EnergyFunc&&  energyFunc,
        floatingpoint lambda,
        floatingpoint LAMBDAMAX,
        floatingpoint FDotFA,
        floatingpoint prevlambda,
        std::vector<floatingpoint>& coordBackup,
        std::vector<floatingpoint>& forceBackup
    ) {
        cout<<"Printing Evslambda information "<<endl;
        cout<<"chosen lambda "<<lambda<<" prev lambda "<<prevlambda<<endl;
    //	for(floatingpoint alpha=0.0;alpha<0.1;alpha=alpha+1e-4    ){
    //		cout<<alpha<<" ";
    //	}
    //	cout<<endl;
        floatingpoint energyzerolambda = energyFunc(coordBackup.data());
        cout<<"Energy zero lambda "<<energyzerolambda<<" "<<endl;
        cout<<"Lambda = [";
        for(floatingpoint alpha=LAMBDAMAX;alpha>=1e-4;alpha=alpha*LAMBDAREDUCE){
            cout<<alpha<<" ";
        }
        cout<<"];"<<endl;
        int count = 0;
        bool exityes = false;
        cout<<"Energy = [";
        for(floatingpoint alpha=LAMBDAMAX;alpha>=1e-4;alpha=alpha*LAMBDAREDUCE){
            //moveBeads
            const auto num = coordBackup.size();
            coordLineSearch.resize(num);
            for(size_t i = 0; i < num; ++i)
                coordLineSearch[i] = coordBackup[i] + alpha * forceBackup[i];
            floatingpoint energy = energyFunc(coordLineSearch.data());
            cout<<energy<<" ";
            if(count > 10)
                exityes = true;
            count++;
        }
        cout<<"];"<<endl;
        cout<<"Armijo = [";
        for(floatingpoint alpha=LAMBDAMAX;alpha>=1e-4;alpha=alpha*LAMBDAREDUCE){
            cout<<energyzerolambda-BACKTRACKSLOPE * alpha * FDotFA<<" ";
        }
        cout<<"];"<<endl;
        if(exityes)
            exit(EXIT_SUCCESS);
    //	exit(EXIT_FAILURE);

    }

	static ofstream _crosscheckdumpMechFile;

};

} // namespace medyan

#endif
