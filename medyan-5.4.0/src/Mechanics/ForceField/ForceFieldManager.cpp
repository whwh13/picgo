
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

#include "ForceFieldManager.h"

#include <numeric> // iota

#include "ForceFieldManagerCUDA.h"

#include <algorithm>
#include <utility>

#include "cross_check.h"
#include "Structure/Bead.h"
#include "Structure/Cylinder.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {






void ForceFieldManager::vectorizeAllForceFields(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {
#ifdef CUDAACCL
    // PT1 Generate single vector of energies from all FF and add them together.
    //@{
    CUDAcommon::cudavars.offset_E=0.0;
    //@}
#endif

    // Cache membrane indices.
    for(auto& m : ps->membranes) {
        medyan::cacheIndicesForFF(m.getMesh(), si);
    }

    for (auto &ff : forceFields) {
        ff->vectorize(si, conf);
    }

#ifdef CUDAACCL
    //reset offset
    if (streamF == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&streamF));
    int nint[1];
    nint[0] = CGMethod::N / 3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_nint, nint, sizeof(int),
                                        cudaMemcpyHostToDevice, streamF));
    CUDAcommon::handleerror(cudaMalloc((void **) &(CUDAcommon::cudavars.gpu_energyvec),
                                       CUDAcommon::cudavars.offset_E * sizeof(floatingpoint)));
    int THREADSPERBLOCK;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    THREADSPERBLOCK = prop.maxThreadsPerBlock;

    blocksnthreads.push_back(CGMethod::N / (3 * THREADSPERBLOCK) + 1);
    if (blocksnthreads[0] == 1) blocksnthreads.push_back(CGMethod::N / 3);
    else blocksnthreads.push_back(THREADSPERBLOCK);

    // PT2 Generate single vector of energies from all FF and add them together.
    //@{
//    std::cout<<"CUDA energy total nint "<<CUDAcommon::cudavars.offset_E<<endl;
    bntaddvec2.clear();
    bntaddvec2 = getaddred2bnt(CUDAcommon::cudavars.offset_E);
    CUDAcommon::handleerror(cudaMalloc((void **) &(CUDAcommon::cudavars.gpu_energyvec), bntaddvec2.at
            (0)*sizeof (floatingpoint)));
    vector<floatingpoint> zerovec(bntaddvec2.at(0));
    fill(zerovec.begin(),zerovec.begin()+bntaddvec2.at(0),0.0);
    CUDAcommon::handleerror(cudaMemcpyAsync(CUDAcommon::cudavars.gpu_energyvec, zerovec.data(),
                            bntaddvec2.at(0) * sizeof(floatingpoint), cudaMemcpyHostToDevice,streamF));
/*    CUDAcommon::handleerror(cudaMemsetAsync(CUDAcommon::cudavars.gpu_energyvec, 0,
                                            bntaddvec2.at(0) * sizeof(floatingpoint), streamF));*/

    params.clear();
    params.push_back(CUDAcommon::cudavars.offset_E);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), sizeof(int),
                                       cudaMemcpyHostToDevice, streamF));
    //@}
#endif
}


floatingpoint ForceFieldManager::computeEnergy(floatingpoint *coord, bool verbose) const {
    chrono::high_resolution_clock::time_point tbegin, tend;
    floatingpoint energy = 0.0;
#ifdef CUDAACCL
#ifdef CUDA_INDIVIDUAL_ESUM
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (floatingpoint)));
    CUDAcommon::handleerror(cudaMemset(gpu_Uvec, 0.0, sizeof (floatingpoint)));
#else
    floatingpoint *gpu_Uvec = CUDAcommon::getCUDAvars().gpu_energy;
    /*floatingpoint *gpu_Uvec;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Uvec, sizeof (floatingpoint)));
    CUDAcommon::handleerror(cudaMemsetAsync(gpu_Uvec, 0, sizeof (floatingpoint),streamF));*/
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaMemset(CUDAcommon::cudavars.gpu_energyvec, 0, bntaddvec2.at(0) * sizeof
            (floatingpoint)));
#endif
/*    auto gU_tot = CUDAcommon::getCUDAvars().gpu_energy;
    setenergytozero << < 1, 1, 0, streamF >> > (gU_tot);*/
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#endif
#ifdef SERIAL_CUDACROSSCHECK
    CUDAcommon::handleerror(cudaDeviceSynchronize());
    floatingpoint cuda_lambda[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));

#endif

    // Update dependent variables using functions defined in force fields.
    computeDependentCoordinates(coord);
    // Compute membrane geometry
    for(auto& m : ps->membranes) {
        updateGeometryValue(m.getMesh(), coord, surfaceGeometrySettings);
    }

    short count = 0;
    CUDAcommon::tmin.computeenergycalls++;
    #ifdef TRACKDIDNOTMINIMIZE
    SysParams::Mininimization().tempEnergyvec.clear();
	#endif
    for (auto &ff : forceFields) {
        tbegin = chrono::high_resolution_clock::now();
        auto tempEnergy = ff->computeEnergy(coord);
        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
//        cout<<ff->getName()<<" "<<tempEnergy<<"pN.nm"<<" ";

#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

        if (verbose) cout << ff->getName() << " energy = " << tempEnergy << endl;
        //if energy is infinity, exit with infinity.
        if (tempEnergy <= -1) {

            cout << "Energy = " << tempEnergy << endl;

            cout
                    << "Energy of system became infinite. Try adjusting minimization parameters."
                    << endl;
            cout << "The culprit was ... " << ff->getName() << endl;

            ff->whoIsCulprit();
            return numeric_limits<floatingpoint>::infinity();
        }
        else energy += tempEnergy;
#ifdef SERIAL_CUDACROSSCHECK
        cudaDeviceSynchronize();
        resetfloatingpointvariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
        addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint)
                                                                   , streamF>>>
        (CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                gpu_Uvec);
        floatingpoint cuda_energyvec[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
        std::cout<<ff->getName()<<" Energy. CUDA "<<cuda_energyvec[0]<<" SERL "
                ""<<energy<<endl;
#endif
    }
#ifdef CUDAACCL
//    std::cout<<"Total nint "<<bntaddvec2.at(0)<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    //Synchronize streams
    for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "computeEnergy",
                                    "ForceFieldManager.cu");
        }
    resetfloatingpointvariableCUDA<<<1,1,0, streamF>>>(gpu_Uvec);
    addvectorred3<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint)
                    , streamF>>>(CUDAcommon::cudavars.gpu_energyvec, gpu_params,
                        gpu_Uvec);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));
#ifdef DETAILEDOUTPUT_ENERGY
    cudaDeviceSynchronize();
    floatingpoint cuda_energyvec[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energyvec, gpu_Uvec, sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"vector energy addition CUDA "<<cuda_energyvec[0]<<" SERL "<<energy<<endl;
#endif
//    CUDAcommon::handleerror(cudaFree(CUDAcommon::cudavars.gpu_energyvec));
#endif
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif

    #ifdef TRACKDIDNOTMINIMIZE
    SysParams::Mininimization().Energyvec.push_back(SysParams::Mininimization().tempEnergyvec);
    #endif
    return energy;
    
}

EnergyReport ForceFieldManager::computeIndividualEnergies(floatingpoint* coord) const {
    // Compute membrane geometry
    for(auto& m : ps->membranes) {
        updateGeometryValue(m.getMesh(), coord, surfaceGeometrySettings);
    }

    EnergyReport result;
    result.total = 0.0;

    // Update dependent variables using functions defined in force fields.
    computeDependentCoordinates(coord);

    for (auto &ff : forceFields) {
        auto tempEnergy = ff->computeEnergy(coord);
        result.individual.push_back({ ff->getName(), tempEnergy });
        result.total += tempEnergy;
    }
    return result;
}

EnergyReport ForceFieldManager::computeEnergyHRMD(floatingpoint *coord) const {
    CUDAcommon::tmin.computeenerycallszero++;
    EnergyReport result = computeIndividualEnergies(coord);
    for(auto& eachResult : result.individual) {
        eachResult.energy /= kT;
    }
    result.total /= kT;
    return result;
}




void ForceFieldManager::computeForces(floatingpoint *coord, floatingpoint* force, int numVar) {
    // Reset forces to zero.
    std::fill(force, force + numVar, 0.0);

#ifdef CUDAACCL
    CUDAvars cvars = CUDAcommon::getCUDAvars();
    if (cross_checkclass::Aux)
        resetForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, streamF >> >
                                                                      (cvars.gpu_forceAux, gpu_nint);
    else
        resetForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, streamF >> >
                                                                      (cvars.gpu_force, gpu_nint);
    CUDAcommon::handleerror(cudaStreamSynchronize(streamF));

    CUDAcommon::handleerror(cudaGetLastError(), "resetForcesCUDA", "ForceFieldManager.cu");
#endif

    // Update dependent variables using functions defined in force fields.
    computeDependentCoordinates(coord);
    // compute membrane geometry with derivatives
    for(auto& m : ps->membranes) {
        updateGeometryValueWithDerivative(m.getMesh(), coord, surfaceGeometrySettings);
    }

    //recompute
//    floatingpoint *F_i = new floatingpoint[CGMethod::N];
    short count = 0;
    CUDAcommon::tmin.computeforcescalls++;
    for (auto &ff : forceFields) {
        tbegin = chrono::high_resolution_clock::now();
        ff->computeForces(coord, force);
        tend = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_energy(tend - tbegin);
        if(CUDAcommon::tmin.individualforces.size() == forceFields.size())
            CUDAcommon::tmin.individualforces[count]+= elapsed_energy.count();
        else
            CUDAcommon::tmin.individualforces.push_back(elapsed_energy.count());
        count++;

#ifdef ALLSYNC
        cudaDeviceSynchronize();
#endif

    }

    // After force computation, propagate all forces accumulated on dependent variables onto the independent variables using the chain rule.
    propagateDependentForces(coord, force);

}

void ForceFieldManager::computeLoadForces() const {

    //reset
    for (auto b: Bead::getBeads()) {
        std::fill(b->loadForcesM.begin(), b->loadForcesM.end(), 0.0);
        std::fill(b->loadForcesP.begin(), b->loadForcesP.end(), 0.0);
//        b->loadForcesP.clear();
//        b->loadForcesM.clear();
    }

    for (auto &f : forceFields)
        f->computeLoadForces();

    //reset lfi as well
    for (auto b: Bead::getBeads()) {
        b->lfip = 0;
        b->lfim = 0;
    }
}
void ForceFieldManager::computeLoadForce(Cylinder* c, ForceField::LoadForceEnd end) const {
    auto  b          = (end == ForceField::LoadForceEnd::Plus ? c->getSecondBead() : c->getFirstBead());
    auto& loadForces = (end == ForceField::LoadForceEnd::Plus ? b->loadForcesP     : b->loadForcesM   );
    auto& lfi        = (end == ForceField::LoadForceEnd::Plus ? b->lfip            : b->lfim          );

    // reset
    std::fill(loadForces.begin(), loadForces.end(), 0.0);

    for(auto& f : forceFields) f->computeLoadForce(*ps, c, end);

    // reset lfi
    lfi = 0;
}


#ifdef CUDAACCL

void ForceFieldManager::CUDAcopyForces(cudaStream_t stream, floatingpoint *fprev, floatingpoint *f) {


//    CUDAcommon::handleerror(cudaFree(CUDAcommon::getCUDAvars().gpu_forceAux));
//    floatingpoint* gpu_forceAux;
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, CGMethod::N * sizeof(floatingpoint)));
//    CUDAvars cvars=CUDAcommon::getCUDAvars();
//    cvars.gpu_forceAux=gpu_forceAux;
//    CUDAcommon::cudavars=cvars;

//    std::cout<<"Copyforces Number of Blocks: "<<blocksnthreads[0]<<endl;
//    std::cout<<"Threads per block: "<<blocksnthreads[1]<<endl;
    copyForcesCUDA << < blocksnthreads[0], blocksnthreads[1], 0, stream >> >
                                                                 (f, fprev, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "copyForcesCUDA", "ForceFieldManager.cu");
}

#endif

void ForceFieldManager::assignallforcemags() {

    for (auto &ff : forceFields)
        ff->assignforcemags();
}


void ForceFieldManager::computeHessian(const std::vector<floatingpoint>& allCoord, int total_DOF, float delta) {
    // store the minimization time and initialize the matrix
    tauVector.push_back(tau());

    chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();

    vector<vector<floatingpoint> > hessianMatrix(total_DOF, vector<floatingpoint>(total_DOF));

    vector<Triplet> tripletList;


    // loop through all the coordinates
    for (auto i = 0; i < total_DOF; i++) {
        cout<<"i "<<i<<" total_DOF "<<total_DOF<<endl;

        // create new vectors for the foorces and coordinates
        vector<floatingpoint> forces_copy_p(allCoord.size());
        vector<floatingpoint> coord_copy_p = allCoord;

        vector<floatingpoint> forces_copy_m(allCoord.size());
        vector<floatingpoint> coord_copy_m = allCoord;

        // perturb the coordinate i
        coord_copy_p[i] += delta;
        coord_copy_m[i] -= delta;

        // calculate the new forces based on perturbation
        computeForces(coord_copy_p.data(), forces_copy_p.data(), forces_copy_p.size());
        computeForces(coord_copy_m.data(), forces_copy_m.data(), forces_copy_p.size());

        for (auto j = 0; j < total_DOF; j++) {

            // store the derivative of the force on each coordinate j
            float h_i_j = -(forces_copy_p[j] - forces_copy_m[j]) / (2 * delta);
            hessianMatrix[i][j] = h_i_j;
            tripletList.push_back(Triplet(i, j, h_i_j));

        }

    }


    // store the full matrix in list
    hessianVector.push_back(hessianMatrix);

    if(SysParams::Mechanics().eigenTracking) {
        // create symmetrized sparse matrix object
        Eigen::SparseMatrix<double> hessMat(total_DOF, total_DOF), hessMatSym;
        hessMat.setFromTriplets(tripletList.begin(), tripletList.end());
        hessMatSym = 0.5 * (Eigen::SparseMatrix<double>(hessMat.transpose()) + hessMat);

        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_vecmat(t1 - t0);

        bool denseEstimationBool = SysParams::Mechanics().denseEstimationBool;

        if (denseEstimationBool) {

            Eigen::MatrixXd denseHessMatSym;
            denseHessMatSym = Eigen::MatrixXd(hessMatSym);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> denseESolver(denseHessMatSym);
            evalues = denseESolver.eigenvalues().real();
            evectors = denseESolver.eigenvectors().real();

        } else {

            Spectra::SparseSymShiftSolve<double> op(hessMatSym);
            //Spectra::SparseSymMatProd<double> op(hessMatSym);
            int numEigs = total_DOF - 1;
            Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op,
                                                                                                                   numEigs,
                                                                                                                   total_DOF,
                                                                                                                   10000);
            //Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double>> eigs(&op, numEigs, numEigs+1);
            /*
            if(evectors.size()!=0){
                //const Eigen::Matrix<double, Eigen::Dynamic, 1> init_vec = evectors.col(0).real();
                const Eigen::Matrix<double, Eigen::Dynamic, 1> init_vec = evectors.real().rowwise().sum();
                if(init_vec.size() == total_DOF){
                    const double * arg = init_vec.data();
                    eigs.init(arg);
                    //eigs.init();
                }else{
                    eigs.init();
                };
            }else{
                eigs.init();
            }*/

            eigs.init();
            int nconv = eigs.compute(Spectra::SortRule::SmallestAlge);
            evalues = eigs.eigenvalues();
            //columns of evectors matrix are the normalized eigenvectors
            evectors = eigs.eigenvectors(numEigs);
        };

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_veceigs(t2 - t1);
        Eigen::VectorXcd IPRI(evectors.cols());
        Eigen::VectorXcd IPRII(evectors.cols());

        // compute participation ratios
        for (auto i = 0; i < evectors.cols(); i++) {
            floatingpoint RI = 0.0;
            Eigen::VectorXd col = evectors.col(i).cwiseAbs2();
            RI = pow(col.norm(), 2);
            floatingpoint RII = 0.0;
            for (auto j = 0; j < evectors.rows() / 3; j++) {
                floatingpoint temp = 0.0;
                for (auto k = 0; k < 3; k++) {
                    temp += pow(evectors(3 * j + k, i).real(), 2);
                }
                RII += pow(temp, 2);
            }
            IPRI(i) = RI;
            IPRII(i) = RII;
        }


        chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_vecPR(t3 - t2);

        cout << "DOF is " << total_DOF << endl;
        std::cout << "Matrix time " << elapsed_vecmat.count() << endl;
        std::cout << "Compute time " << elapsed_veceigs.count() << endl;
        std::cout << "PR time " << elapsed_vecPR.count() << endl;



        // store the eigenvalues in list
        IPRIVector.push_back(IPRI);
        IPRIIVector.push_back(IPRII);
        evaluesVector.push_back(evalues);
    }

}


void ForceFieldManager::setCurrBeadMap(const FFCoordinateStartingIndex& si) {
    prevBeadMap = std::move(currBeadMap);
    currBeadMap.clear();
    for(auto b:Bead::getBeads()){
        currBeadMap[b] = medyan::findBeadCoordIndex(*b, si);
    }
}

void ForceFieldManager::computeProjections(const FFCoordinateStartingIndex& si, const std::vector<floatingpoint>& currCoords) {
    
    // set displacement vector to zeros
    Eigen::VectorXcd disp(evectors.rows());
    for(auto i =0; i<evectors.rows(); i++){
        disp(i) = 0.0;
    }
    
    // set the current beads
    setCurrBeadMap(si);
    

    // loop through the current map
    for(const std::pair<Bead*, int>& c : currBeadMap){
        
        // find the key in the previous map
        auto p = prevBeadMap.find(c.first);
        
        // if it's in the old map get the displacements
        if(p != prevBeadMap.end()){
            int currInd = c.second;
            int prevInd = p->second;
            floatingpoint cx = currCoords[currInd];
            floatingpoint px = prevCoords[prevInd];
            floatingpoint cy = currCoords[currInd + 1];
            floatingpoint py = prevCoords[prevInd + 1];
            floatingpoint cz = currCoords[currInd + 2];
            floatingpoint pz = prevCoords[prevInd + 2];
            disp(prevInd) = cx - px;
            disp(prevInd + 1) = cy - py;
            disp(prevInd + 2) = cz - pz;
        };
    };
    
    // normalize the displacement vector
    if(!disp.isZero(0)) {
        disp = disp.normalized();
    }
    
    // store the projections in a vector
    Eigen::VectorXcd proj(evectors.cols());
    for(auto i = 0; i < evectors.cols(); i++){
        proj(i) = disp.dot(evectors.col(i));
    }
    
    // test that you can reconstruct disp from the projections
    /*
    Eigen::VectorXcd testDisp(evectors.rows());
    for(auto i =0; i<evectors.rows(); i++){
        testDisp(i) = 0.0;
    }
    for(auto i = 0; i < evectors.cols(); i++){
        testDisp += proj(i) * evectors.col(i);
    }
    for(auto i =0; i<evectors.rows(); i++){
        cout<<testDisp(i).real()<<" "<<disp(i).real()<<endl;
    }
    */
    
    projectionsVector.push_back(proj);
    
    // need to reset the prevCoords with current coords
    prevCoords = currCoords;
   
    
}

} // namespace medyan
