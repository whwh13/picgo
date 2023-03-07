
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
#include "SubSystem.h"
#include "BoundaryElement.h"
#include "CompartmentGrid.h"
#include "BindingManager.h"
#include "BindingManagerCUDA.h"
#include "MathFunctions.h"
#include "BoundaryElement.h"
#include "BoundaryElementImpl.h"
#include <vector>
#include "Cylinder.h"
#ifdef SIMDBINDINGEARCH
#include "Util/DistModule/dist_driver.h"
#include "Util/DistModule/dist_coords.h"
#include "Util/DistModule/dist_common.h"
#endif


namespace medyan {

using namespace mathfunc;
void SubSystem::resetNeighborLists() {

#ifdef CUDAACCL_NL
    coord = new floatingpoint[CGMethod::N];
                coord_com = new floatingpoint[3 * Cylinder::getCylinders().size()];
                beadSet = new int[2 * Cylinder::getCylinders().size()];
                cylID = new int[Cylinder::getCylinders().size()];
                filID = new int[Cylinder::getCylinders().size()];
                filType = new int[Cylinder::getCylinders().size()];
                cmpID = new unsigned int[Cylinder::getCylinders().size()];
                fvecpos = new int[Cylinder::getCylinders().size()];

                if(SysParams::Chemistry().numFilaments > 1) {
                    cout << "CUDA NL cannot handle more than one type of filaments." << endl;
                    exit(EXIT_FAILURE);
                }
                int numBindingSites = SysParams::Chemistry().bindingSites[0].size();
                if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                    cmon_state_brancher = new int[ numBindingSites * Cylinder::getCylinders().size()];
                if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                    cmon_state_linker = new int[numBindingSites * Cylinder::getCylinders().size()];
                if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                    cmon_state_motor = new int[numBindingSites * Cylinder::getCylinders().size()];

                int i = 0; //int cID = 0;
                for(auto b: Bead::getBeads()) {
                    //flatten indices
                    int index = 3 * i;
                    coord[index] = b->vcoordinate()[0];
                    coord[index + 1] = b->vcoordinate()[1];
                    coord[index + 2] = b->vcoordinate()[2];
                    i++;
                }
                i = 0;

                 for(auto c:Cylinder::getCylinders()){
                            //flatten indices
                            int index = 3 * i;
                            coord_com[index] = c->coordinate[0];
                            coord_com[index + 1] = c->coordinate[1];
                            coord_com[index + 2] = c->coordinate[2];

                        beadSet[2 * i] = c->getFirstBead()->getStableIndex();
                        beadSet[2 * i + 1] = c->getSecondBead()->getStableIndex();
                        cylID[i] = c->getId();
                        c->_dcIndex = i;
                        fvecpos[i] = c->getPosition();
                        auto fil = dynamic_cast<Filament*>(c->getParent());
                        filID[i] =  fil->getId();
                        cmpID[i] = GController::getCompartmentID(c->getCompartment()->coordinates());
                        filType[i] = fil->getType();
        //                cylstate[i] = c->isFullLength();
                        int j = 0;
                        for(auto it2 = SysParams::Chemistry().bindingSites[fil->getType()].begin();
                            it2 != SysParams::Chemistry().bindingSites[fil->getType()].end(); it2++) {
                            if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                                cmon_state_brancher[numBindingSites * i + j ] = c->getCCylinder()->getCMonomer(*it2)
                                        ->speciesBound(SysParams::Chemistry().brancherBoundIndex[fil->getType()])->getN();
                            if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                                cmon_state_linker[numBindingSites * i + j ] = c->getCCylinder()->getCMonomer(*it2)
                                        ->speciesBound(SysParams::Chemistry().linkerBoundIndex[fil->getType()])->getN();
                            if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                                cmon_state_motor[numBindingSites * i + j ] = c->getCCylinder()->getCMonomer(*it2)
                                        ->speciesBound(SysParams::Chemistry().motorBoundIndex[fil->getType()])->getN();
                            j++;
                        }
                        i++;
                    }
        //        }//Compartment
                //CUDAMALLOC
                //@{
        //        size_t free, total;
        //        CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
        //        fprintf(stdout,"\t### Available VRAM : %g Mo/ %g Mo(total)\n\n",
        //                free/1e6, total/1e6);
        //
        //        cudaFree(0);
        //
        //        CUDAcommon::handleerror(cudaMemGetInfo(&free, &total));
        //        fprintf(stdout,"\t### Available VRAM : %g Mo/ %g Mo(total)\n\n",
        //                free/1e6, total/1e6);
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_coord, CGMethod::N * sizeof(floatingpoint)),"cuda data "
                                        "transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_coord_com, 3 * Cylinder::getCylinders().size() * sizeof
                                        (floatingpoint)),"cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, 2 * Cylinder::getCylinders().size() * sizeof
                                        (int)), "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylID, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_fvecpos, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_filID, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_filType, Cylinder::getCylinders().size() * sizeof(int)),
                                        "cuda data transfer", " SubSystem.h");
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmpID, Cylinder::getCylinders().size() * sizeof(unsigned
                                        int)), "cuda data transfer", " SubSystem.h");
        //        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylstate, Cylinder::getCylinders().size() * sizeof(int)),
        //                                "cuda data transfer", " SubSystem.h");

                if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmon_state_brancher, numBindingSites *
                        Cylinder::getCylinders().size() * sizeof(int)), "cuda data transfer", " SubSystem.h");
                if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmon_state_linker, numBindingSites *
                                            Cylinder::getCylinders().size() * sizeof(int)), "cuda data transfer", " SubSystem.h");
                if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cmon_state_motor, numBindingSites *
                                            Cylinder::getCylinders().size() * sizeof(int)), "cuda data transfer", " SubSystem.h");

                //@}
                //CUDAMEMCPY
                //@{
                CUDAcommon::handleerror(cudaMemcpy(gpu_coord, coord, CGMethod::N *sizeof(floatingpoint), cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_coord_com, coord_com, 3 * Cylinder::getCylinders().size() *sizeof
                                                   (floatingpoint), cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet, 2 * Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_cylID, cylID, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_fvecpos, fvecpos, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_filID, filID, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_filType, filType, Cylinder::getCylinders().size() *sizeof(int),
                                                   cudaMemcpyHostToDevice));
                CUDAcommon::handleerror(cudaMemcpy(gpu_cmpID, cmpID, Cylinder::getCylinders().size() *sizeof(unsigned int),
                                                   cudaMemcpyHostToDevice));
        //        CUDAcommon::handleerror(cudaMemcpy(gpu_cylstate, cylstate, Cylinder::getCylinders().size() *sizeof(int),
        //                                           cudaMemcpyHostToDevice));
                if(SysParams::Chemistry().numBrancherSpecies[0] > 0)
                CUDAcommon::handleerror(cudaMemcpy(gpu_cmon_state_brancher, cmon_state_brancher, numBindingSites *
                                        Cylinder::getCylinders().size() * sizeof(int), cudaMemcpyHostToDevice));
                if(SysParams::Chemistry().numLinkerSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMemcpy(gpu_cmon_state_linker, cmon_state_linker, numBindingSites *
                                        Cylinder::getCylinders().size() *sizeof(int), cudaMemcpyHostToDevice));
                if(SysParams::Chemistry().numMotorSpecies[0] > 0)
                    CUDAcommon::handleerror(cudaMemcpy(gpu_cmon_state_motor, cmon_state_motor, numBindingSites *
                                        Cylinder::getCylinders().size() *sizeof(int), cudaMemcpyHostToDevice));

                CylCylNLvars cylcylnlvars;
                cylcylnlvars.gpu_coord = gpu_coord;
                cylcylnlvars.gpu_coord_com = gpu_coord_com;
                cylcylnlvars.gpu_beadSet = gpu_beadSet;
                cylcylnlvars.gpu_cylID = gpu_cylID;
                cylcylnlvars.gpu_fvecpos = gpu_fvecpos;
                cylcylnlvars.gpu_filID = gpu_filID;
                cylcylnlvars.gpu_filType = gpu_filType;
                cylcylnlvars.gpu_cmpID = gpu_cmpID;
        //        cylcylnlvars.gpu_cylstate = gpu_cylstate;
                cylcylnlvars.gpu_cmon_state_brancher = gpu_cmon_state_brancher;
                cylcylnlvars.gpu_cmon_state_linker = gpu_cmon_state_linker;
                cylcylnlvars.gpu_cmon_state_motor = gpu_cmon_state_motor;
        //        cylcylnlvars.gpu_cylvecpospercmp = gpu_cylvecpospercmp;

                CUDAcommon::cylcylnlvars = cylcylnlvars;
#endif
    chrono::high_resolution_clock::time_point mins, mine;
    mins = chrono::high_resolution_clock::now();

    #if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        _HneighborList->reset();
        mine= chrono::high_resolution_clock::now();
        #ifdef OPTIMOUT
            chrono::duration<floatingpoint> elapsed_H(mine - mins);
            std::cout<<"H NLSTEN reset time "<<elapsed_H.count()<<endl;
            mins = chrono::high_resolution_clock::now();
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_B(mine - mins);
            std::cout<<"H NLSTEN B reset time "<<elapsed_B.count()<<endl;
        #endif
    #endif
    for (auto nl: _neighborLists)
        nl->reset();

    // Other neighbor lists or neighbor list primitives (such as cell lists, BVHs).
    //----------------------------------

    // Bubble neighbor lists.
    if(opBoundaryBubbleNL.has_value()) {
        opBoundaryBubbleNL->reset(*this);
    }
    if(opBubbleBubbleNL.has_value()) {
        opBubbleBubbleNL->reset(*this);
    }
    if(opBubbleBeadNL.has_value()) {
        opBubbleBeadNL->reset(*this);
    }

    // Meshless vertices.
    {
        meshlessSpinVertexCellList.clearElements();
        // Add all vertices to cell lists.
        for(auto& v : meshlessSpinVertices) {
            v.meshlessSpinVertexCellListIndex = meshlessSpinVertexCellList.add(v.sysIndex.value, v.coord);
        }
    }
}
void SubSystem::updateBindingManagers(medyan::SimulConfig& sc) {
#ifdef OPTIMOUT
	chrono::high_resolution_clock::time_point mins, mine, minsinit, mineinit,
			startonetimecost, endonetimecost;
	mins = chrono::high_resolution_clock::now();

#endif
#ifdef CUDAACCL_NL
	if(SysParams::Chemistry().numFilaments > 1) {
		cout << "CUDA Binding Manager cannot handle more than one type of filaments." << endl;
		exit(EXIT_FAILURE);
	}
	initializebindingsitesearchCUDA();

	if(CUDAcommon::getCUDAvars().conservestreams)
		numbindmgrs = 0;
	//Calculate binding sites in CUDA
	Compartment* C0 = _compartmentGrid->getCompartments()[0].get();
	for(auto &manager : C0->getFilamentBindingManagers()) {

		LinkerBindingManager *lManager;
		MotorBindingManager *mManager;
		BranchingManager *bManager;
		auto cylcylnlvars = CUDAcommon::getCylCylNLvars();
		auto coord = cylcylnlvars.gpu_coord;
		auto beadSet = cylcylnlvars.gpu_beadSet;
		auto cylID = cylcylnlvars.gpu_cylID;
		auto filType = cylcylnlvars.gpu_filType;
		auto filID = cylcylnlvars.gpu_filID;
		int *cmpID = cylcylnlvars.gpu_cmpID;
		//Linker
		if ((lManager = dynamic_cast<LinkerBindingManager *>(manager.get()))) {
			//calculate all binding Sites.
			getallpossiblelinkerbindingsitesCUDA(lManager, cylcylnlvars.gpu_cmon_state_linker);
		}
		//Motor
		else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
			//calculate all binding Sites.
			getallpossiblemotorbindingsitesCUDA(mManager, cylcylnlvars
												.gpu_cmon_state_motor);
		}
		//Brancher
		else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {
			//calculate all binding Sites.
			getallpossiblebrancherbindingsitesCUDA(bManager, cylcylnlvars
												   .gpu_cmon_state_brancher);

		}
	}

	//Free vars
	terminatebindingsitesearchCUDA();
	//Assign to respective possible bindings.
	assigntorespectivebindingmanagersCUDA();

	//    for(auto gpb:gpu_possibleBindings_vec)
	//        CUDAcommon::handleerror(cudaFree(gpb),"cudaFree","SubSystem.cu");
	//    for(auto pb:possibleBindings_vec)
	//        CUDAcommon::handleerror(cudaFreeHost(pb), "cudaFree", "SubSystem.cu");
	//    for(auto np:numpairs_vec)
	//        CUDAcommon::handleerror(cudaFreeHost(np),"cudaFree","SubSystem.cu");

	//cudaFree
	endresetCUDA();
#endif
	vectorizeCylinder(sc);
	if(CROSSCHECK_SWITCH)
	    HybridNeighborList::_crosscheckdumpFileNL<<"vectorized Cylinder"<<endl;
	//Version1
	#ifdef NLORIGINAL
	for (auto& C : compartmentGrid->getCompartments()){
		for(auto &manager : C->getFilamentBindingManagers()){
			manager->updateAllPossibleBindings();
		}
	}
	#endif
	//Version2
	#ifdef NLSTENCILLIST
	#if !defined(HYBRID_NLSTENCILLIST) && !defined(SIMDBINDINGSEARCH)
	for (auto& C : compartmentGrid->getCompartments()){
		for(auto &manager : C->getFilamentBindingManagers()){
			manager->updateAllPossibleBindingsstencil();
		}
	}
	#endif
	#endif
	//Version3
	#ifdef HYBRID_NLSTENCILLIST
	for (auto& C : compartmentGrid->getCompartments()) {
		C->getHybridBindingSearchManager()->updateAllPossibleBindingsstencilHYBD();
		for(auto &manager : C->getBranchingManagers()) {
			manager->updateAllPossibleBindingsstencil();
		}
	}
	//UpdateAllBindingReactions
	for (auto& C : compartmentGrid->getCompartments()) {
		C->getHybridBindingSearchManager()->updateAllBindingReactions();
	}
	#endif

	//SIMD cylinder update
#ifdef SIMDBINDINGSEARCH
#ifdef OPTIMOUT
	startonetimecost = chrono::high_resolution_clock::now();
#endif
	if(!initialize) {
        compartmentGrid->getCompartments()[0]->getHybridBindingSearchManager()
			->initializeSIMDvars();
		initialize = true;
	}
#ifdef OPTIMOUT
	endonetimecost = chrono::high_resolution_clock::now();
#endif
	//Generate binding site coordinates in each compartment and seggregate them into
	// different spatial sub-sections.
	for(auto& C : compartmentGrid->getCompartments()) {
		C->SIMDcoordinates4linkersearch_section(1);
		C->SIMDcoordinates4motorsearch_section(1);
	}
    if(CROSSCHECK_SWITCH) {
        HybridNeighborList::_crosscheckdumpFileNL<<"Generated SIMD input files"<<endl;
    }

	//Empty the existing binding pair list
	for (auto& C : compartmentGrid->getCompartments())
		C->getHybridBindingSearchManager()->resetpossibleBindings();
    if(CROSSCHECK_SWITCH) {
        HybridNeighborList::_crosscheckdumpFileNL<<"Completed reset of binding pair maps"<<endl;
    }

	minsSIMD = chrono::high_resolution_clock::now();
	medyan::HybridBindingSearchManager::findtimeV3 = 0.0;
	medyan::HybridBindingSearchManager::SIMDV3appendtime = 0.0;
	// Update binding sites in SIMD
	for (auto& C : compartmentGrid->getCompartments()) {
		C->getHybridBindingSearchManager()->updateAllPossibleBindingsstencilSIMDV3();
        if(CROSSCHECK_SWITCH) {
            HybridNeighborList::_crosscheckdumpFileNL
			<<"L/M Update binding pair map in Cmp "<<C->getId()<<endl;
        }
		for(auto &manager : C->getBranchingManagers()) {
				manager->updateAllPossibleBindingsstencil();
            if(CROSSCHECK_SWITCH) {
                HybridNeighborList::_crosscheckdumpFileNL
				<<"B Update binding pair map in Cmp "<<C->getId()<<endl;
            }
		}
	}
	//UpdateAllBindingReactions
	for (auto& C : compartmentGrid->getCompartments()) {
//        cout<<"Cmp ID "<<C->getID()<<endl;
		C->getHybridBindingSearchManager()->updateAllBindingReactions();
        if(CROSSCHECK_SWITCH) {
            HybridNeighborList::_crosscheckdumpFileNL
		    <<"Updated binding reaction rates in Cmp "<<C->getId()<<endl;
        }
    }

#ifdef OPTIMOUT
	mineSIMD = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_onetimecost(endonetimecost -
	startonetimecost);
	chrono::duration<floatingpoint> elapsed_runSIMDV3(mineSIMD - minsSIMD);
	cout << "SIMDV3 time " << elapsed_runSIMDV3.count() << endl;
	cout << "findV3 time " << medyan::HybridBindingSearchManager::findtimeV3 << endl;
	cout << "Append time " << medyan::HybridBindingSearchManager::SIMDV3appendtime << endl;
#endif
#endif
#ifdef OPTIMOUT
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_orig(mine - mins);
    std::cout<<"BMgr-update time "<<elapsed_orig.count() - elapsed_onetimecost.count()<<endl;
#endif

//free memory
	SysParams::MParams.speciesboundvec.clear();
	#ifdef SIMDBINDINGSEARCH
	for(auto& C : compartmentGrid->getCompartments()) {
		C->deallocateSIMDcoordinates();
	}
	#endif
}

void SubSystem::vectorizeCylinder(medyan::SimulConfig& sc) {
    //vectorize
    SysParams::MParams.speciesboundvec.clear();
//    int cidx = 0;
    vector<int> ncylvec(sc.chemParams.numFilaments);// Number of cylinders
    // corresponding to each filament type.
    auto cylvec = Cylinder::getCylinders();
//    int ncyl = cylvec.size();
    if(cylsqmagnitudevector != nullptr){
        delete [] cylsqmagnitudevector;
    };
    cylsqmagnitudevector = new floatingpoint[Cylinder::rawNumStableElements()];
    unsigned long maxbindingsitespercyl = 0;
    for(auto ftype = 0; ftype < sc.chemParams.numFilaments; ftype++) {
        maxbindingsitespercyl = max<size_t>(maxbindingsitespercyl,SysParams::Chemistry()
                .bindingSites[ftype].size());
    }
    long vectorsize = maxbindingsitespercyl * Cylinder::rawNumStableElements();
    vector<bool> branchspeciesbound(vectorsize);
    vector<bool> linkerspeciesbound(vectorsize);
    vector<bool> motorspeciesbound(vectorsize);//stores species bound corresponding to each
    // cylinder.

    //set the size of each species bound vector
    fill(branchspeciesbound.begin(),branchspeciesbound.begin()+vectorsize, 0);
    fill(linkerspeciesbound.begin(),linkerspeciesbound.begin()+vectorsize, 0);
    fill(motorspeciesbound.begin(),motorspeciesbound.begin()+vectorsize, 0);

    //fill with appropriate values.
    for (auto cyl: cylvec) {
        auto _filamentType = cyl->getType();
        auto x1 = cyl->getFirstBead()->vcoordinate();
        auto x2 = cyl->getSecondBead()->vcoordinate();
        vector<floatingpoint> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
        cylsqmagnitudevector[cyl->getStableIndex()] = sqmagnitude(X1X2);

//        fil->printSelf();
        auto cc = cyl->getCCylinder();
        int idx = 0;
        for (auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
             it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {

            branchspeciesbound[maxbindingsitespercyl * cyl->getStableIndex() + idx] =
                    (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN());
            linkerspeciesbound[maxbindingsitespercyl * cyl->getStableIndex() + idx] =
                    (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN());
            motorspeciesbound[maxbindingsitespercyl * cyl->getStableIndex() + idx] =
                    (cc->getCMonomer(*it1)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[_filamentType])->getN());
            idx++;
        }
    }
    //@}

    SysParams::MParams.speciesboundvec.push_back(branchspeciesbound);
    SysParams::MParams.speciesboundvec.push_back(linkerspeciesbound);
    SysParams::MParams.speciesboundvec.push_back(motorspeciesbound);
    sc.chemParams.maxbindingsitespercylinder = maxbindingsitespercyl;
    SysParams::MParams.cylsqmagnitudevector = cylsqmagnitudevector;
    SysParams::MParams.ncylvec = ncylvec;

	//Enter filament first cylinder position
	for(auto fil:Filament::getFilaments()) {
		int filamentfirstentry = fil->getMinusEndCylinder()->getPosition();
		for(auto cyl:fil->getCylinderVector()){
			auto cindex = cyl->getStableIndex();
			Cylinder::getDbData()[cindex].filamentFirstEntry = filamentfirstentry;
		}
	}

    // Section for backward compatibility only.
    SysParams::CParams.maxbindingsitespercylinder = sc.chemParams.maxbindingsitespercylinder;
}

#ifdef CUDAACCL_NL
void SubSystem::initializebindingsitesearchCUDA() {
    //@{ 1. InitializeBSsearch
    //Reset variables
    numpairs_vec.clear();
    possibleBindings_vec.clear();
    gpu_possibleBindings_vec.clear();
//   auto x = CMonomer::_numBSpecies;
//    auto var = SysParams::Chemistry().bmanagerdistances;
    //Malloc params

    //Copy necessary cylinder data to GPU memory
    //@{
    //        if (gpu_params == NULL) {
    int params[3];
    params[0] = SysParams::Chemistry().numBindingSites[0];//filType dependant
    params[1] = 0;//filType dependant
    params[2] = SysParams::Geometry().cylinderNumMon[0];//filType dependant.
    params[3] = Cylinder::getCylinders().size();
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 4 * sizeof(int)),
                            "cuda data transfer", "SubSystem.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params, 4 * sizeof(int),
                                       cudaMemcpyHostToDevice));
//        }
//        if (gpu_bindingSites == NULL) {
    auto bindingSites = SysParams::Chemistry().bindingSites[0];//filType dependant
    int cpu_bindingSites[bindingSites.size()];
    int iii = 0;
    for (auto bs:bindingSites)
    {cpu_bindingSites[iii] = int(bs); iii++;}
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_bindingSites, bindingSites.size() *
                                       sizeof(int)), "cuda data transfer", "SubSystem.cu");
    CUDAcommon::handleerror(cudaMemcpy(gpu_bindingSites, cpu_bindingSites,
                                       bindingSites.size() *  sizeof(int), cudaMemcpyHostToDevice));
    //        }
    //@}
}

void SubSystem::getallpossiblelinkerbindingsitesCUDA(LinkerBindingManager* lManager,
                                                     int* cmon_state_linker){
    lManager->assigncudavars();
    cudaStream_t  s;
    if(numbindmgrs + 1 > strvec.size() )
    { cudaStreamCreate(&s); strvec.push_back(s);}
    else
        s = strvec.at(numbindmgrs);
    numbindmgrs++;
//    int *cmon_state_linker = cylcylnlvars.gpu_cmon_state_linker;
    //1. Assign optimal blocks and threads
    vector<int> blocksnthreads;
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
    int nint = lManager->getNLsize();
    int *NL = lManager->getNLCUDA();
    int *numNLpairs = lManager->getNLsizeCUDA();
    int *numpairs = lManager->getpossiblebindingssizeCUDA();
    floatingpoint *params2 = lManager->getdistancesCUDA();
    std::cout<<"Total Linker NL size "<<nint<<endl;
//            int *numpairs, test[1];test[0] = 0;
//            CUDAcommon::handleerror(cudaMalloc((void **) &numpairs, sizeof(int)), "cuda data transfer", "SubSystem.cu");
//            CUDAcommon::handleerror(cudaMemcpy(numpairs, test, sizeof(int), cudaMemcpyHostToDevice));
    //2. Calculate binding sites
    if (nint > 0) {
        int *gpu_possibleBindings;

        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_possibleBindings, SysParams::Chemistry()
                                           .numBindingSites[0] * 5 * nint * sizeof(int)));
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                           updateAllPossibleBindingsCUDA, 0, 0);
        blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreads.push_back(blockSize);
        std::cout << "Linker blocks and threads " << blocksnthreads.at(0) << " " << blocksnthreads.at(1)
        << endl;
        resetintvariableCUDA<<<1,1,0,s>>>(numpairs);
        //call CUDA kernel function
        updateAllPossibleBindingsCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                                  (coord, beadSet, cylID, filID, filType, cmpID, NL, numNLpairs, numpairs,
                                                                                          gpu_params, params2, gpu_possibleBindings, cmon_state_linker,
                                                                                          gpu_bindingSites);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
        //copy binding sites back to CPU
        int *cpu_numpairs, *possibleBindings;
        CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                "Subsystem.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                s),"Copy", "Subsystem.cu");
        CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
        CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
        std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
        numpairs_vec.push_back(cpu_numpairs);
        if(cpu_numpairs[0] > 0) {
            CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
                                                  cudaHostAllocDefault), "Copy", "Subsystem.cu");
            CUDAcommon::handleerror(
                    cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] *
                                                                            sizeof(int), cudaMemcpyDeviceToHost,
                                    s), "Copy", "Subsystem.cu");
            possibleBindings_vec.push_back(possibleBindings);
        }
        gpu_possibleBindings_vec.push_back(gpu_possibleBindings);
//                int cpu_numpairs[1];
//                cudaMemcpy(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost);
//                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
//                int possibleBindings[5 * cpu_numpairs[0]];
//                cudaMemcpy(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
//                           cudaMemcpyDeviceToHost);
    }
    lManager->freecudavars();
    //Free NL numpairs
    CUDAcommon::handleerror(cudaFree(numNLpairs),"cudaFree", "SubSystem.cu");
}

void SubSystem::getallpossiblemotorbindingsitesCUDA(MotorBindingManager* mManager, int*
                                                    cmon_state_motor){
    mManager->assigncudavars();
    cudaStream_t  s;
    if(numbindmgrs + 1 > strvec.size() )
    { cudaStreamCreate(&s); strvec.push_back(s);}
    else
        s = strvec.at(numbindmgrs);
    numbindmgrs++;
//    int *cmon_state_motor = cylcylnlvars.gpu_cmon_state_motor;
    //2. Assign optimal blocks and threads
    vector<int> blocksnthreads;
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
    int nint = mManager->getNLsize();
    int *NL = mManager->getNLCUDA();
    int *numNLpairs = mManager->getNLsizeCUDA();
    int *numpairs = mManager->getpossiblebindingssizeCUDA();
    floatingpoint *params2 = mManager->getdistancesCUDA();
    std::cout<<"Total Motor NL size "<<nint<<endl;
    if (nint > 0) {
        int *gpu_possibleBindings;

        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_possibleBindings, SysParams::Chemistry()
                                           .numBindingSites[0] * 5 * nint * sizeof(int)));
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                           updateAllPossibleBindingsCUDA, 0, 0);
        blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreads.push_back(blockSize);
        std::cout << "Motor blocks and threads " << blocksnthreads.at(0) << " " << blocksnthreads.at(1)
        << endl;
        resetintvariableCUDA << < 1, 1, 0, s >> > (numpairs);

        updateAllPossibleBindingsCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                                  (coord, beadSet, cylID, filID, filType, cmpID, NL, numNLpairs, numpairs,
                                                                                          gpu_params, params2,
                                                                                          gpu_possibleBindings, cmon_state_motor,
                                                                                          gpu_bindingSites);

        //copy back to CPU
        int *cpu_numpairs, *possibleBindings;
        CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                "Subsystem.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                s),"Copy", "Subsystem.cu");
        CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
        CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
        std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
        numpairs_vec.push_back(cpu_numpairs);
        if(cpu_numpairs[0] > 0) {
            CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
                                                  cudaHostAllocDefault), "Copy", "Subsystem.cu");
            CUDAcommon::handleerror(
                    cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] *
                                                                            sizeof(int), cudaMemcpyDeviceToHost,
                                    s), "Copy", "Subsystem.cu");
            possibleBindings_vec.push_back(possibleBindings);
        }
        gpu_possibleBindings_vec.push_back(gpu_possibleBindings);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
//                CUDAcommon::handleerror(cudaDeviceSynchronize());
//                //copy back to CPU
//                int cpu_numpairs[1];
//                cudaMemcpy(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost);
//                std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
//                int possibleBindings[5 * cpu_numpairs[0]];
//                cudaMemcpy(possibleBindings, gpu_possibleBindings, 5 * cpu_numpairs[0] * sizeof(int),
//                           cudaMemcpyDeviceToHost);

    }
    mManager->freecudavars();
    //Free NL numpairs
    CUDAcommon::handleerror(cudaFree(numNLpairs),"cudaFree", "SubSystem.cu");
}

void SubSystem::getallpossiblebrancherbindingsitesCUDA(BranchingManager* bManager,
                                                       int* cmon_state_brancher) {
    bManager->assigncudavars();
    cudaStream_t  s;
    if(numbindmgrs + 1 > strvec.size() )
    { cudaStreamCreate(&s); strvec.push_back(s);}
    else
        s = strvec.at(numbindmgrs);
    numbindmgrs++;
//    int *cmon_state_brancher = cylcylnlvars.gpu_cmon_state_brancher;
    //2. Assign optimal blocks and threads
    vector<int> blocksnthreads;
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the maximum occupancy for a full device launch
    int nint = Cylinder::getCylinders().size();
    //int *NL = bManager->getNLCUDA();
    //int *numNLpairs = bManager->getNLsizeCUDA();
    int *numpairs = bManager->getpossiblebindingssizeCUDA();
    floatingpoint *params2 = bManager->getdistancesCUDA();
    int *zone = bManager->getzoneCUDA();
    //std::cout<<"Total Motor NL size "<<nint<<endl;
    if (nint > 0) {
        //Boundary plane
        auto beList = BoundaryElement::getBoundaryElements();
        int nbe = BoundaryElement::getBoundaryElements().size();
        floatingpoint *beListplane = new floatingpoint[4 * nbe];
        floatingpoint *gpu_beListplane;
        for (int i = 0; i < nbe; i++) {

            if(dynamic_cast<PlaneBoundaryElement*>(beList[i])) {
                floatingpoint *x = new floatingpoint[4];
                beList[i]->elementeqn(x);
                beListplane[4 * i] = x[0];
                beListplane[4 * i +1] = x[1];
                beListplane[4 * i +2] = x[2];
                beListplane[4 * i +3] = x[3];
            }
            else{
                cout<<"CUDA cannot handle non-plane type boundaries. Exiting..."<<endl;
                exit(EXIT_FAILURE);
            }
        }
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beListplane, 4 * nbe * sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemcpy(gpu_beListplane, beListplane, 4 * nbe * sizeof(floatingpoint), cudaMemcpyHostToDevice));
        //
        int *gpu_possibleBindings;

        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_possibleBindings, SysParams::Chemistry()
                                           .numBindingSites[0] * 3 * nint * sizeof(int)));
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
                                           updateAllPossibleBindingsCUDA, 0, 0);
        blocksnthreads.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreads.push_back(blockSize);
        std::cout << "Brancher blocks and threads " << blocksnthreads.at(0) << " " << blocksnthreads.at(1)
                  << endl;
        resetintvariableCUDA << < 1, 1, 0, s >> > (numpairs);

        updateAllPossibleBindingsBrancherCUDA << < blocksnthreads[0], blocksnthreads[1],0,s >> >
                                                                                          (coord, beadSet, cylID, filID, filType, cmpID, numpairs,
                                                                                                  gpu_params, params2, zone, gpu_possibleBindings,gpu_bindingSites,
                                                                                                  cmon_state_brancher, gpu_beListplane);

//                CUDAcommon::handleerror(cudaDeviceSynchronize());
        //copy back to CPU
        int *cpu_numpairs, *possibleBindings;
        CUDAcommon::handleerror(cudaHostAlloc(&cpu_numpairs, sizeof(int), cudaHostAllocDefault),"Copy",
                                "Subsystem.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(cpu_numpairs, numpairs, sizeof(int), cudaMemcpyDeviceToHost,
                                                s),"Copy", "Subsystem.cu");
        CUDAcommon::handleerror(cudaStreamSynchronize(s),"Stream Sync","Subsystem.cu");
        //CUDAcommon::handleerror(cudaFree(NL),"cudaFree","NeighborListImpl.cu");
        std::cout << "Number of possibleBindings " << cpu_numpairs[0] << endl;
        numpairs_vec.push_back(cpu_numpairs);
        if(cpu_numpairs[0] > 0) {
            CUDAcommon::handleerror(cudaHostAlloc(&possibleBindings, cpu_numpairs[0] * 3 * sizeof(int),
                                                  cudaHostAllocDefault), "Copy", "Subsystem.cu");
            CUDAcommon::handleerror(
                    cudaMemcpyAsync(possibleBindings, gpu_possibleBindings, cpu_numpairs[0] * 3 *
                                                                            sizeof(int), cudaMemcpyDeviceToHost,
                                    s), "Copy", "Subsystem.cu");
            possibleBindings_vec.push_back(possibleBindings);
        }
        gpu_possibleBindings_vec.push_back(gpu_possibleBindings);
        CUDAcommon::handleerror(cudaFree(gpu_beListplane),"cudaFree", "Subsystem.cu");
    }
    bManager->freecudavars();
}

void SubSystem::terminatebindingsitesearchCUDA(){
    CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "Subsystem.cu");
    CUDAcommon::handleerror(cudaFree(gpu_bindingSites),"cudaFree", "Subsystem.cu");
    //Synchronize streams
    for(auto s:strvec) CUDAcommon::handleerror(cudaStreamSynchronize(s),"stream sync","SubsSystem.cu");
    //Delete sterams
    if(CUDAcommon::getCUDAvars().conservestreams == false) {
        for (auto s:strvec)
            CUDAcommon::handleerror(cudaStreamDestroy(s), "stream destroy", "SubsSystem.cu");
        strvec.clear();
    }
    //clear all possible bindings.
    for(auto& c:_compartmentGrid->getCompartments()){
        for(auto &Mgr : c->getFilamentBindingManagers()) {
            Mgr->clearpossibleBindings();
        }
    }
}

void SubSystem::assigntorespectivebindingmanagersCUDA(){
    Compartment* C0 = _compartmentGrid->getCompartments()[0].get();
    int count = 0;
    for(auto &manager : C0->getFilamentBindingManagers()) {
        LinkerBindingManager *lManager;
        MotorBindingManager *mManager;
        BranchingManager *bManager;
        //Linkers
        if ((lManager = dynamic_cast<LinkerBindingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[5* i];
                int cIndex = possibleBindings[5 * i +1];
                short cbs = short(possibleBindings[5 * i + 2]);
                int cnIndex = possibleBindings[5 * i +3];
                short cnbs = short(possibleBindings[5 * i + 4]);
                auto cylinder = Cylinder::getCylinders()[cIndex];
                auto ncylinder = Cylinder::getCylinders()[cnIndex];
                //get the compartment.
                Compartment* cmp = &GController::getCompartment(cID);
                //get corresponding binding manager
                for(auto &cmanager : cmp->getFilamentBindingManagers()) {
                    if ((lManager = dynamic_cast<LinkerBindingManager *>(cmanager.get()))) {
                        auto t1 = tuple<CCylinder*, short>(cylinder->getCCylinder(), cbs);
                        auto t2 = tuple<CCylinder*, short>(ncylinder->getCCylinder(), cnbs);
                        cmanager->appendpossibleBindings(t1,t2);
                    }
                }
            }
            if(numpairs[0] > 0)
                count++;
        }
            //MOTORS
        else if ((mManager = dynamic_cast<MotorBindingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[5 * i];
                int cIndex = possibleBindings[5 * i +1];
                short cbs = short(possibleBindings[5 * i + 2]);
                int cnIndex = possibleBindings[5 * i +3];
                short cnbs = short(possibleBindings[5 * i + 4]);
                auto cylinder = Cylinder::getCylinders()[cIndex];
                auto ncylinder = Cylinder::getCylinders()[cnIndex];
                //get the compartment
                Compartment* cmp = &GController::getCompartment(cID);
                //get corresponding binding manager
                for(auto &cmanager : cmp->getFilamentBindingManagers()) {
                    if ((mManager = dynamic_cast<MotorBindingManager *>(cmanager.get()))) {
                        auto t1 = tuple<CCylinder*, short>(cylinder->getCCylinder(), cbs);
                        auto t2 = tuple<CCylinder*, short>(ncylinder->getCCylinder(), cnbs);
                        cmanager->appendpossibleBindings(t1,t2);
                    }
                }
            }
            if(numpairs[0] > 0)
                count++;
        }
        else if ((bManager = dynamic_cast<BranchingManager *>(manager.get()))) {
            auto numpairs = numpairs_vec[count];
            auto possibleBindings = possibleBindings_vec[count];
            for(auto i = 0; i < numpairs[0]; i++){
                int cID = possibleBindings[3 * i];
                int cIndex = possibleBindings[3 * i +1];
                short cbs = short(possibleBindings[3 * i + 2]);
                auto cylinder = Cylinder::getCylinders()[cIndex];
                Compartment* cmp = &GController::getCompartment(cID);
                for(auto &cmanager : cmp->getFilamentBindingManagers()) {
                    if ((bManager = dynamic_cast<BranchingManager *>(cmanager.get()))) {
                        auto t1 = tuple<CCylinder*, short>(cylinder->getCCylinder(), cbs);
                        dynamic_cast<BranchingManager *>(cmanager.get())->appendpossibleBindings(t1);
                    }
                }
            }
            if(numpairs[0] > 0)
                count++;
            std::cout<<endl;
        }
    }
}

#endif
bool SubSystem::initialize = false;
floatingpoint SubSystem::SIMDtime  = 0.0;
floatingpoint SubSystem::SIMDtimeV2  = 0.0;
floatingpoint SubSystem::HYBDtime  = 0.0;
floatingpoint SubSystem::timeneighbor  = 0.0;
floatingpoint SubSystem::timedneighbor  = 0.0;
floatingpoint SubSystem::timetrackable  = 0.0;

} // namespace medyan
