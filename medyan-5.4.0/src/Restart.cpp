
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

#include "Restart.h"

namespace medyan {

void Restart::readNetworkSetup() {
    //Flush the filestream pointer
    _inputFile.clear();
    //Go to the first entry in the file.
    _inputFile.seekg(0);

    string line;
    //get first line
    getline(_inputFile, line);
    vector<string> lineVector = split<string>(line);
    _rsystemdata.time = atof(((lineVector[1]).c_str()));
    //get line
    while(getline(_inputFile, line)) {
        //Continue if commented
        if(line.find("#") != string::npos) { continue; }
        vector<string> lineVector = split<string>(line);
        if(line.size()>0) {

            if(lineVector[0] == "NFIL"){
                getline(_inputFile, line);
                //Split string based on space delimiter
                vector<string> lineVector = split<string>(line);
                _rsystemdata.Nfil = atoi(((lineVector[0]).c_str()));
                _rsystemdata.Ncyl = atoi(((lineVector[1]).c_str()));
                _rsystemdata.Nbead = atoi(((lineVector[2]).c_str()));
                _rsystemdata.Nlink = atoi(((lineVector[3]).c_str()));
                _rsystemdata.Nmotor = atoi(((lineVector[4]).c_str()));
                _rsystemdata.Nbranch = atoi(((lineVector[5]).c_str()));
                //resize data structures
                _rBData.bsidvec.resize(_rsystemdata.Nbead);
                _rBData.filidvec.resize(_rsystemdata.Nbead);
                _rBData.filpos.resize(_rsystemdata.Nbead);
                _rBData.coordvec.resize(3*_rsystemdata.Nbead);
                _rBData.forceAuxvec.resize(3*_rsystemdata.Nbead);
            }
            if (lineVector[0] == "BEAD") {
            //Get lines till the line is not empty
                getline(_inputFile, line);
                while (line.size() > 0) {
                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);
                    //Bead stable ID
                    auto beadstableid = atoi(((lineVector[0]).c_str()));
                    _rBData.bsidvec[beadstableid] = beadstableid;
                    //Filament ID
                    _rBData.filidvec[beadstableid] = (atoi(lineVector[1].c_str()));
                    //Filament pos
                    _rBData.filpos[beadstableid] = (atoi(lineVector[2].c_str()));
                    //Copy coords & forces
                    short dim = 0;
                    for (auto it = lineVector.begin() + 3;
                         it != lineVector.begin() + 6; it++) {
                        _rBData.coordvec[3 * beadstableid + dim] = (atof((*it).c_str()));
                        dim++;
                    }
                    dim = 0;
                    for (auto it = lineVector.begin() + 6;
                         it != lineVector.begin() + 9; it++){
                        _rBData.forceAuxvec[3 * beadstableid + dim] = (atof((*it).c_str()));
                        dim++;
                    }
                    getline(_inputFile, line);
                }
            }
            else if (lineVector[0] == "CYLINDER") {
            //Get lines till the line is not empty

                getline(_inputFile, line);
                while (line.size() > 0) {
                    restartCylData _rcyldata;

                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);
                    //Cylinder stable index
                    _rcyldata.cylsid = atoi(((lineVector[0]).c_str()));
                    //Filament ID
                    _rcyldata.filid = atoi((lineVector[1]).c_str());
                    _rcyldata.filtype = atoi((lineVector[2]).c_str());
                    _rcyldata.filpos = atoi((lineVector[3]).c_str());
                    //Filament pos
                    //Bead stable indices
                    _rcyldata.beadsidpairvec[0] = atoi((lineVector[4]).c_str());
                    _rcyldata.beadsidpairvec[1] = atoi((lineVector[5]).c_str());
                    //minus end or plus end?
                    _rcyldata.endstatusvec[0] = atoi((lineVector[6]).c_str());
                    _rcyldata.endstatusvec[1] = atoi((lineVector[7]).c_str());
                    //minus/plus end type
                    _rcyldata.endtypevec[0] = atoi((lineVector[8]).c_str());
                    _rcyldata.endtypevec[1] = atoi((lineVector[9]).c_str());
                    //endmonomerpos
                    _rcyldata.endmonomerpos[0] = atoi((lineVector[10]).c_str());
                    _rcyldata.endmonomerpos[1] = atoi((lineVector[11]).c_str());
                    //totalmonomers
                    _rcyldata.totalmonomers = atoi((lineVector[12]).c_str());
                    //eqlen
                    _rcyldata.eqlen = atof((lineVector[13]).c_str());
                    //append to the vector
                    _rCDatavec.push_back(_rcyldata);

                    //get next cylinder data
                    getline(_inputFile, line);
                }
            }
            else if (lineVector[0] == "FILAMENT") {
            //Get lines till the line is not empty

                getline(_inputFile, line);
                while (line.size() > 0) {
                    restartFilData _rfildata;

                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);
                    //Filament ID
                    _rfildata.filid = atoi((lineVector[0]).c_str());
                    //Filament Type
                    _rfildata.filType = atoi((lineVector[1]).c_str());
                    //Cylinder id vec
                    for (auto it = lineVector.begin() + 2; it != lineVector.end(); it++) {
                        _rfildata.cylsidvec.push_back(atoi((*it).c_str()));
                    }
                    //append to the vector
                    _rFDatavec.push_back(_rfildata);

                    //get next filament data
                    getline(_inputFile, line);
                }
            }

            else if (lineVector[0] == "LINKER") {
                //Get lines till the line is not empty

                getline(_inputFile, line);
                while (line.size() > 0) {
                    restartLinkerData _rldata;
                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);
                    //Motor ID
                    _rldata.linkerid = atoi((lineVector[0]).c_str());
                    //Linker Type
                    _rldata.linkerType = short(atoi((lineVector[1]).c_str()));
                    //Cyl ID1
                    _rldata.cylid1 = atoi((lineVector[2]).c_str());
                    //Cyl ID2
                    _rldata.cylid2 = atoi((lineVector[3]).c_str());
                    //pos1
                    _rldata.pos1 = atoi((lineVector[4]).c_str());
                    //pos2
                    _rldata.pos2 = atoi((lineVector[5]).c_str());
                    //eqlen
                    _rldata.eqlen = atof((lineVector[6]).c_str());
                    //diffusing species name
                    _rldata.diffusingspeciesname = lineVector[7];

                    _rLDatavec.push_back(_rldata);
                    //get next filament data
                    getline(_inputFile, line);
                }
            }

            else if (lineVector[0] == "MOTOR") {
                //Get lines till the line is not empty

                getline(_inputFile, line);
                while (line.size() > 0) {
                    restartMotorData _rmdata;
                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);
                    //Motor ID
                    _rmdata.motorid = atoi((lineVector[0]).c_str());
                    //Motor Type
                    _rmdata.motorType = short(atoi((lineVector[1]).c_str()));
                    //Cyl ID1
                    _rmdata.cylid1 = atoi((lineVector[2]).c_str());
                    //Cyl ID2
                    _rmdata.cylid2 = atoi((lineVector[3]).c_str());
                    //pos1
                    _rmdata.pos1 = atoi((lineVector[4]).c_str());
                    //pos2
                    _rmdata.pos2 = atoi((lineVector[5]).c_str());
                    //eqlen
                    _rmdata.eqlen = atof((lineVector[6]).c_str());
                    //diffusing species name
                    _rmdata.diffusingspeciesname = lineVector[7];
                    if(lineVector.size()>8) {
                        //numHeads
                        _rmdata.numHeads = atoi((lineVector[8]).c_str());
                        //numBoundHeads
                        _rmdata.numBoundHeads = atof((lineVector[9]).c_str());
                    }

                    _rMDatavec.push_back(_rmdata);
                    //get next filament data
                    getline(_inputFile, line);
                }
            }

            else if (lineVector[0] == "BRANCHING") {
                //Get lines till the line is not empty
                getline(_inputFile, line);
                while (line.size() > 0) {
                    restartBrancherData _rbdata;
                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);
                    //Brancher ID
                    _rbdata.branchid = atoi((lineVector[0]).c_str());
                    //Brancher Type
                    _rbdata.branchType = short(atoi((lineVector[1]).c_str()));
                    //Cyl ID1
                    _rbdata.cylid1 = atoi((lineVector[2]).c_str());
                    //Cyl ID2
                    _rbdata.cylid2 = atoi((lineVector[3]).c_str());
                    //pos1
                    _rbdata.pos1 = atoi((lineVector[4]).c_str());
                    //eqlen
                    _rbdata.eqlen = atof((lineVector[5]).c_str());
                    //diffusing species name
                    _rbdata.diffusingspeciesnamebranch = lineVector[6];
                    if(lineVector.size() == 8 )
                        _rbdata.diffusingspeciesnameactin = lineVector[7];
                    else
                        _rbdata.diffusingspeciesnameactin = "A";

                    _rBDatavec.push_back(_rbdata);
                    //get next filament data
                    getline(_inputFile, line);
                }
            }
            /* parse diffusing species copy number in each compartment*/
            if (lineVector[0] == "COMPARTMENT") {
                //Get lines till the line is not empty
                getline(_inputFile, line);
                while (line.size() > 0) {
                    restartCompartmentDiffusingData _rcddata;
                    //Split string based on space delimiter
                    vector<string> lineVector = split<string>(line);

                    _rcddata.id = atoi((lineVector[0]).c_str());
                    //Copy species name
                    for (auto n = 1; n < lineVector.size(); n = n + 2)
                        _rcddata.speciesnamevec.push_back(lineVector[n].c_str());
                    //Copy species copy number
                    for (auto n = 2; n < lineVector.size(); n = n + 2)
                        _rcddata.copynumvec.push_back(atoi(lineVector[n].c_str()));

                    _rCmpDDatavec.push_back(_rcddata);
                    //get next filament data
                    getline(_inputFile, line);
                }
            }

            /* parse bulk species copy number*/
            if (lineVector[0] == "BULKSPECIES") {
                //Get lines till the line is not empty
                getline(_inputFile, line);
                while (line.size() > 0) {
                    vector<string> lineVector = split<string>(line);
                    for (auto n = 0; n < lineVector.size(); n = n + 2)
                        _rbdata.speciesnamevec.push_back(lineVector[n].c_str());
                    //Copy species copy number
                    for (auto n = 1; n < lineVector.size(); n = n + 2)
                        _rbdata.copynumvec.push_back(atoi(lineVector[n].c_str()));
                    //get next bulkspecies data
                    getline(_inputFile, line);
                }
            }

            /* parse total copy number of each species*/
            if (lineVector[0] == "TALLY") {
                //Get lines till the line is not empty
                getline(_inputFile, line);
                string delimiter = ":";
                while (line.size() > 0) {
                    vector<string> lineVector = split<string>(line);
                    auto pos = lineVector[0].find(delimiter);
                    _rbdata.speciesnamevec.push_back(lineVector[0].substr(0,pos));
                    _rbdata.copynumvec.push_back(atoi(lineVector[1].c_str()));
                    //get next bulkspecies data
                    getline(_inputFile, line);
                }
            }
        }
    }
    cout<<"numMotors="<<_rMDatavec.size()<<endl;
    cout<<"numLinkers="<<_rLDatavec.size()<<endl;
    cout<<"numBranchers="<<_rBDatavec.size()<<endl;
}

void Restart::setupInitialNetwork() {

    //Step 1. Create dummy filaments
    //Step 2. Create Beads
    //Step 3. Create cylinders & set plus/minus ends where necessary
    //Step 4. Associate cylinders with each filament by adding it to the cylinder vector
    // in the appropriate order starting from minus end all the way to plus end cylinder.

    map<int, Filament*> filamentmap;//makes it easy to access Filament pointer from
    // Filament ID. Filaments do not follow stable index protocol and hence need not have
    // continuous ID values.

    for(auto &fil : _rFDatavec) {
        //Create dummy filament
        fil.filamentpointer = _subSystem->addTrackable<Filament>(_subSystem, fil.filType);
        //override ID
        fil.filamentpointer->overrideId(fil.filid);
        //add to map
        filamentmap[fil.filid] = fil.filamentpointer;
    }
    cout<<endl;
    cout<<"Num filaments created "<<Filament::getFilaments().size()<<endl;

    for(unsigned int b=0;b<_rBData.bsidvec.size();b++){
        auto bID = _rBData.bsidvec[b];
        auto filptr = filamentmap[_rBData.filidvec[b]];
        //Extract part of the vector.
        vector<floatingpoint> tempcoord(_rBData.coordvec.begin()+3*bID, _rBData.coordvec
        .begin()+3*bID+3);
        //Initialize beads
        auto pBead = _subSystem->addTrackable<Bead>(tempcoord, filptr, _rBData.filpos[b]);
        //Copy Forces
        for(unsigned int dim = 0; dim < 3; dim++)
            pBead->force[dim] = _rBData.forceAuxvec.data()[3*bID+dim];
    }
    cout<<"Num beads created "<<Bead::getBeads().size()<<endl;

    for(auto &cyl : _rCDatavec){
        auto b1 = Bead::getBeads()[cyl.beadsidpairvec[0]];
        auto b2 = Bead::getBeads()[cyl.beadsidpairvec[1]];
        auto filptr = filamentmap[cyl.filid];
        auto _filType = cyl.filtype;
        //initialize cylinder
        Cylinder* c0 = _subSystem->addTrackable<Cylinder> (filptr, b1, b2, _filType,
                                                           cyl.filpos, false, false,
                                                           true, cyl.eqlen);
        cyl.cylinderpointer = c0;

        //set minusend or plusend
        if(cyl.endstatusvec[0])
            c0->setMinusEnd(true);
        else if(cyl.endstatusvec[1])
            c0->setPlusEnd(true);
    }
    cout<<"Num cylinders created "<<Cylinder::getCylinders().size()<<endl;

    for(auto fil : _rFDatavec) {
        vector<Cylinder*> cylvector;
        //Go through cylinder stable indices that should be part of the filament and
        // append the Cylinder pointer to a vector.
        for(auto cylsid:fil.cylsidvec){
            cylvector.push_back(_rCDatavec[cylsid].cylinderpointer);
        }
        fil.filamentpointer->initializerestart(cylvector, _rCDatavec);
    }
}

void Restart::addtoHeapLinkerMotorBrancher(){
//STEP #2. ADD bound Linkers And Motors in inputfile into possible bindings.
    //set total number of reactions to fire.
    _numChemSteps = _rMDatavec.size() + _rLDatavec.size() + _rBDatavec.size();
    //Motors
    for(auto m:_rMDatavec) {
        int cidx1 = m.cylid1;
        int cidx2 = m.cylid2;
        int site1 = m.pos1;
        int site2 = m.pos2;
        Cylinder *c1 = _rCDatavec[cidx1].cylinderpointer;
        Cylinder *c2 = _rCDatavec[cidx2].cylinderpointer;
        if (c1->getId() > c2->getId()) {
            for (auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()) {
                if (dynamic_cast<MotorBindingManager *>(Mgr.get())) {
                    setdiffspeciesnumber(m.diffusingspeciesname,c1);
                    #ifdef NLORIGINAL
                    Mgr->appendpossibleBindings(c1->getCCylinder(),
                                                       c2->getCCylinder(), site1, site2);
                    #else
                    Mgr->appendpossibleBindingsstencil(m.motorType, c1->getCCylinder(),
                                                       c2->getCCylinder(), site1, site2);
                    #endif
                }
            }
        }
        else{
            for (auto &Mgr:c2->getCompartment()->getFilamentBindingManagers()) {
                if (dynamic_cast<MotorBindingManager *>(Mgr.get())) {
                    setdiffspeciesnumber(m.diffusingspeciesname,c2);
                    #ifdef NLORIGINAL
                    Mgr->appendpossibleBindings(c2->getCCylinder(),
                                                       c1->getCCylinder(), site2, site1);
                    #else
                    Mgr->appendpossibleBindingsstencil(m.motorType, c2->getCCylinder(),
                                                       c1->getCCylinder(), site2, site1);
                    #endif
                }
            }
        }
    }
    //Linkers
    for(auto l:_rLDatavec) {
        int cidx1 = l.cylid1;
        int cidx2 = l.cylid2;
        int site1 = l.pos1;
        int site2 = l.pos2;
        Cylinder *c1 = _rCDatavec[cidx1].cylinderpointer;
        Cylinder *c2 = _rCDatavec[cidx2].cylinderpointer;
        if (c1->getId() > c2->getId()) {
            for (auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()) {
                if (dynamic_cast<LinkerBindingManager *>(Mgr.get())) {
                    setdiffspeciesnumber(l.diffusingspeciesname,c1);
                    #ifdef NLORIGINAL
                    Mgr->appendpossibleBindings(c1->getCCylinder(),
                                                       c2->getCCylinder(), site1, site2);
                    #else
                    Mgr->appendpossibleBindingsstencil(l.linkerType, c1->getCCylinder(),
                                                       c2->getCCylinder(), site1, site2);
                    #endif
                }
            }
        }
        else{
            for (auto &Mgr:c2->getCompartment()->getFilamentBindingManagers()) {
                if (dynamic_cast<LinkerBindingManager *>(Mgr.get())) {
                    setdiffspeciesnumber(l.diffusingspeciesname,c2);
                    #ifdef NLORIGINAL
                    Mgr->appendpossibleBindings(c2->getCCylinder(),
                                                       c1->getCCylinder(), site2, site1);
                    #else
                    Mgr->appendpossibleBindingsstencil(l.linkerType, c2->getCCylinder(),
                                                       c1->getCCylinder(), site2, site1);
                    #endif
                }
            }
        }
    }
    //Brancher
    for(auto b:_rBDatavec) {
        int cidx1 = b.cylid1;
        int cidx2 = b.cylid2;
        int site1 = b.pos1;
        Cylinder *c1 = _rCDatavec[cidx1].cylinderpointer;
        Cylinder *c2 = _rCDatavec[cidx2].cylinderpointer;
/*		cout<<"c1 idx "<<c1->getStableIndex()<<endl;
        cout<<"c1 cmp "<<c1->getCompartment()->getId()<<", c2 cmp "<<c2->getCompartment()
        ->getId()<<endl;*/
        for (auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()) {
                if (dynamic_cast<BranchingManager *>(Mgr.get())) {
                    setdiffspeciesnumber(b.diffusingspeciesnamebranch,c1);
                    setdiffspeciesnumber(b.diffusingspeciesnameactin,c1);
                    #ifdef NLORIGINAL
                    Mgr->appendpossibleBindings(c1->getCCylinder(),
                                                       c2->getCCylinder(), site1, site2);
                    #else
                    Mgr->appendpossibleBindingsstencil(b.branchType, c1->getCCylinder(),
                                                       c2->getCCylinder(), site1, 0);
                    #endif
                }
            }
    }
}

void Restart::CBoundinitializerestart(){
    //linkers
    for(auto l:Linker::getLinkers()) {
        int c1 = l->getFirstCylinder()->getStableIndex();
        int c2 = l->getSecondCylinder()->getStableIndex();
        short ftype1 = l->getFirstCylinder()->getType();
        short ftype2 = l->getSecondCylinder()->getType();
        short pos1 = l->getFirstPosition()*SysParams::Geometry().cylinderNumMon[ftype1];
        short pos2 = l->getSecondPosition()*SysParams::Geometry().cylinderNumMon[ftype2];
        for(auto &rldata: _rLDatavec){
            if(!rldata.restartcompletion) {
                int rc1 = rldata.cylid1;
                int rc2 = rldata.cylid2;
                short rpos1 = rldata.pos1;
                short rpos2 = rldata.pos2;
                bool cndn1 = rc1 == c1 && rc2 == c2 && rpos1 == pos1 && rpos2 == pos2;
                bool cndn2 = rc1 == c2 && rc2 == c1 && rpos1 == pos2 && rpos2 == pos1;
                if(cndn1 || cndn2){
                    l->initializerestart(rldata.eqlen);
                    rldata.restartcompletion = true;
                }
            }
        }
    }
    //motors
    for(auto m:MotorGhost::getMotorGhosts()) {
        int c1 = m->getFirstCylinder()->getStableIndex();
        int c2 = m->getSecondCylinder()->getStableIndex();
        short ftype1 = m->getFirstCylinder()->getType();
        short ftype2 = m->getSecondCylinder()->getType();
        short pos1 = m->getFirstPosition()*SysParams::Geometry().cylinderNumMon[ftype1];
        short pos2 = m->getSecondPosition()*SysParams::Geometry().cylinderNumMon[ftype2];
        bool foundstatus = false;
        for(auto &rmdata: _rMDatavec){
            if(!rmdata.restartcompletion) {
                int rc1 = rmdata.cylid1;
                int rc2 = rmdata.cylid2;
                short rpos1 = rmdata.pos1;
                short rpos2 = rmdata.pos2;
                bool cndn1 = rc1 == c1 && rc2 == c2 && rpos1 == pos1 && rpos2 == pos2;
                bool cndn2 = rc1 == c2 && rc2 == c1 && rpos1 == pos2 && rpos2 == pos1;
                if(cndn1 || cndn2){
                    foundstatus = true;
                    m->initializerestart(rmdata.eqlen, rmdata.numHeads, rmdata.numBoundHeads);
                    rmdata.restartcompletion = true;
                }
            }
        }
        if(!foundstatus)
            LOG(ERROR)<<"Motor not found!"<<endl;
    }
    cout<<endl;
    //brancher
    for(auto b:BranchingPoint::getBranchingPoints()) {
        int c1 = b->getFirstCylinder()->getStableIndex();
        int c2 = b->getSecondCylinder()->getStableIndex();
        short ftype1 = b->getFirstCylinder()->getType();
        short pos1 = b->getPosition()*SysParams::Geometry().cylinderNumMon[ftype1];
        for(auto &rbdata: _rLDatavec){
            if(!rbdata.restartcompletion) {
                int rc1 = rbdata.cylid1;
                int rc2 = rbdata.cylid2;
                short rpos1 = rbdata.pos1;
                if(rc1 == c1 && rc2 == c2 && rpos1 == pos1){
                    b->initializerestart(rbdata.eqlen);
                    rbdata.restartcompletion = true;
                }
            }
        }
    }
}

bool Restart::crosscheck(){
    //Filament data
    for(auto &fil : _rFDatavec) {
        auto filptr = fil.filamentpointer;
        short counter = 0;
        auto cylsidvec = fil.cylsidvec;
        bool status = true;
        int startindex = 0;
        for(auto cylinfil:filptr->getCylinderVector()) {
            if(counter == 0){
                startindex = cylinfil->getPosition();
            }
            int currindex = cylinfil->getPosition();
            if(abs(startindex-currindex) != counter ){
                LOG(ERROR) << "Cylinder position does not match" << endl;
                LOG(ERROR) << startindex<<" "<<currindex<<endl;
                status = false;
                break;
            }
            if (cylinfil->getStableIndex() != cylsidvec[counter]) {
                LOG(ERROR) << "Cylinder sidx does not match" << endl;
                status = false;
                break;
            }
            counter++;
        }
        if(!status){
            cout<<"Datadump dictates Fil "<<fil.filid<<" be comprised of cylinder "
                                                       "sidx ";
            for(auto c:cylsidvec)
                cout<<c<<" ";
            cout<<endl;
            cout<<"But Fil "<<filptr->getId()<<" is comprised of cylinder sidx ";
            for(auto c:filptr->getCylinderVector())
                cout<<c->getStableIndex()<<" ";
            cout<<endl;
        }
    }
    //Cylinder data
    return true;
}

} // namespace medyan
