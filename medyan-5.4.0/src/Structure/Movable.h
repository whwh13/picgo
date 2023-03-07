
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

#ifndef MEDYAN_Movable_h
#define MEDYAN_Movable_h

#include "common.h"
#include "Util/DoubleLinkedList.h"

namespace medyan {
/// An abstract base class for a movable element in the SubSystem.

/*! The main function of the Movable class is to implement updatePosition(),
 *  so that the position of any object extending this class can be updated
 *  by the SubSystem.
 */
class Movable {
    
protected:
    Movable() {
        movableList.addNode(this, dllnode);
    }
    
public:
    /// Update the position of this element
    /// @note - this update could be due to an updated force minimization,
    /// a set of chemical steps, or any other event in the SubSystem. This
    /// function will be called by the SubSystem on all Movables.
    virtual void updatePosition() = 0;
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Movable() noexcept {
        movableList.removeNode(dllnode);
    }

    static const linkedlist::DoubleLinkedList<Movable>& getMovableList(){
        return movableList;
    }
private:
    //static dll list here
    inline static linkedlist::DoubleLinkedList<Movable> movableList;
    //Node of this instance
    linkedlist::DLLNode<Movable> dllnode;
};

} // namespace medyan

#endif
