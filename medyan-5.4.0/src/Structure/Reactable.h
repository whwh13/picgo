
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

#ifndef MEDYAN_Reactable_h
#define MEDYAN_Reactable_h

#include "common.h"
#include "Util/DoubleLinkedList.h"

namespace medyan {
/// An abstract base class for a reactable element in the SubSystem.

/*! The main function of the Reactable class is to implement updateReactionRates(),
 *  so that the reactions related to any object extending this class can be updated
 *  by the SubSystem.
 */
class Reactable {
    
protected:
    Reactable() {
    	reactableList.addNode(this, dllnode);
    }

    
public:
    ///Update the reactions in this element
    /// @note - this update could be due to an updated force minimization,
    /// a set of chemical steps, or any other event in the SubSystem. This
    /// function will be called by the SubSystem on all Reactables.
    virtual void updateReactionRates() = 0;
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Reactable() noexcept {
    	reactableList.removeNode(dllnode);
    }

	static const linkedlist::DoubleLinkedList<Reactable>& getReactableList(){
    	return reactableList;
    }

private:
	//static dll list here
	inline static linkedlist::DoubleLinkedList<Reactable> reactableList;
	//Node of this instance
	linkedlist::DLLNode<Reactable> dllnode;
};

} // namespace medyan

#endif
