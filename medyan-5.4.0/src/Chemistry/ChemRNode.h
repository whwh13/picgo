
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

#ifndef MEDYAN_ChemRNode_h
#define MEDYAN_ChemRNode_h

namespace medyan {
/// This is an abstract base class for classes that need to be associated with the
/// given Reaction object.
class RNode{
public:
    /// Dtor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~RNode() noexcept {}
    
    /// This method is called by Reaction::activateReaction(). Its effect depends on the
    /// underlying stochastic simulation algorithm. For example, in the NRM algorithm, a
    /// new tau and a are computed and the heap is updated.
    virtual void activateReaction() = 0;

    /// This method is called by Reaction::passivateReaction(). Its effect depends on
    /// the underlying stochastic simulatin algorithm. For example, in the NRM
    /// algorithm, a tau is set to infinity and the heap is updated.
    virtual void passivateReaction() = 0;
    
    /// Return true if the Reaction is currently passivated
    virtual bool isPassivated() const = 0;

    ///Temporarily adding printSelf to this.
    virtual void printSelf() const = 0;
};

} // namespace medyan

#endif
