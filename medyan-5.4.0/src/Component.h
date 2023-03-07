
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

#ifndef MEDYAN_Component_h
#define MEDYAN_Component_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Composite;
class Visitor;
class SpeciesVisitor;
class ReactionVisitor;

/// The base class for the Composite pattern hieararchy

/*! 
 *  The Composite pattern allows building of complex hieararchical objects, with 
 *  convenient methods for applying a function to all nodes (i.e. the Visitor pattern).
 *  Each node in the hieararchy may have a parent and may contain several children 
 *  nodes. A class that is derived directly from Component and not from Composite,
 *  cannot contain children, i.e. it is a leaf node.
 *  @note Component objects may contain Species and ReactionBase collections, however, 
 *  this is treated seperately from the Composite pattern (i.e. separate methods exist 
 *  for the corresponding access to elements, changes, etc.)
 */
class Component {
private:
    Composite *_parent; ///< The parent of this object. May be a nullptr if this object
                        ///< has no parent.

public:
    /// Default Constructor; Parent is assigned to nullptr
    Component() : _parent(nullptr) {}
    
    /// Virtual Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Component() noexcept {};
    
    /// When this function is applied to a ConditionalVisitor v, the corresponding
    /// v.visit(this) is called, and v is further applied to all children of this node
    /// recursively. However, the ConditionalVisitor allows the visit() function to be
    /// applied only selectively to nodes that conform to specific criteria.
    virtual bool apply (Visitor &v);

    /// Implements the apply_if() method of the Component class by recursively applying
    /// it to itself and all its children that contain Species.
    virtual bool apply (SpeciesVisitor &v) {return apply_impl(v);}

    /// Applies SpeciesVisitor v to every Species* object directly owned by this node.
    /// This method needs to be overriden by descendent classes that contain Species.
    virtual bool apply_impl(SpeciesVisitor &v) {return true;}

    /// Implements the apply_if() method of the Component class by recursively applying
    /// it to itself and all its children that contain ReactionBase.
    virtual bool apply (ReactionVisitor &v) {return apply_impl(v);}
    
    /// Applies ReactionBaseVisitor v to every ReactionBase* object directly owned by
    /// this node. This method needs to be overriden by descendent classes that contain
    /// ReactionBase.
    virtual bool apply_impl(ReactionVisitor &v) {return true;}

    /// Returns the pointer to the parent node. The returned value could be a nullptr if
    /// a parent does not exist.
    Composite* getParent() {return _parent;}
    const Composite* getParent()const { return _parent; }
    
    /// Sets the parent of this node to other.
    void setParent (Composite *other) {_parent=other;}
    
    /// Returns true if this node is the root node. A root node has no parent, and the
    /// corresponding getParent() call would return a nullptr.
    bool isRoot() const {return _parent==nullptr? true : false;}
    
    /// Returns the root node of the hieararchy to which this node belongs to.
    Composite* getRoot();
    
    /// Returns the number of children of this node.
    /// @note Species and ReactionBase objects are not counted.
    virtual size_t numberOfChildren() const {return 0;}
    
    /// Returns true if this node is a Composite node
    virtual bool isComposite() const {return false;}

    /// Return a string indicating the full name of this node (presumably used mainly
    /// for debugging)
    virtual string getFullName() const {return "Component";};
    
    /// Return the total number of nodes contained under this node's hieararchy
    /// @note This is a recursive call, and all nodes under this node are visited.
    virtual size_t countDescendents() const {return 0;}

    /// Return the number of Species contained under this node's hieararchy
    /// @note This is a recursive call, and all nodes under this node are visited.
    virtual size_t countSpecies() const {return 0;}
            
    /// Return the number of ReactionBase objects contained under this node's hieararchy
    /// @note This is a recursive call, and all nodes under this node are visited.
    virtual size_t countReactions() const {return 0;};
    
    /// Prints information about this node. Useful for debugging.
    virtual void printSelf() const = 0;
    
    /// Return the type of this element as an integer index
    virtual int getType() = 0;
};

} // namespace medyan

#endif
