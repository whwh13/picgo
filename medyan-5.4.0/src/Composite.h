
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

#ifndef MEDYAN_Composite_h
#define MEDYAN_Composite_h

#include <string>
#include <vector>
#include <algorithm>
#include <typeinfo>

#include "common.h"

#include "Component.h"

namespace medyan {
//FORWARD DECLARATIONS
class Visitor;
class SpeciesVisitor;
class ReactionVisitor;

/// The aggregating class for the Composite pattern

/*! 
 *  The Composite pattern allows building of complex hieararchical objects, with 
 *  convenient methods for applying a function to all nodes (i.e. the Visitor pattern). 
 *  Each node in the hieararchy may have a parent and may contain several children 
 *  nodes. A class that is derived from Composite can contain children Component 
 *  objects. @note Composite objects may contain other collections,
 *  however, this is treated seperately from the Composite pattern (i.e. separate 
 *  methods exist for the corresponding access to elements, changes, etc.)
 */    
class Composite : public Component {
private:
    vector<unique_ptr<Component>> _children;///< Child node pointers of this
                                            ///< composite node

public:
    /// Default Constructor does nothing.
    Composite() :  Component() {}
    
    /// Virtual Destructor does nothing.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Composite() noexcept {}
    
    /// Implements the apply_if() method of the Component class by recursively applying
    /// it to itself and all its children.
    virtual bool apply (Visitor &v) override;
    
    /// Implements the apply_if() method of the Component class by recursively applying
    /// it to itself and all its children that contain Species.
    virtual bool apply (SpeciesVisitor &v) override;
    
    /// Implements the apply_if() method of the Component class by recursively applying
    /// it to itself and all its children that contain ReactionBase.
    virtual bool apply (ReactionVisitor &v) override;
    
    /// Returns true.
    virtual bool isComposite() const override {return true;}
    
    /// Returns the full name of this node.
    virtual string getFullName() const override {return "Composite";}; 
    
    /// Adds a Component child to this Composite node
    /// @param child - is a unique_ptr<Component>, hence, this node takes the memory
    /// ownership of the corresponding child pointer.
    /// @note the rvalue semantics - the unique_ptr<...> cannot be copied, but only
    /// can be moved
    virtual void addChild (unique_ptr<Component> &&child) {
        _children.push_back(move(child));
        _children.back()->setParent(this);
    }
    
    
    /// Remove *child from this node's children. The child's destructor will be called
    /// and the memory will be freed.
    virtual void removeChild (Component* child) {
        auto child_iter = find_if(_children.begin(),_children.end(),
                    [&child](const unique_ptr<Component> &element)
                     {return element.get()==child ? true : false;});
        if(child_iter!=_children.end())
            _children.erase(child_iter);
        else
            throw out_of_range("Composite::removeChild(): The name child not found");
    }
    
    // Transfer the specified component to another parent.
    //
    // Note:
    //   - If the component is not a child of this, an exception will be thrown.
    //   - The parent of the child will be updated as well.
    void transferChild(Component* pc, Composite& newParent) {
        // First, check whether pc is indeed a child of this.
        auto childIter = std::find_if(
            _children.begin(), _children.end(),
            [pc](const auto& element) { return element.get() == pc; }
        );

        if(childIter == _children.end()) {
            throw std::out_of_range("Composite::transferChild(): child is not found.");
        }
        else {
            // Move the child to new parent.
            newParent._children.push_back(std::move(*childIter));

            // Delete the original child record.
            _children.erase(childIter);

            // Update the parent of the child.
            pc->setParent(&newParent);
        }
    }
    
    /// Returns the number of immediate children of this node.
    /// @note Species and reactions and not included in this tally
    virtual size_t numberOfChildren() const override {return children().size();}

    /// Return the total number of nodes contained under this node's hieararchy
    /// @note This is a recursive call, and all nodes under this node are visited.
    virtual size_t countDescendents() const override {
        size_t sum = numberOfChildren();
        for (auto &c : children())
            sum+=c->numberOfChildren();
        return sum;
    }

    
    /// Returns the number of Species being managed by this node and its
    /// descendent nodes
    virtual size_t countSpecies() const override {
        size_t sum = 0;
        for (auto &c : children())
            sum+=c->countSpecies();
        return sum;
    }
    
    /// Returns the number of ReactionBase objects being managed by this node and its
    /// descendent nodes
    virtual size_t countReactions() const override {
        size_t sum = 0;
        for (auto &c : children())
            sum+=c->countReactions();
        return sum;
    }
    
    /// Returns a reference to the container of Component children of this node
    virtual vector<unique_ptr<Component>>& children () {return _children;}

    /// Returns a const reference to the container of Component children of this node
    virtual const vector<unique_ptr<Component>>& children () const {return _children;}

    /// Returns a pointer to the i-th Component child of this node
    virtual Component* children(size_t i) { return _children[i].get(); }
    virtual const Component* children(size_t i)const { return _children[i].get(); }
    
};

} // namespace medyan

#endif
