#ifndef MEDYAN_Structure_SurfaceMesh_MembraneHierarchy_Hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneHierarchy_Hpp

#include <algorithm> // remove
#include <iostream>
#include <stdexcept> // runtime_error
#include <string>
#include <utility> // move

#include "common.h" // floatingpoint
#include "Component.h"
#include "Composite.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Util/Math/Vec.hpp"
#include "Util/StableVector.hpp"

/******************************************************************************
The closed membranes are not allowed to intersect themselves or each other.
Under this assumption, the orientable closed membranes can have a clear
hierarchy of containing relationship. The 3d space can be divided into 2
regions per membrane regardless of genus, and this hierarchy can clearly
indicate the boundaries of such regions.

This header provides a class for containing hierarchy of closed membranes
represented by a tree structure. The parent node has the membrane directly
containing the membranes of the children; the membranes of the nodes of the
same level does not contain each other.
******************************************************************************/

namespace medyan {
template< typename MemType >
class MembraneHierarchy: public Composite {

public:
    using MembraneIndexType = typename medyan::StableVector<MemType>::Index;
private:

    // -1 is used as sentinel value, indicating no membrane is attached.
    MembraneIndexType mi_ {-1};

    // helper function for printSelf
    void printTree(std::string indent, bool last) const {
        std::cout << indent;
        if (last) {
            std::cout << "\\-";
            indent += "  ";
        }
        else {
            std::cout << "|-";
            indent += "| ";
        }
        
        std::cout << this;
        if(mi_.value != -1)
            std::cout << " (Mem index: " << mi_.value << ")";
        else
            std::cout << " (No membrane attached)";
        
        std::cout << std::endl;

        const auto n = numberOfChildren();
        for (size_t idx = 0; idx < n; ++idx)
            static_cast< const MembraneHierarchy* >(children(idx))->printTree(indent, idx == n - 1);
    }

public:

    /**************************************************************************
    Ctors and Dtors
    **************************************************************************/
    MembraneHierarchy(MembraneIndexType mi): mi_(mi) {}

    /**************************************************************************
    Getters and Setters
    **************************************************************************/
    auto getMembraneIndex() const { return mi_; }

    /**************************************************************************
    Implements Component
    **************************************************************************/
    virtual int getType()override { return 0; }
    virtual void printSelf()const override {
        std::cout << "\nMembraneHierarchy: ptr = " << this << std::endl;

        std::cout << "\nTree structure:\n";
        printTree("", true);

        std::cout << std::endl;
    }

    /**************************************************************************
    Static root
    **************************************************************************/
    static MembraneHierarchy& root() { static MembraneHierarchy root(MembraneIndexType{ -1 }); return root; }

    /**************************************************************************
    Operations on a tree structure
    **************************************************************************/
    // When new membrane is inserted.
    // This function requires that geometry of the membrane has been updated.
    template< typename Context >
    static void addMembrane(Context& sys, MembraneIndexType mi, MembraneHierarchy& root) {
        using MT = typename MemType::MeshType;

        auto& rootChildren = root.children();
        auto& m = sys.membranes[mi];

        // Pick a point on the membrane and check the containing relationship with other membranes
        const medyan::Vec< 3, floatingpoint > p (m.getMesh().attribute(typename MT::VertexIndex {0}).getCoordinate(sys));

        // Recursively search
        for(auto& childPtr: rootChildren) {
            
            MembraneHierarchy* hiePtr = static_cast< MembraneHierarchy* >(childPtr.get());

            if(hiePtr->mi_.value == -1) {
                hiePtr->printSelf();
                throw std::runtime_error("The child node does not point to a specific membrane.");
            }

            if(sys.membranes[hiePtr->mi_].isClosed() && medyan::contains(sys, sys.membranes[hiePtr->mi_].getMesh(), p)) { // Is inside one of the child nodes
                return MembraneHierarchy::addMembrane(sys, mi, *hiePtr); // Search that child
            }
        }

        // Now, the new membrane is outside of every child membrane.

        // First create a new node
        auto newNode = std::make_unique< MembraneHierarchy >(mi);

        // Then check whether any children is inside this membrane.
        // Open membranes cannot contain any children.
        if(m.isClosed()) {

            for(auto& childPtr: rootChildren) {

                MembraneHierarchy* hiePtr = static_cast<MembraneHierarchy*>(childPtr.get());

                const medyan::Vec< 3, floatingpoint > hieP (sys.membranes[hiePtr->mi_].getMesh().attribute(typename MT::VertexIndex {0}).getCoordinate(sys));

                if(medyan::contains(sys, m.getMesh(), hieP)) { // The child membrane is inside new membrane

                    // Add child to the new node
                    newNode->addChild(std::move(childPtr));
                    // Now the content of childPtr is moved and becomes the new child of the new node.
                    // The parent of the child has also been changed.
                    // The childPtr now should be nullptr.
                }
            }

            // Then remove all null children from root using erase-remove idiom.
            // Good news: unique_ptr can be compared to nullptr_t.
            rootChildren.erase(std::remove(rootChildren.begin(), rootChildren.end(), nullptr), rootChildren.end());
        }

        // Finally add the new node to the current root
        root.addChild(std::move(newNode)); // Also manages the deletion of newNode

    } // void addMembrane(MemType* m, MembraneHierarchy& root)
    template< typename Context >
    static void addMembrane(Context& sys, MembraneIndexType mi) { addMembrane(sys, mi, root()); }

    // When a membrane is removed. Must be a closed membrane.
    // Returns whether something is deleted.
    static bool removeMembrane(MembraneIndexType mi, MembraneHierarchy& root) {
        // Find the membrane to be removed by recursive search
        // Only 1 node will be deleted

        MembraneHierarchy* nodeToBeDeleted = nullptr;
        
        for(auto& childPtr: root.children()) {
            
            MembraneHierarchy* hiePtr = static_cast<MembraneHierarchy*>(childPtr.get());

            if(hiePtr->mi_ == mi) { // Found the node to be removed
                nodeToBeDeleted = hiePtr;
            }
            else {
                if(MembraneHierarchy::removeMembrane(mi, *hiePtr)) return true;
                // else the search continues
            }
        }

        if(nodeToBeDeleted) {

            // Bring all its children under its parent.
            for(auto& childHiePtr: nodeToBeDeleted->children()) {
                root.addChild(move(childHiePtr));
            }

            // Remove the node from root and destroy it.
            root.removeChild(nodeToBeDeleted);

            // Exit function
            return true;

        } else {
            return false;
        }
    } // void removeMembrane(MemType* m, MembraneHierarchy& root)
    static bool removeMembrane(MembraneIndexType mi) { return removeMembrane(mi, root()); }

};

} // namespace medyan



#endif
