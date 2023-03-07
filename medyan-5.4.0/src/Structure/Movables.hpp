#ifndef MEDYAN_Structure_Movables_hpp
#define MEDYAN_Structure_Movables_hpp

#include "Chemistry/AdsorptionDesorption.hpp"
#include "Controller/GController.h"
#include "Structure/Bubble.h"
#include "Structure/Compartment.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/Movable.h"
#include "Structure/Special/AFM.h"
#include "Structure/Special/MTOC.h"
#include "Structure/SurfaceMesh/Edge.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "SysParams.h"
#include "utility.h"

namespace medyan::movable {

template< typename T >
constexpr bool isMovable =
    std::is_same_v< T, Vertex > ||
    std::is_same_v< T, Edge > ||
    std::is_same_v< T, Triangle > ||
    std::is_same_v< T, AFM > ||
    std::is_same_v< T, MTOC > ||
    std::is_base_of_v< Movable, T >;
template< typename T >
struct IsMovable : std::integral_constant< bool, isMovable<T> > {};

constexpr Overload updatePosition {
    // Future: replace with concrete types instead of "Movable".
    [](auto&& sys, Movable& m) {
        m.updatePosition();
    },

    [](auto&& sys, Vertex& v) {
        // Update compartment.
        if(sys.getCompartmentGrid()) {
            Compartment* newc = &sys.getCompartmentGrid()->getCompartment(v.coord);

            if(v.cellElement.manager) {
                // Update compartment if needed.
                Compartment* curc = v.cellElement.manager->getHead(v.cellElement);
                if(newc != curc) {
                    v.cellElement.manager->updateElement(v.cellElement, newc->vertexCell);

                    // Update adsorption/desorption reactions.
                    medyan::retarget3DSpeciesToCompartment(v.cVertex, *newc);
                }
            }
            else {
                // Initialize manager and register with compartment.
                v.cellElement.manager = newc->vertexCell.manager;
                v.cellElement.manager->addElement(v.sysIndex, v.cellElement, newc->vertexCell);
            }
        }
    },
    [](auto&& sys, Edge& e) {
        // Update coordinate.
        const auto& mesh = e.getParent(sys).getMesh();
        const auto hei0 = mesh.halfEdge(Membrane::MeshType::edgeIndex(e.getTopoIndex()));
        const auto hei1 = mesh.opposite(hei0);
        const auto v0 = mesh.target(hei0);
        const auto v1 = mesh.target(hei1);

        e.coordinate = 0.5 * (mesh.attribute(v0).getCoordinate(sys) + mesh.attribute(v1).getCoordinate(sys));

        // Update compartment.
        if(sys.getCompartmentGrid()) {
            Compartment* newc = &sys.getCompartmentGrid()->getCompartment(e.coordinate);

            if(e.cellElement.manager) {
                // Update compartment if needed.
                Compartment* curc = e.cellElement.manager->getHead(e.cellElement);
                if(newc != curc) {
                    e.cellElement.manager->updateElement(e.cellElement, newc->edgeCell);
                }
            }
            else {
                // Initialize manager and register with compartment.
                e.cellElement.manager = newc->edgeCell.manager;
                e.cellElement.manager->addElement(e.sysIndex, e.cellElement, newc->edgeCell);
            }
        }
    },
    [](auto&& sys, Triangle& t) {
        // Update coordinate.
        const auto& mesh = t.getParent(sys).getMesh();
        const auto hei0 = mesh.halfEdge(Membrane::MeshType::triangleIndex(t.getTopoIndex()));
        const auto hei1 = mesh.next(hei0);
        const auto hei2 = mesh.next(hei1);
        const auto v0 = mesh.target(hei0);
        const auto v1 = mesh.target(hei1);
        const auto v2 = mesh.target(hei2);

        t.coordinate = (
            mesh.attribute(v0).getCoordinate(sys)
            + mesh.attribute(v1).getCoordinate(sys)
            + mesh.attribute(v2).getCoordinate(sys)
        ) / 3;

        // Update compartment.
        if(sys.getCompartmentGrid()) {
            Compartment* newc = &sys.getCompartmentGrid()->getCompartment(t.coordinate);

            if(t.cellElement.manager) {
                // Update compartment if needed.
                Compartment* curc = t.cellElement.manager->getHead(t.cellElement);
                if(newc != curc) {
                    t.cellElement.manager->updateElement(t.cellElement, newc->triangleCell);
                }
            }
            else {
                // Initialize manager and register with compartment.
                t.cellElement.manager = newc->triangleCell.manager;
                t.cellElement.manager->addElement(t.sysIndex, t.cellElement, newc->triangleCell);
            }
        }
    },

    [](auto&& sys, AFM& obj) {
        obj.updatePosition(sys);
    },
    [](auto&& sys, MTOC& obj) {
        obj.updatePosition(sys);
    },
};


constexpr Overload unregister {
    // Future: replace with concrete types instead of "Movable".
    [](auto&& sys, Movable& m) {
    },

    [](auto&& sys, Vertex& v) {
        if(v.cellElement.manager) {
            v.cellElement.manager->removeElement(v.cellElement);
        }
    },
    [](auto&& sys, Edge& e) {
        if(e.cellElement.manager) {
            e.cellElement.manager->removeElement(e.cellElement);
        }
    },
    [](auto&& sys, Triangle& t) {
        if(t.cellElement.manager) {
            t.cellElement.manager->removeElement(t.cellElement);
        }
    },
    [](auto&& sys, AFM&) {},
    [](auto&& sys, MTOC&) {},
};



} // namespace medyan::movable

#endif
