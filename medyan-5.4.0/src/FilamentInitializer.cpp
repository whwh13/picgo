
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

#include "FilamentInitializer.h"

#include "Boundary.h"
#include "Bubble.h"
#include "Bead.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneRegion.hpp"
#include "SubSystem.h"

#include "MathFunctions.h"
#include "Controller/GController.h"
#include "SysParams.h"
#include "Rand.h"

namespace medyan {
using namespace mathfunc;


FilamentData createFilamentsRandomDist(
    SubSystem&                      sys,
    const MembraneRegion<Membrane>& mr,
    int                             numFilaments,
    int                             filamentType,
    int                             lenFilaments,
    const FilamentSetup&            filamentSetup
) {
    
    FilamentData fd;
    //Create random distribution of filaments
    int filamentCounter = 0;
    
    Boundary *b = mr.getBoundary();
    
    //Qin, if boundary shape is cylinder, create filament in the center of system and vertical to Z axis
    if(b && b->getShape() == BoundaryShape::Cylinder) {
        
        while (filamentCounter < numFilaments) {

            // Create a random filament vector one cylinder long.
            const auto p1 = sys.getCompartmentGrid()->getRandomCenterCoordinates();
            const auto direction = Rand::randUnitVector3();

            const auto p2 = p1 + (lenFilaments * SysParams::Geometry().cylinderSize[filamentType] - 0.01) * direction;

            //check filament orientation
            Vec<3, floatingpoint> midpoint { (p1[0] + p2[0])/2, (p1[1] + p2[1])/2, 0 };
            Vec<3, floatingpoint> rvec { midpoint[0] - SysParams::Boundaries().diameter/2, midpoint[1] - SysParams::Boundaries().diameter/2, 0 };
            Vec<3, floatingpoint> filvec = { p2[0] - p1[0], p2[1] - p1[1], 0 };
            
            auto rdotfil = dot(rvec,filvec);
            auto norm1 =  magnitude(rvec);
            auto norm2 =  magnitude(filvec);
            auto cosrfil = rdotfil / norm1 / norm2;

            if(filamentSetup.restrictCylinderBoundaryAngle) {
                if(abs(cosrfil) > 0.34) // allow angles between -20 and 20 degrees.
                    continue;
            }

            //check if these points are outside bubbles
            bool inBubble = false;
            for(auto& bb : sys.bubbles) {
                auto radius = bb.getRadius();

                if((distance2(bb.coord, p1) < radius * radius) ||
                   (distance2(bb.coord, p2) < radius * radius))
                    inBubble = true;
            }

            //check if within cutoff of boundary
            bool outsideCutoff = false;
            if(b->distance(vec2Vector(p1)) < SysParams::Boundaries().BoundaryCutoff / 4.0 ||
               b->distance(vec2Vector(p2)) < SysParams::Boundaries().BoundaryCutoff / 4.0) {
                outsideCutoff = true;
            }

            // Check if filaments need to form rings.
            if(filamentSetup.formRings) {
                if(b->sidedistance(vec2Vector(p1)) > GController::getSize()[0]/8 ||
                    b->sidedistance(vec2Vector(p2)) > GController::getSize()[0]/8) {
                    outsideCutoff = true;
                }
            }

            if(b->within(vec2Vector(p1)) && b->within(vec2Vector(p2)) && !inBubble && !outsideCutoff) {
                fd.filaments.push_back({ filamentType, { p1, p2 } });
                filamentCounter++;
            }
        }
    }

    else{
        while (filamentCounter < numFilaments) {

            //Create a random filament vector one cylinder long
            const auto p1 = sys.getCompartmentGrid()->getRandomCoordinates();
            const auto direction = Rand::randUnitVector3();

            const auto p2 = p1 + (lenFilaments * SysParams::Geometry().cylinderSize[filamentType] - 0.01) * direction;

            //check if these points are outside bubbles
            bool inBubble = false;
            for(auto& bb : sys.bubbles) {
                auto radius = bb.getRadius();
                
                if((distance2(bb.coord, p1) < radius * radius) ||
                   (distance2(bb.coord, p2) < radius * radius))
                    inBubble = true;
            }

            //check if within cutoff of boundary
            bool outsideCutoff = mr.getBoundary() && (
                mr.getBoundary()->distance(vec2Vector(p1)) < SysParams::Boundaries().BoundaryCutoff / 4.0 ||
                mr.getBoundary()->distance(vec2Vector(p2)) < SysParams::Boundaries().BoundaryCutoff / 4.0
            );
            
            if(mr.contains(sys, p1) && mr.contains(sys, p2) && !inBubble && !outsideCutoff) {
                fd.filaments.push_back({filamentType, { p1, p2 } });
                filamentCounter++;
            }
        }

    }

    return fd;

}

FilamentData createFilamentsConnectedDist(
    SubSystem&                      sys,
    const MembraneRegion<Membrane>& mr,
    int                             numFilaments,
    int                             filamentType,
    int                             lenFilaments,
    const SimulConfig&              conf
) {

    ///SET THIS SPACING PARAMETER
    const FP maxSpacing = 50;

    FilamentData fd;
    
    ///First filament as normal
    //Create a random filament vector one cylinder long
    Vec<3, FP> firstPoint = {500,1000,1000};
    
    Vec<3, FP> direction = {1, 0, 0};
    
    auto secondPoint = firstPoint + (lenFilaments * conf.geoParams.cylinderSize[filamentType] - 0.01) * direction;
    
    FP len = distance(firstPoint, secondPoint);
    fd.filaments.push_back({filamentType, { firstPoint, secondPoint } });
    
    auto prevFirstPoint = firstPoint;
    auto prevSecondPoint = secondPoint;
    
    const FP safeDist = conf.boundParams.BoundaryCutoff;
    
    ///now create properly distanced network
    int filamentCounter = 1;
    while (filamentCounter < numFilaments) {
    
        ///pick a random distance from a random point on the chain
        direction = normalizedVector(secondPoint - firstPoint);
        len = distance(firstPoint, secondPoint);
        const FP randomSeg = Rand::randfloatingpoint(0, len);
        
        const auto randomPoint = firstPoint + randomSeg * direction;
        
        //now pick another random point which is within a certain distance away
        auto randdir = Rand::randUnitVector3();        
        const FP randomDist = Rand::randfloatingpoint(0, maxSpacing);
        const auto nextRandomPoint = randomPoint + randomDist * randdir;
        
        //now pick another random direction for the next filament creation
        randdir = Rand::randUnitVector3();
    
        //random length spacing for new filament
        floatingpoint randomLengthSpacing = Rand::randfloatingpoint(0, lenFilaments * conf.geoParams.cylinderSize[filamentType] - 0.01);
        
        firstPoint = nextRandomPoint + randomLengthSpacing * randdir;
        
        //switch rand direction and create second point
        randdir = -randdir;
        
        secondPoint = firstPoint + (lenFilaments * conf.geoParams.cylinderSize[filamentType] - 0.01) * randdir;
        
        //choose if within boundary
        if(
            mr.contains(sys, firstPoint) && mr.contains(sys, secondPoint) && (
                !mr.getBoundary() || (
                    mr.getBoundary()->distance(vec2Vector(firstPoint)) > safeDist &&
                    mr.getBoundary()->distance(vec2Vector(secondPoint)) > safeDist
                )
            )
        ) {
            fd.filaments.push_back({filamentType, { firstPoint, secondPoint } });
            
            prevFirstPoint = firstPoint;
            prevSecondPoint = secondPoint;
            
            filamentCounter++;
        }
        else { //reset
            firstPoint = prevFirstPoint;
            secondPoint = prevSecondPoint;
        }
        
    }
    
    return fd;
}

FilamentData createFilamentsMTOCDist(
    SubSystem&                      sys,
    const MTOCInit&                 mtocInit,
    floatingpoint                   bubbleRadius,
    const GeoParams&                geoParams
) {
    using namespace std;

    FilamentData fd;
    
    auto theta1 = mtocInit.theta1;
    auto theta2 = mtocInit.theta2;
    auto phi1 = mtocInit.phi1;
    auto phi2 = mtocInit.phi2;

    auto pb = sys.getBoundary();

    for(int fi = 0; fi < mtocInit.numFilaments; ++fi) {
        

        floatingpoint l = Rand::randfloatingpoint(theta1 * 2 * M_PI,theta2 * 2 * M_PI);
        floatingpoint h = Rand::randfloatingpoint(phi1 * M_PI - M_PI/2, phi2 * M_PI - M_PI/2);
        
        
        Vec<3, floatingpoint> point1 {
            mtocInit.bubbleCoord[0] + bubbleRadius * cos(l) * cos(h),
            mtocInit.bubbleCoord[1] + bubbleRadius * sin(h),
        };
                
        // add restrictions to MTOC filament position in Cylinder boundary condition
        if(pb && pb->getShape() == BoundaryShape::Cylinder){
            point1[2] = mtocInit.bubbleCoord[2];
        }
        else{
            point1[2] = mtocInit.bubbleCoord[2] + bubbleRadius * sin(l) * cos(h);
        }
        
        // get projection outward from the MTOC
        auto dir = normalizedVector(point1 - mtocInit.bubbleCoord);
        const auto point2 = point1 + (geoParams.cylinderSize[mtocInit.filamentType]*mtocInit.numCylindersPerFilament - 0.01) * dir;
        
        fd.filaments.push_back({mtocInit.filamentType, { point1, point2 } });
    }
    
    return fd;
}
FilamentData createFilamentsAFMDist(
    SubSystem&                      sys,
    const AFMInit&                  afmInit,
    floatingpoint                   bubbleRadius,
    const GeoParams&                geoParams
) {
    using namespace std;

    FilamentData fd;
    
    auto theta1 = afmInit.theta1;
    auto theta2 = afmInit.theta2;
    auto phi1 = afmInit.phi1;
    auto phi2 = afmInit.phi2;

    auto pb = sys.getBoundary();

    for(int fi = 0; fi < afmInit.numFilaments; ++fi) {
        

        floatingpoint l = Rand::randfloatingpoint(theta1 * 2 * M_PI,theta2 * 2 * M_PI);
        floatingpoint h = Rand::randfloatingpoint(phi1 * M_PI - M_PI/2, phi2 * M_PI - M_PI/2);

        

        Vec<3, floatingpoint> point1 {
            afmInit.bubbleCoord[0] + bubbleRadius * cos(l) * cos(h),
            afmInit.bubbleCoord[1] + bubbleRadius * sin(h),
            afmInit.bubbleCoord[2] + bubbleRadius * sin(l) * cos(h),
        };

        
        // get projection outward from the AFM
        auto dir = normalizedVector(point1 - afmInit.bubbleCoord);
        const auto point2 = point1 + (geoParams.cylinderSize[afmInit.filamentType]*afmInit.numCylindersPerFilament - 0.01) * dir;
        
        fd.filaments.push_back({afmInit.filamentType, { point1, point2 } });
    }
    
    return fd;
}

} // namespace medyan
