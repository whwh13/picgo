
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

#ifndef MEDYAN_FilamentInitializer_h
#define MEDYAN_FilamentInitializer_h

#include "common.h"
#include "SysParams.h"

namespace medyan {

// Forward declarations.
class SubSystem;
template< typename MemType > class MembraneRegion;
class Membrane;

/// An implementation of FilamentInitialzer that creates a completely random
/// Filament distribution within the specified boundary
FilamentData createFilamentsRandomDist(
    SubSystem&                      sys,
    const MembraneRegion<Membrane>& mr,
    int                             numFilaments,
    int                             filamentType,
    int                             lenFilaments,
    const FilamentSetup&            filamentSetup
);

/// An implementation of FilamentInitialzer that creates a sufficiently spaced
/// network of filaments for investigation of small numbers of filaments.
FilamentData createFilamentsConnectedDist(
    SubSystem&                      sys,
    const MembraneRegion<Membrane>& mr,
    int                             numFilaments,
    int                             filamentType,
    int                             lenFilaments,
    const SimulConfig&              conf
);

/// An implementation of FilamentInitialzer that creates a random MTOC configuration
FilamentData createFilamentsMTOCDist(
    SubSystem&                      sys,
    const MTOCInit&                 mtocInit,
    floatingpoint                   bubbleRadius,
    const GeoParams&                geoParams
);

/// An implementation of FilamentInitialzer that creates a random AFM configuration
FilamentData createFilamentsAFMDist(
    SubSystem&                      sys,
    const AFMInit&                  afmInit,
    floatingpoint                   bubbleRadius,
    const GeoParams&                geoParams
);

} // namespace medyan

#endif
