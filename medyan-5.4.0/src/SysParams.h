
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

#ifndef MEDYAN_SysParams_h
#define MEDYAN_SysParams_h

#include <filesystem>
#include <optional>
#include <string>
#include <vector>
#include <list>

#include "common.h"
#include "Mechanics/ForceField/Volume/Types.hpp"
#include "Parameter/Bubble.hpp"
#include "Parameter/Membrane.hpp"
#include "Parameter/Output.hpp"
#include "Structure/FilamentTypes.hpp"
#include "Util/Math/Vec.hpp"

namespace medyan {
//Did not minimize structure
#ifdef TRACKDIDNOTMINIMIZE
struct MinimizationParams{
    vector<floatingpoint> maxF;
    vector<floatingpoint> TotalE;
    vector<vector<floatingpoint>> Energyvec;
    vector<floatingpoint> Lambda;
    vector<floatingpoint> beta;
    vector<bool> safeModeORnot;
    vector<floatingpoint> tempEnergyvec;
    vector<vector<floatingpoint>> gradientvec;
};
#endif

//-------------------------------------
// Parameters read from the system input
//-------------------------------------

/// Struct to hold mechanical parameters for the system
struct MechParams {

    /// Struct to hold the ForceField types
    struct MechanicsFFType {
        
        //@{
        /// FilamentFF type
        string FStretchingType = "";
        string FBendingType    = "";
        string FTwistingType   = "";
        //@}
        
        //@{
        /// LinkerFF type
        string LStretchingType = "";
        string LBendingType    = "";
        string LTwistingType   = "";
        //@}
        
        //@{
        /// MotorFF type
        string MStretchingType = "";
        string MBendingType    = "";
        string MTwistingType   = "";
        //@}
        
        //@{
        /// BranchingFF type
        string BrStretchingType = "";
        string BrBendingType    = "";
        string BrDihedralType   = "";
        string BrPositionType   = "";
        //@}
        
        /// VolumeFF type
        CylinderVolumeExclusionFFType cylinderVolumeExclusionFFType = CylinderVolumeExclusionFFType::none;
        
        /// BoundaryFF Type
        string BoundaryFFType = "";
        
        /// Bubble Type
        string BubbleFFType = "";
        
        /// MTOC Type
        string MTOCFFType = "";
        
        /// AFM Type
        string AFMFFType = "";


        // MembraneFF type
        string memStretchingFFType     = "";
        string memTensionFFType        = "";
        string memBendingFFType        = "";
        string triangleBeadVolumeFFType = "";

        // Volume conservation ff type
        string volumeConservationFFType = "";
    };

    /// Struct to hold mechanics algorithm information
    struct MechanicsAlgorithm {
        string ConjugateGradient = "";
        
        /// Tolerance and cg parameters
        floatingpoint gradientTolerance = 1.0;
        floatingpoint energyChangeRelativeTolerance = 1e-9;
        floatingpoint maxDistance = 1.0;
        floatingpoint lambdaMax = 1.0;
        floatingpoint lambdarunningaverageprobability = 0.0;
        string linesearchalgorithm = "BACKTRACKING";
        bool tryToRecoverInLineSearchError = true;
        
        /// Not yet used
        string MD = "";

        // Other parameters for mechanics algorithm.
        bool adaptMeshDuringMinimization = false;
    };

    // Struct to hold membrane curvature mismatch parameters.
    struct ProteinCurvatureMismatchSetup {
        std::string speciesName;
        double      kBending = 0;
        double      eqCurv = 0;
    };

    struct TriangleBeadVolume {
        floatingpoint k          = 0;
        floatingpoint cutoff     = 0;
        floatingpoint cutoffMech = 0;
    };

    MechanicsFFType       mechanicsFFType;
    MechanicsAlgorithm    mechanicsAlgorithm;

    //@{
    /// Filament parameter
    vector<floatingpoint> FStretchingK  = {};
    vector<floatingpoint> FBendingK     = {};
    vector<floatingpoint> FBendingTheta = {};
    vector<floatingpoint> FTwistingK    = {};
    vector<floatingpoint> FTwistingPhi  = {};
    //@}

    //@{
    /// Linker parameter
    vector<floatingpoint> LStretchingK  = {};
    vector<floatingpoint> LBendingK     = {};
    vector<floatingpoint> LBendingTheta = {};
    vector<floatingpoint> LTwistingK    = {};
    vector<floatingpoint> LTwistingPhi  = {};
    //@}

    //@{
    /// MotorGhost parameter
    vector<floatingpoint> MStretchingK  = {};
    vector<floatingpoint> MBendingK     = {};
    vector<floatingpoint> MBendingTheta = {};
    vector<floatingpoint> MTwistingK    = {};
    vector<floatingpoint> MTwistingPhi  = {};
    //@}

    //@{
    /// BranchingPoint parameter
    vector<floatingpoint> BrStretchingK  = {};
    vector<floatingpoint> BrStretchingL  = {};
    vector<floatingpoint> BrBendingK     = {};
    vector<floatingpoint> BrBendingTheta = {};
    vector<floatingpoint> BrDihedralK    = {};
    vector<floatingpoint> BrPositionK    = {};
    //@}

    //@{
    /// Volume parameter
    vector<floatingpoint> VolumeK                 = {};
    floatingpoint         VolumeCutoff            = 0.0;
    // If monomer volume exclusion is enabled, this value indicates the monomer sampling interval on the filament, indexed by filament type.
    std::vector<Size>     volumeExclusionMonomerInterval;
    TriangleBeadVolume    triangleBeadVolume;
    //@}

    //@{
    /// Bubble parameter
    vector<floatingpoint> BubbleK = {};
    vector<floatingpoint> BubbleRadius = {};
    vector<floatingpoint> BubbleScreenLength = {};
	vector<floatingpoint> MTOCBendingK = {};
    vector<floatingpoint> AFMBendingK = {};

	floatingpoint BubbleCutoff = 0.0;

    ///If using more than one bubble
    short numBubbleTypes = 1;
    //@}

    // Parameters for special force field bead move constraint.
    // Note: bead constraint restricts the range of bead movement of every minimization, using the 3D FENE potential.
    floatingpoint beadConstraintK = 1;
    floatingpoint beadConstraintRmax = 160;

    // Surface protein curvature mismatch parameters.
    bool curvatureMismatchEnableCurvatureSensing = true;
    bool curvatureMismatchEnableCurvatureGeneration = true;
    std::vector<ProteinCurvatureMismatchSetup> proteinCurvatureMismatchSetups;
    
    
    //@{
    /// SPECIAL MECHANICAL PROTOCOLS

    bool pinLowerBoundaryFilaments = false;
    floatingpoint pinFraction = 1.0; //test

    ///To pin filaments on boundary via an attractive potential
    bool pinBoundaryFilaments = false;
    floatingpoint pinDistance = 250; ///< 250nm pinning distance for now
    floatingpoint pinK = 0.0;       ///< Tethered stiffness
    floatingpoint pinTime = 0.0;    ///< Time at which to pin the filaments

    // pin initial filament beads, using pinK
    bool pinInitialFilamentBelowZ = false;
    floatingpoint pinInitialFilamentBelowZValue = 0.0;

    // Membrane equilibrium area expansion.
    struct MembraneEqAreaChange {
        // The equilibrium area of every membrane in the system is increased with this rate (in nm^2/s).
        FP rate = 0;
        // The minimum equilibrium area for every membrane in the system.
        FP minEqArea = 0;
    };
    std::optional<MembraneEqAreaChange> membraneEqAreaChange;

    struct MembraneEqVolumeChange {
        // The equilibrium volume of every membrane in the system is increased with this rate (in nm^3/s).
        FP rate = 0;
        // The minimum equilibrium volume for every membrane in the system.
        FP minEqVolume = 0;
    };
    std::optional<MembraneEqVolumeChange> membraneEqVolumeChange;
    //@}

    //vectorization
    floatingpoint *coord;
    vector<vector<bool>> speciesboundvec;
    floatingpoint* cylsqmagnitudevector;
    vector<int> ncylvec;
    vector<int>bsoffsetvec;
    
    // parameters controlling the calculation of the Hessian matrix
    bool hessTracking = false;
    bool eigenTracking = true;
    bool rockSnapBool = false;
    float hessDelta = 0.0001;
    bool denseEstimationBool = true;
    bool hessMatrixPrintBool = false;
    int hessSkip = 20;

    int cylThresh = 0;

    medyan::FilamentModel globalFilamentModel = medyan::FilamentModel::beadCylinder;
};

/// Struct to hold chemistry parameters for the system
struct ChemParams {

    /// Struct to hold chemistry algorithm information
    struct ChemistryAlgorithm {
        
        string algorithm = "";
        
        //@{
        /// User can specify total run time of the simulation, as well as
        /// frequency of snapshots, neighbor list updates and minimizations.
        floatingpoint runTime = 0.0;
        
        floatingpoint snapshotTime = 0.0;
        floatingpoint datadumpTime = 1000.0;
        
        floatingpoint minimizationTime = 0.0;
        floatingpoint neighborListTime = 0.0;
        floatingpoint initialSlowDownTime = 0.0;
        //@}
        
        //@{
        /// Can also specify a frequency in terms of number of chemical reaction steps
        /// Useful for smaller systems and debugging
        int runSteps = 0;
        
        int snapshotSteps = 0;
        
        int minimizationSteps = 0;
        int neighborListSteps = 0;
        //@}
    };

    /// Struct to hold chem setup information
    struct ChemistrySetup {
        
        std::filesystem::path inputFile; // relative path only
    };

    // Constant parameters.
    //---------------------------------
    // On the same filament, the minimum cylinder position distance that allows pairwise binding.
    static constexpr int minCylinderDistanceSameFilament = 2;

    ChemistryAlgorithm chemistryAlgorithm;
    ChemistrySetup     chemistrySetup;

    //@{
    /// Number of general species
    short numBulkSpecies = 0;
    short numDiffusingSpecies = 0;
    //@}

    /// Number of different filament types
    short numFilaments = 1;

    /// Number of binding sites per cylinder
    /// Vector corresponds to each filament type
    vector<short> numBindingSites = {};

    short maxbindingsitespercylinder = 0;

    //Bindingsites are stored as packed 32bit integers. To ensure that there is adequate
    // memory space to store the binding sites, we need to shift based on the maximum
    // number of binding sites per cylinder.
    uint32_t shiftbybits = 0;
    uint32_t maxStableIndex = 0;

    //@{
    ///Extra motor parameters
    /// Vector corresponds to each motor type
    vector<short> motorNumHeadsMin = {};
    vector<short> motorNumHeadsMax = {};


    vector<floatingpoint> motorStepSize = {};
    //@}

    /// Binding sites on filaments
    /// 2D Vector corresponds to each filament type
    vector<vector<short>> bindingSites = {};

    // Whether binding on the same filament is allowed.
    bool allowPairBindingOnSameFilament = true;

    //@{
    /// Positions of all bound molecules in species vectors
    /// Vector corresponds to each filament type
    vector<short> brancherBoundIndex = vector<short>(MAX_FILAMENT_TYPES);
    vector<short> linkerBoundIndex   = vector<short>(MAX_FILAMENT_TYPES);
    vector<short> motorBoundIndex    = vector<short>(MAX_FILAMENT_TYPES);

    vector<vector<short>> bindingIndices = vector<vector<short>>(MAX_FILAMENT_TYPES);
    //@}
    
    
    //@{
    /// SPECIAL CHEMICAL PROTOCOLS

    /// Option to make Filament objects static after a certain time.
    /// This passivates any polymerization or depolymerization
    /// reactions, resulting in constant length filaments for the rest of simulation.
    bool makeFilamentsStatic = false;
    floatingpoint makeFilamentsStaticTime = 0.0;

    /// Option to make Linker objects static after a certain time.
    /// This passivates any unbinding reactions, resulting in permanently
    /// bound linkers for the rest of simulation.
    bool makeLinkersStatic = false;

    floatingpoint makeLinkersStaticTime = 0.0;

    bool dissTracking = false;
    bool eventTracking = false;
    int linkerbindingskip = 2;
    
    
    /// Make (de)polymerization depends on stress
    bool makeRateDepend = false;
    double makeRateDependTime = 0.0;
    double makeRateDependForce = 0.0;
    
    /// Make (de)polymerization depends on stress
    bool makeAFM = false;
    double AFMStep1 = 0.0;
    double AFMStep2 = 0.0;
    double IterChange = 0.0;
    double StepTotal = 0.0;
    double StepTime = 0.0;
    float originalPolyPlusRate;
    

    //@}
#ifdef CUDAACCL_NL
    string bindingmanagerlist = "";
    vector<floatingpoint> bmanagerdistances = {};
#endif
};

/// Struct to hold geometry parameters for the system
struct GeoParams {

    //@{
    /// Geometry parameter
    short nDim = 3;

    int NX = 0;
    int NY = 0;
    int NZ = 0;

    floatingpoint compartmentSizeX = 0;
    floatingpoint compartmentSizeY = 0;
    floatingpoint compartmentSizeZ = 0;

    // When choosing coordinates randomly in the entire compartment grid, this parameter limits the available grid domain with a certain fraction in each dimension.
    std::array<std::array<floatingpoint, 3>, 2> fracGridSpan = {{ { 0, 0, 0 }, { 1, 1, 1 } }};

    floatingpoint largestCompartmentSide = 0;

    floatingpoint largestCylinderSize = 0;

    vector<floatingpoint> monomerSize = {};

    ///Number of monomers in a cylinder
    vector<int> cylinderNumMon = {};

    vector<floatingpoint> cylinderSize = {};

    vector<floatingpoint> minCylinderSize = {};

    /// Minimum monomer length of a cylinder is preset
    int minCylinderNumMon = 3;
    //@}

    // Surface geometry settings.
    SurfaceGeometrySettings surfaceGeometrySettings {};
    MeshAdapterSettings     meshAdapterSettings {};
};

/// Struct to hold Boundary parameters for the system
struct BoundParams {

    /// Struct to hold the parameters of the Boundary
    struct BoundaryType {
        
        string boundaryShape = "";
        vector<string> boundaryMove = {};
        //string scaleDiffusion = "";
    };

    BoundaryType  boundaryType;

    //@{
    /// Mechanical parameter
    floatingpoint BoundaryK = 0;
    floatingpoint BScreenLength = 0;
    //@}

    /// Cutoff for force calculation
    floatingpoint BoundaryCutoff = 0;
    floatingpoint diameter = 0;

    /// Moving speed (if any)
    floatingpoint moveSpeed = 0;

    //@{
    /// Moving times
    floatingpoint moveStartTime = 0;
    floatingpoint moveEndTime = 0;
    //@}
    int transfershareaxis=-1;       ///Axis along which activate/deactivate protocols should be executed.
    int planestomove = -1; //tracks if both (2), left/bottom/front (1), or
    // right/top/back(0) planes should be moved.

};

/// Struct to hold dynamic rate changing parameters
struct DyRateParams {

    ///Struct to hold dynamic rate changer type
    struct DynamicRateType {
        
        ///Polymerization rate changing
        vector<string> dFPolymerizationType = {};
        
        ///Linker rate changing
        vector<string> dLUnbindingType = {};
        
        //@{
        ///Motor rate changing
        vector<string> dMUnbindingType = {};
        vector<string> dMWalkingType = {};

        //Qin----
        vector<string> dBUnbindingType = {};
        //@}
    };

    DynamicRateType dynamicRateType;

    /// Option for dynamic polymerization rate of filaments
    vector<floatingpoint> dFilPolymerizationCharLength = {};

    /// Option for dynamic unbinding rate of linkers
    vector<floatingpoint> dLinkerUnbindingCharLength = {};
    vector<floatingpoint> dLinkerUnbindingAmplitude = {};

    /// Option for dynamic unbinding rate of motors
    vector<floatingpoint> dMotorUnbindingCharForce = {};

    /// Option for dynamic walking rate of motors
    vector<floatingpoint> dMotorWalkingCharForce = {};

    //Qin
    /// Option for dynamic branching point unbinding rate
    vector<floatingpoint> dBranchUnbindingCharLength = {};

    /// Option for addinh dynamics branching point unbinding rate based on a
    // characteristic Force.
    vector<floatingpoint> dBranchUnbindingCharForce = {};

    /// Option for manual (de)polymerization rate changer
    /// Start time
    double manualCharStartTime = 100000.0;
    /// Plusend Polymerization Rate Ratio
    float manualPlusPolyRate = 1.0;
    /// Plusend Depolymerization Rate Ratio
    float manualPlusDepolyRate = 1.0;
    /// Minusend Polymerization Rate Ratio
    float manualMinusPolyRate = 1.0;
    /// Minusend Depolymerization Rate Ratio
    float manualMinusDepolyRate = 1.0;
};

/// Struct to hold Filament setup information
struct FilamentSetup {
    
    std::filesystem::path inputFile;
    
    ///If want a random distribution, used if inputFile is left blank
    int numFilaments = 0;
    ///Filament length, in number of cylinders
    int filamentLength = 1;
    ///Filament type to create
    short filamentType = 0;
    ///Filament projection type.
    string projectionType="STRAIGHT";

    bool USECHEMCOPYNUM = false; // if set to 0, restart file copy numbers are used. If
        // not, chemistry file copy numbers are used.

    ///For resetting pin positions in restart phase
    string pinRestartFile = "";

    // For cylindrical boundary MTOC simulation only.
    bool restrictCylinderBoundaryAngle = false;
    bool formRings = false;
};





//-------------------------------------
// Parameters read from other input files
//-------------------------------------

/// Struct to hold Species and Reaction information
/// @note - all filament-related reactions and species are 2D vectors corresponding
///         to the filament type specified in the input file.
struct ChemistryData {

    // Species information.
    struct SpeciesGeneralInfo {
        std::string   name;
        std::string   rspeciesType;
    };
    struct SpeciesBulkInfo {
        std::string   name;
        int           initialCopyNumber = 0;
        floatingpoint releaseTime = 0;
        floatingpoint removalTime = 0;
        std::string   rspeciesType;
        // Move boundary specific.
        std::string   copyNumberManipulationType = "NONE";
        // Move boundary specific.
        floatingpoint holdMolarity = 0;
    };
    struct SpeciesDiffusingInfo {
        std::string   name;
        int           initialCopyNumber = 0;
        floatingpoint diffusionCoefficient = 0;
        floatingpoint releaseTime = 0;
        // If removal time is zero, it will not be removed.
        floatingpoint removalTime = 0;
        std::string   rspeciesType;
        // For AVG type of rspecies only.
        int           numEvents = 0;
        // Move boundary specific.
        std::string   copyNumberManipulationType = "NONE";
        // Move boundary specific.
        floatingpoint holdMolarity = 0;
    };
    struct SpeciesMembraneDiffusingInfo {
        std::string name;
        floatingpoint diffusionCoefficient = 0;
        floatingpoint area = 0;
    };

    // Reaction information.
    struct ReactionSurfaceInfo {
        std::vector< std::string > reactants;
        std::vector< std::string > products;
        floatingpoint rate = 0;
    };
    struct ReactionAdsorptionDesorptionInfo {
        std::string speciesName3D;
        std::string speciesName2D;
        // Actual propensity of compartment-based adsorption is
        //   p = (on rate) * (density of 3D diffusing species) * (area of accessible lipid)
        double      onRate = 0;
        // Actual propensity of desorption is
        //   p = (off rate) * (number of 2D diffusing species)
        double      offRate = 0;
    };
    struct ReactionPolymerizationInfo {
        std::string speciesReactantDiffusingBulk;
        std::string speciesReactantPlusEndMinusEnd;
        std::string speciesProductFilament;
        std::string speciesProductPlusEndMinusEnd;
        floatingpoint rate = 0;
        // For dissipation tracking.
        floatingpoint gnum = 0;
        std::string hrcdid;
    };
    struct ReactionDepolymerizationInfo {
        std::string speciesReactantFilament;
        std::string speciesReactantPlusEndMinusEnd;
        std::string speciesProductDiffusingBulk;
        std::string speciesProductPlusEndMinusEnd;
        floatingpoint rate = 0;
        // For dissipation tracking.
        floatingpoint gnum = 0;
        std::string hrcdid;
    };
    struct ReactionAgingInfo {
        std::string speciesReactant;
        std::string speciesProduct;
        floatingpoint rate = 0;
        // For dissipation tracking.
        floatingpoint gnum = 0;
        std::string hrcdid;
    };
    struct ReactionDestructionInfo {
        std::string speciesReactantPlusEnd;
        std::string speciesReactantMinusEnd;
        std::string speciesProduct1;
        std::string speciesProduct2;
        floatingpoint rate = 0;
    };
    struct ReactionLinkerBindingInfo {
        struct ReactantInfo {
            std::string speciesBound1;
            std::string speciesBound2;
            // linkerDiffusing can be a bulk or diffusing species name.
            std::string linkerDiffusing;
        };
        struct ProductInfo {
            std::string linkerBound1;
            std::string linkerBound2;
        };

        ReactantInfo  reactantInfo;
        ProductInfo   productInfo;
        floatingpoint onRate = 0;
        floatingpoint offRate = 0;
        floatingpoint rMin = 0;
        floatingpoint rMax = 0;
        // For dissipation tracking.
        floatingpoint gnum = 0;
        std::string   hrcdid;
    };
    struct ReactionMotorWalkingInfo {
        std::string   speciesReactantMotor;
        std::string   speciesReactantEmptySite;
        std::string   speciesProductMotor;
        std::string   speciesProductEmptySite;
        floatingpoint rate = 0;
        // For dissipation tracking.
        floatingpoint gnum = 0;
        std::string   hrcdid;
    };

    /// Reaction happening between SpeciesBulk and SpeciesDiffusing ONLY
    vector<tuple<vector<string>, vector<string>, floatingpoint, floatingpoint, string>> genReactions = {};
    // Reactions on membrane surfaces.
    std::vector< ReactionSurfaceInfo > reactionsSurface;

    /// Reaction happening between SpeciesBulk ONLY
    vector<tuple<vector<string>, vector<string>, floatingpoint>> bulkReactions = {};

    // Surface adsorption / desorption.
    std::vector< ReactionAdsorptionDesorptionInfo > reactionsAdsorptionDesorption;

    /// Filament nucleation reaction
    vector<vector<tuple<vector<string>, vector<string>, floatingpoint>>> nucleationReactions;

    //@{
    /// Filament reactions
    /*!
        *  All Filament reactions are held using a vector containing a tuple with the 
        *  string of reactants, string of products, and the reaction rate.
        */
    /// Polymerization reactions
    std::vector<std::vector<ReactionPolymerizationInfo>> polymerizationReactions;
    /// Depolymerization reactions
    std::vector<std::vector<ReactionDepolymerizationInfo>> depolymerizationReactions;
    /// Aging reactions
    std::vector<std::vector<ReactionAgingInfo>> agingReactions;
    /// Destruction reactions
    std::vector<std::vector<ReactionDestructionInfo>> destructionReactions;

    /// Branching reactions
    /// This reaction also contains the off rate, and a string
    /// specifying the nucleation zone and relevant distance parameter
    vector<vector<tuple<vector<string>, vector<string>, floatingpoint, floatingpoint, string, floatingpoint>>> branchingReactions;

    /// Severing reactions
    vector<vector<tuple<string, floatingpoint>>> severingReactions;
    //@}

    //@{
    /// Cross Filament binding and unbinding reactions
    /*!
        *  All cross Filament reactions are held using a vector containing a tuple with 
        *  the string of reactants, string of products, the reaction rate, and binding 
        *  range.
        */
    /// Linker reactions
    std::vector< ReactionLinkerBindingInfo > linkerReactions;
    /// MotorGhost reactions
    std::vector< ReactionLinkerBindingInfo > motorReactions;
    //@}

    /// MotorGhost walking reactions
    std::vector< std::vector< ReactionMotorWalkingInfo > > motorWalkingReactions;

    //--------------------------------------------------------------------------
    // Species information.
    //--------------------------------------------------------------------------

    std::vector< SpeciesGeneralInfo > speciesGeneral;

    /// SpeciesBulk parsed, in the form of a tuple which contains the name and
    /// initial copy number, release time, removal time, CONST/REG qualifier, TARGET TYPE
    /// (TOTCONC/FREECONC) and Target CONCENTRATION (needed in move boundary)
    std::vector< SpeciesBulkInfo > speciesBulk;


    /// SpeicesDiffusing parsed, in the form of a tuple which contains name,
    /// initial copy number in reaction volume, the rate of diffusion, release time,
    /// removal time, AVG/REG qualifier, and number of events to average if applicable.
    std::vector< SpeciesDiffusingInfo > speciesDiffusing;
    // 2D diffusing species on membrane surfaces.
    std::vector< SpeciesMembraneDiffusingInfo > speciesMembraneDiffusing;

    //@{
    /// Filament species parsed
    vector<vector<string>> speciesFilament;
    vector<vector<string>> speciesPlusEnd;
    vector<vector<string>> speciesMinusEnd;
    
    vector<vector<string>> speciesBound;
    vector<vector<string>> speciesLinker;
    vector<vector<string>> speciesMotor;
    vector<vector<string>> speciesBrancher;
    //@}
    
    //@{
    /// Binding sites parsed
    vector<string> B_BINDING_INDEX;
    vector<string> L_BINDING_INDEX;
    vector<string> M_BINDING_INDEX;
    //@}
    
    ///Constructor initializes memory of vector members
    ChemistryData()
    
    : nucleationReactions(MAX_FILAMENT_TYPES),
      polymerizationReactions(MAX_FILAMENT_TYPES),
      depolymerizationReactions(MAX_FILAMENT_TYPES),
      agingReactions(MAX_FILAMENT_TYPES),
      destructionReactions(MAX_FILAMENT_TYPES),
    
      branchingReactions(MAX_FILAMENT_TYPES),
      severingReactions(MAX_FILAMENT_TYPES),
      motorWalkingReactions(MAX_FILAMENT_TYPES),
    
      speciesFilament(MAX_FILAMENT_TYPES),
      speciesPlusEnd(MAX_FILAMENT_TYPES),
      speciesMinusEnd(MAX_FILAMENT_TYPES),
      speciesBound(MAX_FILAMENT_TYPES),
      speciesLinker(MAX_FILAMENT_TYPES),
      speciesMotor(MAX_FILAMENT_TYPES),
      speciesBrancher(MAX_FILAMENT_TYPES),
    
      B_BINDING_INDEX(MAX_FILAMENT_TYPES),
      L_BINDING_INDEX(MAX_FILAMENT_TYPES),
      M_BINDING_INDEX(MAX_FILAMENT_TYPES) {}


    // Some auxiliary functions for convenient access.
    //--------------------------------------------------------------------------
    auto numLinkerSpecies() const noexcept { return linkerReactions.size(); }
    auto numMotorSpecies()  const noexcept { return motorReactions.size(); }

    // Find the species index for a given general species name.
    std::optional<Index> findIndexSpeciesGeneral(std::string_view name) const noexcept {
        for(Index i = 0; i < speciesGeneral.size(); ++i) {
            if(speciesGeneral[i].name == name) {
                return i;
            }
        }
        return std::nullopt;
    }
    // Find the species index for a given bulk species name.
    std::optional<Index> findIndexSpeciesBulk(std::string_view name) const noexcept {
        for(Index i = 0; i < speciesBulk.size(); ++i) {
            if(speciesBulk[i].name == name) {
                return i;
            }
        }
        return std::nullopt;
    }
    // Find the species index for a given diffusing species name.
    std::optional<Index> findIndexSpeciesDiffusing(std::string_view name) const noexcept {
        for (Index i = 0; i < speciesDiffusing.size(); ++i) {
            if (speciesDiffusing[i].name == name) {
                return i;
            }
        }
        return std::nullopt;
    }
    // Find the species index for a given membrane diffusing species name.
    std::optional<Index> findIndexSpeciesMembraneDiffusing(std::string_view name) const noexcept {
        for (Index i = 0; i < speciesMembraneDiffusing.size(); ++i) {
            if (speciesMembraneDiffusing[i].name == name) {
                return i;
            }
        }
        return std::nullopt;
    }
};


// Postprocessed membrane mesh chemistry data.
struct MembraneMeshChemistryInfo {
    // Note:
    //   - indices in this struct must correspond to the actual index as is in
    //     the diffusing species names vector

    struct DiffusionInfo {
        int      speciesIndex = 0;

        // Diffusion coeffecient, with dimension L^2 T^-1
        double   diffusionCoeff = 0;

        // Projected surface area of the protein, with dimension L^2.
        double   area = 0;
    };
    struct InternalReactionInfo {
        std::vector< int > reactantSpeciesIndices;
        std::vector< int > productSpeciesIndices;

        // Rate constant, with dimension L^(2*(numReactants - 1)) T^-1
        double rateConstant = 0;
    };

    // Factory that creates mesh chemistry information from ChemistryData.
    static MembraneMeshChemistryInfo fromChemistryData(const ChemistryData& chem) {
        using namespace std;

        MembraneMeshChemistryInfo mmci;

        // Diffusing species and diffusion reactions.
        for(int i = 0; i < chem.speciesMembraneDiffusing.size(); ++i) {
            auto& sinfo = chem.speciesMembraneDiffusing[i];
            mmci.diffusingSpeciesNames.push_back(sinfo.name);
            mmci.diffusion.push_back({
                i,                          // Species index.
                sinfo.diffusionCoefficient, // Diffusion coefficient.
                sinfo.area,                 // Projected surface area.
            });
        }

        // Internal reactions.
        for(auto& rinfo : chem.reactionsSurface) {
            InternalReactionInfo iri;

            const auto findSpeciesIndex = [&](const std::string& sname) {
                const int index = std::find(mmci.diffusingSpeciesNames.begin(), mmci.diffusingSpeciesNames.end(), sname) - mmci.diffusingSpeciesNames.begin();
                if(index >= mmci.diffusingSpeciesNames.size()) {
                    throw runtime_error("Invalid index.");
                }
                return index;
            };

            for(auto& sname : rinfo.reactants) {
                iri.reactantSpeciesIndices.push_back(findSpeciesIndex(sname));
            }
            for(auto& sname : rinfo.products) {
                iri.productSpeciesIndices.push_back(findSpeciesIndex(sname));
            }
            iri.rateConstant = rinfo.rate;

            mmci.internalReactions.push_back(move(iri));
        }

        return mmci;
    }

    // diffusing species
    std::vector< std::string > diffusingSpeciesNames;

    // diffusion reactions
    std::vector< DiffusionInfo > diffusion;

    // Internal reactions involving species
    std::vector< InternalReactionInfo > internalReactions;
};


struct BubbleData {
    struct BubbleInfo {
        int                   type = 0;
        Vec<3, floatingpoint> coord {};
    };

    std::vector<BubbleInfo> bubbles;
};


// Filament initialization data.
struct FilamentData {
    struct Filament {
        int type = 0;
        // Coordinate information to determine initial filament shape.
        std::vector<Vec<3, FP>> coords;
    };

    std::vector<Filament> filaments;
};

// Definition of simulation configuration.
// Note: This data structure should be able to be bijectively mapped to the input file.
struct SimulConfig {

    // The MetaParams class, used to store the source of config file, etc
    struct MetaParams {
        std::filesystem::path systemInputFile;

        // Directory of other input files.
        std::filesystem::path inputDirectory;

        // Whether the system is a restart instead of a fresh start.
        bool isRestart = false;
    };

    MetaParams     metaParams;

    // Parameters from the system input
    GeoParams      geoParams;
    MechParams     mechParams;
    ChemParams     chemParams;
    BoundParams    boundParams;
    DyRateParams   dyRateParams;
    BubbleSetup    bubbleSetup;
    MTOCSettings   mtocSettings;
    AFMSettings    afmSettings;
    MembraneSettings membraneSettings;
    FilamentSetup  filamentSetup;
    OutputParams   outputParams;

    // Parameters from other inputs
    ChemistryData  chemistryData;
    BubbleData     bubbleData;
    FilamentData   filamentData;

    // Postprocessed parameters.
    MembraneMeshChemistryInfo membraneMeshChemistryInfo;
};


//-------------------------------------
// Parameters read from the command line.
//-------------------------------------

enum class MedyanRunMode {
    // Normal simulation.
    simulation,
    // Run GUI only, without any simulation.
    gui,
    // (Deprecated) translate the trajectory file to pdb format.
    analyze,
    // Interactive configuration generation.
    config,
    // Run all tests.
    test,
    // Side simulations and routines.
    side,
};

struct RunAnalyzeParams {
    // Snapshot trajectory analysis specific.
    std::size_t analyzeMembraneBondFrame = 0;
    bool        analyzeMembraneBondAllFrames = false; // If it is set to true, it overrides the specific bond frame info.
    int analyzeFrameInterval = 1; // must be >= 1
};

// All parameters coming from the command line.
// These are not stored in SimulConfig because they do not belong to the input files.
struct CommandLineConfig {
    // Run mode
    MedyanRunMode runMode = MedyanRunMode::simulation;

    // File and directories
    std::filesystem::path inputFile;
    std::filesystem::path inputDirectory;
    std::filesystem::path outputDirectory = std::filesystem::current_path();

    // GUI
    bool guiEnabled = false;

    // Floating point environment.
    bool trapInvalidFP = false;

    // Threadings
    int numThreads = -1;

    // RNG
    bool rngSeedFixed = false;
    unsigned long long rngSeed = 0;

    // side procedure name
    std::string sideProcName;

    // Run mode analyze specific.
    RunAnalyzeParams runAnalyzeParams;
};


//--------------------------------------

// After reading from the input file, some post-processing is needed to compute a few more parameters.
void postprocess(SimulConfig&);



// Forward declaration.
class Controller;
/// Static class that holds all simulation parameters,
/// initialized by the SystemParser
class SysParams {
friend class Controller;
friend class SubSystem;
friend class Cylinder;

public:

    inline static medyan::MechParams MParams;    ///< The mechanical parameters
    static ChemParams CParams;    ///< The chemistry parameters
    static GeoParams GParams;     ///< The geometry parameters
    static BoundParams BParams;   ///< The boundary parameters
    static DyRateParams DRParams; ///< The dynamic rate parameters
    #ifdef TRACKDIDNOTMINIMIZE
    static MinimizationParams MinParams;
	#endif
    
    //@{
#ifdef NLSTENCILLIST
    static short numcylcylNL;//Tracks total number of neighborlists in the system
#endif
    //counter to check excluded volume.
    static int exvolcounter[3]; //positions 0, 1, and 2 correspond to parallel,
    // in-plane and regular cases.
    static long long exvolcounterz[3];

    inline static FilamentSetup filamentSetup;
    static bool RUNSTATE; //0 refers to restart phase and 1 refers to run phase.
    static bool INITIALIZEDSTATUS; // true refers to sucessful initialization. false
    static bool DURINGCHEMISTRY; //true if MEDYAN is running chemistry, false otherwise.
    // corresponds to an on-going initialization state.

    //aravind July11,2016
    static vector<float> BUBBareRate;
    //@
    static const auto& Mechanics() {return MParams;}
    static const ChemParams& Chemistry() {return CParams;}
    static const GeoParams& Geometry() {return GParams;}
    static const BoundParams& Boundaries() {return BParams;}
    static const DyRateParams& DynamicRates() {return DRParams;}


    #ifdef TRACKDIDNOTMINIMIZE
    static MinimizationParams& Mininimization() { return MinParams;}
	#endif

};

} // namespace medyan

#endif
