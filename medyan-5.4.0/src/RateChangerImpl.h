
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

#ifndef MEDYAN_RateChangerImpl_h
#define MEDYAN_RateChangerImpl_h

#include "common.h"

#include "RateChanger.h"

#include "SysParams.h"

namespace medyan {
/// A brownian ratchet implementation of the FilamentRateChanger.
/// Used for filament polymerization when under load force.

/// @note - This function updates polymerization rates based on the
/// Elastic Brownian Ratchet Model (by Peskin et al, Biophys J 1993):
///
///                 k = k_0 * exp(-f * x / kT)
///
/// @note - We note that we have implemented a load force ceiling of 100pN
///         such that recalculated reaction rates are not excessively small.
///         This would produce problems in the chemical simulation algorithm.
///         A 100pN load force ensures that the polymerization rate produced
///         will be small enough such that polymerization events are VERY rare (factor = 1E-29).

class BrownianRatchet : public FilamentRateChanger {
    
private:
    floatingpoint _x; ///< The characteristic length for this function
    const floatingpoint _max_f = 100; ///< 100pN ceiling
    
public:
    BrownianRatchet(floatingpoint charLength) : _x(charLength) {}
    
    virtual float changeRate(float bareRate, floatingpoint force);

    virtual float getRateChangeFactor(floatingpoint force);
};

///A catch-slip bond implementation of the LinkerRateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form of Bell et al, 1978:
///
///  k = k_0 * (a_c * exp(-f * x_c / kT)  + a_s * exp(-f * x_s / kT))
///
///  where x and a are the characteristic lengths and amplitudes
///  of the catch and slip portions of the function, respectively.

class LinkerCatchSlip : public LinkerRateChanger {
    
private:
    floatingpoint _a1; ///< catch bond amplitude
    floatingpoint _a2; ///< slip bond amplitude
    floatingpoint _x1; ///< catch bond characteristic length
    floatingpoint _x2; ///< slip bond characteristic length
    
public:
    LinkerCatchSlip(short linkerType,
                    floatingpoint amplitude1, floatingpoint amplitude2,
                    floatingpoint charLength1, floatingpoint charLength2)
    
    : LinkerRateChanger(linkerType),
    _a1(amplitude1),  _a2(amplitude2),
    _x1(charLength1), _x2(charLength2) {}
    
    virtual float changeRate(float bareRate, floatingpoint force);

    virtual float getRateChangeFactor(floatingpoint force);
};

///A slip bond implementation of the LinkerRateChanger.
///Used for cross-linker unbinding when under stress.

/// @note - This function updates unbinding rates based on the
/// following exponential form of Bell et al, 1978:
///
///                 k = k_0 * exp(f * a / kT)
///
/// So as to exponetially increase the unbinding with more force.

class LinkerSlip : public LinkerRateChanger {
    
private:
    floatingpoint _x; ///< The characteristic length for this function
    
public:
    LinkerSlip(short linkerType, floatingpoint charLength)
    
    : LinkerRateChanger(linkerType), _x(charLength) {}
    
    virtual float changeRate(float bareRate, floatingpoint force);

    virtual float getRateChangeFactor(floatingpoint force);
};

///A slip bond implementation of the BrancherRateChanger.
///The same as crosslinker so far, details is being explored

class BranchSlip : public BranchRateChanger {
    
private:
    floatingpoint _x; ///< The characteristic length for this function
    
public:
    BranchSlip(short branchType, floatingpoint charLength)
    
    : BranchRateChanger(branchType), _x(charLength) {}
    
    virtual float changeRate(float bareRate, floatingpoint force);

    virtual float getRateChangeFactor(floatingpoint force);
};


///A slip bond implementation of the BrancherRateChanger.
/// Uses k = k0 * exp(-F/F0) to change unbinding rate k based on unbinding force at F=0,
// (k0) and characteristic unbinding force F0.
class BranchSlipF : public BranchRateChanger {

private:
	floatingpoint _F0; ///< The characteristic length for this function

public:
	BranchSlipF(short branchType, floatingpoint charForce)

			: BranchRateChanger(branchType), _F0(charForce) {}

	virtual float changeRate(float bareRate, floatingpoint force);

	virtual float getRateChangeFactor(floatingpoint force);
};

///A catch bond implementation of the MotorRateChanger
///Used for a motor unbinding when under stress
///Adopted from the results of Erdmann et al. 2013.

/// @note - This function updates unbinding rates of a
/// Myosin II ensemble based on the following exponential form:
///
///    k_unbinding,eff = k_0 * exp(-F / (N_b(F) * F_0))
///
/// where k_0 is the unbinding rate under zero load,
///
///    k_0 = k_on * (N_t) / (exp(log((k_on + k_off) / k_off) * N_t) - 1)
///
/// F_0 is the characteristic force defining this catch, and
///
///    N_b(F) = rho * N_t + beta * F / N_t
///
/// is the number of bound heads in the ensemble, beta is an emperical parameter
/// which determines the rate of increase for the number of bound heads with
/// respect to applied forces.

// Adding back PLOS version of Mechanochemical feedback based in Medyan3.0.
// TURN it on with macro PLOSFEEDBACK
/// @note - Assuming a duty ratio p = 0.1
/// @note - This function updates unbinding rates of a
/// Myosin II ensemble based on the following exponential form:
///
///      k_unbinding,eff = beta * (k_0 / N_b) * exp(-F / (N_b * F_0))
///
/// where k_0 is the unbinding rate under zero load,
/// F_0 is the characteristic force defining this catch,
/// beta has been chosen to be 0.2,
/// and N_b is the number of bound motor heads in the ensemble,
/// approximated by Erdmann et al. 2013 to be:
///
///             N_b = p * N_t + (F * gamma)
///
/// where gamma has been chosen to be 0.05
/// for a low duty ratio motor (p = 0.1).

class MotorCatch : public MotorRateChanger {
    
private:
    floatingpoint _F0;  ///< characteristic force
    
    //@{
    ///Constant parameters
    floatingpoint _dutyRatio;
    floatingpoint _beta;
    #ifdef PLOSFEEDBACK
    floatingpoint _gamma = 0.05;
    #endif
    //@}
    
public:
    MotorCatch(short motorType, floatingpoint charForce, floatingpoint dutyRatio, floatingpoint beta)
    
    : MotorRateChanger(motorType), _F0(charForce), _dutyRatio(dutyRatio), _beta(beta) {}
    
    /// Set the number of bound heads based on force
    virtual float numBoundHeads(float onRate, float offRate,
                                floatingpoint force, int numHeads);
    
    virtual float changeRate(float onRate, float offRate,
                             floatingpoint numHeads, floatingpoint force);
};

///A low duty catch bond implementation of the MotorRateChanger
///
///  p = 0.1, beta = 2.0
///

class LowDutyMotorCatch : public MotorCatch {
    
public:
    LowDutyMotorCatch(short motorType, floatingpoint charForce)
    
    : MotorCatch(motorType, charForce, 0.1, 2.0){}
};

///A high duty catch bond implementation of the MotorRateChanger
///
///  p = 0.33, alpha = 1.0
///

class HighDutyMotorCatch : public MotorCatch {
    
public:
    HighDutyMotorCatch(short motorType, floatingpoint charForce)
    
    : MotorCatch(motorType, charForce, 0.33, 1.0){}
};


///A stall force implementation of the MotorRateChanger.
///Used for a motor walking when under stress.
///Adopted from Hill et al. 1937, and Erdmann et al. 2013.

/// @note - This function updates walking rates based on the Hill form:
///
///   k_eff = k_0 * (F_0 - F) / (F_0 + (F / (alpha)))
///
/// where F_0 is the characteristic force defining this stall,
/// beta is a dimensionless parameter defining the steepness of the curve,
/// k_0 is the walking rate under zero load, which was approximated
/// by Erdmann et al. 2013 to be:
///
///        k_0 = ((N_t - N_b) / N_b) * k_on
///
/// where k_on is the binding rate of a single motor, d_step is
/// the size of a single motor step, d_total is the total step size
/// of the ensemble in simulation, and N_t is the total number of heads.
///
/// It is noted that the true k_0 is also multipilied by a fractional
/// step size corresponding to the step size in simulation,
/// d_step / d_total where d_total is the total step size in simulation,
/// based on the number of binding sites per cylinder.
///
/// N_b is the number of bound motor heads under zero load,
///
///             N_b = p * Nt
///
class MotorStall : public MotorRateChanger  {
    
private:
    floatingpoint _F0;            ///< characteristic force
    float _stepFrac = 1.0; ///< step size of a single head relative to sim
    
    
    //@{
    ///Constant parameters
    floatingpoint _dutyRatio;
    floatingpoint _alpha;
    //@}
    
    
public:
    MotorStall(short motorType, short filamentType, floatingpoint charForce,
               floatingpoint dutyRatio, floatingpoint alpha)
    
    : MotorRateChanger(motorType), _F0(charForce),
    _dutyRatio(dutyRatio), _alpha(alpha) {
        
        //calculate rate based on step fraction
        floatingpoint d_step = SysParams::Chemistry().motorStepSize[_motorType];
        
        floatingpoint d_total = (floatingpoint)SysParams::Geometry().cylinderSize[filamentType] /
        SysParams::Chemistry().numBindingSites[filamentType];
        
        _stepFrac = d_step / d_total;
    }
    
    virtual float numBoundHeads(float onRate, float offRate,
                                floatingpoint force, int numHeads) {return 0.0;}
    
    virtual float changeRate(float onRate, float offRate,
                             floatingpoint numHeads, floatingpoint force);
};

///A low duty stall force implementation of the MotorRateChanger.
///
///         dutyRatio = 0.1, alpha = 0.2
///
class LowDutyMotorStall : public MotorStall {
    
    
public:
    LowDutyMotorStall(short motorType, short filamentType, floatingpoint charForce)
    
    : MotorStall(motorType, filamentType, charForce, 0.1, 0.2){}
};

///A high duty stall force implementation of the MotorRateChanger.
///
///         dutyRatio = 0.33, alpha = 0.3
///
class HighDutyMotorStall : public MotorStall {
    
    
public:
    HighDutyMotorStall(short motorType, short filamentType, floatingpoint charForce)
    
    : MotorStall(motorType, filamentType, charForce, 0.33, 0.3){}
};

} // namespace medyan

#endif
