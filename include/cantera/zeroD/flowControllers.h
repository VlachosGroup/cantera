//! @file flowControllers.h Some flow devices derived from class FlowDevice.

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWCONTR_H
#define CT_FLOWCONTR_H

#include <iostream>
#include "FlowDevice.h"
#include "ReactorBase.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{
/**
 * A class for mass flow controllers. The mass flow rate is constant or
 * specified as a function of time..
 */
class MassFlowController : public FlowDevice
{
public:
    MassFlowController() : FlowDevice() {
        m_type = MFC_Type;
    }

    virtual bool ready() {
        return FlowDevice::ready() && m_mdot >= 0.0;
    }

    /// If a function of time has been specified for mdot, then update the
    /// stored mass flow rate. Otherwise, mdot is a constant, and does not
    /// need updating.
    virtual void updateMassFlowRate(doublereal time) {
        if (m_func) {
            m_mdot = m_func->eval(time);
        }
        m_mdot = std::max(m_mdot, 0.0);
    }
};

/**
 * A class for flow controllers where the flow rate is equal to the flow rate
 * of a "master" mass flow controller plus a correction proportional to the
 * pressure difference between the inlet and outlet.
 */
class PressureController : public FlowDevice
{
public:
    PressureController() : FlowDevice(), m_master(0) {
        m_type = PressureController_Type;
    }

    virtual bool ready() {
        return FlowDevice::ready() && m_master != 0 && m_coeffs.size() == 1;
    }

    void setMaster(FlowDevice* master) {
        m_master = master;
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = \dot{m}_{master} + c \Delta P \f]
     */
    void setPressureCoeff(double c) {
        m_coeffs = {c};
    }

    virtual void updateMassFlowRate(doublereal time) {
        if (!ready()) {
            throw CanteraError("PressureController::updateMassFlowRate",
                "Device is not ready; some parameters have not been set.");
        }
        m_mdot = m_master->massFlowRate(time)
                 + m_coeffs[0]*(in().pressure() - out().pressure());
        m_mdot = std::max(m_mdot, 0.0);
    }

    virtual double massFlowRateMassDerivative(bool upstream){ // Call it only for non-const Pressure reactors
        double vol, RT, meanW;
        if (upstream){
            vol = in().volume();
            RT = in().contents().RT();
            meanW = in().contents().meanMolecularWeight();
        } else {
            vol = -out().volume();
            RT = out().contents().RT();
            meanW = out().contents().meanMolecularWeight();
        }
        
        auto derivative = m_coeffs[0] * RT / (meanW * vol);
        return derivative;
    }

    virtual double massFlowRateYDerivative(size_t k, bool upstream){ // Call it only for non-const Pressure reactors
        double den, RT, Wk;
        if (upstream){
            den = in().density();
            RT = in().contents().RT();
            Wk = in().contents().molecularWeight(k);
        } else {
            den = -out().density();
            RT = out().contents().RT();
            Wk = out().contents().molecularWeight(k);
        }
        
        auto derivative = m_coeffs[0] * RT * den / Wk;
        return derivative;
    }

    virtual double massFlowRateTDerivative(bool upstream){ 
        double den, meanW;
        if (upstream){
            den = in().density();
            meanW = in().contents().meanMolecularWeight();
        } else {
            den = -out().density();
            meanW = out().contents().meanMolecularWeight();
        }
        auto derivative = m_coeffs[0] * GasConstant * den / meanW;
        return derivative;
    }

protected:
    FlowDevice* m_master;
};

//! Supply a mass flow rate that is a function of the pressure drop across the
//! valve.
/*!
 * The default behavior is a linearly proportional to the pressure difference.
 * Note that real valves do not have this behavior, so this class does not
 * model real, physical valves.
 */
class Valve : public FlowDevice
{
public:
    Valve() : FlowDevice() {
        m_type = Valve_Type;
    }

    virtual bool ready() {
        return FlowDevice::ready() && (m_coeffs.size() == 1 || m_func);
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = c \Delta P \f]
     */
    void setPressureCoeff(double c) {
        m_coeffs = {c};
    }

    /// Compute the currrent mass flow rate, based on the pressure difference.
    virtual void updateMassFlowRate(doublereal time) {
        if (!ready()) {
            throw CanteraError("Valve::updateMassFlowRate",
                "Device is not ready; some parameters have not been set.");
        }
        double delta_P = in().pressure() - out().pressure();
        if (m_func) {
            m_mdot = m_func->eval(delta_P);
        } else {
            m_mdot = m_coeffs[0]*delta_P;
        }
        m_mdot = std::max(m_mdot, 0.0);
    }

    virtual double massFlowRateMassDerivative(bool upstream){
        double der_P;
        double vol, RT, meanW;
        if (upstream){
            vol = in().volume();
            RT = in().contents().RT();
            meanW = in().contents().meanMolecularWeight();
            der_P = RT / (vol * meanW); 
        } else {
            vol = out().volume();
            RT = out().contents().RT();
            meanW = out().contents().meanMolecularWeight();
            der_P = -RT / (vol * meanW); 
        }

        if (m_func){
            double delta_P = in().pressure() - out().pressure();
            return m_func->derivative()(delta_P) * der_P;
        }
        return der_P * m_coeffs[0];
    }

    virtual double massFlowRateTDerivative(bool upstream){
        double der_P;
        double den, meanW;
        if (upstream){
            den = in().density();
            meanW = in().contents().meanMolecularWeight();
            der_P = GasConstant * den /  meanW; 
        } else {
            den = out().density();
            meanW = out().contents().meanMolecularWeight();
            der_P = -GasConstant * den / meanW; 
        }

        if (m_func){
            double delta_P = in().pressure() - out().pressure();
            return m_func->derivative()(delta_P) * der_P;
        }
        return der_P * m_coeffs[0];
    }

    virtual double massFlowRateYDerivative(size_t k, bool upstream){
        double der_P;
        double den, RT, Wk;
        if (upstream){
            den = in().density();
            RT = in().contents().RT();
            Wk = in().contents().molecularWeight(k);
            der_P = RT * den / Wk; 
        } else {
            den = out().density();
            RT = out().contents().RT();
            Wk = out().contents().molecularWeight(k);
            der_P = -RT * den / Wk; 
        }

        if (m_func){
            double delta_P = in().pressure() - out().pressure();
            return m_func->derivative()(delta_P) * der_P;
        } 
        return der_P * m_coeffs[0];
    }
};

}
#endif
