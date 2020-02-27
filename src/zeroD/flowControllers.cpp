//! @file flowControllers.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/flowControllers.h"
#include "cantera/zeroD/ReactorBase.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

MassFlowController::MassFlowController() : FlowDevice() {
    m_type = MFC_Type;
}

void MassFlowController::setMassFlowRate(double mdot)
{
    if (m_tfunc) {
        delete m_tfunc;
    }
    m_coeff = mdot;
}

void MassFlowController::updateMassFlowRate(double time)
{
    if (!ready()) {
        throw CanteraError("MassFlowController::updateMassFlowRate",
                           "Device is not ready; some parameters have not been set.");
    }
    double mdot = m_coeff;
    if (m_tfunc) {
        mdot *= m_tfunc->eval(time);
    }
    m_mdot = std::max(mdot, 0.0);
}

PressureController::PressureController() : FlowDevice(), m_master(0) {
    m_type = PressureController_Type;
}

void PressureController::updateMassFlowRate(double time)
{
    if (!ready()) {
        throw CanteraError("PressureController::updateMassFlowRate",
                           "Device is not ready; some parameters have not been set.");
    }
    double mdot = m_coeff;
    double delta_P = in().pressure() - out().pressure();
    if (m_pfunc) {
        mdot *= m_pfunc->eval(delta_P);
    } else {
        mdot *= delta_P;
    }
    mdot += m_master->massFlowRate(time);
    m_mdot = std::max(mdot, 0.0);
}

double PressureController::massFlowRateMassDerivative(bool upstream)
{ 
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
    
    auto derivative = m_master->massFlowRateMassDerivative(false) + m_coeff * RT / (meanW * vol);
    return derivative;
}

double PressureController::massFlowRateYDerivative(size_t k, bool upstream)
{ 
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
    
    auto derivative = m_coeff * RT * den / Wk;
    return derivative;
}

double PressureController::massFlowRateTDerivative(bool upstream)
{
    double den, meanW;
    if (upstream){
        den = in().density();
        meanW = in().contents().meanMolecularWeight();
    } else {
        den = -out().density();
        meanW = out().contents().meanMolecularWeight();
    }
    auto derivative = m_coeff * GasConstant * den / meanW;
    return derivative;
}

Valve::Valve() : FlowDevice() {
    m_type = Valve_Type;
}

void Valve::updateMassFlowRate(double time)
{
    if (!ready()) {
        throw CanteraError("Valve::updateMassFlowRate",
                           "Device is not ready; some parameters have not been set.");
    }
    double mdot = m_coeff;
    if (m_tfunc) {
        mdot *= m_tfunc->eval(time);
    }
    double delta_P = in().pressure() - out().pressure();
    if (m_pfunc) {
        mdot *= m_pfunc->eval(delta_P);
    } else {
        mdot *= delta_P;
    }
    m_mdot = std::max(mdot, 0.0);
}

double Valve::massFlowRateMassDerivative(bool upstream)
{
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

    if (m_pfunc){
        double delta_P = in().pressure() - out().pressure();
        return m_pfunc->derivative()(delta_P) * der_P;
    }
    return der_P * m_coeff;
}

double Valve::massFlowRateTDerivative(bool upstream)
{
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

    if (m_pfunc){
        double delta_P = in().pressure() - out().pressure();
        return m_pfunc->derivative()(delta_P) * der_P;
    }
    return der_P * m_coeff;
}

double Valve::massFlowRateYDerivative(size_t k, bool upstream)
{
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

    if (m_pfunc){
        double delta_P = in().pressure() - out().pressure();
        return m_pfunc->derivative()(delta_P) * der_P;
    } 
    return der_P * m_coeff;
}

}
