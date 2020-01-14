//! @file GasSurfReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/Array.h"
#include "cantera/zeroD/GasSurfReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

using namespace std;

namespace Cantera
{

void GasSurfReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "IdealGas") {
        throw CanteraError("GasSurfReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

/* Original Code
void GasSurfReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // Set the third component to the temperature
    y[2] = m_thermo->temperature();

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+3);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 3);
}
*/

void GasSurfReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    size_t mf_start;
    mf_start = 1;
    // Set the third component to the temperature

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+mf_start);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + mf_start);
}


void GasSurfReactor::initialize(doublereal t0)
{
    Reactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

/* Original Code
void GasSurfReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];
    m_thermo->setMassFractions_NoNorm(y+3);
    m_thermo->setState_TR(y[2], m_mass / m_vol);
    updateSurfaceState(y + m_nsp + 3);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}
*/
void GasSurfReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    size_t mf_start;
    double T;
    mf_start = 1;
    T = m_thermo->temperature();
    
    m_thermo->setMassFractions_NoNorm(y+mf_start);
    m_thermo->setState_TR(T, m_mass / m_vol);
    updateSurfaceState(y + m_nsp + mf_start);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}


/* Original code 
void GasSurfReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double mcvdTdt = 0.0; // m * c_v * dT/dt
    double* dYdt = ydot + 3;

    m_thermo->restoreState(m_state);
    applySensitivity(params);
    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    evalWalls(time);
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + 3);
    dmdt += mdot_surf;

    // compression work and external heat transfer
    mcvdTdt += - m_pressure * m_vdot - m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // production in gas phase and from surfaces
        dYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n] / m_mass;
        // dilution by net surface mass flux
        dYdt[n] -= Y[n] * mdot_surf / m_mass;
    }

    // add terms for outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        double mdot_out = m_outlet[i]->massFlowRate(time);
        dmdt -= mdot_out; // mass flow out of system
        mcvdTdt -= mdot_out * m_pressure * m_vol / m_mass; // flow work
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        double mdot_in = m_inlet[i]->massFlowRate(time);
        dmdt += mdot_in; // mass flow into system
        mcvdTdt += m_inlet[i]->enthalpy_mass() * mdot_in;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;

            // In combination with h_in*mdot_in, flow work plus thermal
            // energy carried with the species
            mcvdTdt -= m_uk[n] / mw[n] * mdot_spec;
        }
    }

    ydot[0] = dmdt;
    ydot[1] = m_vdot;
    if (m_energy) {
        ydot[2] = mcvdTdt / (m_mass * m_thermo->cv_mass());
    } else {
        ydot[2] = 0;
    }

    resetSensitivity(params);
}
*/

void GasSurfReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    size_t ystart_no;
    ystart_no = 1;

    double* dYdt = ydot + ystart_no;

    m_thermo->restoreState(m_state);
    applySensitivity(params);
    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    evalWalls(time);
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + ystart_no);
    dmdt += mdot_surf;


    for (size_t n = 0; n < m_nsp; n++) {
        // production in gas phase and from surfaces
        dYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n] / m_mass;
        // dilution by net surface mass flux
        dYdt[n] -= Y[n] * mdot_surf / m_mass;
    }

    // add terms for outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        double mdot_out = m_outlet[i]->massFlowRate(time);
        dmdt -= mdot_out; // mass flow out of system
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        double mdot_in = m_inlet[i]->massFlowRate(time);
        dmdt += mdot_in; // mass flow into system
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;
        }
    }

    ydot[0] = dmdt;
    resetSensitivity(params);
}



void GasSurfReactor::evalJacEqs(doublereal time, doublereal* y, doublereal* ydot,
                                 Array2D* jac)
{
    Array2D &J = *jac;

    const vector_fp& mw = m_thermo->molecularWeights();
    m_work.resize(m_nsp);
    m_kin->updateROPDerivatives();
    auto yloc = 1;
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + yloc);
    for (size_t j = 0; j < m_nsp; j++){             // Eq. (49) of pyjac
        m_kin->getNetProductionRateYDerivatives(m_work.data(), j);
        auto j1 = j + yloc;
        J(j1, 0) = (mdot_surf * y[j1] - m_wdot[j] * m_vol -  m_sdot[j] * mw[j]) / (m_mass * m_mass);
        double wt_frac = m_thermo->meanMolecularWeight() / mw[j];
        for (size_t k = 0; k < m_nsp; k++) {
            auto k1 = k + yloc;
            //J(j,k) = ydot[j] * wt_frac;
            J(k1,j1) = mw[k] / m_thermo->density() * m_work[k];
            if (k1 == j1){
                J(k1, j1) -= mdot_surf/m_mass;
            }
        }
    }
    evalSurfaceDerivatives(time, y, ydot, jac);

    // Outlets 
    for (size_t i = 0; i < m_outlet.size(); i++) {
        J(0,0) -= m_outlet[i]->massFlowRateMassDerivative(true);

        for (size_t j = 0; j < m_nsp; j++) {
            auto j1 = j + yloc;
            J(0, j1) -= m_outlet[i]->massFlowRateYDerivative(j, true);
        }
    }

    // Inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        auto inlet_mass_der = m_inlet[i]->massFlowRateMassDerivative(false);
        J(0,0) += inlet_mass_der; 

        double mdot_in = m_inlet[i]->massFlowRate(time);
        for (size_t j = 0; j < m_nsp; j++) {
            auto j1 = j + yloc;
            J(0, j1) += m_inlet[i]->massFlowRateYDerivative(j, false);

            J(j1,j1) -= mdot_in / m_mass;
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(j);
            double mdot_spec_mass_der = m_inlet[i]->outletSpeciesMassFlowRateMassDerivative(j);
            J(j1, 0) += (inlet_mass_der * y[j1] - mdot_spec_mass_der)/m_mass;
            J(j1, 0) -= (mdot_in * y[j1] - mdot_spec)/(m_mass*m_mass);

        }
    }


}

size_t GasSurfReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 3;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else if (nm == "temperature") {
        return 2;
    } else {
        return npos;
    }
}

std::string GasSurfReactor::componentName(size_t k) {
    if (k == 2) {
        return "temperature";
    } else {
        return Reactor::componentName(k);
    }
}


}
