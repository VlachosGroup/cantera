//! @file GasOnlyReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/Array.h"
#include "cantera/zeroD/GasOnlyReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

using namespace std;

namespace Cantera
{

void GasOnlyReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "IdealGas") {
        throw CanteraError("GasOnlyReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void GasOnlyReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y);

}


void GasOnlyReactor::initialize(doublereal t0)
{
    Reactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

void GasOnlyReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    double T = m_thermo->temperature();
    m_thermo->setMassFractions_NoNorm(y);
    //m_thermo->setState_TR(T, m_mass / m_vol);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

void GasOnlyReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    double* dYdt = ydot;

    m_thermo->restoreState(m_state);
    applySensitivity(params);
    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }


    for (size_t n = 0; n < m_nsp; n++) {
        // production in gas phase 
        dYdt[n] = (m_wdot[n] * m_vol) * mw[n] / m_mass;
    }


    resetSensitivity(params);
}



void GasOnlyReactor::evalJacEqs(doublereal time, doublereal* y, doublereal* ydot,
                                 Array2D* jac)
{
    Array2D &J = *jac;

    //cout << "Before getCp_R" << endl;
    //m_thermo->getCp_R(m_work.data());                   // C_p/R
    //cout << "After getCp_R" << endl;
    //m_thermo->getPartialMolarIntEnergies(m_uk.data());
    //cout << "After getPartialMolarIntEnergies" << endl;
    //m_kin->getNetProductionRates(m_wdot.data());              
    //cout << "Before getNetProductionRateTDerivatives" << endl;
    //vector_fp dwdotdT(m_nsp);
    //m_kin->getNetProductionRateTDerivatives(dwdotdT.data());
    //cout << "After getNetProductionRateTDerivatives" << endl;
    //double mdot_surf = evalSurfaces(time, ydot + m_nsp + 3); 

    /*
    double df1dT {0}, df1dT_2t{0};
    for (size_t i = 0; i < m_nsp; i++) { 
        auto CvR = m_work[i] - 1;
        CvR -= m_uk[i] * dcvRdT / m_thermo->cv_mass();  // C_v(k)/R - u(k) * dc_v/R/dT
        CvR -= m_uk[i]/RT;          //TODO: Get U/RT as well
        CvR *= (m_wdot[i] * m_vol + m_sdot[i]);         // End of first term Eq (46)
        df1dT += CvR;
    }
    df1dT *= inv_mcv * GasConstant; // First term of Eq. (46) of pyjac
    cout << "After df1dT" << endl;

    for (size_t i = 0; i < m_nsp; i++) { 
        auto prod_rate = m_vol * dwdotdT[i]; //TODO: + dsdotdT 
        df1dT_2t +=  m_uk[i] *  prod_rate;
    }
    df1dT_2t *= inv_mcv;
    cout << "After df1dT_2t" << endl;
    
    */

    
    // J(j, k) = $\dho Y_j / dho Y_k$
    const vector_fp& mw = m_thermo->molecularWeights();
    m_work.resize(m_nsp);
    m_kin->updateROPDerivatives();
    for (size_t j = 0; j < m_nsp; j++){             // Eq. (49) of pyjac
        m_kin->getNetProductionRateYDerivatives(m_work.data(), j);
        double wt_frac = m_thermo->meanMolecularWeight() / mw[j];
        for (size_t k = 0; k < m_nsp; k++) {
            //J(j,k) = ydot[j] * wt_frac;
            J(k,j) = mw[k] / m_thermo->density() * m_work[k];
        }
    }
}

size_t GasOnlyReactor::componentIndex(const string& nm) const
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

std::string GasOnlyReactor::componentName(size_t k) {
    if (k == 2) {
        return "temperature";
    } else {
        return Reactor::componentName(k);
    }
}


}
