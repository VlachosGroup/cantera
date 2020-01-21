//! @file IdealGasReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/Array.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

using namespace std;

namespace Cantera
{

void IdealGasReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

/* Original Code
void IdealGasReactor::getState(double* y)
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

void IdealGasReactor::getState(double* y)
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
    if (m_energy){
        mf_start = 2;
    } else {
        mf_start = 1;
    }
    // Set the third component to the temperature
    if (m_energy){
        y[1] = m_thermo->temperature();
    }

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+mf_start);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + mf_start);
}


void IdealGasReactor::initialize(doublereal t0)
{
    Reactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

/* Original Code
void IdealGasReactor::updateState(doublereal* y)
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
void IdealGasReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    size_t mf_start;
    double T;
    if (m_energy) {
        mf_start = 2;
        T = y[1];
    } else {
        mf_start = 1;
        T = m_thermo->temperature();
    }
    
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
void IdealGasReactor::evalEqs(doublereal time, doublereal* y,
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

void IdealGasReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double mcvdTdt = 0.0; // m * c_v * dT/dt
    size_t ystart_no;
    if (m_energy)
        ystart_no = 2;
    else
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

    // compression work and external heat transfer
    //mcvdTdt += - m_pressure * m_vdot - m_Q;
    mcvdTdt -= m_Q;  // Const volume no moving walls

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
    //ydot[1] = m_vdot; Volume term removed
    if (m_energy) {
        ydot[1] = mcvdTdt / (m_mass * m_thermo->cv_mass());
    } 

    resetSensitivity(params);
}



void IdealGasReactor::evalJacEqs(doublereal time, doublereal* y, doublereal* ydot,
                                 Array2D* jac)
{
    Array2D &J = *jac;
    
    const size_t m_ind = 0;
    //const size_t V_ind = 1;
    const size_t T_ind = 1;
    size_t y_ind;
    if (m_energy)
        y_ind = 2;
    else
        y_ind = 1;

    /* Volume related terms. 
     * For Idea Gas Constant volume Reactor, these are trivial w.r.t. reactor */
    /*
    J(V_ind, m_ind) = 0;
    J(V_ind, V_ind) = 1;
    J(V_ind, T_ind) = 0;
    for (size_t i = 0; i < m_nsp; i++) J(V_ind, y_ind + i) =0;
    */

    // Temperature Derivatives
    const vector_fp& mw = m_thermo->molecularWeights(); 
    auto inv_mcv = 1.0/(m_mass * m_thermo->cv_mass());

    if (m_energy){
        vector_fp dwdotdT(m_nsp);
        m_kin->getNetProductionRateTDerivatives(dwdotdT.data());
        for (size_t j = 0; j < m_nsp; j++){
            auto j1 = j + y_ind;
            J(j1, T_ind) = mw[j] * dwdotdT[j] / m_thermo->density();
        }

        double dcvRdT = 0;                                  // Small c/R
        vector<double> dCvRdT(m_nsp);                       // Capital C/R
        //double RT = GasConstant * m_thermo->temperature();
        double RT = m_thermo->RT();
        m_thermo->getdCp_RdT(dCvRdT.data());                // dC_p/R/dT=dC_v/R/dT
        for (size_t i = 0; i < m_nsp; i++) {                // Compute dcvdT
            dcvRdT += dCvRdT[i] / mw[i] * y[i];
        }
        m_work.resize(m_nsp);
        m_thermo->getCp_R(m_work.data());                   // C_p/R

        m_thermo->getPartialMolarIntEnergies(m_uk.data());  // U
        m_kin->getNetProductionRates(m_wdot.data());              

        double df1dT {0}, df1dT_2t{0};
        for (size_t j = 0; j < m_nsp; j++) { 
            auto CvR = m_work[j] - 1;
            CvR -= m_uk[j] * dcvRdT / m_thermo->cv_mass();  // C_v(k)/R - u(k) * dc_v/R/dT
            CvR -= m_uk[j]/RT;          
            CvR *= (m_wdot[j] * m_vol + m_sdot[j]);         // End of first term Eq (46)
            df1dT += CvR;
        }
        df1dT *= inv_mcv * GasConstant;                     // First term of Eq. (46) of pyjac

        for (size_t i = 0; i < m_nsp; i++) { 
            auto prod_rate = m_vol * dwdotdT[i]; //TODO: + dsdotdT 
            df1dT_2t +=  m_uk[i] *  prod_rate;
        }
        df1dT_2t *= inv_mcv;

        J(T_ind, T_ind) = df1dT - df1dT_2t;

        double dfTdm_1t{0}; 
        for (size_t j = 0; j < m_nsp; j++) { 
            dfTdm_1t  += m_uk[j] * (m_wdot[j] * m_vol + m_sdot[j]);  
        }
        J(T_ind, m_ind) = dfTdm_1t * inv_mcv / m_mass; 
    }

    m_work.resize(m_nsp);
    m_kin->updateROPDerivatives();
    vector_fp Cv(m_nsp);
    m_thermo->getCp_R(Cv.data());                   // C_p/R
    double dfTdg_t1 = 0;
    if (m_energy){
        for (size_t j = 0; j < m_nsp; j++) {
            dfTdg_t1 += (m_vol * m_wdot[j] -  m_sdot[j]) * m_uk[j];
        }
    }

    double mdot_surf = evalSurfaces(time, ydot + m_nsp + y_ind);
    for (size_t j = 0; j < m_nsp; j++){             // Eq. (49) of pyjac
        m_kin->getNetProductionRateYDerivatives(m_work.data(), j);
        auto j1 = j + y_ind;
        J(j1, m_ind) = (mdot_surf * y[j1] - (m_wdot[j] * m_vol -  m_sdot[j]) * mw[j]) / (m_mass * m_mass);
        double wt_frac = m_thermo->meanMolecularWeight() / mw[j];
        for (size_t k = 0; k < m_nsp; k++) {
            auto k1 = k + y_ind;
            J(k1,j1) = mw[k] / m_thermo->density() * m_work[k];
            if (k1 == j1){
                J(k1, j1) -= mdot_surf/m_mass;
            }
        }

        if (m_energy){
            double dfTdg_t2 = 0;
            for (size_t k = 0; k < m_nsp; k++) {
                dfTdg_t2 += m_work[k] * m_uk[k];
            }

            J(T_ind, j1) = inv_mcv * (
                    dfTdg_t1 * (Cv[j]-1) * GasConstant / (mw[j] * m_thermo->cv_mass()) -
                    dfTdg_t2 * m_vol);
        }
    }
    evalSurfaceDerivatives(time, y, ydot, jac);

    // Outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        J(m_ind, m_ind) -= m_outlet[i]->massFlowRateMassDerivative(true);

        for (size_t j = 0; j < m_nsp; j++) {
            auto j1 = j + y_ind;
            J(m_ind, j1) -= m_outlet[i]->massFlowRateYDerivative(j, true);
        }
    }

    // Inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        auto inlet_mass_der = m_inlet[i]->massFlowRateMassDerivative(false);
        J(m_ind,m_ind) += inlet_mass_der;

        double mdot_in = m_inlet[i]->massFlowRate(time);
        for (size_t j = 0; j < m_nsp; j++) {
            auto j1 = j + y_ind;
            J(m_ind, j1) += m_inlet[i]->massFlowRateYDerivative(j, false);

            J(j1, j1) -= mdot_in / m_mass;
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(j);
            double mdot_spec_mass_der = m_inlet[i]->outletSpeciesMassFlowRateMassDerivative(j);
            J(j1, m_ind) += (inlet_mass_der * y[j1] - mdot_spec_mass_der)/m_mass;
            J(j1, m_ind) -= (mdot_in * y[j1] - mdot_spec)/(m_mass*m_mass);
        }
    }


    /*
    // Mass related terms 
    J(m_ind, m_ind) = 1;
    //J(m_ind, V_ind) = 0;



    // Mass fraction derivatives 
    // TODO: Implement A\sum_j W_j ds_j/dY_k 
    for (size_t i = 0; i < m_nsp; i++) J(m_ind, y_ind + i) =0;

     Temperature related terms 
    // J(T_ind, T_ind) = $\dho Tdot / dho T$
    cout << "After dcvrdt" << endl;
    //vector_fp dsdotdT(m_nsp);

    cout << "After getCp_R" << endl;
    cout << "After getPartialMolarIntEnergies" << endl;
    cout << "Before getNetProductionRateTDerivatives" << endl;
    cout << "After getNetProductionRateTDerivatives" << endl;
    //double mdot_surf = evalSurfaces(time, ydot + m_nsp + 3); 

    cout << "After df1dT" << endl;

        cout << J(T_ind, T_ind) << endl;

    // J(k, T_ind) = $\dho Y_k/ dho T$
    double  T = m_thermo->temperature();
    for (size_t i = 0; i < m_nsp; i++) {
        J(y_ind + i, T_ind) = mw[i] /m_mass * ( 
                (m_vol*m_wdot[i]+m_sdot[i])/T  + dwdotdT[i]); //TODO: + dsdotdT
    }
    cout << "After yind" << endl;

    // Mass fraction derivatives
    // J(T_ind, k) = $\dho Tdot / dho Y_k$

    
    // J(j, k) = $\dho Y_j / dho Y_k$
    */
    

}

size_t IdealGasReactor::componentIndex(const string& nm) const
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

std::string IdealGasReactor::componentName(size_t k) {
    if (k == 2) {
        return "temperature";
    } else {
        return Reactor::componentName(k);
    }
}


}
