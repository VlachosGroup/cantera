/**
 *  @file SurfLatIntPhase.cpp
 *  Definitions for thermodynamic model of a surface phase
 *  derived from ThermoPhase, where lateral interactions of surface species are accounted for
 *  (see \ref thermoprops and class
 *  \link Cantera::SurfIntPhase SurfIntPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/MultiSpeciesInterThermo.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
SurfLatIntPhase::SurfLatIntPhase(doublereal n0):
    SurfPhase(n0)
{
}

SurfLatIntPhase::SurfLatIntPhase(const std::string& infile, const std::string& id_) :
    SurfPhase(infile, id_)
{
    
}

SurfLatIntPhase::SurfLatIntPhase(XML_Node& xmlphase) :
    SurfPhase(xmlphase)
{
}

/*
doublereal SurfLatIntPhase::enthalpy_mole() const
{
    if (m_n0 <= 0.0) {
        return 0.0;
    }
    _updateThermo();
    return mean_X(m_h0);
}

void SurfLatIntPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }
}

void SurfLatIntPhase::_updateThermo(bool force) const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow || force) {
        m_spthermo.update(tnow, m_cp0.data(), m_h0.data(), m_s0.data());
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_h0[k] *= GasConstant * tnow;
            m_s0[k] *= GasConstant;
            m_cp0[k] *= GasConstant;
            m_mu0[k] = m_h0[k] - tnow*m_s0[k];
        }
        _updateInterThermo();
        m_tlast = tnow;
    }
}

void SurfLatIntPhase::_updateThermo(bool force, ) const
{
    doublereal tnow = temperature();
    m_spInterThermo.update(tnow, m_h0.data());
}

void SurfLatIntPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","SurfaceCoverage");
    doublereal n = getFloat(eosdata, "site_density", "toSI");
    setSiteDensity(n);
}
*/
}
