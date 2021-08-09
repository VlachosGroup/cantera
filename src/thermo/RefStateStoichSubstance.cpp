/**
 *  @file RefStateStoichSubstance.cpp
 * Definition file for the RefStateStoichSubstance class, which represents a fixed-composition
 * incompressible substance and is used in conjuction with surface species, where the effect of
 * pressure on enthalpy is eliminated to conform with the surface species (see \ref thermoprops and
 * class \link Cantera::StoichSubstance StoichSubstance\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/RefStateStoichSubstance.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/ctml.h"

namespace Cantera
{

// ----  Constructors -------

RefStateStoichSubstance::RefStateStoichSubstance(const std::string& infile, const std::string& id_)
{
    initThermoFile(infile, id_);
}

RefStateStoichSubstance::RefStateStoichSubstance(XML_Node& xmlphase, const std::string& id_)
{
    importPhase(xmlphase, this);
}

// ----- Mechanical Equation of State ------


// ---- Chemical Potentials and Activities ----

// Properties of the Standard State of the Species in the Solution

void RefStateStoichSubstance::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
}

// ---- Thermodynamic Values for the Species Reference States ----


// ---- Initialization and Internal functions

/*
void RefStateStoichSubstance::initThermo()
{
    // Make sure there is one and only one species in this phase.
    if (m_kk != 1) {
        throw CanteraError("RefStateStoichSubstance::initThermo",
                           "stoichiometric substances may only contain one species.");
    }

    if (species(0)->input.hasKey("equation-of-state")) {
        auto& eos = species(0)->input["equation-of-state"].getMapWhere(
            "model", "constant-volume");
        if (eos.hasKey("density")) {
            assignDensity(eos.convert("density", "kg/m^3"));
        } else if (eos.hasKey("molar-density")) {
            assignDensity(meanMolecularWeight() *
                            eos.convert("molar-density", "kmol/m^3"));
        } else if (eos.hasKey("molar-volume")) {
            assignDensity(meanMolecularWeight() /
                            eos.convert("molar-volume", "m^3/kmol"));
        } else {
            throw InputFileError("StoichSubstance::initThermo", eos,
                "equation-of-state entry for species '{}' is missing 'density',"
                " 'molar-volume' or 'molar-density' specification",
                speciesName(0));
        }
    } else if (m_input.hasKey("density")) {
        assignDensity(m_input.convert("density", "kg/m^3"));
    }

    // Store the reference pressure in the variables for the class.
    m_p0 = refPressure();

    // Call the base class thermo initializer
    SingleSpeciesTP::initThermo();
}
*/
void RefStateStoichSubstance::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("RefStateStoichSubstance::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    std::string model = tnode["model"];
    if (model != "RefStateStoichSubstance") {
        throw CanteraError("RefStateStoichSubstance::initThermoXML",
                           "thermo model attribute must be RefStateStoichSubstance");
    }
    double dens = getFloat(tnode, "density", "toSI");
    assignDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

void RefStateStoichSubstance::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "RefStateStoichSubstance") {
        throw CanteraError("RefStateStoichSubstance::setParametersFromXML",
                           "thermo model attribute must be RefStateStoichSubstance");
    }
    assignDensity(getFloat(eosdata, "density", "toSI"));
}

}
