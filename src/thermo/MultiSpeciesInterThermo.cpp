/**
 *  @file MultiSpeciesInterThermo.cpp
 *  Declarations for a thermodynamic property manager for multiple species
 *  in a phase (see \ref spthermo and
 * \link Cantera::MultiSpeciesInterThermo MultiSpeciesInterThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MultiSpeciesInterThermo.h"
#include "cantera/thermo/LateralInteraction.h"
//#include "cantera/thermo/MultiSpeciesThermo.h"
//#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/eigen_dense.h"
#include <utility>
#include <iostream>

using namespace std;

namespace Cantera
{

MultiSpeciesInterThermo::MultiSpeciesInterThermo() :
    m_nSpecies(0)
{
}

MultiSpeciesInterThermo::MultiSpeciesInterThermo(
        vector<shared_ptr<LateralInteraction> > interactions) :
    m_nSpecies(0),
    m_interactions(interactions)
{
}

void MultiSpeciesInterThermo::addInteraction(
        shared_ptr<LateralInteraction> interaction)
{
    m_interactions.push_back(interaction); 
}

int get_index(vector<string> species, string name) {
    auto it = find(species.begin(), species.end(), name);
    if (it == species.end())
    {
        // name not in vector
        throw CanteraError("MultiSpeciesInterThermo::get_index",
                "Name of species in interaction not found in supplied species");
    }
    return distance(species.begin(), it);
}

void MultiSpeciesInterThermo::buildSpeciesInterMap(std::vector<std::string> species)
{
    auto cnt = 0;
    for (auto inter : m_interactions) {
        string sp1 = inter->species1Name();
        string sp2 = inter->species2Name();

        auto ind1 = get_index(species, sp1);
        auto ind2 = get_index(species, sp2);
        m_species_intrxn_map[make_pair(ind1, ind2)] = cnt;
        cnt++;
    }

    if (m_int_strengths.size() == 0) { // Initialize to zero values. 
        m_int_strengths = Eigen::MatrixXd::Zero(species.size(), species.size());
    }

    m_nSpecies = species.size();

    //if (m_coverages.size() == 0) { // Not yet initialized, initialize to zero values
    //    m_coverages = Eigen::VectorXd::Zero(species.size());
    //}

}

void MultiSpeciesInterThermo::update(doublereal t, double* coverages, 
                                     double* h_RT) const
{
    for (int i = 0; i <  m_nSpecies; i++){
        //if (coverages[i] != m_coverages[i]){
            //m_coverages(i) = coverages[i];
            //for (auto j=0; j < m_coverages.size(); j++){
                //auto inter = m_interactions[m_specie_inter_map[std::make_pair(j, i)]];
                //m_int_strengths(j,i) = inter->strength(m_coverages[i]);
            //}
        //}
        for (int j=0; j < m_nSpecies; j++){
            auto it = m_species_intrxn_map.find(make_pair(j,i));
            if (it != m_species_intrxn_map.end()) {
                auto ind = it->second;
                m_int_strengths(j,i) = m_interactions[ind]->strength(coverages[i]);
            }
        }
    }

    Eigen::VectorXd inter_h0_RT = m_int_strengths.rowwise().sum()/(GasConstant * t);
    for (size_t i = 0; i < m_nSpecies; i++) {
        h_RT[i] = inter_h0_RT[i];
    }
}

/*
void MultiSpeciesInterThermo::reportInteractionParams(size_t index, int& type,
        doublereal* const c, doublereal& covIntercept1, doublereal& covIntercept2,
        doublereal& refPressure_) const
{
    const SpeciesThermoInterpType* sp = provideSTIT(index);
    size_t n;
    if (sp) {
        sp->reportParameters(n, type, minTemp_, maxTemp_,
                             refPressure_, c);
    } else {
        type = -1;
    }
}
*/



bool MultiSpeciesInterThermo::ready(size_t nSpecies) {
    if (m_interactions.size() < nSpecies * nSpecies) {
        return false;
    }
    return true;
}


}
