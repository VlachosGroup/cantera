/**
 * @file MultiSpeciesInterThermo.h
 *  Header for a general species thermodynamic property manager for a phase where 
 *  species have interactions (see
 * \link Cantera::MultiSpeciesInterThermo MultiSpeciesInterThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTISPECIESINTERTHERMO_H
#define CT_MULTISPECIESINTERTHERMO_H

//#include "MultiSpeciesThermo.h"
#include <memory>
#include <vector>
#include <utility>
#include <map>
#include "cantera/numerics/eigen_dense.h"

namespace Cantera
{

class LateralInteraction;

//! A species thermodynamic property manager for a phase that explicitly accounts 
//! for interactions between species.
/*!
 * This is a derived class of MultiSpeciesThermo, which is a general manager that 
 * can handle a wide variety of species * thermodynamic polynomials for individual 
 * species and compute their * nondimensional, reference-state thermodynamic 
 * properties (i.e. as a function * of temperature only). Whereas MultiSpeciesThermo
 * does not account for interaction between species, MultiSpeciesInteractionThermo
 * does based on the interaction strength defined. 
 *
 * Developed to account for lateral interactions between surface species
 *
 * The ThermoPhase object relies on MultiSpeciesInterThermo to calculate the
 * contributions of interactions to the thermodynamic properties of the 
 * reference state for all of the species in the
 * phase, for a range of temperatures. 
 *
 * The most important member function for the MultiSpeciesInterThermo class is the
 * member function MultiSpeciesInterThermo::update(). The function calculates the
 * values of Cp/R, H/RT, and S/R coming from lateral interactions for all of the 
 * surface species at once at the specified coverage.
 *
 * Usually, all of the species in a interacting surface phase are installed into a
 * MultiSpeciesInterThermo object. However, there is no requirement that a
 * MultiSpeciesInterThermo object handles all of the species in a phase. The member
 * function
 * \link MultiSpeciesInterThermo::install_STIT() install_STIT()\endlink
 * is called to install each species into the MultiSpeciesSurfThermo object.
 *
 * @ingroup spthermo
 */
class MultiSpeciesInterThermo 
{
public:
    //! Constructor
    MultiSpeciesInterThermo() {}

    MultiSpeciesInterThermo(
            std::vector<std::shared_ptr<LateralInteraction> > interactions);

    // MultiSpeciesInterThermo objects are not copyable or assignable
    MultiSpeciesInterThermo(const MultiSpeciesInterThermo& b) = delete;
    MultiSpeciesInterThermo& operator=(const MultiSpeciesInterThermo& b) = delete;
    virtual ~MultiSpeciesInterThermo() {}

    //! Install interaction parameterization of species
    /*!
     *  @param index           Index of the species being installed
     *  @param interParameters interaction parameters of the species
     */
    void addInteraction(std::shared_ptr<LateralInteraction> interaction);

    //! Build the map between the indices of species and interactions
    /*!
     *  @param species   Ordered species strings in a phase 
     */
    void buildSpeciesInterMap(std::vector<std::string> species);

    //! Compute the coverage dependent lateral interactions for all species.
    /*!
     * Given temperature T in K, this method updates the values of -
     * enthalpy for the given coverages, Pref of each of the standard states.
     *
     * @param T         Temperature (Kelvin)
     * @param coverages Coverages (vector of floats)
     * @param h_RT      Vector of Dimensionless enthalpies. (length m_kk).
     */
    void update(doublereal T, double* coverages, 
                double* h_rt) const;

    int nSpecies() {return m_nSpecies;}

    //! This utility function reports back the interaction
    //! parameters for the species with index number *index*.
    /*!
     * @param index         Species index
     * @param type          Integer type of the Interaction type.
     * @param c             Vector of interaction coefficients.
     * @param covIntercept1 output -  1st intercept of coverage. One means doesn't exist.
     * @param covIntercept2 output -  2nd intercept of coverage. One means doesn't exist.
     */
    /*void reportInteractionParams(size_t index, int& type,
                                 doublereal* const c,
                                 doublereal& covIntercept1,
                                 doublereal& covIntercept2) const;
    */
    //! Check if data for all species (0 through nSpecies-1) has been installed.
    bool ready(size_t nSpecies);

protected:

    //! Store the interaction objects to obtain the interactions strengths
    std::vector<std::shared_ptr<LateralInteraction> > m_interactions;
    //! Map the interaction index to species indices
    std::map< std::pair<int, int>, int > m_species_intrxn_map;
    //! Interaction strength parameters  stored in matrix form
    mutable Eigen::MatrixXd m_int_strengths;

    int m_nSpecies;
    //! Vector of coverages 
    //Eigen::VectorXd m_coverages;

};

}

#endif
