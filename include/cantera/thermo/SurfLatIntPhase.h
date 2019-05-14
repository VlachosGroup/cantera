/**
 *  @file SurfLatIntPhase.h
 *  Header for a thermodynamics model of a surface phase where laterial interactions
 *  between different species types are accounted for
 *  derived from SurfPhase,
 *  (see \ref thermoprops and class \link Cantera::SurLatIntfPhase SurfLatIntPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_SURFLATINTPHASE_H
#define CT_SURFLATINTPHASE_H

#include <memory>
#include "SurfPhase.h"
#include "MultiSpeciesInterThermo.h"

namespace Cantera
{


//class MultiSpeciesInterThermo;
//! A thermodynamic model for a surface phase, where lateral interactions between
//! different species types is accounted for 
//! model.
/*!
 * The surface consists of a grid of equivalent sites. Surface species may be
 * defined to occupy one or more sites. The surface species are assumed to be
 * interacting and the interaction strength depends on their coverage. At 
 * dilute limit, ideal solution limit should be recovered.
 *
 * The density of surface sites is given by the variable \f$ n_0 \f$,
 * which has SI units of kmol m-2.
 *
 * ## Specification of Species Standard State Properties
 *
 * It is assumed that the reference state thermodynamics may be obtained by a
 * pointer to a populated species thermodynamic property manager class (see
 * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
 * state thermodynamics is resolved at this level. How to relate coverage 
 * changes to the reference state thermodynamics is resolved at this level.
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * Therefore, The standard state internal energy for species *k* is equal to the
 * enthalpy for species *k*.
 *
 *       \f[
 *            u^o_k = h^o_k
 *       \f]
 *
 * Also, the standard state chemical potentials, entropy, and heat capacities
 * are independent of pressure. The standard state Gibbs free energy is obtained
 * from the enthalpy and entropy functions.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The activity of species defined in the phase is given by
 *       \f[
 *            a_k = \theta_k
 *       \f]
 *
 * The chemical potential for species *k* is equal to
 *       \f[
 *            \mu_k(T,P) = \mu^o_k(T) + R T \log(\theta_k)
 *       \f]
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * The internal energy for species k is equal to the enthalpy for species *k*
 *       \f[
 *            u_k = h_k
 *       \f]
 *
 * The entropy for the phase is given by the following relation, which is
 * independent of the pressure:
 *
 *       \f[
 *            s_k(T,P) = s^o_k(T) - R \log(\theta_k)
 *       \f]
 *
 * ## Specification of lateral interactions based on coverages
 * 
 * Depending on the coverage, the enthalpy for species *k* is modified by
 *       \f[
 *            h_k -= \frac{\sum_j{max(\theta_j - \theta_0, 0)*f_{j,k}}}{R T},
 *       \f]
 * where \f$\theta_0\f$ is the threshold coverage and \f$f_{j,k}\f$ is the interaction
 * strength between the species j and k.
 *
 * ## %Application within Kinetics Managers
 *
 * The activity concentration,\f$  C^a_k \f$, used by the kinetics manager, is equal to
 * the actual concentration, \f$ C^s_k \f$, and is given by the following
 * expression.
 *       \f[
 *            C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}
 *       \f]
 *
 * The standard concentration for species *k* is:
 *        \f[
 *            C^0_k = \frac{n_0}{s_k}
 *        \f]
 *
 * ## Instantiation of the Class
 *
 * The constructor for this phase is located in the default ThermoFactory
 * for %Cantera. A new SurfLatIntPhase may be created by the following code snippet:
 *
 * @code
 *    XML_Node *xc = get_XML_File("diamond.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "diamond_100");
 *    ThermoPhase *diamond100TP_tp = newPhase(*xs);
 *    SurfLatIntPhase *diamond100TP = dynamic_cast <SurfLatIntPhase *>(diamond100TP_tp);
 * @endcode
 *
 * or by the following constructor:
 *
 * @code
 *    XML_Node *xc = get_XML_File("diamond.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "diamond_100");
 *    SurfLatIntPhase *diamond100TP = new SurfLatIntPhase(*xs);
 * @endcode
 *
 * ## XML Example
 *
 * An example of an XML Element named phase setting up a SurfLatIntPhase object named
 * diamond_100 is given below.
 *
 * @code
 * <phase dim="2" id="diamond_100">
 *    <elementArray datasrc="elements.xml">H C</elementArray>
 *    <speciesArray datasrc="#species_data">c6HH c6H* c6*H c6** c6HM c6HM* c6*M c6B </speciesArray>
 *    <reactionArray datasrc="#reaction_data"/>
 *    <state>
 *       <temperature units="K">1200.0</temperature>
 *       <coverages>c6H*:0.1, c6HH:0.9</coverages>
 *    </state>
 *    <thermo model="SurfaceCoverage">
 *       <site_density units="mol/cm2">3e-09</site_density>
 *       <define lateral interaction parameters here>
 *    </thermo>
 *    <kinetics model="Interface"/>
 *    <transport model="None"/>
 *    <phaseArray>
 *         gas_phase diamond_bulk
 *    </phaseArray>
 * </phase>
 * @endcode
 *
 * The model attribute, "SurfaceCoverage", on the thermo element identifies the phase as being
 * a SurfLatIntPhase object.
 *
 * @ingroup thermoprops
 */
class SurfLatIntPhase : public SurfPhase
{
public:
    //! Constructor.
    /*!
     *  @param n0 Site Density of the Surface Phase
     *            Units: kmol m-2.
     */
    SurfLatIntPhase(doublereal n0 = 1.0);

    //! Construct and initialize a SurfLatIntPhase ThermoPhase object directly from an
    //! ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    SurfLatIntPhase(const std::string& infile, const std::string& id);

    //! Construct and initialize a SurfLatIntPhase ThermoPhase object directly from an
    //! XML database
    /*!
     *  @param xmlphase XML node pointing to a SurfLatIntPhase description
     */
    SurfLatIntPhase(XML_Node& xmlphase);

    virtual std::string type() const {
        return "Surf";
    }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     * \f[
     * \hat h(T,P) = \sum_k X_k \hat h^0_k(T),
     * \f]
     * and is a function only of temperature. The standard-state pure-species
     * Enthalpies \f$ \hat h^0_k(T) \f$ are computed by the species
     * thermodynamic property manager.
     *
     * \see MultiSpeciesThermo
     */
    //virtual doublereal enthalpy_mole() const;

    //! Set the Equation-of-State parameters by reading an XML Node Input
    /*!
     * The Equation-of-State data consists of one item, the site density.
     *
     * @param thermoData   Reference to an XML_Node named thermo containing the
     *                     equation-of-state data. The XML_Node is within the
     *                     phase XML_Node describing the SurfPhase object.
     *
     * An example of the contents of the thermoData XML_Node is provided below.
     * The units attribute is used to supply the units of the site density in
     * any convenient form. Internally it is changed into MKS form.
     *
     * @code
     *    <thermo model="SurfaceCoverage">
     *       <site_density units="mol/cm2"> 3e-09 </site_density>
     *    </thermo>
     * @endcode
     */
    virtual void setParametersFromXML(const XML_Node& thermoData);

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    bool addInteraction(std::shared_ptr<LateralInteraction> intrxn);
    
    virtual bool addSpecies(shared_ptr<Species> spec);

    bool installInteractionArrays(const XML_Node& p, bool check_for_duplicates);

    void checkDuplicates() {}

    std::vector<std::string> getAffectedInteractions(std::string speciesName);

    std::vector<std::string> getAffectingInteractions(std::string speciesName);

    std::shared_ptr<LateralInteraction> getInteractionfromID(std::string id);

    int nInteractions() const { return m_interactions.size(); }

protected:
    //! Vector of lateral interaction parameters (number of species squared). length m_kk**2.
    MultiSpeciesInterThermo m_spInterThermo;

    std::map<std::string, std::shared_ptr<LateralInteraction> > m_interactions;

    mutable vector_fp m_coverages;  // Working array for the interactions

    mutable vector_fp m_h0_inter;   // Temporarily store the enthalpy delta from interactions


//private:
    //! Update the species reference state thermodynamic functions
    /*!
     * The polynomials for the standard state functions are only reevaluated if
     * the temperature has changed.
     *
     * @param force   Boolean, which if true, forces a reevaluation of the thermo
     *                polynomials. default = false.
     * @param lat_int Boolean, which if true, modifies the enthalpy to account for
     *                lateral interactions of surface species. default = false.
     */
    virtual void _updateThermo(bool force=false) const;
};
}

#endif
