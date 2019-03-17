/**
 *  @file LateralInteraction.h
 *  Header for a piecewise linear model of lateral interactions of surface
 *  species (see \ref thermoprops and calss 
 *  \link Cantera::LateralInteraction LateralInteraction\endlink).
 */

// This file is part of Cantera/hetero_ct. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef HTRCT_LATINT_H
#define HTRCT_LATINT_H

#include <set>
#include <memory>

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class XML_Node;
class Species;

//! A generic piecewise linear model for coverage dependent lateral 
//! interactions between surface species.  
/*!
 * Coverage dependent lateral interactions of surface species are accounted 
 * for by their contribution to the formation enthalpies of the surface species. 
 * In this model, the contribution of lateral interactions to formation 
 * enthalpy of a surface species \f$ i \f$ is given as
 *
 *           \f[
 *               \Delta H_i = \sum_j{LI_{i,j}\theta_j},
 *           \f]
 * 
 * where \f$ LI_{i,j} \f$ represents the lateral interaction parameter between 
 * species \f$ i \f$ and \f$ j \f$, and \f$ \theta_j \f$ represent the coverage 
 * of species \f$ j \f$. \f$ LI_{i,j} \f$ need not be symmetric w.r.t. \f$ i \f$
 * and \f$ j \f$.  Further, the interaction strength parameters are given as 
 * piece wise linear functions of surface coverage. The piecewise linear 
 * functions are defined by their slopes and the coverage boundaries. 0 and 1
 * are the default lower and upper coverage boundaries. A two piecewise linear function, 
 * \f$ L_{i,j} \f$, given as
 *
 *           \f[
 *              LI'_{i,j} &= \alpha_{i,j}\;\;\; \forall\; \theta_j < \theta_0 \\
 *                        &= \beta_{i,j}\;\;\; \forall\; \theta_j > \theta_0,
 *           \f]
 *
 * is specified by three parameters \f$ \alpha \f$, \f$ \beta \f$, and 
 * \f$ \theta_0 \f$. This model is also called 3 parameter model. Similarly, one 
 * can specify a one-parameter model which is essentially a linear function in 
 * between the coverages 0 and 1, and 5-parameter model with 3 piecewise linear 
 * functions in between 0 and 1 coverage limits. The implmentation allows for
 * arbitrary number of piecewise linear functions. The lateral interaction is
 * then given by integrating the function slope values till the given coverage.
 *
 */
class LateralInteraction
{
public:
    LateralInteraction();
    
    LateralInteraction(std::string species1, std::string species2,
                       const vector_fp& interaction_strengths,
                       const vector_fp& coverage_intercepts, 
                       std::string name);
                       
    ~LateralInteraction();

    //! Check that the specified coverage bounds are valid 
    bool validate();

    //! The species affected by the interaction
    std::string species1Name();     

    //! The species influencing the interaction
    std::string species2Name();    

    //! Id of the interaction
    std::string name() { return m_id; }

    //! Get the interaction strength at the given coverage
    double strength(const double coverage) const;

protected:
    std::pair<std::string, std::string> m_species; 
    vector_fp m_strengths;
    vector_fp m_cov_thresholds;
    std::string m_id;

};

//! Create a new LateralInteraction object for the lateral interaction 
//! defined in `intrxn_node`
shared_ptr<LateralInteraction> newLateralInteraction(const XML_Node& intrxn_node);

//! Create LateralInteraction objects for all `<interaction>` nodes in an 
//! XML document.
//!
//! The `<interaction>` nodes are assumed to be children of the 
//! `<interactionData>` node in an XML document with a `<ctml>` root node, 
//! as in the case of XML files produced by conversion from CTI files.
//!
//! This function can be used in combination with get_XML_File() and
//! get_XML_from_string() to get Interaction objects from either a file or a
//! string, respectively, where the string or file is formatted as either CTI
//! or XML.
//!
//! If Interaction objects are being created from a CTI definition that does not
//! contain corresponding phase definitions, then one of the following must be
//! true, or the resulting interaction parameters will be incorrect:
//!
//!   - The interaction strength constants are expressed in (Joules/mol) units
//!   - A `units` directive is included 
std::vector<shared_ptr<LateralInteraction> > getInteractions(const XML_Node& node);


}
#endif
