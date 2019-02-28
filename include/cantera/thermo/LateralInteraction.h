/**
 *  @file LateralInteraction.h
 */

// This file is part of hetero_ct.

#ifndef HTRCT_LATINT_H
#define HTRCT_LATINT_H

#include <set>
#include <memory>

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class XML_Node;
class Species;

//! Class which stores data about the lateral interactions between the 
//! surface species.
class LateralInteraction
{
public:
    LateralInteraction();
    /*
    LateralInteraction(const Species* const species, 
                       const double* const interaction_strengths,
                       const double* const coverage_intercepts);
                       */
    ~LateralInteraction();

    bool validate();

    //! Get the interaction strength at the given coverage
    double strength(double coverage);


protected:
    std::set<std::shared_ptr<Species> > m_species; 
    vector_fp m_strengths;
    vector_fp m_cov_thresholds;
    std::string m_id;

};

std::vector<shared_ptr<LateralInteraction> > getInteractions(const XML_Node& node);


}
#endif
