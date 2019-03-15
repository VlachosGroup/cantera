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
    
    LateralInteraction(std::string species1, std::string species2,
                       const vector_fp& interaction_strengths,
                       const vector_fp& coverage_intercepts, 
                       std::string name);
                       
    ~LateralInteraction();

    bool validate();

    std::string species1Name();     //Affected Species

    std::string species2Name();    //Affecting Species

    std::string name() { return m_id; }

    //! Get the interaction strength at the given coverage
    double strength(const double coverage) const;


protected:
    //std::pair<std::shared_ptr<Species>, std::shared_ptr<Species> > m_species; 
    std::pair<std::string, std::string> m_species; 
    vector_fp m_strengths;
    vector_fp m_cov_thresholds;
    std::string m_id;

};

shared_ptr<LateralInteraction> newLateralInteraction(const XML_Node& interaction_node);


std::vector<shared_ptr<LateralInteraction> > getInteractions(const XML_Node& node);


}
#endif
