

#include "cantera/thermo/LateralInteraction.h"
#include "cantera/thermo/Species.h"
//#include "cantera/thermo/LatIntThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include <iostream>
#include <limits>

namespace Cantera {

LateralInteraction::LateralInteraction()
{
}

/*
LateralInteraction::LateralInteraction(const Species* const species, 
                                       const double* const strengths,   
                                       const double* const intercepts)
    : m_spcies(species), m_strengths(strengths), m_intercepts(intercepts)
{
}
*/

LateralInteraction::~LateralInteraction()
{
}

bool LateralInteraction::validate()
{
    if (m_strengths.size() == m_cov_thresholds.size()+1)
        return true;
    else
        return false;
}


std::string LateralInteraction::species1Name() { return m_species.first->name; }

std::string LateralInteraction::species2Name() { return m_species.second->name; }


shared_ptr<LateralInteraction> newLateralInteraction(const XML_Node& interaction_node)
{
    if (interaction_node.hasAttrib("id")){
        std::string id = interaction_node["id"];
    }
    const XML_Node& sp_array = interaction_node.child("speciesArray");
    XML_Node* db = get_XML_Node(sp_array["datasrc"], &interaction_node.root());
    /*if (db == 0) {
        throw CanteraError(" Can not find the XML node for species databases: ");
                           //+ sp_array["datasrc"]);
    }*/
    
    std::vector<XML_Node*> fas = interaction_node.getChildren("floatArray");
    vector_fp m_strengths, m_cov_thresholds;
    for (auto fa: fas){
        if (fa->name() == "strength")
            getFloatArray(*fa, m_strengths);//, nodeName = "strength")
        if (fa->name() == "coverage_threshold")
            getFloatArray(*fa, m_cov_thresholds);//, nodeName = "coverage_threshold")
    }

    
    //auto interaction = make_shared<LateralInteraction>(species, strengths, intercepts);
    auto interaction = make_shared<LateralInteraction>();

    return interaction;
}


std::vector<shared_ptr<LateralInteraction> > getInteractions(const XML_Node& node)
{
    std::vector<shared_ptr<LateralInteraction> > interactions;
    for (const auto& intnode : node.child("interactionData").getChildren("interaction")) {
        interactions.push_back(newLateralInteraction(*intnode));
    }
    return interactions;
}

}
