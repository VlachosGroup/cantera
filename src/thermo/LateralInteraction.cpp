

#include "cantera/thermo/LateralInteraction.h"
#include "cantera/thermo/Species.h"
//#include "cantera/thermo/LatIntThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <utility>


using namespace std;

namespace Cantera {

LateralInteraction::LateralInteraction()
{
}


LateralInteraction::LateralInteraction(std::string species1, std::string species2,
                                       const vector_fp& strengths,   
                                       const vector_fp& intercepts, 
                                       std::string name) :
    m_strengths(strengths), m_cov_thresholds(intercepts), m_id(name)
{
    m_species  = make_pair(species1, species2);
}


LateralInteraction::~LateralInteraction()
{
}

bool LateralInteraction::validate()
{
    return (m_strengths.size() + 1 == m_cov_thresholds.size()) ? true : false;
}


std::string LateralInteraction::species1Name() const { 
    //return m_species.first->name; 
    return m_species.first; 
}

std::string LateralInteraction::species2Name() const { 
    //return m_species.second->name; 
    return m_species.second; 
}

double LateralInteraction::strength(const double coverage) const 
{
    /*
    if (coverage > 1) {
        throw CanteraError ("Cantera::LateralInteraction", 
                "Coverage '{}' greater than 1", coverage);
    }
    else if (coverage < 0){
        throw CanteraError ("Cantera::LateralInteraction", 
                "Coverage '{}' less than 0", coverage);
    }
    */
    double cov = (coverage < 0) ? 0.0 : coverage;
    cov = (cov > 1) ? 1.0 : coverage;
    doublereal val = 0.0;
    for (size_t i = 0; i < m_strengths.size(); i++) {
        auto cov_low_thr = m_cov_thresholds[i];
        auto cov_up_thr = m_cov_thresholds[i+1];
        if (cov_up_thr < cov) {
            val += (cov_up_thr - cov_low_thr) * m_strengths[i];
        } 
        else {
            val += (cov- cov_low_thr) * m_strengths[i];
            break;
        }
    }

    return  val;
}



shared_ptr<LateralInteraction> newLateralInteraction(const XML_Node& interaction_node)
{
    std::string id = interaction_node["id"];
    const XML_Node& sp_array = interaction_node.child("speciesArray");
    std::vector<std::string> species;
    getStringArray(sp_array, species);
    if (species.size() != 2)
        throw CanteraError("Cantera::newLateralInteraction", 
                "The size of the species array: '{}' is different from 2",  
                species.size());
                           //+ sp_array["datasrc"]);
    //XML_Node* db = get_XML_Node(sp_array["datasrc"], &interaction_node.root());
    /*if (db == 0) {
        throw CanteraError(" Can not find the XML node for species databases: ");
                           //+ sp_array["datasrc"]);
    }*/
    
    std::vector<XML_Node*> fas = interaction_node.getChildren("floatArray");
    vector_fp strengths, cov_thresholds;
    for (auto fa: fas){
        if (fa->attrib("name") == "strength")
            getFloatArray(*fa, strengths, true, "actEnergy");
        if (fa->attrib("name") == "coverage_threshold")
            getFloatArray(*fa, cov_thresholds, false);
    }

    auto interaction = make_shared<LateralInteraction>(species[0], species[1], 
                                                       strengths, cov_thresholds, id);
    return interaction;
}

std::vector<shared_ptr<LateralInteraction>> getInteractions(const XML_Node& node)
{
    std::vector<shared_ptr<LateralInteraction>> interactions;
    for (const auto& intnode : node.child("interactionData").getChildren("interaction")) {
        interactions.push_back(newLateralInteraction(*intnode));
    }
    return interactions;
}

shared_ptr<LateralInteraction> newLateralInteraction(const AnyMap& intrxnNode)
{
    
    std::string id = intrxnNode["id"].asString();
    auto species = intrxnNode["species"].asVector<string>();
    if (species.size() != 2)
        throw CanteraError("Cantera::newLateralInteraction", 
                "The size of the species array: '{}' is different from 2",  
                species.size());
    auto units = intrxnNode.units();
    auto strengths = intrxnNode["strength"].asVector<AnyValue>(); // Update to account for units
    //auto strengths = intrxnNode.convertVector("strength", "K", 
    vector<double> slopes;
    for (auto& strength : strengths) {
        slopes.push_back(units.convertActivationEnergy(strength, "J/kmol"));
    }
    auto covs = intrxnNode["coverage-threshold"].asVector<double>();
    

    auto interaction = make_shared<LateralInteraction>(species[0], species[1], 
                                                       slopes, covs, id);
    return interaction;
}

std::vector<shared_ptr<LateralInteraction>> getInteractions(const AnyMap& intrxnNode, const AnyMap& rootNode)
{
    std::vector<shared_ptr<LateralInteraction>> interactions;
    /*for (const auto& intnode : node.child("interactionData").getChildren("interaction")) {
        interactions.push_back(newLateralInteraction(*intnode));
    }*/
    if (intrxnNode["interactions"].is<vector<AnyMap>>()) {
        for (const auto& intrxnNd : intrxnNode["interactions"].asVector<AnyMap>()){
            interactions.push_back(newLateralInteraction(intrxnNd));
        }
    }
    return interactions;
}


}
