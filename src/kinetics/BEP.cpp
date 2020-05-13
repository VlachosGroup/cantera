

#include "cantera/kinetics/BEP.h"
#include "cantera/kinetics/Reaction.h"
//#include "cantera/thermo/LatIntThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include "cantera/base/AnyMap.h"
//#include <iostream>
//#include <limits>
#include <vector>
#include <string>
//#include <utility>


using namespace std;

namespace Cantera {

void BEP::installCleaveReactions(const vector<string>& rxnIds)
{
    installReactions(rxnIds, true);
}

void BEP::installSynthesisReactions(const vector<string>& rxnIds)
{
    installReactions(rxnIds, false);
}

void BEP::installReactions(const vector<string>& rxnIds, bool cleave)
{
    m_reactionIds.insert(m_reactionIds.end(),  rxnIds.begin(), rxnIds.end());
    m_rxnIdIsClvg.insert(m_rxnIdIsClvg.end(), rxnIds.size(), cleave);
}


void BEP::addReaction(string reactionId, bool cleave) 
{
    m_reactionIds.push_back(reactionId);
    m_rxnIdIsClvg.push_back(cleave);
}

void BEP::addReaction(size_t reactionIndex, bool cleave) 
{
    m_reactionIndices.push_back(reactionIndex);
    m_rxnIsClvg.push_back(cleave);
}

bool BEP::reactionInBEP(size_t i) const
{
    auto it = find (m_reactionIndices.begin(), m_reactionIndices.end(), i);
    if (it == m_reactionIndices.end()){
        return false;
    } else {
        return true;
    }
}

bool BEP::reactionInBEP(string id) const
{
    auto it = find (m_reactionIds.begin(), m_reactionIds.end(), id);
    if (it == m_reactionIds.end()){
        return false;
    } else {
        return true;
    }
}

double BEP::actEnergy(double delH, bool rxnIsClvg) const
{
    if (m_BEPIsClvg  == rxnIsClvg) {
        return m_m * delH + m_b;
    } else {
        return (m_m - 1) * delH + m_b;
    }
}


void BEP::computeActivationEnergies(const double* deltaH, double* Ea_R) const
{
    for (size_t i = 0; i < nReactions(); i++){
        if (m_BEPIsClvg == m_rxnIsClvg[i]) {
            Ea_R[i] = m_m * deltaH[m_reactionIndices[i]] + m_b;
        } else {
            Ea_R[i] = (m_m - 1) * deltaH[m_reactionIndices[i]] + m_b;
        }
    }
}

int getBEPReactionArrayIds(const XML_Node& p, vector<string>& ids)
{
    int itot = 0;
    ids.clear();

    // Search the children of the phase element for the XML element named
    // reactionArray. If we can't find it, then return signaling having not
    // found any reactions. Apparently, we allow multiple reactionArray elements
    // here Each one will be processed sequentially, with the end result being
    // purely additive.
    vector<XML_Node*> rarrays = p.getChildren("reactionArray");
    if (rarrays.empty()) {
        return 0;
    }
    for (size_t n = 0; n < rarrays.size(); n++) {
        // Go get a reference to the current XML element, reactionArray. We will
        // process this element now.
        const XML_Node& rxns = *rarrays[n];

        // The reactionArray element has an attribute called, datasrc. The value
        // of the attribute is the XML element comprising the top of the tree of
        // reactions for the phase. Find this datasrc element starting with the
        // root of the current XML node.
        const XML_Node* rdata = get_XML_Node(rxns["datasrc"], &rxns.root());

        // Search for child elements called include. We only include a reaction
        // if it's tagged by one of the include fields. Or, we include all
        // reactions if there are no include fields.
        vector<XML_Node*> incl = rxns.getChildren("include");
        vector<XML_Node*> allrxns = rdata->getChildren("reaction");
        // if no 'include' directive, then include all reactions
        if (incl.empty()) {
            for (size_t i = 0; i < allrxns.size(); i++) {
                ids.push_back(allrxns[i]->attrib("id"));
                ++itot;
            }
        } else {
            for (size_t nii = 0; nii < incl.size(); nii++) {
                const XML_Node& ii = *incl[nii];
                string imin = ii["min"];
                string imax = ii["max"];

                string::size_type iwild = string::npos;
                if (imax == imin) {
                    iwild = imin.find("*");
                    if (iwild != string::npos) {
                        imin = imin.substr(0,iwild);
                        imax = imin;
                    }
                }

                for (size_t i = 0; i < allrxns.size(); i++) {
                    const XML_Node* r = allrxns[i];
                    string rxid;
                    if (r) {
                        rxid = r->attrib("id");
                        if (iwild != string::npos) {
                            rxid = rxid.substr(0,iwild);
                        }

                        // To decide whether the reaction is included or not we
                        // do a lexical min max and operation. This sometimes
                        // has surprising results.
                        if ((rxid >= imin) && (rxid <= imax)) {
                            ids.push_back(rxid);
                            ++itot;
                        }
                    }
                }
            }
        }
    }

    return itot;
}

vector<string> getBEPReactionIds(const AnyMap& bepNode, const AnyMap& rootNode)
{
    // Find sections containing reactions to add
    vector<string> sections, rules;
    vector<string> ids;

    if (bepNode.hasKey("reactions")) {
        const auto& reactionsNode = bepNode.at("reactions");
        if (reactionsNode.is<string>()) {
            if (rootNode.hasKey("reactions")) {
                // Specification of the rule for adding species from the default
                // 'reactions' section, if it exists
                sections.push_back("reactions");
                rules.push_back(reactionsNode.asString());
            } else if (reactionsNode.asString() != "none") {
                throw InputFileError("getBEPReactionIds", reactionsNode,
                    "BEP entry implies existence of 'reactions' section "
                    "which does not exist in the current input file.");
            }
        } else if (reactionsNode.is<vector<string>>()) {
            // List of sections from which all species should be added
            for (const auto& item : reactionsNode.as<vector<string>>()) {
                sections.push_back(item);
                rules.push_back("all");
            }
        } else if (reactionsNode.is<vector<AnyMap>>()) {
            // Mapping of rules to apply for each specified section containing
            // reactions
            for (const auto& item : reactionsNode.as<vector<AnyMap>>()) {
                sections.push_back(item.begin()->first);
                rules.push_back(item.begin()->second.asString());
            }
        }
    } 

    // Add reactions from each section
    for (size_t i = 0; i < sections.size(); i++) {
        
        const auto& slash = boost::ifind_last(sections[i], "/");
        if (slash) {
            // specified section is in a different file
            string fileName (sections[i].begin(), slash.begin());
            string node(slash.end(), sections[i].end());
            AnyMap reactions = AnyMap::fromYamlFile(fileName,
                rootNode.getString("__file__", ""));
            for (const auto& R : reactions[node].asVector<AnyMap>()) {
                ids.push_back(R["id"].asString());
            }
        } else {
            // specified section is in the current file
            for (const auto& R : rootNode.at(sections[i]).asVector<AnyMap>()) {
                ids.push_back(R["id"].asString());
            }
        }
    }

    return ids;
}


shared_ptr<BEP> newBEP(const XML_Node& bep_node)
{
    string id = bep_node["id"];
    double slope = getFloat(bep_node, "alpha");
    double intercept = getFloat(bep_node, "beta", "actEnergy");
    string dir = bep_node.child("direction").value();

    bool isCleave;
    if (dir == "cleavage") {
        isCleave = true;
    } else if (dir == "synthesis") {
        isCleave = false;
    } else {
        throw CanteraError("newBEP", "BEP direction takes only one of "
                "'cleavage', 'synthesis', but {} is given", dir);
    }
    auto bep = make_shared<BEP>(slope, intercept, isCleave, id);

    const XML_Node& clv_rxn_node = bep_node.child("cleavage_reactions");
    vector<string> rxnids;
    auto clv_rxn_no = getBEPReactionArrayIds(clv_rxn_node, rxnids);
    bep->installCleaveReactions(rxnids);

    const XML_Node& syn_rxn_node = bep_node.child("synthesis_reactions");
    auto syn_rxn_no = getBEPReactionArrayIds(syn_rxn_node, rxnids);
    bep->installSynthesisReactions(rxnids);
    return bep;
}


std::vector<shared_ptr<BEP> > getBEPs(const XML_Node& node)
{
    std::vector<shared_ptr<BEP> > beps;
    for (const auto bepnode : node.child("bepData").getChildren("bep")) {
        beps.push_back(newBEP(*bepnode));
    }
    return beps;
}

shared_ptr<BEP> newBEP(const AnyMap& bep_node, const AnyMap& rootNode)
{
    string id = bep_node["id"].asString();
    double slope = bep_node["alpha"].asDouble();
    auto units = bep_node.units();
    double intercept = units.convertActivationEnergy(bep_node["beta"], "K");
    string dir = bep_node["direction"].asString();

    bool isCleave;
    if (dir == "cleavage") {
        isCleave = true;
    } else if (dir == "synthesis") {
        isCleave = false;
    } else {
        throw CanteraError("newBEP", "BEP direction takes only one of "
                "'cleavage', 'synthesis', but {} is given", dir);
    }
    auto bep = make_shared<BEP>(slope, intercept, isCleave, id);

    //auto clv_rxn_node = bep_node["cleavage_reactions"].as<AnyMap>();
    //auto rxnids = getBEPReactionIds(clv_rxn_node, rootNode);
    if (bep_node.hasKey("cleavage_reactions")) {
        auto rxnids = bep_node["cleavage_reactions"].as<vector<string>>();
        bep->installCleaveReactions(rxnids);
    }

    //auto syn_rxn_node = bep_node["synthesis_reactions"].as<AnyMap>();
    //rxnids = getBEPReactionIds(syn_rxn_node, rootNode);
    if (bep_node.hasKey("synthesis_reactions")) {
        auto rxnids = bep_node["synthesis_reactions"].as<vector<string>>();
        bep->installSynthesisReactions(rxnids);
    }
    return bep;
}


}
