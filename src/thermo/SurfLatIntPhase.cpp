/**
 *  @file SurfLatIntPhase.cpp
 *  Definitions for thermodynamic model of a surface phase
 *  derived from ThermoPhase, where lateral interactions of surface species are accounted for
 *  (see \ref thermoprops and class
 *  \link Cantera::SurfLatIntPhase SurfLatIntPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/MultiSpeciesInterThermo.h"
#include "cantera/thermo/LateralInteraction.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
SurfLatIntPhase::SurfLatIntPhase(doublereal n0):
    SurfPhase(n0)
{
}

SurfLatIntPhase::SurfLatIntPhase(const string& infile, const string& id_) :
    SurfPhase(infile, id_)
{
    
}

SurfLatIntPhase::SurfLatIntPhase(XML_Node& xmlphase) :
    SurfPhase(xmlphase)
{
}

bool SurfLatIntPhase::addSpecies(shared_ptr<Species> spec) 
{

    bool added = SurfPhase::addSpecies(spec);
    if (added){
        m_h0_inter.push_back(0.0);
        m_coverages.push_back(0.0);
    }

    return added;
}

//!  Gather a vector of pointers to XML_Nodes for a phase
/*!
 *   @param intrxnDataNodeList   Output vector of pointer to XML_Nodes which 
 *      contain the species XML_Nodes for the interactions in the current phase.
 *   @param intrxnNamesList      Output Vector of strings, which contain 
 *      the names/ids of the interactions in the phase
 *   @param intrxnRuleList       Output Vector of ints, which contain the value 
 *      of intrxnrule for each interaction in the phase
 *   @param intrxnArray_names    Vector of pointers to the XML_Nodes which 
 *      contains the names of the interactions in the phase
 *   @param intrxnArray_dbases   Input vector of pointers to interaction 
 *      databases. We search each data base for the required interaction names
 *   @param  intrxnrule          Input vector of intrxnrule values
 */
static void formInteractionXMLNodeList(vector<XML_Node*> &intrxnDataNodeList,
                                       vector<string> &intrxnNamesList,
                                       vector_int &intrxnRuleList,
                                       const vector<XML_Node*> intrxnArray_names,
                                       const vector<XML_Node*> intrxnArray_dbases,
                                       const vector_int intrxnrule)
{
    // used to check that each interaction is declared only once
    map<string, bool> declared;

    for (size_t jintrxn = 0; jintrxn < intrxnArray_dbases.size(); jintrxn++) {
        const XML_Node& intrxnArray = *intrxnArray_names[jintrxn];

        // Get the top XML for the database
        const XML_Node* db = intrxnArray_dbases[jintrxn];

        // Get the array of intrxn name strings and then count them
        vector<string> intrxnnames;
        getStringArray(intrxnArray, intrxnnames);
        size_t nintrxn = intrxnnames.size();

        // if 'all' is specified as the one and only intrxn in the
        // intrxnArray_names field, then add all intrxn defined in the
        // corresponding database to the phase
        if (nintrxn == 1 && intrxnnames[0] == "all") {
            vector<XML_Node*> allintrxn = db->getChildren("interaction");
            nintrxn = allintrxn.size();
            intrxnnames.resize(nintrxn);
            for (size_t nn = 0; nn < nintrxn; nn++) {
                string stemp = (*allintrxn[nn])["name"];
                if (!declared[stemp] || intrxnrule[jintrxn] < 10) {
                    declared[stemp] = true;
                    intrxnNamesList.push_back(stemp);
                    intrxnDataNodeList.push_back(allintrxn[nn]);
                    intrxnRuleList.push_back(intrxnrule[jintrxn]);
                }
            }
        } else if (nintrxn == 1 && intrxnnames[0] == "unique") {
            vector<XML_Node*> allintrxn = db->getChildren("interaction");
            nintrxn = allintrxn.size();
            intrxnnames.resize(nintrxn);
            for (size_t nn = 0; nn < nintrxn; nn++) {
                string stemp = (*allintrxn[nn])["name"];
                if (!declared[stemp]) {
                    declared[stemp] = true;
                    intrxnNamesList.push_back(stemp);
                    intrxnDataNodeList.push_back(allintrxn[nn]);
                    intrxnRuleList.push_back(intrxnrule[jintrxn]);
                }
            }
        } else {
            map<string, XML_Node*> intrxnNodes;
            for (size_t k = 0; k < db->nChildren(); k++) {
                XML_Node& child = db->child(k);
                intrxnNodes[child["name"]] = &child;
            }
            for (size_t k = 0; k < nintrxn; k++) {
                string stemp = intrxnnames[k];
                if (!declared[stemp] || intrxnrule[jintrxn] < 10) {
                    declared[stemp] = true;
                    // Find the intrxn in the database by name.
                    auto iter = intrxnNodes.find(stemp);
                    if (iter == intrxnNodes.end()) {
                        throw CanteraError("importPhase","no data for interaction, \""
                                           + stemp + "\"");
                    }
                    intrxnNamesList.push_back(stemp);
                    intrxnDataNodeList.push_back(iter->second);
                    intrxnRuleList.push_back(intrxnrule[jintrxn]);
                }
            }
        }
    }
}

bool SurfLatIntPhase::addInteraction(shared_ptr<LateralInteraction> intrxn) {
    if (m_interactions.find(toLowerCopy(intrxn->name())) != m_interactions.end()) {
        throw CanteraError("Phase::addInteraction",
            "Phase '{}' already contains a interaction named '{}'.",
            name(), intrxn->name());
    }

    m_interactions[toLowerCopy(intrxn->name())] = intrxn;
} 

vector<string> SurfLatIntPhase::getAffectedInteractions(string speciesName) 
{
    vector<string> intrxnNames;
    intrxnNames.reserve(nSpecies());
    for (auto const& intrxn: m_interactions) 
        if (intrxn.second->species2Name() == speciesName)
            intrxnNames.push_back(intrxn.first);

    return intrxnNames;
}

vector<string> SurfLatIntPhase::getAffectingInteractions(string speciesName) 
{
    vector<string> intrxnNames;
    intrxnNames.reserve(nSpecies());
    for (auto const& intrxn: m_interactions) 
        if (intrxn.second->species1Name() == speciesName)
            intrxnNames.push_back(intrxn.first);

    return intrxnNames;
}

shared_ptr<LateralInteraction> SurfLatIntPhase::getInteractionfromID(string id) 
{
    return m_interactions[toLowerCopy(id)];
}

bool SurfLatIntPhase::installInteractionArrays(const XML_Node& p, 
                                               bool check_for_duplicates)
{
    int itot = 0;

    vector<XML_Node*> iarrays = p.getChildren("interactionArray");
    if (iarrays.empty()) {
        return false;
    }
    for (size_t n = 0; n < iarrays.size(); n++) {
        // Go get a reference to the current XML element, interactionArray. We will
        // process this element now.
        const XML_Node& intrxns = *iarrays[n];

        // The interactionArray element has an attribute called, datasrc. The value
        // of the attribute is the XML element comprising the top of the tree of
        // reactions for the phase. Find this datasrc element starting with the
        // root of the current XML node.
        const XML_Node* idata = get_XML_Node(intrxns["datasrc"], &intrxns.root());

        // We will set intrxnrule.skipUndeclaredSpecies to 'true'. intrxnrule is
        // passed to the routine that parses each individual interaction so that
        // the parser will skip all interactions containing a species not part
        // of the current phase without throwing an error.
        //kin.skipUndeclaredSpecies(true);

        // Search for child elements called include. We only include a reaction
        // if it's tagged by one of the include fields. Or, we include all
        // reactions if there are no include fields.
        vector<XML_Node*> incl = intrxns.getChildren("include");
        vector<XML_Node*> allintrxns = idata->getChildren("interaction");
        // if no 'include' directive, then include all reactions
        if (incl.empty()) {
            for (size_t i = 0; i < allintrxns.size(); i++) {
                string intrxnid = allintrxns[i]->attrib("id");
                addInteraction(newLateralInteraction(*allintrxns[i]));
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

                for (size_t i = 0; i < allintrxns.size(); i++) {
                    const XML_Node* intrxn = allintrxns[i];
                    string intrxnid;
                    if (intrxn) {
                        intrxnid = intrxn->attrib("id");
                        if (iwild != string::npos) {
                            intrxnid = intrxnid.substr(0,iwild);
                        }

                        // To decide whether the interaction is included or not we
                        // do a lexical min max and operation. This sometimes
                        // has surprising results.
                        if ((intrxnid >= imin) && (intrxnid <= imax)) {
                            addInteraction(newLateralInteraction(*intrxn));
                            ++itot;
                        }
                    }
                }
            }
        }
    }

    /* Add the interactions to MultiSpeciesInterThermo */
    for (const auto & intrx_id_pair : m_interactions) {
        m_spInterThermo.addInteraction(intrx_id_pair.second);
    }
    m_spInterThermo.buildSpeciesInterMap(speciesNames()); 


    if (check_for_duplicates) { 
        checkDuplicates(); //Not yet implemented
    }

    return true;

}

void SurfLatIntPhase::initThermoXML(XML_Node& phaseNode, const string& id)
{
    ThermoPhase::initThermoXML(phaseNode, id);

    // Read the interaction data here. This could be put in importPhase, but
    // interactions are not generic for all phases. So as a stopgap measure
    // reading the interaction data as remnant data for SurfLatIntPhase.
    // Lot of code is duplicated from species importing part in the importPhase 
    // function. Could be simplified.
    installInteractionArrays(phaseNode, true);
    return;

    if (phaseNode.hasChild("interactionArray")){
        vector<XML_Node*> intrxnArrays = phaseNode.getChildren("interactionArray");

        vector<XML_Node*> dbases;
        vector_int intrxnrule(intrxnArrays.size(), 0);

        // Default behavior when importing interactions containing unknown species
        // is to discard them

        for (size_t i = 0; i < intrxnArrays.size(); i++) {
            const XML_Node& intrxnArray = *intrxnArrays[i];

            XML_Node* db = get_XML_Node(intrxnArray["datasrc"], &phaseNode.root());
            if (db == 0) {
                throw CanteraError("importPhase() Can not find XML node for ", 
                                   "interaction database: "
                                   + intrxnArray["datasrc"]);
            }
            dbases.push_back(db);
        }

        // Now collect all the interxns and the XML_Node * pointer for those
        // interxns in a single vector. This is where we decide what interxns
        // are to be included in the phase. 
        vector<XML_Node*> intrxnDataNodeList;
        vector<string> intrxnNamesList;
        vector_int intrxnRuleList;
        formInteractionXMLNodeList(intrxnDataNodeList, intrxnNamesList, 
                                   intrxnRuleList, intrxnArrays, 
                                   dbases, intrxnrule);

        size_t nIntrxn = intrxnDataNodeList.size();
        for (size_t k = 0; k < nIntrxn; k++) {
            XML_Node* s = intrxnDataNodeList[k];
            AssertTrace(s != 0);
            //if (intrxnRuleList[k]) {
            //   ignoreUndefinedSpecies();
            //}
            this->addInteraction(newLateralInteraction(*s));
        }
    }

}


doublereal SurfLatIntPhase::enthalpy_mole() const
{
    if (m_n0 <= 0.0) {
        return 0.0;
    }
    _updateThermo();
    return mean_X(m_h0);
}

/*
void SurfLatIntPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }
}
*/

doublereal SurfLatIntPhase::entropy_mole() const
{
    _updateThermo();
    doublereal s = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        s += moleFraction(k) * (m_s0[k] -
            GasConstant * log(std::max(concentration(k) * size(k)/m_n0, SmallNumber)));
    }
    return s;
}

doublereal SurfLatIntPhase::cp_mole() const
{
    _updateThermo();
    return mean_X(m_cp0);
}

void SurfLatIntPhase::getStandardChemPotentials(doublereal* mu0) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu0);
}

void SurfLatIntPhase::getChemPotentials(doublereal* mu) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu);
    getActivityConcentrations(m_work.data());
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += RT() * (log(m_work[k]) - logStandardConc(k));
    }
}

void SurfLatIntPhase::getPureGibbs(doublereal* g) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), g);
}

void SurfLatIntPhase::getGibbs_RT(doublereal* grt) const
{
    _updateThermo();
    scale(m_mu0.begin(), m_mu0.end(), grt, 1.0/RT());
}

void SurfLatIntPhase::getEnthalpy_RT(doublereal* hrt) const
{
    _updateThermo();
    scale(m_h0.begin(), m_h0.end(), hrt, 1.0/RT());
}

void SurfLatIntPhase::getEntropy_R(doublereal* sr) const
{
    _updateThermo();
    scale(m_s0.begin(), m_s0.end(), sr, 1.0/GasConstant);
}

void SurfLatIntPhase::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    scale(m_cp0.begin(), m_cp0.end(), cpr, 1.0/GasConstant);
}

void SurfLatIntPhase::getStandardVolumes(doublereal* vol) const
{
    _updateThermo();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = 1.0/standardConcentration(k);
    }
}



void SurfLatIntPhase::_updateThermo(bool force) const
{
    doublereal tnow = temperature();
    //if (m_tlast != tnow || force) { //TODO: This check is needed 
        m_spthermo.update(tnow, m_cp0.data(), m_h0.data(), m_s0.data());

        getCoverages(m_coverages.data());
        m_spInterThermo.update(tnow, m_coverages.data(), m_h0_inter.data());

        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_h0[k] += m_h0_inter[k];
            m_h0[k] *= GasConstant * tnow;
            m_s0[k] *= GasConstant;
            m_cp0[k] *= GasConstant;
            m_mu0[k] = m_h0[k] - tnow*m_s0[k];
        }
        m_tlast = tnow;
    //}
}

/*
void SurfLatIntPhase::_updateInterThermo() 
{
    doublereal tnow = temperature();
    m_spInterThermo.update(tnow, m_h0.data());
}
*/

void SurfLatIntPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","SurfaceCoverage");
    doublereal n = getFloat(eosdata, "site_density", "toSI");
    setSiteDensity(n);
}

}
