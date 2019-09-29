/**
 *  @file BEP.h
 *  Header for defining activation energy in terms of reaction enthalpies 
 *  through Brønsted–Evans–Polanyi (BEP) principle.
 *  \link Cantera::BEP BEP\endlink).
 */

// This file is part of Cantera/openmkm. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef HTRCT_BEP_H
#define HTRCT_BEP_H

#include <set>
#include <string>
#include <memory>

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class XML_Node;
class Reaction;

//! Definition of a BEP for a class of reactions.  
/*!
 * Activation energy of a reaction is defined through a BEP relation when 
 * direct data is not avaialable. This is a popular approach in heterogenous
 * catalysis. BEP relation formulates activation energy as a linear function 
 * of Gibbs reaction energy. The slope and the intercept of the linear function
 * are generally identical for a class of reactions.
 * In this model, the activation energy \f$ E_a \f$ is given as
 *
 *           \f[
 *               E_a = m \Delta H_{rxn} + b,
 *           \f]
 * 
 * where \f$ m \f$ and \f$ b \f$ represent the slope and intercep of the BEP
 * relation.
 *
 */
class BEP
{
public:
    BEP() : m_m(0), m_b(0), m_BEPIsClvg(true), m_id(0) {}

    BEP(double slope, double intercept, bool cleave, std::string id) : 
        m_m(slope), m_b(intercept), m_BEPIsClvg(cleave), m_id(id) {}
                       
    ~BEP() {}

    //! Add cleavage reactions 
    void installCleaveReactions(const std::vector<std::string>& rxnIds); 

    //! Add cleavage reactions 
    void installSynthesisReactions(const std::vector<std::string>& rxnIds); 

    //! Add reaction ids in bulk
    void installReactions(const std::vector<std::string>& rxnIds, bool cleave); 

    //! Check that the specified coverage bounds are valid 
    void addReaction(std::string reactionId, bool cleave);

    //! Check that the specified coverage bounds are valid 
    void addReaction(size_t reactionIndex, bool cleave);

    //! Indices of the reactions in the kinetics object
    const std::vector<size_t>& getReactionIndices() const 
    {
        return m_reactionIndices; 
    }

    //! Ids of the reactions
    const std::vector<std::string>& getReactionIds() const 
    {
        return m_reactionIds; 
    }

    size_t nReactions() const 
    {
        return m_reactionIndices.size();
    }

    size_t reactionIndexSize() const 
    {
        return m_reactionIndices.size();
    }

    size_t reactionIdSize() const 
    {
        return m_reactionIds.size();
    }

    bool rxnIdIsCleave(size_t i) const 
    {
        return m_rxnIdIsClvg[i];
    }

    bool rxnIsCleave(size_t i) const 
    {
        return m_rxnIsClvg[i];
    }

    bool reactionInBEP(size_t i) const; 

    bool reactionInBEP(std::string id) const;

    double actEnergy(double delH, bool rxnIsClvg) const;

    void computeActivationEnergies(const double* deltaH_R, double* Ea_R) const;

    //! Get the slope of the BEP
    double slope() const { return m_m; }

    //! Get the intercept of the BEP
    double intercept() const { return m_b; }

    //! Get the id of the BEP
    std::string id() const { return m_id; }

protected:
    double m_m;
    double m_b;
    std::string m_id;
    bool m_BEPIsClvg;
    std::vector<size_t> m_reactionIndices = {};     // To identify the reactions in kinetics objects associated with this BEP
    std::vector<std::string> m_reactionIds = {}; // Used to read the reactions from xml file associated with this BEP
    std::vector<bool> m_rxnIsClvg = {};
    std::vector<bool> m_rxnIdIsClvg = {};

};

//! Create a new BEP object for a given BEP relation 
//! defined in `bep_node`
shared_ptr<BEP> newBEP(const XML_Node& bep_node);

//! Create BEP objects for all `<bep>` nodes in an 
//! XML document.
//!
//! The `<bep>` nodes are assumed to be children of the 
//! `<bepData>` node in an XML document with a `<ctml>` root node, 
//! as in the case of XML files produced by conversion from CTI files.
//!
//! This function can be used in combination with get_XML_File() and
//! get_XML_from_string() to get BEP objects from either a file or a
//! string, respectively, where the string or file is formatted as either CTI
//! or XML.
//!
std::vector<shared_ptr<BEP> > getBEPs(const XML_Node& node);

}
#endif
