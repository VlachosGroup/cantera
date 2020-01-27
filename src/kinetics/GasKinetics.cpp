/**
 *  @file GasKinetics.cpp Homogeneous kinetics in ideal gases
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/GasKinetics.h"

using namespace std;

namespace Cantera
{
GasKinetics::GasKinetics(thermo_t* thermo) :
    BulkKinetics(thermo),
    m_logp_ref(0.0),
    m_logc_ref(0.0),
    m_logStandConc(0.0),
    m_pres(0.0)
{
}

void GasKinetics::update_rates_T()
{
    doublereal T = thermo().temperature();
    doublereal P = thermo().pressure();
    m_logStandConc = log(thermo().standardConcentration());
    doublereal logT = log(T);

    if (T != m_temp) {
        if (!m_rfn.empty()) {
            m_rates.update(T, logT, m_rfn.data());
        }

        if (!m_rfn_low.empty()) {
            m_falloff_low_rates.update(T, logT, m_rfn_low.data());
            m_falloff_high_rates.update(T, logT, m_rfn_high.data());
        }
        if (!falloff_work.empty()) {
            m_falloffn.updateTemp(T, falloff_work.data());
        }
        updateKc();
        m_ROP_ok = false;
    }

    if (T != m_temp || P != m_pres) {
        if (m_plog_rates.nReactions()) {
            m_plog_rates.update(T, logT, m_rfn.data());
            m_ROP_ok = false;
        }

        if (m_cheb_rates.nReactions()) {
            m_cheb_rates.update(T, logT, m_rfn.data());
            m_ROP_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}

void GasKinetics::updateTDerivativeFactors()
{
    doublereal T = thermo().temperature();
    doublereal logT = log(T);

    if (!m_rfn_dTMult.empty()) {
        m_rates.update_TDerivative(T, logT, m_rfn_dTMult.data());
    }

    if (!m_rfn_low.empty()) {
        m_falloff_low_rates.update_TDerivative(T, logT, m_rfn_low_dT.data());
        m_falloff_high_rates.update_TDerivative(T, logT, m_rfn_high_dT.data());

        // Update the m_rfn_dTMult for falloff reactions
        for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
            if (reactionType(m_fallindx[i]) == FALLOFF_RXN) {
                m_rfn_dTMult[m_fallindx[i]] = m_rfn_high_dT[i];
            } else { // CHEMACT_RXN
                m_rfn_dTMult[m_fallindx[i]] = m_rfn_high_dT[i];
            }
        }
    }

    if (!falloff_work_dT.empty()) {
        m_falloffn.updateTempDerivative(T, falloff_work_dT.data());
    }

    vector_fp dbdt(m_kk);
    thermo().getdBdT(dbdt.data());
    for (size_t i = 0; i < m_kk; i++){
        dbdt[i] -= 1/T;
    }
    getRevReactionDelta(dbdt.data(),  m_dBdT.data());
}

void GasKinetics::updateYDerivativeFactors()
{
    //update_rates_C();
}

void GasKinetics::update_rates_C()
{
    thermo().getActivityConcentrations(m_conc.data());
    doublereal ctot = thermo().molarDensity();

    // 3-body reactions
    if (!concm_3b_values.empty()) {
        m_3b_concm.update(m_conc, ctot, concm_3b_values.data());
    }

    // Falloff reactions
    if (!concm_falloff_values.empty()) {
        m_falloff_concm.update(m_conc, ctot, concm_falloff_values.data());
    }

    // P-log reactions
    if (m_plog_rates.nReactions()) {
        double logP = log(thermo().pressure());
        m_plog_rates.update_C(&logP);
    }

    // Chebyshev reactions
    if (m_cheb_rates.nReactions()) {
        double log10P = log10(thermo().pressure());
        m_cheb_rates.update_C(&log10P);
    }

    m_ROP_ok = false;
}

void GasKinetics::updateKc()
{
    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        m_rkcn[irxn] = std::min(exp(m_rkcn[irxn]*rrt - m_dn[irxn]*m_logStandConc),
                                BigNumber);
    }

    for (size_t i = 0; i != m_irrev.size(); ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

void GasKinetics::getEquilibriumConstants(doublereal* kc)
{
    update_rates_T();
    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reactions
    getReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-m_rkcn[i]*rrt + m_dn[i]*m_logStandConc);
    }

    // force an update of T-dependent properties, so that m_rkcn will
    // be updated before it is used next.
    m_temp = 0.0;
}

void GasKinetics::computeFalloffPrefactors()
{
    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        falloff_pr[i] = concm_falloff_values[i] * m_rfn_low[i] / (m_rfn_high[i] + SmallNumber);
        AssertFinite(falloff_pr[i], "GasKinetics::processFalloffReactions",
                     "pr[{}] is not finite.", i);
    }

    m_falloffn.pr_to_falloff(falloff_pr.data(), falloff_work.data());
}
void GasKinetics::processFalloffReactions()
{
    computeFalloffPrefactors();
    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        if (reactionType(m_fallindx[i]) == FALLOFF_RXN) {
            m_ropf[m_fallindx[i]] = falloff_pr[i] * m_rfn_high[i];
        } else { // CHEMACT_RXN
            m_ropf[m_fallindx[i]] = falloff_pr[i] * m_rfn_low[i];
        }
    }
}

void GasKinetics::processFalloffReactionDerivatives()
{
    // use m_ropr for temporary storage of reduced pressure

    double invT = 1.0/thermo().temperature();
    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        falloff_pr[i] = concm_falloff_values[i] * m_rfn_low[i] / (m_rfn_high[i] + SmallNumber);
        AssertFinite(falloff_pr[i], "GasKinetics::processFalloffReactions",
                     "pr[{}] is not finite.", i);
        m_dPr_dT[i] = (m_rfn_low_dT[i] - m_rfn_high_dT[i] - invT) * falloff_pr[i];
    }

    auto W = thermo().meanMolecularWeight();
    const auto weights = thermo().molecularWeights();
    auto density = thermo().density();
    vector_fp work(m_falloff_low_rates.nReactions());
    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        work[i] = m_rfn_low[i] / (m_rfn_high[i] + SmallNumber);
    }

    if (!m_dPr_dY.data().size() && m_falloff_low_rates.nReactions()) 
        m_dPr_dY.resize(m_falloff_low_rates.nReactions(), m_kk);
    for (size_t j = 0; j < m_kk; j++) {
        m_falloff_concm.getEfficiencies(j, m_dPr_dY.ptrColumn(j)); 
        auto den_by_wt = density/weights[j];
        for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
            // TODO: For non-const vol reactors add the first term of pyjac (85)
            m_dPr_dY(i, j) *= work[i] * den_by_wt; 
        }
    }
    
    m_falloffn.dc_dPr(falloff_pr.data(), m_dc_dPr.data());

    m_falloffn.dF_dPr(falloff_pr.data(), falloff_work.data(), work.data());
    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        m_dc_dPr[i] += work[i];
    }

    if (!m_dc_dY.data().size() && m_falloff_low_rates.nReactions()) 
        m_dc_dY.resize(m_falloff_low_rates.nReactions(), m_kk);
    for (size_t j = 0; j < m_kk; j++) {
        for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
            m_dc_dY(i, j) = m_dc_dPr[i] *  m_dPr_dY(i,j);
        }
    }
}

void GasKinetics::updateROP()
{
    update_rates_C();
    update_rates_T();
    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    m_ropf = m_rfn;

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm.multiply(m_ropf.data(), concm_3b_values.data());
    }

    if (m_falloff_high_rates.nReactions()) {
        processFalloffReactions();
    }

    for (size_t i = 0; i < nReactions(); i++) {
        // Scale the forward rate coefficient by the perturbation factor
        m_ropf[i] *= m_perturb[i];
        // For reverse rates computed from thermochemistry, multiply the forward
        // rate coefficients by the reciprocals of the equilibrium constants
        m_ropr[i] = m_ropf[i] * m_rkcn[i];
    }

    // multiply ropf by concentration products
    m_reactantStoich.multiply(m_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    m_revProductStoich.multiply(m_conc.data(), m_ropr.data());

    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    for (size_t i = 0; i < m_rfn.size(); i++) {
        AssertFinite(m_rfn[i], "GasKinetics::updateROP",
                     "m_rfn[{}] is not finite.", i);
        AssertFinite(m_ropf[i], "GasKinetics::updateROP",
                     "m_ropf[{}] is not finite.", i);
        AssertFinite(m_ropr[i], "GasKinetics::updateROP",
                     "m_ropr[{}] is not finite.", i);
    }
    m_ROP_ok = true;
}

void GasKinetics::updateROPDerivatives(bool constPressure)
{
    updateROP();
    updateTDerivativeFactors();
    processFalloffReactionDerivatives();

    // Temperature derivatives
    doublereal invT = 1.0/m_temp;
    if (constPressure){
        for (size_t i = 0; i < nReactions(); i++){
            m_dNetROPdT[i] = m_ropf[i] * 
                (m_rfn_dTMult[i] - m_reactant_stoichsum[i] * invT);
            m_dNetROPdT[i] -= m_ropr[i] * 
                ((m_rfn_dTMult[i] - m_product_stoichsum[i] * invT) - m_dBdT[i]);
        }
    } else {
        for (size_t i = 0; i < nReactions(); i++){
            m_dNetROPdT[i] = m_ropf[i] * m_rfn_dTMult[i];// * invT;
            m_dNetROPdT[i] -= m_ropr[i] * (m_rfn_dTMult[i]  - m_dBdT[i]);
        }
    }

    // Process 3b reactions 
    if (!concm_3b_values.empty()) {
        m_3b_concm.add(m_dNetROPdT.data(), m_ropnet.data(), -invT);
    }
    
    // TODO: Process T derivatives of falloff reactions


    // Mass fraction derivatives
    updateYDerivativeFactors();

    if (!m_dRevROPdY.data().size()){   //TODO: Push these to class initialization
        m_dRevROPdY.resize(nReactions(), m_kk);
    }
    if (!m_dNetROPdY.data().size()){
        m_dNetROPdY.resize(nReactions(), m_kk);
    }

    auto W = thermo().meanMolecularWeight();
    const auto weights = thermo().molecularWeights();
    auto den = thermo().density();

    vector_fp rfn_tmp(m_rfn);
    // Process Fallof reaction rate constants
    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        if (reactionType(m_fallindx[i]) == FALLOFF_RXN) {
            rfn_tmp[m_fallindx[i]] = m_rfn_high[i];
        } else { // CHEMACT_RXN
            rfn_tmp[m_fallindx[i]] = m_rfn_low[i];
        }
    }
    
    computeFalloffPrefactors();
    for (size_t j = 0; j < m_kk; j++){    // Fwd part of Eq. 76 of pyjac fwd  
        m_dNetROPdY.setColumn(j, rfn_tmp.data());
        m_reactantStoich.derivative_multiply(
                m_conc.data(), m_dNetROPdY.ptrColumn(j), j);
        auto den_by_wght = den / weights[j];
        auto muwght_by_wght = W / weights[j] ;
        for (size_t i = 0; i < nReactions(); i++){
            m_dNetROPdY(i, j) *=  den_by_wght;
        }
        if (!concm_3b_values.empty()) {     // Eq. 58 of pyjac net
            m_3b_concm.multiply(m_dNetROPdY.ptrColumn(j), concm_3b_values.data());
        }
        // Process Falloff prefactors
        for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
            m_dNetROPdY(m_fallindx[i], j) *= falloff_pr[i];
        }

        if (constPressure){
            for (size_t i = 0; i < nReactions(); i++){
                m_dNetROPdY(i, j) -= m_ropf[i] * m_reactant_stoichsum[i] * 
                                     muwght_by_wght;
            }
        }
    }

    for (size_t j = 0; j < m_kk; j++){    // Rev part of Eq. 76 of pyjac 
        m_dRevROPdY.setColumn(j, rfn_tmp.data());
        m_revProductStoich.derivative_multiply(
                m_conc.data(), m_dRevROPdY.ptrColumn(j), j);
        auto den_by_wght = den / weights[j];
        auto muwght_by_wght = W / weights[j] ;
        for (size_t i = 0; i < nReactions(); i++){
            m_dRevROPdY(i, j) *=  den_by_wght * m_rkcn[i];
        }
        // Eq. 58 of pyjac net
        if (!concm_3b_values.empty()) {  // Process 3b reactions
            m_3b_concm.multiply(m_dRevROPdY.ptrColumn(j), concm_3b_values.data());
        }
        // Process Falloff prefactors
        for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
            m_dRevROPdY(m_fallindx[i], j) *= falloff_pr[i];
        }
        if (constPressure){
            for (size_t i = 0; i < nReactions(); i++){
                m_dRevROPdY(i, j) -= m_ropr[i] * m_product_stoichsum[i] * 
                                     muwght_by_wght;
            }
        }
    }

    vector_fp work2;  // Temporary arrays
    if (!concm_3b_values.empty()) 
        work2.resize(concm_3b_values.size());

    for (size_t j = 0; j < m_kk; j++){    
        for (size_t i = 0; i < nReactions(); i++){
            m_dNetROPdY(i, j) -= m_dRevROPdY(i, j);
        }
    }

    if (!concm_3b_values.empty()) {     // Eq. 58 of pyjac net
        for (size_t j = 0; j < m_kk; j++){    
            //m_3b_concm.multiply(m_dNetROPdY.ptrColumn(j), concm_3b_values.data());

            m_3b_concm.getEfficiencies(j, work2.data()); 
            auto den_by_wtj = den / weights[j];
            for (size_t i = 0; i < concm_3b_values.size(); i++) {
                work2[i] *= den_by_wtj / concm_3b_values[i];     // 2nd term in Eq.(80) of pyjac
            }
            if (constPressure){
                auto W_by_Wj = W / weights[j] ;
                for (size_t i = 0; i < concm_3b_values.size(); i++) {
                    work2[i] -= W_by_Wj; //* concm_3b_values[i];   // 1st term in Eq.(80) of pyjac
                }
            }

            m_3b_concm.add(m_dNetROPdY.ptrColumn(j), m_ropnet.data(), work2.data());
        }
    }

        // Process Falloff reaction values
    for (size_t j = 0; j < m_kk; j++){    
        for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
            m_dNetROPdY(m_fallindx[i], j) += m_dc_dY(i, j) * m_ropnet[m_fallindx[i]];
        }
    }

    /*if (m_falloff_high_rates.nReactions()) {
        processFalloffReactionDerivatives();
    }*/

}

void GasKinetics::getFwdRateConstants(doublereal* kfwd)
{
    update_rates_C();
    update_rates_T();

    // copy rate coefficients into ropf
    m_ropf = m_rfn;

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm.multiply(m_ropf.data(), concm_3b_values.data());
    }

    if (m_falloff_high_rates.nReactions()) {
        processFalloffReactions();
    }

    for (size_t i = 0; i < nReactions(); i++) {
        // multiply by perturbation factor
        kfwd[i] = m_ropf[i] * m_perturb[i];
    }
}

bool GasKinetics::addReaction(shared_ptr<Reaction> r)
{
    // operations common to all reaction types
    bool added = BulkKinetics::addReaction(r);
    if (!added) {
        return false;
    }

    switch (r->reaction_type) {
    case ELEMENTARY_RXN:
        addElementaryReaction(dynamic_cast<ElementaryReaction&>(*r));
        break;
    case THREE_BODY_RXN:
        addThreeBodyReaction(dynamic_cast<ThreeBodyReaction&>(*r));
        break;
    case FALLOFF_RXN:
    case CHEMACT_RXN:
        addFalloffReaction(dynamic_cast<FalloffReaction&>(*r));
        break;
    case PLOG_RXN:
        addPlogReaction(dynamic_cast<PlogReaction&>(*r));
        break;
    case CHEBYSHEV_RXN:
        addChebyshevReaction(dynamic_cast<ChebyshevReaction&>(*r));
        break;
    default:
        throw CanteraError("GasKinetics::addReaction",
            "Unknown reaction type specified: {}", r->reaction_type);
    }
    return true;
}

void GasKinetics::addFalloffReaction(FalloffReaction& r)
{
    // install high and low rate coeff calculators and extend the high and low
    // rate coeff value vectors
    size_t nfall = m_falloff_high_rates.nReactions();
    m_falloff_high_rates.install(nfall, r.high_rate);
    m_rfn_high.push_back(0.0);
    m_falloff_low_rates.install(nfall, r.low_rate);
    m_rfn_low.push_back(0.0);
    falloff_pr.push_back(0.0);

    // add this reaction number to the list of falloff reactions
    m_fallindx.push_back(nReactions()-1);
    m_rfallindx[nReactions()-1] = nfall;

    // Extend the derivative vectors associated with falloff reactions
    m_rfn_low_dT.push_back(0.0);
    m_rfn_high_dT.push_back(0.0);
    falloff_work_dT.push_back(0.0);
    m_dPr_dT.push_back(0.0);
    m_dc_dT.push_back(0.0);
    m_dc_dPr.push_back(0.0);

    // install the enhanced third-body concentration calculator
    map<size_t, double> efficiencies;
    for (const auto& eff : r.third_body.efficiencies) {
        size_t k = kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            efficiencies[k] = eff.second;
        } else if (!m_skipUndeclaredThirdBodies) {
            throw CanteraError("GasKinetics::addFalloffReaction", "Found "
                "third-body efficiency for undefined species '" + eff.first +
                "' while adding reaction '" + r.equation() + "'");
        }
    }
    m_falloff_concm.install(nfall, efficiencies,
                            r.third_body.default_efficiency);
    concm_falloff_values.resize(m_falloff_concm.workSize());

    // install the falloff function calculator for this reaction
    m_falloffn.install(nfall, r.reaction_type, r.falloff);
    falloff_work.resize(m_falloffn.workSize());
}

void GasKinetics::addThreeBodyReaction(ThreeBodyReaction& r)
{
    m_rates.install(nReactions()-1, r.rate);
    map<size_t, double> efficiencies;
    for (const auto& eff : r.third_body.efficiencies) {
        size_t k = kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            efficiencies[k] = eff.second;
        } else if (!m_skipUndeclaredThirdBodies) {
            throw CanteraError("GasKinetics::addThreeBodyReaction", "Found "
                "third-body efficiency for undefined species '" + eff.first +
                "' while adding reaction '" + r.equation() + "'");
        }
    }
    m_3b_concm.install(nReactions()-1, efficiencies,
                       r.third_body.default_efficiency);
    concm_3b_values.resize(m_3b_concm.workSize());
}

void GasKinetics::addPlogReaction(PlogReaction& r)
{
    m_plog_rates.install(nReactions()-1, r.rate);
}

void GasKinetics::addChebyshevReaction(ChebyshevReaction& r)
{
    m_cheb_rates.install(nReactions()-1, r.rate);
}

void GasKinetics::modifyReaction(size_t i, shared_ptr<Reaction> rNew)
{
    // operations common to all reaction types
    BulkKinetics::modifyReaction(i, rNew);

    switch (rNew->reaction_type) {
    case ELEMENTARY_RXN:
        modifyElementaryReaction(i, dynamic_cast<ElementaryReaction&>(*rNew));
        break;
    case THREE_BODY_RXN:
        modifyThreeBodyReaction(i, dynamic_cast<ThreeBodyReaction&>(*rNew));
        break;
    case FALLOFF_RXN:
    case CHEMACT_RXN:
        modifyFalloffReaction(i, dynamic_cast<FalloffReaction&>(*rNew));
        break;
    case PLOG_RXN:
        modifyPlogReaction(i, dynamic_cast<PlogReaction&>(*rNew));
        break;
    case CHEBYSHEV_RXN:
        modifyChebyshevReaction(i, dynamic_cast<ChebyshevReaction&>(*rNew));
        break;
    default:
        throw CanteraError("GasKinetics::modifyReaction",
            "Unknown reaction type specified: {}", rNew->reaction_type);
    }

    // invalidate all cached data
    m_ROP_ok = false;
    m_temp += 0.1234;
    m_pres += 0.1234;
}

void GasKinetics::modifyThreeBodyReaction(size_t i, ThreeBodyReaction& r)
{
    m_rates.replace(i, r.rate);
}

void GasKinetics::modifyFalloffReaction(size_t i, FalloffReaction& r)
{
    size_t iFall = m_rfallindx[i];
    m_falloff_high_rates.replace(iFall, r.high_rate);
    m_falloff_low_rates.replace(iFall, r.low_rate);
    m_falloffn.replace(iFall, r.falloff);
}

void GasKinetics::modifyPlogReaction(size_t i, PlogReaction& r)
{
    m_plog_rates.replace(i, r.rate);
}

void GasKinetics::modifyChebyshevReaction(size_t i, ChebyshevReaction& r)
{
    m_cheb_rates.replace(i, r.rate);
}

void GasKinetics::init()
{
    BulkKinetics::init();
    m_logp_ref = log(thermo().refPressure()) - log(GasConstant);
}

void GasKinetics::invalidateCache()
{
    BulkKinetics::invalidateCache();
    m_pres += 0.13579;
}

}
