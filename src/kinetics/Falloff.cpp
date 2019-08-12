/**
 *  @file Falloff.cpp Definitions for member functions of classes derived from
 *      Falloff
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/kinetics/Falloff.h"

using namespace std;

namespace Cantera
{

void Falloff::init(const vector_fp& c)
{
    if (c.size() != 0) {
        throw CanteraError("Falloff::init",
            "Incorrect number of parameters. 0 required. Received {}.",
            c.size());
    }
}

void Troe::init(const vector_fp& c)
{
    if (c.size() != 3 && c.size() != 4) {
        throw CanteraError("Troe::init",
            "Incorrect number of parameters. 3 or 4 required. Received {}.",
            c.size());
    }
    m_a = c[0];
    if (std::abs(c[1]) < SmallNumber) {
        m_rt3 = std::numeric_limits<double>::infinity();
    } else {
        m_rt3 = 1.0 / c[1];
    }

    if (std::abs(c[2]) < SmallNumber) {
        m_rt1 = std::numeric_limits<double>::infinity();
    } else {
        m_rt1 = 1.0 / c[2];
    }

    if (c.size() == 4) {
        m_t2 = c[3];
    }
}

void Troe::updateTemp(double T, double* work) const
{
    double Fcent = (1.0 - m_a) * exp(-T*m_rt3) + m_a * exp(-T*m_rt1);
    if (m_t2) {
        Fcent += exp(- m_t2 / T);
    }
    *work = log10(std::max(Fcent, SmallNumber));
}

void Troe::updateTempDerivative(double T, double* work) const
{
    double dFcent = -(1.0 - m_a) * m_rt3 * exp(-T*m_rt3) 
                   - m_a * m_rt1 * exp(-T*m_rt1);
    if (m_t2) {
        auto invT = 1.0 / T;
        auto t2byT = m_t2 * invT;
        dFcent += t2byT * invT * exp(-t2byT);
    }
    *work = std::max(dFcent, SmallNumber);
}

double Troe::F(double pr, const double* work) const
{
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

//TODO: Change the name (get rid of Fcent and use generic one applicable for all falloff classes
double Troe::dF_dFcent(double pr, const double* work) const
{
    auto log10Fcent = *work;
    auto Fcent = pow(10.0, log10Fcent);
    auto logpr = log10(std::max(pr, SmallNumber)); 
    double invFcent = 1.0/(Fcent);
    double A = logpr - 0.67 * log10Fcent - 0.4;
    double B = 0.806 - 1.1762 * log10Fcent - 0.14 * logpr;
    double Bcube = B * B * B;
    double AbyBsq = A * A / (B * B);
    double One_1pAbyBsq = 1.0/(1 +  AbyBsq);
    return One_1pAbyBsq * (            // Fi removed
                1.0 / Fcent  - log10Fcent * 2*A/Bcube * One_1pAbyBsq * 
                               invFcent * (-0.67 * B + 1.1762 * A));
}

// Eq.(94) of pyjac.
// In the final value F_i is omitted because it gets cancelled in Eq.(87)
double Troe::dF_dPr(double pr, const double* work) const
{
    auto log10Fcent = *work;
    auto logpr = log10(std::max(pr, SmallNumber)); 
    double A = logpr - 0.67 * log10Fcent - 0.4;
    double B = 0.806 - 1.1762 * log10Fcent - 0.14 * logpr;
    double Bcube = B * B * B;
    double AbyBsq = A * A / (B * B);
    double One_1pAbyBsq = 1.0/(1 +  AbyBsq);
    return -One_1pAbyBsq * One_1pAbyBsq * log10Fcent * 2*A *  // Fi is removed
                 (B + 0.14 * A) / (Bcube * pr) ;
}

/*
double Troe::dFdT(double pr, double dprdT, double logFcent, double dFcent) const
{
    throw CanteraError("Troe::dFdT", "Not implemented."); 
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}

double Troe::dFdY(double pr, const double* work, size_t j) const
{
    throw CanteraError("Troe::dFdY", "Not implemented."); 
    double lpr = log10(std::max(pr,SmallNumber));
    double cc = -0.4 - 0.67 * (*work);
    double nn = 0.75 - 1.27 * (*work);
    double f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
    double lgf = (*work) / (1.0 + f1 * f1);
    return pow(10.0, lgf);
}
*/
void Troe::getParameters(double* params) const {
    params[0] = m_a;
    params[1] = 1.0/m_rt3;
    params[2] = 1.0/m_rt1;
    params[3] = m_t2;
}

void SRI::init(const vector_fp& c)
{
    if (c.size() != 3 && c.size() != 5) {
        throw CanteraError("SRI::init",
            "Incorrect number of parameters. 3 or 5 required. Received {}.",
            c.size());
    }

    if (c[2] < 0.0) {
        throw CanteraError("SRI::init()",
                           "m_c parameter is less than zero: {}", c[2]);
    }
    m_a = c[0];
    m_b = c[1];
    m_c = c[2];

    if (c.size() == 5) {
        if (c[3] < 0.0) {
            throw CanteraError("SRI::init()",
                               "m_d parameter is less than zero: {}", c[3]);
        }
        m_d = c[3];
        m_e = c[4];
    } else {
        m_d = 1.0;
        m_e = 0.0;
    }
}

void SRI::updateTemp(double T, double* work) const
{
    *work = m_a * exp(- m_b / T);
    if (m_c != 0.0) {
        *work += exp(- T/m_c);
    }
    work[1] = m_d * pow(T,m_e);
}

void SRI::updateTempDerivative(double T, double* work) const
{
    throw CanteraError("SRI::updateTempDerivative", "Not implemented."); 
}

double SRI::F(double pr, const double* work) const
{
    double lpr = log10(std::max(pr,SmallNumber));
    double xx = 1.0/(1.0 + lpr*lpr);
    return pow(*work, xx) * work[1];
}

double SRI::dF_dFcent(double pr, const double* work) const
{
    throw CanteraError("SRI::dF_dFcent", "Not implemented."); 
}

double SRI::dF_dPr(double pr, const double* work) const
{
    double lpr = log10(max(pr,SmallNumber));
    double xx = 1.0/(1.0 + lpr*lpr);
    return -xx * xx * log(max(*work, SmallNumber)) * 2 * lpr / (pr  * log(10));
}

/*
double SRI::dFdT(double pr, double dprdT, double logFcent, double dFcent) const
{
    throw CanteraError("SRI::dFdT", "Not implemented."); 
}

double SRI::dFdY(double pr, const double* work, size_t j) const
{
    throw CanteraError("SRI::dFdY", "Not implemented."); 
}
*/

void SRI::getParameters(double* params) const
{
    params[0] = m_a;
    params[1] = m_b;
    params[2] = m_c;
    params[3] = m_d;
    params[4] = m_e;
}

}
