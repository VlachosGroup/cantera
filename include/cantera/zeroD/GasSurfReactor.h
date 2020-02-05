//! @file GasSurfReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_GASSURFREACTOR_H
#define CT_GASSURFREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class GasSurfReactor is a class for stirred reactors that is specifically
 * optimized for ideal gases. In this formulation, temperature replaces the
 * total internal energy as a state variable.
 */
class GasSurfReactor: public Reactor 
{
public:
    GasSurfReactor () {}

    virtual int type() const {
        return GasSurfReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void evalJacEqs(doublereal t, doublereal* y, doublereal* ydot,
                            Array2D* jac);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "temperature", the name of a homogeneous phase species, or the
    //! name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);

//protected:
//    vector_fp m_uk; //!< Species molar internal energies
};

}

#endif
