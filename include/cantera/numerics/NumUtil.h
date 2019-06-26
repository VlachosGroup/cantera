//! @file NumUtil.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_NUMUTIL_H
#define CT_NUMUTIL_H

namespace Cantera
{

const int c_NONE = 0;
const int c_GE_ZERO = 1;
const int c_GT_ZERO = 2;
const int c_LE_ZERO = -1;
const int c_LT_ZERO = -2;

const int DIAG = 1;
const int DENSE = 2;
const int NOJAC = 4;
const int JAC = 8;
const int GMRES = 16;
const int BAND = 32;



inline bool checkFlag(const int constraintFlag)
 {
     auto cflag = constraintFlag;
     bool valid_cflag = false;
     if (cflag == c_NONE || cflag == c_GE_ZERO || cflag == c_GT_ZERO ||
         cflag == c_LE_ZERO || cflag == c_LT_ZERO)
         valid_cflag = true;
     return valid_cflag;
 }
 

}

#endif
