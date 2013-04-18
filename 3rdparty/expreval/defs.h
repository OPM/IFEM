// File:    defs.h
// Author:  Brian Vanderburg II
// Purpose: Definitions for ExprEval
//------------------------------------------------------------------------------

#ifndef __EXPREVAL_DEFS_H
#define __EXPREVAL_DEFS_H

namespace ExprEval
{
    // constants
    const double EXPREVAL_PI = 3.14159265358979323846;
    const double EXPREVAL_E = 2.7182818284590452354;
    
} // namespace ExprEval

// CHANGED from original code: due to auto_ptr deprecation in c++0x
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && !defined(__INTEL_COMPILER)
#define aptr(x) std::unique_ptr<x>
#else
#define aptr(x) std::auto_ptr<x>
#endif


#endif // __EXPREVAL_DEFS_H

