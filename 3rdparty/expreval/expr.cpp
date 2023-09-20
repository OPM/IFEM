// File:    expr.cpp
// Author:  Brian Vanderburg II
// Purpose: Expression object
//------------------------------------------------------------------------------

// Includes
#include <new>
#include <memory>

#include "defs.h"
#include "expr.h"
#include "parser.h"
#include "node.h"
#include "except.h"

using namespace std;
using namespace ExprEval;

#include <autodiff/reverse/var.hpp>

// Expression object
//------------------------------------------------------------------------------

// Constructor
template<class Value>
Expression<Value>::Expression() : m_vlist(0), m_flist(0), m_expr(0)
{
    m_abortcount = 200000;
    m_abortreset = 200000;
}

// Destructor
template<class Value>
Expression<Value>::~Expression()
{
    // Delete expression nodes
    delete m_expr;
}

// Set value list
template<class Value>
void Expression<Value>::SetValueList(ValueList<Value>* vlist)
{
    m_vlist = vlist;
}

// Get value list
template<class Value>
ValueList<Value>* Expression<Value>::GetValueList() const
{
    return m_vlist;
}

// Set function list
template<class Value>
void Expression<Value>::SetFunctionList(FunctionList<Value> *flist)
{
    m_flist = flist;
}

// Get function list
template<class Value>
FunctionList<Value>* Expression<Value>::GetFunctionList() const
{
    return m_flist;
}

// Test for an abort
template<class Value>
bool Expression<Value>::DoTestAbort()
{
    // Derive a class to test abort
    return false;
}

// Test for an abort
template<class Value>
void Expression<Value>::TestAbort(bool force)
{
    if(force)
    {
        // Test for an abort now
        if(DoTestAbort())
        {
            throw(AbortException());
        }
    }
    else
    {
        // Test only if abort count is 0
        if(m_abortcount == 0)
        {
            // Reset count
            m_abortcount = m_abortreset;

            // Test abort
            if(DoTestAbort())
            {
                throw(AbortException());
            }
        }
        else
        {
            // Decrease abort count
            m_abortcount--;
        }
    }
}

// Set test abort count
template<class Value>
void Expression<Value>::SetTestAbortCount(unsigned long count)
{
    m_abortreset = count;
    if(m_abortcount > count)
        m_abortcount = count;
}

// Parse expression
template<class Value>
void Expression<Value>::Parse(const string &exstr)
{
    // Clear the expression if needed
    if(m_expr)
        Clear();

    // Create parser
    aptr(Parser<Value>) p(new Parser<Value>(this));

    // Parse the expression
    m_expr = p->Parse(exstr);
}

// Clear the expression
template<class Value>
void Expression<Value>::Clear()
{
    delete m_expr;
    m_expr = 0;
}

// Evaluate an expression
template<class Value>
Value Expression<Value>::Evaluate()
{
    if(m_expr)
    {
        return m_expr->Evaluate();
    }
    else
    {
        throw(EmptyExpressionException());
    }
}

namespace ExprEval {
template class Expression<double>;
template class Expression<autodiff::var>;
}
