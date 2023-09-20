// File:    funclist.cpp
// Author:  Brian Vanderburg II
// Purpose: Function list
//------------------------------------------------------------------------------


// Includes
#include <new>

#include "funclist.h"
#include "except.h"
#include "node.h"

#include <autodiff/reverse/var.hpp>

using namespace std;
using namespace ExprEval;

// Function factory
//------------------------------------------------------------------------------

// Constructor
template<class Value>
FunctionFactory<Value>::FunctionFactory()
{
}
    
// Destructor
template<class Value>
FunctionFactory<Value>::~FunctionFactory()
{
}
    
// Create
template<class Value>
FunctionNode<Value> *FunctionFactory<Value>::Create(Expression<Value> *expr)
{
    FunctionNode<Value> *n = DoCreate(expr);
    if(n)
        n->m_factory = this;
        
    return n;
}
        

// Function list
//------------------------------------------------------------------------------

// Constructor
template<class Value>
FunctionList<Value>::FunctionList()
{
}
                    
// Destructor
template<class Value>
FunctionList<Value>::~FunctionList()
{
    // Free function factories
    Clear();
}
            
// Add factory to list
template<class Value>
void FunctionList<Value>::Add(FunctionFactory<Value> *factory)
{
    // Check it
    if(factory == 0)
        throw(NullPointerException("FunctionList::Add"));
    
    // Make sure it does not exist
    size_type pos;
    
    for(pos  = 0; pos < m_functions.size(); pos++)
    {
        if(m_functions[pos]->GetName() == factory->GetName())
            throw(AlreadyExistsException(factory->GetName()));
    }
    
    m_functions.push_back(factory);
}
    
// Create a node for a function
template<class Value>
FunctionNode<Value> *FunctionList<Value>::Create(const string &name, Expression<Value> *expr)
{
    // Make sure pointer exists
    if(expr == 0)
        throw(NullPointerException("FunctionList::Create"));
    
    size_type pos;
    
    for(pos = 0; pos < m_functions.size(); pos++)
    {
        if(m_functions[pos]->GetName() == name)
        {
            // Found it
            return m_functions[pos]->Create(expr);
        }
    }
        
    // Not found
    return 0;
}
            
// FunctionList::AddDefaultFunctions is located in func.cpp
// along with the default function factories            
          
// Free function list
template<class Value>
void FunctionList<Value>::Clear()
{
    size_type pos;
    
    for(pos = 0; pos < m_functions.size(); pos++)
    {
        delete m_functions[pos];
    }
}

#define INSTANCE(...) \
template class FunctionFactory<__VA_ARGS__>; \
template class FunctionList<__VA_ARGS__>;

namespace ExprEval {
INSTANCE(double)
INSTANCE(autodiff::var)
}
