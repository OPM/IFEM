// File:    vallist.cpp
// Author:  Brian Vanderburg II
// Purpose: Value list used for variable and constants
//------------------------------------------------------------------------------


// Includes
#include <new>
#include <memory>

#include "autodiff/reverse/var/var.hpp"
#include "defs.h"
#include "vallist.h"
#include "except.h"

using namespace std;
using namespace ExprEval;

// ValueListItem
//------------------------------------------------------------------------------

// Constructor for internal value
template<class Value>
ValueListItem<Value>::ValueListItem(const string &name, Value def, bool constant)
{
    m_name = name;
    m_constant = constant;
    m_value = m_def = def;
    m_ptr = 0;
}

// Constructor for external value
template<class Value>
ValueListItem<Value>::ValueListItem(const string &name, Value *ptr, Value def, bool constant)
{
    m_name = name;
    m_constant = constant;
    m_value = m_def = def;    
    m_ptr = ptr;
    
    if(m_ptr)
        *m_ptr = def;
    else
        throw(NullPointerException("ValueListItem::ValueListItem"));        
}    
    
// Get the name
template<class Value>
const string& ValueListItem<Value>::GetName() const
{
    return m_name;
}
        
// Return if it is constant
template<class Value>
bool ValueListItem<Value>::IsConstant() const
{
    return m_constant;
}
    
// Get value address
template<class Value>
Value* ValueListItem<Value>::GetAddress()
{
    return m_ptr ? m_ptr : &m_value;
} 
    
// Reset to default value
template<class Value>
void ValueListItem<Value>::Reset()
{
    if(m_ptr)
        *m_ptr = m_def;
    else
        m_value = m_def;
}    

    
// ValueList
//------------------------------------------------------------------------------

// Constructor
template<class Value>
ValueList<Value>::ValueList()
{
}
    
// Destructor
template<class Value>
ValueList<Value>::~ValueList()
{
    Clear();
}    
    
// Add value to list
template<class Value>
void ValueList<Value>::Add(const string &name, Value def, bool constant)
{
    // Ensure value does not already exist
    if(GetAddress(name))
    {
        throw(AlreadyExistsException(name));
    }
    else
    {
        // Create value
        aptr(ValueListItem<Value>) i(new ValueListItem<Value>(name, def, constant));
        
        // Add value to list
        m_values.push_back(i.get());
        i.release();
    }
}
    
// Add an external value to the list
template<class Value>
void ValueList<Value>::AddAddress(const string &name, Value* ptr, Value def, bool constant)
{
    if(GetAddress(name))
    {
        throw(AlreadyExistsException(name));
    }
    else if(ptr == 0)
    {
        throw(NullPointerException("ValueList::AddAddress"));
    }
    else
    {
        // Create value
        aptr(ValueListItem<Value>) i(new ValueListItem<Value>(name, ptr, def, constant));
        
        // Add value to list
        m_values.push_back(i.get());
        i.release();
    }
}    
    
// Get the address of the value, internal or external
template<class Value>
Value* ValueList<Value>::GetAddress(const string &name) const
{
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
    {
        if(m_values[pos]->GetName() == name)
        {
            // The name matches
            return m_values[pos]->GetAddress();
        }
    }
        
    // No item found
    return 0;
}    
    
// Is the value a constant
template<class Value>
bool ValueList<Value>::IsConstant(const string &name) const
{
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
    {
        if(m_values[pos]->GetName() == name && m_values[pos]->IsConstant())
        {
            return true;
        }
    }
        
    return false;
}   
    
// Number of values in the list
template<class Value>
typename ValueList<Value>::size_type ValueList<Value>::Count() const
{
    return m_values.size();
}
    
// Get an item
template<class Value>
void ValueList<Value>::Item(size_type pos, string *name, Value* value) const
{
    if(name)
        *name = m_values[pos]->GetName();
        
    if(value)
        *value = *(m_values[pos]->GetAddress());
}
    
// Add some default values
template<class Value>
void ValueList<Value>::AddDefaultValues()
{
    // Math constant 'e'
    Add("E", EXPREVAL_E, true);
    
    // Math constant PI
    Add("PI", EXPREVAL_PI, true);
}
    
// Reset values
template<class Value>
void ValueList<Value>::Reset()
{
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
    {
        m_values[pos]->Reset();
    }
}
    
// Free values
template<class Value>
void ValueList<Value>::Clear()
{
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
    {
        delete m_values[pos];
    }
}

namespace ExprEval {
template class ValueList<double>;
template class ValueList<autodiff::var>;
}
