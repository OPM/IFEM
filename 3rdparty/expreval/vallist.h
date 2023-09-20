// File:    vallist.h
// Author:  Brian Vanderburg II
// Purpose: Value list used for variable and constants
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_VALLIST_H
#define __EXPREVAL_VALLIST_H

// Includes
#include <string>
#include <vector>

// Part of expreval namespace
namespace ExprEval
{
// Value list item
//--------------------------------------------------------------------------
template<class Value>
    class ValueListItem
    {
    public:
        ValueListItem(const ::std::string &name, Value def = 0.0, bool constant = false);
        ValueListItem(const ::std::string &name, Value *ptr, Value def = 0.0, bool constant = false);
        
        const ::std::string &GetName() const;
        bool IsConstant() const;
        
        Value *GetAddress();
        void Reset();
        
    private:
        ::std::string m_name; // Name of value
        bool m_constant; // Value is constant
        
        Value m_value; // Internal value (if ptr == 0)
        Value *m_ptr; // Pointer to extern value if not 0
            
        Value m_def; // Default value when reset
    };
        
    
        
    // Value list
    //--------------------------------------------------------------------------
    template<class Value>
    class ValueList
    {
    public:
        using ValueVector = std::vector<ValueListItem<Value>*>;
        using size_type = typename ValueVector::size_type;
                
        ValueList();
        ~ValueList();
        
        // Add variable or constant to the list
        void Add(const ::std::string &name, Value def = 0.0, bool constant = false);
        
        // Add an external variable or constant to the list
        void AddAddress(const ::std::string &name, Value *ptr, Value def = 0.0, bool constant = false);
        
        // Get the address of a variable or constant, internal or external
        Value *GetAddress(const ::std::string &name) const;
        
        // Is the value constant
        bool IsConstant(const ::std::string &name) const;
        
        // Enumerate values
        size_type Count() const;
        void Item(size_type pos, ::std::string *name = 0, Value *value = 0) const;
        
        // Initialize some default values (math constants)
        void AddDefaultValues();
        
        // Reset items to default values (constants are not changed)
        void Reset();
        
        // Free items and clear memory
        void Clear();
    
    private:
        ::std::vector<ValueListItem<Value>*> m_values;
    };
};

#endif // __EXPREVAL_VALLIST_H

