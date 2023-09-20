// File:    funclist.h
// Author:  Brian Vanderburg II
// Purpose: Function list
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_FUNCLIST_H
#define __EXPREVAL_FUNCLIST_H

// Includes
#include <string>
#include <vector>

// Part of expreval namespace
namespace ExprEval
{
    // Forward declarations
    template<class Value> class FunctionNode;
    template<class Value> class Expression;
    
    // Function factory
    //--------------------------------------------------------------------------
    template<class Value>
    class FunctionFactory
    {
    public:
        FunctionFactory();
        virtual ~FunctionFactory();
        
        virtual ::std::string GetName() const = 0;
        virtual FunctionNode<Value> *DoCreate(Expression<Value> *expr) = 0;
        
        FunctionNode<Value> *Create(Expression<Value> *expr);
    };
        
    // Function list
    //--------------------------------------------------------------------------
    template<class Value>
    class FunctionList
    {
    public:
        using FactoryVec = std::vector<FunctionFactory<Value>*>;
        using size_type = typename FactoryVec::size_type;
                
        FunctionList();
        ~FunctionList();
        
        // Add a function factory to the list
        void Add(FunctionFactory<Value> *factory);
        
        // Create a node for a function
        FunctionNode<Value> *Create(const ::std::string &name, Expression<Value> *expr);
        
        // Initialize default functions
        void AddDefaultFunctions();
        
        // Free items and clear memory
        void Clear();
    
    private:
        ::std::vector<FunctionFactory<Value>*> m_functions;
    };
        
}

#endif // __EXPREVAL_FUNCLIST_H

