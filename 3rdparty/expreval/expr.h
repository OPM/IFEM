// File:    expr.h
// Author:  Brian Vanderburg II
// Purpose: Expression object
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_EXPR_H
#define __EXPREVAL_EXPR_H

// Includes
#include <string>

// Part of expreval namespace
namespace ExprEval
{
    // Forward declarations
    template<class Value> class ValueList;
    template<class Value> class FunctionList;
    template<class Value> class Node;

    // Expression class
    //--------------------------------------------------------------------------
    template<class Value>
    class Expression
    {
    public:
        Expression();
        virtual ~Expression();

        // Variable list
        void SetValueList(ValueList<Value> *vlist);
        ValueList<Value> *GetValueList() const;

        // Function list
        void SetFunctionList(FunctionList<Value> *flist);
        FunctionList<Value> *GetFunctionList() const;

        // Abort control
        virtual bool DoTestAbort();
        void TestAbort(bool force = false);
        void SetTestAbortCount(unsigned long count);

        // Parse an expression
        void Parse(const ::std::string &exstr);

        // Clear an expression
        void Clear();

        // Evaluate expression
        Value Evaluate();

    protected:
        ValueList<Value> *m_vlist;
        FunctionList<Value> *m_flist;
        Node<Value> *m_expr;
        unsigned long m_abortcount;
        unsigned long m_abortreset;
    };



} // namespace ExprEval

#endif // __EXPREVAL_EXPR_H
