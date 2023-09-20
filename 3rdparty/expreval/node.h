// File:    node.h
// Author:  Brian Vanderburg II
// Purpose: Expression node
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_NODE_H
#define __EXPREVAL_NODE_H

// Includes
#include <vector>

#include "parser.h"

// Part of expreval namespace
namespace ExprEval
{
    // Forward declarations
    template<class Value> class Expression;
    template<class Value> class FunctionFactory;
    
    // Node class
    //--------------------------------------------------------------------------
    template<class Value>
    class Node
    {
    public:
        explicit Node(Expression<Value>* expr);
        virtual ~Node();
        
        virtual Value DoEvaluate() = 0;
        virtual void Parse(Parser<Value> &parser,
                           typename Parser<Value>::size_type start,
                           typename Parser<Value>::size_type end,
                           typename Parser<Value>::size_type v1 = 0) = 0;
                
        Value Evaluate(); // Calls Expression::TestAbort, then DoEvaluate
        
    protected:
        Expression<Value> *m_expr;
    };
        
    // General function node class
    //--------------------------------------------------------------------------
    template<class Value>
    class FunctionNode : public Node<Value>
    {
    public:
        explicit FunctionNode(Expression<Value> *expr);
        ~FunctionNode();
        
        // Parse nodes and references
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
                
    private:
        // Function factory
        FunctionFactory<Value> *m_factory;
        
        // Argument count
        long m_argMin;
        long m_argMax;
        long m_refMin;
        long m_refMax;

    protected:
        // Set argument count (called in derived constructors)
        void SetArgumentCount(long argMin = 0, long argMax = 0,
                long refMin = 0, long refMax = 0);

        // Function name (using factory)
        ::std::string GetName() const;
    
        // Normal, reference, and data parameters
        ::std::vector<Node<Value>*> m_nodes;
        ::std::vector<Value*> m_refs;

        friend class FunctionFactory<Value>;
    };
        
    // Mulit-expression node
    //--------------------------------------------------------------------------
    template<class Value>
    class MultiNode : public Node<Value>
    {
    public:
        explicit MultiNode(Expression<Value> *expr);
        ~MultiNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value>& parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        ::std::vector<Node<Value>*> m_nodes;
    };
        
    // Assign node
    //--------------------------------------------------------------------------
    template<class Value>
    class AssignNode : public Node<Value>
    {
    public:
        explicit AssignNode(Expression<Value> *expr);
        ~AssignNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Value *m_var;
        Node<Value> *m_rhs;
    };
        
    // Add node
    //--------------------------------------------------------------------------
    template<class Value>
    class AddNode : public Node<Value>
    {
    public:
        explicit AddNode(Expression<Value> *expr);
        ~AddNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Node<Value> *m_lhs;
        Node<Value> *m_rhs;
    };
        
    // Subtract node
    //--------------------------------------------------------------------------
    template<class Value>
    class SubtractNode : public Node<Value>
    {
    public:
        explicit SubtractNode(Expression<Value> *expr);
        ~SubtractNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Node<Value> *m_lhs;
        Node<Value> *m_rhs;
    };
        
    // Multiply node
    //--------------------------------------------------------------------------
    template<class Value>
    class MultiplyNode : public Node<Value>
    {
    public:
        explicit MultiplyNode(Expression<Value> *expr);
        ~MultiplyNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Node<Value> *m_lhs;
        Node<Value> *m_rhs;
    };
        
    // Divide node
    //--------------------------------------------------------------------------
    template<class Value>
    class DivideNode : public Node<Value>
    {
    public:
        explicit DivideNode(Expression<Value> *expr);
        ~DivideNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Node<Value> *m_lhs;
        Node<Value> *m_rhs;
    };
        
    // Negate node
    //--------------------------------------------------------------------------
    template<class Value>
    class NegateNode : public Node<Value>
    {
    public:
        explicit NegateNode(Expression<Value> *expr);
        ~NegateNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Node<Value> *m_rhs;
    };
        
    // Exponent node
    //--------------------------------------------------------------------------
    template<class Value>
    class ExponentNode : public Node<Value>
    {
    public:
        explicit ExponentNode(Expression<Value> *expr);
        ~ExponentNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Node<Value> *m_lhs;
        Node<Value> *m_rhs;
    };
        
    // Variable node (also used for constants)
    //--------------------------------------------------------------------------
    template<class Value>
    class VariableNode : public Node<Value>
    {
    public:
        explicit VariableNode(Expression<Value> *expr);
        ~VariableNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Value *m_var;
    };
        
    // Value node
    //--------------------------------------------------------------------------
    template<class Value>
    class ValueNode : public Node<Value>
    {
    public:
        explicit ValueNode(Expression<Value> *expr);
        ~ValueNode();
        
        Value DoEvaluate() override;
        void Parse(Parser<Value> &parser,
                   typename Parser<Value>::size_type start,
                   typename Parser<Value>::size_type end,
                   typename Parser<Value>::size_type v1 = 0) override;
                
    private:
        Value m_val;
    };        
        
} // namespace ExprEval
    
#endif // __EXPREVAL_NODE_H  

