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
    class Expression;
    class FunctionFactory;

    
    // Node class
    //--------------------------------------------------------------------------
    class Node
    {
    public:
        explicit Node(Expression *expr);
        virtual ~Node();
        
        virtual double DoEvaluate() = 0;
        virtual void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) = 0;
                
        double Evaluate(); // Calls Expression::TestAbort, then DoEvaluate
        
    protected:
        Expression *m_expr;    
    };
        
    // General function node class
    //--------------------------------------------------------------------------
    class FunctionNode : public Node
    {
    public:
        explicit FunctionNode(Expression *expr);
        ~FunctionNode();
        
        // Parse nodes and references
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
                
    private:
        // Function factory
        FunctionFactory *m_factory;
        
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
        ::std::vector<Node*> m_nodes;
        ::std::vector<double*> m_refs;

    friend class FunctionFactory;
    };
        
    // Mulit-expression node
    //--------------------------------------------------------------------------
    class MultiNode : public Node
    {
    public:
        explicit MultiNode(Expression *expr);
        ~MultiNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        ::std::vector<Node*> m_nodes;
    };
        
    // Assign node
    //--------------------------------------------------------------------------
    class AssignNode : public Node
    {
    public:
        explicit AssignNode(Expression *expr);
        ~AssignNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        double *m_var;
        Node *m_rhs;
    };
        
    // Add node
    //--------------------------------------------------------------------------
    class AddNode : public Node
    {
    public:
        explicit AddNode(Expression *expr);
        ~AddNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        Node *m_lhs;
        Node *m_rhs;
    };
        
    // Subtract node
    //--------------------------------------------------------------------------
    class SubtractNode : public Node
    {
    public:
        explicit SubtractNode(Expression *expr);
        ~SubtractNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        Node *m_lhs;
        Node *m_rhs;
    };
        
    // Multiply node
    //--------------------------------------------------------------------------
    class MultiplyNode : public Node
    {
    public:
        explicit MultiplyNode(Expression *expr);
        ~MultiplyNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        Node *m_lhs;
        Node *m_rhs;
    };
        
    // Divide node
    //--------------------------------------------------------------------------
    class DivideNode : public Node
    {
    public:
        explicit DivideNode(Expression *expr);
        ~DivideNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        Node *m_lhs;
        Node *m_rhs;
    };
        
    // Negate node
    //--------------------------------------------------------------------------
    class NegateNode : public Node
    {
    public:
        explicit NegateNode(Expression *expr);
        ~NegateNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        Node *m_rhs;
    };
        
    // Exponent node
    //--------------------------------------------------------------------------
    class ExponentNode : public Node
    {
    public:
        explicit ExponentNode(Expression *expr);
        ~ExponentNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        Node *m_lhs;
        Node *m_rhs;
    };
        
    // Variable node (also used for constants)
    //--------------------------------------------------------------------------
    class VariableNode : public Node
    {
    public:
        explicit VariableNode(Expression *expr);
        ~VariableNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        double *m_var;
    };
        
    // Value node
    //--------------------------------------------------------------------------
    class ValueNode : public Node
    {
    public:
        explicit ValueNode(Expression *expr);
        ~ValueNode();
        
        double DoEvaluate() override;
        void Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
                Parser::size_type v1 = 0) override;
                
    private:
        double m_val;
    };        
        
} // namespace ExprEval
    
#endif // __EXPREVAL_NODE_H  

