// File:    node.cpp
// Author:  Brian Vanderburg II
// Purpose: Expression node
//------------------------------------------------------------------------------


// Includes
#include <new>
#include <memory>
#include <cmath>
#include <cerrno>

#include "defs.h"
#include "node.h"
#include "expr.h"
#include "vallist.h"
#include "funclist.h"
#include "except.h"

#include <autodiff/reverse/var.hpp>

using namespace std;
using namespace ExprEval;

// Node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
Node<Value>::Node(Expression<Value>* expr) : m_expr(expr)
{
    if(expr == 0)
        throw(NullPointerException("Node::Node"));
}

// Destructor
template<class Value>
Node<Value>::~Node()
{
}
    
// Evaluate
template<class Value>
Value Node<Value>::Evaluate()
{
    m_expr->TestAbort();
    
    return DoEvaluate();
}
    
// Function node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
FunctionNode<Value>::FunctionNode(Expression<Value> *expr) : Node<Value>(expr),
        m_argMin(0), m_argMax(0), m_refMin(0), m_refMax(0)
{
}
    
// Destructor
template<class Value>
FunctionNode<Value>::~FunctionNode()
{
    for (auto& node : m_nodes)
        delete node;
}
    
// Get name
template<class Value>
string FunctionNode<Value>::GetName() const
{
    return m_factory->GetName();
}
    
// Set argument count
template<class Value>
void FunctionNode<Value>::SetArgumentCount(long argMin, long argMax, long refMin, long refMax)
{
    m_argMin = argMin;
    m_argMax = argMax;
    m_refMin = refMin;
    m_refMax = refMax;
}
    
// Parse expression
template<class Value>
void FunctionNode<Value>::Parse(Parser<Value> &parser,
                                typename Parser<Value>::size_type start,
                                typename Parser<Value>::size_type end,
                                typename Parser<Value>::size_type v1)
{
    typename Parser<Value>::size_type pos, last;
    int plevel = 0;
    
    // Sanity check (start/end are function parenthesis)
    if(start >= end)
        throw(SyntaxException());
        
    // Look
    last = start + 1;
    for(pos = start + 1; pos <= end && start + 1 != end; pos++)
    {
        switch(parser[pos].GetType())
        {
            case Token::TypeOpenParenthesis:
            {
                plevel++;
                break;
            };
                
            case Token::TypeComma:
            case Token::TypeCloseParenthesis:            
            {
                // Handle Close parenthesis for all but the last one at the end
                if(parser[pos].GetType() == Token::TypeCloseParenthesis && pos != end)
                {
                    plevel--;
                
                    if(plevel < 0)
                    {
                        UnmatchedParenthesisException e;
                
                        e.SetStart(parser[pos].GetStart());
                        e.SetEnd(parser[pos].GetEnd());
            
                        throw(e);
                    }
                    
                    break;                    
                }
                    
                // Handle comma, or if it was the ending parenthesis treat it like comma
                if(plevel == 0)
                {
                    if(pos > last)
                    {
                        // reference parameter?
                        if(parser[last].GetType() == Token::TypeAmpersand)
                        {
                            // Reference parameter, check position and type of next parameter
                            if(last == pos - 2 && parser[last + 1].GetType() == Token::TypeIdentifier)
                            {
                                // Get value list
                                ValueList<Value> *vlist = this->m_expr->GetValueList();
                                if(vlist == 0)
                                {
                                    NoValueListException e;
                                    
                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());
                                    throw(e);
                                }
                                    
                                // Get name
                                string ident = parser[last + 1].GetIdentifier();
                                
                                // Make sure it is not a constant
                                if(vlist->IsConstant(ident))
                                {
                                    ConstantReferenceException e(ident);
                                    
                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());
                                    throw(e);
                                }
                                    
                                // Get address
                                Value *vaddr = vlist->GetAddress(ident);
                                if(vaddr == 0)
                                {
                                    // Try to add it and get again
                                    vlist->Add(ident);
                                    vaddr = vlist->GetAddress(ident);
                                }
                                    
                                if(vaddr == 0)
                                {
                                    NotFoundException e(ident);
                                    
                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());
                                    
                                    throw(e);
                                }
                                    
                                // Add it
                                m_refs.push_back(vaddr);
                            }
                            else
                            {
                                SyntaxException e;
                                
                                e.SetStart(parser[last].GetStart());
                                e.SetEnd(parser[pos].GetEnd());
                                throw(e);
                            }
                        } // TypeAmpersand
                        else
                        {
                            // Create node
                            aptr(Node<Value>) n(parser.ParseRegion(last, pos - 1));
                            m_nodes.push_back(n.get());
                            n.release();
                        }
                    }
                    else
                    {
                        SyntaxException e;
                        
                        e.SetStart(parser[pos].GetStart());
                        e.SetEnd(parser[pos].GetEnd());
                        throw(e);
                    }
                    
                    last = pos + 1;
                }
                    
                break;
            }

            default:
                break;
        }
    }
        
    // plevel should be zero
    if(plevel != 0)
    {
        UnmatchedParenthesisException e;
        
        e.SetStart(parser[end].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);        
    }


    // Check argument count
    if(m_argMin != -1 && m_nodes.size() < static_cast<size_t>(m_argMin))
    {
        InvalidArgumentCountException e(GetName());
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }

    if(m_argMax != -1 && m_nodes.size() > static_cast<size_t>(m_argMax))
    {
        InvalidArgumentCountException e(GetName());
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }

    if(m_refMin != -1 && m_refs.size() < static_cast<size_t>(m_refMin))
    {
        InvalidArgumentCountException e(GetName());
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }

    if(m_refMax != -1 && m_refs.size() > static_cast<size_t>(m_refMax))
    {
        InvalidArgumentCountException e(GetName());
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
}    

// Multi node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
MultiNode<Value>::MultiNode(Expression<Value> *expr) : Node<Value>(expr)
{
}
    
// Destructor
template<class Value>
MultiNode<Value>::~MultiNode()
{
    for(auto& node : m_nodes)
    {
        delete node;
    }
}
        
// Evaluate
template<class Value>
Value MultiNode<Value>::DoEvaluate()
{
    Value result = 0.0;
    
    for(auto& node : m_nodes)
    {
        result = node->Evaluate();
    }
        
    return result;
}
    
// Parse
template<class Value>
void MultiNode<Value>::Parse(Parser<Value> &parser,
                             typename Parser<Value>::size_type start,
                             typename Parser<Value>::size_type end,
                             typename Parser<Value>::size_type v1)
{
    typename Parser<Value>::size_type pos, last;
    int plevel = 0;
    
    // Sanity check
    if(start >= end)
        throw(SyntaxException());
        
    // Look
    last = start;
    for(pos = start; pos <= end; pos++)
    {
        switch(parser[pos].GetType())
        {
            case Token::TypeOpenParenthesis:
            {
                plevel++;
                break;
            };
                
            case Token::TypeCloseParenthesis:
            {
                plevel--;
                
                if(plevel < 0)
                {
                    UnmatchedParenthesisException e;
                    
                    e.SetStart(parser[pos].GetStart());
                    e.SetEnd(parser[pos].GetEnd());
                    throw(e);
                }
                    
                break;
            }
                
            case Token::TypeSemicolon:
            {
                if(plevel == 0)
                {
                    if(pos > last)
                    {
                        // Everything from last to pos - 1
                        aptr(Node<Value>) n(parser.ParseRegion(last, pos - 1));
                        m_nodes.push_back(n.get());
                        n.release();
                    }
                    else
                    {
                        SyntaxException e;
                        
                        e.SetStart(parser[last].GetStart());
                        e.SetEnd(parser[pos].GetEnd());
                        throw(e);
                    }
                    
                    last = pos + 1;
                }
                    
                break;
            }
          default:
                break;
        }
    }
        
    // plevel should be zero
    if(plevel != 0)
    {
        UnmatchedParenthesisException e;
        
        e.SetStart(parser[pos].GetStart());
        e.SetEnd(parser[pos].GetEnd());
        throw(e);
    }
        
    // If the end was not a semicolon, test it as well
    if(last < end + 1)
    {
        aptr(Node<Value>) n(parser.ParseRegion(last, end));
        m_nodes.push_back(n.get());
        n.release();
    }        
}
    
// Assign node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
AssignNode<Value>::AssignNode(Expression<Value> *expr) : Node<Value>(expr), m_var(0), m_rhs(0)
{
}
    
// Destructor
template<class Value>
AssignNode<Value>::~AssignNode()
{
    // Free child node
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value AssignNode<Value>::DoEvaluate()
{
    return (*m_var = m_rhs->Evaluate());        
}
    
// Parse
template<class Value>
void AssignNode<Value>::Parse(Parser<Value> &parser,
                              typename Parser<Value>::size_type start,
                              typename Parser<Value>::size_type end,
                              typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(v1 != start + 1 || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Get the value list
    ValueList<Value> *vlist = this->m_expr->GetValueList();
    if(vlist == 0)
    {
        NoValueListException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
    }
        
    // Determine variable name
    string ident = parser[start].GetIdentifier();
    
    // Make sure it is not a constant
    if(vlist->IsConstant(ident))
    {
        ConstantAssignException e(ident);
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[v1].GetEnd());
        throw(e);
    }
        
    // Get address
    Value *vaddr = vlist->GetAddress(ident);
    
    if(vaddr == 0)
    {
        // If it does not already exist, try to create it
        vlist->Add(ident);
        vaddr = vlist->GetAddress(ident);
    }
        
    if(vaddr == 0)
    {
        NotFoundException e(ident);
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
    }
        
    // Parse the node (will throw if it can not parse)
    aptr(Node<Value>) n(parser.ParseRegion(v1 + 1, end));
    
    // Set data
    m_var = vaddr;
    m_rhs = n.release();
}    
    

// Add node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
AddNode<Value>::AddNode(Expression<Value> *expr) : Node<Value>(expr), m_lhs(0), m_rhs(0)
{
}
    
// Destructor
template<class Value>
AddNode<Value>::~AddNode()
{
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value AddNode<Value>::DoEvaluate()
{
    return m_lhs->Evaluate() + m_rhs->Evaluate();        
}
    
// Parse
template<class Value>
void AddNode<Value>::Parse(Parser<Value> &parser,
                           typename Parser<Value>::size_type start,
                           typename Parser<Value>::size_type end,
                           typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Parse sides
    aptr(Node<Value>) left(parser.ParseRegion(start, v1 - 1));
    aptr(Node<Value>) right(parser.ParseRegion(v1 + 1, end));
    
    m_lhs = left.release();
    m_rhs = right.release();
}
    
// Subtract node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
SubtractNode<Value>::SubtractNode(Expression<Value> *expr) : Node<Value>(expr), m_lhs(0), m_rhs(0)
{
}
    
// Destructor
template<class Value>
SubtractNode<Value>::~SubtractNode()
{
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value SubtractNode<Value>::DoEvaluate()
{
    return m_lhs->Evaluate() - m_rhs->Evaluate();        
}
    
// Parse
template<class Value>
void SubtractNode<Value>::Parse(Parser<Value> &parser,
                                typename Parser<Value>::size_type start,
                                typename Parser<Value>::size_type end,
                                typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Parse sides
    aptr(Node<Value>) left(parser.ParseRegion(start, v1 - 1));
    aptr(Node<Value>) right(parser.ParseRegion(v1 + 1, end));
    
    m_lhs = left.release();
    m_rhs = right.release();
}    
    
// Multiply node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
MultiplyNode<Value>::MultiplyNode(Expression<Value> *expr) : Node<Value>(expr), m_lhs(0), m_rhs(0)
{
}
    
// Destructor
template<class Value>
MultiplyNode<Value>::~MultiplyNode()
{
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value MultiplyNode<Value>::DoEvaluate()
{
    return m_lhs->Evaluate() * m_rhs->Evaluate();        
}
    
// Parse
template<class Value>
void MultiplyNode<Value>::Parse(Parser<Value> &parser,
                                typename Parser<Value>::size_type start,
                                typename Parser<Value>::size_type end,
                                typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Parse sides
    aptr(Node<Value>) left(parser.ParseRegion(start, v1 - 1));
    aptr(Node<Value>) right(parser.ParseRegion(v1 + 1, end));
    
    m_lhs = left.release();
    m_rhs = right.release();
}    
    
// Divide node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
DivideNode<Value>::DivideNode(Expression<Value> *expr) : Node<Value>(expr), m_lhs(0), m_rhs(0)
{
}
    
// Destructor
template<class Value>
DivideNode<Value>::~DivideNode()
{
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value DivideNode<Value>::DoEvaluate()
{
    Value r2 = m_rhs->Evaluate();
    
    if(r2 != Value(0.0))
    {
        return m_lhs->Evaluate() / r2;
    }
    else
    {
        throw(DivideByZeroException());
    }
}
    
// Parse
template<class Value>
void DivideNode<Value>::Parse(Parser<Value> &parser,
                              typename Parser<Value>::size_type start,
                              typename Parser<Value>::size_type end,
                              typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Parse sides
    aptr(Node<Value>) left(parser.ParseRegion(start, v1 - 1));
    aptr(Node<Value>) right(parser.ParseRegion(v1 + 1, end));
    
    m_lhs = left.release();
    m_rhs = right.release();
}    
    
// Negate node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
NegateNode<Value>::NegateNode(Expression<Value> *expr) : Node<Value>(expr), m_rhs(0)
{
}
    
// Destructor
template<class Value>
NegateNode<Value>::~NegateNode()
{
    // Free child nodes
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value NegateNode<Value>::DoEvaluate()
{
    return -(m_rhs->Evaluate());        
}
    
// Parse
template<class Value>
void NegateNode<Value>::Parse(Parser<Value> &parser,
                              typename Parser<Value>::size_type start,
                              typename Parser<Value>::size_type end,
                              typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(start != v1 || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Parse sides
    aptr(Node<Value>) right(parser.ParseRegion(v1 + 1, end));
    
    m_rhs = right.release();
}
    
// Exponent node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
ExponentNode<Value>::ExponentNode(Expression<Value> *expr) : Node<Value>(expr), m_lhs(0), m_rhs(0)
{
}
    
// Destructor
template<class Value>
ExponentNode<Value>::~ExponentNode()
{
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
}
        
// Evaluate
template<class Value>
Value ExponentNode<Value>::DoEvaluate()
{
    errno = 0;

    Value result = pow(m_lhs->Evaluate(), m_rhs->Evaluate());

    if constexpr (std::is_same_v<Value, double>) {
        if(errno)
            throw(MathException("^"));
    }
        
    return result;        
}
    
// Parse
template<class Value>
void ExponentNode<Value>::Parse(Parser<Value> &parser,
                                typename Parser<Value>::size_type start,
                                typename Parser<Value>::size_type end,
                                typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Parse sides
    aptr(Node<Value>) left(parser.ParseRegion(start, v1 - 1));
    aptr(Node<Value>) right(parser.ParseRegion(v1 + 1, end));
    
    m_lhs = left.release();
    m_rhs = right.release();
}    
    
// Variable node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
VariableNode<Value>::VariableNode(Expression<Value> *expr) : Node<Value>(expr), m_var(0)
{
}
    
// Destructor
template<class Value>
VariableNode<Value>::~VariableNode()
{
}
        
// Evaluate
template<class Value>
Value VariableNode<Value>::DoEvaluate()
{
    return *m_var;        
}
    
// Parse
template<class Value>
void VariableNode<Value>::Parse(Parser<Value> &parser,
                                typename Parser<Value>::size_type start,
                                typename Parser<Value>::size_type end,
                                typename Parser<Value>::size_type v1)
{
    // Check some basic syntax
    if(start != end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Get value list
    ValueList<Value> *vlist = this->m_expr->GetValueList();
    if(vlist == 0)
    {
        NoValueListException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
    }
    
    // Get name
    string ident = parser[start].GetIdentifier();
    
    // Get address
    Value *vaddr = vlist->GetAddress(ident);
    
    if(vaddr == 0)
    {
        // If it does not already exist, try to create it
        vlist->Add(ident);
        vaddr = vlist->GetAddress(ident);
    }
        
    if(vaddr == 0)
    {
        NotFoundException e(ident);
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
    }
        
    // Set information
    m_var = vaddr;
}            
    
// Value node
//------------------------------------------------------------------------------

// Constructor
template<class Value>
ValueNode<Value>::ValueNode(Expression<Value> *expr) : Node<Value>(expr), m_val(0)
{
}
    
// Destructor
template<class Value>
ValueNode<Value>::~ValueNode()
{
}
        
// Evaluate
template<class Value>
Value ValueNode<Value>::DoEvaluate()
{
    return m_val;        
}
    
// Parse
template<class Value>
void ValueNode<Value>::Parse(Parser<Value> &parser,
                             typename Parser<Value>::size_type start,
                             typename Parser<Value>::size_type end,
                             typename Parser<Value>::size_type v1)
{
    // Check basic syntax
    if(start != end)
    {
        SyntaxException e;
        
        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
    }
        
    // Set info
    m_val = parser[start].GetValue();
}

#define INSTANCE(T) \
  template class T<double>; \
  template class T<autodiff::var>;

namespace ExprEval {
INSTANCE(AddNode)
INSTANCE(AssignNode)
INSTANCE(DivideNode)
INSTANCE(ExponentNode)
INSTANCE(FunctionNode)
INSTANCE(MultiNode)
INSTANCE(MultiplyNode)
INSTANCE(NegateNode)
INSTANCE(Node)
INSTANCE(SubtractNode)
INSTANCE(ValueNode)
INSTANCE(VariableNode)
}
