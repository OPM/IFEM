// File:    except.h
// Author:  Brian Vanderburg II
// Purpose: ExprEval exceptions
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_EXCEPT_H
#define __EXPREVAL_EXCEPT_H

// Includes
#include <exception>
#include <string>

// Part of expreval namespace
namespace ExprEval
{
    // Forward declarations
    class Expression;

    // Exception class
    //--------------------------------------------------------------------------
    class Exception : public ::std::exception
    {
    public:
        // Types of exceptions (to simplify some error handling)
        // Each exception must set m_type in the construct
        // Each exception should have a different type
        enum Type
        {
            Type_Exception = 1,
            Type_NotFoundException,
            Type_AlreadyExistsException,
            Type_NullPointerException,
            Type_MathException,
            Type_DivideByZeroException,
            Type_NoValueListException,
            Type_NoFunctionListException,
            Type_AbortException,
            Type_EmptyExpressionException,
            Type_UnknownTokenException,
            Type_InvalidArgumentCountException,
            Type_ConstantAssignException,
            Type_ConstantReferenceException,
            Type_SyntaxException,
            Type_UnmatchedParenthesisException
        };

        Exception();
        virtual ~Exception() throw();

        Type GetType() const { return m_type; }
        const ::std::string &GetValue() const { return m_value; }
        const ::std::string &GetError() const { return m_error; }

        void SetStart(::std::string::size_type start) { m_start = start; }
        void SetEnd(::std::string::size_type end) { m_end = end; }

        ::std::string::size_type GetStart() const { return m_start; }
        ::std::string::size_type GetEnd() const { return m_end; }

    protected:
        Type m_type;
        ::std::string m_value, m_error;
        ::std::string::size_type m_start, m_end;
    };

    // Not found exception (for functions)
    //--------------------------------------------------------------------------
    class NotFoundException : public Exception
    {
    public:
        explicit NotFoundException(const ::std::string &name);
    };

    // Already exists exception (function or value already exists)
    //--------------------------------------------------------------------------
    class AlreadyExistsException : public Exception
    {
    public:
        explicit AlreadyExistsException(const ::std::string &name);
    };

    // A null pointer was passed
    //--------------------------------------------------------------------------
    class NullPointerException : public Exception
    {
    public:
        explicit NullPointerException(const ::std::string &method);
    };

    // A bad math error occured
    //--------------------------------------------------------------------------
    class MathException : public Exception
    {
    public:
        explicit MathException(const ::std::string &function);
    };

    // Divide by zero exception
    //--------------------------------------------------------------------------
    class DivideByZeroException : public Exception
    {
    public:
        DivideByZeroException();
    };

    // No value list exception (if expression uses variables or constants
    // but list does not exist) during parsing
    //--------------------------------------------------------------------------
    class NoValueListException : public Exception
    {
    public:
        NoValueListException();
    };

    // No function list exception (if expression uses functions) during parsing
    //--------------------------------------------------------------------------
    class NoFunctionListException : public Exception
    {
    public:
        NoFunctionListException();
    };

    // Abort exception (if DoTestAbort returns true)
    //--------------------------------------------------------------------------
    class AbortException : public Exception
    {
    public:
        AbortException();
    };

    // Empty expression (for parsing or evaluation)
    //--------------------------------------------------------------------------
    class EmptyExpressionException : public Exception
    {
    public:
        EmptyExpressionException();
    };

    // Unknown token in expression string
    //--------------------------------------------------------------------------
    class UnknownTokenException : public Exception
    {
    public:
        UnknownTokenException();
    };

    // Invalid argument count to function
    //--------------------------------------------------------------------------
    class InvalidArgumentCountException : public Exception
    {
    public:
        explicit InvalidArgumentCountException(const ::std::string &function);
    };

    // Assign to a constant
    //--------------------------------------------------------------------------
    class ConstantAssignException : public Exception
    {
    public:
        explicit ConstantAssignException(const ::std::string &value);
    };

    // Pass constant by reference
    //--------------------------------------------------------------------------
    class ConstantReferenceException : public Exception
    {
    public:
        explicit ConstantReferenceException(const ::std::string &value);
    };

    // A general syntax exception
    //--------------------------------------------------------------------------
    class SyntaxException : public Exception
    {
    public:
        SyntaxException();
    };

    // Unmatched parenthesis
    //--------------------------------------------------------------------------
    class UnmatchedParenthesisException : public Exception
    {
    public:
        UnmatchedParenthesisException();
    };

} // namespace ExprEval

#endif // __EXPREVAL_EXCEPT_H
