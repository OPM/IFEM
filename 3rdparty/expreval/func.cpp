// File:    func.cpp
// Author:  Brian Allen Vanderburg II
// Purpose: ExprEval internal functions
//------------------------------------------------------------------------------

// Includes
#include <new>
#include <memory>

#include <cmath>
#include <cerrno>
#include <ctime>

#include "autodiff/reverse/var/var.hpp"
#include "defs.h"
#include "funclist.h"
#include "node.h"
#include "except.h"

using namespace ExprEval;

// Anonymous namespace for items
namespace
{
    // Absolute value
    //--------------------------------------------------------------------------
    template<class Value>
    class abs_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit abs_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            if constexpr (std::is_same_v<Value, Real>)
                return fabs(this->m_nodes[0]->Evaluate());
            else
                return abs(this->m_nodes[0]->Evaluate());
        }
    };

    template<class Value>
    class abs_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "abs";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new abs_FunctionNode<Value>(expr);
        }
    };

    // Modulus
    //--------------------------------------------------------------------------
    template<class Value>
    class mod_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit mod_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = fmod(this->m_nodes[0]->Evaluate(),
                                this->m_nodes[1]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class mod_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "mod";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new mod_FunctionNode<Value>(expr);
        }
    };

    // Integer part
    //--------------------------------------------------------------------------
    template<class Value>
    class ipart_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit ipart_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value result;

            modf(this->m_nodes[0]->Evaluate(), &result);

            return result;
        }
    };

    template<class Value>
    class ipart_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "ipart";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new ipart_FunctionNode<Value>(expr);
        }
    };

    // Fraction part
    //--------------------------------------------------------------------------
    template<class Value>
    class fpart_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit fpart_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value dummy;

            return modf(this->m_nodes[0]->Evaluate(), &dummy);
        }
    };

    template<class Value>
    class fpart_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "fpart";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new fpart_FunctionNode<Value>(expr);
        }
    };

    // Minimum
    //--------------------------------------------------------------------------
    template<class Value>
    class min_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit min_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, -1, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value result = 0;
            size_t pos = 0;
            for (const auto& node : this->m_nodes) {
                Value tmp = node->Evaluate();
                if (pos == 0 || tmp < result)
                  result = tmp;
                ++pos;
            }

            return result;
        }
    };

    template<class Value>
    class min_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "min";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new min_FunctionNode<Value>(expr);
        }
    };

    // Maximum
    //--------------------------------------------------------------------------
    template<class Value>
    class max_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit max_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, -1, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value result = 0;
            size_t pos = 0;
            for (const auto& node : this->m_nodes)
            {
                Value tmp = node->Evaluate();
                if (pos == 0 || tmp > result)
                  result = tmp;
                ++pos;
            }

            return result;
        }
    };

    template<class Value>
    class max_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "max";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new max_FunctionNode<Value>(expr);
        }
    };

    // Square root
    //--------------------------------------------------------------------------
    template<class Value>
    class sqrt_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit sqrt_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = sqrt(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class sqrt_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "sqrt";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new sqrt_FunctionNode<Value>(expr);
        }
    };

    // Sine
    //--------------------------------------------------------------------------
    template<class Value>
    class sin_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit sin_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = sin(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class sin_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "sin";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new sin_FunctionNode<Value>(expr);
        }
    };

    // Cosine
    //--------------------------------------------------------------------------
    template<class Value>
    class cos_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit cos_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = cos(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class cos_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "cos";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new cos_FunctionNode<Value>(expr);
        }
    };

    // Tangent
    //--------------------------------------------------------------------------
    template<class Value>
    class tan_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit tan_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = tan(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class tan_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "tan";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new tan_FunctionNode<Value>(expr);
        }
    };

    // Hyperbolic Sine
    //--------------------------------------------------------------------------
    template<class Value>
    class sinh_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit sinh_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            auto result = sinh(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class sinh_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "sinh";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new sinh_FunctionNode<Value>(expr);
        }
    };

    // Hyperbolic Cosine
    //--------------------------------------------------------------------------
    template<class Value>
    class cosh_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit cosh_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = cosh(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class cosh_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "cosh";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new cosh_FunctionNode<Value>(expr);
        }
    };

    // Hyperbolic Tangent
    //--------------------------------------------------------------------------
    template<class Value>
    class tanh_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit tanh_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = tanh(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class tanh_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "tanh";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new tanh_FunctionNode<Value>(expr);
        }
    };

    // Arc Sine
    //--------------------------------------------------------------------------
    template<class Value>
    class asin_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit asin_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = asin(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class asin_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "asin";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new asin_FunctionNode<Value>(expr);
        }
    };

    // Arc Cosine
    //--------------------------------------------------------------------------
    template<class Value>
    class acos_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit acos_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = acos(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class acos_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "acos";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new acos_FunctionNode<Value>(expr);
        }
    };

    // Arc Tangent
    //--------------------------------------------------------------------------
    template<class Value>
    class atan_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit atan_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = atan(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class atan_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "atan";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new atan_FunctionNode<Value>(expr);
        }
    };

    // Arc Tangent 2: atan2(y, x)
    //--------------------------------------------------------------------------
    template<class Value>
    class atan2_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit atan2_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = atan2(this->m_nodes[0]->Evaluate(), this->m_nodes[1]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class atan2_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "atan2";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new atan2_FunctionNode<Value>(expr);
        }
    };

    // Log
    //--------------------------------------------------------------------------
    template<class Value>
    class log_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit log_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = log10(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class log_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "log";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new log_FunctionNode<Value>(expr);
        }
    };

    // Ln
    //--------------------------------------------------------------------------
    template<class Value>
    class ln_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit ln_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = log(this->m_nodes[0]->Evaluate());

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class ln_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "ln";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new ln_FunctionNode<Value>(expr);
        }
    };

    // Exp
    //--------------------------------------------------------------------------
    template<class Value>
    class exp_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit exp_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;
            Value x = this->m_nodes[0]->Evaluate();
            Value result = exp(x);

            if (errno && x > Value(0.0)) // Fixed: No exception on underflow
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class exp_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "exp";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new exp_FunctionNode<Value>(expr);
        }
    };

    // Logn
    //--------------------------------------------------------------------------
    template<class Value>
    class logn_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit logn_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            // Check for division by zero
            Value tmp = log(this->m_nodes[1]->Evaluate());

            if(tmp == 0.0)
                throw(MathException(this->GetName()));

            // Calculate result
            Value result = log(this->m_nodes[0]->Evaluate()) / tmp;

            if(errno)
                throw(MathException(this->GetName()));

            return result;
        }
    };

    template<class Value>
    class logn_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "logn";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new logn_FunctionNode<Value>(expr);
        }
    };

    template<class Value>
    class pow_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit pow_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value result = pow(this->m_nodes[0]->Evaluate(), this->m_nodes[1]->Evaluate());

            if constexpr (std::is_same_v<Value, double>) { // disable exception with autodiff
                if(errno)
                    throw(MathException(this->GetName()));
            }

            return result;
        }
    };

    template<class Value>
    class pow_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "pow";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new pow_FunctionNode<Value>(expr);
        }
    };

    // Ceil
    //--------------------------------------------------------------------------
    template<class Value>
    class ceil_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit ceil_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            return ceil(this->m_nodes[0]->Evaluate());
        }
    };

    template<class Value>
    class ceil_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "ceil";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new ceil_FunctionNode<Value>(expr);
        }
    };

    // Floor
    //--------------------------------------------------------------------------
    template<class Value>
    class floor_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit floor_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            return floor(this->m_nodes[0]->Evaluate());
        }
    };

    template<class Value>
    class floor_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "floor";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new floor_FunctionNode<Value>(expr);
        }
    };

    // Rand
    //--------------------------------------------------------------------------

    // Get next random (0,1)
    inline double NextRandom(double *seed)
    {
        long a = (long)(*seed) * 214013L + 2531011L;
        *seed = (double)a;

        return (double)((a >> 16) & 0x7FFF) / 32767.0;
    }

    template<class Value>
    class rand_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit rand_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(0, 0, 1, 1);
        }

        Value DoEvaluate() override
        {
            return NextRandom(this->m_refs[0]);
        }
    };

    template<class Value>
    class rand_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "rand";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new rand_FunctionNode<Value>(expr);
        }
    };

    // Random
    //--------------------------------------------------------------------------
    template<class Value>
    class random_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit random_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 1, 1);
        }

        Value DoEvaluate() override
        {
            Value a = this->m_nodes[0]->Evaluate();
            Value b = this->m_nodes[1]->Evaluate();

            return NextRandom(this->m_refs[0]) * (b - a) + a;
        }
    };

    template<class Value>
    class random_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "random";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new random_FunctionNode<Value>(expr);
        }
    };

    // Randomize
    //--------------------------------------------------------------------------
    template<class Value>
    class randomize_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit randomize_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(0, 0, 1, 1);
        }

        Value DoEvaluate() override
        {
            static long curcall = 1;

            *this->m_refs[0] = (Value)(((clock() + 1024u + curcall) * time(NULL)) % 2176971487u);
            curcall++;

            return 0.0;
        }
    };

    template<class Value>
    class randomize_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "randomize";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new randomize_FunctionNode<Value>(expr);
        }
    };

    // Radians to degrees
    //--------------------------------------------------------------------------
    template<class Value>
    class deg_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit deg_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            return (this->m_nodes[0]->Evaluate() * 180.0) / M_PI;
        }
    };

    template<class Value>
    class deg_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "deg";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new deg_FunctionNode<Value>(expr);
        }
    };

    // Degrees to radians
    //--------------------------------------------------------------------------
    template<class Value>
    class rad_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit rad_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            return (this->m_nodes[0]->Evaluate() * M_PI) / 180.0;
        }
    };

    template<class Value>
    class rad_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "rad";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new rad_FunctionNode<Value>(expr);
        }
    };

    // Rectangular to polar: rect2pol(x, y, &distance, &angle)
    //--------------------------------------------------------------------------
    template<class Value>
    class rect2pol_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit rect2pol_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 2, 2);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value x = this->m_nodes[0]->Evaluate();
            Value y = this->m_nodes[1]->Evaluate();

            Value d = hypot(x, y);
            Value a = atan2(y, x);

            if(errno)
                throw(MathException(this->GetName()));

            *this->m_refs[0] = d;
            if(a < 0.0)
                *this->m_refs[1] = a + (2.0 * M_PI);
            else
                *this->m_refs[1] = a;

            return d;
        }
    };

    template<class Value>
    class rect2pol_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "rect2pol";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new rect2pol_FunctionNode<Value>(expr);
        }
    };

    // Polar to rectangular: pol2rect(distance, angle, &x, &y)
    //--------------------------------------------------------------------------
    template<class Value>
    class pol2rect_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit pol2rect_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 2, 2);
        }

        Value DoEvaluate() override
        {
            errno = 0;

            Value d = this->m_nodes[0]->Evaluate();
            Value a = this->m_nodes[1]->Evaluate();

            Value x = d * cos(a);
            Value y = d * sin(a);

            if(errno)
                throw(MathException(this->GetName()));

            *this->m_refs[0] = x;
            *this->m_refs[1] = y;

            return x;
        }
    };

    template<class Value>
    class pol2rect_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "pol2rect";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new pol2rect_FunctionNode<Value>(expr);
        }
    };

    // If
    //--------------------------------------------------------------------------
    template<class Value>
    class if_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit if_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(3, 3, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value c = this->m_nodes[0]->Evaluate();

            if(c == 0.0)
                return this->m_nodes[2]->Evaluate();
            else
                return this->m_nodes[1]->Evaluate();
        }
    };

    template<class Value>
    class if_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "if";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new if_FunctionNode<Value>(expr);
        }
    };

    // Select
    //--------------------------------------------------------------------------
    template<class Value>
    class select_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit select_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(3, 4, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value c = this->m_nodes[0]->Evaluate();

            if(c < 0.0)
                return this->m_nodes[1]->Evaluate();
            else if(c == 0.0)
                return this->m_nodes[2]->Evaluate();
            else
            {
                if(this->m_nodes.size() == 3)
                    return this->m_nodes[2]->Evaluate();
                else
                    return this->m_nodes[3]->Evaluate();
            }
        }
    };

    template<class Value>
    class select_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "select";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new select_FunctionNode<Value>(expr);
        }
    };

    // Equal
    //--------------------------------------------------------------------------
    template<class Value>
    class equal_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit equal_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            if(this->m_nodes[0]->Evaluate() == this->m_nodes[1]->Evaluate())
                return 1.0;
            else
                return 0.0;
        }
    };

    template<class Value>
    class equal_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "equal";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new equal_FunctionNode<Value>(expr);
        }
    };

    // Above
    //--------------------------------------------------------------------------
    template<class Value>
    class above_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit above_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            if(this->m_nodes[0]->Evaluate() > this->m_nodes[1]->Evaluate())
                return 1.0;
            else
                return 0.0;
        }
    };

    template<class Value>
    class above_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "above";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new above_FunctionNode<Value>(expr);
        }
    };

    // Below
    //--------------------------------------------------------------------------
    template<class Value>
    class below_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit below_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            if(this->m_nodes[0]->Evaluate() < this->m_nodes[1]->Evaluate())
                return 1.0;
            else
                return 0.0;
        }
    };

    template<class Value>
    class below_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "below";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new below_FunctionNode<Value>(expr);
        }
    };

    // Clip
    //--------------------------------------------------------------------------
    template<class Value>
    class clip_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit clip_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(3, 3, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value v = this->m_nodes[0]->Evaluate();
            Value a = this->m_nodes[1]->Evaluate();
            Value b = this->m_nodes[2]->Evaluate();

            if(v < a)
                return a;
            else if(v > b)
                return b;
            else
                return v;
        }
    };

    template<class Value>
    class clip_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "clip";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new clip_FunctionNode<Value>(expr);
        }
    };

    // Clamp
    //--------------------------------------------------------------------------
    template<class Value>
    class clamp_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit clamp_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(3, 3, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value v = this->m_nodes[0]->Evaluate();
            Value a = this->m_nodes[1]->Evaluate();
            Value b = this->m_nodes[2]->Evaluate();

            if(a == b)
                return a;
            else
            {
                Value tmp = fmod(v - a, b - a);

                if(tmp < 0)
                    return tmp + b;
                else
                    return tmp + a;
            }
        }
    };

    template<class Value>
    class clamp_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "clamp";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new clamp_FunctionNode<Value>(expr);
        }
    };

    // Rescale
    //--------------------------------------------------------------------------
    template<class Value>
    class rescale_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit rescale_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(5, 5, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value pnt = this->m_nodes[0]->Evaluate();
            Value o1 = this->m_nodes[1]->Evaluate();
            Value o2 = this->m_nodes[2]->Evaluate();
            Value n1 = this->m_nodes[3]->Evaluate();
            Value n2 = this->m_nodes[4]->Evaluate();

            Value odiff = o2 - o1;
            if(odiff == 0.0)
                return n1;
            else
            {
                return (pnt - o1) * (n2 - n1) / odiff + n1;
            }
        }
    };

    template<class Value>
    class rescale_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "rescale";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new rescale_FunctionNode<Value>(expr);
        }
    };

    // Poly: poly(x, c3, c2, c1, c0): c3*x^3 + c2*x^2 + c1*x + c0
    //--------------------------------------------------------------------------
    template<class Value>
    class poly_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit poly_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, -1, 0, 0);
        }

        Value DoEvaluate() override
        {
            Value total = 0.0;
            double curpow;

            curpow = (double)(this->m_nodes.size() - 2);

            // Value of x
            Value x = this->m_nodes[0]->Evaluate();

            errno = 0;

            for (const auto& node : this->m_nodes)
            {
                total += node->Evaluate() * pow(x, curpow);
                curpow -= 1.0;
            }

            if(errno)
                throw(MathException(this->GetName()));

            return total;
        }
    };

    template<class Value>
    class poly_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "poly";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new poly_FunctionNode<Value>(expr);
        }
    };

    // And
    //--------------------------------------------------------------------------
    template<class Value>
    class and_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit and_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            if(this->m_nodes[0]->Evaluate() == 0.0 || this->m_nodes[1]->Evaluate() == 0.0)
                return 0.0;
            else
                return 1.0;
        }
    };

    template<class Value>
    class and_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "and";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new and_FunctionNode<Value>(expr);
        }
    };

    // Or
    //--------------------------------------------------------------------------
    template<class Value>
    class or_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit or_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(2, 2, 0, 0);
        }

        Value DoEvaluate() override
        {
            if(this->m_nodes[0]->Evaluate() == 0.0 && this->m_nodes[1]->Evaluate() == 0.0)
                return 0.0;
            else
                return 1.0;
        }
    };

    template<class Value>
    class or_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "or";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new or_FunctionNode<Value>(expr);
        }
    };

    // Not
    //--------------------------------------------------------------------------
    template<class Value>
    class not_FunctionNode : public FunctionNode<Value>
    {
    public:
        explicit not_FunctionNode(Expression<Value> *expr) : FunctionNode<Value>(expr)
        {
            this->SetArgumentCount(1, 1, 0, 0);
        }

        Value DoEvaluate() override
        {
            if(this->m_nodes[0]->Evaluate() == 0.0)
                return 1.0;
            else
                return 0.0;
        }
    };

    template<class Value>
    class not_FunctionFactory : public FunctionFactory<Value>
    {
    public:
        std::string GetName() const override
        {
            return "not";
        }

        FunctionNode<Value> *DoCreate(Expression<Value> *expr) override
        {
            return new not_FunctionNode<Value>(expr);
        }
    };


} // namespace


// Initialize default functions
template<class Value>
void FunctionList<Value>::AddDefaultFunctions()
{
    #define ADDFUNCTION(name) \
    aptr(FunctionFactory<Value>) name ## _func(new name ## _FunctionFactory<Value>()); \
        m_functions.push_back(name ## _func.get()); \
        name ## _func.release();

    ADDFUNCTION(abs);

    ADDFUNCTION(sqrt);

    ADDFUNCTION(sin);
    ADDFUNCTION(cos);
    ADDFUNCTION(tan);

    ADDFUNCTION(sinh);
    ADDFUNCTION(cosh);
    ADDFUNCTION(tanh);

    ADDFUNCTION(asin);
    ADDFUNCTION(acos);
    ADDFUNCTION(atan);
    ADDFUNCTION(atan2);

    ADDFUNCTION(log);
    ADDFUNCTION(ln);
    ADDFUNCTION(exp);
    ADDFUNCTION(logn);
    ADDFUNCTION(pow);

    ADDFUNCTION(deg);
    ADDFUNCTION(rad);

    ADDFUNCTION(rect2pol);
    ADDFUNCTION(pol2rect);

    ADDFUNCTION(if); // Preprocess will take care of this beforehand
    ADDFUNCTION(select);

    ADDFUNCTION(equal);
    ADDFUNCTION(above);
    ADDFUNCTION(below);

    ADDFUNCTION(clip);
    ADDFUNCTION(rescale);

    ADDFUNCTION(poly);

    ADDFUNCTION(and);
    ADDFUNCTION(or);
    ADDFUNCTION(not);

    if constexpr (std::is_same_v<Value,double>) {
      ADDFUNCTION(mod);

      ADDFUNCTION(ipart);
      ADDFUNCTION(fpart);

      ADDFUNCTION(min);
      ADDFUNCTION(max);

      ADDFUNCTION(ceil);
      ADDFUNCTION(floor);

      ADDFUNCTION(rand);
      ADDFUNCTION(random);
      ADDFUNCTION(randomize);

      ADDFUNCTION(clamp);
    }
}

namespace ExprEval {
template void FunctionList<double>::AddDefaultFunctions();
template class FunctionList<autodiff::var>;
}
