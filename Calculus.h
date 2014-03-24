#ifndef __CALCULUS_H_
#define __CALCULUS_H_

#include <QThread>
#include "params.h"

struct AnswerItem;
struct Params;

class Calculus : public QThread
{
    Q_OBJECT
  protected:
    Params p;
    virtual void run();
    virtual void emit_x(double t, size_t row, size_t n_steps) = 0;
    virtual void solve_newton(double t, size_t& n_steps) = 0;
    virtual const char* name() const {  return "";  }
  signals:
    void a_step_done(AnswerItem);
    void newton_does_not_converge(QString name, double t, unsigned n_steps);
  public:
    Calculus(const Params& _p) : p(_p) {}
    virtual ~Calculus() {}
};

class CalculusEiler : public Calculus
{
    Q_OBJECT
    struct Impl;
    Impl* pimpl;
    virtual void emit_x(double t, size_t row, size_t n_steps);
    virtual void solve_newton(double t, size_t& n_steps);
    virtual const char* name() const {  return "Eiler";  }
 public:
    CalculusEiler(const Params& _p);
    virtual ~CalculusEiler() {}

    static double calculate_Pg(double Eqprime, double U, double Xdprime, double Xd, double s_d_v, double c_d_v);
    static double calculate_Qg(double Eqprime, double U, double Xdprime, double Xd, double s_d_v, double c_d_v);

    static double calculate_Pc(double U, double V, double Y11, double Y12, double A11, double A12, double Uc);
    static double calculate_Qc(double U, double V, double Y11, double Y12, double A11, double A12, double Uc);
};

class CalculusTrapeze : public Calculus
{
    Q_OBJECT
    struct Impl;
    Impl* pimpl;
    virtual void emit_x(double t, size_t row, size_t n_steps);
    virtual void solve_newton(double t, size_t& n_steps);
    virtual const char* name() const {  return "Trapeze";  }
 public:
    CalculusTrapeze(const Params& _p);
    virtual ~CalculusTrapeze() {}
};

class CalculusSequensive : public Calculus
{
    Q_OBJECT
    struct Impl;
    Impl* pimpl;
    virtual void emit_x(double t, size_t row, size_t n_steps);
    virtual void solve_newton(double t, size_t& n_steps);
    virtual const char* name() const {  return "Sequensive";  }
 public:
    CalculusSequensive(const Params& _p);
    virtual ~CalculusSequensive() {}
    void set_X(double x1, double x2, double x4, double x5, double x6, double x7, double x9, double x0);
};

class CalculusParallel : public Calculus
{
    Q_OBJECT
    struct Impl;
    Impl* pimpl;
    virtual void emit_x(double t, size_t row, size_t n_steps);
    virtual void solve_newton(double t, size_t& n_steps);
    virtual const char* name() const {  return "Parallel";  }
 public:
    CalculusParallel(const Params& _p);
    virtual ~CalculusParallel() {}
};

struct Equiv
{
    double Y11, Y12, A11, A12;
    double Pd, Upphi;
    double Tsw_low, Tsw_high;
    Equiv() : Y11(0), Y12(0), A11(0), A12(0), Pd(0), Upphi(0), Tsw_low(-1), Tsw_high(-1) {}
};

Equiv recalculate_equiv_params(double t, const double U, const Params& p, Equiv e);

#endif // __CALCULUS_H_
