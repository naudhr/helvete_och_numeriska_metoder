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
    void newton_does_not_converge(QString name, double t, size_t n_steps);
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

#endif // __CALCULUS_H_
