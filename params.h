#ifndef __PARAMS_H_
#define __PARAMS_H_

#include <QMetaType>

struct Params
{
    struct Consts {
        double K0f, K1f, Ur0, K0u, K1u;
        double T, Ty, Tu, Tg, Te, Tf, Tphi, Td0;
        double Pt0, omega_nom, Xdprime, Xd, Eqenom, Uc;
    } reg;
    struct Repl {
        double Y11, Y12, A11, A12, Y11em, Y12em, A11em, A12em, Pd;
    } repl;
    struct Start {
        double Delta0, Eqe0, Eqprime0, U0, V0;
    } start;
    double Tstart, Tstop, dt;
    double eps;
    size_t max_iterations;
    bool dirty_hack;
};

struct AnswerItem
{
    size_t row, set_no, n_steps;
    double time, delta, omega, Eqe, Eqprime, V, U;
};

Q_DECLARE_METATYPE(AnswerItem);

#endif // __PARAMS_H_
