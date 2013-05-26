#ifndef __PARAMS_H_
#define __PARAMS_H_

struct Params
{
    struct Reg {
        double K0f, K1f, Ur0, K0u, K1u;
    } reg;
    struct Repl {
        double Y11, Y12, A11, A12, Y11em, Y12em, A11em, A12em, P0f;
    } repl;
    struct Start {
        double Delta0, Eqe0, Eqprime0, U0, V0;
    } start;
    double Tstart, Tstop, dt;
    double eps;
    size_t max_iterations;
};

struct AnswerItem
{
    double time, delta, omega, Eqe, Eqprime, V, U;
};

#endif // __PARAMS_H_
