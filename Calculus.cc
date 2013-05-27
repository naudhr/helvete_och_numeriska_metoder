#include "Calculus.h"
#include "params.h"
#include <QDebug>
#include <cassert>
#include <cmath>

// the Trapeze and Eiler matrix and vector statements are fully and trully
// Chich-Pych's intellectual property and thus are subjects of some kinds
// of copyright.

//-----------------------------------------------------------------------
// Ok, just to spill ublas off. Can't imagine how to build that ublas stuff on a winbox...
class IMatrix
{
    double data[12][12];
  public:
    IMatrix() {
        for(size_t r=0; r<12; r++)
            for(size_t c=0; c<12; c++)
                data[r][c] = 0.;
    }
    double operator()(size_t r, size_t c) const {
        assert(r < 12 and c < 12);
        return data[r][c];
    }
    double& operator()(size_t r, size_t c) {
        assert(r < 12 and c < 12);
        return data[r][c];
    }
};

class X
{
    double data[12];
  public:
    X() {
        for(size_t i=0; i<12; i++)
            data[i] = 0.;
    }
    double operator()(size_t i) const {
        assert(i< 12);
        return data[i];
    }
    double& operator()(size_t i) {
        assert(i< 12);
        return data[i];
    }
    X& operator+=(const X& d) {
        for(size_t i=0; i<12; i++)
            data[i] += d(i);
        return *this;
    }
    bool less_then_eps(double eps) const {
        for(size_t i=0; i<12; i++)
            if(std::abs(data[i]) < eps)
                ;
            else
                return false;
        return true;
    }
};

class AMatrix
{
    double data[12][13];
  public:
    AMatrix(const IMatrix& I, const X& W) {
        for(size_t r=0; r<12; r++) {
            for(size_t c=0; c<12; c++)
                data[r][c] = I(r,c);
            data[r][12] = -W(r);
        }
    }
    double operator()(size_t r, size_t c) const {
        assert(r < 12 and c < 13);
        return data[r][c];
    }
    void swap_rows(size_t i, size_t j) {
        assert(i < 12 and j < 12);
        if(i != j)
            for(size_t c=0; c<13; c++)
                std::swap(data[i][c], data[j][c]);
    }
    void normalise_row(size_t i) {
        assert(i < 12);
        const double denominator = data[i][i];
        for(size_t c=0; c<13; c++)
            data[i][c] /= denominator;
    }
    void subst_row_mult(size_t i, size_t j) {
        assert(i < 12 and j < 12);
        const double multtiplier = data[i][j];
        for(size_t c=0; c<13; c++)
            data[i][c] -= data[j][c]*multtiplier;
    }
    X last_column() const {
        X x;
        for(size_t r=0; r<12; r++)
            x(r) = data[r][12];
        return x;
    }
};

//-----------------------------------------------------------------------

static X make_x0(const Params& p)
{
    X x;
    x(0) = 0.; // delta_omega_0
    x(1) = p.start.Delta0;
    x(2) = p.start.Eqprime0;
    x(3) = p.start.Eqe0;
    x(4) = 0.; // x3_0
    x(5) = 0.; // x4_0
    x(6) = 0.; // x5_0
    x(7) = 0.; // x6_0
    x(8) = 0.; // x8_0
    x(9) = 0.; // x9_0
    x(10) = p.start.V0;
    x(11) = p.start.U0;
    return x;
}

//-----------------------------------------------------------------------

static inline bool non_zero(const AMatrix& A, size_t i, size_t j) {  return std::abs(A(i,j)) > 1e-100;  }

static void ensure_non_zero_diagonal_elem(AMatrix& A, size_t n)
{
    if(non_zero(A,n,n))
        return;
    for(size_t i=n+1; i<12/*A.size1()*/; i++)
        if(non_zero(A,i,n))
        {
            A.swap_rows(i,n);
            return;
        }
}

//-----------------------------------------------------------------------

X solve_gauss(const IMatrix& I, const X& W)
{
    AMatrix A(I,W);

    // Well, actually, we know a lot about the A guts and could use some
    // tricks to reduce computations. But here we're just showing how to
    // brutegauss LASes.
    for(size_t i=0; i<12; i++)
    {
        ensure_non_zero_diagonal_elem(A,i);
        A.normalise_row(i);
        for(size_t j=i+1; j<12; j++)
            A.subst_row_mult(j,i);
    }
    for(size_t k=0; k<12; k++)
    {
        const size_t i = 11-k;
        for(size_t j=0; j<i; j++)
            A.subst_row_mult(j,i);
    }
    return A.last_column();
}

//-----------------------------------------------------------------------

struct Equiv
{
    double Y11, Y12, A11, A12;
    double Pd, Upphi;
};

struct NewtonDoesNotConverge {};

template <typename Method> X solve_newton(const X& x_k_1, const Params& p, const Equiv& e)
{
    // I*delta_x + W = 0
    X x_i_1 = x_k_1;
    for(size_t max_iterations=p.max_iterations; max_iterations; max_iterations--)
    {
        const X W = Method::calculate_W(x_k_1, x_i_1, p.dt, p.reg, e);
        if(W.less_then_eps(p.eps))
            return x_i_1;
        const IMatrix I = Method::calculate_I(x_i_1, p.dt, p.reg, e);
        const X delta_x_i = solve_gauss(I, W);
        x_i_1 += delta_x_i;
    }
    throw NewtonDoesNotConverge();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

const static double Pt0 = 1.;
const static double T = 1.;
const static double Ty = 5.1;
const static double Tu = 0.05;
const static double Tg = 0.02;
const static double Te = 0.03;
const static double Tf = 0.05;
const static double Tphi = 0.05;
const static double Td0 = 5.87;
const static double omega_nom = 18000.;
const static double Xdprime = 0.207;
const static double Xd = 1.364;
const static double Xdp = (Xd-Xdprime)/Xd/Xdprime;
const static double Eqenom = 1.87;
const static double Uc = 1.;

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

Equiv equiv_params(double t, const X& x, const Params& p, Equiv e)
{
    e.Y11 = t<0.12 ? p.repl.Y11em : p.repl.Y11;
    e.Y12 = t<0.12 ? p.repl.Y12em : p.repl.Y12;
    e.A11 = t<0.12 ? p.repl.A11em : p.repl.A11;
    e.A12 = t<0.12 ? p.repl.A12em : p.repl.A12;
    e.Pd = p.repl.Pd;
    const double U = x(11);
    if(U < 0.85) e.Upphi = 2.*Eqenom - (2.*Eqenom - p.start.Eqe0)*std::exp(-p.dt/Te);
    if(U > 0.9)  e.Upphi = 0.;
    return e;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct Eiler {

    static X calculate_W(const X& x_k_1, const X& x_i_1, double dt, const Params::Reg& reg, const Equiv& e)
    {
        const double domega_k = x_k_1(0);   const double domega = x_i_1(0);
        const double delta_k = x_k_1(1);    const double delta = x_i_1(1);
        const double Eqprime_k = x_k_1(2);  const double Eqprime = x_i_1(2);
        const double Eqe_k = x_k_1(3);      const double Eqe = x_i_1(3);
        const double X3_k = x_k_1(4);       const double X3 = x_i_1(4);
        const double X4_k = x_k_1(5);       const double X4 = x_i_1(5);
        const double X5_k = x_k_1(6);       const double X5 = x_i_1(6);
        const double X6_k = x_k_1(7);       const double X6 = x_i_1(7);
        const double X8_k = x_k_1(8);       const double X8 = x_i_1(8);
        const double X9_k = x_k_1(9);       const double X9 = x_i_1(9);
                                            const double V = x_i_1(10);
                                            const double U = x_i_1(11);
        const double d_v = delta-V;

        X W;
        W(0) = delta - delta_k - dt*domega;
        W(1) = domega - domega_k - (dt/Ty)*omega_nom*(Pt0 - e.Pd/omega_nom*domega - Eqprime*U/Xdprime*sin(d_v) + U*U*Xdp*sin(d_v)*cos(d_v));
        W(2) = Eqprime - Eqprime_k - (dt/Td0)*(Eqe - Eqprime*Xd/Xdprime + U*Xd*Xdp*cos(d_v));
        W(3) = X3 - X3_k - (dt/Tu)*(reg.K0u*(reg.Ur0 - U) - X3);
        W(4) = X4 - X4_k - (dt/Tg)*(reg.K1u/Tu*(reg.Ur0 - U - X3/reg.K0u) - X4);
        W(5) = X5 - X5_k - (dt/Tphi)*(V - X5);
        W(6) = X6 - X6_k - (dt/Tf)*((V - X5)/Tphi - X6);
        W(7) = X8 - X8_k - (dt/Tf)*reg.K0f*((V - X5)/Tphi - X6) - dt/T*X8;
        W(8) = X9 - X9_k - (dt/Tg)*(reg.K1f/Tf*((V - X5)/Tphi - X6) - reg.K1f/reg.K0f/T*X8 - X9);
        W(9) = Eqe - Eqe_k - (dt/Te)*(X3 + X4 + X8 + X9 + e.Upphi - Eqe);
        W(10) = Eqprime*U/Xdprime*sin(d_v) - U*U*Xdp*sin(d_v)*cos(d_v) - U*U*e.Y11*sin(e.A11) - U*Uc*e.Y12*sin(V-e.A12);
        W(11) = Eqprime*U/Xdprime*cos(d_v) - U*U/Xdprime + U*U*Xdp*sin(d_v)*sin(d_v) - U*U*e.Y11*cos(e.A11) + U*Uc*e.Y12*cos(V-e.A12);
        return W;
    }
    static IMatrix calculate_I(const X& x_i_1, double dt, const Params::Reg& reg, const Equiv& e)
    {
        IMatrix I;
        //I = boost::numeric::ublas::zero_matrix<double>(12,12);

        const double delta = x_i_1(1);
        const double Eqprime = x_i_1(2);
        const double V = x_i_1(10);
        const double U = x_i_1(11);
        const double d_v = delta-V;

        I(0,0) = -dt;
        I(0,1) = 1.;

        I(1,0) = (dt/Ty)*e.Pd + 1.;
        I(1,1) = (dt/Ty)*omega_nom*(Eqprime*U/Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v));
        I(1,2) = (dt/Ty)*omega_nom*U/Xdprime*sin(d_v);
        I(1,10) = -I(1,1);
        I(1,11) = (dt/Ty)*omega_nom*(Eqprime/Xdprime*sin(d_v) - U*Xdp*sin(2*d_v));

        I(2,1) = (dt/Td0)*U*Xd*Xdp*sin(d_v);
        I(2,2) = (dt/Td0)*Xd/Xdprime + 1.;
        I(2,3) = -(dt/Td0);
        I(2,10) = -I(2,1);
        I(2,11) = (dt/Td0)*Xd*Xdp*cos(d_v);

        I(3,4) = (dt/Tu) + 1.;
        I(3,11) = (dt/Tu)*reg.K0u;

        I(4,4) = (dt/Tg)*reg.K1u/reg.K0u/Tu;
        I(4,5) = (dt/Tg) + 1.;
        I(4,11) = I(4,4)*reg.K0u;

        I(5,6) = (dt/Tphi) + 1.;
        I(5,10) = -(dt/Tphi);

        I(6,6) = (dt/Tphi)/Tf;
        I(6,7) = (dt/Tf) + 1.;
        I(6,10) = -I(6,6);

        I(7,6) = (dt/Tphi)/Tf*reg.K0f;
        I(7,7) = I(7,6)*Tphi;
        I(7,8) = 1. - dt/T;
        I(7,10) = -I(7,6);

        I(8,6) = (dt/Tg)*reg.K1f/Tf/Tphi;
        I(8,7) = I(8,6)*Tphi;
        I(8,8) = (dt/Tg)*reg.K1f/reg.K0f/T;
        I(8,9) = (dt/Tg) + 1.;
        I(8,10) = -I(8,6);

        I(9,3) = (dt/Te) + 1.;
        I(9,4) = -(dt/Te);
        I(9,5) = -(dt/Te);
        I(9,8) = -(dt/Te);
        I(9,9) = -(dt/Te);

        I(10,1) = Eqprime*U/Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v);
        I(10,2) = U/Xdprime*sin(d_v);
        I(10,10) = -I(10,1) - U*Uc*e.Y12*cos(V-e.A12);
        I(10,11) = -Uc*e.Y12*sin(V-e.A12) - 2.*U*e.Y11*sin(e.A11) + Eqprime/Xdprime*sin(d_v) - U*Xdp*sin(2.*d_v);

        I(11,1) = U*U*Xdp*sin(2*d_v) - Eqprime*U/Xdprime*sin(d_v);
        I(11,2) = U/Xdprime*cos(d_v);
        I(11,10) = -I(11,1) - U*Uc*e.Y12*sin(V-e.A12);
        I(11,11) = Eqprime/Xdprime*cos(d_v) - 2.*U/Xdprime + 2.*U*Xdp*sin(d_v)*sin(d_v) + Uc*e.Y12*cos(V-e.A12) - 2.*U*e.Y11*cos(e.A11);
        return I;
    }
};

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

static AnswerItem make_answer_item(double t, const X& x)
{
    AnswerItem i = { 0., 0., 0., 0., 0., 0., 0. };
    i.time=t, i.delta=x(1), i.omega=x(0), i.Eqe=x(3), i.Eqprime=x(2), i.V=x(10), i.U=x(11);
    return i;
}

struct RatherSuspiciousEqe
{
    double value;
    RatherSuspiciousEqe(double v) : value(v) {}
};

QVector<AnswerItem> CalculusEiler::doWork(const Params& p)
{
    QVector<AnswerItem> a;
    a.reserve((p.Tstop-p.Tstart)/p.dt+1);

    X x = make_x0(p);
    a.push_back( make_answer_item(p.Tstart,x) );

    const Equiv e = { 0., 0., 0., 0., 0., 0. };
    for(double t=p.Tstart+p.dt; t<p.Tstop; t+=p.dt) try
    {
        x = solve_newton<Eiler>(x, p, equiv_params(t,x,p,e));
        a.push_back( make_answer_item(t,x) );
        emit a_step_done();
        if(false and std::abs(a.back().Eqe) > 1e2)
            throw RatherSuspiciousEqe(a.back().Eqe);
    }
    catch(NewtonDoesNotConverge)
    {
        qCritical() << "Eiler" << a.size() << "points: aaaaand NewtonDoesNotConverge";
        break;
        //throw;
    }
    return a;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct Trapeze {

    static X calculate_W(const X& x_k_1, const X& x_i_1, double dt, const Params::Reg& reg, const Equiv& e)
    {
        const double domega_k = x_k_1(0);   const double domega = x_i_1(0);
        const double delta_k = x_k_1(1);    const double delta = x_i_1(1);
        const double Eqprime_k = x_k_1(2);  const double Eqprime = x_i_1(2);
        const double Eqe_k = x_k_1(3);      const double Eqe = x_i_1(3);
        const double X3_k = x_k_1(4);       const double X3 = x_i_1(4);
        const double X4_k = x_k_1(5);       const double X4 = x_i_1(5);
        const double X5_k = x_k_1(6);       const double X5 = x_i_1(6);
        const double X6_k = x_k_1(7);       const double X6 = x_i_1(7);
        const double X8_k = x_k_1(8);       const double X8 = x_i_1(8);
        const double X9_k = x_k_1(9);       const double X9 = x_i_1(9);
        const double V_k = x_k_1(10);       const double V = x_i_1(10);
        const double U_k = x_k_1(11);       const double U = x_i_1(11);
        const double d_v_k = delta_k-V_k;   const double d_v = delta-V;

        X W;
        W(0) = delta - delta_k - dt/2.*(domega + domega_k);
        W(1) = domega - domega_k - dt/Ty/2.*omega_nom*(Pt0 + Pt0 - e.Pd/omega_nom*(domega + domega_k)
                - Eqprime*U/Xdprime*sin(d_v) - Eqprime_k*U_k/Xdprime*sin(d_v_k)
                + U*U*Xdp*sin(d_v)*cos(d_v) + U_k*U_k*Xdp*sin(d_v_k)*cos(d_v_k));
        W(2) = Eqprime - Eqprime_k - dt/Td0/2.*(Eqe + Eqe_k - (Eqprime + Eqprime_k)*Xd/Xdprime
                + U*Xd*Xdp*cos(d_v) + U_k*Xd*Xdp*cos(d_v_k));
        W(3) = X3 - X3_k - dt/Tu/2.*(2.*reg.K0u*reg.Ur0 - reg.K0u*(U + U_k) - X3 - X3_k);
        W(4) = X4 - X4_k - dt/Tg/2.*(2.*reg.K1u/Tu*reg.Ur0 - reg.K1u/Tu*(U + U_k) - reg.K1u/reg.K0u/Tu*(X3 + X3_k) - X4 - X4_k);
        W(5) = X5 - X5_k - dt/Tphi/2.*(V + V_k - X5 - X5_k);
        W(6) = X6 - X6_k - dt/Tf/2.*((V + V_k - X5 - X5_k)/Tphi - X6 - X6_k);
        W(7) = X8 - X8_k - dt/Tf/2.*reg.K0f*((V + V_k - X5 - X5_k)/Tphi - X6 - X6_k) - dt/T*(X8 + X8_k);
        W(8) = X9 - X9_k - dt/Tg/2.*(reg.K1f/Tf*((V + V_k - X5 - X5_k)/Tphi - X6 - X6_k)
                - reg.K1f/reg.K0f/T*(X8 + X8_k) - X9 - X9_k);
        W(9) = Eqe - Eqe_k - dt/Te/2.*(X3 + X3_k + X4 + X4_k + X8 + X8_k + X9 + X9_k + 2.*e.Upphi - Eqe - Eqe_k);
        W(10) = Eqprime*U/Xdprime*sin(d_v) - U*U*Xdp*sin(d_v)*cos(d_v) - U*U*e.Y11*sin(e.A11) - U*Uc*e.Y12*sin(V-e.A12);
        W(11) = Eqprime*U/Xdprime*cos(d_v) - U*U/Xdprime + U*U*Xdp*sin(d_v)*sin(d_v) - U*U*e.Y11*cos(e.A11) + U*Uc*e.Y12*cos(V-e.A12);
        return W;
    }
    static IMatrix calculate_I(const X& x_i_1, double dt, const Params::Reg& reg, const Equiv& e)
    {
        IMatrix I;
        //I = boost::numeric::ublas::zero_matrix<double>(12,12);

        const double delta = x_i_1(1);
        const double Eqprime = x_i_1(2);
        const double V = x_i_1(10);
        const double U = x_i_1(11);
        const double d_v = delta-V;

        I(0,0) = -dt/2.;
        I(0,1) = 1.;

        I(1,0) = dt/Ty/2.*e.Pd + 1.;
        I(1,1) = dt/Ty/2.*omega_nom*(Eqprime*U/Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v));
        I(1,2) = dt/Ty/2.*omega_nom*U/Xdprime*sin(d_v);
        I(1,10) = -I(1,1);
        I(1,11) = dt/Ty/2.*omega_nom*(Eqprime/Xdprime*sin(d_v) - U*Xdp*sin(2*d_v));

        I(2,1) = dt/Td0/2.*U*Xd*Xdp*sin(d_v);
        I(2,2) = dt/Td0/2.*Xd/Xdprime + 1.;
        I(2,3) = -dt/Td0/2.;
        I(2,10) = -I(2,1);
        I(2,11) = dt/Td0/2.*Xd*Xdp*cos(d_v);

        I(3,4) = dt/Tu/2. + 1.;
        I(3,11) = dt/Tu/2.*reg.K0u;

        I(4,4) = dt/Tg/2.*reg.K1u/reg.K0u/Tu;
        I(4,5) = dt/Tg/2. + 1.;
        I(4,11) = I(4,4)*reg.K0u;

        I(5,6) = dt/Tphi/2. + 1.;
        I(5,10) = -dt/Tphi/2.;

        I(6,6) = dt/Tf/2./Tphi;
        I(6,7) = dt/Tf/2. + 1.;
        I(6,10) = -I(6,6);

        I(7,6) = dt/Tphi/Tf/2.*reg.K0f;
        I(7,7) = I(7,6)*Tphi;
        I(7,8) = 1. - dt/T/2.;
        I(7,10) = -I(7,6);

        I(8,6) = dt/Tg/Tf/Tphi/2.*reg.K1f;
        I(8,7) = I(8,6)*Tphi;
        I(8,8) = dt/Tg/2.*reg.K1f/reg.K0f/T;
        I(8,9) = dt/Tg/2. + 1.;
        I(8,10) = -I(8,6);

        I(9,3) = dt/Te/2. + 1.;
        I(9,4) = -dt/Te/2.;
        I(9,5) = -dt/Te/2.;
        I(9,8) = -dt/Te/2.;
        I(9,9) = -dt/Te/2.;

        I(10,1) = Eqprime*U/Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v);
        I(10,2) = U/Xdprime*sin(d_v);
        I(10,10) = -I(10,1) - U*Uc*e.Y12*cos(V-e.A12);
        I(10,11) = -Uc*e.Y12*sin(V-e.A12) - 2.*U*e.Y11*sin(e.A11) + Eqprime/Xdprime*sin(d_v) - U*Xdp*sin(2.*d_v);

        I(11,1) = U*U*Xdp*sin(2*d_v) - Eqprime*U/Xdprime*sin(d_v);
        I(11,2) = U/Xdprime*cos(d_v);
        I(11,10) = -I(11,1) - U*Uc*e.Y12*sin(V-e.A12);
        I(11,11) = Eqprime/Xdprime*cos(d_v) - 2.*U/Xdprime + 2.*U*Xdp*sin(d_v)*sin(d_v) + Uc*e.Y12*cos(V-e.A12) - 2.*U*e.Y11*cos(e.A11);
        return I;
    }
};

//-----------------------------------------------------------------------

QVector<AnswerItem> CalculusTrapeze::doWork(const Params& p)
{
    QVector<AnswerItem> a;
    a.reserve((p.Tstop-p.Tstart)/p.dt+1);

    X x = make_x0(p);
    a.push_back( make_answer_item(p.Tstart,x) );

    const Equiv e = { 0., 0., 0., 0., 0., 0. };
    for(double t=p.Tstart+p.dt; t<p.Tstop; t+=p.dt) try
    {
        x = solve_newton<Trapeze>(x, p, equiv_params(t,x,p,e));
        a.push_back( make_answer_item(t,x) );
        emit a_step_done();
    }
    catch(NewtonDoesNotConverge)
    {
        qDebug() << "Trapeze" << a.size() << "points: NewtonDoesNotConverge";
        break;
        //throw;
    }
    return a;
}

