#include "Calculus.h"
#include "params.h"
#include <QDebug>
#include <cassert>
#include <cmath>

// the Trapeze and Eiler matrix and vector statements are fully and trully
// Lyuba Chich-Pych's intellectual property and thus are subjects of some
// kinds of copyright.

//-----------------------------------------------------------------------
// Ok, just to spill ublas off. Can't imagine how to build that ublas stuff on a winbox...
template <size_t R, size_t C> class Matrix_N
{
    double data[R][C];
  public:
    Matrix_N() {
        for(size_t r=0; r<R; r++)
            for(size_t c=0; c<C; c++)
                data[r][c] = 0.;
    }
    double operator()(size_t r, size_t c) const {
        assert(r < R and c < C);
        return data[r][c];
    }
    double& operator()(size_t r, size_t c) {
        assert(r < R);
        assert(c < C);
        assert(r < R and c < C);
        return data[r][c];
    }
};

template <size_t N> struct IMatrix_N : public Matrix_N<N,N> {};

template <size_t N> class X_N : public Matrix_N<N,1>
{
  public:
    double operator()(size_t i) const {
        return Matrix_N<N,1>::operator()(i,0);
    }
    double& operator()(size_t i) {
        return Matrix_N<N,1>::operator()(i,0);
    }
    X_N& operator+=(const X_N& d) {
        for(size_t i=0; i<N; i++)
            (*this)(i) += d(i);
        return *this;
    }
    bool less_then_eps(double eps) const {
        for(size_t i=0; i<N; i++)
        {
            if(std::isnan((*this)(i)) or std::abs((*this)(i)) > eps)
                return false;
        }
        return true;
    }
};

template <size_t N> class AMatrix_N : public Matrix_N<N,N+1>
{
  public:
    AMatrix_N(const IMatrix_N<N>& I, const X_N<N>& W) {
        for(size_t r=0; r<N; r++) {
            for(size_t c=0; c<N; c++)
                (*this)(r,c) = I(r,c);
            (*this)(r,N) = -W(r);
        }
    }
    void ensure_non_zero_diagonal_elem(size_t n, const double eps = 1e-100) {
        if(std::abs((*this)(n,n)) > eps)
            return;
        for(size_t i=n+1; i<N; i++)
            if(std::abs((*this)(i,n)) > eps)
            {
                for(size_t c=0; c<=N; c++)
                    std::swap((*this)(i,c), (*this)(n,c));
                return;
            }
    }
    void normalise_row(size_t i) {
        const double denominator = (*this)(i,i);
        for(size_t c=0; c<=N; c++)
            (*this)(i,c) /= denominator;
    }
    void subst_row_mult(size_t i, size_t j) {
        const double multtiplier = (*this)(i,j);
        for(size_t c=0; c<=N; c++)
            (*this)(i,c) -= (*this)(j,c)*multtiplier;
    }
    X_N<N> last_column() const {
        X_N<N> x;
        for(size_t r=0; r<N; r++)
            x(r) = (*this)(r,N);
        return x;
    }
};

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

// gaussforcing a NLAS
template<size_t N> X_N<N> solve_gauss(const IMatrix_N<N>& I, const X_N<N>& W)
{
    // the AMatrix type designed specially to operate with the newton source
    // data in a convenient way. After construction, it manipulates its rows
    // consistently and returns the answer with properly applied signs.
    AMatrix_N<N> A(I,W);

    // Well, actually, we know a lot about the A guts and could use some
    // tricks to reduce computations. But here we're just showing how to
    // brutegauss LASes. A pretty straightforward gauss elimination
    // algorithm is implemented.

    // Rolling down, converting the matrix into an upper-triangle one
    for(size_t i=0; i<N; i++)
    {
        A.ensure_non_zero_diagonal_elem(i);
        A.normalise_row(i);
        for(size_t j=i+1; j<N; j++)
            A.subst_row_mult(j,i);
    }
    // Rolling up, making the matrix diagonal
    for(size_t k=0; k<N; k++)
    {
        const size_t i = N-1-k;
        for(size_t j=0; j<i; j++)
            A.subst_row_mult(j,i);
    }
    return A.last_column();
}

//-----------------------------------------------------------------------

struct NewtonDoesNotConverge {
    NewtonDoesNotConverge(size_t n) : n_steps(n) {}
    const size_t n_steps;
};

// a classical Newton solver
// templated to provide the same framework for either methods
template <class Method> typename Method::X solve_newton_impl(const Method* pimpl, size_t& n_steps)
{
    // incrementally solving "I*delta_x + W = 0 (where I = W')" until both W and delta_x > eps
    typename Method::X x_i_1 = pimpl->x;
    n_steps = 0;
    for(size_t max_iterations=pimpl->p.max_iterations; max_iterations; max_iterations--)
    {
        const typename Method::X W = pimpl->calculate_W(x_i_1);
        const typename Method::IMatrix I = pimpl->calculate_I(x_i_1);
        // converts our NLAS into a LAS and solves it
        const typename Method::X delta_x_i = solve_gauss(I, W);
        x_i_1 += delta_x_i;
        n_steps++;
        if(W.less_then_eps(pimpl->p.eps))// and delta_x_i.less_then_eps(pimpl->p.eps))
            return x_i_1; // luckily, the system depends on the x_k_1 so we can return the new
                          // step value, omitting (x_i_1 - x_k_1) -> x_k_1 + delta_x
    }
    throw NewtonDoesNotConverge(n_steps);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

Equiv recalculate_equiv_params(double t, const double U, const Params& p, Equiv e)
{
    e.Y11 = t<0.12 ? p.repl.Y11em : p.repl.Y11;
    e.Y12 = t<0.12 ? p.repl.Y12em : p.repl.Y12;
    e.A11 = t<0.12 ? p.repl.A11em : p.repl.A11;
    e.A12 = t<0.12 ? p.repl.A12em : p.repl.A12;
    e.Pd = p.repl.Pd;
    if(U < 0.85) {
        if(e.Tsw_low < 0) {
            e.Tsw_low = t + 0.0;
        } else if(t > e.Tsw_low) {
            e.Upphi = 2.*p.reg.Eqenom - (2.*p.reg.Eqenom - p.start.Eqe0)*std::exp(-t/p.reg.Te);
            e.Tsw_low = -1.;
        }
    }
    if(U > 0.9) {
        if(e.Tsw_high < 0) {
            e.Tsw_high = t + 0.1;
        } else if(t > e.Tsw_high) {
            e.Upphi = 0.;
            e.Tsw_high = -1.;
        }
    }
    return e;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

// does the real work and emits signal on every time step to update a progress bar
void Calculus::run()
{
    size_t n=0, n_steps=0;
    emit_x(p.Tstart, n++, n_steps); // initial vector according to the blueprint

    for(double t=p.Tstart+p.dt; t<p.Tstop; t+=p.dt, n++) try
    {
        solve_newton(t, n_steps);
        emit_x(t, n, n_steps);
    }
    catch(const NewtonDoesNotConverge& exc)
    {
        emit newton_does_not_converge(name(), t, static_cast<unsigned>(exc.n_steps));
        return;
    }
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct CalculusEiler::Impl
{
    typedef X_N<12> X;
    typedef IMatrix_N<12> IMatrix;
    typedef AMatrix_N<12> AMatrix;

    X calculate_W(const X& x_i_1) const;
    IMatrix calculate_I(const X& x_i_1) const;
    void dirty_hack();

    Impl(const Params& _p) : p(_p)
    {
        x(0) = 0.; // delta_omega_0
        x(1) = p.start.Delta0;
        x(2) = p.start.Eqprime0;
        x(3) = p.start.Eqe0;
        x(4) = p.reg.K0u*(p.reg.Ur0 - p.start.U0); // x3(0)
        x(6) = p.start.V0; // x5(0)
	 //   x(6) = 0.508;
        x(10) = p.start.V0;
        x(11) = p.start.U0;
    }

    X x;
    Equiv e;
    Params p;
    size_t n_steps;
};

void CalculusEiler::Impl::dirty_hack()
{
    if(not p.dirty_hack)
        return;
    x(3) = qMin(x(3), 2*p.reg.Eqenom);
    x(3) = qMax(x(3), 0.);
}

double CalculusEiler::calculate_Pg(double Eqprime, double U,
                                   double Xdprime, double Xd,
                                   double s_d_v, double c_d_v)
{
    const double Xdp = (Xd-Xdprime)/Xd/Xdprime;
    return Eqprime*U/Xdprime*s_d_v - U*U*Xdp*s_d_v*c_d_v;
}

double CalculusEiler::calculate_Pc(double U, double V,
                                   double Y11, double Y12,
                                   double A11, double A12,
                                   double Uc)
{
    return U*U*Y11*sin(A11) + U*Uc*Y12*sin(V-A12);
}

double CalculusEiler::calculate_Qg(double Eqprime, double U,
                                   double Xdprime, double Xd,
                                   double s_d_v, double c_d_v)
{
    return Eqprime*U/Xdprime*c_d_v - U*U*(c_d_v*c_d_v/Xdprime + s_d_v*s_d_v/Xd);
}

double CalculusEiler::calculate_Qc(double U, double V,
                                   double Y11, double Y12,
                                   double A11, double A12,
                                   double Uc)
{
    return U*U*Y11*cos(A11) - U*Uc*Y12*cos(V-A12);
}

CalculusEiler::Impl::X CalculusEiler::Impl::calculate_W(const X& x_i_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;
    const X& x_k_1 = x;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

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
    const double s_d_v = sin(d_v);
    const double c_d_v = cos(d_v);

    X W;
    W(0) = delta - delta_k - dt*domega;
    W(1) = domega - domega_k - (dt/r.Ty)*r.omega_nom*(r.Pt0 - e.Pd/r.omega_nom*domega - Eqprime*U/r.Xdprime*sin(d_v) + U*U*Xdp*sin(d_v)*cos(d_v));
    W(2) = Eqprime - Eqprime_k - (dt/r.Td0)*(Eqe - Eqprime*r.Xd/r.Xdprime + U*r.Xd*Xdp*cos(d_v));
    W(3) = X3 - X3_k - (dt/r.Tu)*(r.K0u*(r.Ur0 - U) - X3);
    W(4) = X4 - X4_k - (dt/r.Tg)*(r.K1u/r.Tu*(r.Ur0 - U - X3/r.K0u) - X4);
    W(5) = X5 - X5_k - (dt/r.Tphi)*(V - X5);
    W(6) = X6 - X6_k - (dt/r.Tf)*((V - X5)/r.Tphi - X6);
    W(7) = X8 - X8_k - (dt/r.Tf)*r.K0f*((V - X5)/r.Tphi - X6) + dt/r.T*X8;
    W(8) = X9 - X9_k - (dt/r.Tg)*(r.K1f/r.Tf*((V - X5)/r.Tphi - X6) - r.K1f/r.K0f/r.T*X8 - X9);
    W(9) = Eqe - Eqe_k - (dt/r.Te)*(X3 + X4 + X8 + X9 + p.start.Eqe0 + e.Upphi - Eqe);
    W(10) = calculate_Pg(Eqprime, U, r.Xdprime, r.Xd, s_d_v, c_d_v) -
            calculate_Pc(U, V, e.Y11, e.Y12, e.A11, e.A12, r.Uc);
    //qDebug()<<"W(10) ="<<calculate_Pg(Eqprime, U, r.Xdprime, r.Xd, s_d_v, c_d_v)<<"-"
    //                   <<calculate_Pc(U, V, e.Y11, e.Y12, e.A11, e.A12, r.Uc)<<"="<<W(10);
    W(11) = calculate_Qg(Eqprime, U, r.Xdprime, r.Xd, s_d_v, c_d_v) -
            calculate_Qc(U, V, e.Y11, e.Y12, e.A11, e.A12, r.Uc);
    //qDebug()<<"W(11) ="<<calculate_Qg(Eqprime, U, r.Xdprime, r.Xd, s_d_v, c_d_v)<<"-"
    //                   <<calculate_Qc(U, V, e.Y11, e.Y12, e.A11, e.A12, r.Uc)<<"="<<W(11);
    return W;
}

CalculusEiler::Impl::IMatrix CalculusEiler::Impl::calculate_I(const X& x_i_1) const
{
    const double dt = p.dt;
    const Params::Consts& r = p.reg;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

    const double delta = x_i_1(1);
    const double Eqprime = x_i_1(2);
    const double V = x_i_1(10);
    const double U = x_i_1(11);
    const double d_v = delta-V;

    IMatrix I;

    I(0,0) = -dt;
    I(0,1) = 1.;

    I(1,0) = (dt/r.Ty)*e.Pd + 1.;
    I(1,1) = (dt/r.Ty)*r.omega_nom*(Eqprime*U/r.Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v));
    I(1,2) = (dt/r.Ty)*r.omega_nom*U/r.Xdprime*sin(d_v);
    I(1,10) = -I(1,1);
    I(1,11) = (dt/r.Ty)*r.omega_nom*(Eqprime/r.Xdprime*sin(d_v) - U*Xdp*sin(2*d_v));

    I(2,1) = (dt/r.Td0)*U*r.Xd*Xdp*sin(d_v);
    I(2,2) = (dt/r.Td0)*r.Xd/r.Xdprime + 1.;
    I(2,3) = -(dt/r.Td0);
    I(2,10) = -I(2,1);
    I(2,11) = -(dt/r.Td0)*r.Xd*Xdp*cos(d_v);

    I(3,4) = (dt/r.Tu) + 1.;
    I(3,11) = (dt/r.Tu)*r.K0u;

    I(4,4) = (dt/r.Tg)*r.K1u/r.K0u/r.Tu;
    I(4,5) = (dt/r.Tg) + 1.;
    I(4,11) = (dt/r.Tg)*r.K1u/r.Tu;

    I(5,6) = (dt/r.Tphi) + 1.;
    I(5,10) = -(dt/r.Tphi);

    I(6,6) = (dt/r.Tphi)/r.Tf;
    I(6,7) = (dt/r.Tf) + 1.;
    I(6,10) = -I(6,6);

    I(7,6) = (dt/r.Tphi)/r.Tf*r.K0f;
    I(7,7) = I(7,6)*r.Tphi;
    I(7,8) = 1. + dt/r.T;
    I(7,10) = -I(7,6);

    I(8,6) = (dt/r.Tg)*r.K1f/r.Tf/r.Tphi;
    I(8,7) = I(8,6)*r.Tphi;
    I(8,8) = (dt/r.Tg)*r.K1f/r.K0f/r.T;
    I(8,9) = (dt/r.Tg) + 1.;
    I(8,10) = -I(8,6);

    I(9,3) = (dt/r.Te) + 1.;
    I(9,4) = -(dt/r.Te);
    I(9,5) = -(dt/r.Te);
    I(9,8) = -(dt/r.Te);
    I(9,9) = -(dt/r.Te);

    I(10,1) = Eqprime*U/r.Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v);
    I(10,2) = U/r.Xdprime*sin(d_v);
    I(10,10) = -I(10,1) - U*r.Uc*e.Y12*cos(V-e.A12);
    I(10,11) = -r.Uc*e.Y12*sin(V-e.A12) - 2.*U*e.Y11*sin(e.A11) + Eqprime/r.Xdprime*sin(d_v) - U*Xdp*sin(2.*d_v);

    I(11,1) = U*U*Xdp*sin(2*d_v) - Eqprime*U/r.Xdprime*sin(d_v);
    I(11,2) = U/r.Xdprime*cos(d_v);
    I(11,10) = -I(11,1) - U*r.Uc*e.Y12*sin(V-e.A12);
    I(11,11) = Eqprime/r.Xdprime*cos(d_v) - 2*U*(cos(d_v)*cos(d_v)/r.Xdprime + sin(d_v)*sin(d_v)/r.Xd) - 2*U*e.Y11*cos(e.A11) + r.Uc*e.Y12*cos(V-e.A12);
    return I;
}

//-----------------------------------------------------------------------

CalculusEiler::CalculusEiler(const Params& p) : Calculus(p), pimpl(new Impl(p)) {}

void CalculusEiler::emit_x(double t, size_t row, size_t n_steps)
{
    AnswerItem i;
    i.time=t, i.row=row, i.set_no=0, i.n_steps=n_steps;
    i.delta=pimpl->x(1), i.omega=pimpl->x(0), i.Eqe=pimpl->x(3);
    i.Eqprime=pimpl->x(2), i.V=pimpl->x(10), i.U=pimpl->x(11);
    emit a_step_done(i);
}

void CalculusEiler::solve_newton(double t, size_t& n_steps)
{
    // recalculating and reassigning 'e' at every step is essential to get some hysteresis behaviour
    pimpl->e = recalculate_equiv_params(t, pimpl->x(11), pimpl->p, pimpl->e);
    // as our newton implementation is able to return a new step value, omit some arithmetics
    pimpl->x = solve_newton_impl(pimpl, n_steps);
    pimpl->dirty_hack();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct CalculusTrapeze::Impl
{
    typedef X_N<12> X;
    typedef IMatrix_N<12> IMatrix;
    typedef AMatrix_N<12> AMatrix;

    X calculate_W(const X& x_i_1) const;
    IMatrix calculate_I(const X& x_i_1) const;
    void dirty_hack();

    Impl(const Params& _p) : p(_p)
    {
        x(0) = 0.; // delta_omega_0
        x(1) = p.start.Delta0;
        x(2) = p.start.Eqprime0;
        x(3) = p.start.Eqe0;
        x(4) = p.reg.K0u*(p.reg.Ur0 - p.start.U0); // x3(0)
        x(6) = p.start.V0; // x5(0)
        x(10) = p.start.V0;
        x(11) = p.start.U0;
    }

    X x;
    Equiv e;
    Params p;
    size_t n_steps;
};

void CalculusTrapeze::Impl::dirty_hack()
{
    if(not p.dirty_hack)
        return;
    x(3) = qMin(x(3), 2*p.reg.Eqenom);
    x(3) = qMax(x(3), 0.);
}

CalculusTrapeze::Impl::X CalculusTrapeze::Impl::calculate_W(const X& x_i_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;
    const X& x_k_1 = x;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

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
    W(1) = domega - domega_k - dt/r.Ty/2.*r.omega_nom*(r.Pt0 + r.Pt0 - e.Pd/r.omega_nom*(domega + domega_k)
            - Eqprime*U/r.Xdprime*sin(d_v) - Eqprime_k*U_k/r.Xdprime*sin(d_v_k)
            + U*U*Xdp*sin(d_v)*cos(d_v) + U_k*U_k*Xdp*sin(d_v_k)*cos(d_v_k));
    W(2) = Eqprime - Eqprime_k - dt/r.Td0/2.*(Eqe + Eqe_k - (Eqprime + Eqprime_k)*r.Xd/r.Xdprime
            + U*r.Xd*Xdp*cos(d_v) + U_k*r.Xd*Xdp*cos(d_v_k));
    W(3) = X3 - X3_k - dt/r.Tu/2.*(2.*r.K0u*r.Ur0 - r.K0u*(U + U_k) - X3 - X3_k);
    W(4) = X4 - X4_k - dt/r.Tg/2.*(2.*r.K1u/r.Tu*r.Ur0 - r.K1u/r.Tu*(U + U_k) - r.K1u/r.K0u/r.Tu*(X3 + X3_k) - X4 - X4_k);
    W(5) = X5 - X5_k - dt/r.Tphi/2.*(V + V_k - X5 - X5_k);
    W(6) = X6 - X6_k - dt/r.Tf/2.*((V + V_k - X5 - X5_k)/r.Tphi - X6 - X6_k);
    W(7) = X8 - X8_k - dt/r.Tf/2.*r.K0f*((V + V_k - X5 - X5_k)/r.Tphi - X6 - X6_k) + dt/r.T*(X8 + X8_k);
    W(8) = X9 - X9_k - dt/r.Tg/2.*(r.K1f/r.Tf*((V + V_k - X5 - X5_k)/r.Tphi - X6 - X6_k)
            - r.K1f/r.K0f/r.T*(X8 + X8_k) - X9 - X9_k);
    W(9) = Eqe - Eqe_k - dt/r.Te/2.*(X3 + X3_k + X4 + X4_k + X8 + X8_k + X9 + X9_k + 2.*e.Upphi + 2.*p.start.Eqe0 - Eqe - Eqe_k);
    W(10) = Eqprime*U/r.Xdprime*sin(d_v) - U*U*Xdp*sin(2*d_v)/2 - U*U*e.Y11*sin(e.A11) - U*r.Uc*e.Y12*sin(V-e.A12);
    W(11) = Eqprime*U/r.Xdprime*cos(d_v) - U*U/r.Xdprime*cos(d_v)*cos(d_v) - U*U/r.Xd*sin(d_v)*sin(d_v) - U*U*e.Y11*cos(e.A11) + U*r.Uc*e.Y12*cos(V-e.A12);
    return W;
}

CalculusTrapeze::Impl::IMatrix CalculusTrapeze::Impl::calculate_I(const X& x_i_1) const
{
    const double dt = p.dt;
    const Params::Consts& r = p.reg;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

    const double delta = x_i_1(1);
    const double Eqprime = x_i_1(2);
    const double V = x_i_1(10);
    const double U = x_i_1(11);
    const double d_v = delta-V;

    IMatrix I;

    I(0,0) = -dt/2.;
    I(0,1) = 1.;

    I(1,0) = dt/r.Ty/2.*e.Pd + 1.;
    I(1,1) = dt/r.Ty/2.*r.omega_nom*(Eqprime*U/r.Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v));
    I(1,2) = dt/r.Ty/2.*r.omega_nom*U/r.Xdprime*sin(d_v);
    I(1,10) = -I(1,1);
    I(1,11) = dt/r.Ty/2.*r.omega_nom*(Eqprime/r.Xdprime*sin(d_v) - U*Xdp*sin(2*d_v));

    I(2,1) = dt/r.Td0/2.*U*r.Xd*Xdp*sin(d_v);
    I(2,2) = dt/r.Td0/2.*r.Xd/r.Xdprime + 1.;
    I(2,3) = -dt/r.Td0/2.;
    I(2,10) = -I(2,1);
    I(2,11) = dt/r.Td0/2.*r.Xd*Xdp*cos(d_v);

    I(3,4) = dt/r.Tu/2. + 1.;
    I(3,11) = dt/r.Tu/2.*r.K0u;

    I(4,4) = dt/r.Tg/2.*r.K1u/r.K0u/r.Tu;
    I(4,5) = dt/r.Tg/2. + 1.;
    I(4,11) = I(4,4)*r.K0u;

    I(5,6) = dt/r.Tphi/2. + 1.;
    I(5,10) = -dt/r.Tphi/2.;

    I(6,6) = dt/r.Tf/2./r.Tphi;
    I(6,7) = dt/r.Tf/2. + 1.;
    I(6,10) = -I(6,6);

    I(7,6) = dt/r.Tphi/r.Tf/2.*r.K0f;
    I(7,7) = I(7,6)*r.Tphi;
    I(7,8) = 1. + dt/r.T/2.;
    I(7,10) = -I(7,6);

    I(8,6) = dt/r.Tg/r.Tf/r.Tphi/2.*r.K1f;
    I(8,7) = I(8,6)*r.Tphi;
    I(8,8) = dt/r.Tg/2.*r.K1f/r.K0f/r.T;
    I(8,9) = dt/r.Tg/2. + 1.;
    I(8,10) = -I(8,6);

    I(9,3) = dt/r.Te/2. + 1.;
    I(9,4) = -dt/r.Te/2.;
    I(9,5) = -dt/r.Te/2.;
    I(9,8) = -dt/r.Te/2.;
    I(9,9) = -dt/r.Te/2.;

    I(10,1) = Eqprime*U/r.Xdprime*cos(d_v) - U*U*Xdp*cos(2*d_v);
    I(10,2) = U/r.Xdprime*sin(d_v);
    I(10,10) = -I(10,1) - U*r.Uc*e.Y12*cos(V-e.A12);
    I(10,11) = -r.Uc*e.Y12*sin(V-e.A12) - 2.*U*e.Y11*sin(e.A11) + Eqprime/r.Xdprime*sin(d_v) - U*Xdp*sin(2.*d_v);

    I(11,1) = U*U*Xdp*sin(2*d_v) - Eqprime*U/r.Xdprime*sin(d_v);
    I(11,2) = U/r.Xdprime*cos(d_v);
    I(11,10) = -I(11,1) - U*r.Uc*e.Y12*sin(V-e.A12);
    I(11,11) = Eqprime/r.Xdprime*cos(d_v) - 2*U/r.Xdprime*cos(d_v)*cos(d_v) - 2*U/r.Xd*sin(d_v)*sin(d_v) + r.Uc*e.Y12*cos(V-e.A12) - 2*U*e.Y11*cos(e.A11);
    return I;
}

//-----------------------------------------------------------------------

CalculusTrapeze::CalculusTrapeze(const Params& p) : Calculus(p), pimpl(new Impl(p)) {}

void CalculusTrapeze::emit_x(double t, size_t row, size_t n_steps)
{
    AnswerItem i;
    i.time=t, i.row=row, i.set_no=1, i.n_steps=n_steps;
    i.delta=pimpl->x(1), i.omega=pimpl->x(0), i.Eqe=pimpl->x(3);
    i.Eqprime=pimpl->x(2), i.V=pimpl->x(10), i.U=pimpl->x(11);
    emit a_step_done(i);
}

void CalculusTrapeze::solve_newton(double t, size_t& n_steps)
{
    // recalculating and reassigning 'e' at every step is essential to get some hysteresis behaviour
    pimpl->e = recalculate_equiv_params(t, pimpl->x(11), pimpl->p, pimpl->e);
    // as our newton implementation is able to return a new step value, omit some arithmetics
    pimpl->x = solve_newton_impl(pimpl, n_steps);
    pimpl->dirty_hack();
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct CalculusSequensive::Impl
{
    typedef X_N<5> X;
    typedef IMatrix_N<5> IMatrix;
    typedef AMatrix_N<5> AMatrix;
    typedef X_N<9> Xj;

    X calculate_W(const X& x_i_1) const;
    IMatrix calculate_I(const X& x_i_1) const;

    Xj calculate_xj(const X& x_n, const X& x_n_1, const Xj& xj_n_1) const;
    void dirty_hack();

    Impl(const Params& _p) : p(_p)
    {
        x(0) = 0.; // delta_omega_0
        x(1) = p.start.Delta0;
        x(2) = p.start.Eqprime0;
        x(3) = p.start.V0;
        x(4) = p.start.U0;

        xj(0) = p.start.U0;
        xj(1) = p.reg.Ur0 - p.start.U0;
        xj(3) = p.start.V0;
        xj(7) = p.start.Eqe0;
        xj(8) = p.start.Eqe0;
    }

    X x;
    Xj xj;
    Equiv e;
    Params p;
    size_t n_steps;
};

void CalculusSequensive::Impl::dirty_hack()
{
    if(not p.dirty_hack)
        return;
    xj(8) = qMin(xj(8), 2*p.reg.Eqenom);
    xj(8) = qMax(xj(8), 0.);
}

#define _exp(T) std::exp(-dt/(T))
#define n_exp(T) (1. - _exp(T))

CalculusSequensive::Impl::Xj CalculusSequensive::Impl::calculate_xj(const X& x_n, const X& x_n_1, const Xj& xj_n_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;

    const double X1_k = xj_n_1(0), X2_k = xj_n_1(1), X4_k = xj_n_1(2), X5_k = xj_n_1(3),
                 X6_k = xj_n_1(4), X7_k = xj_n_1(5), X9_k = xj_n_1(6), X0_k = xj_n_1(7),
                 Eqe_k = xj_n_1(8);

    const double V_k = x_n_1(3);        const double V = x_n(3);
    const double U_k = x_n_1(4);        const double U = x_n(4);

    const double X1 = X1_k * _exp(r.Tu) + U - U_k + n_exp(r.Tu) * (U_k - r.Tu/dt * (U - U_k));
    const double X2 = r.Ur0 - X1;
    const double X3 = r.K0u * X2;
    const double X4 = X4_k * _exp(r.Tg) + r.K1u/dt*(X2 - X2_k) * n_exp(r.Tg);
    const double X5 = X5_k * _exp(r.Tphi) + V - V_k + n_exp(r.Tphi) * (V_k - r.Tphi/dt * (V - V_k));
    const double X6 = X6_k * _exp(r.Tf) +  1./dt*(X5 - X5_k) * n_exp(r.Tf);
    const double X7 = X7_k * _exp(r.T ) + r.T/dt*(X6 - X6_k) * n_exp(r.T );
    const double X8 = r.K0f * X7;
    const double X9 = X9_k * _exp(r.Tg) + r.K1f/dt*(X7 - X7_k) * n_exp(r.Tg);
    const double X0 = X3 + X4 + X8 + X9 + p.start.Eqe0 + e.Upphi;
    const double Eqe = Eqe_k* _exp(r.Te) + X0 - X0_k + n_exp(r.Te) * (X0_k - r.Te/dt * (X0 - X0_k));

    Xj xj_n;
    xj_n(0) = X1;
    xj_n(1) = X2;
    xj_n(2) = X4;
    xj_n(3) = X5;
    xj_n(4) = X6;
    xj_n(5) = X7;
    xj_n(6) = X9;
    xj_n(7) = X0;
    xj_n(8) = Eqe;

    return xj_n;
}

CalculusSequensive::Impl::X CalculusSequensive::Impl::calculate_W(const X& x_i_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;
    const X& x_k_1 = x;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

    const double domega_k = x_k_1(0);   const double domega = x_i_1(0);
    const double delta_k = x_k_1(1);    const double delta = x_i_1(1);
    const double Eqprime_k = x_k_1(2);  const double Eqprime = x_i_1(2);
                                        const double V = x_i_1(3);
                                        const double U = x_i_1(4);

    const double d_v = delta-V;
    const double s_d_v = sin(d_v);
    const double c_d_v = cos(d_v);

    const Xj xj_n = calculate_xj(x_i_1, x_k_1, xj);
    const double Eqe = xj_n(8);

    X W;
    W(0) = delta - delta_k - dt*domega;
    W(1) = domega - domega_k - (dt/r.Ty)*r.omega_nom*(r.Pt0 - e.Pd/r.omega_nom*domega - Eqprime*U/r.Xdprime*s_d_v + U*U*Xdp*s_d_v*c_d_v);
    W(2) = Eqprime - Eqprime_k - (dt/r.Td0)*(Eqe - Eqprime*r.Xd/r.Xdprime + U*r.Xd*Xdp*c_d_v);
    W(3) = Eqprime*U/r.Xdprime*s_d_v - U*U*Xdp*s_d_v*c_d_v - U*U*e.Y11*sin(e.A11) - U*r.Uc*e.Y12*sin(V-e.A12);
    W(4) = Eqprime*U/r.Xdprime*c_d_v - U*U*(c_d_v*c_d_v/r.Xdprime + s_d_v*s_d_v/r.Xd) - U*U*e.Y11*cos(e.A11) + U*r.Uc*e.Y12*cos(V-e.A12);

    return W;
}

CalculusSequensive::Impl::IMatrix CalculusSequensive::Impl::calculate_I(const X& x_i_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

    const double delta = x_i_1(1);
    const double Eqprime = x_i_1(2);
    const double V = x_i_1(3);
    const double U = x_i_1(4);
    const double d_v = delta-V;
    const double s_d_v = sin(d_v);
    const double c_d_v = cos(d_v);

    const double dEqe_Te = 1. - n_exp(r.Te)*r.Te/dt;

    const double dEqe_dV = dEqe_Te * (1. - r.Tphi/dt*n_exp(r.Tphi)) * n_exp(r.Tf) * n_exp(r.T)
                                   * r.T/dt/dt * ( r.K0f + r.K1f/dt * n_exp(r.Tg) );

    const double dEqe_dU = dEqe_Te * (1. - r.Tu/dt*n_exp(r.Tu)) * ( -r.K0u - r.K1u/dt * n_exp(r.Tg) );

    IMatrix I;

    I(0,0) = -dt;
    I(0,1) = 1.;

    I(1,0) = (dt/r.Ty)*e.Pd + 1.;
    I(1,1) = (dt/r.Ty)*r.omega_nom*(Eqprime*U/r.Xdprime*c_d_v - U*U*Xdp*cos(2*d_v));
    I(1,2) = (dt/r.Ty)*r.omega_nom*U/r.Xdprime*s_d_v;
    I(1,3) = -I(1,1);
    I(1,4) = (dt/r.Ty)*r.omega_nom*(Eqprime/r.Xdprime*s_d_v - 2*U*Xdp*s_d_v*c_d_v);

    I(2,1) = (dt/r.Td0)*U*r.Xd*Xdp*s_d_v;
    I(2,2) = (dt/r.Td0)*r.Xd/r.Xdprime + 1.;
    I(2,3) = -I(2,1) - (dt/r.Td0)*dEqe_dV;
    I(2,4) = -(dt/r.Td0)*(r.Xd*Xdp*c_d_v + dEqe_dU);

    I(3,1) = Eqprime*U/r.Xdprime*c_d_v - U*U*Xdp*cos(2*d_v);
    I(3,2) = U/r.Xdprime*s_d_v;
    I(3,3) = -I(3,1) - U*r.Uc*e.Y12*cos(V-e.A12);
    I(3,4) = Eqprime/r.Xdprime*s_d_v - 2.*U*Xdp*s_d_v*c_d_v - 2.*U*e.Y11*sin(e.A11) - r.Uc*e.Y12*sin(V-e.A12);

    I(4,1) = U*U*Xdp*2*s_d_v*c_d_v - Eqprime*U/r.Xdprime*s_d_v;
    I(4,2) = U/r.Xdprime*c_d_v;
    I(4,3) = -I(4,1) - U*r.Uc*e.Y12*sin(V-e.A12);
    I(4,4) = Eqprime/r.Xdprime*c_d_v - 2*U*(c_d_v*c_d_v/r.Xdprime + s_d_v*s_d_v/r.Xd) - 2*U*e.Y11*cos(e.A11) + r.Uc*e.Y12*cos(V-e.A12);
    return I;
}

//-----------------------------------------------------------------------

CalculusSequensive::CalculusSequensive(const Params& p) : Calculus(p), pimpl(new Impl(p)) {}

void CalculusSequensive::emit_x(double t, size_t row, size_t n_steps)
{
    AnswerItem i;
    i.time=t, i.row=row, i.set_no=2, i.n_steps=n_steps;
    i.delta=pimpl->x(1), i.omega=pimpl->x(0), i.Eqe=pimpl->xj(8);
    i.Eqprime=pimpl->x(2), i.V=pimpl->x(3), i.U=pimpl->x(4);
    emit a_step_done(i);
}

void CalculusSequensive::solve_newton(double t, size_t& n_steps)
{
    // recalculating and reassigning 'e' at every step is essential to get some hysteresis behaviour
    pimpl->e = recalculate_equiv_params(t, pimpl->x(4), pimpl->p, pimpl->e);
    // as our newton implementation is able to return a new step value, omit some arithmetics
    Impl::X x = solve_newton_impl(pimpl, n_steps);
    pimpl->xj = pimpl->calculate_xj(x, pimpl->x, pimpl->xj);
    pimpl->dirty_hack();
    pimpl->x = x;
}

void CalculusSequensive::set_X(double x1, double x2, double x4, double x5, double x6, double x7, double x9, double x0)
{
    pimpl->xj(0) = x1;
    pimpl->xj(1) = x2;
    pimpl->xj(2) = x4;
    pimpl->xj(3) = x5;
    pimpl->xj(4) = x6;
    pimpl->xj(5) = x7;
    pimpl->xj(6) = x9;
    pimpl->xj(7) = x0;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct CalculusParallel::Impl
{
    typedef X_N<5> X;
    typedef IMatrix_N<5> IMatrix;
    typedef AMatrix_N<5> AMatrix;
    typedef X_N<12> Xj;

    X calculate_W(const X& ) const;
    IMatrix calculate_I(const X& ) const;

    Xj calculate_xj(const X& x_n, const X& x_n_1, const Xj& xj_n_1) const;
    void dirty_hack();

    Impl(const Params& _p) : p(_p)
    {
        x(0) = 0.; // delta_omega_0
        x(1) = p.start.Delta0;
        x(2) = p.start.Eqprime0;
        x(3) = p.start.V0;
        x(4) = p.start.U0;

        xj(0) = p.start.U0;
        xj(1) = p.reg.Ur0 - xj(0);
        xj(2) = p.reg.K1U * xj(1);
        xj(3) = p.reg.K2U * xj(1);
        xj(4) = p.reg.K1V * p.start.V0;
        xj(5) = p.reg.K2V * p.start.V0;
        xj(6) = p.reg.K3V * p.start.V0;
        xj(7) = p.reg.K4V * p.start.V0;
        xj(8) = p.reg.K5V * p.start.V0;
        xj(9) = p.start.Eqe0;
        xj(10) = p.start.Eqe0;
    }

    X x;
    Xj xj;
    Equiv e;
    Params p;
    size_t n_steps;
};

void CalculusParallel::Impl::dirty_hack()
{
    if(not p.dirty_hack)
        return;
}

CalculusParallel::Impl::Xj CalculusParallel::Impl::calculate_xj(const X& x_n, const X& x_n_1, const Xj& xj_n_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;

    const double X1_k = xj_n_1(0), X2_k = xj_n_1(1), X3_k = xj_n_1(2), X4_k = xj_n_1(3), X5_k = xj_n_1(4),
                 X6_k = xj_n_1(5), X7_k = xj_n_1(6), X8_k = xj_n_1(7), X9_k = xj_n_1(8), X10_k = xj_n_1(9),
                 X11_k = xj_n_1(10), Upphi_k = xj_n_1(11);

    const double V_k = x_n_1(3);        const double V = x_n(3);
    const double U_k = x_n_1(4);        const double U = x_n(4);

    const double X1 = X1_k * _exp(r.Tu) + U - U_k + n_exp(r.Tu) * (U_k - r.Tu/dt * (U - U_k));
    const double X2 = r.Ur0 - X1;
    const double X3 = X3_k * _exp(r.Te) + r.K1U * (X2 - X2_k + n_exp(r.Te) * (X2_k - r.Te/dt * (X2 - X2_k)));
    const double X4 = X4_k * _exp(r.Tg) + r.K2U * (X2 - X2_k + n_exp(r.Tg) * (X2_k - r.Tg/dt * (X2 - X2_k)));
    const double X5 = X5_k * _exp(r.Tphi) + r.K1V * (V - V_k + n_exp(r.Tphi) * (V_k - r.Tphi/dt * (V - V_k)));
    const double X6 = X6_k * _exp(r.Tf  ) + r.K2V * (V - V_k + n_exp(r.Tf  ) * (V_k - r.Tf  /dt * (V - V_k)));
    const double X7 = X7_k * _exp(r.T   ) + r.K3V * (V - V_k + n_exp(r.T   ) * (V_k - r.T   /dt * (V - V_k)));
    const double X8 = X8_k * _exp(r.Tg  ) + r.K4V * (V - V_k + n_exp(r.Tg  ) * (V_k - r.Tg  /dt * (V - V_k)));
    const double X9 = X9_k * _exp(r.Te  ) + r.K5V * (V - V_k + n_exp(r.Te  ) * (V_k - r.Te  /dt * (V - V_k)));
    const double X10 = X10_k * _exp(r.Te) + p.start.Eqe0 * n_exp(r.Te);
    const double X11 = X11_k * _exp(r.Te) + e.Upphi - Upphi_k + n_exp(r.Te) * (Upphi_k - r.Te/dt * (e.Upphi - Upphi_k));
    //const double Eqe = X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11;

    Xj xj_n;
    xj_n(0) = X1;
    xj_n(1) = X2;
    xj_n(2) = X3;
    xj_n(3) = X4;
    xj_n(4) = X5;
    xj_n(5) = X6;
    xj_n(6) = X7;
    xj_n(7) = X8;
    xj_n(8) = X9;
    xj_n(9) = X10;
    xj_n(10) = X11;
    xj_n(11) = e.Upphi;

    return xj_n;
}

CalculusParallel::Impl::X CalculusParallel::Impl::calculate_W(const X& x_i_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;
    const X& x_k_1 = x;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

    const double domega_k = x_k_1(0);   const double domega = x_i_1(0);
    const double delta_k = x_k_1(1);    const double delta = x_i_1(1);
    const double Eqprime_k = x_k_1(2);  const double Eqprime = x_i_1(2);
                                        const double V = x_i_1(3);
                                        const double U = x_i_1(4);

    const double d_v = delta-V;
    const double s_d_v = sin(d_v);
    const double c_d_v = cos(d_v);

    const Xj xj_n = calculate_xj(x_i_1, x_k_1, xj);
    //const double Eqe = X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11;
    const double Eqe = xj_n(2) + xj_n(3) + xj_n(4) + xj_n(5) + xj_n(6) + xj_n(7) + xj_n(8) + xj_n(9) + xj_n(10);
    const double Eqe_dirty_hacked = p.dirty_hack ? qMax(qMin(Eqe, 2*p.reg.Eqenom), 0.) : Eqe;

    X W;
    W(0) = delta - delta_k - dt*domega;
    W(1) = domega - domega_k - (dt/r.Ty)*r.omega_nom*(r.Pt0 - e.Pd/r.omega_nom*domega - Eqprime*U/r.Xdprime*s_d_v + U*U*Xdp*s_d_v*c_d_v);
    W(2) = Eqprime - Eqprime_k - (dt/r.Td0)*(Eqe_dirty_hacked - Eqprime*r.Xd/r.Xdprime + U*r.Xd*Xdp*c_d_v);
    W(3) = Eqprime*U/r.Xdprime*s_d_v - U*U*Xdp*s_d_v*c_d_v - U*U*e.Y11*sin(e.A11) - U*r.Uc*e.Y12*sin(V-e.A12);
    W(4) = Eqprime*U/r.Xdprime*c_d_v - U*U*(c_d_v*c_d_v/r.Xdprime + s_d_v*s_d_v/r.Xd) - U*U*e.Y11*cos(e.A11) + U*r.Uc*e.Y12*cos(V-e.A12);

    return W;
}

CalculusParallel::Impl::IMatrix CalculusParallel::Impl::calculate_I(const X& x_i_1) const
{
    const double dt =p.dt;
    const Params::Consts& r = p.reg;

    const double Xdp = (r.Xd-r.Xdprime)/r.Xd/r.Xdprime;

    const double delta = x_i_1(1);
    const double Eqprime = x_i_1(2);
    const double V = x_i_1(3);
    const double U = x_i_1(4);
    const double d_v = delta-V;
    const double s_d_v = sin(d_v);
    const double c_d_v = cos(d_v);

    const double dEqe_dV = r.K1V * (1 - r.Tphi/dt * n_exp(r.Tphi)) +
                           r.K2V * (1 - r.Tf  /dt * n_exp(r.Tf  )) +
                           r.K3V * (1 - r.T   /dt * n_exp(r.T   )) +
                           r.K4V * (1 - r.Tg  /dt * n_exp(r.Tg  )) +
                           r.K5V * (1 - r.Te  /dt * n_exp(r.Te  ));

    const double dEqe_dU = (r.Tu/dt * n_exp(r.Tu) - 1) * (r.K1U * (1 - r.Te/dt * n_exp(r.Te)) + r.K2U * (1 - r.Tg/dt * n_exp(r.Tg)));

    IMatrix I;

    I(0,0) = -dt;
    I(0,1) = 1.;

    I(1,0) = (dt/r.Ty)*e.Pd + 1.;
    I(1,1) = (dt/r.Ty)*r.omega_nom*(Eqprime*U/r.Xdprime*c_d_v - U*U*Xdp*cos(2*d_v));
    I(1,2) = (dt/r.Ty)*r.omega_nom*U/r.Xdprime*s_d_v;
    I(1,3) = -I(1,1);
    I(1,4) = (dt/r.Ty)*r.omega_nom*(Eqprime/r.Xdprime*s_d_v - 2*U*Xdp*s_d_v*c_d_v);

    I(2,1) = (dt/r.Td0)*U*r.Xd*Xdp*s_d_v;
    I(2,2) = (dt/r.Td0)*r.Xd/r.Xdprime + 1.;
    I(2,3) = -I(2,1) - (dt/r.Td0)*dEqe_dV;
    I(2,4) = -(dt/r.Td0)*(r.Xd*Xdp*c_d_v + dEqe_dU);

    I(3,1) = Eqprime*U/r.Xdprime*c_d_v - U*U*Xdp*cos(2*d_v);
    I(3,2) = U/r.Xdprime*s_d_v;
    I(3,3) = -I(3,1) - U*r.Uc*e.Y12*cos(V-e.A12);
    I(3,4) = Eqprime/r.Xdprime*s_d_v - 2.*U*Xdp*s_d_v*c_d_v - 2.*U*e.Y11*sin(e.A11) - r.Uc*e.Y12*sin(V-e.A12);

    I(4,1) = U*U*Xdp*2*s_d_v*c_d_v - Eqprime*U/r.Xdprime*s_d_v;
    I(4,2) = U/r.Xdprime*c_d_v;
    I(4,3) = -I(4,1) - U*r.Uc*e.Y12*sin(V-e.A12);
    I(4,4) = Eqprime/r.Xdprime*c_d_v - 2*U*(c_d_v*c_d_v/r.Xdprime + s_d_v*s_d_v/r.Xd) - 2*U*e.Y11*cos(e.A11) + r.Uc*e.Y12*cos(V-e.A12);
    return I;
}

//-----------------------------------------------------------------------

CalculusParallel::CalculusParallel(const Params& p) : Calculus(p), pimpl(new CalculusParallel::Impl(p)) {}

void CalculusParallel::emit_x(double t, size_t row, size_t n_steps)
{
    //const double Eqe = X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11;
    const double Eqe = pimpl->xj(2) + pimpl->xj(3) + pimpl->xj(4) + pimpl->xj(5) + pimpl->xj(6) + pimpl->xj(7) + pimpl->xj(8) + pimpl->xj(9) + pimpl->xj(10);
    const double Eqe_dirty_hacked = p.dirty_hack ? qMax(qMin(Eqe, 2*p.reg.Eqenom), 0.) : Eqe;

    AnswerItem i;
    i.time=t, i.row=row, i.set_no=3, i.n_steps=n_steps;
    i.delta=pimpl->x(1), i.omega=pimpl->x(0), i.Eqe=Eqe_dirty_hacked;
    i.Eqprime=pimpl->x(2), i.V=pimpl->x(3), i.U=pimpl->x(4);
    emit a_step_done(i);
}

void CalculusParallel::solve_newton(double t, size_t& n_steps)
{
    // recalculating and reassigning 'e' at every step is essential to get some hysteresis behaviour
    pimpl->e = recalculate_equiv_params(t, pimpl->x(4), pimpl->p, pimpl->e);
    // as our newton implementation is able to return a new step value, omit some arithmetics
    Impl::X x = solve_newton_impl(pimpl, n_steps);
    pimpl->xj = pimpl->calculate_xj(x, pimpl->x, pimpl->xj);
    pimpl->dirty_hack();
    pimpl->x = x;
}

