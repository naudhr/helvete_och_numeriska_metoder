#include "StartingPage.h"
#include "params.h"
#include "Calculus.h"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QProgressBar>
#include <QPushButton>
#include <QLineEdit>
#include <QThread>
#include <qwt6/qwt_plot.h>
#include <qwt6/qwt_plot_curve.h>
#include <qwt6/qwt_legend.h>
#include <qwt6/qwt_plot_zoomer.h>
#include <qwt6/qwt_plot_grid.h>
//#include <qwt6/qwt_plot_magnifier.h>


StartingPage::StartingPage() : QWidget()
{
    auto reg = new QHBoxLayout;
    auto eqv = new QHBoxLayout;
    auto srt = new QHBoxLayout;
    auto cpt = new QHBoxLayout;

    progress_bar_eiler = new QProgressBar(this);
    progress_bar_eiler->setFormat("Eiler: %v");
    progress_bar_eiler->setMinimum(0);
    progress_bar_eiler->setMaximum(1);
    progress_bar_eiler->setValue(0);
    progress_bar_trapeze = new QProgressBar(this);
    progress_bar_trapeze->setFormat("Trapeze: %v");
    progress_bar_trapeze->setMinimum(0);
    progress_bar_trapeze->setMaximum(1);
    progress_bar_trapeze->setValue(0);
    auto prb = new QHBoxLayout;
    prb->addWidget(progress_bar_eiler);
    prb->addWidget(progress_bar_trapeze);

    start_button = new QPushButton("Start",this);

    plot = new QwtPlot(this);
    curve_eiler_Delta = new QwtPlotCurve("Delta Eiler");
    curve_eiler_Delta->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Delta->setPen(QPen(QColor(160,0,210)));
    curve_eiler_Delta->setVisible(false);
    curve_eiler_Delta->attach(plot);
    curve_eiler_Omega = new QwtPlotCurve("Omega Eiler");
    curve_eiler_Omega->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Omega->setPen(QPen(QColor(210,0,0)));
    curve_eiler_Omega->setVisible(false);
    curve_eiler_Omega->attach(plot);
    curve_eiler_Eqe = new QwtPlotCurve("Eqe Eiler");
    curve_eiler_Eqe->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Eqe->setPen(QPen(QColor(5,205,0)));
    curve_eiler_Eqe->setVisible(false);
    curve_eiler_Eqe->attach(plot);
    curve_eiler_Eqprime = new QwtPlotCurve("Eqprime Eiler");
    curve_eiler_Eqprime->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Eqprime->setPen(QPen(QColor(210,120,0)));
    curve_eiler_Eqprime->setVisible(false);
    curve_eiler_Eqprime->attach(plot);
    curve_eiler_U = new QwtPlotCurve("U Eiler");
    curve_eiler_U->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_U->setPen(QPen(Qt::black));
    curve_eiler_U->setVisible(false);
    curve_eiler_U->attach(plot);

    curve_trapeze_Delta = new QwtPlotCurve("Delta Trapeze");
    curve_trapeze_Delta->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Delta->setPen(QPen(QColor(210,50,255)));
    curve_trapeze_Delta->setVisible(false);
    curve_trapeze_Delta->attach(plot);
    curve_trapeze_Omega = new QwtPlotCurve("Omega Trapeze");
    curve_trapeze_Omega->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Omega->setPen(QPen(QColor(255,45,45)));
    curve_trapeze_Omega->setVisible(false);
    curve_trapeze_Omega->attach(plot);
    curve_trapeze_Eqe = new QwtPlotCurve("Eqe Trapeze");
    curve_trapeze_Eqe->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Eqe->setPen(QPen(QColor(50,255,50)));
    curve_trapeze_Eqe->setVisible(false);
    curve_trapeze_Eqe->attach(plot);
    curve_trapeze_Eqprime = new QwtPlotCurve("Eqprime Trapeze");
    curve_trapeze_Eqprime->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Eqprime->setPen(QPen(QColor(255,160,50)));
    curve_trapeze_Eqprime->setVisible(false);
    curve_trapeze_Eqprime->attach(plot);

    auto legend = new QwtLegend;
    legend->setItemMode(QwtLegend::ReadOnlyItem);
    plot->insertLegend(legend,QwtPlot::TopLegend);

    zoom = new QwtPlotZoomer(plot->canvas());
    zoom->setRubberBandPen(QPen(Qt::white));

    auto grid = new QwtPlotGrid;
    grid->enableXMin(true);
    grid->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid->setMinPen(QPen(Qt::gray,0,Qt::DotLine));
    grid->attach(plot);

    auto l = new QVBoxLayout(this);
    l->addLayout(reg);
    l->addLayout(eqv);
    l->addLayout(srt);
    l->addLayout(cpt);
    l->addLayout(prb);
    l->addWidget(plot,2);
    l->addWidget(start_button);

    plot->setVisible(true);
    progress_bar_eiler->setVisible(true);
    progress_bar_trapeze->setVisible(true);
    connect(start_button, SIGNAL(clicked()), SLOT(start()));
    //start();
}

void StartingPage::start()
{
    start_button->setEnabled(false);

    auto p = collect_params();
    auto n = (p.Tstop - p.Tstart)/p.dt;

    plot->setAxisScale(QwtPlot::xBottom,p.Tstart-p.dt,p.Tstop+p.dt);

    progress_bar_eiler->setMaximum(n);
    progress_bar_eiler->setValue(0);
    progress_bar_eiler->setVisible(true);

    CalculusEiler e;//(p);
    connect(&e, SIGNAL(a_step_done()), SLOT(increment_eiler()));
    //connect(&e, SIGNAL(answer(Answer)), this, SLOT(plot_answer(Answer)));
    plot_answer_eiler(e.doWork(p));
/*
    auto e = new CalculusEiler(p);
    auto t = new QThread(this);
    connect(t, SIGNAL(started()), e, SLOT(doWork()), Qt::QueuedConnection);
    connect(t, SIGNAL(finished()), e, SLOT(deleteLater()), Qt::QueuedConnection);
    connect(e, SIGNAL(a_step_done()), this, SLOT(increment_eiler()), Qt::QueuedConnection);
    connect(e, SIGNAL(answer(Answer)), this, SLOT(plot_answer(Answer)), Qt::QueuedConnection);
    e->moveToThread(t);
    t->start();
*/
    progress_bar_trapeze->setMaximum(n);
    progress_bar_trapeze->setValue(0);
    progress_bar_trapeze->setVisible(true);

    CalculusTrapeze t;
    connect(&t, SIGNAL(a_step_done()), SLOT(increment_trapeze()));
    plot_answer_trapeze(t.doWork(p));

    plot->setAxisAutoScale(QwtPlot::xBottom);
    plot->setAxisAutoScale(QwtPlot::yLeft);
    plot->setVisible(true);
    plot->replot();
    zoom->setZoomBase();
    start_button->setEnabled(true);
}

bool StartingPage::if_all_fields_are_filled()
{
    return true;
}

void StartingPage::change_params_reactor()
{
    start_button->setEnabled( if_all_fields_are_filled() );
}

Params StartingPage::collect_params()
{
    Params p = {};
    p.reg.K0f=10.;
    p.reg.K1f=0.5;
    p.reg.Ur0=1.05;
    p.reg.K0u=50.;
    p.reg.K1u=3.6;
    p.repl.Y11=0.496;
    p.repl.Y12=0.532;
    p.repl.A11=0.958;
    p.repl.A12=0.715;
    p.repl.Y11em=0.35;
    p.repl.Y12em=0.45;
    p.repl.A11em=0.9;
    p.repl.A12em=0.3;
    p.start.Delta0=1.365;
    p.start.Eqe0=1.87;
    p.start.Eqprime0=1.1;
    p.start.U0=1.059;
    p.start.V0=0.597;
    p.Tstart=0.;
    p.Tstop=3.;
    p.dt=0.005;
    p.eps=0.005;
    p.max_iterations=10;
    return p;
}

void StartingPage::increment_eiler()
{
    progress_bar_eiler->setValue(progress_bar_eiler->value()+1);
}

void StartingPage::increment_trapeze()
{
    progress_bar_trapeze->setValue(progress_bar_trapeze->value()+1);
}

void StartingPage::plot_answer_eiler(const QVector<AnswerItem> ans)
{
    qDebug()<<__FUNCTION__<<ans.size();

    QVector<QPointF> data_delta(ans.size());
    QVector<QPointF> data_omega(ans.size());
    QVector<QPointF> data_eqe(ans.size());
    QVector<QPointF> data_eqp(ans.size());
    QVector<QPointF> data_u(ans.size());
    double max_delta = 0., max_omega = 0., max_Eqprime = 0., max_Eqe = 0.;
    double min_delta = 0., min_omega = 0., min_Eqprime = 0., min_Eqe = 0.;
    for(int i=0; i<ans.size(); i++)
    {
        const auto a = ans.at(i);
        data_delta[i] = QPointF(a.time, a.delta);
        data_omega[i] = QPointF(a.time, a.omega);
        data_eqp[i] = QPointF(a.time, a.Eqprime);
        data_eqe[i] = QPointF(a.time, a.Eqe);
        data_u[i] = QPointF(a.time, a.U);
        max_delta = std::max(max_delta, a.delta);
        max_omega = std::max(max_omega, a.omega);
        max_Eqprime = std::max(max_Eqprime, a.Eqprime);
        max_Eqe = std::max(max_Eqe, a.Eqe);
        min_delta = std::min(min_delta, a.delta);
        min_omega = std::min(min_omega, a.omega);
        min_Eqprime = std::min(min_Eqprime, a.Eqprime);
        min_Eqe = std::min(min_Eqe, a.Eqe);
    }
    curve_eiler_U->setSamples(data_u);
    curve_eiler_Delta->setSamples(data_delta);
    curve_eiler_Omega->setSamples(data_omega);
    curve_eiler_Eqprime->setSamples(data_eqp);
    curve_eiler_Eqe->setSamples(data_eqe);
    curve_eiler_U->setVisible(true);
    curve_eiler_Delta->setVisible(true);
    curve_eiler_Omega->setVisible(true);
    curve_eiler_Eqprime->setVisible(true);
    curve_eiler_Eqe->setVisible(true);

    double min = std::min(std::min(std::min(min_delta,min_omega),min_Eqprime),min_Eqe);
    double max = std::max(std::max(std::max(max_delta,max_omega),max_Eqprime),max_Eqe);
    plot->setAxisScale(QwtPlot::yLeft,min,max);
    qDebug()<<"plot->setVisible(true): max"<<max<<"max_delta"<<max_delta<<"max_omega"<<max_omega<<"max_Eqprime"<<max_Eqprime<<"max_Eqe"<<max_Eqe;
    qDebug()<<"plot->setVisible(true): min"<<min<<"min_delta"<<min_delta<<"min_omega"<<min_omega<<"min_Eqprime"<<min_Eqprime<<"min_Eqe"<<min_Eqe;
}

void StartingPage::plot_answer_trapeze(const QVector<AnswerItem> ans)
{
    qDebug()<<__FUNCTION__<<ans.size();

    QVector<QPointF> data_delta(ans.size());
    QVector<QPointF> data_omega(ans.size());
    QVector<QPointF> data_eqe(ans.size());
    QVector<QPointF> data_eqp(ans.size());
    double max_delta = 0., max_omega = 0., max_Eqprime = 0., max_Eqe = 0.;
    double min_delta = 0., min_omega = 0., min_Eqprime = 0., min_Eqe = 0.;
    for(int i=0; i<ans.size(); i++)
    {
        const auto a = ans.at(i);
        data_delta[i] = QPointF(a.time, a.delta);
        data_omega[i] = QPointF(a.time, a.omega);
        data_eqp[i] = QPointF(a.time, a.Eqprime);
        data_eqe[i] = QPointF(a.time, a.Eqe);
        max_delta = std::max(max_delta, a.delta);
        max_omega = std::max(max_omega, a.omega);
        max_Eqprime = std::max(max_Eqprime, a.Eqprime);
        max_Eqe = std::max(max_Eqe, a.Eqe);
        min_delta = std::min(min_delta, a.delta);
        min_omega = std::min(min_omega, a.omega);
        min_Eqprime = std::min(min_Eqprime, a.Eqprime);
        min_Eqe = std::min(min_Eqe, a.Eqe);
    }
    curve_trapeze_Delta->setSamples(data_delta);
    curve_trapeze_Omega->setSamples(data_omega);
    curve_trapeze_Eqprime->setSamples(data_eqp);
    curve_trapeze_Eqe->setSamples(data_eqe);
    curve_trapeze_Delta->setVisible(true);
    curve_trapeze_Omega->setVisible(true);
    curve_trapeze_Eqprime->setVisible(true);
    curve_trapeze_Eqe->setVisible(true);

    double min = std::min(std::min(std::min(min_delta,min_omega),min_Eqprime),min_Eqe);
    double max = std::max(std::max(std::max(max_delta,max_omega),max_Eqprime),max_Eqe);
    //plot->setAxisScale(QwtPlot::yLeft,min,max);
    qDebug()<<"plot->setVisible(true): max"<<max<<"max_delta"<<max_delta<<"max_omega"<<max_omega<<"max_Eqprime"<<max_Eqprime<<"max_Eqe"<<max_Eqe;
    qDebug()<<"plot->setVisible(true): min"<<min<<"min_delta"<<min_delta<<"min_omega"<<min_omega<<"min_Eqprime"<<min_Eqprime<<"min_Eqe"<<min_Eqe;
}

