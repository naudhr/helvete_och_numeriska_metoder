#include "StartingPage.h"
#include "params.h"
#include "Calculus.h"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QProgressBar>
#include <QPushButton>
#include <QLineEdit>
#include <QCheckBox>
#include <QLabel>
#include <QDoubleValidator>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_grid.h>


StartingPage::StartingPage() : QWidget()
{
    validator = new QDoubleValidator(this);

    QHBoxLayout* reg = new QHBoxLayout;
    reg->addWidget(new QLabel("Regulator parameters:",this));
    reg->addStretch(1);
    reg->addWidget(new QLabel("K0f:",this));
    reg->addWidget(K0f = new QLineEdit("10",this));
    K0f->setValidator(validator);
    reg->addStretch(1);
    reg->addWidget(new QLabel("K1f:",this));
    reg->addWidget(K1f = new QLineEdit("0.5",this));
    K1f->setValidator(validator);
    reg->addStretch(1);
    reg->addWidget(new QLabel("Ur0:",this));
    reg->addWidget(Ur0 = new QLineEdit("1.05",this));
    Ur0->setValidator(validator);
    reg->addStretch(1);
    reg->addWidget(new QLabel("K0u:",this));
    reg->addWidget(K0u = new QLineEdit("50",this));
    K0u->setValidator(validator);
    reg->addStretch(1);
    reg->addWidget(new QLabel("K1u:",this));
    reg->addWidget(K1u = new QLineEdit("3.6",this));
    K1u->setValidator(validator);

    QHBoxLayout* eqv = new QHBoxLayout;
    eqv->addWidget(new QLabel("Equivalent circuit parameters:",this));
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("Y11:",this));
    eqv->addWidget(Y11 = new QLineEdit("0.496",this));
    Y11->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("Y12:",this));
    eqv->addWidget(Y12 = new QLineEdit("0.532",this));
    Y12->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("A11 (degrees):",this));
    eqv->addWidget(A11 = new QLineEdit("0.958",this));
    A11->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("A12 (degrees):",this));
    eqv->addWidget(A12 = new QLineEdit("0.715",this));
    A12->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("Pd:",this));
    eqv->addWidget(Pd = new QLineEdit("5",this));
    Pd->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("Y11em:",this));
    eqv->addWidget(Y11em = new QLineEdit("0.35",this));
    Y11em->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("Y12em:",this));
    eqv->addWidget(Y12em = new QLineEdit("0.45",this));
    Y12em->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("A11em (degrees):",this));
    eqv->addWidget(A11em = new QLineEdit("0.9",this));
    A11em->setValidator(validator);
    eqv->addStretch(1);
    eqv->addWidget(new QLabel("A12em (degrees):",this));
    eqv->addWidget(A12em = new QLineEdit("0.3",this));
    A12em->setValidator(validator);

    QHBoxLayout* srt = new QHBoxLayout;
    srt->addWidget(new QLabel("Initial values:",this));
    srt->addStretch(1);
    srt->addWidget(new QLabel("Delta0:",this));
    srt->addWidget(Delta0 = new QLineEdit("1.365",this));
    Delta0->setValidator(validator);
    srt->addStretch(1);
    srt->addWidget(new QLabel("Eqe0:",this));
    srt->addWidget(Eqe0 = new QLineEdit("1.87",this));
    Eqe0->setValidator(validator);
    srt->addStretch(1);
    srt->addWidget(new QLabel("Eqprime0:",this));
    srt->addWidget(Eqprime0 = new QLineEdit("1.1",this));
    Eqprime0->setValidator(validator);
    srt->addStretch(1);
    srt->addWidget(new QLabel("U0:",this));
    srt->addWidget(U0 = new QLineEdit("1.059",this));
    U0->setValidator(validator);
    srt->addStretch(1);
    srt->addWidget(new QLabel("V0:",this));
    srt->addWidget(V0 = new QLineEdit("0.597",this));
    V0->setValidator(validator);

    QHBoxLayout* cpt = new QHBoxLayout;
    cpt->addWidget(new QLabel("Computational parameters:",this));
    cpt->addStretch(1);
    cpt->addWidget(new QLabel("Tstart:",this));
    cpt->addWidget(Tstart = new QLineEdit("0",this));
    Tstart->setValidator(validator);
    cpt->addStretch(1);
    cpt->addWidget(new QLabel("Tstop:",this));
    cpt->addWidget(Tstop = new QLineEdit("3",this));
    Tstop->setValidator(validator);
    cpt->addStretch(1);
    cpt->addWidget(new QLabel("dt",this));
    cpt->addWidget(dt = new QLineEdit("0.005",this));
    dt->setValidator(validator);
    cpt->addStretch(1);
    cpt->addWidget(new QLabel("eps:",this));
    cpt->addWidget(eps = new QLineEdit("0.005",this));
    eps->setValidator(validator);
    cpt->addStretch(1);
    cpt->addWidget(new QLabel("max_iterations:",this));
    cpt->addWidget(max_iterations = new QLineEdit("10",this));
    max_iterations->setValidator(new QIntValidator(5,500,this));

    enable_eiler = new QCheckBox(this);
    enable_eiler->setChecked(false);
    enable_trapeze = new QCheckBox(this);
    enable_trapeze->setChecked(false);
    progress_bar_eiler = new QProgressBar(this);
    progress_bar_eiler->setFormat("Eiler");
    progress_bar_eiler->setMinimum(0);
    progress_bar_eiler->setMaximum(1);
    progress_bar_eiler->setValue(0);
    progress_bar_eiler->setEnabled(false);
    progress_bar_trapeze = new QProgressBar(this);
    progress_bar_trapeze->setFormat("Trapeze");
    progress_bar_trapeze->setMinimum(0);
    progress_bar_trapeze->setMaximum(1);
    progress_bar_trapeze->setValue(0);
    progress_bar_trapeze->setEnabled(false);
    QHBoxLayout* prb = new QHBoxLayout;
    prb->addWidget(enable_eiler);
    prb->addWidget(progress_bar_eiler);
    prb->addSpacing(30);
    prb->addWidget(enable_trapeze);
    prb->addWidget(progress_bar_trapeze);

    start_button = new QPushButton("Start",this);
    start_button->setEnabled(false);

    plot = new QwtPlot(this);
    curve_eiler_Delta = new QwtPlotCurve("Delta Eiler");
    curve_eiler_Delta->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Delta->setPen(QPen(QColor(160,0,210)));
    curve_eiler_Delta->setVisible(false);
    //curve_eiler_Delta->attach(plot);
    curve_eiler_Omega = new QwtPlotCurve("Omega Eiler");
    curve_eiler_Omega->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Omega->setPen(QPen(QColor(210,0,0)));
    curve_eiler_Omega->setVisible(false);
    //curve_eiler_Omega->attach(plot);
    curve_eiler_Eqe = new QwtPlotCurve("Eqe Eiler");
    curve_eiler_Eqe->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Eqe->setPen(QPen(QColor(5,205,0)));
    curve_eiler_Eqe->setVisible(false);
    //curve_eiler_Eqe->attach(plot);
    curve_eiler_Eqprime = new QwtPlotCurve("Eqprime Eiler");
    curve_eiler_Eqprime->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_Eqprime->setPen(QPen(QColor(210,120,0)));
    curve_eiler_Eqprime->setVisible(false);
    //curve_eiler_Eqprime->attach(plot);
    curve_eiler_U = new QwtPlotCurve("U Eiler");
    curve_eiler_U->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_eiler_U->setPen(QPen(Qt::black));
    curve_eiler_U->setVisible(false);
    //curve_eiler_U->attach(plot);

    curve_trapeze_Delta = new QwtPlotCurve("Delta Trapeze");
    curve_trapeze_Delta->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Delta->setPen(QPen(QColor(210,50,255)));
    curve_trapeze_Delta->setVisible(false);
    //curve_trapeze_Delta->attach(plot);
    curve_trapeze_Omega = new QwtPlotCurve("Omega Trapeze");
    curve_trapeze_Omega->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Omega->setPen(QPen(QColor(255,45,45)));
    curve_trapeze_Omega->setVisible(false);
    //curve_trapeze_Omega->attach(plot);
    curve_trapeze_Eqe = new QwtPlotCurve("Eqe Trapeze");
    curve_trapeze_Eqe->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Eqe->setPen(QPen(QColor(50,255,50)));
    curve_trapeze_Eqe->setVisible(false);
    //curve_trapeze_Eqe->attach(plot);
    curve_trapeze_Eqprime = new QwtPlotCurve("Eqprime Trapeze");
    curve_trapeze_Eqprime->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve_trapeze_Eqprime->setPen(QPen(QColor(255,160,50)));
    curve_trapeze_Eqprime->setVisible(false);
    //curve_trapeze_Eqprime->attach(plot);

    QwtLegend* legend = new QwtLegend;
    legend->setItemMode(QwtLegend::ReadOnlyItem);
    plot->insertLegend(legend,QwtPlot::TopLegend);

    zoom = new QwtPlotZoomer(plot->canvas());
    zoom->setRubberBandPen(QPen(Qt::white));

    QwtPlotGrid* grid = new QwtPlotGrid;
    grid->enableXMin(true);
    grid->setMajPen(QPen(Qt::black,0,Qt::DotLine));
    grid->setMinPen(QPen(Qt::gray,0,Qt::DotLine));
    grid->attach(plot);

    QVBoxLayout* l = new QVBoxLayout(this);
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

    connect(enable_eiler, SIGNAL(clicked()), SLOT(wanna_eiler()));
    connect(enable_trapeze, SIGNAL(clicked()), SLOT(wanna_trapeze()));
    //start_button->setFocus(Qt::TabFocusReason);
}

StartingPage::~StartingPage()
{
    if(not curve_eiler_U->plot() )  delete curve_eiler_U;
    if(not curve_eiler_Delta->plot() )  delete curve_eiler_Delta;
    if(not curve_eiler_Omega->plot() )  delete curve_eiler_Omega;
    if(not curve_eiler_Eqe->plot() )  delete curve_eiler_Eqe;
    if(not curve_eiler_Eqprime->plot() )  delete curve_eiler_Eqprime;
    if(not curve_trapeze_Delta->plot() )  delete curve_trapeze_Delta;
    if(not curve_trapeze_Omega->plot() )  delete curve_trapeze_Omega;
    if(not curve_trapeze_Eqe->plot() )  delete curve_trapeze_Eqe;
    if(not curve_trapeze_Eqprime->plot() )  delete curve_trapeze_Eqprime;
}

void StartingPage::wanna_eiler()
{
    progress_bar_eiler->setEnabled( enable_eiler->isChecked() );
    curve_eiler_U->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Delta->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Omega->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Eqe->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Eqprime->attach( enable_eiler->isChecked() ? plot : NULL );
    start_button->setEnabled( enable_eiler->isChecked() or enable_trapeze->isChecked() );
}

void StartingPage::wanna_trapeze()
{
    progress_bar_trapeze->setEnabled( enable_trapeze->isChecked() );
    curve_trapeze_Delta->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_Omega->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_Eqe->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_Eqprime->attach( enable_trapeze->isChecked() ? plot : NULL );
    start_button->setEnabled( enable_eiler->isChecked() or enable_trapeze->isChecked() );
}

void StartingPage::start()
{
    enable_everything(false);

    const Params p = collect_params();
    const size_t n = (p.Tstop - p.Tstart)/p.dt;

    plot->setAxisScale(QwtPlot::xBottom,p.Tstart-p.dt,p.Tstop+p.dt);

    if( enable_eiler->isChecked() )
    {
        progress_bar_eiler->setMaximum(n);
        progress_bar_eiler->setValue(0);
        progress_bar_eiler->setFormat("Eiler: %v");

        CalculusEiler e;
        connect(&e, SIGNAL(a_step_done()), SLOT(increment_eiler()));
        plot_answer_eiler(e.doWork(p));
    }
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
    if( enable_trapeze->isChecked() )
    {
        progress_bar_trapeze->setMaximum(n);
        progress_bar_trapeze->setValue(0);
        progress_bar_trapeze->setFormat("Trapeze: %v");

        CalculusTrapeze t;
        connect(&t, SIGNAL(a_step_done()), SLOT(increment_trapeze()));
        plot_answer_trapeze(t.doWork(p));
    }

    plot->setAxisAutoScale(QwtPlot::xBottom);
    plot->setAxisAutoScale(QwtPlot::yLeft);
    plot->setVisible(true);
    plot->replot();
    zoom->setZoomBase();

    enable_everything(true);
}

void StartingPage::enable_everything(bool e)
{
    K0f->setEnabled(e); K1f->setEnabled(e); Ur0->setEnabled(e); K0u->setEnabled(e); K1u->setEnabled(e);
    Y11->setEnabled(e); Y12->setEnabled(e); A11->setEnabled(e); A12->setEnabled(e);
    Y11em->setEnabled(e); Y12em->setEnabled(e); A11em->setEnabled(e); A12em->setEnabled(e); Pd->setEnabled(e);
    Delta0->setEnabled(e); Eqe0->setEnabled(e); Eqprime0->setEnabled(e); U0->setEnabled(e); V0->setEnabled(e);
    Tstart->setEnabled(e); Tstop->setEnabled(e); dt->setEnabled(e);
    eps->setEnabled(e); max_iterations->setEnabled(e);
    start_button->setEnabled(e);
}

Params StartingPage::collect_params()
{
    Params p = {};
    p.reg.K0f = K0f->text().toDouble();
    p.reg.K1f = K1f->text().toDouble();
    p.reg.Ur0 = Ur0->text().toDouble();
    p.reg.K0u = K0u->text().toDouble();
    p.reg.K1u = K1u->text().toDouble();
    p.repl.Pd = Pd->text().toDouble();
    p.repl.Y11 = Y11->text().toDouble();
    p.repl.Y12 = Y12->text().toDouble();
    p.repl.A11 = A11->text().toDouble() * 3.1415926 / 180;
    p.repl.A12 = A12->text().toDouble() * 3.1415926 / 180;
    p.repl.Y11em = Y11em->text().toDouble();
    p.repl.Y12em = Y12em->text().toDouble();
    p.repl.A11em = A11em->text().toDouble() * 3.1415926 / 180;
    p.repl.A12em = A12em->text().toDouble() * 3.1415926 / 180;
    p.start.Delta0 = Delta0->text().toDouble();
    p.start.Eqe0 = Eqe0->text().toDouble();
    p.start.Eqprime0 = Eqprime0->text().toDouble();
    p.start.U0 = U0->text().toDouble();
    p.start.V0 = V0->text().toDouble();
    p.Tstart = Tstart->text().toDouble();
    p.Tstop = Tstop->text().toDouble();
    p.dt = dt->text().toDouble();
    p.eps = eps->text().toDouble();
    p.max_iterations = max_iterations->text().toUInt();
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
        const AnswerItem& a = ans.at(i);
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
        const AnswerItem& a = ans.at(i);
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

