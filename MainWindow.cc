#include "MainWindow.h"
#include "Calculus.h"
#include "params.h"

#include <QApplication>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QProgressBar>
#include <QPushButton>
#include <QLineEdit>
#include <QCheckBox>
#include <QGroupBox>
#include <QLabel>
#include <QDoubleValidator>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_grid.h>


int main(int argc, char *argv[])
{
  QApplication qapp(argc,argv);

  MainWindow w;
  w.show();
  return qapp.exec();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

MainWindow::MainWindow() : QMainWindow(NULL)
{
    setCentralWidget(new CentralWidget);
    setFixedSize(1600,900);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

CentralWidget::CentralWidget() : QWidget(NULL)
{
    QVBoxLayout* l = new QVBoxLayout(this);
    tabs = new QTabWidget(this);
    tabs->setTabsClosable(false);
    tabs->setMovable(false);
    l->addWidget(tabs,1);

    start_button = new QPushButton("Start",this);
    start_button->setEnabled(false);
    l->addWidget(start_button);

    tabs->addTab(calculus = new CalculusWidget(this), "Calculus");
    tabs->addTab(sysparams = new SystemParamsWidget(this), "System Parameters");

    connect(start_button, SIGNAL(clicked()), SLOT(start()));
    connect(calculus, SIGNAL(enable_start_button(bool)), start_button, SLOT(setEnabled(bool)));
}

void CentralWidget::start()
{
    tabs->setCurrentIndex(0);
    start_button->setEnabled(false);
    calculus->start( sysparams->collect_params() );
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

static QLineEdit* add_double_input(QBoxLayout* layout, const char* label, const char* value, QValidator* validator)
{
    layout->addStretch(1);
    layout->addWidget(new QLabel(label));
    QLineEdit* edit = new QLineEdit(value);
    edit->setValidator(validator);
    layout->addWidget(edit);
    return edit;
}

static QwtPlotCurve* add_plot_curve(const char* label, int r, int g, int b)
{
    QwtPlotCurve* curve = new QwtPlotCurve(label);
    curve->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve->setPen(QPen(QColor(r,g,b)));
    curve->setVisible(false);
    return curve;
}

CalculusWidget::CalculusWidget(QWidget* p) : QWidget(p)
{
    QDoubleValidator* validator = new QDoubleValidator(this);

    QHBoxLayout* eqv = new QHBoxLayout;
    eqv->addWidget(new QLabel("Equivalent circuit parameters:",this));
    Y11 = add_double_input(eqv, "Y11:", "0.496", validator);
    Y12 = add_double_input(eqv, "Y12:", "0.532", validator);
    A11 = add_double_input(eqv, "A11:", "0.958", validator);
    A12 = add_double_input(eqv, "A12:", "0.715", validator);
    Pd = add_double_input(eqv, "Pd:", "5", validator);
    Y11em = add_double_input(eqv, "Y11em:", "0.35", validator);
    Y12em = add_double_input(eqv, "Y12em:", "0.45", validator);
    A11em = add_double_input(eqv, "A11em:", "0.9", validator);
    A12em = add_double_input(eqv, "A12em:", "0.3", validator);

    QHBoxLayout* srt = new QHBoxLayout;
    srt->addWidget(new QLabel("Initial values:",this));
    Delta0 = add_double_input(srt, "Delta0:", "1.365", validator);
    Eqe0 = add_double_input(srt, "Eqe0:", "1.87", validator);
    Eqprime0 = add_double_input(srt, "Eqprime0:", "1.1", validator);
    U0 = add_double_input(srt, "U0:", "1.059", validator);
    V0 = add_double_input(srt, "V0:", "0.597", validator);

    QHBoxLayout* cpt = new QHBoxLayout;
    cpt->addWidget(new QLabel("Computational parameters:",this));
    Tstart = add_double_input(cpt, "Tstart:", "0", validator);
    Tstop = add_double_input(cpt, "Tstop:", "3", validator);
    dt = add_double_input(cpt, "dt:", "0.0001", validator);
    eps = add_double_input(cpt, "eps:", "0.005", validator);
    max_iterations = add_double_input(cpt, "max_iterations:", "50", validator);
    max_iterations->setValidator(new QIntValidator(5,5000,this));

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

    plot = new QwtPlot(this);

    curve_eiler_Delta = add_plot_curve("Delta Eiler",160,0,210);
    curve_eiler_Omega = add_plot_curve("Omega Eiler",210,0,0);
    curve_eiler_Eqe = add_plot_curve("Eqe Eiler",5,205,0);
    curve_eiler_Eqprime = add_plot_curve("Eqprime Eiler",210,120,0);
    curve_eiler_U = add_plot_curve("U Eiler",0,0,0);

    curve_trapeze_Delta = add_plot_curve("Delta Trapeze",210,50,255);
    curve_trapeze_Omega = add_plot_curve("Omega Trapeze",255,45,45);
    curve_trapeze_Eqe = add_plot_curve("Eqe Trapeze",50,255,50);
    curve_trapeze_Eqprime = add_plot_curve("Eqprime Trapeze",255,160,50);
    curve_trapeze_U = add_plot_curve("U Trapeze",60,60,60);

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
    l->addLayout(eqv);
    l->addLayout(srt);
    l->addLayout(cpt);
    l->addLayout(prb);
    l->addWidget(plot,2);

    plot->setVisible(true);
    progress_bar_eiler->setVisible(true);
    progress_bar_trapeze->setVisible(true);

    connect(enable_eiler, SIGNAL(clicked()), SLOT(wanna_eiler()));
    connect(enable_trapeze, SIGNAL(clicked()), SLOT(wanna_trapeze()));
}

CalculusWidget::~CalculusWidget()
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
    if(not curve_trapeze_U->plot() )  delete curve_trapeze_U;
}

void CalculusWidget::wanna_eiler()
{
    progress_bar_eiler->setEnabled( enable_eiler->isChecked() );
    curve_eiler_U->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Delta->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Omega->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Eqe->attach( enable_eiler->isChecked() ? plot : NULL );
    curve_eiler_Eqprime->attach( enable_eiler->isChecked() ? plot : NULL );
    emit enable_start_button( enable_eiler->isChecked() or enable_trapeze->isChecked() );
}

void CalculusWidget::wanna_trapeze()
{
    progress_bar_trapeze->setEnabled( enable_trapeze->isChecked() );
    curve_trapeze_Delta->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_Omega->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_Eqe->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_Eqprime->attach( enable_trapeze->isChecked() ? plot : NULL );
    curve_trapeze_U->attach( enable_trapeze->isChecked() ? plot : NULL );
    emit enable_start_button( enable_eiler->isChecked() or enable_trapeze->isChecked() );
}

void CalculusWidget::start(const Params::Consts& reg)
{
    enable_everything(false);

    Params p = collect_params();
    p.reg = reg;
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
    emit enable_start_button(true);
}

void CalculusWidget::enable_everything(bool e)
{
    Y11->setEnabled(e); Y12->setEnabled(e); A11->setEnabled(e); A12->setEnabled(e);
    Y11em->setEnabled(e); Y12em->setEnabled(e); A11em->setEnabled(e); A12em->setEnabled(e); Pd->setEnabled(e);
    Delta0->setEnabled(e); Eqe0->setEnabled(e); Eqprime0->setEnabled(e); U0->setEnabled(e); V0->setEnabled(e);
    Tstart->setEnabled(e); Tstop->setEnabled(e); dt->setEnabled(e);
    eps->setEnabled(e); max_iterations->setEnabled(e);
}

Params CalculusWidget::collect_params()
{
    Params p = {};
    p.repl.Pd = Pd->text().toDouble();
    p.repl.Y11 = Y11->text().toDouble();
    p.repl.Y12 = Y12->text().toDouble();
    p.repl.A11 = A11->text().toDouble();
    p.repl.A12 = A12->text().toDouble();
    p.repl.Y11em = Y11em->text().toDouble();
    p.repl.Y12em = Y12em->text().toDouble();
    p.repl.A11em = A11em->text().toDouble();
    p.repl.A12em = A12em->text().toDouble();
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

void CalculusWidget::increment_eiler()
{
    progress_bar_eiler->setValue(progress_bar_eiler->value()+1);
}

void CalculusWidget::increment_trapeze()
{
    progress_bar_trapeze->setValue(progress_bar_trapeze->value()+1);
}

void CalculusWidget::plot_answer_eiler(const QVector<AnswerItem> ans)
{
    qDebug()<<__FUNCTION__<<ans.size();

    QVector<QPointF> data_delta(ans.size());
    QVector<QPointF> data_omega(ans.size());
    QVector<QPointF> data_eqe(ans.size());
    QVector<QPointF> data_eqp(ans.size());
    QVector<QPointF> data_u(ans.size());
    for(int i=0; i<ans.size(); i++)
    {
        const AnswerItem& a = ans.at(i);
        data_delta[i] = QPointF(a.time, a.delta);
        data_omega[i] = QPointF(a.time, a.omega);
        data_eqp[i] = QPointF(a.time, a.Eqprime);
        data_eqe[i] = QPointF(a.time, a.Eqe);
        data_u[i] = QPointF(a.time, a.U);
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
/*
    double min = std::min(std::min(std::min(min_delta,min_omega),min_Eqprime),min_Eqe);
    double max = std::max(std::max(std::max(max_delta,max_omega),max_Eqprime),max_Eqe);
    plot->setAxisScale(QwtPlot::yLeft,min,max);
    qDebug()<<"plot->setVisible(true): max"<<max<<"max_delta"<<max_delta<<"max_omega"<<max_omega<<"max_Eqprime"<<max_Eqprime<<"max_Eqe"<<max_Eqe;
    qDebug()<<"plot->setVisible(true): min"<<min<<"min_delta"<<min_delta<<"min_omega"<<min_omega<<"min_Eqprime"<<min_Eqprime<<"min_Eqe"<<min_Eqe;
    */
}

void CalculusWidget::plot_answer_trapeze(const QVector<AnswerItem> ans)
{
    qDebug()<<__FUNCTION__<<ans.size();

    QVector<QPointF> data_delta(ans.size());
    QVector<QPointF> data_omega(ans.size());
    QVector<QPointF> data_eqe(ans.size());
    QVector<QPointF> data_eqp(ans.size());
    QVector<QPointF> data_u(ans.size());
    for(int i=0; i<ans.size(); i++)
    {
        const AnswerItem& a = ans.at(i);
        data_delta[i] = QPointF(a.time, a.delta);
        data_omega[i] = QPointF(a.time, a.omega);
        data_eqp[i] = QPointF(a.time, a.Eqprime);
        data_eqe[i] = QPointF(a.time, a.Eqe);
        data_u[i] = QPointF(a.time, a.U);
    }
    curve_trapeze_Delta->setSamples(data_delta);
    curve_trapeze_Omega->setSamples(data_omega);
    curve_trapeze_Eqprime->setSamples(data_eqp);
    curve_trapeze_Eqe->setSamples(data_eqe);
    curve_trapeze_U->setSamples(data_u);
    curve_trapeze_Delta->setVisible(true);
    curve_trapeze_Omega->setVisible(true);
    curve_trapeze_Eqprime->setVisible(true);
    curve_trapeze_Eqe->setVisible(true);
    curve_trapeze_U->setVisible(true);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

SystemParamsWidget::SystemParamsWidget(QWidget* p) : QWidget(p)
{
    QDoubleValidator* validator = new QDoubleValidator(this);

    QGroupBox *reg = new QGroupBox("Regulator parameters:", this);
    QHBoxLayout* reg_l = new QHBoxLayout(reg);
    K0f = add_double_input(reg_l, "K0f:", "10", validator);
    K1f = add_double_input(reg_l, "K1f:", "0.5", validator);
    Ur0 = add_double_input(reg_l, "Ur0:", "1.05", validator);
    K0u = add_double_input(reg_l, "K0u:", "50", validator);
    K1u = add_double_input(reg_l, "K1u:", "3.6", validator);

    QGroupBox* tms = new QGroupBox("Time scales:",this);
    QHBoxLayout* tms_l = new QHBoxLayout(tms);
    T = add_double_input(tms_l, "T:", "1", validator);
    Ty = add_double_input(tms_l, "Tj:", "5.1", validator);
    Tu = add_double_input(tms_l, "Tu:", "0.05", validator);
    Tg = add_double_input(tms_l, "Tg:", "0.02", validator);
    Te = add_double_input(tms_l, "Te:", "0.03", validator);
    Tf = add_double_input(tms_l, "Tf:", "0.05", validator);
    Tphi = add_double_input(tms_l, "Tphi:", "0.05", validator);
    Td0 = add_double_input(tms_l, "Td0:", "5.87", validator);

    QGroupBox* oth = new QGroupBox("Other constants:",this);
    QHBoxLayout* oth_l = new QHBoxLayout(oth);
    Pt0 = add_double_input(oth_l, "Pt0:", "1", validator);
    omega_nom = add_double_input(oth_l, "omega_nom:", "18000", validator);
    Xdprime = add_double_input(oth_l, "Xd':", "0.207", validator);
    Xd = add_double_input(oth_l, "Xd:", "1.364", validator);
    Eqenom = add_double_input(oth_l, "Eqenom:", "1.87", validator);
    Uc = add_double_input(oth_l, "Uc:", "1", validator);

    QVBoxLayout* l = new QVBoxLayout(this);
    l->addWidget(reg);
    l->addWidget(tms);
    l->addWidget(oth);
}

Params::Consts SystemParamsWidget::collect_params()
{
    Params::Consts c;
    c.K0f = K0f->text().toDouble();
    c.K1f = K1f->text().toDouble();
    c.Ur0 = Ur0->text().toDouble();
    c.K0u = K0u->text().toDouble();
    c.K1u = K1u->text().toDouble();
    c.T = T->text().toDouble();
    c.Ty = Ty->text().toDouble();
    c.Tu = Tu->text().toDouble();
    c.Tg = Tg->text().toDouble();
    c.Te = Te->text().toDouble();
    c.Tf = Tf->text().toDouble();
    c.Tphi = Tphi->text().toDouble();
    c.Td0 = Td0->text().toDouble();
    c.Pt0 = Pt0->text().toDouble();
    c.omega_nom = omega_nom->text().toDouble();
    c.Xdprime = Xdprime->text().toDouble();
    c.Xd = Xd->text().toDouble();
    c.Eqenom = Eqenom->text().toDouble();
    c.Uc = Uc->text().toDouble();
    return c;
}

