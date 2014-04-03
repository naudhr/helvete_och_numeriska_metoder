#include "MainWindow.h"
#include "Calculus.h"
#include "params.h"
#include "NoQwt.h"

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
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QShortcut>
#include <QHeaderView>
#include <QFileDialog>
#include <QFile>

#include <QDebug>
#include <QMessageBox>
#include <cmath>

int main(int argc, char *argv[])
{
  QApplication qapp(argc,argv);

  qRegisterMetaType<AnswerItem>("AnswerItem");
  qRegisterMetaType<AnswerItem>("AnswerItem");

  MainWindow w;
  QObject::connect(&w,SIGNAL(quit()),&qapp,SLOT(quit()));

  w.showMaximized();
  return qapp.exec();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

MainWindow::MainWindow() : QMainWindow(NULL)
{
    setCentralWidget(new CentralWidget);
    /*QShortcut* q =*/ new QShortcut(Qt::CTRL + Qt::Key_Q, this, SIGNAL(quit()));
    //setFixedSize(1600,900);
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
    tabs->addTab(to_excel = new ToExcel(this), "Table Data");

    connect(start_button, SIGNAL(clicked()), SLOT(start()));
    connect(calculus, SIGNAL(enable_start_button(bool)), start_button, SLOT(setEnabled(bool)));
    connect(calculus, SIGNAL(to_excel_populate(AnswerItem)), to_excel, SLOT(populate(AnswerItem)), Qt::QueuedConnection);
    connect(calculus, SIGNAL(enable_start_button(bool)), to_excel, SLOT(make_up(bool)), Qt::QueuedConnection);
    connect(sysparams, SIGNAL(params_changed(Params::Consts)), calculus, SLOT(sysparams_changed(Params::Consts)));

    calculus->sysparams_changed( sysparams->collect_params() );
}

void CentralWidget::start()
{
    tabs->setCurrentIndex(0);
    start_button->setEnabled(false);
    calculus->start( sysparams->collect_params() );
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

static QLineEdit* add_double_input(QBoxLayout* layout, const QString& label, const char* value, QValidator* validator)
{
    layout->addStretch(1);
    layout->addWidget(new QLabel(label));
    QLineEdit* edit = new QLineEdit(value);
    edit->setValidator(validator);
    layout->addWidget(edit);
    return edit;
}

static NoQwtPlotCurve* add_plot_curve(NoQwtPlot* plot, const QString& label, int r, int g, int b, const QString& tag)
{
    NoQwtPlotCurve* curve = new NoQwtPlotCurve(plot, label, QPen(QColor(r,g,b)), tag);
    curve->setVisible(true);
    return curve;
}

CalculusWidget::CalculusWidget(QWidget* p) : QWidget(p)
{
    QDoubleValidator* validator = new QDoubleValidator(this);

    QHBoxLayout* eqv = new QHBoxLayout;
    eqv->addWidget(new QLabel("Equivalent circuit parameters:",this));
    Y11 = add_double_input(eqv, "Y11:", "2.523", validator);
    Y12 = add_double_input(eqv, "Y12:", "2.726", validator);
    A11 = add_double_input(eqv, QChar(0x03b1)+QLatin1String("11:"), "0.089", validator);
    A12 = add_double_input(eqv, QChar(0x03b1)+QLatin1String("12:"), "0.08", validator);
    Pd = add_double_input(eqv, "Pd:", "5", validator);
    Y11em = add_double_input(eqv, "Y11av:", "3.085", validator);
    Y12em = add_double_input(eqv, "Y12av:", "1.862", validator);
    A11em = add_double_input(eqv, QChar(0x03b1)+QLatin1String("11av:"), "0.081", validator);
    A12em = add_double_input(eqv, QChar(0x03b1)+QLatin1String("12av:"), "0.078", validator);

    QHBoxLayout* cpt = new QHBoxLayout;
    cpt->addWidget(new QLabel("Computational parameters:",this));
    Tstart = add_double_input(cpt, "Tstart:", "0", validator);
    Tstop = add_double_input(cpt, "Tstop:", "3", validator);
    dt = add_double_input(cpt, "dt:", "0.001", validator);
    eps = add_double_input(cpt, "eps:", "0.005", validator);
    max_iterations = add_double_input(cpt, "max_iterations:", "50", validator);
    max_iterations->setValidator(new QIntValidator(5,5000,this));
    cpt->addSpacing(30);
    cpt->addStretch(2);
    cpt->addWidget(new QLabel("Initial values:",this));
    Delta0 = add_double_input(cpt, QChar(0x03b4)+QLatin1String("(0):"), "1.294", validator);
    Eqe0 = add_double_input(cpt, "Eqe(0):", "1.698", validator);
    Eqprime0 = add_double_input(cpt, "E'q(0):", "0.759", validator);
    U0 = add_double_input(cpt, "U(0):", "1.035", validator);
    V0 = add_double_input(cpt, QChar(0x03bd)+QLatin1String("(0):"), "0.353", validator);

    connect(U0, SIGNAL(editingFinished()), SLOT(sysparams_changed()));
    connect(V0, SIGNAL(editingFinished()), SLOT(sysparams_changed()));
    connect(Eqe0, SIGNAL(editingFinished()), SLOT(sysparams_changed()));

    seq_params = new QWidget(this);
    QHBoxLayout* srt = new QHBoxLayout(seq_params);
    srt->addWidget(new QLabel("Sequensive related"));
    seq_x1_0 = add_double_input(srt, "x1(0)", "0", validator);
    seq_x2_0 = add_double_input(srt, "x2(0)", "0", validator);
    seq_x4_0 = add_double_input(srt, "x4(0)", "0", validator);
    seq_x5_0 = add_double_input(srt, "x5(0)", "0", validator);
    seq_x6_0 = add_double_input(srt, "x6(0)", "0", validator);
    seq_x7_0 = add_double_input(srt, "x7(0)", "0", validator);
    seq_x9_0 = add_double_input(srt, "x9(0)", "0", validator);
    seq_x0_0 = add_double_input(srt, "x10(0)", "0", validator);

    power_part = new QPushButton(QString::fromWCharArray(L"W11, W12, Pг, Pс, Qг, Qс"), this);
    connect(power_part, SIGNAL(clicked()), SLOT(popup_power_widget()));

    QHBoxLayout* prb = new QHBoxLayout;
    online_plotting = new QCheckBox("Online",this);
    prb->addWidget(online_plotting);
    dirty_hack = new QCheckBox("Dirty hack",this);
    prb->addWidget(dirty_hack);
    prb->addSpacing(30);

    enable_eiler = new QCheckBox(this);
    enable_eiler->setChecked(false);
    progress_bar_eiler = new QProgressBar(this);
    progress_bar_eiler->setFormat("Eiler");
    progress_bar_eiler->setMinimum(0);
    progress_bar_eiler->setMaximum(1);
    progress_bar_eiler->setValue(0);
    progress_bar_eiler->setEnabled(false);
    prb->addWidget(enable_eiler);
    prb->addWidget(progress_bar_eiler);
    prb->addWidget(power_part);
    prb->addSpacing(30);
    enable_trapeze = new QCheckBox(this);
    enable_trapeze->setChecked(false);
    progress_bar_trapeze = new QProgressBar(this);
    progress_bar_trapeze->setFormat("Trapeze");
    progress_bar_trapeze->setMinimum(0);
    progress_bar_trapeze->setMaximum(1);
    progress_bar_trapeze->setValue(0);
    progress_bar_trapeze->setEnabled(false);
    prb->addWidget(enable_trapeze);
    prb->addWidget(progress_bar_trapeze);
    prb->addSpacing(30);
    enable_sequensive = new QCheckBox(this);
    enable_sequensive->setChecked(false);
    progress_bar_sequensive = new QProgressBar(this);
    progress_bar_sequensive->setFormat("Sequensive");
    progress_bar_sequensive->setMinimum(0);
    progress_bar_sequensive->setMaximum(1);
    progress_bar_sequensive->setValue(0);
    progress_bar_sequensive->setEnabled(false);
    prb->addWidget(enable_sequensive);
    prb->addWidget(progress_bar_sequensive);
    prb->addSpacing(30);
    enable_parallel = new QCheckBox(this);
    enable_parallel->setChecked(false);
    enable_parallel->setEnabled(false);
    progress_bar_parallel = new QProgressBar(this);
    progress_bar_parallel->setFormat("Parallel");
    progress_bar_parallel->setMinimum(0);
    progress_bar_parallel->setMaximum(1);
    progress_bar_parallel->setValue(0);
    progress_bar_parallel->setEnabled(false);
    prb->addWidget(enable_parallel);
    prb->addWidget(progress_bar_parallel);

    view = new NoQwtGraphicsView(this);
    plot = new NoQwtPlot(view, "", "", "");

    add_plot_curve(plot, QChar(0x03b4)+QLatin1String(" Eiler"),160,0,210, "delta E");
    add_plot_curve(plot, QChar(0x0394)+QString(0x03c9)+QLatin1String(" Eiler"),210,0,0, "omega E");
    add_plot_curve(plot, "Eqe Eiler",5,205,0, "eqe E");
    add_plot_curve(plot, "E'q Eiler",210,120,0, "eqp E");
    add_plot_curve(plot, "U Eiler",0,0,0, "u E");
    add_plot_curve(plot, QChar(0x03bd)+QLatin1String(" Eiler"),0,100,100, "v E");

    add_plot_curve(plot, QChar(0x03b4)+QLatin1String(" Trapeze"),210,50,255, "delta T");
    add_plot_curve(plot, QChar(0x0394)+QString(0x03c9)+QLatin1String(" Trapeze"),255,45,45, "omega T");
    add_plot_curve(plot, "Eqe Trapeze",50,255,50, "eqe T");
    add_plot_curve(plot, "E'q Trapeze",255,160,50, "eqp T");
    add_plot_curve(plot, "U Trapeze",60,60,60, "u T");
    add_plot_curve(plot, QChar(0x03bd)+QLatin1String(" Trapeze"),0,180,180, "v T");

    add_plot_curve(plot, QChar(0x03b4)+QLatin1String(" Sequensive"),210,50,255, "delta S");
    add_plot_curve(plot, QChar(0x0394)+QString(0x03c9)+QLatin1String(" Sequensive"),255,45,45, "omega S");
    add_plot_curve(plot, "Eqe Sequensive",50,255,50, "eqe S");
    add_plot_curve(plot, "E'q Sequensive",255,160,50, "eqp S");
    add_plot_curve(plot, "U Sequensive",60,60,60, "u S");
    add_plot_curve(plot, QChar(0x03bd)+QLatin1String(" Sequensive"),0,180,180, "v S");

    /*NoQwtPlotCurve* w1 =*/ add_plot_curve(plot, "W11 Power", 255, 0, 0, "W11");
    /*NoQwtPlotCurve* w2 =*/ add_plot_curve(plot, "W12 Power", 0, 255, 0, "W12");
    /*NoQwtPlotCurve* pg =*/ add_plot_curve(plot, QString::fromWCharArray(L"Pг Power"), 0, 0, 255, "pg");
    /*NoQwtPlotCurve* pc =*/ add_plot_curve(plot, QString::fromWCharArray(L"Pс Power"), 255, 255, 0, "pc");
    /*NoQwtPlotCurve* qg =*/ add_plot_curve(plot, QString::fromWCharArray(L"Qг Power"), 255, 0, 255, "qg");
    /*NoQwtPlotCurve* qc =*/ add_plot_curve(plot, QString::fromWCharArray(L"Qс Power"), 0, 255, 255, "qc");

    QVBoxLayout* l = new QVBoxLayout(this);
    l->addLayout(eqv);
    l->addLayout(cpt);
    l->addWidget(seq_params);
    l->addLayout(prb);
    l->addWidget(view,2);

    plot->setVisible(false);
    //progress_bar_eiler->setVisible(true);
    //progress_bar_trapeze->setVisible(true);

    connect(enable_eiler, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    connect(enable_trapeze, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    connect(enable_sequensive, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    connect(enable_parallel, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    jobs = 0;

    some_calc_enabled();
}

CalculusWidget::~CalculusWidget()
{
}

void CalculusWidget::some_calc_enabled()
{
    progress_bar_eiler->setEnabled( enable_eiler->isChecked() );
    progress_bar_trapeze->setEnabled( enable_trapeze->isChecked() );
    progress_bar_sequensive->setEnabled( enable_sequensive->isChecked() );
    progress_bar_parallel->setEnabled( enable_parallel->isChecked() );
    emit enable_start_button( enable_eiler->isChecked() or
                              enable_trapeze->isChecked() or
                              enable_sequensive->isChecked() or
                              enable_parallel->isChecked() );

    seq_params->setVisible(enable_sequensive->isChecked());
    power_part->setVisible(enable_eiler->isChecked() and not enable_trapeze->isChecked() and not enable_sequensive->isChecked() and not enable_parallel->isChecked());

    view->legend()->setVisibleSection("Eiler", enable_eiler->isChecked());
    view->legend()->setVisibleSection("Power", enable_eiler->isChecked());
    view->legend()->setVisibleSection("Trapeze", enable_trapeze->isChecked());
    view->legend()->setVisibleSection("Sequensive", enable_sequensive->isChecked());
    view->legend()->setVisibleSection("Parallel", enable_parallel->isChecked());
    view->legend()->setVisible(true);
}

void CalculusWidget::popup_power_widget()
{
    const Params::Consts& r = collected_params.reg;
    Equiv e;

    NoQwtPlotCurve* w1 = plot->curve("W11");
    NoQwtPlotCurve* w2 = plot->curve("W12");
    NoQwtPlotCurve* pg = plot->curve("pg");
    NoQwtPlotCurve* pc = plot->curve("pc");
    NoQwtPlotCurve* qg = plot->curve("qg");
    NoQwtPlotCurve* qc = plot->curve("qc");

    w1->reset();
    w2->reset();
    pg->reset();
    pc->reset();
    qg->reset();
    qc->reset();

    const NoQwtPlotCurve* d_E = plot->curve("delta E");
    const NoQwtPlotCurve* p_E = plot->curve("eqp E");
    const NoQwtPlotCurve* u_E = plot->curve("u E");
    const NoQwtPlotCurve* v_E = plot->curve("v E");

    for(int i=0, z=1; i<z; i++)
    {
        const QVector<QPointF> d = d_E->points();
        z = d.size();
        if(i >= z)
            break;

        const double t = d[i].x();
        const double delta = d[i].y();
        const double Eqprime = p_E->points()[i].y();
        const double U = u_E->points()[i].y();
        const double V = v_E->points()[i].y();
        const double s_d_v = std::sin(delta - V);
        const double c_d_v = std::cos(delta - V);
        e = recalculate_equiv_params(t, U, collected_params, e);

        const double Pg = CalculusEiler::calculate_Pg(Eqprime, U, r.Xdprime, r.Xd, s_d_v, c_d_v);
        const double Pc = CalculusEiler::calculate_Pc(U, V, e.Y11, e.Y12, e.A11, e.A12, r.Uc);
        const double Qg = CalculusEiler::calculate_Qg(Eqprime, U, r.Xdprime, r.Xd, s_d_v, c_d_v);
        const double Qc = CalculusEiler::calculate_Qc(U, V, e.Y11, e.Y12, e.A11, e.A12, r.Uc);
        const double W1 = Pg - Pc;
        const double W2 = Qg - Qc;

        w1->addData(t, W1);
        w2->addData(t, W2);
        pg->addData(t, Pg);
        pc->addData(t, Pc);
        qg->addData(t, Qg);
        qc->addData(t, Qc);
    }
}

void CalculusWidget::sysparams_changed()
{
    this->sysparams_changed(collected_params.reg);
}

void CalculusWidget::sysparams_changed(Params::Consts sc)
{
    collect_params();
    collected_params.reg = sc;

    const Params& p = collected_params;

    seq_x1_0->setEnabled(false);
    seq_x1_0->setText(U0->text());
    seq_x2_0->setEnabled(false);
    seq_x2_0->setText(QString::number(p.reg.Ur0 - p.start.U0));
    //seq_x4_0->setEditable(false);
    //seq_x4_0->setValue(QString::number(p.reg.Ur0 - p.start.U0));
    seq_x5_0->setEnabled(false);
    seq_x5_0->setText(V0->text());
    seq_x9_0->setEnabled(false);
    seq_x9_0->setText(Eqe0->text());
    seq_x0_0->setEnabled(false);
    seq_x0_0->setText(Eqe0->text());
}

void CalculusWidget::start(const Params::Consts& reg)
{
    enable_everything(false);

    collect_params();
    collected_params.reg = reg;
    sysparams_changed(reg);

    const Params& p = collected_params;

    const size_t n = (p.Tstop - p.Tstart)/p.dt;
    answer_buffer.clear();
    answer_buffer.reserve(n * 4);

    if( enable_eiler->isChecked() )
    {
        progress_bar_eiler->setMaximum(n);
        progress_bar_eiler->setValue(0);
        progress_bar_eiler->setFormat("Eiler: %v");

        QThread* c = new CalculusEiler(p);
        connect(c, SIGNAL(a_step_done(AnswerItem)), this, SLOT(eiler_step(AnswerItem)), Qt::QueuedConnection);
        connect(c, SIGNAL(finished()), this, SLOT(a_part_of_the_plot_done()), Qt::QueuedConnection);
        connect(c, SIGNAL(newton_does_not_converge(QString,double,unsigned)), SLOT(ndnc(QString,double,unsigned)));
        connect(c, SIGNAL(finished()), c, SLOT(deleteLater()));
        c->start();
        jobs ++;
    }

    if( enable_trapeze->isChecked() )
    {
        progress_bar_trapeze->setMaximum(n);
        progress_bar_trapeze->setValue(0);
        progress_bar_trapeze->setFormat("Trapeze: %v");

        QThread* c = new CalculusTrapeze(p);
        connect(c, SIGNAL(a_step_done(AnswerItem)), this, SLOT(trapeze_step(AnswerItem)), Qt::QueuedConnection);
        connect(c, SIGNAL(finished()), this, SLOT(a_part_of_the_plot_done()), Qt::QueuedConnection);
        connect(c, SIGNAL(newton_does_not_converge(QString,double,unsigned)), SLOT(ndnc(QString,double,unsigned)));
        connect(c, SIGNAL(finished()), c, SLOT(deleteLater()));
        c->start();
        jobs ++;
    }

    if( enable_sequensive->isChecked() )
    {
        progress_bar_sequensive->setMaximum(n);
        progress_bar_sequensive->setValue(0);
        progress_bar_sequensive->setFormat("Sequensive: %v");

        CalculusSequensive* c = new CalculusSequensive(p);

        c->set_X(seq_x1_0->text().toDouble(), seq_x2_0->text().toDouble(),
                 seq_x4_0->text().toDouble(), seq_x5_0->text().toDouble(),
                 seq_x6_0->text().toDouble(), seq_x7_0->text().toDouble(),
                 seq_x9_0->text().toDouble(), seq_x0_0->text().toDouble());

        connect(c, SIGNAL(a_step_done(AnswerItem)), this, SLOT(sequensive_step(AnswerItem)), Qt::QueuedConnection);
        connect(c, SIGNAL(finished()), this, SLOT(a_part_of_the_plot_done()), Qt::QueuedConnection);
        connect(c, SIGNAL(newton_does_not_converge(QString,double,unsigned)), SLOT(ndnc(QString,double,unsigned)));
        connect(c, SIGNAL(finished()), c, SLOT(deleteLater()));
        c->start();
        jobs ++;
    }

    if( enable_parallel->isChecked() )
    {
        progress_bar_parallel->setMaximum(n);
        progress_bar_parallel->setValue(0);
        progress_bar_parallel->setFormat("Parallel: %v");

        QThread* c = new CalculusParallel(p);
        connect(c, SIGNAL(a_step_done(AnswerItem)), this, SLOT(parallel_step(AnswerItem)), Qt::QueuedConnection);
        connect(c, SIGNAL(finished()), this, SLOT(a_part_of_the_plot_done()), Qt::QueuedConnection);
        connect(c, SIGNAL(newton_does_not_converge(QString,double,unsigned)), SLOT(ndnc(QString,double,unsigned)));
        connect(c, SIGNAL(finished()), c, SLOT(deleteLater()));
        c->start();
        jobs ++;
    }

    plot->reset();
    plot->setVisible(online_plotting->isChecked());

    enable_everything(jobs == 0);
    emit enable_start_button(jobs == 0);
}

void CalculusWidget::enable_everything(bool e)
{
    online_plotting->setEnabled(e);
    Y11->setEnabled(e); Y12->setEnabled(e); A11->setEnabled(e); A12->setEnabled(e);
    Y11em->setEnabled(e); Y12em->setEnabled(e); A11em->setEnabled(e); A12em->setEnabled(e); Pd->setEnabled(e);
    Delta0->setEnabled(e); Eqe0->setEnabled(e); Eqprime0->setEnabled(e); U0->setEnabled(e); V0->setEnabled(e);
    Tstart->setEnabled(e); Tstop->setEnabled(e); dt->setEnabled(e);
    eps->setEnabled(e); max_iterations->setEnabled(e);
}

void CalculusWidget::collect_params()
{
    Params& p = collected_params;
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
    p.dirty_hack = dirty_hack->isChecked();
}

static void plot_answer_step(NoQwtPlot* plot, const QLatin1String& suffix, const AnswerItem& a)
{
    plot->curve(QLatin1String("delta ")+suffix)->addData(a.time, a.delta);
    plot->curve(QLatin1String("omega ")+suffix)->addData(a.time, a.omega);
    plot->curve(QLatin1String("eqp ")+suffix)->addData(a.time, a.Eqprime);
    plot->curve(QLatin1String("eqe ")+suffix)->addData(a.time, a.Eqe);
    plot->curve(QLatin1String("u ")+suffix)->addData(a.time, a.U);
    plot->curve(QLatin1String("v ")+suffix)->addData(a.time, a.V);
}

void CalculusWidget::sequensive_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("S"), ans);
    progress_bar_sequensive->setValue(progress_bar_sequensive->value()+1);
    answer_buffer.append(ans);
    if(online_plotting->isChecked())
        view->fitInView(plot);
}

void CalculusWidget::parallel_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("P"), ans);
    progress_bar_parallel->setValue(progress_bar_parallel->value()+1);
    answer_buffer.append(ans);
    if(online_plotting->isChecked())
        view->fitInView(plot);
}

void CalculusWidget::trapeze_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("T"), ans);
    progress_bar_trapeze->setValue(progress_bar_trapeze->value()+1);
    answer_buffer.append(ans);
    if(online_plotting->isChecked())
        view->fitInView(plot);
}

void CalculusWidget::eiler_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("E"), ans);
    progress_bar_eiler->setValue(progress_bar_eiler->value()+1);
    answer_buffer.append(ans);
    if(online_plotting->isChecked())
        view->fitInView(plot);
}

void CalculusWidget::a_part_of_the_plot_done()
{
    jobs --;
    view->fitInView(plot);
    plot->scaled();
    plot->setVisible(true);
    view->update();

    if(jobs == 0)
    {
        foreach(const AnswerItem& a, answer_buffer)
            emit to_excel_populate(a);
        enable_everything(true);
        emit enable_start_button(true);
    }
}

void CalculusWidget::ndnc(QString name, double t, unsigned n_steps)
{
    QMessageBox::warning(this, "Newton does not converge",
                         QString("The Newton algorithm for %1 does not converge\n"
                                 "at the time %2 after %3 steps")
                             .arg(name).arg(t).arg(n_steps));
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

SystemParamsWidget::SystemParamsWidget(QWidget* p) : QWidget(p)
{
    QDoubleValidator* validator = new QDoubleValidator(this);

    QGroupBox *reg = new QGroupBox("Regulator parameters:", this);
    QHBoxLayout* reg_l = new QHBoxLayout(reg);
    K0f = add_double_input(reg_l, "K0f:", "1", validator);
    K1f = add_double_input(reg_l, "K1f:", "0.1", validator);
    Ur0 = add_double_input(reg_l, "Ur0:", "1.04", validator);
    K0u = add_double_input(reg_l, "K0u:", "50", validator);
    K1u = add_double_input(reg_l, "K1u:", "3.6", validator);

    QGroupBox* tms = new QGroupBox("Time scales:",this);
    QHBoxLayout* tms_l = new QHBoxLayout(tms);
    T = add_double_input(tms_l, "T:", "1", validator);
    Ty = add_double_input(tms_l, "Tj:", "5.59", validator);
    Tu = add_double_input(tms_l, "Tu:", "0.05", validator);
    Tg = add_double_input(tms_l, "Tg:", "0.02", validator);
    Te = add_double_input(tms_l, "Te:", "0.03", validator);
    Tf = add_double_input(tms_l, "Tf:", "0.05", validator);
    Tphi = add_double_input(tms_l, "Tphi:", "0.05", validator);
    Td0 = add_double_input(tms_l, "Td0:", "7", validator);

    QGroupBox* oth = new QGroupBox("Other constants:",this);
    QHBoxLayout* oth_l = new QHBoxLayout(oth);
    Pt0 = add_double_input(oth_l, "Pt0:", "1", validator);
    omega_nom = add_double_input(oth_l, "omega_nom:", "314.1592653589793", validator);
    Xdprime = add_double_input(oth_l, "Xd':", "0.194", validator);
    Xd = add_double_input(oth_l, "Xd:", "1.42", validator);
    Eqenom = add_double_input(oth_l, "Eqenom:", "2.79", validator);
    Uc = add_double_input(oth_l, "Uc:", "1", validator);

    connect(K0f, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(K1f, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Ur0, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(K0u, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(K1u, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(T  , SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Ty , SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Tu , SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Tg , SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Te , SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Tf , SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Tphi, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Td0, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Pt0, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(omega_nom, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Xdprime, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Xd, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Eqenom, SIGNAL(editingFinished()), SLOT(emit_params_changed()));
    connect(Uc , SIGNAL(editingFinished()), SLOT(emit_params_changed()));

    QVBoxLayout* l = new QVBoxLayout(this);
    l->addWidget(reg);
    l->addWidget(tms);
    l->addWidget(oth);
}

void SystemParamsWidget::emit_params_changed()
{
    emit params_changed(collect_params());
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

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
const size_t set_size = 6;

ToExcel::ToExcel(QWidget* p) : QWidget(p), table(NULL)
{
    table = new QTableWidget(this);
    table->setColumnCount(1 + set_size*4);
    QStringList vh;
    vh << "Time" << "n" << (QChar(0x03b4)+QLatin1String(" Eiler")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Eiler"))
                 << "Eqe Eiler" << (QChar(0x03bd)+QLatin1String(" Eiler")) << "U Eiler" /*<< "E'q Eiler"*/
                 << "n" << (QChar(0x03b4)+QLatin1String(" Trapeze")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Trapeze"))
                 << "Eqe Trapeze" << (QChar(0x03bd)+QLatin1String(" Trapeze")) << "U Trapeze" /*<< "E'q  Trapeze"*/
                 << "n" << (QChar(0x03b4)+QLatin1String(" Sequensive")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Sequensive"))
                 << "Eqe Sequensive" << (QChar(0x03bd)+QLatin1String(" Sequensive")) << "U Sequensive" /*<< "E'q  Sequensive"*/
                 << "n" << (QChar(0x03b4)+QLatin1String(" Parallel")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Parallel"))
                 << "Eqe Parallel" << (QChar(0x03bd)+QLatin1String(" Parallel")) << "U Parallel" /*<< "E'q  Parallel"*/;
    table->setHorizontalHeaderLabels(vh);
    table->verticalHeader()->setVisible(false);
    table->setSortingEnabled(false);

    QPushButton* button = new QPushButton("Save to file (then you can open it with Excel)",this);
    connect(button, SIGNAL(clicked()), SLOT(export_to_excel()));
    //button->setVisible(false);

    QVBoxLayout* l = new QVBoxLayout(this);
    l->addWidget(table);
    l->addWidget(button);
}

void ToExcel::populate(const AnswerItem& data)
{
    const int set_off = set_size * data.set_no;

    if(data.row >= table->rowCount())
        table->setRowCount(data.row + 1);
    table->setItem(data.row, 0, new QTableWidgetItem(QString::number(data.time)));

    table->setItem(data.row, set_off+1, new QTableWidgetItem(QString::number(data.n_steps)));
    table->setItem(data.row, set_off+2, new QTableWidgetItem(QString::number(data.delta)));
    table->setItem(data.row, set_off+3, new QTableWidgetItem(QString::number(data.omega)));
    table->setItem(data.row, set_off+4, new QTableWidgetItem(QString::number(data.Eqe)));
    table->setItem(data.row, set_off+5, new QTableWidgetItem(QString::number(data.V)));
    table->setItem(data.row, set_off+6, new QTableWidgetItem(QString::number(data.U)));
    //table->setItem(row, set_off+6, new QTableWidgetItem(""));//QString::number(data.Eqeprime)));
}

void ToExcel::make_up(bool filled)
{
    if(filled)
        table->resizeColumnsToContents();
    else
        table->setRowCount(0);
}

void ToExcel::export_to_excel()
{
    const QString fn = QFileDialog::getSaveFileName(this,"Save to file","","Capable to import with Excel (*.csv)");
    QFile file(fn);
    if(not file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file);
    out.setGenerateByteOrderMark(true);

    const char delimiter = ';';

    out << table->horizontalHeaderItem(0)->text();
    for(int c=1; c<table->columnCount(); c++)
        out << delimiter << table->horizontalHeaderItem(c)->text();
    out << '\n';

    for(int r=0; r<table->rowCount(); r++)
    {
        out << table->item(r,0)->text();
        for(int c=1; c<table->columnCount(); c++)
            if(const QTableWidgetItem* item = table->item(r,c))
                out << delimiter << item->text().replace('.',',');
            else
                out << delimiter << ' ';
        out << '\n';
    }
}

