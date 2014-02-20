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

int main(int argc, char *argv[])
{
  QApplication qapp(argc,argv);

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
    connect(calculus, SIGNAL(to_excel_populate(AnswerItem,int)), to_excel, SLOT(populate(AnswerItem,int)));
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

static void add_plot_curve(NoQwtPlot* plot, const QString& label, int r, int g, int b, const QString& tag)
{
    NoQwtPlotCurve* curve = new NoQwtPlotCurve(plot, label, QPen(QColor(r,g,b)), QBrush(Qt::black), tag);
    curve->setVisible(true);
    //return curve;
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

    QHBoxLayout* srt = new QHBoxLayout;
    srt->addWidget(new QLabel("Initial values:",this));
    Delta0 = add_double_input(srt, "Delta0:", "1.294", validator);
    Eqe0 = add_double_input(srt, "Eqe0:", "1.698", validator);
    Eqprime0 = add_double_input(srt, "Eqprime0:", "0.759", validator);
    U0 = add_double_input(srt, "U0:", "1.035", validator);
    V0 = add_double_input(srt, "V0:", "0.353", validator);

    QHBoxLayout* cpt = new QHBoxLayout;
    cpt->addWidget(new QLabel("Computational parameters:",this));
    Tstart = add_double_input(cpt, "Tstart:", "0", validator);
    Tstop = add_double_input(cpt, "Tstop:", "3", validator);
    dt = add_double_input(cpt, "dt:", "0.001", validator);
    eps = add_double_input(cpt, "eps:", "0.005", validator);
    max_iterations = add_double_input(cpt, "max_iterations:", "50", validator);
    max_iterations->setValidator(new QIntValidator(5,5000,this));

    QHBoxLayout* prb = new QHBoxLayout;
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
    progress_bar_parallel = new QProgressBar(this);
    progress_bar_parallel->setFormat("Parallel");
    progress_bar_parallel->setMinimum(0);
    progress_bar_parallel->setMaximum(1);
    progress_bar_parallel->setValue(0);
    progress_bar_parallel->setEnabled(false);
    prb->addWidget(enable_parallel);
    prb->addWidget(progress_bar_parallel);

    view = new NoQwtGraphicsView(this);
    plot = new NoQwtPlot(view->scene(), "", "", "");

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

    QVBoxLayout* l = new QVBoxLayout(this);
    l->addLayout(eqv);
    l->addLayout(srt);
    l->addLayout(cpt);
    l->addLayout(prb);
    l->addWidget(view,2);

    plot->setVisible(true);
    progress_bar_eiler->setVisible(true);
    progress_bar_trapeze->setVisible(true);

    connect(enable_eiler, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    connect(enable_trapeze, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    connect(enable_sequensive, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    connect(enable_parallel, SIGNAL(clicked()), SLOT(some_calc_enabled()));
    jobs = 0;
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
}

void CalculusWidget::start(const Params::Consts& reg)
{
    enable_everything(false);

    Params p = collect_params();
    p.reg = reg;
    const int n = (p.Tstop - p.Tstart)/p.dt;

    if( enable_eiler->isChecked() )
    {
        progress_bar_eiler->setMaximum(n);
        progress_bar_eiler->setValue(0);
        progress_bar_eiler->setFormat("Eiler: %v");

        QThread* c = new CalculusEiler(p);
        connect(c, SIGNAL(a_step_done(AnswerItem)), this, SLOT(eiler_step(AnswerItem)), Qt::QueuedConnection);
        connect(c, SIGNAL(finished()), this, SLOT(a_part_of_the_plot_done()), Qt::QueuedConnection);
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
        connect(c, SIGNAL(finished()), c, SLOT(deleteLater()));
        c->start();
        jobs ++;
    }

    if( enable_sequensive->isChecked() )
    {
        progress_bar_sequensive->setMaximum(n);
        progress_bar_sequensive->setValue(0);
        progress_bar_sequensive->setFormat("Sequensive: %v");

        QThread* c = new CalculusSequensive(p);
        connect(c, SIGNAL(a_step_done(AnswerItem)), this, SLOT(sequensive_step(AnswerItem)), Qt::QueuedConnection);
        connect(c, SIGNAL(finished()), this, SLOT(a_part_of_the_plot_done()), Qt::QueuedConnection);
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
        connect(c, SIGNAL(finished()), c, SLOT(deleteLater()));
        c->start();
        jobs ++;
    }

    plot->reset();
    plot->setVisible(true);

    enable_everything(jobs == 0);
    emit enable_start_button(jobs == 0);
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
    Params p;// = {};
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
    emit to_excel_populate(ans,2);
}

void CalculusWidget::parallel_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("P"), ans);
    progress_bar_parallel->setValue(progress_bar_parallel->value()+1);
    emit to_excel_populate(ans,3);
}

void CalculusWidget::trapeze_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("T"), ans);
    progress_bar_trapeze->setValue(progress_bar_trapeze->value()+1);
    emit to_excel_populate(ans,1);
}

void CalculusWidget::eiler_step(const AnswerItem& ans)
{
    plot_answer_step(plot, QLatin1String("E"), ans);
    progress_bar_eiler->setValue(progress_bar_eiler->value()+1);
    emit to_excel_populate(ans,0);
}

void CalculusWidget::a_part_of_the_plot_done()
{
    jobs --;
    view->fitInView(plot);
    //scene->setSceneRect(plot->boundingRect());
    view->update();

    enable_everything(jobs == 0);
    emit enable_start_button(jobs == 0);
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

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

ToExcel::ToExcel(QWidget* p) : QWidget(p), table(NULL)
{
    table = new QTableWidget(this);
    table->setColumnCount(21/*25*/);
    QStringList vh;
    vh << "Time" << (QChar(0x03b4)+QLatin1String(" Eiler")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Eiler")) << "Eqe Eiler"
                 << (QChar(0x03bd)+QLatin1String(" Eiler")) << "U Eiler" /*<< "E'q Eiler"*/
                 << (QChar(0x03b4)+QLatin1String(" Trapeze")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Trapeze")) << "Eqe Trapeze"
                 << (QChar(0x03bd)+QLatin1String(" Trapeze")) << "U Trapeze" /*<< "E'q  Trapeze"*/
                 << (QChar(0x03b4)+QLatin1String(" Sequensive")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Sequensive")) << "Eqe Sequensive"
                 << (QChar(0x03bd)+QLatin1String(" Sequensive")) << "U Sequensive" /*<< "E'q  Sequensive"*/
                 << (QChar(0x03b4)+QLatin1String(" Parallel")) << (QChar(0x0394)+QString(0x03c9)+QLatin1String(" Parallel")) << "Eqe Parallel"
                 << (QChar(0x03bd)+QLatin1String(" Parallel")) << "U Parallel" /*<< "E'q  Parallel"*/;
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

void ToExcel::populate(const AnswerItem& data, int set_no)
{
    const int set_size = (table->columnCount()-1)/4;
    const int set_off = set_size * set_no;

    int row = 0;
    for( ; row<table->rowCount(); row++)
        if(table->item(row,1+set_off) == NULL/*->text().isEmpty()*/)
            break;
    if(row >= table->rowCount())
    {
        table->setRowCount(row + 1);
        table->setItem(row, 0, new QTableWidgetItem(QString::number(data.time)));
    }
    table->setItem(row, set_off+1, new QTableWidgetItem(QString::number(data.delta)));
    table->setItem(row, set_off+2, new QTableWidgetItem(QString::number(data.omega)));
    table->setItem(row, set_off+3, new QTableWidgetItem(QString::number(data.Eqe)));
    table->setItem(row, set_off+4, new QTableWidgetItem(QString::number(data.V)));
    table->setItem(row, set_off+5, new QTableWidgetItem(QString::number(data.U)));
    //table->setItem(row, set_off+6, new QTableWidgetItem(""));//QString::number(data.Eqeprime)));
}

void ToExcel::clear()
{
    table->setRowCount(0);
}

void ToExcel::export_to_excel()
{
    const QString fn = QFileDialog::getSaveFileName(this,"Save to file","","Capable to import with Excel (*.csv)");
    QFile file(fn);
    if(not file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file);
    const char delimiter = ';';

    out << table->horizontalHeaderItem(0)->text() << delimiter
        << table->horizontalHeaderItem(1)->text() << delimiter
        << table->horizontalHeaderItem(2)->text() << delimiter
        << table->horizontalHeaderItem(3)->text() << delimiter
        << table->horizontalHeaderItem(4)->text() << delimiter
        << table->horizontalHeaderItem(5)->text() << delimiter
        << table->horizontalHeaderItem(6)->text() << '\n';

    for(int r=0; r<table->rowCount(); r++)
        out << table->item(r,0)->text() << delimiter
            << table->item(r,1)->text() << delimiter
            << table->item(r,2)->text() << delimiter
            << table->item(r,3)->text() << delimiter
            << table->item(r,4)->text() << delimiter
            << table->item(r,5)->text() << delimiter
            << table->item(r,6)->text() << '\n';
}

