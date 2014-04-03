#ifndef __MAIN_WINDOW_H_
#define __MAIN_WINDOW_H_

#include <QMainWindow>
#include "params.h"


class QTabWidget;
class QPushButton;
class QTableWidget;
class CalculusWidget;
class SystemParamsWidget;
class ToExcel;

class MainWindow : public QMainWindow
{
    Q_OBJECT
  public:
    MainWindow();
  signals:
    void quit();
};

class CentralWidget : public QWidget
{
    Q_OBJECT

    QTabWidget* tabs;
    ToExcel* to_excel;
    CalculusWidget* calculus;
    SystemParamsWidget* sysparams;
    QPushButton* start_button;

  private slots:
    void start();

  public:
    CentralWidget();
};

class QDoubleValidator;
class QProgressBar;
class QLineEdit;
class QCheckBox;
class NoQwtPlot;
class NoQwtGraphicsView;

class CalculusWidget : public QWidget
{
    Q_OBJECT

    NoQwtGraphicsView* view;
    NoQwtPlot* plot;

    QCheckBox* dirty_hack;
    QCheckBox* online_plotting;
    QCheckBox* enable_eiler;
    QCheckBox* enable_trapeze;
    QCheckBox* enable_sequensive;
    QCheckBox* enable_parallel;
    QProgressBar* progress_bar_eiler;
    QProgressBar* progress_bar_trapeze;
    QProgressBar* progress_bar_sequensive;
    QProgressBar* progress_bar_parallel;

    QLineEdit *Y11, *Y12, *A11, *A12, *Y11em, *Y12em, *A11em, *A12em, *Pd;
    QLineEdit *Delta0, *Eqe0, *Eqprime0, *U0, *V0;
    QLineEdit *Tstart, *Tstop, *dt;
    QLineEdit *eps;
    QLineEdit *max_iterations;
    QLineEdit *seq_x1_0, *seq_x2_0, *seq_x4_0, *seq_x5_0, *seq_x6_0, *seq_x7_0, *seq_x9_0, *seq_x0_0;
    QWidget* seq_params;
    QPushButton* power_part;

    Params collected_params;
    void collect_params();

    size_t jobs;
    QVector<AnswerItem> answer_buffer;

  signals:
    void enable_start_button(bool);
    void to_excel_populate(const AnswerItem& ans);

  private slots:
    void enable_everything(bool);
    void eiler_step(const AnswerItem& );
    void trapeze_step(const AnswerItem& );
    void sequensive_step(const AnswerItem& );
    void parallel_step(const AnswerItem& );
    void some_calc_enabled();
    void a_part_of_the_plot_done();
    void ndnc(QString name, double t, unsigned n_steps);
    void popup_power_widget();
    void sysparams_changed();

  public slots:
    void sysparams_changed(Params::Consts);

  public:
    CalculusWidget(QWidget* );
    ~CalculusWidget();
    void start(const Params::Consts& );
};

class SystemParamsWidget : public QWidget
{
    Q_OBJECT
    
    QLineEdit *T, *Ty, *Tu, *Tg, *Te, *Tf, *Tphi, *Td0;
    QLineEdit *Pt0, *omega_nom, *Xdprime, *Xd, *Eqenom, *Uc;
    QLineEdit *K0f, *K1f, *Ur0, *K0u, *K1u;
  private slots:
    void emit_params_changed();
  signals:
    void params_changed(Params::Consts);
  public:
    SystemParamsWidget(QWidget* );
    Params::Consts collect_params();
};

class ToExcel : public QWidget
{
    Q_OBJECT

    QTableWidget* table;
  public:
    ToExcel(QWidget* );
  private slots:
    void export_to_excel();
  public slots:
    void populate(const AnswerItem& data);
    void make_up(bool filled);
};

#endif // __MAIN_WINDOW_H_
