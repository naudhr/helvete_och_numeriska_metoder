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
class QwtPlot;
class QwtPlotCurve;
class QwtPlotZoomer;

class CalculusWidget : public QWidget
{
    Q_OBJECT

    QwtPlotZoomer* zoom;

    QwtPlot* plot;
    QwtPlotCurve* curve_eiler_U;
    QwtPlotCurve* curve_eiler_Delta;
    QwtPlotCurve* curve_eiler_Omega;
    QwtPlotCurve* curve_eiler_Eqe;
    QwtPlotCurve* curve_eiler_Eqprime;
    QwtPlotCurve* curve_trapeze_Delta;
    QwtPlotCurve* curve_trapeze_Omega;
    QwtPlotCurve* curve_trapeze_Eqe;
    QwtPlotCurve* curve_trapeze_U;
    QwtPlotCurve* curve_trapeze_Eqprime;

    QCheckBox* enable_eiler;
    QCheckBox* enable_trapeze;
    QProgressBar* progress_bar_eiler;
    QProgressBar* progress_bar_trapeze;

    QLineEdit *Y11, *Y12, *A11, *A12, *Y11em, *Y12em, *A11em, *A12em, *Pd;
    QLineEdit *Delta0, *Eqe0, *Eqprime0, *U0, *V0;
    QLineEdit *Tstart, *Tstop, *dt;
    QLineEdit *eps;
    QLineEdit *max_iterations;

    Params collect_params();

  private slots:
    void enable_everything(bool);
    void increment_eiler();
    void increment_trapeze();
    void wanna_eiler();
    void wanna_trapeze();
    void plot_answer_eiler(const QVector<AnswerItem> ans);
    void plot_answer_trapeze(const QVector<AnswerItem> ans);
  signals:
    void enable_start_button(bool);
    void to_excel_populate(double, QVector<x2_U_D_E>);
  public:
    CalculusWidget(QWidget* );
    ~CalculusWidget();
    void start(const Params::Consts& );
};

class SystemParamsWidget : public QWidget
{
    QLineEdit *T, *Ty, *Tu, *Tg, *Te, *Tf, *Tphi, *Td0;
    QLineEdit *Pt0, *omega_nom, *Xdprime, *Xd, *Eqenom, *Uc;
    QLineEdit *K0f, *K1f, *Ur0, *K0u, *K1u;
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
    void populate(double dt, const QVector<x2_U_D_E>& data);
};

#endif // __MAIN_WINDOW_H_
