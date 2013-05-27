#ifndef __STARTING_PAGE_H_
#define __STARTING_PAGE_H_

#include <QWidget>
#include "params.h"

class QDoubleValidator;
class QProgressBar;
class QPushButton;
class QLineEdit;
class QCheckBox;
class QwtPlot;
class QwtPlotCurve;
class QwtPlotZoomer;

class StartingPage : public QWidget
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
    QwtPlotCurve* curve_trapeze_Eqprime;
    QPushButton* start_button;

    QCheckBox* enable_eiler;
    QCheckBox* enable_trapeze;
    QProgressBar* progress_bar_eiler;
    QProgressBar* progress_bar_trapeze;

    QLineEdit *K0f, *K1f, *Ur0, *K0u, *K1u;
    QLineEdit *Y11, *Y12, *A11, *A12, *Y11em, *Y12em, *A11em, *A12em, *Pd;
    QLineEdit *Delta0, *Eqe0, *Eqprime0, *U0, *V0;
    QLineEdit *Tstart, *Tstop, *dt;
    QLineEdit *eps;
    QLineEdit *max_iterations;
    QDoubleValidator* validator;

    Params collect_params();

  private slots:
    void start();
    void enable_everything(bool);
    void increment_eiler();
    void increment_trapeze();
    void wanna_eiler();
    void wanna_trapeze();
    void plot_answer_eiler(const QVector<AnswerItem> ans);
    void plot_answer_trapeze(const QVector<AnswerItem> ans);
  signals:
  public:
    StartingPage();
    ~StartingPage();
};

#endif // __STARTING_PAGE_H_
