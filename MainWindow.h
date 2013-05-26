#ifndef __MAIN_WINDOW_H_
#define __MAIN_WINDOW_H_

#include <QMainWindow>

class QTabWidget;
struct Params;

class MainWindow : public QMainWindow
{
    Q_OBJECT

    //QTabWidget* tabs;
  private slots:
    //void starting_page();
    //void show_answer(const Answer& );
      /*
    void options();
    void echoMail();
    void echoState(int,int,int,int);
    void echoConnected(bool);
    void echoReceived();
    void echoSent();
    void s_check_mail();
    void show_all();
    */

  public slots:
    //void applySettings();

  signals:
    void quit();
    //void how_you_doin();

  public:
    MainWindow();
};

#endif // __MAIN_WINDOW_H_
