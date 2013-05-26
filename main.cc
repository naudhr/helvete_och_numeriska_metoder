#include <QApplication>
#include "MainWindow.h"

int main(int argc, char *argv[])
{
  QApplication qapp(argc,argv);

  MainWindow w;
  QObject::connect(&w,SIGNAL(quit()),&qapp,SLOT(quit()));

  w.show();
  //w.showMaximized();
  return qapp.exec();
}

