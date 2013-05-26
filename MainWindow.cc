#include "MainWindow.h"
#include "StartingPage.h"
#include "params.h"

MainWindow::MainWindow() : QMainWindow(NULL)
{
    setCentralWidget(new StartingPage);
    setFixedSize(1600,900);
}

