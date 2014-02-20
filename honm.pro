QT       += core gui

TEMPLATE = app
TARGET = honm
DEPENDPATH += .
INCLUDEPATH += .

CONFIG += qwt release qaxcontainer
CONFIG -= debug


# Input
HEADERS += MainWindow.h params.h Calculus.h NoQwt.h
SOURCES += MainWindow.cc Calculus.cc NoQwt.cc
