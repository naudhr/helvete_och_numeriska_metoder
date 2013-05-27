QT       += core gui

TEMPLATE = app
TARGET = honm
DEPENDPATH += .
INCLUDEPATH += .

CONFIG += qwt

# Input
HEADERS += MainWindow.h StartingPage.h params.h Calculus.h
SOURCES += main.cc MainWindow.cc StartingPage.cc Calculus.cc
