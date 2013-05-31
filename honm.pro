QT       += core gui

TEMPLATE = app
TARGET = honm
DEPENDPATH += .
INCLUDEPATH += .

CONFIG += qwt release
CONFIG -= debug


# Input
HEADERS += MainWindow.h params.h Calculus.h
SOURCES += MainWindow.cc Calculus.cc
