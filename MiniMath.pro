#-------------------------------------------------
#
# Project created by QtCreator 2014-12-23T14:07:44
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = MiniMath
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    menus.cpp \
    functions.cpp \
    algorithms.cpp \
    codewindow.cpp

HEADERS  += mainwindow.h \
    functions.h \
    algorithms.h \
    menus.h \
    codewindow.h

FORMS    += mainwindow.ui

CONFIG += mobility
MOBILITY = 

OTHER_FILES +=

RESOURCES += \
    text.qrc

