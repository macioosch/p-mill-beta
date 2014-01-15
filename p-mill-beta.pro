#-------------------------------------------------
#
# Project created by QtCreator 2013-12-16T15:35:06
#
#-------------------------------------------------
TARGET = p-mill-beta

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TEMPLATE = app

SOURCES += main.cpp \
    parse_cli_opts.cpp \
    ball_mill.cpp

LIBS += -lboost_program_options

HEADERS += \
    parse_cli_opts.h \
    parameters.h \
    ball_mill.h

QMAKE_CXXFLAGS += -std=c++0x
