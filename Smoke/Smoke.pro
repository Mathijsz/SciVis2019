#-------------------------------------------------
#
# Project created by QtCreator 2019-02-11T18:47:59
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Smoke
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

# Needed for freeglut
DEFINES += TARGET_HOST_MS_WINDOWS X_DISPLAY_MISSING FREEGLUT_STATIC

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    ../fluids.c

HEADERS += \
        mainwindow.h

FORMS += \
        mainwindow.ui

LIBS += -lm

INCLUDEPATH += "$$PWD/../freeglut-2.8.1/include"

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

# FFTW

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/release/ -lrfftw
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/debug/ -lrfftw
else:unix: LIBS += -L$$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/ -lrfftw

INCLUDEPATH += $$PWD/../fftw-2.1.5/sourceAndDoc/rfftw
DEPENDPATH += $$PWD/../fftw-2.1.5/sourceAndDoc/rfftw

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/release/librfftw.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/debug/librfftw.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/release/rfftw.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/debug/rfftw.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/rfftw/librfftw.a

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/release/ -lfftw
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/debug/ -lfftw
else:unix: LIBS += -L$$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/ -lfftw

INCLUDEPATH += $$PWD/../fftw-2.1.5/sourceAndDoc/fftw
DEPENDPATH += $$PWD/../fftw-2.1.5/sourceAndDoc/fftw

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/release/libfftw.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/debug/libfftw.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/release/fftw.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/debug/fftw.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../fftw-2.1.5/sourceAndDoc/fftw/libfftw.a

# GLUT

win32 {
    CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../freeglut-2.8.1/release/ -lfreeglut-2.8.1
    else:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../freeglut-2.8.1/debug/ -lfreeglut-2.8.1

    INCLUDEPATH += $$PWD/../freeglut-2.8.1
    DEPENDPATH += $$PWD/../freeglut-2.8.1

    win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../freeglut-2.8.1/release/libfreeglut-2.8.1.a
    else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../freeglut-2.8.1/debug/libfreeglut-2.8.1.a
    else:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../freeglut-2.8.1/release/freeglut-2.8.1.lib
    else:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../freeglut-2.8.1/debug/freeglut-2.8.1.lib

    LIBS += -lopengl32 -lglu32 -lgdi32 -lwinmm

}
else:LIBS += -lglut -lGL -lGLU -lGLEW
