QT       -= core gui

TARGET = fftw
TEMPLATE = lib
CONFIG += staticlib

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

HEADERS += \
    fftw-int.h \
    fftw.h \
    config.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}

SOURCES += \
    wisdomio.c \
    wisdom.c \
    twiddle.c \
    timer.c \
    rader.c \
    putils.c \
    planner.c \
    malloc.c \
    generic.c \
    ftwi_64.c \
    ftwi_32.c \
    ftwi_16.c \
    ftwi_10.c \
    ftwi_9.c \
    ftwi_8.c \
    ftwi_7.c \
    ftwi_6.c \
    ftwi_5.c \
    ftwi_4.c \
    ftwi_3.c \
    ftwi_2.c \
    ftw_64.c \
    ftw_32.c \
    ftw_16.c \
    ftw_10.c \
    ftw_9.c \
    ftw_8.c \
    ftw_7.c \
    ftw_6.c \
    ftw_5.c \
    ftw_4.c \
    ftw_3.c \
    ftw_2.c \
    fni_64.c \
    fni_32.c \
    fni_16.c \
    fni_15.c \
    fni_14.c \
    fni_13.c \
    fni_12.c \
    fni_11.c \
    fni_10.c \
    fni_9.c \
    fni_8.c \
    fni_7.c \
    fni_6.c \
    fni_5.c \
    fni_4.c \
    fni_3.c \
    fni_2.c \
    fni_1.c \
    fn_64.c \
    fn_32.c \
    fn_16.c \
    fn_15.c \
    fn_14.c \
    fn_13.c \
    fn_12.c \
    fn_11.c \
    fn_10.c \
    fn_9.c \
    fn_8.c \
    fn_7.c \
    fn_6.c \
    fn_5.c \
    fn_4.c \
    fn_3.c \
    fn_2.c \
    fn_1.c \
    fftwnd.c \
    fftwf77.c \
    executor.c \
    config.c
