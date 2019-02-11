QT       -= core gui

TARGET = rfftw
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

INCLUDEPATH += ../fftw

HEADERS += \
    rfftw.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}

SOURCES += \
    fcr_1.c \
    fcr_2.c \
    fcr_3.c \
    fcr_4.c \
    fcr_5.c \
    fcr_6.c \
    fcr_7.c \
    fcr_8.c \
    fcr_9.c \
    fcr_10.c \
    fcr_11.c \
    fcr_12.c \
    fcr_13.c \
    fcr_14.c \
    fcr_15.c \
    fcr_16.c \
    fcr_32.c \
    fcr_64.c \
    fcr_128.c \
    fhb_2.c \
    fhb_3.c \
    fhb_4.c \
    fhb_5.c \
    fhb_6.c \
    fhb_7.c \
    fhb_8.c \
    fhb_9.c \
    fhb_10.c \
    fhb_16.c \
    fhb_32.c \
    fhf_2.c \
    fhf_3.c \
    fhf_4.c \
    fhf_5.c \
    fhf_6.c \
    fhf_7.c \
    fhf_8.c \
    fhf_9.c \
    fhf_10.c \
    fhf_16.c \
    fhf_32.c \
    frc_1.c \
    frc_2.c \
    frc_3.c \
    frc_4.c \
    frc_5.c \
    frc_6.c \
    frc_7.c \
    frc_8.c \
    frc_9.c \
    frc_10.c \
    frc_11.c \
    frc_12.c \
    frc_13.c \
    frc_14.c \
    frc_15.c \
    frc_16.c \
    frc_32.c \
    frc_64.c \
    frc_128.c \
    rconfig.c \
    rexec.c \
    rexec2.c \
    rfftwf77.c \
    rfftwnd.c \
    rgeneric.c \
    rplanner.c
