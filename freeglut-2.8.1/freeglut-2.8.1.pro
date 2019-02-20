QT       -= core gui

TARGET = freeglut-2.8.1
TEMPLATE = lib
CONFIG += staticlib warn_off

INCLUDEPATH += include/

DEFINES += FREEGLUT_STATIC

win32: DEFINES += TARGET_HOST_MS_WINDOWS X_DISPLAY_MISSING

# Prehistoric code in fluid.c does not know unicode
DEFINES -= UNICODE

QMAKE_CFLAGS =
QMAKE_CFLAGS_RELEASE = -O3
DEF_FILE = src/freeglutdll.def

LIBS += -lopengl32 -lgdi32 -lwinmm
QMAKE_LFLAGS = -flto -O3

HEADERS += \
    include/GL/freeglut.h \
    include/GL/freeglut_ext.h \
    include/GL/freeglut_std.h \
    include/GL/glut.h \
    src/freeglut_teapot_data.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

SOURCES += \
    src/freeglut_callbacks.c \
    src/freeglut_cursor.c \
    src/freeglut_display.c \
    src/freeglut_ext.c \
    src/freeglut_font.c \
    src/freeglut_font_data.c \
    src/freeglut_gamemode.c \
    src/freeglut_geometry.c \
    src/freeglut_glutfont_definitions.c \
    src/freeglut_init.c \
    src/freeglut_input_devices.c \
    src/freeglut_joystick.c \
    src/freeglut_main.c \
    src/freeglut_menu.c \
    src/freeglut_misc.c \
    src/freeglut_overlay.c \
    src/freeglut_spaceball.c \
    src/freeglut_state.c \
    src/freeglut_stroke_mono_roman.c \
    src/freeglut_stroke_roman.c \
    src/freeglut_structure.c \
    src/freeglut_teapot.c \
    src/freeglut_videoresize.c \
    src/freeglut_window.c

!win32: SOURCES += src/freeglut_xinput.c
