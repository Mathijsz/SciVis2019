TEMPLATE = subdirs

CONFIG += ordered

SUBDIRS = fftw-2.1.5

Smoke.depends = fftw-2.1.5

win32 {
    SUBDIRS += freeglut-2.8.1
    Smoke.depends += freeglut-2.8.1
}

SUBDIRS += Smoke
