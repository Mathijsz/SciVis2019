TEMPLATE = subdirs

CONFIG += ordered

SUBDIRS = fftw-2.1.5 \
    Smoke

Smoke.depends = fftw-2.1.5
