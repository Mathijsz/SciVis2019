#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication smoke(argc, argv);
    MainWindow w;
    w.show();

    return smoke.exec();
}
