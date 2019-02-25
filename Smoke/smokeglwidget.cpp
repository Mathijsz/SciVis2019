#include "smokeglwidget.h"
#include "fluids.h"

#include <QTimer>
#include <QMouseEvent>
#include <QDebug>

SmokeGLWidget::SmokeGLWidget(QWidget *parent)
 : QOpenGLWidget(parent), timer(new QTimer(nullptr))
{
    fluids::init_simulation(50);
    connect(timer, SIGNAL(timeout()), this, SLOT(step()));
    timer->start(17);
}

SmokeGLWidget::~SmokeGLWidget()
{
    delete timer;
}

void SmokeGLWidget::initializeGL()
{
    initializeOpenGLFunctions();
}

void SmokeGLWidget::resizeGL(int w, int h)
{
    fluids::reshape(w, h);
}

void SmokeGLWidget::paintGL()
{
    fluids::display();
}

void SmokeGLWidget::step()
{
    fluids::do_one_simulation_step();
    update();
}

void SmokeGLWidget::mouseMoveEvent(QMouseEvent *e)
{
    if (e->buttons() == Qt::LeftButton)
        fluids::drag(e->x(), e->y());
}

void SmokeGLWidget::keyPressEvent(QKeyEvent *e)
{
    if (e->text().length() <= 0)
        return;
    // keyboard function in the original code only deals with ascii
    qDebug() << e->text();
    QString t = e->text();
    char c = t.at(0).toLatin1();
    fluids::keyboard(static_cast<unsigned char>(c), 0, 0);
}

void SmokeGLWidget::set_color(int status)
{
    fluids::color_dir = status;
}