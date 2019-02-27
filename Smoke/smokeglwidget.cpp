#include "smokeglwidget.h"
#include "fluids.h"

#include <QTimer>
#include <QMouseEvent>
#include <QDebug>
#include <QRadioButton>
#include <QString>
#include <QRect>
#include <QPainter>

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
    qDebug() << "Using OpenGL"
             << reinterpret_cast<const char *>(glGetString(GL_VERSION))
             << "on"
             << reinterpret_cast<const char *>(glGetString(GL_RENDERER));
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
    QString t = e->text();
    char c = t.at(0).toLatin1();
    fluids::keyboard(static_cast<unsigned char>(c), 0, 0);
}

void SmokeGLWidget::set_color(int status)
{
    fluids::color_dir = status;
}

void SmokeGLWidget::change_color(bool toggle)
{
    if (toggle) {
        QRadioButton *button = static_cast<QRadioButton *>(sender());
        fluids::scalar_col = (button->objectName().end()-1)->digitValue();
        trigger_colormap();
    }
}

void SmokeGLWidget::set_hedgehogs(int status)
{
    fluids::draw_vecs = status;
}

void SmokeGLWidget::set_smoke(int status)
{
    fluids::draw_smoke = status;
}
