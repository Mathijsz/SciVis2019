#include "smokeglwidget.h"
#include "fluids.h"

#include <QTimer>
#include <QMouseEvent>
#include <QDebug>
#include <QRadioButton>
#include <QString>
#include <QRect>
#include <QPainter>
#include <QSpinBox>

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
    fluids::destroy_simulation();
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
    if (fluids::autoscale_colormaps) {
        update_max_box(fluids::max_col);
        update_min_box(fluids::min_col);
        trigger_colormap();
    }
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
        int digit = (button->objectName().end()-1)->digitValue();
        if (digit >= 0)
            fluids::scalar_col = digit;
        else
            qDebug() << "Invalid radio button sender for setting color map";
        trigger_colormap();
    }
}

void SmokeGLWidget::set_hedgehogs(int status)
{
    fluids::draw_vecs = status;
}

void SmokeGLWidget::set_smoke(int status)
{
    fluids::enable_smoke = status;
}

void SmokeGLWidget::set_bands(int bands)
{
    fluids::bands = bands;
}

void SmokeGLWidget::enable_bands(int status)
{
    fluids::enable_bands = (bool)status;
}

void SmokeGLWidget::set_minmax(double value)
{
    QSpinBox *sb = static_cast<QSpinBox *>(sender());
    if (sb->objectName().startsWith("min"))
        fluids::min_col = (float)value;
    else
        fluids::max_col = (float)value;
    trigger_colormap();
}

void SmokeGLWidget::set_color_data(bool toggle)
{
    if (toggle) {
        QRadioButton *button = static_cast<QRadioButton *>(sender());
        int digit = (button->objectName().end()-1)->digitValue();
        if (digit >= 0) {
            fluids::color_data_type = (vis_data_type)digit;
        } else
            qDebug() << "Invalid radio button sender for setting data type";
    }
}

void SmokeGLWidget::set_vector_data(bool toggle)
{
    if (toggle) {
        QRadioButton *button = static_cast<QRadioButton *>(sender());
        int digit = (button->objectName().end()-1)->digitValue();
        if (digit >= 0)
            fluids::vector_data_type = (vis_data_type)digit;
        else
            qDebug() << "Invalid radio button sender for setting data type";
    }
}

void SmokeGLWidget::set_dim_x(int x)
{
    fluids::DIM_X = x;
}

void SmokeGLWidget::set_dim_y(int y)
{
    fluids::DIM_Y = y;
}

void SmokeGLWidget::set_interpol_type(bool toggle)
{
    if (toggle) {
        QRadioButton *button = static_cast<QRadioButton *>(sender());
        int digit = (button->objectName().end()-1)->digitValue();
        if (digit >= 0)
            fluids::interpolation = (interpol_type)digit;
        else
            qDebug() << "Invalid radio button sender for setting interpolation type";
    }
}

void SmokeGLWidget::set_glyph_scale(int n)
{
    fluids::vec_scale = (float)n;
}

void SmokeGLWidget::set_glyph_shape(bool toggle)
{
    if (toggle) {
        QRadioButton *button = static_cast<QRadioButton *>(sender());
        int digit = (button->objectName().end()-1)->digitValue();
        if (digit >= 0)
            fluids::glyph_shape = (glyph_type)digit;
        else
            qDebug() << "Invalid radio button sender for setting glyph shape";
    }
}

void SmokeGLWidget::set_autoscale_colormap(int status)
{
    fluids::autoscale_colormaps = (bool)status;
}

void SmokeGLWidget::enable_isolines(int status)
{
    fluids::enable_isolines = (bool)status;
}

void SmokeGLWidget::set_isolines(double value)
{
    fluids::isoline = (float)value;
}
