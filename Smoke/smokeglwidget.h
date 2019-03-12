#ifndef SMOKEGLWIDGET_H
#define SMOKEGLWIDGET_H

#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>


class SmokeGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
Q_OBJECT

public:
    SmokeGLWidget(QWidget *parent);
    ~SmokeGLWidget();

public slots:
    void step();

    void set_color(int status);
    void change_color(bool toggle);
    void set_hedgehogs(int status);
    void set_smoke(int status);
    void set_bands(int bands);
    void enable_bands(int status);
    void set_minmax(double value);
    void set_color_data(bool toggle);
    void set_vector_data(bool toggle);

    void set_dim_x(int x);
    void set_dim_y(int y);
    void set_interpol_type(bool toggle);

signals:
    void trigger_colormap();

protected:
    void initializeGL() override;

    void resizeGL(int w, int h) override;

    void paintGL() override;

    void mouseMoveEvent(QMouseEvent *event) override;

    void keyPressEvent(QKeyEvent *event) override;

private:
    QTimer *timer;

};

#endif // SMOKEGLWIDGET_H
