#ifndef SMOKEGLWIDGET_H
#define SMOKEGLWIDGET_H

#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QVector2D>

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
    void set_repeat_levels(int repeats);
    void enable_repeats(int status);
    void set_minmax(double value);
    void set_color_data(bool toggle);
    void set_vector_data(bool toggle);

    void set_dim_x(int x);
    void set_dim_y(int y);
    void set_interpol_type(bool toggle);

    void set_glyph_scale(int n);
    void set_glyph_shape(bool toggle);
    void set_autoscale_colormap(int status);
    void enable_isolines(int status);
    void set_isolines(double value);

    void set_upper_isoline(double value);
    void set_lower_isoline(double value);
    void set_isoline_count(int value);
    void enable_bounded_isolines(int status);

    void enable_heightmap(int status);
    void set_height_data(bool toggle);

    void toggle_shading(int status);
    void reset_view();

    void set_height(int h);

    void enable_streamtubes(int status);
    void reset_streamtubes();

signals:
    void trigger_colormap();
    void update_min_box(double min);
    void update_max_box(double max);
    void set_streamtube_count(int count);

protected:
    void initializeGL() override;

    void resizeGL(int w, int h) override;

    void paintGL() override;

    void mouseMoveEvent(QMouseEvent *event) override;

    void mousePressEvent(QMouseEvent *event) override;

    void wheelEvent(QWheelEvent *event) override;

    void keyPressEvent(QKeyEvent *event) override;

private:
    QTimer *timer;
    QVector2D last_mouse_pos;
    bool new_input;

};

#endif // SMOKEGLWIDGET_H
