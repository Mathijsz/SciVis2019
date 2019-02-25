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
