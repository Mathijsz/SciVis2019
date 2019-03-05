#ifndef COLORLEGENDWIDGET_H
#define COLORLEGENDWIDGET_H

#include <QWidget>

class ColorLegendWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ColorLegendWidget(QWidget *parent = 0);

protected:
    void paintEvent(QPaintEvent *event);

signals:

private:
    QImage constructLegend(void (*f)(float, float*, float*, float*), int banding_levels = 0);

};

#endif // COLORLEGENDWIDGET_H
