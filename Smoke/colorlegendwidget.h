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

};

#endif // COLORLEGENDWIDGET_H
