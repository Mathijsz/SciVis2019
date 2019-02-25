#include "colorlegendwidget.h"

#include <QPainter>
#include "fluids.h"

ColorLegendWidget::ColorLegendWidget(QWidget *parent) : QWidget(parent)
{

}

void ColorLegendWidget::paintEvent(QPaintEvent *event)
{
    QRect bar(0, 0, width(), height());
    QLinearGradient gradient(bar.topLeft(), bar.bottomLeft());

    switch (fluids::scalar_col) {
        case COLOR_RAINBOW:
            gradient.setColorAt(0, Qt::blue);
            gradient.setColorAt(0.5, Qt::green);
            gradient.setColorAt(1, Qt::red);
            break;
        case COLOR_BLUE_TO_YELLOW:
            gradient.setColorAt(0, Qt::blue);
            gradient.setColorAt(1, Qt::yellow);
            break;
        case COLOR_RED_TO_WHITE:
            gradient.setColorAt(0, Qt::red);
            gradient.setColorAt(1, Qt::white);
            break;
        case COLOR_BLUE_TO_RED_VIA_WHITE:
            gradient.setColorAt(0, Qt::blue);
            gradient.setColorAt(0.5, Qt::white);
            gradient.setColorAt(1, Qt::red);
            break;
        case COLOR_BANDS:

            break;
        default:

            break;
    }

    QPainter painter(this);
    painter.fillRect(bar, gradient);
}
