#include "colorlegendwidget.h"

#include <QPainter>
#include <QImage>
#include <QDebug>
#include "fluids.h"

ColorLegendWidget::ColorLegendWidget(QWidget *parent) : QWidget(parent)
{

}

void ColorLegendWidget::paintEvent(QPaintEvent *event)
{
    QRect bar(0, 0, width(), height());
    QImage image;

    switch (fluids::scalar_col) {
        case COLOR_RAINBOW:
            image = constructLegend(&fluids::rainbow);
            break;
        case COLOR_BLUE_TO_YELLOW:
            image = constructLegend(&fluids::blue_to_yellow);
            break;
        case COLOR_RED_TO_WHITE:
            image = constructLegend(&fluids::red_to_white);
            break;
        case COLOR_BLUE_TO_RED_VIA_WHITE:
            image = constructLegend(&fluids::blue_to_red_via_white);
            break;
        case COLOR_BANDS:
            image = constructLegend(&fluids::rainbow, 10);
            break;
        default:
            image = constructLegend(&fluids::white_to_black);
            break;
    }
    image = image.scaled(bar.width(), bar.height(), Qt::IgnoreAspectRatio);
    QPainter painter(this);
    painter.fillRect(bar, image);
}


QImage ColorLegendWidget::constructLegend(void (*f)(float, float*, float*, float*), int banding_levels)
{
    float r, g, b;
    QImage image(1, 256, QImage::Format_RGB32);
    for (int i = 0; i < 256; i++) {
        fluids::with_banding(f, (float)i / 255, &r, &g, &b, banding_levels);
        image.setPixel(0, i, qRgb(255*r, 255*g, 255*b));
    }
    return image;
}
