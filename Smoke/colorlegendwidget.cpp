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

    using namespace fluids;
    image = constructLegend(get_color_func((colormap)scalar_col), bands);

    image = image.scaled(bar.width(), bar.height(), Qt::IgnoreAspectRatio);
    QPainter painter(this);
    painter.fillRect(bar, image);
}


QImage ColorLegendWidget::constructLegend(color_func f, int banding_levels)
{
    float r, g, b;
    QImage image(1, 256, QImage::Format_RGB32);
    for (int i = 0; i < 256; i++) {
        // Also display max value in the color legend
        int div = 256;
        if (fluids::enable_bands) {
            div -= 255 / banding_levels;
            fluids::with_banding(f, (float)i / div, &r, &g, &b, banding_levels);
        } else {
            f((float)i / div, &r, &g, &b);
        }
        image.setPixel(0, 255 - i, qRgb(255*r, 255*g, 255*b));
    }
    return image;
}
