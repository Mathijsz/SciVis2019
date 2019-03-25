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
    using namespace fluids;

    float r, g, b;
    int range = 256;
    QImage image(1, range, QImage::Format_RGB32);
    for (int i = 0; i < range; i++) {
        // Also display max value in the color legend
        int div = range;
        if (enable_bands)
            div -= 255 / banding_levels;
        float scaled_col = (max_col - min_col) * ((float)i/div) + min_col;
        fluids::with_banding(f, scaled_col, &r, &g, &b, banding_levels);
        image.setPixel(0, range - 1 - i, qRgb(255*r, 255*g, 255*b));
    }
    return image;
}
