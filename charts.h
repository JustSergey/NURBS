#ifndef CHARTS_H
#define CHARTS_H

#include "functions.h"
#include "ui_widget.h"

void plot_deriv_for_point(const QVector<Point_curve>& data_NURBS, Ui::Widget* ui)
{
    for(const auto& point: data_NURBS)
    {
        QCPItemLine *line = new QCPItemLine(ui->graph_function);
        line->setHead(QCPLineEnding::esFlatArrow);
        line->start->setCoords(point.curve.first, point.curve.second);
        line->end->setCoords(point.derivative_1.first + point.curve.first, point.derivative_1.second + point.curve.second);
    }

    ui->graph_function->replot();
}

#endif // CHARTS_H
