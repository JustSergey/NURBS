#ifndef CHART_H
#define CHART_H

#include "nurbs.h"
#include "ui_widget.h"

class Chart
{
public:
    explicit Chart(QCustomPlot *canvas, const QString &title = "");

    // Рисует линию (касательную), начинающуюся от точки кривой (point)
    void paint_tangent(const Point_curve &point_curve) const;

    // Рисует определяющий многоугольник с вершинами
    void paint_polygon(const QVector<QVector<double>> &polygon_points, const QString &label, const QColor &color = QColor(0, 0, 0, 255), double width = 1) const;

    // Рисует точку
    void paint_point(double x, double y, double width = 5, const QString &text = "", const QColor &color = QColor(0, 0, 0, 255)) const;

    // Рисует линию по двум точкам
    void paint_line(double x_1, double y_1, double x_2, double y_2, const QColor& color = QColor(0, 0, 0), double width = 3.5) const;

    // Рисует двойную стрелку "<-->" по двум точкам
    void paint_double_arrow(double x_1, double y_1, double x_2, double y_2) const;

    // Рисует стрелку "-->" по двум точкам
    void paint_arrow(double x_1, double y_1, double x_2, double y_2) const;

    // Рисует линию (касательную) c центром в точке кривой (point)
    void paint_tangent_centred(const Point_curve &point, const QColor &color = QColor(0, 100, 0), double width = 2.8) const;

    // Рисует кривую NURBS
    void paint_curve(const NURBS &points_NURBS, const QString &label, const Qt::PenStyle &penStyle = Qt::PenStyle::SolidLine, const QColor &color = QColor(0, 0, 0, 255), double width = 1.5) const;

    // Рисует кривую NURBS
    void paint_curve(const QVector<Point_curve> &points_NURBS, const QString &label, const Qt::PenStyle &penStyle = Qt::PenStyle::SolidLine, const QColor &color = QColor(0, 0, 0, 255), double width = 1.5) const;

    // Рисует надпись
    void paint_lable(double x, double y, const QString &text, double font_size = 10) const;

    // Рисует надпись со стрелкой
    void paint_lable_with_arrow(double x_start, double y_start, double x_end, double y_end, const QString &text) const;

    // Рисует вектор первой производной кривой
    void paint_first_deriv(const NURBS &points_NURBS, const QString &label, const Qt::PenStyle &penStyle = Qt::PenStyle::SolidLine, const QColor &color = QColor(0, 0, 0, 255), double width = 1.5) const;

    // Рисует вектор второй производной кривой
    void paint_second_deriv(const NURBS &points_NURBS, const QString &label, const Qt::PenStyle &penStyle = Qt::PenStyle::SolidLine, const QColor &color = QColor(0, 0, 0, 255), double width = 1.5) const;

    // Рисует аппроксимирующие кривые, на расстоянии эпсилон от исходной кривой
    void paint_curves_shift(const QPair<QVector<Point_curve>, QVector<Point_curve>> &points_NURBS) const;

private:
    QCustomPlot* _canvas;
    QString _title;
};

#endif // CHART_H
