#include "widget.h"
#include "ui_widget.h"
#include "nurbs.h"
#include "chart.h"

Widget::Widget(QWidget *parent)
    : QWidget(parent), ui(new Ui::Widget)
{
    ui->setupUi(this);

    const QVector<QVector<double>> control_points // Точки определяющего многоугольника
    {
        {1, 1},
        {3.5, 1.25},
        {5, 4.5},
        {7.5, 2},
        {10, 4},
        {12, 1.5}
    };

    const QVector<double> weight {1, 1, 1, 1, 1, 1}; // Весовые коэффициенты
    const uint degree = 5;      // Степень аппроксимирующих полиномов
    const int num_splits = 60;  // Кол-во разбиений (точек) в реальной части узлов. вектора

    NURBS curve(control_points, weight, degree, num_splits);
    curve.fill_nodal_vector();
    curve.calc_points_derivs_NURBS();

    Chart canvas_1(ui->canvas_1, "Название");
    canvas_1.paint_polygon(control_points, "Многоугольник");
    canvas_1.paint_curve(curve, "Исходный NURBS", Qt::PenStyle::SolidLine, QColor(30, 144, 255));
    canvas_1.paint_curves_shift(curve.calc_curves_shift(1));
}

Widget::~Widget()
{
    delete ui;
}

