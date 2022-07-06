#include "widget.h"
#include "ui_widget.h"
#include "functions.h"
#include "charts.h"
#include <QDebug>

Widget::Widget(QWidget *parent)
    : QWidget(parent), ui(new Ui::Widget)
{
    ui->setupUi(this);

    using namespace std;

    const QVector<QVector<double>> b // Массив точек многоугольника
    {
        {1.25, 1.3},
        {2.5, 3.9},
        {5.6, 3.9},
        {6.25, 1.3},
        {7.5, 2.6}
    };

    const vector<double> h {1, 1, 1, 1, 1};              // Весовые коэффициенты
    const vector<double> u {0, 0, 0, 0.4, 0.6, 1, 1, 1}; // Узловой вектор

    const uint p = 2;            // Степень аппроксимирующих полиномов
    const uint m = u.size() - 1; // n_kn - количество узлов (длина) в узловом векторе
    const uint n = m - p - 1;    // n_real - количество узлов (длина) реальной части узлового вектора

    // Реальный диапазон
    const double u_start = u[p];
    const double u_stop = u[n + 1];

    vector<vector<double>> nders(p + 1, vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и  1-я и 2-я производные
    vector<QPair<double, double>> c2(p + 1); // Индекс 2 для 2D задачи

    const int n_u = 60; // Кол-во разбиений (точек) в реальной части узлов. вектора

    QVector<Point_curve> data_NURBS(n_u + 1); // Содержит точки кривой, 1-ой, 2-ой производной

    for(int i = 0; i < n_u + 1; ++i)
    {
        double u_i = (i / static_cast<double>(n_u)) * (u_stop - u_start);
        curve_point_and_deriv_NURBS(data_NURBS[i], n, p, u, b, h, u_i, c2, nders);

        data_NURBS[i].u = u_i;
        data_NURBS[i].curve = c2[0];
        data_NURBS[i].derivative_1 = c2[1];
        data_NURBS[i].derivative_2 = c2[2];
    }



    QString title = "B-сплайн 3-го порядка";
    QString labels_legend_1 = "Контур. многоуг.";
    QString labels_legend_2 = "B-сплайн";

    int x_min = -10, x_max = 20;
    int y_min = -10, y_max = 20;

    curve_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    title = "1-я прoизвдная B-сплайна 3-го порядка";
    x_min = -10, x_max = 30;
    y_min = -25, y_max = 30;

    first_derivative_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    title = "2-я прoизвдная B-сплайна 3-го порядка";
    x_min = -95, x_max = 60;
    y_min = -110, y_max = 120;

    second_derivative_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);




    QVector<double> point_u { 0, 0.2, 0.4, 0.6, 1 }; // Массив, хранящий точки u, от которых пойдёт производная на графике
    QVector<Point_curve> derivs_curve(point_u.size());

    for(int i = 0; i < point_u.size(); ++i) // Считаем координаты и производные для точек u
    {
        curve_point_and_deriv_NURBS(derivs_curve[i], n, p, u, b, h, point_u[i], c2, nders);

        derivs_curve[i].u = point_u[i];
        derivs_curve[i].curve = c2[0];
        derivs_curve[i].derivative_1 = c2[1];
        derivs_curve[i].derivative_2 = c2[2];
    }

    derivative_point_line(derivs_curve, ui); // Рисуем линию производной в точке (касательную)



    QPair<double, double> point(5.5, 4.5); // Точка, к которой мы будем проводить перпендикуляр

    Point_curve u_perpendicular = finding_perpendicular(n, p, u, b, h, point); // Ближайшая точка для перпендикуляра

    plot_tangent(u_perpendicular, ui); // Рисуем касательную к точке
    plot_line(point.first, point.second, u_perpendicular.curve.first, u_perpendicular.curve.second, ui); // Рисуем перпендикуляр между точкой и кривой



    QVector<double> real_spans = real_span_calc(p, n, u); // Вектор с точками реального диапазона спанов
    QVector<Point_curve> u_real_spans(real_spans.size());

    for(int i = 0; i < u_real_spans.size(); ++i)
    {
        curve_point_and_deriv_NURBS(u_real_spans[i], n, p, u, b, h, real_spans[i], c2, nders);
        u_real_spans[i].curve = c2[0];
        u_real_spans[i].derivative_1 = c2[1];
        u_real_spans[i].derivative_2 = c2[2];
        plot_point(u_real_spans[i].curve.first, u_real_spans[i].curve.second, ui, "", QColor(123, 104, 238)); // Рисуем точки на графике (точки границ реального диапазона спанов)
    }

}

Widget::~Widget()
{
    delete ui;
}
