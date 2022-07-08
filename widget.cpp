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
    const int n_u = 60;          // Кол-во разбиений (точек) в реальной части узлов. вектора

    // Реальный диапазон
    const double u_start = u[p];
    const double u_stop = u[n + 1];

    vector<vector<double>> nders(p + 1, vector<double>(p + 1)); // Содержит для заданного "u" массив BASIS функций и 1-ую и 2-ую производную
    vector<QPair<double, double>> c2(p + 1);  // Индекс 2 для 2D задачи
    QVector<Point_curve> data_NURBS(n_u + 1); // Содержит точки кривой, 1-ой и 2-ой производной

    for(int i = 0; i < n_u + 1; ++i)
    {
        double u_i = (i / static_cast<double>(n_u)) * (u_stop - u_start);
        curve_point_and_deriv_NURBS(data_NURBS[i], n, p, u, b, h, u_i, c2, nders);
    }



    QString title = "Кратчайшее расстояние от точки до интервалов";
    QString labels_legend_1 = "Определяющий многоуг.";
    QString labels_legend_2 = "NURBS";

    int x_min = -10, x_max = 10;
    int y_min = -10, y_max = 10;

    curve_plot(ui->graph_function, b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);

    title = "1-я прoизвдная NURBS 2-ой степени";
    x_min = -15, x_max = 15;
    y_min = -15, y_max = 15;

    first_derivative_plot(ui->graph_first_derivative, b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);

    title = "2-я прoизвдная NURBS 2-ой степени";
    x_min = -100, x_max = 100;
    y_min = -100, y_max = 100;

    second_derivative_plot(ui->graph_second_derivative, b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);




    QVector<double> real_spans = real_span_calc(p, n, u); // Вектор с точками реального диапазона спанов
    QVector<Point_curve> u_real_spans(real_spans.size());

    for(int i = 0; i < u_real_spans.size(); ++i)
    {
        curve_point_and_deriv_NURBS(u_real_spans[i], n, p, u, b, h, real_spans[i], c2, nders);
        plot_point(ui->graph_function, u_real_spans[i].curve.first, u_real_spans[i].curve.second, 9, "", QColor(147, 112, 219)); // Рисуем точки на графике (точки границ реального диапазона спанов)
    }

    // *ДОБАВЛЯЕТ В ЛЕГЕНДУ* //
    ui->graph_function->addGraph();
    ui->graph_function->graph()->setPen(QColor(147, 112, 219));
    ui->graph_function->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 9)); // Формируем вид точек
    ui->graph_function->graph()->setLineStyle(QCPGraph::lsNone);
    ui->graph_function->graph()->addData(u_real_spans[0].curve.first, u_real_spans[0].curve.second);
    ui->graph_function->graph()->setName("Точка смены"); // Точка смены полинома/элемента узлового вектора

    ui->graph_function->replot();


    /*
    QVector<double> point_u { 0, 0.2, 0.4, 0.6, 0.8, 1 }; // Массив, хранящий точки u, от которых пойдёт производная на графике
    QVector<Point_curve> derivs_curve(point_u.size());

    for(int i = 0; i < point_u.size(); ++i) // Считаем координаты и производные для точек u
    {
        curve_point_and_deriv_NURBS(derivs_curve[i], n, p, u, b, h, point_u[i], c2, nders);

        derivs_curve[i].u = point_u[i];
        derivs_curve[i].curve = c2[0];
        derivs_curve[i].derivative_1 = c2[1];
        derivs_curve[i].derivative_2 = c2[2];
    }

    derivative_point_line(ui->graph_function, derivs_curve); // Рисуем линию производной в точке (касательную)

    plot_lable(ui->graph_function, derivs_curve[0].derivative_1.first + derivs_curve[0].curve.first + 1, derivs_curve[0].derivative_1.second + derivs_curve[0].curve.second + 0.2, "C'(0)");
    plot_lable(ui->graph_function,derivs_curve[1].derivative_1.first + derivs_curve[1].curve.first + 1.4, derivs_curve[1].derivative_1.second + derivs_curve[1].curve.second + 0.4, "C'(1/5)");
    plot_lable(ui->graph_function,derivs_curve[2].derivative_1.first + derivs_curve[2].curve.first + 1, derivs_curve[2].derivative_1.second + derivs_curve[2].curve.second - 0.5, "C'(2/5)");
    plot_lable(ui->graph_function,derivs_curve[3].derivative_1.first + derivs_curve[3].curve.first + 1, derivs_curve[3].derivative_1.second + derivs_curve[3].curve.second - 0.5, "C'(3/5)");
    plot_lable(ui->graph_function,derivs_curve[4].derivative_1.first + derivs_curve[4].curve.first + 1.3, derivs_curve[4].derivative_1.second + derivs_curve[4].curve.second - 0.5, "C'(4/5)");
    plot_lable(ui->graph_function,derivs_curve[5].derivative_1.first + derivs_curve[5].curve.first + 1, derivs_curve[5].derivative_1.second + derivs_curve[5].curve.second + 0.1, "C'(1)");
    derivative_point_line(ui->graph_first_derivative, derivs_curve); // Рисуем линию производной в точке (касательную)

    plot_lable(ui->graph_first_derivative, derivs_curve[0].derivative_1.first + 1, derivs_curve[0].derivative_1.second + 0.5, "C'(0)");
    plot_lable(ui->graph_first_derivative, derivs_curve[1].derivative_1.first + 1.5, derivs_curve[1].derivative_1.second + 0.3, "C'(1/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[2].derivative_1.first + 1.4, derivs_curve[2].derivative_1.second - 0.5, "C'(2/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[3].derivative_1.first + 0.5, derivs_curve[3].derivative_1.second - 0.8, "C'(3/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[4].derivative_1.first + 1.4, derivs_curve[4].derivative_1.second - 0.5, "C'(4/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[5].derivative_1.first + 0.5, derivs_curve[5].derivative_1.second + 0.8, "C'(1)");

    for(const auto& p: derivs_curve)
    {
        plot_point(ui->graph_function, p.curve.first, p.curve.second);
    }
*/


    /*
    QPair<double, double> point(2, 6); // Точка, к которой мы будем проводить перпендикуляр
    Point_curve u_perpendicular = finding_perpendicular(n, p, u, b, h, point, ui->graph_function); // Ближайшая точка для перпендикуляра

    //plot_line(ui->graph_function, point.first, point.second, u_perpendicular.curve.first, u_perpendicular.curve.second, QColor(178, 34, 34)); // Рисуем перпендикуляр между точкой и кривой
    //plot_tangent(ui->graph_function, u_perpendicular); // Рисуем касательную к точке

    // *ДОБАВЛЯЕТ В ЛЕГЕНДУ* //
    ui->graph_function->addGraph();
    ui->graph_function->graph()->setPen(QColor(244, 164, 96));
    ui->graph_function->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 11)); // Формируем вид точек
    ui->graph_function->graph()->setLineStyle(QCPGraph::lsNone);
    ui->graph_function->graph()->addData(point.first, point.second);
    ui->graph_function->graph()->setName("Заданная точка"); // Точка смены полинома/элемента узлового вектора



    QCPCurve *curve = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    QPen pen;
    pen.setColor(QColor(0, 100, 0));
    pen.setWidthF(2.8);
    curve->setPen(pen);
    curve->setName("Касательная"); // Обзываем полигон в легенде графика
*/
    /*
    QCPCurve *curve1 = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    QPen pen1;
    pen1.setWidthF(2);
    pen1.setColor(QColor(0, 0, 0));
    pen1.setStyle(Qt::PenStyle::DashLine);
    curve1->setPen(pen1);
    curve1->setName("Перпендикуляр"); // Обзываем полигон в легенде графика

    plot_lable_with_arrow(ui->graph_function, 6.35, 6.35, 4.8, 4.27, "Предельные\nслучаи");
    plot_lable_with_arrow(ui->graph_function, 6.35, 6.35, 4, 4.36, "");
    plot_lable_with_arrow(ui->graph_function, 2, 3.5, 2.75, 4.5, "");
    plot_lable(ui->graph_function, 0.8, 3.1, "Перпендикуляр\nс мин. длиной");
    */
}

Widget::~Widget()
{
    delete ui;
}
