#include "widget.h"
#include "ui_widget.h"
#include "functions.h"
#include "charts.h"
#include <QDebug>

// Проецирует точку в пространстве на кривую
void pointProjection(const int& n, const int& p, const std::vector<double>& u_vector, const QVector<QVector<double>>& polygon, const std::vector<double>& h,
                     const QPair<double, double>& point, Ui::Widget* ui)
{
    QVector<Point_curve> point_u(n - 1); // Массив точек - перпендикуляров

    std::vector<std::vector<double>> nders(p + 1, std::vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и  1-я и 2-я производные
    std::vector<QPair<double, double>> c2(p + 1); // Индекс 2 для 2D задачи

    // Реальный диапазон
    double u_start = u_vector[p];
    const double u_end = u_vector[n + 1];

    QVector<double> u_real_span; // Спаны реального диапазона узлового вектора

    for(int i = 1; u_start < u_end; ++i)
    {
        u_real_span.push_back(u_start);
        u_start = u_vector[p + i];
    }

    u_real_span.push_back(u_end);

    for(int i = 0; i < u_real_span.size() - 1; ++i)
    {
        point_u[i].u = (u_real_span[i + 1] - u_real_span[i]) / 2 + u_real_span[i]; // Берём среднее спана

        for(int k = 0; k < 15; ++k)
        {
            if(point_u[i].u < u_real_span[i]) // Если точка вышла из спана
            {
                point_u[i].u = u_real_span[i];
                curve_point_and_deriv_NURBS(point_u, n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
                point_u[i].curve = c2[0];
                point_u[i].derivative_1 = c2[1];
                point_u[i].derivative_2 = c2[2];
                break;
            }
            else if(point_u[i].u > u_real_span[i + 1])
            {
                point_u[i].u = u_real_span[i + 1];
                curve_point_and_deriv_NURBS(point_u, n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
                point_u[i].curve = c2[0];
                point_u[i].derivative_1 = c2[1];
                point_u[i].derivative_2 = c2[2];
                break;
            }

            curve_point_and_deriv_NURBS(point_u, n, p, u_vector, polygon, h, point_u[i].u, c2, nders);

            point_u[i].curve = c2[0];
            point_u[i].derivative_1 = c2[1];
            point_u[i].derivative_2 = c2[2];

            qDebug() << point_u[i].u;

            double x = point_u[i].curve.first - point.first;
            double y = point_u[i].curve.second - point.second;
            double numerator = x * point_u[i].derivative_1.first + y * point_u[i].derivative_1.second;
            double denominator = x * point_u[i].derivative_2.first + y * point_u[i].derivative_2.second + pow(vector_len(point_u[i].derivative_1), 2);
            point_u[i].u = point_u[i].u - numerator / denominator; // Новая точка кривой
        }
    }

    Point_curve min_point;
    min_point = point_u[0];
    double min_len = vector_len(point, point_u[0].curve);

    for(int i = 1; i < point_u.size(); ++i)
    {
        double temp_len = vector_len(point, point_u[i].curve);

        if(min_len > temp_len)
        {
            min_len = temp_len;
            min_point = point_u[i];
        }
    }

    QCPItemLine *line = new QCPItemLine(ui->graph_function);
    line->start->setCoords(min_point.curve.first, min_point.curve.second);
    line->end->setCoords(point.first, point.second);

    ui->graph_function->replot();
}

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
        curve_point_and_deriv_NURBS(data_NURBS, n, p, u, b, h, u_i, c2, nders);

        data_NURBS[i].u = u_i;
        data_NURBS[i].curve = c2[0];
        data_NURBS[i].derivative_1 = c2[1];
        data_NURBS[i].derivative_2 = c2[2];
    }

    QString title = "B-сплайн 4-го порядка";
    QString labels_legend_1 = "Контур. многоуг.";
    QString labels_legend_2 = "B-сплайн";

    int x_min = 0, x_max = 20;
    int y_min = -5, y_max = 15;

    curve_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    title = "1-я прoизвдная B-сплайна 4-го порядка";
    x_min = -10, x_max = 30;
    y_min = -25, y_max = 30;

    first_derivative_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    title = "2-я прoизвдная B-сплайна 4-го порядка";
    x_min = -95, x_max = 60;
    y_min = -110, y_max = 120;

    second_derivative_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    QVector<double> point_u { 0, 0.2, 0.4, 0.6, 1 }; // Массив, хранящий точки u, от которых пойдёт производная на графике
    QVector<Point_curve> derivs_curve(point_u.size());

    for(int i = 0; i < point_u.size(); ++i) // Считаем координаты и производные для точек u
    {
        curve_point_and_deriv_NURBS(derivs_curve, n, p, u, b, h, point_u[i], c2, nders);

        derivs_curve[i].u = point_u[i];
        derivs_curve[i].curve = c2[0];
        derivs_curve[i].derivative_1 = c2[1];
        derivs_curve[i].derivative_2 = c2[2];
    }

    derivative_point_line(derivs_curve, ui); // Рисует линию производной в точке (касательную)

    QPair<double, double> point(4.5, 4.5);

    pointProjection(n, p, u, b, h, point, ui);
}

Widget::~Widget()
{
    delete ui;
}
