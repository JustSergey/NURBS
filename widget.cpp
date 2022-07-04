#include "widget.h"
#include "ui_widget.h"
#include "functions.h"
#include "charts.h"
#include <QDebug>

void point_distance(const int& n, const int& p, const std::vector<double>& u_vector, const QVector<QVector<double>>& b,
                    const std::vector<double>& h, const QPair<double, double>& point, Ui::Widget* ui)
{
    QVector<Point_curve> point_u(1); // Точка кривой
    point_u[0].u = 0.2; // Тестовое

    std::vector<std::vector<double>> nders(p + 1, std::vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и  1-я и 2-я производные
    std::vector<QPair<double, double>> c2(p + 1); // Индекс 2 для 2D задачи

    for(int i = 0; i < 15; ++i)
    {
        curve_point_and_deriv_NURBS(point_u, n, p, u_vector, b, h, point_u[0].u, c2, nders);

        point_u[0].curve = c2[0];
        point_u[0].derivative_1 = c2[1];
        point_u[0].derivative_2 = c2[2];

        qDebug() << point_u[0].u;

        for(const auto& p: point_u)
        {
            double x = p.curve.first - point.first;
            double y = p.curve.second - point.second;
            double numerator = x * p.derivative_1.first + y * p.derivative_1.second;
            double denominator = x * p.derivative_2.first + y * p.derivative_2.second + pow(vector_len(p.derivative_1), 2);
            point_u[0].u = p.u - numerator / denominator;
        }
    }

    QCPItemLine *line = new QCPItemLine(ui->graph_function);
    line->start->setCoords(point_u[0].curve.first, point_u[0].curve.second);
    line->end->setCoords(point.first, point.second);

    ui->graph_function->replot();
}




void function_plot(const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS,
                   const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                   const QString& title, const QString& labels_legend_1, const QString& labels_legend_2, Ui::Widget* ui)
{
    ui->graph_function->clearGraphs(); // Очищаем все графики

    QCPCurve *curve_point = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);

    curve_point->setPen(QColor(0, 0, 0, 255)); // Задаем чёрный цвет
    curve_point->setLineStyle(QCPCurve::lsNone); // Убираем линии
    curve_point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    QPen pen;
    pen.setWidth(2); // Устанавливаем ширину
    curve_point->setPen(pen);

    for(const auto& el: data_point) // Рисуем точки
        curve_point->addData(el[0], el[1]);

    curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    QCPCurve *curve_spline = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    pen.setColor(QColor(30, 144, 255));
    curve_spline->setPen(pen);

    for(const auto& el: data_NURBS) // Рисуем сплайн
        curve_spline->addData(el.curve.first, el.curve.second);

    ui->graph_function->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    ui->graph_function->legend->setVisible(true); // Включаем Легенду графика
    curve_point->setName(labels_legend_1);
    curve_spline->setName(labels_legend_2);

    // Подписываем оси Ox и Oy
    ui->graph_function->xAxis->setLabel("Ось X");
    ui->graph_function->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    ui->graph_function->xAxis->setRange(x_min, x_max); // Для оси Ox
    ui->graph_function->yAxis->setRange(y_min, y_max); // Для оси Oy

    ui->graph_function->plotLayout()->insertRow(0); // Вставляем строку
    ui->graph_function->plotLayout()->addElement (0, 0, new QCPTextElement(ui->graph_function, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец заглавие

    ui->graph_function->replot();
}

void first_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS,
                           const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                           const QString& title, const QString& labels_legend_1, const QString& labels_legend_2, Ui::Widget* ui)
{
    ui->graph_first_derivative->clearGraphs(); // Очищаем все графики

    QCPCurve *curve_point = new QCPCurve(ui->graph_first_derivative->xAxis, ui->graph_first_derivative->yAxis);

    curve_point->setPen(QColor(0, 0, 0, 255)); // Задаем чёрный цвет
    curve_point->setLineStyle(QCPCurve::lsNone); // Убираем линии
    curve_point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    for(const auto& el: data_point) // Рисуем точки
        curve_point->addData(el[0], el[1]);

    curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    QCPCurve *curve_deriv { nullptr };
    curve_deriv = new QCPCurve(ui->graph_first_derivative->xAxis, ui->graph_first_derivative->yAxis);

    for(const auto& el: data_NURBS) // Рисуем сплайн
        curve_deriv->addData(el.derivative_1.first, el.derivative_1.second);

    ui->graph_first_derivative->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    ui->graph_first_derivative->legend->setVisible(true); // Включаем Легенду графика
    curve_point->setName(labels_legend_1);
    curve_deriv->setName(labels_legend_2);

    // Подписываем оси Ox и Oy
    ui->graph_first_derivative->xAxis->setLabel("Ось X");
    ui->graph_first_derivative->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    ui->graph_first_derivative->xAxis->setRange(x_min, x_max); // Для оси Ox
    ui->graph_first_derivative->yAxis->setRange(y_min, y_max); // Для оси Oy

    ui->graph_first_derivative->plotLayout()->insertRow(0); // Вставляем строку
    ui->graph_first_derivative->plotLayout()->addElement (0, 0, new QCPTextElement(ui->graph_first_derivative, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец заглавие

    ui->graph_first_derivative->replot();
}

void second_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS,
                           const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                           const QString& title, const QString& labels_legend_1, const QString& labels_legend_2, Ui::Widget* ui)
{
    ui->graph_second_derivative->clearGraphs(); // Очищаем все графики

    QCPCurve *curve_point = new QCPCurve(ui->graph_second_derivative->xAxis, ui->graph_second_derivative->yAxis);

    curve_point->setPen(QColor(0, 0, 0, 255)); // Задаем чёрный цвет
    curve_point->setLineStyle(QCPCurve::lsNone); // Убираем линии
    curve_point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    for(const auto& el: data_point) // Рисуем точки
        curve_point->addData(el[0], el[1]);

    curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    QCPCurve *curve_deriv { nullptr };
    curve_deriv = new QCPCurve(ui->graph_second_derivative->xAxis, ui->graph_second_derivative->yAxis);

    for(const auto& el: data_NURBS) // Рисуем сплайн
        curve_deriv->addData(el.derivative_2.first, el.derivative_2.second);

    ui->graph_second_derivative->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    ui->graph_second_derivative->legend->setVisible(true); // Включаем Легенду графика
    curve_point->setName(labels_legend_1);
    curve_deriv->setName(labels_legend_2);

    // Подписываем оси Ox и Oy
    ui->graph_second_derivative->xAxis->setLabel("Ось X");
    ui->graph_second_derivative->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    ui->graph_second_derivative->xAxis->setRange(x_min, x_max); // Для оси Ox
    ui->graph_second_derivative->yAxis->setRange(y_min, y_max); // Для оси Oy

    ui->graph_second_derivative->plotLayout()->insertRow(0); // Вставляем строку
    ui->graph_second_derivative->plotLayout()->addElement (0, 0, new QCPTextElement(ui->graph_second_derivative, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец заглавие

    ui->graph_second_derivative->replot();
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

    function_plot(b, data_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

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

    plot_deriv_for_point(derivs_curve, ui);

    QPair<double, double> point(1, 4);

    point_distance(n, p, u, b, h, point, ui);
}


Widget::~Widget()
{
    delete ui;
}
