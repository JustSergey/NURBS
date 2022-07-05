#ifndef CHARTS_H
#define CHARTS_H

#include "functions.h"
#include "ui_widget.h"

// Рисует линию (касательную) производной в точке
void derivative_point_line(const QVector<Point_curve>& data_NURBS, Ui::Widget* ui)
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

// Рисует кривую и многоугольник по заданным точкам
void curve_plot(const QVector<QVector<double>>& b, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
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

    for(const auto& el: b) // Рисуем точки
        curve_point->addData(el[0], el[1]);

    curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    QCPCurve *curve_spline = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    pen.setColor(QColor(30, 144, 255));
    curve_spline->setPen(pen);

    for(const auto& el: data_NURBS) // Рисуем сплайн
        curve_spline->addData(el.curve.first, el.curve.second);

    ui->graph_function->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    ui->graph_function->legend->setVisible(true); // Включаем легенду графика
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

// Рисует первую производную и многоугольник по заданным точкам
void first_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
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

// Рисует вторую производную и многоугольник по заданным точкам
void second_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
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

// Рисует линию между двумя точками
void plot_line(const QPair<double, double>& point, const Point_curve& u_perpendicula, Ui::Widget* ui)
{
    QCPItemLine *line = new QCPItemLine(ui->graph_function);
    line->start->setCoords(u_perpendicula.curve.first, u_perpendicula.curve.second);
    line->end->setCoords(point.first, point.second);
    ui->graph_function->replot();
}

// Рисует касательную к точке
void plot_tangent(const Point_curve& point, Ui::Widget* ui)
{
    QCPItemLine *line = new QCPItemLine(ui->graph_function);
    line->setPen(QColor(0, 128, 0)); // Задаем зелёный цвет
    line->start->setCoords(point.curve.first - point.derivative_1.first / 10, point.curve.second - point.derivative_1.second / 10);
    line->end->setCoords(point.curve.first + point.derivative_1.first / 10, point.curve.second + point.derivative_1.second / 10);
    ui->graph_function->replot();
}

// Рисует точки на графике
void plot_point(const QVector<Point_curve>& points, const QColor& color, Ui::Widget* ui)
{
    for(const auto& p: points)
    {
        ui->graph_function->addGraph();
        ui->graph_function->graph()->setPen(color); // Задаем чёрный цвет)
        ui->graph_function->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 8)); // Формируем вид точек
        ui->graph_function->graph()->setLineStyle(QCPGraph::lsNone);
        ui->graph_function->graph()->addData(p.curve.first, p.curve.second);
        ui->graph_function->replot();
    }
}

#endif // CHARTS_H
