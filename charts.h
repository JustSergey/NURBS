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

// Рисует многоугольник с вершинами
void plot_polygon(const QVector<QVector<double>>& polygon, Ui::Widget* ui, const QString& label, const QColor& color = QColor(0, 0, 0, 255), const double& width = 1)
{
    QCPCurve *shape = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);

    shape->setPen(color);
    shape->setLineStyle(QCPCurve::lsNone); // Убираем линии
    shape->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    QPen pen;
    pen.setWidth(width); // Устанавливаем ширину
    shape->setPen(pen);

    for(const auto& point: polygon) // Рисуем точки
        shape->addData(point[0], point[1]);

    shape->setLineStyle(QCPCurve::lsLine); // Добавляем линии
    shape->setName(label); // Обзываем полигон в легенде графика

    ui->graph_second_derivative->replot();
}

// Рисует точки на графике
void plot_point(const double& x, const double& y, Ui::Widget* ui, const QString& text = "", const QColor& color = QColor(0, 0, 0, 255))
{
    ui->graph_function->addGraph();
    ui->graph_function->graph()->setPen(color);
    ui->graph_function->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 8)); // Формируем вид точек
    ui->graph_function->graph()->setLineStyle(QCPGraph::lsNone);
    ui->graph_function->graph()->addData(x, y);

    QCPItemText *textLabel = new QCPItemText(ui->graph_function);
    textLabel->position->setCoords(x + 0.1, y - 0.1);
    textLabel->setText(text);

    ui->graph_function->replot();
}

// Рисует линию между двумя точками
void plot_line(const double& point_x_1, const double& point_y_1, const double& point_x_2, const double& point_y_2, Ui::Widget* ui)
{
    QCPItemLine *line = new QCPItemLine(ui->graph_function);
    line->start->setCoords(point_x_1, point_y_1);
    line->end->setCoords(point_x_2, point_y_2);

    ui->graph_function->replot();
}

void plot_curve(const QVector<Point_curve>& data_NURBS, Ui::Widget* ui,  const QString& label, const QColor& color = QColor(0, 0, 0, 255))
{
    QCPCurve *curve = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    curve->setPen(color);

    for(const auto& point: data_NURBS) // Рисуем сплайн
        curve->addData(point.curve.first, point.curve.second);

    curve->setName(label); // Обзываем кривую в легенде графика

    ui->graph_function->replot();
}

// Рисует кривую и многоугольник по заданным точкам
void curve_plot(const QVector<QVector<double>>& b, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                const QString& title, const QString& labels_legend_1, const QString& labels_legend_2, Ui::Widget* ui)
{
    ui->graph_function->clearGraphs(); // Очищаем все графики
    ui->graph_function->legend->setVisible(true); // Включаем легенду графика

    plot_polygon(b, ui, labels_legend_1); // Рисуем многоугольник с вершинами

    plot_curve(data_NURBS, ui, labels_legend_2, QColor(30, 144, 255)); // Рисуем сплайн

    ui->graph_function->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Делаем график перетаскиваемым + масштабирование колеса прокрутки

    // Подписываем оси Ox и Oy
    ui->graph_function->xAxis->setLabel("Ось X");
    ui->graph_function->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    ui->graph_function->xAxis->setRange(x_min, x_max);
    ui->graph_function->yAxis->setRange(y_min, y_max);

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

// Рисует касательную к точке
void plot_tangent(const Point_curve& point, Ui::Widget* ui)
{
    QCPItemLine *line = new QCPItemLine(ui->graph_function);
    line->setPen(QColor(0, 128, 0)); // Задаем зелёный цвет
    line->start->setCoords(point.curve.first - point.derivative_1.first / 10, point.curve.second - point.derivative_1.second / 10);
    line->end->setCoords(point.curve.first + point.derivative_1.first / 10, point.curve.second + point.derivative_1.second / 10);
    ui->graph_function->replot();
}

#endif // CHARTS_H
