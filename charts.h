#ifndef CHARTS_H
#define CHARTS_H

#include "functions.h"
#include "ui_widget.h"

// Рисует линию (касательную) производной в точке
void plot_tangent_point(QCustomPlot* canvas, const Point_curve& point)
{
    QCPItemLine *line = new QCPItemLine(canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->start->setCoords(point.curve.first, point.curve.second);
    line->end->setCoords(point.derivative_1.first + point.curve.first, point.derivative_1.second + point.curve.second);
    canvas->replot();
}

// Рисует многоугольник с вершинами
void plot_polygon(QCustomPlot* canvas, const QVector<QVector<double>>& polygon_points, const QString& label, const QColor& color = QColor(0, 0, 0, 255), const double& width = 1)
{
    QCPCurve *shape = new QCPCurve(canvas->xAxis, canvas->yAxis);

    shape->setPen(color);
    shape->setLineStyle(QCPCurve::lsNone); // Убираем линии
    shape->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    QPen pen;
    pen.setWidth(width); // Устанавливаем ширину
    shape->setPen(pen);

    // uint counter = 0;

    for(const auto& point: polygon_points) // Рисуем точки
    {
        /*
        // Делаем подписи к каждой вершине многоугольника
        QCPItemText *label = new QCPItemText(canvas);
        label->position->setCoords(point[0] + 0.35, point[1] - 0.2);
        label->setFont(QFont("sans", 10));
        label->setText(QString("P%1").arg(counter++));
        */

        shape->addData(point[0], point[1]);
    }

    //shape->setLineStyle(QCPCurve::lsLine); // Добавляем линии
    shape->setName(label); // Обзываем полигон в легенде графика
    canvas->replot();
}

// Рисует точку на графике
void plot_point(QCustomPlot* canvas, const double& x, const double& y, const double& width = 5, const QString& text = "", const QColor& color = QColor(0, 0, 0, 255))
{
    canvas->addGraph();
    canvas->graph()->setPen(color);
    canvas->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, width)); // Формируем вид точек
    canvas->graph()->setLineStyle(QCPGraph::lsNone);
    canvas->graph()->addData(x, y);
    canvas->legend->removeItem(canvas->legend->itemCount() - 1); // Удаляем точку из легенды

    if(!text.isEmpty()) // Если есть текст для подписи к точке
    {
        QCPItemText *label = new QCPItemText(canvas); // Подпись к точке
        label->position->setCoords(x + 0.2, y - 0.2);
        label->setText(text);
    }

    canvas->replot();
}

// Рисует линию между двумя точками
void plot_line(QCustomPlot* canvas, const double& x_1, const double& y_1, const double& x_2, const double& y_2, const QColor& color = QColor(0, 0, 0), const double& width = 3.5)
{
    QCPItemLine *line = new QCPItemLine(canvas);
    QPen pen;
    pen.setStyle(Qt::PenStyle::DashLine);
    pen.setColor(color);
    pen.setWidth(width);
    line->setPen(pen);
    line->start->setCoords(x_1, y_1);
    line->end->setCoords(x_2, y_2);
    canvas->replot();
}

// Рисует касательную к точке
void plot_tangent(QCustomPlot* canvas, const Point_curve& point, const QColor& color = QColor(0, 100, 0), const double& width = 2.8)
{
    QCPItemLine *line = new QCPItemLine(canvas);
    QPen pen;
    pen.setColor(color);
    pen.setWidth(width);
    line->setPen(pen);
    line->start->setCoords(point.curve.first - point.derivative_1.first / 10, point.curve.second - point.derivative_1.second / 10);
    line->end->setCoords(point.curve.first + point.derivative_1.first / 10, point.curve.second + point.derivative_1.second / 10);
    canvas->replot();
}

// Рисует кривую NURBS
void plot_curve(QCustomPlot* canvas, const QVector<Point_curve>& data_NURBS, const QString& label, const Qt::PenStyle& penStyle = Qt::PenStyle::SolidLine,
                const QColor& color = QColor(0, 0, 0, 255), const double& width = 1.5)
{
    QCPCurve *curve = new QCPCurve(canvas->xAxis, canvas->yAxis);
    QPen pen;
    pen.setColor(QColor(color));
    pen.setStyle(penStyle);
    pen.setWidthF(width);
    curve->setPen(pen);

    for(const auto& point: data_NURBS) // Рисуем сплайн
        curve->addData(point.curve.first, point.curve.second);

    curve->setName(label); // Обзываем кривую в легенде графика
    canvas->replot();
}

// Рисует надпись
void plot_lable(QCustomPlot* canvas, const double& x, const double& y, const QString& text, const double& font_size = 10)
{
    QCPItemText *label = new QCPItemText(canvas);
    label->setFont(QFont("sans", font_size));
    label->position->setCoords(x, y);
    label->setText(text);
    canvas->replot();
}

// Рисует надпись со стрелкой
void plot_lable_with_arrow(QCustomPlot* canvas, const double& x_1, const double& y_1, const double& x_2, const double& y_2, const QString& text)
{
    QCPItemText *label = new QCPItemText(canvas);
    label->setFont(QFont("sans", 10));
    label->position->setCoords(x_1 + 0.3, y_1 + 0.7);
    label->setText(text);

    QCPItemLine *line = new QCPItemLine(canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->start->setCoords(x_1, y_1);
    line->end->setCoords(x_2, y_2);
    canvas->replot();
}

// Рисует кривую и многоугольник по заданным точкам
void curve_plot(QCustomPlot* canvas, const QVector<QVector<double>>& control_points, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                const QString& title, const QString& labels_legend_1, const QString& labels_legend_2)
{
    canvas->clearGraphs(); // Очищаем все графики
    canvas->legend->setVisible(true); // Включаем легенду графика

    plot_polygon(canvas, control_points, labels_legend_1); // Рисуем многоугольник с вершинами
    plot_curve(canvas, data_NURBS, labels_legend_2, Qt::PenStyle::SolidLine, QColor(30, 144, 255)); // Рисуем сплайн

/*
    // Рисуем подписи к спанам реального диапазон (убрать при необходимости)
    plot_lable_with_arrow(canvas, 1.5, 4.8, 2.58, 3.14, "u∈[0, 1/5)");

    QCPItemText *label = new QCPItemText(canvas);
    label->setFont(QFont("sans", 10));
    label->position->setCoords(4, 1.2);
    label->setText("u∈[1/5, 3/5)");
    QCPItemLine *line = new QCPItemLine(canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->start->setCoords(4, 1.2 + 0.2);
    line->end->setCoords(5.46, 3.64);
    canvas->replot();

    plot_lable_with_arrow(canvas, 7.4, 4.3, 6.28, 2.14, "u∈[3/5, 1)");
*/

   canvas->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Делаем график перетаскиваемым + масштабирование колеса прокрутки

    // Подписываем оси Ox и Oy
    canvas->xAxis->setLabel("Ось X");
    canvas->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    canvas->xAxis->setRange(x_min, x_max);
    canvas->yAxis->setRange(y_min, y_max);

    canvas->plotLayout()->insertRow(0); // Вставляем строку
    canvas->plotLayout()->addElement (0, 0, new QCPTextElement(canvas, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец оглавление

    canvas->replot();
}

// Рисует первую производную и многоугольник по заданным точкам
void first_derivative_plot(QCustomPlot* canvas, const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                           const QString& title, const QString& labels_legend_1, const QString& labels_legend_2)
{
    canvas->clearGraphs(); // Очищаем все графики

    //QCPCurve *curve_point = new QCPCurve(canvas->xAxis, canvas->yAxis);

    //curve_point->setPen(QColor(0, 0, 0, 255)); // Задаем чёрный цвет
    //curve_point->setLineStyle(QCPCurve::lsNone); // Убираем линии
    //curve_point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    //for(const auto& el: data_point) // Рисуем точки
    //    curve_point->addData(el[0], el[1]);

    //curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    QCPCurve *curve_deriv = new QCPCurve(canvas->xAxis, canvas->yAxis);

    for(const auto& el: data_NURBS) // Рисуем сплайн
        curve_deriv->addData(el.derivative_1.first, el.derivative_1.second);

    canvas->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    canvas->legend->setVisible(true); // Включаем Легенду графика
    //curve_point->setName(labels_legend_1);
    curve_deriv->setName(labels_legend_2);

    // Подписываем оси Ox и Oy
    canvas->xAxis->setLabel("Ось X");
    canvas->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    canvas->xAxis->setRange(x_min, x_max); // Для оси Ox
    canvas->yAxis->setRange(y_min, y_max); // Для оси Oy

    canvas->plotLayout()->insertRow(0); // Вставляем строку
    canvas->plotLayout()->addElement (0, 0, new QCPTextElement(canvas, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец заглавие

    canvas->replot();
}

// Рисует вторую производную и многоугольник по заданным точкам
void second_derivative_plot(QCustomPlot* canvas, const QVector<QVector<double>>& data_point, const QVector<Point_curve>& data_NURBS, const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                            const QString& title, const QString& labels_legend_1, const QString& labels_legend_2)
{
    canvas->clearGraphs(); // Очищаем все графики

    QCPCurve *curve_point = new QCPCurve(canvas->xAxis, canvas->yAxis);

    curve_point->setPen(QColor(0, 0, 0, 255)); // Задаем чёрный цвет
    curve_point->setLineStyle(QCPCurve::lsNone); // Убираем линии
    curve_point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    for(const auto& el: data_point) // Рисуем точки
        curve_point->addData(el[0], el[1]);

    curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    QCPCurve *curve_deriv = new QCPCurve(canvas->xAxis, canvas->yAxis);

    for(const auto& el: data_NURBS) // Рисуем сплайн
        curve_deriv->addData(el.derivative_2.first, el.derivative_2.second);

    canvas->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    canvas->legend->setVisible(true); // Включаем Легенду графика
    curve_point->setName(labels_legend_1);
    curve_deriv->setName(labels_legend_2);

    // Подписываем оси Ox и Oy
    canvas->xAxis->setLabel("Ось X");
    canvas->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    canvas->xAxis->setRange(x_min, x_max); // Для оси Ox
    canvas->yAxis->setRange(y_min, y_max); // Для оси Oy

    canvas->plotLayout()->insertRow(0); // Вставляем строку
    canvas->plotLayout()->addElement (0, 0, new QCPTextElement(canvas, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец заглавие

    canvas->replot();
}

#endif // CHARTS_H
