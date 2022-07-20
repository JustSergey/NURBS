#include "chart.h"

Chart::Chart(QCustomPlot *canvas, const QString &title)
    : m_canvas { canvas }, m_title { title }
{
    m_canvas->clearGraphs(); // Очищаем все графики
    m_canvas->legend->setVisible(true); // Включаем легенду графика
    m_canvas->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Делаем график перетаскиваемым + масштабирование колеса прокрутки
    m_canvas->plotLayout()->insertRow(0); // Вставляем строку
    m_canvas->plotLayout()->addElement (0, 0, new QCPTextElement(m_canvas, m_title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец оглавление
    m_canvas->xAxis->setRange(-15, 15);
    m_canvas->yAxis->setRange(-15, 15);
}

void Chart::set_title(const QString &title)
{
    m_title = title;
}

// Рисует линию (касательную), начинающуюся от точки кривой (point)
void Chart::paint_tangent(const Point_curve& point_curve) const
{
    QCPItemLine *line = new QCPItemLine(m_canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->start->setCoords(point_curve.curve.first, point_curve.curve.second);
    line->end->setCoords(point_curve.deriv_1.first + point_curve.curve.first, point_curve.deriv_1.second + point_curve.curve.second);
    m_canvas->replot();
}

// Рисует определяющий многоугольник с вершинами
void Chart::paint_polygon(const QVector<QVector<double>> &polygon_points, const QString &label, const QColor &color, double width) const
{
    QCPCurve *shape = new QCPCurve(m_canvas->xAxis, m_canvas->yAxis);
    shape->setPen(color);
    shape->setLineStyle(QCPCurve::lsNone); // Убираем линии
    shape->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    QPen pen;
    pen.setWidth(width);
    shape->setPen(pen);

    // uint counter = 0;

    for(const auto& point: polygon_points) // Рисуем точки
    {
        shape->addData(point[0], point[1]);

        /* Делает подписи к каждой вершине многоугольника
        QCPItemText *label = new QCPItemText(canvas);
        label->position->setCoords(point[0] + 0.35, point[1] - 0.2);
        label->setFont(QFont("sans", 10));
        label->setText(QString("P%1").arg(counter++));
        */
    }

    shape->setLineStyle(QCPCurve::lsLine); // Добавляем линии
    shape->setName(label); // Устанавливаем название полигона в легенде графика
    m_canvas->replot();
}

// Рисует точку
void Chart::paint_point(double x, double y, double width, const QString &text, const QColor &color) const
{
    m_canvas->addGraph();
    m_canvas->graph()->setPen(color);
    m_canvas->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, width)); // Формируем вид точек
    m_canvas->graph()->setLineStyle(QCPGraph::lsNone);
    m_canvas->graph()->addData(x, y);
    m_canvas->legend->removeItem(m_canvas->legend->itemCount() - 1); // Удаляем точку из легенды, если надо

    if(!text.isEmpty()) // Если есть текст для подписи к точке
    {
        QCPItemText *label = new QCPItemText(m_canvas); // Подпись к точке
        label->position->setCoords(x + 0.2, y - 0.2);
        label->setText(text);
    }

    m_canvas->replot();
}

// Рисует линию по двум точкам
void Chart::paint_line(double x_1, double y_1, double x_2, double y_2, const QColor &color, double width) const
{
    QCPItemLine *line = new QCPItemLine(m_canvas);
    QPen pen;
    pen.setStyle(Qt::PenStyle::DashLine);
    pen.setColor(color);
    pen.setWidth(width);
    line->setPen(pen);
    line->start->setCoords(x_1, y_1);
    line->end->setCoords(x_2, y_2);
    m_canvas->replot();
}

// Рисует двойную стрелку "<-->" по двум точкам
void Chart::paint_double_arrow(double x_1, double y_1, double x_2, double y_2) const
{
    QCPItemLine *line = new QCPItemLine(m_canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->setTail(QCPLineEnding::esFlatArrow);
    line->start->setCoords(x_1, y_1);
    line->end->setCoords(x_2, y_2);
    m_canvas->replot();
}

// Рисует стрелку "-->" по двум точкам
void Chart::paint_arrow(double x_1, double y_1, double x_2, double y_2) const
{
    QCPItemLine *line = new QCPItemLine(m_canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->start->setCoords(x_1, y_1);
    line->end->setCoords(x_2, y_2);
    m_canvas->replot();
}

// Рисует линию (касательную) c центром в точке кривой (point)
void Chart::paint_tangent_centred(const Point_curve &point, const QColor &color, double width) const
{
    QCPItemLine *line = new QCPItemLine(m_canvas);
    QPen pen;
    pen.setColor(color);
    pen.setWidth(width);
    line->setPen(pen);
    line->start->setCoords(point.curve.first - point.deriv_1.first / 10, point.curve.second - point.deriv_1.second / 10);
    line->end->setCoords(point.curve.first + point.deriv_1.first / 10, point.curve.second + point.deriv_1.second / 10);
    m_canvas->replot();
}

// Рисует кривую NURBS
void Chart::paint_curve(const NURBS &points_NURBS, const QString &label, const Qt::PenStyle &penStyle, const QColor &color, double width) const
{
    QCPCurve *curve = new QCPCurve(m_canvas->xAxis, m_canvas->yAxis);
    QPen pen;
    pen.setColor(QColor(color));
    pen.setStyle(penStyle);
    pen.setWidthF(width);
    curve->setPen(pen);


    for(const auto& point: points_NURBS.get_point_curve())
        curve->addData(point.curve.first, point.curve.second);

    curve->setName(label);
    m_canvas->replot();
}

// Рисует кривую NURBS
void Chart::paint_curve(const QVector<Point_curve> &points_NURBS, const QString &label, const Qt::PenStyle &penStyle, const QColor &color, double width) const
{
    QCPCurve *curve = new QCPCurve(m_canvas->xAxis, m_canvas->yAxis);
    QPen pen;
    pen.setColor(QColor(color));
    pen.setStyle(penStyle);
    pen.setWidthF(width);
    curve->setPen(pen);

    for(const auto& point: points_NURBS)
        curve->addData(point.curve.first, point.curve.second);

    curve->setName(label);
    m_canvas->replot();
}

// Рисует надпись
void Chart::paint_lable(double x, double y, const QString &text, double font_size) const
{
    QCPItemText *label = new QCPItemText(m_canvas);
    label->setFont(QFont("sans", font_size));
    label->position->setCoords(x, y);
    label->setText(text);
    m_canvas->replot();
}

// Рисует надпись со стрелкой
void Chart::paint_lable_with_arrow(double x_start, double y_start, double x_end, double y_end, const QString &text) const
{
    QCPItemText *label = new QCPItemText(m_canvas);
    label->setFont(QFont("sans", 10));
    label->position->setCoords(x_start + 0.3, y_start + 0.7);
    label->setText(text);

    QCPItemLine *line = new QCPItemLine(m_canvas);
    line->setHead(QCPLineEnding::esFlatArrow);
    line->start->setCoords(x_start, y_start);
    line->end->setCoords(x_end, y_end);
    m_canvas->replot();
}

// Рисует вектор первой производной кривой
void Chart::paint_first_deriv(const NURBS &points_NURBS, const QString &label, const Qt::PenStyle &penStyle, const QColor &color, double width) const
{
    QCPCurve *curve = new QCPCurve(m_canvas->xAxis, m_canvas->yAxis);
    QPen pen;
    pen.setColor(QColor(color));
    pen.setStyle(penStyle);
    pen.setWidthF(width);
    curve->setPen(pen);

    for(const auto& point: points_NURBS.get_point_curve())
        curve->addData(point.deriv_1.first, point.deriv_1.second);

    curve->setName(label);
    m_canvas->replot();
}

// Рисует вектор второй производной кривой
void Chart::paint_second_deriv(const NURBS &points_NURBS, const QString &label, const Qt::PenStyle &penStyle, const QColor &color, double width) const
{
    QCPCurve *curve = new QCPCurve(m_canvas->xAxis, m_canvas->yAxis);
    QPen pen;
    pen.setColor(QColor(color));
    pen.setStyle(penStyle);
    pen.setWidthF(width);
    curve->setPen(pen);

    for(const auto& point: points_NURBS.get_point_curve())
        curve->addData(point.deriv_2.first, point.deriv_2.second);

    curve->setName(label);
    m_canvas->replot();
}

// Рисует аппроксимирующие кривые, на расстоянии эпсилон от исходной кривой
void Chart::paint_curves_shift(const QPair<QVector<Point_curve>, QVector<Point_curve>> &points_NURBS) const
{
    paint_curve(points_NURBS.first, "Граница", Qt::PenStyle::DashLine);
    paint_curve(points_NURBS.second, "", Qt::PenStyle::DashLine);
}
