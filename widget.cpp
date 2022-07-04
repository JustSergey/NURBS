#include "widget.h"
#include "ui_widget.h"
#include "functions.h"
#include <QDebug>

void curve_point_and_deriv_NURBS(QVector<curve>& data_CurvePoin_and_Deriv_NURBS, const int& n, const int& p, const std::vector<double>& u,
                               const QVector<QVector<double>>& b, const std::vector<double>& h,
                               const double& u_i, std::vector<QPair<double, double>>& c2,
                               std::vector<std::vector<double>> nders)
/* Функция расчитывает для заданного "u" одну точку на В-сплайне и 1-ю и 2-ю проиизв. для этой точки
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * b - контрольные точки (control polygon)
 * u_i - точка внутри РЕАЛЬНОГО диатазона в узловом векторе */
{
    double span = findSpan(data_CurvePoin_and_Deriv_NURBS, n, p, u, u_i); // Диапазон узлового веткора
    qDebug() << "Span =" << span << "\tu =" << u_i;

    if ((b.size() - 1) != n)
        qDebug() << "** Сообщение из curvePoin_and_Deriv_NURBS -- (b[0].size() - 1) != n";

    dersBasisFuns(span, u_i, p, u, nders);

    //c2.resize(p + 1, std::vector<QPair<double, double>>);
    //c2.assign(p + 1);
     //vector<vector<QPair<double, double>>> c2(p + 1)

    double d  = 0; // Знаменатель формулы NURBS (формула 5-122, Роджерс (рус.) стр 360)
    std::vector<double> n0(2); // Числитель формулы NURBS (формула 5-122, Роджерс (рус.) стр 360)
    std::vector<double> n1(2); // Числитель Первого слагаемого формулы 1-ой Производ. NURBS (формула 5-126, Роджерс (рус.) стр 372)
    std::vector<double> n2(2); // Множитель в Числителе Второго слагаемого формулы 1-ой Производ. NURBS (формула 5-126, Роджерс (рус.) стр 372)
    std::vector<double> n3(2); // Множитель в Числителе при расчёте 2-ой Производ. NURBS ((см. мои листы)
    std::vector<double> n4(2); // Мночитель в Числителе при расчёте 2-ой Производ. NURBS ((см. мои листы)

    int j = 0; // кривая (нулевая производная)

    for(int i = 0; i < p + 1; ++i)
    {
        qDebug() << "----------";
        qDebug() << "j =" << j << " i =" << i << " span - p + i =" << span - p + i;
        qDebug() << "nders[j][i] =" << nders[j][i] << " b[span - p + i] =" << b[span - p + i] <<
                    " h[span - p + i] =" << h[span - p + i];

        for(int k = 0; k < b[0].size(); ++k)
            n0[k] += b[span - p + i][k] * h[span - p + i] * nders[j][i];

        d += nders[j][i] * h[span - p + i];
    }

    c2[0].first = n0[0] / d;
    c2[0].second = n0[1] / d;

    qDebug() << "j = 0 - кривая \nc2 =" << c2;

    if(p == 1)
        return;

    // 1-я производная

    for(int i = 0; i < p + 1; ++i)
    {
        for(size_t k = 0; k < n1.size(); ++k)
        {
            n1[k] += b[span - p + i][k] * h[span - p + i] * nders[1][i];
            n2[k] += h[span - p +i] * nders[1][i];
        }
    }

    c2[1].first = n1[0] / d - (n0[0] * n2[0]) / (d * d);
    c2[1].second = n1[1] / d - (n0[1] * n2[1]) / (d * d);

    qDebug() << "j = 1 - 1-я производная \nc2 =" << c2;

    // 2-я производная
    
    for(int i = 0; i < p + 1; ++i)
    {
        for(size_t k = 0; k < n1.size(); ++k)
        {
            n3[k] += b[span - p + i][k] * h[span - p + i] * nders[2][i];
            n4[k] += h[span - p +i] * nders[2][i];
        }
    }
    
    std::vector<double> s1(2);
    
    for(size_t i = 0; i < s1.size(); ++i)
         s1[i] = n3[i] / d - (n1[i] * n2[i]) / (d * d);
    
    std::vector<double> nn(2);

    for(size_t i = 0; i < nn.size(); ++i)
        nn[i] = n1[i] * n2[i];

    std::vector<double> nn_deriv(2);

    for(size_t i = 0; i < nn_deriv.size(); ++i)
        nn_deriv[i] = n1[i] * n2[i] + n0[i] * n4[i];
    
    std::vector<double> s2(2);

    for(size_t i = 0; i < s2.size(); ++i)
        s2[i] = nn_deriv[i] / (d * d) - (nn[i] * 2 * n2[i]) / (d * d * d * d);

    c2[2].first = s1[0] - s2[0];
    c2[2].second = s1[1] - s2[1];

    qDebug() << "j = 2 - 2-я производная \nc2 =" << c2 << "\n";

    return;
}

void plot_deriv_for_point(const QVector<curve>& pointDeriv, Ui::Widget* ui)
{
    for(const auto& point: pointDeriv)
    {
        QCPItemLine *line = new QCPItemLine(ui->graph_function);
        line->setHead(QCPLineEnding::esFlatArrow);
        line->start->setCoords(point.curve.first, point.curve.second);
        line->end->setCoords(point.derivative_1.first + point.curve.first, point.derivative_1.second + point.curve.second);
    }

    ui->graph_function->replot();
}


void function_plot(const QVector<QVector<double>>& data_point, const QVector<curve>& data_spline,
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

    for(const auto& el: data_spline) // Рисуем сплайн
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

void first_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<curve>& data_deriv,
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

    for(const auto& el: data_deriv) // Рисуем сплайн
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

void second_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<curve>& data_deriv,
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

    for(const auto& el: data_deriv) // Рисуем сплайн
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

    QVector<curve> data_CurvePoin_and_Deriv_NURBS(n_u + 1);

    for(int i = 0; i < n_u + 1; ++i)
    {
        double u_i = (i / static_cast<double>(n_u)) * (u_stop - u_start);
        curve_point_and_deriv_NURBS(data_CurvePoin_and_Deriv_NURBS, n, p, u, b, h, u_i, c2, nders);

        data_CurvePoin_and_Deriv_NURBS[i].u = u_i;
        data_CurvePoin_and_Deriv_NURBS[i].curve = c2[0];
        data_CurvePoin_and_Deriv_NURBS[i].derivative_1 = c2[1];
        data_CurvePoin_and_Deriv_NURBS[i].derivative_2 = c2[2];
    }

    QString title = "B-сплайн 4-го порядка";
    QString labels_legend_1 = "Контур. многоуг.";
    QString labels_legend_2 = "B-сплайн";

    int x_min = 0, x_max = 20;
    int y_min = -5, y_max = 15;

    function_plot(b, data_CurvePoin_and_Deriv_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    title = "1-я прoизвдная B-сплайна 4-го порядка";
    x_min = -10, x_max = 30;
    y_min = -25, y_max = 30;

    first_derivative_plot(b, data_CurvePoin_and_Deriv_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    title = "2-я прoизвдная B-сплайна 4-го порядка";
    x_min = -95, x_max = 60;
    y_min = -110, y_max = 120;

    second_derivative_plot(b, data_CurvePoin_and_Deriv_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

    QVector<double> point_u(b.size()); // Массив, хранящий u, от которых пойдёт производная на графике

    for(double i = b.size() - 1; i >= 1; --i)
        point_u[b.size() - i] = 1 / i;

    QVector<curve> data;
    for(uint i = 0; i < point_u.size() - 1;)
    {
        for(const auto& point: data_CurvePoin_and_Deriv_NURBS)
        {
            if(point_u[i] == point.u)
            {
                data.push_back(point);
                ++i;
            }
        }
    }

    plot_deriv_for_point(data, ui);
}


Widget::~Widget()
{
    delete ui;
}
