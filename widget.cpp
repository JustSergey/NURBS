#include "widget.h"
#include "ui_widget.h"
#include <QDebug>

void printError(const QString& error_message)
{
    qDebug() << error_message;
    return;
}

// Определяет индекс узлового промежутка
uint findSpan(const int& n, const int& p, const std::vector<double>& u, const double& u_i)
/*
 * Вход: n,p,u,u_i. Выход: индекс  узлового  промежутка
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * u_i - точка внутри реального диатазона в узловом векторе
*/
{
    for(uint k = 0; k < u.size() - 1; ++k)
    {
        if(u[k] > u[k + 1])
            printError("** Сообщение из findSpan -- Узловой вектор убывает (u[k] > u[k + 1])");
    }

    if((static_cast<int>(u.size()) - 1) != (n + 1 + p)) // Кастуем в int, чтобы не было предупреждения
        printError("** Сообщение из findSpan -- (u.size() - 1) не равен (n + 1 + p)");

    if(u_i < u[p] || u_i > u[n + 1])
        printError("** Сообщение из FindSpan -- u вышел за реальный диапазон -- (u_i < u[p] || u_i > u[n + 1])");

    if(u_i == u[n + 1]) // Последний диапазон (начало диапазона), в котором может находиться u
        return n;

    // Далее идёт двочиный поиск
    uint low = p, high = n + 1;
    uint middle = static_cast<uint>((low + high) / 2);

    while((u_i < u[middle]) || (u_i >= u[middle + 1]))
    {
        if(u_i < u[middle])
            high = middle;
        else
            low = middle;

        middle = (low + high) / 2;
    }

    return middle;
}

void dersBasisFuns(const double& i, const double& u_i, const int& p, const std::vector<double>& u,
                   std::vector<std::vector<double>>& nders)
/* Функция расчитывает все ненулевые базис. функции (НЕ Нулевые - в смысле только, например, 4-ре для кубич. полиномов)
 * для заданного аргумента ("u")  и все производные (от 0 до p) помещает их в 2D массив ders
 *
 * i - номер диапазона для которого расчитывется N (N_i-p,p ... N_i,p)
 * u_i - значение u
 * p - степень полинома (=degree)
 * u - узловой вектор
 * nders - массив basis функций - N[0],...,N[p] и их производных */
{
    using namespace std;

    vector<double> left(p + 1), right(p + 1);
    vector<vector<double>> ndu(p + 1, vector<double>(p + 1)), a(p + 1, vector<double>(p + 1));

    ndu[0][0] = 1.0;

    for(int j = 1; j < p + 1; ++j)
    {
        left[j] = u_i - u[i + 1 - j];
        right[j] = u[i + j] - u_i;
        double saved = 0;

        for(int r = 0; r < j; ++r)
        {
            // Lower triangle
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = ndu[r][j - 1] / ndu[j][r];
            // Upper triangle
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        ndu[j][j] = saved;
    }

    // Load the basis functions

    for(int j = 0; j < p + 1; ++j)
        nders[0][j] = ndu[j][p];

    // This section computes the derivatives (Eq. [2.9])

    for(int r = 0; r < p + 1; ++r)
    {
        int s1 { 0 }, s2 { 1 }; // Alternate rows in array "a"
        a[0][0] = 1.0;

        for(int k = 1; k < p + 1; ++k)
        {
            double d { 0.0 };
            double rk = r - k;
            double pk = p - k;

            if(r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            double j1 { 0 };

            if(rk >= -1)
                j1 = 1;
            else
                j1 = -rk;

            double j2 { 0 };

            if(r - 1 <= pk)
                j2 = k - 1;
            else
                j2 = p - r;

            for(uint j = j1; j < j2 + 1; ++j)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if(r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }

            nders[k][r] = d;
            double temp = s1;
            s1 = s2;
            s2 = temp; // Switch rows
        }
    }

    // Multiply through by the correct factors

    double r = p;

    for(int k = 1; k < p + 1; ++k)
    {
        for(int j = 0; j < p + 1; ++j)
        {
            nders[k][j] *= r;
        }

        r *= p - k;
    }

    // Для контроля Суммируем значения Базис. Функций в точке "u".
    // Если все верно, то Сумма долж. быть = 1

    double sum = 0;
    for(uint i = 0; i < nders.size(); ++i)
        sum += nders[0][i];

    if((sum < (1 - 1e-10)) || (sum > 1 + 1e-10))
        printError("** Сообщение из DersBasisFuns -- Сумма Базис. Функций НЕ РАВНА 1");
}

void curve_point_and_deriv_NURBS(const int& n, const int& p, const std::vector<double>& u,
                               const QVector<QVector<double>>& b, const std::vector<double>& h,
                               const double& u_i, std::vector<std::vector<double>>& c2,
                               std::vector<std::vector<double>> nders, const QVector<double>& point_u, QVector<int>& index_u)
/* Функция расчитывает для заданного "u" одну точку на В-сплайне и 1-ю и 2-ю проиизв. для этой точки
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * b - контрольные точки (control polygon)
 * u_i - точка внутри РЕАЛЬНОГО диатазона в узловом векторе */
{
    double span = findSpan(n, p, u, u_i); // Диапазон узлового веткора
    qDebug() << "Span =" << span << "\tu =" << u_i;

    static uint counter; // Для отсчёта индекса в массиве u

    for(const auto& el: point_u)
    {
        if(u_i == el)
        {
            index_u.push_back(counter);
            break;
        }
    }

    ++counter;

    if ((b.size() - 1) != n)
        printError("** Сообщение из curvePoin_and_Deriv_NURBS -- (b[0].size() - 1) != n");

    dersBasisFuns(span, u_i, p, u, nders);

    c2.assign(p + 1, std::vector<double>(2));

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

    for(size_t i = 0; i < c2[0].size(); ++i)
        c2[0][i] = n0[i] / d;

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

    for(size_t i = 0; i < c2[0].size(); ++i)
        c2[1][i] = n1[i] / d - (n0[i] * n2[i]) / (d * d);

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

    for(size_t i = 0; i < c2[2].size(); ++i)
        c2[2][i] = s1[i] - s2[i];

    qDebug() << "j = 2 - 2-я производная \nc2 =" << c2 << "\n";

    return;
}

void plot_deriv_for_point(const QVector<QVector<QVector<double>>>& pointDeriv, Ui::Widget* ui)
{
    for(const auto& point: pointDeriv)
    {
        QCPItemLine *line = new QCPItemLine(ui->graph_function);
        line->setHead(QCPLineEnding::esFlatArrow);
        line->start->setCoords(point[0][0], point[0][1]);
        line->end->setCoords(point[1][0] + point[0][0], point[1][1] + point[0][1]);
    }

    ui->graph_function->replot();

}

void function_plot(const QVector<QVector<double>>& data_point, const QVector<QVector<QVector<double>>>& data_spline,
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
        curve_spline->addData(el[0][0], el[0][1]);

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

void first_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<QVector<QVector<double>>>& data_deriv,
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
        curve_deriv->addData(el[1][0], el[1][1]);

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

void second_derivative_plot(const QVector<QVector<double>>& data_point, const QVector<QVector<QVector<double>>>& data_deriv,
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
        curve_deriv->addData(el[2][0], el[2][1]);

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

    vector<double> h {1, 1, 1, 1, 1}; // Весовые коэффициенты
    vector<double> u {0, 0, 0, 0.4, 0.6, 1, 1, 1}; // Узловой вектор

    const uint p = 2; // Степень аппроксимирующих полиномов
    const uint m = u.size() - 1; // n_kn - количество узлов (длина) в узловом векторе
    const uint n = m - p - 1; // n_real - количество узлов (длина) реальной части узлового вектора

    // Реальный диапазон
    const double u_start = u[p];
    const double u_stop = u[n + 1];

    vector<vector<double>> nders(p + 1, vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и  1-я и 2-я производные
    vector<vector<double>> c2(p + 1, vector<double>(2)); // Индекс 2 для 2D задачи

    const int n_u = 60; // Кол-во разбиений (точек) в реальной части узлов. вектора

    QVector<QVector<QVector<double>>> data_CurvePoin_and_Deriv_NURBS(n_u + 1, QVector<QVector<double>>(p + 1, QVector<double>(2)));

    QVector<int> index_u; // Массив, хранящий индексы u в
    QVector<double> point_u(b.size()); // Массив, хранящий u, от которых пойдёт производная

    for(double i = b.size() - 1; i >= 1; --i)
        point_u[b.size() - i] = 1 / i;

    for(int i = 0; i < n_u + 1; ++i)
    {
        double u_i = (i / static_cast<double>(n_u)) * (u_stop - u_start);
        curve_point_and_deriv_NURBS(n, p, u, b, h, u_i, c2, nders, point_u, index_u);

        for(size_t k = 0; k < c2.size(); ++k)
        {
            for(size_t l = 0; l < c2[k].size(); ++l)
                data_CurvePoin_and_Deriv_NURBS[i][k][l] = c2[k][l];
        }
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

    QVector<QVector<QVector<double>>> pointDeriv;

    for(const auto& index: index_u)
        pointDeriv.push_back(data_CurvePoin_and_Deriv_NURBS[index]);

    plot_deriv_for_point(pointDeriv, ui);
}


Widget::~Widget()
{
    delete ui;
}
