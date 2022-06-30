#include "widget.h"
#include "ui_widget.h"
#include <QDebug>

void printError(const QString& error_message)
{
    qDebug() << error_message;
    return;
}

double findSpan(const int& n, const int& p, const std::vector<double>& u, const double& u_i)
/* n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * u_i - точка внутри реального диатазона в узловом векторе */
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
    int low {p}, high {n + 1};
    uint middle = static_cast<uint>((low + high) / 2);

    while((u_i < u[middle]) || (u_i >= u[middle + 1]))
    {
        if(u_i < u[middle])
            high = middle;
        else
            low = middle;

        middle = static_cast<uint>((low + high) / 2);
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

    double sum { 0 };
    for(uint i = 0; i < nders.size(); ++i)
        sum += nders[0][i];

    if((sum < (1 - 1e-10)) || (sum > 1 + 1e-10))
        printError("** Сообщение из DersBasisFuns -- Сумма Базис. Функций НЕ РАВНА 1");
}

void curve_point_and_deriv_NURBS(const int& n, const int& p, const std::vector<double>& u,
                               const QVector<QVector<double>>& b, const std::vector<double>& h,
                               const double& u_i, std::vector<std::vector<double>>& c2,
                               std::vector<std::vector<double>> nders)
/* Функция расчитывает для заданного "u" одну точку на В-сплайне и 1-ю и 2-ю проиизв. для этой точки
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * b - контрольные точки (control polygon)
 * u_i - точка внутри РЕАЛЬНОГО диатазона в узловом векторе */
{
    double span = findSpan(n, p, u, u_i);
    qDebug() << "Span =" << span << "\tu =" << u_i;

    if ((b.size() - 1) != n)
        printError("** Сообщение из curvePoin_and_Deriv_NURBS -- (b[0].size() - 1) != n");

    dersBasisFuns(span, u_i, p, u, nders);

    c2.assign(p + 1, std::vector<double>(2));

    double d  { 0 }; // Знаменатель формулы NURBS (формула 5-122, Роджерс (рус.) стр 360)
    std::vector<double> n0(2); // Числитель формулы NURBS (формула 5-122, Роджерс (рус.) стр 360)
    std::vector<double> n1(2); // Числитель Первого слагаемого формулы 1-ой Производ. NURBS (формула 5-126, Роджерс (рус.) стр 372)
    std::vector<double> n2(2); // Множитель в Числителе Второго слагаемого формулы 1-ой Производ. NURBS (формула 5-126, Роджерс (рус.) стр 372)
    std::vector<double> n3(2); // Множитель в Числителе при расчёте 2-ой Производ. NURBS ((см. мои листы)
    std::vector<double> n4(2); // Мночитель в Числителе при расчёте 2-ой Производ. NURBS ((см. мои листы)

    int j { 0 }; // кривая (нулевая производная)

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

    j = 1; // 1-я производная

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

    j = 2; // 2-я производная
    
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

void plot_trace(const QVector<QVector<double>>& data_point, const QVector<QVector<QVector<double>>>& data_spline,
                const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                const QString& title, const QString& labels_legend_1, const QString& labels_legend_2, Ui::Widget* ui)
{
    ui->widget->clearGraphs(); // Очищаем все графики

    QCPCurve *curve_point { nullptr };
    curve_point = new QCPCurve(ui->widget->xAxis, ui->widget->yAxis);

    curve_point->setPen(QColor(0, 0, 0, 255)); // Задаем чёрный цвет
    curve_point->setLineStyle(QCPCurve::lsNone); // Убираем линии
    curve_point->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 6)); // Формируем вид точек

    for(const auto& el: data_point) // Рисуем точки
        curve_point->addData(el[0], el[1]);

    curve_point->setLineStyle(QCPCurve::lsLine); // Добавляем линии

    //ui->widget->addGraph(); // Добавляем ещё один график в widget

    QCPCurve *curve_spline { nullptr };
    curve_spline = new QCPCurve(ui->widget->xAxis, ui->widget->yAxis);

    for(const auto& el: data_spline) // Рисуем сплайн
        curve_spline->addData(el[0][0], el[0][1]);

    ui->widget->setInteractions(QCP :: iRangeDrag | QCP :: iRangeZoom); // Перетаскиваемый + масштабирование колеса прокрутки

    ui->widget->legend->setVisible(true); // Включаем Легенду графика
    curve_point->setName(labels_legend_1);
    curve_spline->setName(labels_legend_2);

    // Подписываем оси Ox и Oy
    ui->widget->xAxis->setLabel("Ось X");
    ui->widget->yAxis->setLabel("Ось Y");

    // Установим область, которая будет показываться на графике
    ui->widget->xAxis->setRange(x_min, x_max); // Для оси Ox
    ui->widget->yAxis->setRange(y_min, y_max); // Для оси Oy

    ui->widget->plotLayout()->insertRow(0); // Вставляем строку
    ui->widget->plotLayout()->addElement (0, 0, new QCPTextElement(ui->widget, title, QFont("sans", 12))); // Добавляем в первую строку и первый столбец заглавие

    ui->widget->replot();
}

Widget::Widget(QWidget *parent)
    : QWidget(parent), ui(new Ui::Widget)
{
    ui->setupUi(this);

    using namespace std;

    const QVector<QVector<double>> b // Массив точек многоугольника
    {
        {1.0, 1.0},
        {2.0, 3.0},
        {4.0, 3.0},
        {3.0, 1.0}
    };

    vector<double> h {1, 1, 1, 1, 1, 1}; // Весовые коэффициенты
    vector<double> u {0, 0, 0, 0, 1, 1, 1, 1}; // Узловой вектор

    const uint p { 3 }; // Степень аппроксимирующих полиномов
    const uint m { 7 }; // n_kn - количество узлов (длина) в узловом векторе
    const uint n = m - p - 1; // n_real - количество узлов (длина) реальной части узлового вектора

    // Реальный диапазон
    const double u_start = u[p];
    const double u_stop = u[n + 1];

    vector<vector<double>> nders(p + 1, vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и  1-я и 2-я производные
    vector<vector<double>> c2(p + 1, vector<double>(2)); // Индекс 2 для 2D задачи

    const int n_u = 60; // Кол-во разбиений (точек) в реальной части узлов. вектора

    QVector<QVector<QVector<double>>> data_CurvePoin_and_Deriv_NURBS(n_u + 1, QVector<QVector<double>>(p + 1, QVector<double>(2)));

    for(int i = 0; i < n_u + 1; ++i)
    {
        double u_i = (i / static_cast<double>(n_u)) * (u_stop - u_start);
        curve_point_and_deriv_NURBS(n, p, u, b, h, u_i, c2, nders);

        for(size_t k = 0; k < c2.size(); ++k)
        {
            for(size_t l = 0; l < c2[k].size(); ++l)
                data_CurvePoin_and_Deriv_NURBS[i][k][l] = c2[k][l];
        }
    }

    QString title = "NURBS 4-го порядка (кубич. полиномы)";
    QString labels_legend_1 = "Контур. многоуг.";
    QString labels_legend_2 = "B-сплайн";

    int x_min = -1, x_max = 6;
    int y_min = -1, y_max = 6;

    plot_trace(b, data_CurvePoin_and_Deriv_NURBS, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2, ui);

}


Widget::~Widget()
{
    delete ui;
}
