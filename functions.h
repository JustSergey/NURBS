#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <QDebug>

struct Point_curve // Хранит точку кривой, её производную (1 и 2), интервал (span) и параметр u
{
    QPair<double, double> curve;
    QPair<double, double> derivative_1;
    QPair<double, double> derivative_2;
    int span;
    double u;
};

// Определяет индекс узлового промежутка (интервал)
uint findSpan(const uint& n, const int& p, const std::vector<double>& u, const double& u_i)
/*
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m = n + 1 + p
 * u_i - точка внутри реального диатазона в узловом векторе
*/
{
    for(uint k = 0; k < u.size() - 1; ++k)
    {
        if(u[k] > u[k + 1])
            qDebug() << "Сообщение из findSpan - Узловой вектор убывает u[k] > u[k + 1]";
    }

    if((u.size() - 1) != (n + 1 + p))
        qDebug() << "Сообщение из findSpan - (u.size() - 1) != (n + 1 + p)";

    if(u_i < u[p] || u_i > u[n + 1])
        qDebug() << "Сообщение из FindSpan - u вышел за реальный диапазон u_i < u[p] || u_i > u[n + 1]";

    if(u_i == u[n + 1]) // Последний диапазон (начало диапазона), в котором может находиться u
        return n;

    uint low = p, high = n + 1, middle = (low + high) / 2;

    // Выполняем  двоичный  поиск
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

// Вычисляет все ненулевые базисные функции и производные (от 0 до p)
void dersBasisFuns(const double& i, const double& u_i, const int& p, const std::vector<double>& u, std::vector<std::vector<double>>& nders)
/*
 * i - номер диапазона для которого расчитывется N (N_i-p,p ... N_i,p)
 * u_i - значение u
 * p - степень полинома (=degree)
 * u - узловой вектор
 * nders - массив basis функций - N[0],...,N[p] и их производных
*/
{
    QVector<double> left(p + 1), right(p + 1);
    QVector<QVector<double>> ndu(p + 1, QVector<double>(p + 1)); // Для хранения базисных функций и узлов различия
    QVector<QVector<double>> a(2, QVector<double>(p + 1)); // Хранит два наиболее недавно вычисленных ряда

    ndu[0][0] = 1.0;

    for(int j = 1; j < p + 1; ++j)
    {
        left[j] = u_i - u[i + 1 - j];
        right[j] = u[i + j] - u_i;
        double saved = 0;

        for(int r = 0; r < j; ++r)
        {
            // Нижний треугольник
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = ndu[r][j - 1] / ndu[j][r];
            // Верхний треугольник
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        ndu[j][j] = saved;
    }

    for(int j = 0; j <= p; ++j) // Загрузим базисные функции
        nders[0][j] = ndu[j][p];

    // В этом разделе вычисляем производные

    for(int r = 0; r < p + 1; ++r) // Цикл по индексу функции
    {
        int s1 = 0, s2 = 1; // Альтернативные строки в массиве
        a[0][0] = 1.0;

        for(int k = 1; k <= p; ++k)
        {
            double d = 0;
            double rk = r - k;
            double pk = p - k;

            if(r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            double j1 = 0;

            if(rk >= -1)
                j1 = 1;
            else
                j1 = -rk;

            double j2 = 0;

            if(r - 1 <= pk)
                j2 = k - 1;
            else
                j2 = p - r;

            for(uint j = j1; j <= j2; ++j)
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

            // Меняем строки местами
            double temp = s1;
            s1 = s2;
            s2 = temp;
        }
    }

    // Умножаем на правильные коэффициенты

    double r = p;

    for(int k = 1; k <= p; ++k)
    {
        for(int j = 0; j < p + 1; ++j)
        {
            nders[k][j] *= r;
        }

        r *= p - k;
    }

    // Для контроля суммируем значения базисных Функций в точке "u".
    // Если все верно, то сумма должна быть = 1

    double sum = 0;
    for(uint i = 0; i < nders.size(); ++i)
        sum += nders[0][i];

    if((sum < (1 - 1e-10)) || (sum > 1 + 1e-10))
        qDebug() << "Сообщение из DersBasisFuns - Сумма базисных Функций != 1";
}

void curve_point_and_deriv_NURBS(QVector<Point_curve>& data_NURBS, const int& n, const int& p, const std::vector<double>& u, const QVector<QVector<double>>& b,
                                 const std::vector<double>& h, const double& u_i, std::vector<QPair<double, double>>& c2,  std::vector<std::vector<double>> nders)
/*
 * Функция расчитывает для заданного "u" одну точку на В-сплайне и 1-ю и 2-ю проиизв. для этой точки
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * b - контрольные точки (control polygon)
 * u_i - точка внутри РЕАЛЬНОГО диатазона в узловом векторе
*/
{
    double span = findSpan(n, p, u, u_i); // Диапазон узлового веткора

    static int counter; // Счётчик для присваивания нужного span
    data_NURBS[counter].span = span;

    if(counter == data_NURBS.size() - 1) // Если мы дошли до конца массива
        counter = 0; // Обнуляем счётчик
    else
        ++counter;

    qDebug() << "Span =" << span << "\tu =" << u_i;

    if ((b.size() - 1) != n)
        qDebug() << "** Сообщение из curvePoin_and_Deriv_NURBS -- (b[0].size() - 1) != n";

    dersBasisFuns(span, u_i, p, u, nders);

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
         s1[i] = n3[i] / d - (n1[i] * n2[i]) / (pow(d, 2));

    std::vector<double> nn(2);

    for(size_t i = 0; i < nn.size(); ++i)
        nn[i] = n1[i] * n2[i];

    std::vector<double> nn_deriv(2);

    for(size_t i = 0; i < nn_deriv.size(); ++i)
        nn_deriv[i] = n1[i] * n2[i] + n0[i] * n4[i];

    std::vector<double> s2(2);

    for(size_t i = 0; i < s2.size(); ++i)
        s2[i] = nn_deriv[i] / (d * d) - (nn[i] * 2 * n2[i]) / (pow(d, 4));

    c2[2].first = s1[0] - s2[0];
    c2[2].second = s1[1] - s2[1];

    qDebug() << "j = 2 - 2-я производная \nc2 =" << c2 << "\n";

    return;
}

// Вычисляет длину вектора (2D)
double vector_len(const QPair<double, double>& point)
{
    return sqrt(pow(point.first, 2) + pow(point.second, 2));
}

#endif // FUNCTIONS_H
