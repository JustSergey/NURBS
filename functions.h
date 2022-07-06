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

// Вычисляет длину для радиус вектора
double vector_len(const QPair<double, double>& point)
{
    return sqrt(pow(point.first, 2) + pow(point.second, 2));
}

// Вычисляет длину для радиус вектора
double vector_len(const double& x, const double& y)
{
    return sqrt(pow(x, 2) + pow(y, 2));
}

// Вычисляет длину для вектора по координатам
double vector_len(const QPair<double, double>& p1, const QPair<double, double>& p2)
{
    return sqrt(pow(p2.first - p1.first, 2) + pow(p2.second - p1.second, 2));
}

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

void curve_point_and_deriv_NURBS(Point_curve& data_NURBS, const int& n, const int& p, const std::vector<double>& u, const QVector<QVector<double>>& b, const std::vector<double>& h,
                                 const double& u_i, std::vector<QPair<double, double>>& c2,  std::vector<std::vector<double>> nders)
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
    data_NURBS.span = span;

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
        qDebug() << "nders[j][i] =" << nders[j][i] << " b[span - p + i] =" << b[span - p + i] << " h[span - p + i] =" << h[span - p + i];

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

// Рассчитывает спаны реального диапазона узлового вектора
QVector<double> real_span_calc(const uint& p, const uint& n, const std::vector<double>& u)
{ // Возвращает вектор с точками спанов реального узлового вектора
    // Реальный диапазон
    double u_start = u[p];
    const double u_end = u[n + 1];

    QVector<double> u_real_span;

    for(int i = 1; u_start < u_end; ++i)
    {
        u_real_span.push_back(u_start);
        u_start = u[p + i];
    }

    u_real_span.push_back(u_end);

    return u_real_span;
}

/*
// Возвращает точку кривой, перпендикулярной точке на плоскости
Point_curve finding_perpendicular(const int& n, const int& p, const std::vector<double>& u_vector, const QVector<QVector<double>>& polygon, const std::vector<double>& h, const QPair<double, double>& point)
{
    QVector<Point_curve> point_u(n - 1); // Массив точек - перпендикуляров
    QVector<double> u_real_span = real_span_calc(p, n, u_vector); // Спаны реального диапазона узлового вектора

    for(int i = 0; i < u_real_span.size() - 1; ++i)
    {
        point_u[i].u = (u_real_span[i + 1] - u_real_span[i]) / 2 + u_real_span[i]; // Берём среднее спана
        std::vector<std::vector<double>> nders(p + 1, std::vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и 1-я и 2-я производные
        std::vector<QPair<double, double>> c2(p + 1); // Индекс 2 для 2D задачи

        for(int k = 0; k < 1000; ++k) // ПОКА 35 ИТЕРАЦИЙ
        {

            if(point_u[i].u < u_real_span[i]) // Если точка вышла из спана
            {
                point_u[i].u = u_real_span[i];
                curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
                point_u[i].curve = c2[0];
                point_u[i].derivative_1 = c2[1];
                point_u[i].derivative_2 = c2[2];
                continue;
            }
            else if(point_u[i].u > u_real_span[i + 1])
            {
                point_u[i].u = u_real_span[i + 1];
                curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
                point_u[i].curve = c2[0];
                point_u[i].derivative_1 = c2[1];
                point_u[i].derivative_2 = c2[2];
                continue;
            }

            curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);

            point_u[i].curve = c2[0];
            point_u[i].derivative_1 = c2[1];
            point_u[i].derivative_2 = c2[2];

            double x = point_u[i].curve.first - point.first;
            double y = point_u[i].curve.second - point.second;
            double numerator = x * point_u[i].derivative_1.first + y * point_u[i].derivative_1.second;
            double denominator = x * point_u[i].derivative_2.first + y * point_u[i].derivative_2.second + pow(vector_len(point_u[i].derivative_1), 2);
            point_u[i].u = point_u[i].u - numerator / denominator * 0.01; // Новая, приближённая точка кривой
        }
    }

    Point_curve point_min_len; // Точка кривой с минимальным расстоянием до точки на плоскости
    point_min_len = point_u[0]; // Присваиваем первую точку для дальнейшего сравнения
    double min_len = vector_len(point, point_u[0].curve); // Минимальная длина вектора

    for(int i = 1; i < point_u.size(); ++i) // Ищем вектор с минимальной длиной
    {
        double temp_len = vector_len(point, point_u[i].curve);

        if(min_len > temp_len)
        {
            min_len = temp_len;
            point_min_len = point_u[i];
        }
    }

    return point_min_len;
}
*/

// Возвращает точку кривой, перпендикулярной точке на плоскости
Point_curve finding_perpendicular(const int& n, const int& p, const std::vector<double>& u_vector, const QVector<QVector<double>>& polygon, const std::vector<double>& h, const QPair<double, double>& point)
{
    QVector<Point_curve> point_u(n - 1); // Массив точек - перпендикуляров
    QVector<double> u_real_span = real_span_calc(p, n, u_vector); // Спаны реального диапазона узлового вектора

    for(int i = 0; i < u_real_span.size() - 1; ++i) // Итерируемся по спанам, начиная с нулевого
    {
        point_u[i].u = (u_real_span[i + 1] - u_real_span[i]) / 2 + u_real_span[i]; // Берём среднее спана
        std::vector<std::vector<double>> nders(p + 1, std::vector<double>(p + 1)); // nders - для заданного "u" массив BASIS функций и 1-я и 2-я производные
        std::vector<QPair<double, double>> c2(p + 1); // Индекс 2 для 2D задачи

        curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
        point_u[i].curve = c2[0];
        point_u[i].derivative_1 = c2[1];
        point_u[i].derivative_2 = c2[2];

        double cosine_check = 1;
        const double EPSILON_ANGLE = 0.001; // Эпсилон для косинуса прямого угла

        do
        {
            double x = point_u[i].curve.first - point.first;
            double y = point_u[i].curve.second - point.second;
            double numerator = x * point_u[i].derivative_1.first + y * point_u[i].derivative_1.second;
            double denominator = x * point_u[i].derivative_2.first + y * point_u[i].derivative_2.second + pow(vector_len(point_u[i].derivative_1), 2);
            point_u[i].u = point_u[i].u - numerator / denominator; // Новая, приближённая точка кривой

            curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);

            point_u[i].curve = c2[0];
            point_u[i].derivative_1 = c2[1];
            point_u[i].derivative_2 = c2[2];

            if(point_u[i].u < u_real_span[i]) // Если точка вышла из спана
            {
                point_u[i].u = u_real_span[i];
                curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
                point_u[i].curve = c2[0];
                point_u[i].derivative_1 = c2[1];
                point_u[i].derivative_2 = c2[2];
                break;
            }
            else if(point_u[i].u > u_real_span[i + 1])
            {
                point_u[i].u = u_real_span[i + 1];
                curve_point_and_deriv_NURBS(point_u[i], n, p, u_vector, polygon, h, point_u[i].u, c2, nders);
                point_u[i].curve = c2[0];
                point_u[i].derivative_1 = c2[1];
                point_u[i].derivative_2 = c2[2];
                break;
            }

            // Проверка на нулевой косинус
            numerator = abs(numerator);
            denominator = vector_len(point_u[i].derivative_1) * vector_len(x, y);
            cosine_check = numerator / denominator; // Угол между точкой и u

        } while (cosine_check > EPSILON_ANGLE);
    }

    Point_curve point_min_len; // Точка кривой с минимальным расстоянием до точки на плоскости
    point_min_len = point_u[0]; // Присваиваем первую точку для дальнейшего сравнения
    double min_len = vector_len(point, point_u[0].curve); // Минимальная длина вектора

    for(int i = 1; i < point_u.size(); ++i) // Ищем вектор с минимальной длиной
    {
        double temp_len = vector_len(point, point_u[i].curve);

        if(min_len > temp_len)
        {
            min_len = temp_len;
            point_min_len = point_u[i];
        }
    }

    return point_min_len;
}

#endif // FUNCTIONS_H
