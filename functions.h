#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <QDebug>

// Хранит точку кривой, её производную (1 и 2), интервал (span) и параметр u
struct Point_curve
{
    QPair<double, double> curve;
    QPair<double, double> derivative_1;
    QPair<double, double> derivative_2;
    int span;
    double u;
};

// Возвращает заполненный узловой вектор с параметрами
QVector<double> u_fill(const QVector<QVector<double>>& control_points, const double& degree)
{
    const uint n_ver = control_points.size(); // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
    const uint n_kn = n_ver + degree + 1;     // Количество узлов (длина) в узловом векторе (n_knots)
    const uint n_real = n_ver - degree + 1;   // Количество узлов (длина) реальной части узлового вектора
    QVector<double> u(n_kn);                  // Узловой вектор

    // Реальный диапазон
    const double u_start = degree;
    const double u_stop = n_kn - degree - 1;
    const double step = 1 / static_cast<double>(n_real - 1); // Шаг в реальном диапазоне

    for(int i = u_start + 1; i < u_stop; ++i) // Заполняем реальный диапазон
        u[i] = u[i - 1] + step;

    for(int i = u_stop; i < u.size(); ++i)    // Заполняем последние параметры единицами
        u[i] = 1;

    return u;
}

// Возвращает длину радиус-вектора
double radius_vector_len(const QPair<double, double>& point)
{
    return sqrt(pow(point.first, 2) + pow(point.second, 2));
}

// Возвращает длину радиус-вектора
double radius_vector_len(const double& x, const double& y)
{
    return sqrt(pow(x, 2) + pow(y, 2));
}

// Возвращает длину вектора
double vector_len(const QPair<double, double>& point_1, const QPair<double, double>& point_2)
{
    return sqrt(pow(point_2.first - point_1.first, 2) + pow(point_2.second - point_1.second, 2));
}

// Возвращает длину вектора
double vector_len(const QPair<double, double>& point_1, const Point_curve& point_2)
{
    return sqrt(pow(point_2.curve.first - point_1.first, 2) + pow(point_2.curve.second - point_1.second, 2));
}

// Возвращает угол между двумя радиус-векторами
double vector_angle(const QPair<double, double>& vec_1, const QPair<double, double>& vec_2, const Point_curve& start_point)
{
    double a_x = vec_1.first - start_point.curve.first;
    double a_y = vec_1.second - start_point.curve.second;
    double b_x = vec_2.first - start_point.curve.first;
    double b_y = vec_2.second - start_point.curve.second;
    double numerator = a_x * b_x + a_y * b_y;
    double v1_len = vector_len(vec_1, start_point);
    double v2_len = vector_len(vec_2, start_point);
    double cos = numerator / (v2_len * v1_len);
    return acos(cos) * 180 / M_PI;
}

// Возвращает индекс узлового промежутка (интервал)
uint findSpan(const uint& n_kn, const uint& n_ver, const int& degree, const QVector<double>& u, const double& u_i)
{
    for(int k = 0; k < u.size() - 1; ++k)
    {
        if(u[k] > u[k + 1])
            qDebug() << "Сообщение из findSpan - Узловой вектор убывает u[k] > u[k + 1]";
    }

    if(static_cast<uint>(u.size()) != (n_ver + degree + 1))
        qDebug() << "Сообщение из findSpan - (u.size()) != (n_ver + degree + 1)";

    if(u_i < u[degree] || u_i > u[n_kn - degree - 1])
        qDebug() << "Сообщение из FindSpan - u вышел за реальный диапазон u_i < u[degree] || u_i > u[n_kn - degree - 1]";

    if(u_i == u[n_kn - degree - 1]) // Если точка дошла до последнего спана
        return n_kn - degree - 2;

    uint low = degree, high = n_kn - degree - 1, middle = (low + high) / 2;

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

// Вычисляет все ненулевые базисные функции и производные
void dersBasisFuns(const double& span, const double& u_i, const int& degree, const QVector<double>& u, QVector<QVector<double>>& nders)
{
    QVector<double> left(degree + 1), right(degree + 1);
    QVector<QVector<double>> ndu(degree + 1, QVector<double>(degree + 1)); // Для хранения базисных функций и узлов различия
    QVector<QVector<double>> a(2, QVector<double>(degree + 1)); // Хранит два наиболее недавно вычисленных ряда

    ndu[0][0] = 1.0;

    for(int j = 1; j < degree + 1; ++j)
    {
        left[j] = u_i - u[span + 1 - j];
        right[j] = u[span + j] - u_i;
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

    for(int j = 0; j <= degree; ++j) // Загрузим базисные функции
        nders[0][j] = ndu[j][degree];

    // В этом разделе вычисляем производные

    for(int r = 0; r < degree + 1; ++r) // Цикл по индексу функции
    {
        int s1 = 0, s2 = 1; // Альтернативные строки в массиве
        a[0][0] = 1.0;

        for(int k = 1; k <= degree; ++k)
        {
            double d = 0;
            double rk = r - k;
            double pk = degree - k;

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
                j2 = degree - r;

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

    double r = degree;

    for(int k = 1; k <= degree; ++k)
    {
        for(int j = 0; j < degree + 1; ++j)
        {
            nders[k][j] *= r;
        }

        r *= degree - k;
    }

    // Для контроля суммируем значения базисных Функций в точке "u"
    double sum = 0;

    for(int i = 0; i < nders.size(); ++i)
        sum += nders[0][i];

    if((sum < (1 - 1e-10)) || (sum > 1 + 1e-10)) // Если все верно, то сумма должна быть = 1
        qDebug() << "Сообщение из DersBasisFuns - Сумма базисных Функций != 1";
}

// Расчитывает для заданного "u" одну точку на В-сплайне и 1-ю и 2-ю проиизв. для этой точки
void curve_point_and_deriv_NURBS(Point_curve& data_NURBS, const uint& n_kn, const int& degree, const QVector<double>& u, const QVector<QVector<double>>& control_points, const QVector<double>& w, const double& u_i)
{
    const double n_ver = control_points.size();

    double span = findSpan(n_kn, n_ver, degree, u, u_i); // Диапазон узлового веткора
    data_NURBS.span = span;

    //qDebug() << "Span =" << span << "\tu =" << u_i;

    QVector<QVector<double>> nders(degree + 1, QVector<double>(degree + 1)); // Содержит для заданного "u" массив BASIS функций и 1-ую и 2-ую производную

    dersBasisFuns(span, u_i, degree, u, nders);

    double d  = 0; // Знаменатель формулы NURBS (формула 5-122, Роджерс (рус.) стр 360)
    std::vector<double> n0(2); // Числитель формулы NURBS (формула 5-122, Роджерс (рус.) стр 360)
    std::vector<double> n1(2); // Числитель Первого слагаемого формулы 1-ой Производ. NURBS (формула 5-126, Роджерс (рус.) стр 372)
    std::vector<double> n2(2); // Множитель в Числителе Второго слагаемого формулы 1-ой Производ. NURBS (формула 5-126, Роджерс (рус.) стр 372)
    std::vector<double> n3(2); // Множитель в Числителе при расчёте 2-ой Производ. NURBS
    std::vector<double> n4(2); // Мночитель в Числителе при расчёте 2-ой Производ. NURBS

    QVector<QPair<double, double>> c2(degree + 1);  // Индекс 2 для 2D задачи

    int j = 0; // кривая (нулевая производная)

    for(int i = 0; i < degree + 1; ++i)
    {
        //qDebug() << "----------";
        // qDebug() << "j =" << j << " i =" << i << " span - p + i =" << span - degree + i;
        //qDebug() << "nders[j][i] =" << nders[j][i] << " b[span - p + i] =" <<  control_points[span - degree + i] << " h[span - p + i] =" << w[span - degree + i];

        for(int k = 0; k < control_points[0].size(); ++k)
            n0[k] += control_points[span - degree + i][k] * w[span - degree + i] * nders[j][i];

        d += nders[j][i] * w[span - degree + i];
    }

    c2[0].first = n0[0] / d;
    c2[0].second = n0[1] / d;

    //qDebug() << "j = 0 - кривая \nc2 =" << c2;

    data_NURBS.curve = c2[0];
    data_NURBS.u = u_i;

    if(degree == 1)
        return;

    // 1-я производная

    for(int i = 0; i < degree + 1; ++i)
    {
        for(size_t k = 0; k < n1.size(); ++k)
        {
            n1[k] += control_points[span - degree + i][k] * w[span - degree + i] * nders[1][i];
            n2[k] += w[span - degree +i] * nders[1][i];
        }
    }

    c2[1].first = n1[0] / d - (n0[0] * n2[0]) / (d * d);
    c2[1].second = n1[1] / d - (n0[1] * n2[1]) / (d * d);

    //qDebug() << "j = 1 - 1-я производная \nc2 =" << c2;
    data_NURBS.derivative_1 = c2[1];

    // 2-я производная

    for(int i = 0; i < degree + 1; ++i)
    {
        for(size_t k = 0; k < n1.size(); ++k)
        {
            n3[k] += control_points[span - degree + i][k] * w[span - degree + i] * nders[2][i];
            n4[k] += w[span - degree +i] * nders[2][i];
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

    //qDebug() << "j = 2 - 2-я производная \nc2 =" << c2 << "\n";

    // Присваиваем координаты точки на кривой, 1-ю и 2-ю производную
    data_NURBS.derivative_2 = c2[2];
}

// Возвращает точку, повернутую на angle
QPair<double, double> rotate_point(const Point_curve& point, const double& angle = M_PI / 2)
{
    double rotatedX = point.derivative_1.first * cos(angle) - point.derivative_1.second * sin(angle);
    double rotatedY = point.derivative_1.first * sin(angle) + point.derivative_1.second * cos(angle);
    QPair<double, double> perpendicular {rotatedX + point.curve.first, rotatedY + point.curve.second};
    return perpendicular;
}

// Возвращает точки с длинной + epsilon
QPair<double, double> epsilon_point(const QPair<double, double>& rotated_point, const Point_curve& curve_point, const double& eps)
{
    double x = rotated_point.first - curve_point.curve.first;
    double y = rotated_point.second - curve_point.curve.second;

    double len = sqrt(x * x + y * y);
    x *= eps / len;
    y *= eps / len;

    return QPair<double, double> { x + curve_point.curve.first, y + curve_point.curve.second};
}

// Возвращает вектор с точками спанов реального узлового вектора
QVector<double> real_span_calc(const uint& degree, const uint& n_kn, const QVector<double>& u)
{
    // Реальный диапазон
    double u_start = u[degree];
    const double u_end = u[n_kn - degree - 1];
    QVector<double> u_real_span;

    for(int i = 1; u_start < u_end; ++i)
    {
        u_real_span.push_back(u_start);
        u_start = u[degree + i];
    }

    u_real_span.push_back(u_end);
    return u_real_span;
}

// Возвращает косинус между двумя точками
double cos_calc(const Point_curve& point_1, const QPair<double, double>& point_2)
{
    double x = point_1.curve.first - point_2.first;
    double y = point_1.curve.second - point_2.second;
    double numerator = abs(x * point_1.derivative_1.first + y * point_1.derivative_1.second);
    double denominator = radius_vector_len(point_1.derivative_1) * radius_vector_len(x, y);
    if(denominator == 0)
        return 0;
    else
        return numerator / denominator;
}

// Возвращает точку кривой, перпендикулярной точке на плоскости
Point_curve finding_perpendicular(const int& n_real, const uint& n_kn, const int& degree, const QVector<double>& u_vector,
                                  const QVector<QVector<double>>& polygon, const QVector<double>& w, const QPair<double, double>& point)
{
    QVector<Point_curve> point_u(n_real - 1); // Массив точек - перпендикуляров
    QVector<double> u_real_span = real_span_calc(degree, n_kn, u_vector); // Спаны реального диапазона узлового вектора
    QVector<double> cos_data(u_real_span.size() - 1);

    for(int i = 0; i < u_real_span.size() - 1; ++i) // Итерируемся по спанам, начиная с нулевого
    {
        point_u[i].u = (u_real_span[i + 1] - u_real_span[i]) / 2 + u_real_span[i]; // Берём среднее спана

        curve_point_and_deriv_NURBS(point_u[i], n_kn, degree, u_vector, polygon, w, point_u[i].u);

        const double EPSILON_ANGLE = 0.01; // Эпсилон для косинуса прямого угла

        do
        {
            double x = point_u[i].curve.first - point.first;
            double y = point_u[i].curve.second - point.second;
            double numerator = x * point_u[i].derivative_1.first + y * point_u[i].derivative_1.second;
            double denominator = x * point_u[i].derivative_2.first + y * point_u[i].derivative_2.second + pow(radius_vector_len(point_u[i].derivative_1), 2);

            if(denominator == 0)
                point_u[i].u = point_u[i].u;
            else
                point_u[i].u = point_u[i].u - numerator / denominator * 0.01; // Новая, приближённая точка кривой

            if(point_u[i].u < u_real_span[i]) // Если точка вышла из спана
            {
                point_u[i].u = u_real_span[i];
                curve_point_and_deriv_NURBS(point_u[i], n_kn, degree, u_vector, polygon, w, point_u[i].u);
                cos_data[i] = cos_calc(point_u[i], point);
                break;
            }
            else if(point_u[i].u > u_real_span[i + 1])
            {
                point_u[i].u = u_real_span[i + 1];
                curve_point_and_deriv_NURBS(point_u[i], n_kn, degree, u_vector, polygon, w, point_u[i].u);
                cos_data[i] = cos_calc(point_u[i], point);
                break;
            }

            curve_point_and_deriv_NURBS(point_u[i], n_kn, degree, u_vector, polygon, w, point_u[i].u);
            cos_data[i] = cos_calc(point_u[i], point);
        } while (cos_data[i] > EPSILON_ANGLE);
    }

    Point_curve point_min_len; // Точка кривой с минимальным расстоянием до точки на плоскости
    point_min_len = point_u[0]; // Присваиваем первую точку для дальнейшего сравнения
    double min_len = vector_len(point, point_u[0].curve);
    double min_cos = cos_data[0];

    for(int i = 1; i < point_u.size(); ++i) // Ищем точку с максимально нулевым косинусом
    {
        double temp_len = vector_len(point, point_u[i].curve);
        double temp_cos = cos_data[i];

        if(min_cos > temp_cos) // Сравниваем косинус двух точек
        {
            min_len = temp_len;
            min_cos = temp_cos;
            point_min_len = point_u[i];
        }
        else if(min_cos == temp_cos && min_len > temp_len) // Если косинусы равны - сравниваем длину векторов
        {
            min_len = temp_len;
            min_cos = temp_cos;
            point_min_len = point_u[i];
        }
    }

    /*
    plot_tangent(canvas, point_u[0]); // Рисуем касательную к точке
    plot_tangent(canvas, point_u[1]); // Рисуем касательную к точке
    plot_tangent(canvas, point_u[2]); // Рисуем касательную к точке

    plot_line(canvas, point.first, point.second, point_u[0].curve.first, point_u[0].curve.second, QColor(0, 0, 0), 2); // Рисуем перпендикуляр между точкой и кривой
    plot_line(canvas, point.first, point.second, point_u[1].curve.first, point_u[1].curve.second, QColor(0, 0, 0), 2); // Рисуем перпендикуляр между точкой и кривой
    plot_line(canvas, point.first, point.second, point_u[2].curve.first, point_u[2].curve.second, QColor(0, 0, 0), 2); // Рисуем перпендикуляр между точкой и кривой
    */

    return point_min_len;
}

// Возвращает кривую с пониженной степенью
QVector<Point_curve> declining_degree_curve(const QVector<QVector<double>>& control_points, const double& degree, const QVector<double> w, const int& number_u, const double& epsilon)
{
    // Рассчитываем точки исходной кривой
    const QVector<double> u = u_fill(control_points, degree); // Узловой вектор
    QVector<Point_curve> data_NURBS(number_u + 1);  // Содержит точки кривой, 1-ой и 2-ой производной

    const uint n_ver = control_points.size(); // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
    const uint n_kn = n_ver + degree + 1;   // Количество узлов (длина) в узловом векторе (n_knots)
    const uint n_real = n_ver - degree + 1; // Количество узлов (длина) реальной части узлового вектора

    for(int i = 0; i < number_u + 1; ++i)
    {
        double u_i = (i / static_cast<double>(number_u));
        curve_point_and_deriv_NURBS(data_NURBS[i], n_kn, degree, u, control_points, w, u_i);
    }

    uint degree_new = 1; // Начинаем с 1-ой степени кривой
    double max_perpendicular = 0;
    QVector<Point_curve> data_NURBS_new(number_u + 1);  // Содержит точки кривой, 1-ой и 2-ой производной

    do
    {
        const QVector<double> u_new = u_fill(control_points, degree_new); // Узловой вектор
        const uint n_ver_new = control_points.size(); // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
        const uint n_kn_new = n_ver_new + degree_new + 1;     // Количество узлов (длина) в узловом векторе (n_knots)

        for(int i = 0; i < number_u + 1; ++i)
        {
            double u_i = (i / static_cast<double>(number_u));
            curve_point_and_deriv_NURBS(data_NURBS_new[i], n_kn_new, degree_new, u_new, control_points, w, u_i);
        }

        Point_curve max_p1, max_p2;
        max_perpendicular = 0;
        for(int i = 0; i < data_NURBS_new.size() - 1; ++i)
        {
            Point_curve u_perpendicular = finding_perpendicular(n_real, n_kn, degree, u, control_points, w, data_NURBS_new[i].curve);
            double temp_max_perpendicular = vector_len(data_NURBS_new[i].curve, u_perpendicular);

            if(max_perpendicular < temp_max_perpendicular)
            {
                max_perpendicular = temp_max_perpendicular;
                max_p1 = u_perpendicular;
                max_p2 = data_NURBS_new[i];
            }
        }

        QVector<QPair<double, double>> rotated_points;
        QVector<Point_curve> epsilon_point_DATA (number_u + 1);

        QVector<QPair<double, double>> rotated_points1;
        QVector<Point_curve> reverse_epsilon_point_DATA (number_u + 1);

        for(int i = 0; i < data_NURBS.size(); ++i)
        {
            rotated_points.push_back(rotate_point(data_NURBS[i]));
            epsilon_point_DATA[i].curve = epsilon_point(rotated_points[i], data_NURBS[i], epsilon);

            rotated_points1.push_back(rotate_point(data_NURBS[i], - M_PI / 2));
            reverse_epsilon_point_DATA[i].curve = epsilon_point(rotated_points1[i], data_NURBS[i], epsilon);
        }

        if(epsilon > max_perpendicular)
            continue;
        else
            ++degree_new;

    }while(epsilon < max_perpendicular);

    qDebug() << "Степень новой кривой = " << degree_new;
    return data_NURBS_new;
}

#endif // FUNCTIONS_H
