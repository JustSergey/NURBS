#include "nurbs.h"
#include "mathutils.h"
#include <cmath>

class NURBS::Impl
{

};

NURBS::NURBS()
{
}

NURBS::NURBS(const QVector<QVector<double>> &control_points, const QVector<double> &weight, const uint &degree, const uint &num_splits)
    : m_control_points { control_points }, m_weight { weight }, m_degree { degree }, m_num_splits { num_splits }
{
    m_points_NURBS.resize(m_num_splits + 1);

    m_num_ver = control_points.size();
    m_num_kn = m_num_ver + degree + 1;
    m_num_real_kn = m_num_ver - degree + 1;

    /* Рассчитываем реальный диапазон */
    m_real_kn_start = m_degree;
    m_real_kn_stop = m_num_kn - m_degree - 1;
};

void NURBS::set_control_point(const QVector<QVector<double>> &control_points)
{
    m_control_points = control_points;
}

void NURBS::set_weight(const QVector<double> &weight)
{
    m_weight = weight;
}

void NURBS::set_degree(const uint &degree)
{
    m_degree = degree;
}

void NURBS::set_num_splits(const uint &num_splits)
{
    m_num_splits = num_splits;
}

QVector<Point_curve> NURBS::get_point_curve() const
{
    return m_points_NURBS;
}

uint NURBS::get_num_splits() const
{
    return m_num_splits;
}

// Равномерно заполняет узловой вектор
void NURBS::fill_nodal_vector()
{
    m_nodal_vector.resize(m_num_kn);

    const double step = 1 / static_cast<double>(m_num_real_kn - 1); // Шаг в реальном диапазоне

    for(uint i = m_real_kn_start + 1; i < m_real_kn_stop; ++i)   // Заполняем реальный диапазон
        m_nodal_vector[i] = m_nodal_vector[i - 1] + step;

    for(int i = m_real_kn_stop; i < m_nodal_vector.size(); ++i)  // Заполняем последние параметры единицами
        m_nodal_vector[i] = 1;
}

// Возвращает индекс узлового промежутка для заданной точки
uint NURBS::find_span(const double &real_point) const
{
    for(int i = 0; i < m_nodal_vector.size() - 1; ++i)
    {
        if(m_nodal_vector[i] > m_nodal_vector[i + 1])
            qDebug() << "Error! find_span: узловой вектор убывает m_nodal_vector[i] > m_nodal_vector[i + 1]";
    }

    if(static_cast<uint>(m_nodal_vector.size()) != (m_num_ver + m_degree + 1))
        qDebug() << "Error! find_span: m_nodal_vector.size()) != (m_num_ver + m_degree + 1)";

    if(real_point < m_nodal_vector[m_degree] || real_point > m_nodal_vector[m_num_kn - m_degree - 1])
        qDebug() << "Error! find_span: m_nodal_vector вышел за реальный диапазон";

    if(real_point == m_nodal_vector[m_num_kn - m_degree - 1]) // Если дошли до конца реального диапазона
        return m_num_kn - m_degree - 2;

    uint low = m_degree, high = m_num_kn - m_degree - 1, middle = (low + high) / 2;

    /* Выполняем двоичный поиск */
    while((real_point < m_nodal_vector[middle]) || (real_point >= m_nodal_vector[middle + 1]))
    {
        if(real_point < m_nodal_vector[middle])
            high = middle;
        else
            low = middle;

        middle = (low + high) / 2;
    }

    return middle;
}

// Вычисляет ненулевые базисные функции и их производные
void NURBS::calc_derivs_basis_func(QVector<QVector<double>> &derivs_basis, const double &real_point, const double &span)
{
    QVector<double> left(m_degree + 1), right(m_degree + 1);
    QVector<QVector<double>> ndu(m_degree + 1, QVector<double>(m_degree + 1)); // Для хранения базисных функций и узлов различия
    ndu[0][0] = 1;

    for(uint i = 1; i < m_degree + 1; ++i)
    {
        left[i] = real_point - m_nodal_vector[span + 1 - i];
        right[i] = m_nodal_vector[span + i] - real_point;
        double saved = 0;

        for(uint j = 0; j < i; ++j)
        {
            // Нижний треугольник
            ndu[i][j] = right[j + 1] + left[i - j];
            double temp = ndu[j][i - 1] / ndu[i][j];
            // Верхний треугольник
            ndu[j][i] = saved + right[j + 1] * temp;
            saved = left[i - j] * temp;
        }

        ndu[i][i] = saved;
    }

    for(uint i = 0; i <= m_degree; ++i) // Загружаем базисные функции
        derivs_basis[0][i] = ndu[i][m_degree];

    /* Вычисляем производные */

    QVector<QVector<double>> a(2, QVector<double>(m_degree + 1)); // Хранит два вычесленных ряда

    for(int i = 0; i < static_cast<int>(m_degree + 1); ++i) // Цикл по индексу функции
    {
        int s1 = 0, s2 = 1; // Альтернативные строки в массиве
        a[0][0] = 1.0;

        for(int k = 1; k <= static_cast<int>(m_degree); ++k)
        {
            double d = 0;
            double rk = i - k;
            double pk = m_degree - k;

            if(i >= k)
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

            if(i - 1 <= pk)
                j2 = k - 1;
            else
                j2 = m_degree - i;

            for(uint j = j1; j <= j2; ++j)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if(i <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][i];
                d += a[s2][k] * ndu[i][pk];
            }

            derivs_basis[k][i] = d;

            /* Меняем строки местами */
            double temp = s1;
            s1 = s2;
            s2 = temp;
        }
    }

    /* Умножаем на правильные коэффициенты */

    double r = m_degree;

    for(uint i = 1; i <= m_degree; ++i)
    {
        for(uint j = 0; j < m_degree + 1; ++j)
            derivs_basis[i][j] *= r;

        r *= m_degree - i;
    }

    /* Для контроля суммируем значения базисных Функций в точке */
    double sum = 0;

    for(int i = 0; i < derivs_basis.size(); ++i)
        sum += derivs_basis[0][i];

    if((sum < (1 - 1e-10)) || (sum > 1 + 1e-10)) // Если все верно, то сумма должна = 1
        qDebug() << "Error! calc_derivs_basis_func: Сумма базисных Функций != 1";
}

// Рассчитывает точку кривой и её 1-ую и 2-ую производную
void NURBS::calc_point_derivs_NURBS(Point_curve &point_NURBS, const double &real_point)
{
    double span = find_span(real_point);
    point_NURBS.span = span;

    // Содержит для заданного узлового вектора ненулевые базис. функции и их производные
    QVector<QVector<double>> derivs_basis(m_degree + 1, QVector<double>(m_degree + 1));

    calc_derivs_basis_func(derivs_basis, real_point, span);

    /* Рассчитываем точку кривой */

    double denominator = 0;
    QVector<double> n0(2);

    for(uint i = 0; i < m_degree + 1; ++i)
    {
        for(int k = 0; k < m_control_points[0].size(); ++k)
            n0[k] += m_control_points[span - m_degree + i][k] * m_weight[span - m_degree + i] * derivs_basis[0][i];

        denominator += derivs_basis[0][i] * m_weight[span - m_degree + i];
    }

    QVector<QPair<double, double>> point(m_degree + 1); // Содержит точку кривой, её 1-ую и 2-ую производную
    point[0].first = n0[0] / denominator;
    point[0].second = n0[1] / denominator;

    point_NURBS.curve = point[0];
    point_NURBS.real_point = real_point;

    /* Рассчитываем 1-ую производную */

    QVector<double> n1(2);
    QVector<double> n2(2);

    for(uint i = 0; i < m_degree + 1; ++i)
    {
        for(int k = 0; k < n1.size(); ++k)
        {
            n1[k] += m_control_points[span - m_degree + i][k] * m_weight[span - m_degree + i] * derivs_basis[1][i];
            n2[k] += m_weight[span - m_degree +i] * derivs_basis[1][i];
        }
    }

    point[1].first = n1[0] / denominator - (n0[0] * n2[0]) / (denominator * denominator);
    point[1].second = n1[1] / denominator - (n0[1] * n2[1]) / (denominator * denominator);

    point_NURBS.deriv_1 = point[1];

    if(m_degree == 1) // Если степень полинома = 1, то заканчиваем считать
        return;

    /* Рассчитываем 2-ую производную */

    QVector<double> n3(2);
    QVector<double> n4(2);

    for(uint i = 0; i < m_degree + 1; ++i)
    {
        for(int k = 0; k < n1.size(); ++k)
        {
            n3[k] += m_control_points[span - m_degree + i][k] * m_weight[span - m_degree + i] * derivs_basis[2][i];
            n4[k] += m_weight[span - m_degree +i] * derivs_basis[2][i];
        }
    }

    QVector<double> s1(2);

    for(int i = 0; i < s1.size(); ++i)
         s1[i] = n3[i] / denominator - (n1[i] * n2[i]) / (denominator * denominator);

    QVector<double> nn(2);

    for(int i = 0; i < nn.size(); ++i)
        nn[i] = n1[i] * n2[i];

    QVector<double> nn_deriv(2);

    for(int i = 0; i < nn_deriv.size(); ++i)
        nn_deriv[i] = n1[i] * n2[i] + n0[i] * n4[i];

    QVector<double> s2(2);

    for(int i = 0; i < s2.size(); ++i)
        s2[i] = nn_deriv[i] / (denominator * denominator) - (nn[i] * 2 * n2[i]) / (pow(denominator, 4));

    point[2].first = s1[0] - s2[0];
    point[2].second = s1[1] - s2[1];

    // Присваиваем координаты точки кривой, 1-ю и 2-ю производную
    point_NURBS.deriv_2 = point[2];
}

// Рассчитывает все точки кривой и их 1-ую и 2-ую производную
void NURBS::calc_points_derivs_NURBS()
{
    for(uint i = 0; i < m_num_splits + 1; ++i)
    {
        double real_point = (i / static_cast<double>(m_num_splits)); // Точка реальной части узлового вектора
        calc_point_derivs_NURBS(m_points_NURBS[i], real_point);
    }
}

// Возвращает повёрнутую точку на определённый угол (angle)
QPair<double, double> NURBS::rotate_point(const Point_curve &point, const double &angle) const
{
    double rotated_X = point.deriv_1.first * cos(angle) - point.deriv_1.second * sin(angle);
    double rotated_Y = point.deriv_1.first * sin(angle) + point.deriv_1.second * cos(angle);
    return QPair<double, double> {rotated_X + point.curve.first, rotated_Y + point.curve.second};
}

// Возвращает точку вектора со сдвигом в длину (len)
QPair<double, double> NURBS::point_shift(const Point_curve &point_1, const QPair<double, double> &point_2, const double &len) const
{
    double x = point_2.first - point_1.curve.first;
    double y = point_2.second - point_1.curve.second;
    double len_vec = MathUtils::radiusVectorLength(x, y);
    x *= len / len_vec;
    y *= len / len_vec;
    return QPair<double, double> {x + point_1.curve.first, y + point_1.curve.second};
}

// Возвращает вектор с параметрами спанов реального узлового вектора
QVector<double> NURBS::calc_points_real_span() const
{
    double p_start = m_nodal_vector[m_degree];
    const double p_end = m_nodal_vector[m_num_kn - m_degree - 1];
    QVector<double> points_real_span;

    for(int i = 1; p_start < p_end; ++i)
    {
        points_real_span.push_back(p_start);
        p_start = m_nodal_vector[m_degree + i];
    }

    points_real_span.push_back(p_end);
    return points_real_span;
}

// Возвращает косинус между производной в точке кривой и вектором, начинающимся от точки кривой
double NURBS::cos_between_vectors(const Point_curve &point_curve, const QPair<double, double> &point_end_vec) const
{
    double x = point_curve.curve.first - point_end_vec.first;
    double y = point_curve.curve.second - point_end_vec.second;
    double numerator = abs(x * point_curve.deriv_1.first + y * point_curve.deriv_1.second);
    double denominator = MathUtils::radiusVectorLength(point_curve.deriv_1) * MathUtils::radiusVectorLength(x, y);

    if(denominator == 0)
        return 0;
    else
        return numerator / denominator;
}

// Возвращает точку кривой, перпендикулярной точке на плоскости (point)
Point_curve NURBS::find_perpendicular_curve(const QPair<double, double> &point)
{
    QVector<Point_curve> perpend_points(m_num_real_kn - 1); // Массив точек - перпендикуляров
    QVector<double> real_spans = calc_points_real_span();   // Спаны реального диапазона узлового вектора
    QVector<double> cos_perpendicular(real_spans.size() - 1);

    for(int i = 0; i < real_spans.size() - 1; ++i) // Итерируемся по спанам, начиная с нулевого
    {
        perpend_points[i].real_point = (real_spans[i + 1] - real_spans[i]) / 2 + real_spans[i]; // Берём среднее спана

        calc_point_derivs_NURBS(perpend_points[i], perpend_points[i].real_point);

        const double EPSILON_ANGLE = 0.01; // Эпсилон для косинуса прямого угла

        do
        {
            double x = perpend_points[i].curve.first - point.first;
            double y = perpend_points[i].curve.second - point.second;
            double numerator = x * perpend_points[i].deriv_1.first + y * perpend_points[i].deriv_1.second;
            double denominator = x * perpend_points[i].deriv_2.first + y * perpend_points[i].deriv_2.second + pow(MathUtils::radiusVectorLength(perpend_points[i].deriv_1), 2);

            if(denominator == 0)
                perpend_points[i].real_point = perpend_points[i].real_point;
            else
                perpend_points[i].real_point = perpend_points[i].real_point - numerator / denominator * 0.01; // Новая, приближённая точка кривой

            if(perpend_points[i].real_point < real_spans[i]) // Если точка вышла из спана
            {
                perpend_points[i].real_point = real_spans[i];
                calc_point_derivs_NURBS(perpend_points[i], perpend_points[i].real_point);
                cos_perpendicular[i] = cos_between_vectors(perpend_points[i], point);
                break;
            }
            else if(perpend_points[i].real_point > real_spans[i + 1])
            {
                perpend_points[i].real_point = real_spans[i + 1];
                calc_point_derivs_NURBS(perpend_points[i], perpend_points[i].real_point);
                cos_perpendicular[i] = cos_between_vectors(perpend_points[i], point);
                break;
            }

            calc_point_derivs_NURBS(perpend_points[i], perpend_points[i].real_point);
            cos_perpendicular[i] = cos_between_vectors(perpend_points[i], point);
        } while (cos_perpendicular[i] > EPSILON_ANGLE);
    }

    Point_curve point_min_len; // Точка кривой с минимальным расстоянием до точки на плоскости
    point_min_len = perpend_points[0]; // Присваиваем первую точку для дальнейшего сравнения
    double min_len = MathUtils::vectorLenght(point, perpend_points[0].curve);
    double min_cos = cos_perpendicular[0];

    for(int i = 1; i < perpend_points.size(); ++i)
    {
        double temp_len = MathUtils::vectorLenght(point, perpend_points[i].curve);
        double temp_cos = cos_perpendicular[i];

        if(min_cos > temp_cos) // Сравниваем косинус двух точек
        {
            min_len = temp_len;
            min_cos = temp_cos;
            point_min_len = perpend_points[i];
        }
        else if(min_cos == temp_cos && min_len > temp_len) // Если косинусы равны - сравниваем длину векторов
        {
            min_len = temp_len;
            min_cos = temp_cos;
            point_min_len = perpend_points[i];
        }
    }

    return point_min_len;
}

// Возвращает кривую с пониженной степенью
NURBS NURBS::decrease_degree_curve(const double &epsilon)
{
    uint degree_new = 1; // Начинаем с 1-ой степени кривой
    NURBS curve_new;

    do
    {
        NURBS curve_new(m_control_points, m_weight, degree_new, m_num_splits);
        curve_new.fill_nodal_vector();
        curve_new.calc_points_derivs_NURBS();

        Point_curve max_p1, max_p2;
        double max_perpendicular = 0;

        for(int i = 0; i < curve_new.m_points_NURBS.size() - 1; ++i)
        {
            Point_curve u_perpendicular = find_perpendicular_curve(curve_new.m_points_NURBS[i].curve);
            double temp_max_perpendicular = MathUtils::vectorLenght(curve_new.m_points_NURBS[i].curve, u_perpendicular);

            if(max_perpendicular < temp_max_perpendicular)
            {
                max_perpendicular = temp_max_perpendicular;
                max_p1 = u_perpendicular;
                max_p2 = curve_new.m_points_NURBS[i];
            }
        }

        if(epsilon > max_perpendicular)
            break;

        ++degree_new;
    }while(this->m_degree > degree_new);

    return curve_new;
}

// Возвращает пару векторов с точками кривых, отдалённых от исходной кривой на длину (len)
QPair<QVector<Point_curve>, QVector<Point_curve> > NURBS::calc_curves_shift(const double &len) const
{
    QVector<QPair<double, double>> rotated_points;
    QVector<Point_curve> epsilon_point_DATA (m_num_splits + 1);

    QVector<QPair<double, double>> rotated_points1;
    QVector<Point_curve> reverse_epsilon_point_DATA (m_num_splits + 1);

    for(int i = 0; i < m_points_NURBS.size(); ++i)
    {
        rotated_points.push_back(rotate_point(m_points_NURBS[i], M_PI / 2));
        epsilon_point_DATA[i].curve = point_shift(m_points_NURBS[i], rotated_points[i], len);

        rotated_points1.push_back(rotate_point(m_points_NURBS[i], - M_PI / 2));
        reverse_epsilon_point_DATA[i].curve = point_shift(m_points_NURBS[i], rotated_points1[i], len);
    }

    return QPair<QVector<Point_curve>, QVector<Point_curve>> (epsilon_point_DATA, reverse_epsilon_point_DATA);
}
