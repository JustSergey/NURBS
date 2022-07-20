#ifndef NURBS_H
#define NURBS_H

#include <math.h>
#include <QDebug>

// Хранит точку кривой, её 1-ую и 2-ую производную, интервал (span) и параметр
struct Point_curve
{
    QPair<double, double> curve;    // Точка кривой
    QPair<double, double> deriv_1;  // Точка 1-ой производной
    QPair<double, double> deriv_2;  // Точка 2-ой производной
    double real_point;              // Точка реальной части узл. вектора "u"
    uint span;                      // Интервал
};

class NURBS
{
public:
    explicit NURBS();
    explicit NURBS(const QVector<QVector<double>> &control_points, const QVector<double> &weight, const uint &degree, const uint &num_point);

    void set_control_point(const QVector<QVector<double>> &control_points);
    void set_weight(const QVector<double> &weight);
    void set_degree(const uint &degree);
    void set_num_splits(const uint &num_point);

    QVector<Point_curve> get_point_curve() const;
    uint get_num_splits() const;

    void fill_nodal_vector();        // Равномерно заполняет узловой вектор
    void calc_points_derivs_NURBS(); // Рассчитывает все точки кривой и их 1-ую и 2-ую производную
    Point_curve find_perpendicular_curve(const QPair<double, double>& point); // Возвращает точку кривой, перпендикулярной точке на плоскости (point)
    NURBS decrease_degree_curve(const double &epsilon); // Возвращает кривую с пониженной степенью
    QPair<QVector<Point_curve>, QVector<Point_curve>> calc_curves_shift(const double &len) const; // Возвращает пару векторов с точками кривых, отдалённых от исходной кривой на длину (len)

private:
    QVector<Point_curve> m_points_NURBS;       // Точки NURBS
    QVector<QVector<double>> m_control_points; // Точки определяющего многоугольника
    QVector<double> m_nodal_vector;            // Узловой вектор
    double m_real_kn_start;   // Начало реального диапазона узл. вектора
    double m_real_kn_stop;    // Конец реального диапазона узл. вектора
    QVector<double> m_weight; // Весовые коэффициенты точек опр. многоуг.
    uint m_degree;            // Степень аппроксимирующих полиномов
    uint m_num_splits;        // Кол-во разбиений (точек) в реальной части узл. вектора "u"
    uint m_num_ver;           // Кол-во вершин в определяющем многоугольнике (num_vertices)
    uint m_num_kn;            // Кол-во узлов (длина) в узловом векторе (num_knots)
    uint m_num_real_kn;       // Кол-во узлов (длина) реальной части узлового вектора

    uint find_span(const double &u_i) const; // Возвращает индекс узлового промежутка для заданной точки
    void calc_derivs_basis_func(QVector<QVector<double>> &derivs_basis, const double &real_point, const double &span); // Вычисляет ненулевые базисные функции и их производные
    void calc_point_derivs_NURBS(Point_curve &point_NURBS, const double &real_point); // Рассчитывает точку кривой и её 1-ую и 2-ую производную
    QVector<double> calc_points_real_span() const; // Возвращает вектор с параметрами спанов реального узлового вектора
    QPair<double, double> rotate_point(const Point_curve &point, const double &angle) const; // Возвращает повёрнутую точку на определённый угол (angle)
    QPair<double, double> point_shift(const Point_curve &point_1, const QPair<double, double> &point_2, const double &len) const; // Возвращает точку вектора со сдвигом в длину (len)
    double cos_between_vectors(const Point_curve &point_curve, const QPair<double, double> &point_end_vec) const; // Возвращает косинус между производной в точке кривой и вектором, начинающимся от точки кривой
    class Impl;
    Impl* _impl;
};

#endif // NURBS_H
