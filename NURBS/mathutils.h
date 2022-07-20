#ifndef MATHUTILS_H
#define MATHUTILS_H

#include "nurbs.h"

class MathUtils
{
public:
    // Возвращают длину радиус-вектора
    static double radiusVectorLength(double x, double y);
    static double radiusVectorLength(const QPair<double, double> &point);

    // Возвращают длину вектора
    static double vectorLenght(double x_1, double y_1, double x_2, double y_2);
    static double vectorLenght(const QPair<double, double> &point_1, const QPair<double, double> &point_2);
    static double vectorLenght(const QPair<double, double> &point_1, const Point_curve &point_2);

    // Возвращает угол между двумя векторами
    static double angleBetweenVectors(const QPair<double, double> &vec_1_start, const QPair<double, double> &vec_1_end, const QPair<double, double> &vec_2_start, const QPair<double, double> &vec_2_end);
};

#endif // MATHUTILS_H
