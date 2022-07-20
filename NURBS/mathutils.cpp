#include "mathutils.h"
#include <cmath>

double MathUtils::radiusVectorLength(double x, double y)
{
    return sqrt(x * x + y * y);
}

double MathUtils::radiusVectorLength(const QPair<double, double> &point)
{
    return radiusVectorLength(point.first, point.second);
}

double MathUtils::vectorLenght(double x_1, double y_1, double x_2, double y_2)
{
    return sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2));
}

double MathUtils::vectorLenght(const QPair<double, double> &point_1, const QPair<double, double> &point_2)
{
    return vectorLenght(point_1.first, point_1.second, point_2.first, point_2.second);
}

double MathUtils::vectorLenght(const QPair<double, double> &point_1, const Point_curve &point_2)
{
    return vectorLenght(point_1.first, point_1.second, point_2.curve.first, point_2.curve.second);
}

double MathUtils::angleBetweenVectors(const QPair<double, double> &vec_1_start, const QPair<double, double> &vec_1_end, const QPair<double, double> &vec_2_start, const QPair<double, double> &vec_2_end)
{
    double vec_x_1 = vec_1_end.first - vec_1_start.first;
    double vec_y_1 = vec_1_end.second - vec_1_start.second;
    double vec_x_2 = vec_2_end.first - vec_2_start.first;
    double vec_y_2 = vec_2_end.second - vec_2_start.second;
    double numerator = vec_x_1 * vec_x_2 + vec_y_1 * vec_y_2;
    double len_vec_1 = MathUtils::vectorLenght(vec_1_start, vec_1_end);
    double len_vec_2 = MathUtils::vectorLenght(vec_2_start, vec_2_end);
    double cos = numerator / (len_vec_1 * len_vec_2);
    return acos(cos) * 180 / M_PI;
}
