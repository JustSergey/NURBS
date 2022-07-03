#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "widget.h"
#include <QDebug>

struct curve // Хранит точку кривой, её производную (1 и 2), интервал (span) и точку u
{
    QPair<double, double> curve;
    QPair<double, double> derivative_1;
    QPair<double, double> derivative_2;
    int span;
    double u;
};

// Определяет индекс узлового промежутка (интервал)
uint findSpan(QVector<curve>& data_CurvePoin_and_Deriv_NURBS, const uint& n, const int& p, const std::vector<double>& u, const double& u_i)
/*
 * Вход: n, p, u, u_i
 * Выход: индекс  узлового  промежутка
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m=n+1+p
 * u_i - точка внутри реального диатазона в узловом векторе
*/
{
    static uint counter; // Для индексирования нужной точки нужного span

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

    data_CurvePoin_and_Deriv_NURBS[counter++].span = middle;

    return middle;
}

#endif // FUNCTIONS_H