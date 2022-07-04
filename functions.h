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
 * n - кол-во Control Points (счёт от нуля)
 * p - степень полинома(=degree)
 * u - узловой вектор - мах индекс в нем m = n + 1 + p
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

#endif // FUNCTIONS_H
