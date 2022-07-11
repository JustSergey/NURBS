#include "widget.h"
#include "ui_widget.h"
#include "functions.h"
#include "charts.h"
#include <QDebug>

Widget::Widget(QWidget *parent)
    : QWidget(parent), ui(new Ui::Widget)
{
    ui->setupUi(this);

/*
    const QVector<QVector<double>> control_points_1 // Точки определяющего многоугольника
    {
        {1.25, 1.3},
        {2.5, 3.9},
        {5.6, 3.9},
        {6.25, 1.3},
        {7.5, 2.6}
    };

    const QVector<double> w_1 {1, 1, 1, 1, 1}; // Весовые коэффициенты
    const uint degree_1 = 2;                   // Степень аппроксимирующих полиномов
    const int number_u_1 = 60;                 // Кол-во разбиений (точек) в реальной части узлов. вектора

    const QVector<double> u_1 = u_fill(control_points_1, degree_1); // Узловой вектор

    QVector<QVector<double>> nders_1(degree_1 + 1, QVector<double>(degree_1 + 1)); // Содержит для заданного "u" массив BASIS функций и 1-ую и 2-ую производную
    QVector<QPair<double, double>> c2_1(degree_1 + 1);  // Индекс 2 для 2D задачи
    QVector<Point_curve> data_NURBS_1(number_u_1 + 1);    // Содержит точки кривой, 1-ой и 2-ой производной

    const uint n_ver_1 = control_points_1.size(); // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
    const uint n_real_1 = n_ver_1 - degree_1 + 1;   // Количество узлов (длина) реальной части узлового вектора
    const uint n_kn_1 = n_ver_1 + degree_1 + 1;     // Количество узлов (длина) в узловом векторе (n_knots)

    for(int i = 0; i < number_u_1 + 1; ++i)
    {
        double u_i = (i / static_cast<double>(number_u_1));
        curve_point_and_deriv_NURBS(data_NURBS_1[i], n_real_1, n_kn_1, degree_1, u_1, control_points_1, w_1, u_i, c2_1, nders_1);
    }



    QString title = "Метрики разности между двух кривых";
    QString labels_legend_1 = "Определяющий многоуг.";
    QString labels_legend_2 = "Исходный NURBS";

    int x_min = -10, x_max = 10;
    int y_min = -10, y_max = 10;
    curve_plot(ui->graph_function, control_points_1, data_NURBS_1, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);

    title = "1-я прoизвдная NURBS 2-ой степени";
    x_min = -15, x_max = 15;
    y_min = -15, y_max = 15;
    first_derivative_plot(ui->graph_first_derivative, control_points_1, data_NURBS_1, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);

    title = "2-я прoизвдная NURBS 2-ой степени";
    x_min = -100, x_max = 100;
    y_min = -100, y_max = 100;
    second_derivative_plot(ui->graph_second_derivative, control_points_1, data_NURBS_1, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);




    QVector<double> real_spans = real_span_calc(degree_1, n_kn_1, u_1); // Вектор с точками реального диапазона спанов
    QVector<Point_curve> u_real_spans(real_spans.size());

    for(int i = 0; i < u_real_spans.size(); ++i)
    {
        curve_point_and_deriv_NURBS(u_real_spans[i], n_real_1, n_kn_1, degree_1, u_1, control_points_1, w_1, real_spans[i], c2_1, nders_1);
        plot_point(ui->graph_function, u_real_spans[i].curve.first, u_real_spans[i].curve.second, 9, "", QColor(147, 112, 219)); // Рисуем точки на графике (точки границ реального диапазона спанов)
    }

    // *ДОБАВЛЯЕТ В ЛЕГЕНДУ* //
    ui->graph_function->addGraph();
    ui->graph_function->graph()->setPen(QColor(147, 112, 219));
    ui->graph_function->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 9)); // Формируем вид точек
    ui->graph_function->graph()->setLineStyle(QCPGraph::lsNone);
    ui->graph_function->graph()->addData(u_real_spans[0].curve.first, u_real_spans[0].curve.second);
    ui->graph_function->graph()->setName("Точка смены"); // Точка смены полинома/элемента узлового вектора

    ui->graph_function->replot();



    QVector<double> point_u { 0, 0.2, 2/5.0, 0.6, 0.8, 1 }; // Массив, хранящий точки u, от которых пойдёт производная на графике
    QVector<Point_curve> derivs_curve(point_u.size());

    for(int i = 0; i < point_u.size(); ++i) // Считаем координаты и производные для точек u
    {
        curve_point_and_deriv_NURBS(derivs_curve[i], n_real_1, n_kn_1, degree_1, u_1, control_points_1, w_1, point_u[i], c2_1, nders_1);
    }

    for(const auto& el: derivs_curve)
    {
        derivative_point_line(ui->graph_function, el); // Рисуем линию производной в точке (касательную)
    }


    plot_lable(ui->graph_function, derivs_curve[0].derivative_1.first + derivs_curve[0].curve.first + 1, derivs_curve[0].derivative_1.second + derivs_curve[0].curve.second + 0.2, "C'(0)");
    plot_lable(ui->graph_function,derivs_curve[1].derivative_1.first + derivs_curve[1].curve.first + 1.4, derivs_curve[1].derivative_1.second + derivs_curve[1].curve.second + 0.4, "C'(1/5)");
    plot_lable(ui->graph_function,derivs_curve[2].derivative_1.first + derivs_curve[2].curve.first + 1, derivs_curve[2].derivative_1.second + derivs_curve[2].curve.second - 0.5, "C'(2/5)");
    plot_lable(ui->graph_function,derivs_curve[3].derivative_1.first + derivs_curve[3].curve.first + 1, derivs_curve[3].derivative_1.second + derivs_curve[3].curve.second - 0.5, "C'(3/5)");
    plot_lable(ui->graph_function,derivs_curve[4].derivative_1.first + derivs_curve[4].curve.first + 1.3, derivs_curve[4].derivative_1.second + derivs_curve[4].curve.second - 0.5, "C'(4/5)");
    plot_lable(ui->graph_function,derivs_curve[5].derivative_1.first + derivs_curve[5].curve.first + 1, derivs_curve[5].derivative_1.second + derivs_curve[5].curve.second + 0.1, "C'(1)");
    //derivative_point_line(ui->graph_first_derivative, derivs_curve); // Рисуем линию производной в точке (касательную)

    plot_lable(ui->graph_first_derivative, derivs_curve[0].derivative_1.first + 1, derivs_curve[0].derivative_1.second + 0.5, "C'(0)");
    plot_lable(ui->graph_first_derivative, derivs_curve[1].derivative_1.first + 1.5, derivs_curve[1].derivative_1.second + 0.3, "C'(1/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[2].derivative_1.first + 1.4, derivs_curve[2].derivative_1.second - 0.5, "C'(2/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[3].derivative_1.first + 0.5, derivs_curve[3].derivative_1.second - 0.8, "C'(3/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[4].derivative_1.first + 1.4, derivs_curve[4].derivative_1.second - 0.5, "C'(4/5)");
    plot_lable(ui->graph_first_derivative, derivs_curve[5].derivative_1.first + 0.5, derivs_curve[5].derivative_1.second + 0.8, "C'(1)");

    for(const auto& p: derivs_curve)
    {
        plot_point(ui->graph_function, p.curve.first, p.curve.second);
    }



    /*ц
    QPair<double, double> point(2, 6); // Точка, к которой мы будем проводить перпендикуляр
    Point_curve u_perpendicular = finding_perpendicular(n, p, u, b, h, point, ui->graph_function); // Ближайшая точка для перпендикуляра

    //plot_line(ui->graph_function, point.first, point.second, u_perpendicular.curve.first, u_perpendicular.curve.second, QColor(178, 34, 34)); // Рисуем перпендикуляр между точкой и кривой
    //plot_tangent(ui->graph_function, u_perpendicular); // Рисуем касательную к точке

    // *ДОБАВЛЯЕТ В ЛЕГЕНДУ* //
    ui->graph_function->addGraph();
    ui->graph_function->graph()->setPen(QColor(244, 164, 96));
    ui->graph_function->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 11)); // Формируем вид точек
    ui->graph_function->graph()->setLineStyle(QCPGraph::lsNone);
    ui->graph_function->graph()->addData(point.first, point.second);
    ui->graph_function->graph()->setName("Заданная точка"); // Точка смены полинома/элемента узлового вектора



    QCPCurve *curve = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    QPen pen;
    pen.setColor(QColor(0, 100, 0));
    pen.setWidthF(2.8);
    curve->setPen(pen);
    curve->setName("Касательная"); // Обзываем полигон в легенде графика
    */

    /*
    QCPCurve *curve1 = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    QPen pen1;
    pen1.setWidthF(2);
    pen1.setColor(QColor(0, 0, 0));
    pen1.setStyle(Qt::PenStyle::DashLine);
    curve1->setPen(pen1);
    curve1->setName("Перпендикуляр"); // Обзываем полигон в легенде графика

    plot_lable_with_arrow(ui->graph_function, 6.35, 6.35, 4.8, 4.27, "Предельные\nслучаи");
    plot_lable_with_arrow(ui->graph_function, 6.35, 6.35, 4, 4.36, "");
    plot_lable_with_arrow(ui->graph_function, 2, 3.5, 2.75, 4.5, "");
    plot_lable(ui->graph_function, 0.8, 3.1, "Перпендикуляр\nс мин. длиной");
    */


/*
    const QVector<QVector<double>> control_points_2 // Массив точек многоугольника
    {
        {0.8, 2},
        {4, 3},
        {7.6, 2}
    };

    const QVector<double> w_2 {1, 1, 1}; // Весовые коэффициенты
    const uint degree_2 = 2;          // Степень аппроксимирующих полиномов
    const int number_u_2 = 60;        // Кол-во разбиений (точек) в реальной части узлов. вектора

    const QVector<double> u_2 = u_fill(control_points_2, degree_2); // Узловой вектор

    QVector<QVector<double>> nders_2(degree_2 + 1, QVector<double>(degree_2 + 1)); // Содержит для заданного "u" массив BASIS функций и 1-ую и 2-ую производную
    QVector<QPair<double, double>> c2_2(degree_2 + 1);  // Индекс 2 для 2D задачи
    QVector<Point_curve> data_NURBS_2(number_u_2 + 1);      // Содержит точки кривой, 1-ой и 2-ой производной

    const uint n_ver_2 = control_points_2.size();  // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
    const uint n_real_2 = n_ver_2 - degree_2 + 1;  // Количество узлов (длина) реальной части узлового вектора

    for(int i = 0; i < number_u_2 + 1; ++i)
    {
        double u_i = (i / static_cast<double>(number_u_2));
        curve_point_and_deriv_NURBS(data_NURBS_2[i], n_real_2, degree_2, u_2, control_points_2, w_2, u_i, c2_2, nders_2);
    }

    plot_curve(ui->graph_function, data_NURBS_2, "Аппроксимирующий NURBS", QColor(34, 139, 34)); // Рисуем сплайн

/*
    for(int i = 0; i < data_NURBS_2.size();)
    {
        Point_curve u_perpendicular = finding_perpendicular(n_real_1, degree_1, u_1, control_points_1, w_1, data_NURBS_2[i].curve);
        plot_line(ui->graph_function, data_NURBS_2[i].curve.first, data_NURBS_2[i].curve.second, u_perpendicular.curve.first, u_perpendicular.curve.second, QColor(178, 34, 34)); // Рисуем перпендикуляр между точкой и кривой
        i += 2;
    }
*/
    /*
    Point_curve u_perpendicular;
    for(int i = 0; i < data_NURBS_1.size();)
    {
        u_perpendicular = finding_perpendicular(n_real_2, degree_2, u_2, control_points_2, w_2, data_NURBS_1[i].curve);
        plot_line(ui->graph_function, data_NURBS_1[i].curve.first, data_NURBS_1[i].curve.second, u_perpendicular.curve.first, u_perpendicular.curve.second); // Рисуем перпендикуляр между точкой и кривой
        i += 3;
    }

    QCPCurve *curve = new QCPCurve(ui->graph_function->xAxis, ui->graph_function->yAxis);
    QPen pen;
    pen.setColor(QColor(0, 0, 0));
    pen.setStyle(Qt::PenStyle::DashLine);
    pen.setWidthF(2.5);
    curve->setPen(pen);
    curve->setName("Метрика (перпендикуляр)"); // Обзываем полигон в легенде графика
    */


    ///---------------////


   /* const QVector<QVector<double>> control_points_1 // Точки определяющего многоугольника
    {
        {1.2, 1},
        {4.5, 2,5},
        {7, 6.5},
        {10, 4},
        {13, 6},
        {16, 3},
    };
*/

    const QVector<QVector<double>> control_points_1 // Точки определяющего многоугольника
    {
        {1, 1},
        {3, 1.5},
        {5, 4.5},
        {7.5, 2},
        {10, 4},
        {12, 1.5}
    };

    const QVector<double> w_1 {1, 1, 1, 1, 1, 1}; // Весовые коэффициенты
    const uint degree_1 = 5;                   // Степень аппроксимирующих полиномов
    const int number_u_1 = 60;                 // Кол-во разбиений (точек) в реальной части узлов. вектора

    const QVector<double> u_1 = u_fill(control_points_1, degree_1); // Узловой вектор

    QVector<QVector<double>> nders_1(degree_1 + 1, QVector<double>(degree_1 + 1)); // Содержит для заданного "u" массив BASIS функций и 1-ую и 2-ую производную
    QVector<QPair<double, double>> c2_1(degree_1 + 1);  // Индекс 2 для 2D задачи
    QVector<Point_curve> data_NURBS_1(number_u_1 + 1);    // Содержит точки кривой, 1-ой и 2-ой производной

    const uint n_ver_1 = control_points_1.size(); // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
    const uint n_kn_1 = n_ver_1 + degree_1 + 1;     // Количество узлов (длина) в узловом векторе (n_knots)
    const uint n_real_1 = n_ver_1 - degree_1 + 1;   // Количество узлов (длина) реальной части узлового вектора

    for(int i = 0; i < number_u_1 + 1; ++i)
    {
        double u_i = (i / static_cast<double>(number_u_1));
        curve_point_and_deriv_NURBS(data_NURBS_1[i], n_real_1, n_kn_1, degree_1, u_1, control_points_1, w_1, u_i, c2_1, nders_1);
    }

    QString title = "Наибольшая метрика разности между 2-мя кривыми";
    QString labels_legend_1 = "Контрольные точки";
    QString labels_legend_2 = "Исходный NURBS";

    int x_min = -10, x_max = 10;
    int y_min = -10, y_max = 10;
    curve_plot(ui->graph_function, control_points_1, data_NURBS_1, x_min, x_max, y_min, y_max, title, labels_legend_1, labels_legend_2);

    //plot_polygon(ui->graph_function, control_points_1, labels_legend_1);

    const QVector<QVector<double>> control_points_2 // Точки определяющего многоугольника
    {
        {1, 1},
        {2.5, 3},
        {5.2, 5},
        {7, 1},
        {9.5, 5},
        {11, 3},
        {12, 1.5}
    };

    const QVector<double> w_2 {1, 1, 1, 1, 1, 1, 1}; // Весовые коэффициенты
    const uint degree_2 = 2;                   // Степень аппроксимирующих полиномов
    const int number_u_2 = 60;                 // Кол-во разбиений (точек) в реальной части узлов. вектора

    const QVector<double> u_2 = u_fill(control_points_2, degree_2); // Узловой вектор

    QVector<QVector<double>> nders_2(degree_2 + 1, QVector<double>(degree_2 + 1)); // Содержит для заданного "u" массив BASIS функций и 1-ую и 2-ую производную
    QVector<QPair<double, double>> c2_2(degree_2 + 1);  // Индекс 2 для 2D задачи
    QVector<Point_curve> data_NURBS_2(number_u_2 + 1);    // Содержит точки кривой, 1-ой и 2-ой производной

    const uint n_ver_2 = control_points_2.size(); // Количество вершин в определяющем многоугольнике (n_vertices) (отсчёт с 1)
    const uint n_kn_2 = n_ver_2 + degree_2 + 1;     // Количество узлов (длина) в узловом векторе (n_knots)
    const uint n_real_2 = n_ver_2 - degree_2 + 1;   // Количество узлов (длина) реальной части узлового вектора

    for(int i = 0; i < number_u_2 + 1; ++i)
    {
        double u_i = (i / static_cast<double>(number_u_2));
        curve_point_and_deriv_NURBS(data_NURBS_2[i], n_real_2, n_kn_2, degree_2, u_2, control_points_2, w_2, u_i, c2_2, nders_2);
    }

    labels_legend_2 = "Аппроксимирующий NURBS";
    plot_curve(ui->graph_function, data_NURBS_2, labels_legend_2, QColor(28, 172, 120));

    Point_curve max_p1, max_p2;
    double max_perpendicular = 0;
    for(int i = 0; i < data_NURBS_2.size() - 1; ++i)
    {
        Point_curve u_perpendicular = finding_perpendicular(n_real_1, n_kn_1, degree_1, u_1, control_points_1, w_1, data_NURBS_2[i].curve);
        //plot_line(ui->graph_function, data_NURBS_2[i].curve.first, data_NURBS_2[i].curve.second, u_perpendicular.curve.first, u_perpendicular.curve.second, QColor(178, 34, 34)); // Рисуем перпендикуляр между точкой и кривой
        double temp_max_perpendicular = vector_len(data_NURBS_2[i].curve, u_perpendicular);
        if(max_perpendicular < temp_max_perpendicular)
        {
            max_perpendicular = temp_max_perpendicular;
            max_p1 = u_perpendicular;
            max_p2 = data_NURBS_2[i];
        }
    }

    plot_line(ui->graph_function, max_p1.curve.first, max_p1.curve.second, max_p2.curve.first, max_p2.curve.second, QColor(178, 34, 34)); // Рисуем перпендикуляр между точкой и кривой
    qDebug() << "MAX_LEN: " << max_perpendicular;

    plot_lable_with_arrow(ui->graph_function, 2, 5, 4.812, 3.124, "Наибольшее расстояние\nмежду кривыми");




    /*
    QVector<QPair<double, double>> epsilon;
    for(const auto& p: data_NURBS_1)
        epsilon.push_back(calc_epsilon(p));

    for(const auto& p: epsilon)
        plot_point(ui->graph_function, p.first, p.second);
    */



    //QPair<double, double> unit_vector {data_NURBS_1[0].curve.first + 1, data_NURBS_1[0].curve.second};
    //double angle = vector_angle(rotated_points[0], data_NURBS_1[0], unit_vector);


    const double eps = 3;
    QVector<QPair<double, double>> rotated_points;
    QVector<Point_curve> epsilon_point_DATA (number_u_1 + 1);

    QVector<QPair<double, double>> rotated_points1;
    QVector<Point_curve> epsilon_point_DATA1 (number_u_1 + 1);

    for(int i = 0; i < data_NURBS_1.size(); ++i)
    {
        rotated_points.push_back(point_rotated(data_NURBS_1[i]));
        epsilon_point_DATA[i].curve = epsilon_point(rotated_points[i], data_NURBS_1[i], eps);

        rotated_points1.push_back(point_rotated(data_NURBS_1[i], - M_PI / 2));
        epsilon_point_DATA1[i].curve = epsilon_point(rotated_points1[i], data_NURBS_1[i], eps);
    }

    plot_curve(ui->graph_function, epsilon_point_DATA, "");
    plot_curve(ui->graph_function, epsilon_point_DATA1, "");

 //   plot_line(ui->graph_function, data_NURBS_1[0].curve.first,  data_NURBS_1[0].curve.second,
 //           data_NURBS_1[0].curve.first + x, data_NURBS_1[0].curve.second + y);

/*
    QVector<QPair<double, double>> epsilon2;
    epsilon2.push_back(point_rotated(data_NURBS_1[0], M_PI / -2));
    plot_point(ui->graph_function, epsilon2[0].first, epsilon2[0].second);
    derivative_point_line(ui->graph_function, data_NURBS_1[0]);
*/
}

Widget::~Widget()
{
    delete ui;
}
