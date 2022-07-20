/********************************************************************************
** Form generated from reading UI file 'widget.ui'
**
** Created by: Qt User Interface Compiler version 5.15.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_WIDGET_H
#define UI_WIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_Widget
{
public:
    QGridLayout *gridLayout_2;
    QGridLayout *gridLayout;
    QCustomPlot *graph_second_derivative;
    QCustomPlot *graph_function;
    QCustomPlot *graph_first_derivative;

    void setupUi(QWidget *Widget)
    {
        if (Widget->objectName().isEmpty())
            Widget->setObjectName(QString::fromUtf8("Widget"));
        Widget->resize(1800, 600);
        gridLayout_2 = new QGridLayout(Widget);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        graph_second_derivative = new QCustomPlot(Widget);
        graph_second_derivative->setObjectName(QString::fromUtf8("graph_second_derivative"));

        gridLayout->addWidget(graph_second_derivative, 0, 2, 1, 1);

        graph_function = new QCustomPlot(Widget);
        graph_function->setObjectName(QString::fromUtf8("graph_function"));

        gridLayout->addWidget(graph_function, 0, 0, 1, 1);

        graph_first_derivative = new QCustomPlot(Widget);
        graph_first_derivative->setObjectName(QString::fromUtf8("graph_first_derivative"));

        gridLayout->addWidget(graph_first_derivative, 0, 1, 1, 1);


        gridLayout_2->addLayout(gridLayout, 0, 0, 1, 1);


        retranslateUi(Widget);

        QMetaObject::connectSlotsByName(Widget);
    } // setupUi

    void retranslateUi(QWidget *Widget)
    {
        Widget->setWindowTitle(QCoreApplication::translate("Widget", "NURBS", nullptr));
    } // retranslateUi

};

namespace Ui {
    class Widget: public Ui_Widget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_WIDGET_H
