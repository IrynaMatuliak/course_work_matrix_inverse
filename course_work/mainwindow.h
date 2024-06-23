#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "matrix.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class InverseOfMatrix;
}
QT_END_NAMESPACE




class InverseOfMatrix : public QMainWindow
{
    Q_OBJECT

public:
    InverseOfMatrix(QWidget *parent = nullptr);
    ~InverseOfMatrix();

private slots:
    void on_pushButton_clicked();

    void on_SetSize_currentIndexChanged(int index);

    void on_CalculateButton_clicked();

    void on_DeleteButton_clicked();

    void on_HelpButton_clicked();

    void on_GenerateRandomMatrix_clicked();

private:
    Ui::InverseOfMatrix *ui;
    SquareMatrix* _matrix;
    QString formResults(const Results& res);
};
#endif // MAINWINDOW_H
