#include <random>
#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "resultdialog.h"
#include "helpdialog.h"
#include <QMessageBox>
#include <QFileDialog>
#include <cmath>
std::vector<const char*> methods {"LUP-decomposition", "Shults Method"};
std::vector<int> sizeOfMatrix {2, 3, 4, 5, 6, 7, 8, 9, 10};
std::vector<double> precision {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};

QString output_matrix(const QString& caption, const SquareMatrix& m)
{
    QString result = caption + "\n";
    for (size_t j = 0; j < m.size(); ++j) {
        for (size_t i = 0; i < m.size(); ++i) {
            double d = m(i, j);
            if ((d > -1e5 && d < -1e-5) || (d > 1e-5 && d < 1e5) || d == 0 || std::isnan(d)) {
                result += QString().asprintf("%15.4Lf", m[i][j]);
            } else {
                if ((d > -1e11 && d < -1e-11) || (d > 1e-11 && d < 1e11)) {
                    result += QString().asprintf("%15.4Le", m[i][j]);
                } else {
                    result += QString().asprintf("%15.4lf", 0.);
                }
            }
        }
        result += "\n";
    }
    result += "\n";
    return result;
}

InverseOfMatrix::InverseOfMatrix(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::InverseOfMatrix)
    , _matrix(new SquareMatrix(2))
{

    ui->setupUi(this);
    for (int i = 0; i < methods.size(); ++i) {
        ui->SetMethod->addItem(methods[i]);
    }
    ui->SetMethod->setCurrentIndex(0);

    for (int i = 0; i < sizeOfMatrix.size(); ++i) {
        ui->SetSize->addItem(QString().asprintf("%d", sizeOfMatrix[i]));
    }
    ui->SetSize->setCurrentIndex(0);

    for (int i = 0; i < precision.size(); ++i) {
        ui->SetAccuracy->addItem(QString().asprintf("%.0e", precision[i]));
    }
    ui->SetAccuracy->setCurrentIndex(0);

    srand(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
}
InverseOfMatrix::~InverseOfMatrix()
{
    delete ui;
    delete _matrix;
}

void InverseOfMatrix::on_pushButton_clicked()
{
    close();
}


void InverseOfMatrix::on_SetSize_currentIndexChanged(int index)
{
    delete _matrix;
    _matrix = new SquareMatrix(sizeOfMatrix[index]);
    ui->Matrix->setColumnCount(sizeOfMatrix[index]);
    ui->Matrix->setRowCount(sizeOfMatrix[index]);
}


void InverseOfMatrix::on_CalculateButton_clicked()
{
    ResultDialog dlg(this);
    QString result;

    QString prefix = "Expected double value ";
    for (size_t j = 0; j < _matrix->size(); ++j) {
        for (size_t i = 0; i < _matrix->size(); ++i) {
            double d = 0;
            bool ok = false;
            QString text;
            if (ui->Matrix->item(i, j)) {
                d = ui->Matrix->item(i, j)->text().toDouble(&ok);
                if (!((d >= -1e10 && d <= -1e-10) || (d >= 1e-10 && d <= 1e10) ||  d == 0)) {
                    ok = false;
                    prefix = "Value is out of range [-1e10; -1e-10] U [1e-10; 1e10] U {0}.\nInvalid double value ";
                }
                text = ui->Matrix->item(i, j)->text();
            }
            if (!ok) {
                QMessageBox msg = QMessageBox();
                msg.setWindowTitle("Error");
                msg.setText(prefix + "in cell (" + QString::number(i + 1) + ", "
                            + QString::number(j + 1) + "). Received '" + text  + "'");
                msg.setIcon(QMessageBox::Critical);
                msg.exec();
                return;
            }
            (*_matrix)(j, i) = d;
        }
    }

    if (_matrix->isSingular()) {
        QMessageBox msg = QMessageBox();
        msg.setWindowTitle("Error");
        msg.setText("Matrix is singular. Unable find inverse matrix.");
        msg.setIcon(QMessageBox::Critical);
        msg.exec();
        return;
    }

    int idx = ui->SetMethod->currentIndex();
    double prec = precision[ui->SetAccuracy->currentIndex()];
    Results r = _matrix->inverse(idx, prec);
    if (std::isnan(r.X[0][0])) {
        QMessageBox msg = QMessageBox();
        msg.setWindowTitle("Error");
        msg.setText("Method doesn't converge. Unable find inverse matrix.");
        msg.setIcon(QMessageBox::Critical);
        msg.exec();
        return;
    }
    result += formResults(r);
    dlg.setText(result);
    dlg.exec();
}

void InverseOfMatrix::on_DeleteButton_clicked()
{
    for (int i = 0; i < ui->Matrix->rowCount(); ++i) {
        for (int j = 0; j < ui->Matrix->columnCount(); ++j) {
            if (ui->Matrix->item(i, j)) {
                ui->Matrix->item(i, j)->setText("");
            }
        }
    }

}

void InverseOfMatrix::on_HelpButton_clicked()
{
    HelpDialog dlg(this);
    dlg.exec();
}

void InverseOfMatrix::on_GenerateRandomMatrix_clicked()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(-1000000, 1000000);
    for (int j = 0; j < ui->Matrix->rowCount(); ++j) {
        for (int i = 0; i < ui->Matrix->columnCount(); ++i) {
            double sign = (static_cast<double>(rand()) / RAND_MAX) > 0.5 ? 1. : -1.;
            //double val = (static_cast<double>(rand()) / RAND_MAX) * 1000. * sign;
            double value = static_cast<double>(distr(gen)) / 10007. * sign;
            ui->Matrix->setItem(i, j, new QTableWidgetItem(QString().asprintf("%.4lf", value)));
        }
    }
}

QString InverseOfMatrix::formResults(const Results& res)
{
    QString result;
    result += output_matrix("Original matrix:", res.A);
    switch(res.method) {
    case 0:
        result += output_matrix("LU matrix:", res.LU);
        result += output_matrix("P matrix:", res.P);
        break;
    case 1:
        result += "Precision: " + QString().asprintf("%Le\n", res.precision);
        break;
    }
    result += "Determinant: " + QString().asprintf("%15.4Lf\n", res.determinant);
    result += output_matrix("Inverted matrix:", res.X);
    result += output_matrix("Check A*A^(-1) matrix:", res.A * res.X);
    result += "Complexity of calculations: " + QString().asprintf("%ld\n", res.complexity);
    result += "\n";
    result += "\n";
    return result;
}
