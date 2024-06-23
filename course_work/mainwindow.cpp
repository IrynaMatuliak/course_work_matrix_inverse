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

std::vector<std::vector<double>> initialMatrixValues = {
    { 12.8538, -14.0461,  -3.4774, -14.3966, -12.8699,  55.9413,  36.1487, -25.7588,  41.1328, -96.9730 },
    { -95.9739,  16.5915,  61.9532, -49.8691,  78.3523,  66.6649,  58.5234,   4.9418, -25.4568, -39.3329 },
    {  64.5612, -80.4812,  51.4682,  28.4054,  22.4191, -78.3781, -55.7360,  91.0555, -79.8147,  -5.0118 },
    { -81.5095, -62.2833,  68.7113, -87.0333,  30.8999, -36.7704,   6.8185, -45.9398,  66.7764,  77.7809 },
    { -42.2875, -59.1954,  82.8187,  72.7472,  -7.1253,  88.4588, -87.4979,  70.1233,  24.1458, -96.0872 },
    { -50.6628,  61.8374, -29.5576, -99.2569,  82.9759,  12.7770, -84.3253,  83.6093,  39.9766,  73.7821 },
    {  95.1122, -16.7203,  82.2234,  76.0680,  14.7936, -33.9456,  25.7464,  55.3051, -83.7945, -28.4071 },
    { -41.2227, -42.8707, -10.4209, -35.9189, -84.8385,  36.8642, -13.3561, -11.7313, -65.4956,  94.4538 },
    { -38.4501,  37.8566,  -6.5625, -23.2112, -29.2655,  15.2509,  20.9572, -53.2138, -11.8873,  36.9784 },
    {  75.3726,  69.7440, -57.1748, -97.2084,  -9.5201,  19.3610,  89.1597, -19.8904, -34.4917,  55.9287 }
};

void InverseOfMatrix::initializeMatrixWithValues(const std::vector<std::vector<double>>& values)
{
    int n = values.size();
    ui->Matrix->setColumnCount(n);
    ui->Matrix->setRowCount(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = values[i][j];
            (*_matrix)(j, i) = value;
            ui->Matrix->setItem(i, j, new QTableWidgetItem(QString().asprintf("%.4lf", value)));
        }
    }
}

InverseOfMatrix::InverseOfMatrix(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::InverseOfMatrix)
    , _matrix(new SquareMatrix(initialMatrixValues.size()))
{
    ui->setupUi(this);
    for (int i = 0; i < methods.size(); ++i) {
        ui->SetMethod->addItem(methods[i]);
    }
    ui->SetMethod->setCurrentIndex(0);

    for (int i = 0; i < sizeOfMatrix.size(); ++i) {
        ui->SetSize->addItem(QString().asprintf("%d", sizeOfMatrix[i]));
    }
    ui->SetSize->setCurrentIndex(sizeOfMatrix.size() - 1);

    for (int i = 0; i < precision.size(); ++i) {
        ui->SetAccuracy->addItem(QString().asprintf("%.0e", precision[i]));
    }
    ui->SetAccuracy->setCurrentIndex(0);

    srand(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));

    initializeMatrixWithValues(initialMatrixValues);
}

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

/*InverseOfMatrix::InverseOfMatrix(QWidget *parent)
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
}*/

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
