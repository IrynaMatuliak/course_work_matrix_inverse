#include "resultdialog.h"
#include "ui_resultdialog.h"
#include <QFileDialog>
#include <QMessageBox>

ResultDialog::ResultDialog(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::ResultDialog)
{
    ui->setupUi(this);
}

ResultDialog::~ResultDialog()
{
    delete ui;
}
void ResultDialog::setText(const QString& txt) {
    ui->ResultOutput->setText(txt);
}

void ResultDialog::on_SaveToFile_clicked()
{
    QString filename = QFileDialog::getSaveFileName(
        this,
        tr("Save matrix"),
        QDir::currentPath(),
        tr("Text files (*.txt);;All files (*.*)"));
    if (!filename.isNull())
    {
        out_stream wrap_stream(filename.toStdString());

        if (!wrap_stream->is_open()) {
            QMessageBox msg = QMessageBox();
            msg.setWindowTitle("Error");
            msg.setText(QString("Unable to open ") + filename);
            msg.setIcon(QMessageBox::Critical);
            msg.exec();
            return;
        }

        wrap_stream.getStream() << ui->ResultOutput->toPlainText().toStdString() << std::endl;
    }
}

