#ifndef RESULTDIALOG_H
#define RESULTDIALOG_H

#include <QDialog>
#include <fstream>
namespace Ui {
class ResultDialog;
}

class out_stream {
public:
    out_stream(const std::string& fname) : outputFile(fname) { }
    ~out_stream() {
        outputFile.close();
    }
    std::ofstream& getStream() {
        return outputFile;
    }
    std::ofstream* operator->() {
        return &outputFile;
    }
private:
    std::ofstream outputFile;
};

class ResultDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ResultDialog(QWidget *parent = nullptr);
    ~ResultDialog();

    void setText(const QString& txt);
private slots:
    void on_SaveToFile_clicked();

private:
    Ui::ResultDialog *ui;
};

#endif // RESULTDIALOG_H
