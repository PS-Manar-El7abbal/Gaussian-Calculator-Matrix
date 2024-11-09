#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_solveButton_clicked();
    void updateTableSize();
    void clearMatrix();

private:
    Ui::MainWindow *ui;
    void printMatrixToUI(const QVector<QVector<double>>& matrix, const QString& comment);
    void swapRows(QVector<QVector<double>>& matrix, int row1, int row2);
    void normalizeRow(QVector<QVector<double>>& matrix, int row, double pivot);
    void eliminateBelow(QVector<QVector<double>>& matrix, int row);
    void gaussianElimination(QVector<QVector<double>>& matrix, QVector<double>& solution);
    void backSubstitution(QVector<QVector<double>>& matrix, QVector<double>& solution);
    int findPivotRow(const QVector<QVector<double>>& matrix, int row, int col);

    bool has_NO_Solution(const QVector<QVector<double>>& matrix);
    bool has_Infinite_Solution(const QVector<QVector<double>>& matrix);
};

#endif // MAINWINDOW_H
