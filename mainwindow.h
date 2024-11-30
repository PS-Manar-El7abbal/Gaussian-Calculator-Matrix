#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QString>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void updateTableSize();                // Slot to update table size based on spin box values
    void on_solveButton_clicked();         // Slot to handle the solve button click event
    void clearMatrix();  // Declare the clearMatrix slot

private:
    Ui::MainWindow *ui;                    // Pointer to the UI object
    void printMatrixToUI(const QVector<QVector<double>>& matrix, const QString& comment); // Function to display matrix in UI
    void displaySpecialMessage(const QString& message); // Display messages for special cases
    void swapRows(QVector<QVector<double>>& matrix, int row1, int row2); // Swap two rows in the matrix
    void normalizeRow(QVector<QVector<double>>& matrix, int row, double pivot); // Normalize a row by its pivot
    void eliminateBelow(QVector<QVector<double>>& matrix, int row); // Eliminate elements below the pivot
    void gaussianElimination(QVector<QVector<double>>& matrix, QVector<double>& solution); // Perform Gaussian elimination
    void gaussJordanElimination(QVector<QVector<double>>& matrix, QVector<double>& solution);
    void backSubstitution(QVector<QVector<double>>& matrix, QVector<double>& solution); // Perform back substitution
    int findPivotRow(const QVector<QVector<double>>& matrix, int row, int col); // Find the pivot row for Gaussian elimination
    bool inverseMatrix(const QVector<QVector<double>>& matrix, QVector<QVector<double>>& inverse);
    void solvebyInverse(QVector<QVector<double>>& matrix, QVector<double>& sol);

};

#endif // MAINWINDOW_H
