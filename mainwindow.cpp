#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QTableWidget>
#include <QVector>
#include <QString>
#include <QDebug>
#include <QMessageBox>
#include <QRegularExpression>

// Constructor for MainWindow
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Connect the spin boxes to the slot to update table dimensions
    connect(ui->spinBox, SIGNAL(valueChanged(int)), this, SLOT(updateTableSize()));
    connect(ui->spinBox_2, SIGNAL(valueChanged(int)), this, SLOT(updateTableSize()));
    connect(ui->clearmatrix, SIGNAL(clicked()), this, SLOT(clearMatrix()));
    connect(ui->solveButton, SIGNAL(clicked()), this, SLOT(solveSystem()));

    // Set initial size for tableWidget (2x2 to avoid smaller matrix)
    ui->spinBox->setValue(2);   // Default rows
    ui->spinBox_2->setValue(3);  // Default columns

    // Add constraints to ensure the rows and columns can't be smaller than 2
    ui->spinBox->setMinimum(2);  // Minimum rows = 2
    ui->spinBox_2->setMinimum(3); // Minimum columns = 2

    updateTableSize();
}

//clear matrix and solve
void MainWindow::clearMatrix()
{
    int rows = ui->tableWidget->rowCount();
    int cols = ui->tableWidget->columnCount();


    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            QTableWidgetItem* item = ui->tableWidget->item(i, j);
            if (item) {
                item->setText("0");  //clear each cell
            } else {
                ui->tableWidget->setItem(i, j, new QTableWidgetItem(""));

            }
        }
    }
    ui->outputTextEdit->clear();
}

// Slot function to update the table size
void MainWindow::updateTableSize()
{
    int rows = ui->spinBox->value();
    int cols = ui->spinBox_2->value();

    ui->tableWidget->setRowCount(rows);
    ui->tableWidget->setColumnCount(cols);

    // Initialize each cell with zero as default if it's empty
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (!ui->tableWidget->item(i, j)) {
                ui->tableWidget->setItem(i, j, new QTableWidgetItem("0"));
            }
        }
    }
}

// Function to print the matrix to QTableWidget with a given comment, with enhanced formatting
void MainWindow::printMatrixToUI(const QVector<QVector<double>>& matrix, const QString& comment) {
    ui->outputTextEdit->append("<b>" + comment + "</b>\n");

    // Print matrix with bold formatting for readability
    for (const auto& row : matrix) {
        QString rowStr = "<i>(";
        int lastColumn = row.size() - 1;
        for (int i = 0; i < row.size(); ++i) {
            QString numberStr = QString::number((std::abs(row[i]) < 1e-10 ? 0.0 : row[i]), 'f', 2); // 2 decimal precision and avoid "-0"
            rowStr += numberStr;
            if (i == lastColumn - 1) {
                rowStr += " | ";
            } else if (i < lastColumn) {
                rowStr += " ";
            }
        }
        rowStr += ")</i>";
        ui->outputTextEdit->append(rowStr);
    }
    ui->outputTextEdit->append("\n");
}

// Display special messages like no solution or infinite solutions
void MainWindow::displaySpecialMessage(const QString& message) {
    ui->outputTextEdit->append("<b style='color: #e74c3c;'>" + message + "</b>");
}

// Swap two rows and log the operation
void MainWindow::swapRows(QVector<QVector<double>>& matrix, int row1, int row2) {
    if (row1 != row2) {
        std::swap(matrix[row1], matrix[row2]);
        printMatrixToUI(matrix, QString("Swapped R%1 with R%2").arg(row1 + 1).arg(row2 + 1));
    }
}

// Normalize a specific row and log the operation, with check to avoid dividing by zero
void MainWindow::normalizeRow(QVector<QVector<double>>& matrix, int row, double pivot) {
    if (qAbs(pivot) < 1e-10) {
        displaySpecialMessage("No solution for this system (division by near-zero pivot).");
        return;
    }
    int cols = matrix[row].size();
    for (int i = 0; i < cols; ++i) {
        matrix[row][i] /= pivot;
    }
    printMatrixToUI(matrix, QString("Normalized R%1 by dividing by %2").arg(row + 1).arg(QString::number(pivot, 'f', 2)));
}

// Eliminate elements below the pivot in a given column and log the operation
void MainWindow::eliminateBelow(QVector<QVector<double>>& matrix, int row) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int i = row + 1; i < rows; ++i) {
        double factor = matrix[i][row] / matrix[row][row];
        for (int j = 0; j < cols; ++j)
            matrix[i][j] -= factor * matrix[row][j];
        printMatrixToUI(matrix, QString("Eliminated below pivot in R%1 by factor %2").arg(i + 1).arg(factor, 0, 'f', 2));
    }
}

// Perform Gaussian elimination and log each step until reaching row echelon form
void MainWindow::gaussianElimination(QVector<QVector<double>>& matrix, QVector<double>& solution) {
    int rows = matrix.size();
    bool infiniteSolution = false;
    bool noSolution = false;

    for (int i = 0; i < rows; ++i) {

        int pivotRow = findPivotRow(matrix, i, i); // Find pivot row
        swapRows(matrix, i, pivotRow);             // Swap current row with pivot row


        // Check if pivot element is effectively zero (using a small epsilon value to handle floating-point precision)
        if (qAbs(matrix[i][i]) < 1e-10) {
            bool isZeroRow = true;
            for (int j = 0; j < matrix[i].size() - 1; ++j) {
                if (qAbs(matrix[i][j]) > 1e-10) {
                    isZeroRow = false;
                    break;
                }
            }
            if (isZeroRow && qAbs(matrix[i].back()) > 1e-10) {
                noSolution = true;
                break;                // Inconsistent row, meaning no solution
            } else if (isZeroRow) {
                infiniteSolution = true;
                break;                // Consistent zero row indicates infinite solutions
            }
        }

        normalizeRow(matrix, i, matrix[i][i]); // Normalize the current row
        eliminateBelow(matrix, i);           // Eliminate elements below the pivot
    }

    // Print the matrix after reaching row echelon form
    printMatrixToUI(matrix, "Row Echelon Form:");

    if (noSolution) {
        displaySpecialMessage("No solution for this system.");
    } else if (infiniteSolution) {
        displaySpecialMessage("Infinite solutions for this system.");
    } else {
        // Proceed to back substitution if a unique solution exists
        backSubstitution(matrix, solution);
    }
}

// Perform Gauss-Jordan elimination and log each step
void MainWindow::gaussJordanElimination(QVector<QVector<double>>& matrix, QVector<double>& solution) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    bool noSolution = false;
    bool infiniteSolutions = false;

    for (int i = 0; i < rows; ++i) {
        // Find the pivot row and swap
        int pivotRow = findPivotRow(matrix, i, i);
        swapRows(matrix, i, pivotRow);

        // Check if the pivot is effectively zero
        if (qAbs(matrix[i][i]) < 1e-10) {
            bool isZeroRow = true;
            for (int j = 0; j < cols - 1; ++j) {
                if (qAbs(matrix[i][j]) > 1e-10) {
                    isZeroRow = false;
                    break;
                }
            }
            if (isZeroRow && qAbs(matrix[i].back()) > 1e-10) {
                noSolution = true;
                break; // Inconsistent row
            } else if (isZeroRow) {
                infiniteSolutions = true;
                continue; // Skip to the next row
            }
        }

        // Normalize the pivot row
        normalizeRow(matrix, i, matrix[i][i]);

        // Eliminate all other rows in the current column
        for (int j = 0; j < rows; ++j) {
            if (j != i) {
                double factor = matrix[j][i];
                for (int k = 0; k < cols; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
                printMatrixToUI(matrix, QString("Eliminated R%1 using R%2").arg(j + 1).arg(i + 1));
            }
        }
    }

    // Print the matrix after RREF
    printMatrixToUI(matrix, "Reduced Row Echelon Form:");

    if (noSolution) {
        displaySpecialMessage("No solution for this system.");
    } else if (infiniteSolutions) {
        displaySpecialMessage("Infinite solutions for this system.");
    } else {
        // Extract solution from the last column of the matrix
        solution.resize(rows);
        for (int i = 0; i < rows; ++i) {
            solution[i] = matrix[i].back();
        }

        // Output the solution with bold formatting
        ui->outputTextEdit->append("<b>Solution:</b>");
        for (int i = 0; i < rows; ++i) {
            QString solutionStr = QString("x%1 = %2").arg(i + 1).arg((std::abs(solution[i]) < 1e-10 ? 0.0 : solution[i]), 0, 'f', 2);
            ui->outputTextEdit->append("<i>" + solutionStr + "</i>");
        }
    }
}


// Back substitution to find solution after Gaussian elimination with formatted solution output
void MainWindow::backSubstitution(QVector<QVector<double>>& matrix, QVector<double>& solution) {
    int n = matrix.size();
    solution.resize(n, 0);

    for (int i = n - 1; i >= 0; --i) {
        solution[i] = matrix[i].back();
        for (int j = i + 1; j < n; ++j) {
            solution[i] -= matrix[i][j] * solution[j];
        }
        solution[i] /= matrix[i][i];
    }

    // Output the solution with bold formatting and precision
    ui->outputTextEdit->append("<b>Solution:</b>");
    for (int i = 0; i < n; ++i) {
        QString solutionStr = QString("x%1 = %2").arg(i + 1).arg((std::abs(solution[i]) < 1e-10 ? 0.0 : solution[i]), 0, 'f', 2);
        ui->outputTextEdit->append("<i>" + solutionStr + "</i>");
    }
}

// Function to calculate the inverse of matrix A
bool MainWindow::inverseMatrix(const QVector<QVector<double>>& matrix, QVector<QVector<double>>& inverse) {
    int n = matrix.size();
    inverse = matrix;  // Start with a copy of the input matrix

    // Create an identity matrix of size n
    QVector<QVector<double>> identity(n, QVector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
    }

    // Gaussian elimination to create row echelon form
    for (int i = 0; i < n; ++i) {
        int pivotRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (qAbs(inverse[j][i]) > qAbs(inverse[pivotRow][i])) {
                pivotRow = j;
            }
        }

        if (qAbs(inverse[pivotRow][i]) < 1e-10) {
            return false;  // Singular matrix (not invertible)
        }

        // Swap rows
        std::swap(inverse[i], inverse[pivotRow]);
        std::swap(identity[i], identity[pivotRow]);

        // Normalize pivot row
        double pivot = inverse[i][i];
        for (int j = 0; j < n; ++j) {
            inverse[i][j] /= pivot;
            identity[i][j] /= pivot;
        }

        // Eliminate below pivot
        for (int j = i + 1; j < n; ++j) {
            double factor = inverse[j][i];
            for (int k = 0; k < n; ++k) {
                inverse[j][k] -= factor * inverse[i][k];
                identity[j][k] -= factor * identity[i][k];
            }
        }
    }

    // Back substitution to create reduced row echelon form
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            double factor = inverse[j][i];
            for (int k = 0; k < n; ++k) {
                inverse[j][k] -= factor * inverse[i][k];
                identity[j][k] -= factor * identity[i][k];
            }
        }
    }

    // Copy the result into inverse matrix
    inverse = identity;
    return true;
}


// Function to solve the system using the inverse of A
void MainWindow::solvebyInverse(QVector<QVector<double>>& inverse, QVector<double>& sol) {
    int rows = ui->tableWidget->rowCount();
    int cols = ui->tableWidget->columnCount();

    // Ensure that the matrix is square (number of rows == number of columns)
    if (rows != cols - 1) {
        displaySpecialMessage("The system is not a square system. Inverse can't be calculated.");
        return;
    }

    QVector<QVector<double>> matrix(rows, QVector<double>(cols - 1));
    QVector<double> constants(rows);

    // Populate the matrix and the constants vector from the table
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols - 1; ++j) {
            matrix[i][j] = ui->tableWidget->item(i, j)->text().toDouble();
        }
        constants[i] = ui->tableWidget->item(i, cols - 1)->text().toDouble();
    }

    QVector<QVector<double>> inverse;
    bool invertible = inverseMatrix(matrix,inverse);  // Compute the inverse of the matrix

    if (!invertible) {
        displaySpecialMessage("Matrix is singular, so the system has no unique solution.");
        return;
    }

    // Multiply the inverse matrix with the constants vector to get the solution
    QVector<double> solution(rows);
    for (int i = 0; i < rows; ++i) {
        solution[i] = 0;
        for (int j = 0; j < rows; ++j) {
            solution[i] += inverse[i][j] * constants[j];
        }
    }

    // Output the solution
    ui->outputTextEdit->append("<b>Solution:</b>");
    for (int i = 0; i < rows; ++i) {
        QString solutionStr = QString("x%1 = %2").arg(i + 1).arg((std::abs(solution[i]) < 1e-10 ? 0.0 : solution[i]), 0, 'f', 2);
        ui->outputTextEdit->append("<i>" + solutionStr + "</i>");
    }
}


// Destructor
MainWindow::~MainWindow() {
    delete ui;
}
