#include <iostream>
#include <stdexcept>
#include <vector>

class Matrix {
private:
    int rows_, cols_;
    std::vector<std::vector<double>> matrix_;

public:
    Matrix(){
        rows_ = 2;
        cols_ = 2;
        matrix_.resize(rows_, std::vector<double>(cols_, 0));
    }

    Matrix(const Matrix& other){
        rows_ = other.rows_;
        cols_ = other.cols_;
        matrix_ = other.matrix_;
    }

    Matrix(Matrix&& other) noexcept : matrix_(std::move(other.matrix_)) {
        rows_ = other.rows_;
        cols_ = other.cols_;
        other.rows_ = 0;
        other.cols_ = 0;
    }

    Matrix(int rows, int cols){
        rows_ = rows;
        cols_ = cols;
        matrix_.resize(rows_, std::vector<double>(cols_, 0));
    }

    ~Matrix() = default;

    int getCols() const {
        return cols_;
    }

    int getRows() const {
        return rows_;
    }

    void show(){
        for(auto & i : matrix_){
            for(double j : i){
                std::cout<<j<<" ";
            }
            std::cout<<"\n";
        }
    }


    bool EqMatrix(const Matrix& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            return false;
        }

        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_; j++) {
                if (matrix_[i][j] != other.matrix_[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    void SumMatrix(const Matrix& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("can't sum matrix");
        }
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                matrix_[i][j] += other.matrix_[i][j];
            }
        }
    }

    void SubMatrix(const Matrix& other) {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("matrices cannot be subtracted");
        }
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                matrix_[i][j] -= other.matrix_[i][j];
            }
        }
    }

    void MulNumber(const double num) {
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                matrix_[i][j] *= num;
            }
        }
    }

    void MulMatrix(const Matrix& other) {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("first matrix cols != second matrix rows");
        }
        Matrix result(rows_, other.cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < other.cols_; ++j) {
                for (int k = 0; k < cols_; ++k) {
                    result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
                }
            }
        }
        *this = result;
    }

    void Transpose() {
        Matrix result(cols_, rows_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                result.matrix_[j][i] = matrix_[i][j];
            }
        }
        *this = result;
    }

    double Determinant() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("the matrix must be square");
        }

        // Копируем матрицу, чтобы не изменять оригинальную
        std::vector<std::vector<double>> matrix = matrix_;
        double det = 1.0;

        for (int i = 0; i < rows_; ++i) {
            // Поиск ведущего элемента (ved)
            int ved = i;
            for (int j = i + 1; j < rows_; ++j) {
                if (std::abs(matrix[j][i]) > std::abs(matrix[ved][i])) {
                    ved = j;
                }
            }

            if (ved != i) {
                std::swap(matrix[i], matrix[ved]);
                det *= -1; // Меняем знак определителя при перестановке строк
            }

            if (matrix[i][i] == 0) {
                return 0; // Если ведущий элемент нулевой, определитель равен 0
            }

            det *= matrix[i][i]; // Умножаем определитель на ведущий элемент

            for (int j = i + 1; j < rows_; ++j) {
                double factor = matrix[j][i] / matrix[i][i];
                for (int k = i + 1; k < rows_; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }

        return det;
    }

    Matrix CalcComplements() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("the matrix must be square");
        }
        Matrix result(rows_, cols_);
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                Matrix subMat(rows_ - 1, cols_ - 1);
                int subRow = 0;
                for (int k = 0; k < rows_; ++k) {
                    if (k == i) continue;
                    int subCol = 0;
                    for (int l = 0; l < cols_; ++l) {
                        if (l == j) continue;
                        subMat(subRow, subCol) = matrix_[k][l];
                        subCol++;
                    }
                    subRow++;
                }
                double det = subMat.Determinant();
                result(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * det;
            }
        }
        return result;
    }

    Matrix InverseMatrix() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("the matrix must be square");
        }

        double det = Determinant();
        if (det == 0) {
            throw std::invalid_argument("det cant be = 0");
        }

        // Вычисляем матрицу алгебраических дополнений
        Matrix complements = CalcComplements();
        // Транспонируем её (получаем присоединённую матрицу)
        complements.Transpose();
        // Делим каждый элемент на определитель
        complements.MulNumber(1.0 / det);

        return complements;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result(*this);
        result.SumMatrix(other);
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result(*this);
        result.SubMatrix(other);
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result(*this);
        result.MulMatrix(other);
        return result;
    }

    Matrix operator*(const double num) const {
        Matrix result(*this);
        result.MulNumber(num);
        return result;
    }

    bool operator==(const Matrix& other) const {
        return EqMatrix(other);
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows_ = other.rows_;
            cols_ = other.cols_;
            matrix_ = other.matrix_;
        }
        return *this;
    }

    Matrix& operator+=(const Matrix& other) {
        SumMatrix(other);
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        SubMatrix(other);
        return *this;
    }

    Matrix& operator*=(const Matrix& other) {
        MulMatrix(other);
        return *this;
    }

    Matrix& operator*=(const double num) {
        MulNumber(num);
        return *this;
    }

    double& operator()(int i, int j) {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
            throw std::out_of_range("Index incorrect");
        }
        return matrix_[i][j];
    }

    const double& operator()(int i, int j) const {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
            throw std::out_of_range("Index incorrect");
        }
        return matrix_[i][j];
    }
};


int main() {
    // Создание матриц
    Matrix mat1(2, 2); // Матрица 2x2, заполненная нулями
    Matrix mat2(2, 2);

    // Заполнение матриц значениями
    mat1(0, 0) = 1; mat1(0, 1) = 2;
    mat1(1, 0) = 3; mat1(1, 1) = 4;

    mat2(0, 0) = 5; mat2(0, 1) = 6;
    mat2(1, 0) = 7; mat2(1, 1) = 8;

    // Вывод матриц
    std::cout << "Matrix 1:" << std::endl;
    mat1.show();
    std::cout << "\nMatrix 2:" << std::endl;
    mat2.show();

    // Сравнение матриц
    std::cout << "\nMatrix 1 == Matrix 2: " << (mat1 == mat2) << std::endl;

    // Сложение матриц
    Matrix matSum = mat1 + mat2;
    std::cout << "\nMatrix 1 + Matrix 2:" << std::endl;
    matSum.show();

    // Вычитание матриц
    Matrix matSub = mat1 - mat2;
    std::cout << "\nMatrix 1 - Matrix 2:" << std::endl;
    matSub.show();

    // Умножение матриц
    Matrix matMul = mat1 * mat2;
    std::cout << "\nMatrix 1 * Matrix 2:" << std::endl;
    matMul.show();

    // Умножение матрицы на число
    Matrix matMulNum = mat1 * 2.5;
    std::cout << "\nMatrix 1 * 2.5:" << std::endl;
    matMulNum.show();

    // Транспонирование матрицы
    Matrix matTransposed = mat1;
    matTransposed.Transpose();
    std::cout << "\nTransposed Matrix 1:" << std::endl;
    matTransposed.show();

    // Определитель матрицы
    double det = mat1.Determinant();
    std::cout << "\nDeterminant of Matrix 1: " << det << std::endl;

    // Матрица алгебраических дополнений
    Matrix matComplements = mat1.CalcComplements();
    std::cout << "\nCofactor Matrix of Matrix 1:" << std::endl;
    matComplements.show();

    // Обратная матрица
    Matrix matInverse = mat1.InverseMatrix();
    std::cout << "\nInverse of Matrix 1:" << std::endl;
    matInverse.show();

    // Проверка оператора присваивания
    Matrix mat3 = mat1;
    std::cout << "\nMatrix 3 (copy of Matrix 1):" << std::endl;
    mat3.show();

    // Проверка оператора +=
    mat3 += mat2;
    std::cout << "\nMatrix 3 += Matrix 2:" << std::endl;
    mat3.show();

    // Проверка оператора -=
    mat3 -= mat2;
    std::cout << "\nMatrix 3 -= Matrix 2:" << std::endl;
    mat3.show();

    // Проверка оператора *=
    mat3 *= mat2;
    std::cout << "\nMatrix 3 *= Matrix 2:" << std::endl;
    mat3.show();

    // Проверка оператора *=
    mat3 *= 2.0;
    std::cout << "\nMatrix 3 *= 2.0:" << std::endl;
    mat3.show();

    return 0;
}
