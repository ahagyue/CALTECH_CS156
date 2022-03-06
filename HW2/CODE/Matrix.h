#ifndef __MATRIX_H
#define __MATRIX_H

#include <iostream>
#include <vector>
#include <string>
#include <random>

using namespace std;

template<typename V> class Matrix {
private:
    vector<vector<V>> mat;
    unsigned row, col;

public:
    // 1. Constructor / Destructor
    Matrix();
    Matrix(unsigned, unsigned, const V&);
    Matrix(vector<vector<V>>);
    Matrix(const Matrix<V>&);
    ~Matrix();


    // 2-1. Arithmetic Operations 
    Matrix<V>& operator= (const Matrix<V>&);
    Matrix<V> operator> (const V&) const;
    Matrix<V> operator< (const V&) const;

    void operator+= (const Matrix<V>&);
    void operator*= (const Matrix<V>&);
    void operator-= (const Matrix<V>&);

    void operator+= (const V);
    void operator*= (const V);
    void operator-= (const V);

    Matrix<V>& operator- ();

    // 2-2 Arithmetic Operators 
    Matrix<V> operator+ (const Matrix<V>&)const;
    Matrix<V> operator- (const Matrix<V>&)const;
    Matrix<V> operator* (const Matrix<V>&)const;

    template<typename U> friend Matrix<U> operator+(const U, const Matrix<U>&);
    template<typename U> friend Matrix<U> operator-(const U, const Matrix<U>&);
    template<typename U> friend Matrix<U> operator*(const U, const Matrix<U>&);

    template<typename U> friend Matrix<U> operator+(const Matrix<U>&, const U);
    template<typename U> friend Matrix<U> operator-(const Matrix<U>&, const U);
    template<typename U> friend Matrix<U> operator*(const Matrix<U>&, const U);

    // 3. Mathematic Tools

    Matrix<V> T() const;
    Matrix<double> inverse() const;
    Matrix<double> pseudo_inverse() const;
    Matrix<double>* LUP_decomposition() const;
    static Matrix<V> eye(unsigned);

    Matrix<V> diagonal() const;
    Matrix<V> get_line(const unsigned&, int) const;

    V sum(const unsigned&, int) const;
    double determinant() const;

    void concatenate(const Matrix<V>&, int);
    void random(V, V, string);
    void switching(const unsigned&, const unsigned&, int);

    //4. Utils
    template<typename U> friend ostream& operator<<(ostream&, const Matrix<U>&);
    template<typename U> friend istream& operator>>(istream&, Matrix<U>&);

    template<typename U> Matrix<U> typecast() const;

    V& operator()(const unsigned&, const unsigned&);
    const V& operator()(const unsigned&, const unsigned&) const;

    unsigned get_row() const;
    unsigned get_col() const;

    void rebuild_mat(unsigned, unsigned, V);
    void clear();
};


// --------------------- EXCEPTION DEFINITION ---------------------

class MatrixSizeException : public exception {
private:
    string message;
public:
    MatrixSizeException(string _m) : message("Matrix Size Error : " + _m) {}
    virtual const char* what() const throw() { return message.c_str(); }
};

class MatrixInputException : public exception {
private:
    string message;
public:
    MatrixInputException(string _m) : message("Matrix Input Error : " + _m) {}
    virtual const char* what() const throw() { return message.c_str(); }
};


// --------------------- METHOD DEFINITION ---------------------

// 1. Constructor / Destructor

template<typename V>
inline Matrix<V>::Matrix()
{
    row = 0;
    col = 0;
}

template<typename V>
Matrix<V>::Matrix(unsigned row, unsigned col, const V& initial) {

    this->row = row;
    this->col = col;

    mat.resize(row);
    for (unsigned i = 0; i < mat.size(); i++) {
        mat[i].resize(col, initial);
    }
}


template<typename V>
Matrix<V>::Matrix(vector<vector<V>> dMat) {
    row = unsigned(dMat.size());
    col = unsigned(dMat[0].size());
    mat = dMat;
}

template<typename V>
Matrix<V>::Matrix(const Matrix<V>& nMat) {
    row = nMat.row;
    col = nMat.col;
    mat = nMat.mat;
}


template<typename V>
Matrix<V>::~Matrix() {
    vector<vector<V>>().swap(mat);
}

// 2-1. Arithmetic Operations

template<typename V>
Matrix<V>& Matrix<V>::operator= (const Matrix<V>& nMat) {
    row = nMat.row;
    col = nMat.col;
    mat = nMat.mat;

    return *this;
}

template<typename V>
Matrix<V> Matrix<V>::operator> (const V& dead_point) const {
    Matrix<V> bool_mat(row, col, (V)(false));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            bool_mat(i, j) = (V)(mat[i][j] > dead_point);
        }
    }
    return bool_mat;
}

template<typename V>
Matrix<V> Matrix<V>::operator< (const V& dead_point) const {
    Matrix<V> bool_mat(row, col, (V)(false));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            bool_mat(i, j) = (V)(mat[i][j] < dead_point);
        }
    }
    return bool_mat;
}

template<typename V>
void Matrix<V>::operator+= (const Matrix<V>& nMat) {
    *this = *this + nMat;
}

template<typename V>
void Matrix<V>::operator*= (const Matrix<V>& nMat) {
    *this = *this * nMat;
}

template<typename V>
void Matrix<V>::operator-= (const Matrix<V>& nMat) {
    *this = *this - nMat;
}

template<typename V>
void Matrix<V>::operator+= (const V x) {
    *this = *this + x;
}

template<typename V>
void Matrix<V>::operator*= (const V x) {
    *this = *this * x;
}

template<typename V>
void Matrix<V>::operator-= (const V x) {
    *this = *this - x;
}

// 2-2 Arithmetic Operators

template<typename V>
Matrix<V> Matrix<V>::operator+ (const Matrix<V>& nMat)const {
    unsigned m_row = max(row, nMat.row);
    unsigned m_col = max(col, nMat.col);

    vector<vector<V>> m(m_row, vector<V>(m_col, 0));

    for (unsigned i = 0; i < m_row; i++) {
        for (unsigned j = 0; j < m_col; j++) {
            m[i][j] += (mat[i % row][j % col] + nMat.mat[i % nMat.row][j % nMat.col]);
        }
    }
    Matrix<V> add_mat(m);
    return add_mat;
}

template<typename V>
Matrix<V> Matrix<V>::operator- (const Matrix<V>& nMat)const {
    unsigned m_row = max(row, nMat.row);
    unsigned m_col = max(col, nMat.col);

    vector<vector<V>> m(m_row, vector<V>(m_col, 0));

    for (unsigned i = 0; i < m_row; i++) {
        for (unsigned j = 0; j < m_col; j++) {
            m[i][j] += (mat[i % row][j % col] - nMat.mat[i % nMat.row][j % nMat.col]);
        }
    }
    Matrix<V> sub_mat(m);
    return sub_mat;
}

template<typename V>
Matrix<V> Matrix<V>::operator* (const Matrix<V>& nMat)const {
    if (col != nMat.row) {
        throw MatrixSizeException("The number of columns in first matrix is " + to_string(col) +
            ", but the number of rows in second matrix is " + to_string(row));
    }
    else {
        vector<vector<V>> m(row, vector<V>(nMat.col, 0));

        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j < nMat.col; j++) {
                V element = 0;
                for (unsigned k = 0; k < col; k++)
                    element += mat[i][k] * nMat.mat[k][j];
                m[i][j] = element;
            }
        }
        Matrix<V> mul_mat(m);
        return mul_mat;
    }
}

template<typename V>
Matrix<V>& Matrix<V>::operator- () {
    for (unsigned i = 0; i < row; i++) {
        for (unsigned j = 0; j < col; j++) {
            mat[i][j] = -mat[i][j];
        }
    }

    return *this;
}

template<typename U>
Matrix<U> operator+(const U x, const Matrix<U>& nMat) {
    Matrix<U> add_mat(nMat);

    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            add_mat.mat[i][j] += x;
        }
    }
    return add_mat;
}

template<typename U>
Matrix<U> operator-(const U x, const Matrix<U>& nMat) {
    Matrix<U> sub_mat(nMat);

    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            sub_mat.mat[i][j] -= x;
        }
    }
    return sub_mat;
}

template<typename U>
Matrix<U> operator*(const U x, const Matrix<U>& nMat) {
    Matrix<U> mul_mat(nMat);

    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            mul_mat.mat[i][j] *= x;
        }
    }
    return mul_mat;
}

template<typename U>
Matrix<U> operator+(const Matrix<U>& nMat, const U x) {
    Matrix<U> add_mat(nMat);

    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            add_mat.mat[i][j] += x;
        }
    }
    return add_mat;
}

template<typename U>
Matrix<U> operator-(const Matrix<U>& nMat, const U x) {
    Matrix<U> sub_mat(nMat);

    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            sub_mat.mat[i][j] -= x;
        }
    }
    return sub_mat;
}

template<typename U>
Matrix<U> operator*(const Matrix<U>& nMat, const U x) {
    Matrix<U> mul_mat(nMat);

    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            mul_mat.mat[i][j] *= x;
        }
    }
    return mul_mat;
}


// 3. Mathematics Tools

template<typename V>
Matrix<V> Matrix<V>::T() const {
    Matrix<V> transposed(col, row, (V)0);

    for (unsigned i = 0; i < row; i++) {
        for (unsigned j = 0; j < col; j++) {
            transposed.mat[j][i] = mat[i][j];
        }
    }
    return transposed;
}

template<typename V>
inline Matrix<double> Matrix<V>::inverse() const
{
    if (col != row) {
        return (*this).pseudo_inverse();
    }
    else {
        Matrix<double> L_inv = Matrix<double>::eye(row);
        Matrix<double> U_inv = Matrix<double>::eye(row);

        Matrix<double>* LU = (*this).LUP_decomposition();
        Matrix<double> L = LU[0];
        Matrix<double> U = LU[1];

        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j <= i; j++) {
                L_inv(i, j) /= L(i, i);
                U_inv(row - i - 1, row - j - 1) /= U(row - i - 1, row - i - 1);
            }
           
            for (unsigned j = i+1; j < row; j++) {
                for (unsigned k = 0; k <= i; k++) {
                    L_inv(j, k) -= L_inv(i, k) * L(j, i);
                    U_inv(row - j - 1, row - k - 1) -= U_inv(row - i - 1, row - k - 1) * U(row - j - 1, row - i - 1);
                }
            }
        }

        return U_inv * L_inv;
    }
}

template<typename V>
inline Matrix<double> Matrix<V>::pseudo_inverse() const
{
    Matrix<V> transversed = (*this).T();
    Matrix<double> d_transversed = transversed.typecast<double>();
    return (transversed * (*this)).inverse() * d_transversed;
}

template<typename V>
inline Matrix<double>* Matrix<V>::LUP_decomposition() const{
    if (row != col) {
        throw MatrixSizeException("This is not square matrix. You can not use LUP decomposition function.");
    }
    else {
        Matrix<V> PA(*this);
        Matrix<double> L(row, col, 0);
        Matrix<double> U(row, col, 0);

        if (mat[0][0] == 0) {
            for (unsigned i = 0; i < row; i++) {
                if (mat[i][0] != 0) {
                    PA.switching(0, i, 0);
                    cout << "Switched 1st row and " + to_string(i + 1) + "th row. (Because LU decomposition is impossible.)";
                    break;
                }
            }
        }

        for (unsigned i = 0; i < row; i++) {
            L(i, 0) = double(PA(i, 0)) / PA(0, 0);
            L(i, i) = 1;
            U(0, i) = PA(0, i);
        }

        for (unsigned i = 1; i < row; i++) {
            for (unsigned j = i; j < row; j++) {
                double u = 0;
                double l = 0;
                for (unsigned k = 0; k < j; k++) {
                    u += L(i, k) * U(k, j);
                    l += L(j, k) * U(k, i);
                }
                U(i, j) = PA(i, j) - u;
                L(j, i) = double(PA(j, i) - l) / U(i, i);
            }
        }
        Matrix<double> *LU = new Matrix<double>[2];
        LU[0] = L;
        LU[1] = U;
        return LU;
    }
}

template<typename V>
inline Matrix<V> Matrix<V>::eye(unsigned n)
{
    Matrix<V> I(n, n, V(0));
    for (unsigned i = 0; i < n; i++)
        I(i, i) = V(1);
    return I;
}

template<typename V>
inline Matrix<V> Matrix<V>::diagonal() const
{
    Matrix<V> diag(1, min(row, col), 0);

    for (unsigned i = 0; i < min(row, col); i++)
        diag(0, i) = mat[i][i];

    return diag;
}

template<typename V>
inline Matrix<V> Matrix<V>::get_line(const unsigned& n, int axis) const{
    if (axis == 0) {
        Matrix<V> line(1, col, 0);
        for (unsigned i = 0; i < col; i++)
            line(0, i) = mat[n][i];
        return line;
    }
    else if(axis == 1) {
        Matrix<V> line(1, row, 0);
        for (unsigned i = 0; i < row; i++)
            line(0, i) = mat[i][n];
        return line;
    }
    else {
        throw MatrixSizeException("Matrix is 2D. You need to reset axis");
    }
}

template<typename V>
inline V Matrix<V>::sum(const unsigned&, int) const
{
    return V();
}

template<typename V>
inline double Matrix<V>::determinant() const {
    Matrix<double>* LU = (*this).LUP_decomposition();
    Matrix<double> diag = LU[1].diagonal();
    double det = 1;
    for (unsigned i = 0; i < row; i++)
        det *= diag(0, i);
    return det;
}

template<typename V>
inline void Matrix<V>::concatenate(const Matrix<V>& nMat, int axis){
    if(axis == 0) {
        if(row != nMat.get_row())
            throw MatrixSizeException("Row number should be same for row concatenation.");
        
        Matrix<V> newMat(row, col + nMat.get_col(), 0);
        for(unsigned i = 0; i < row; i++) {
            for(unsigned j = 0; j < col; j++) 
                newMat(i, j) = mat[i][j];
            for(unsigned j = col; j < col + nMat.get_col(); j++)
                newMat(i, j) = nMat(i, j - col);
        }
        (*this) = newMat;
    } 
    else if(axis == 1) {
        if(col != nMat.get_col())
            throw MatrixSizeException("Column number should be same for column concatenation.");
        
        Matrix<V> newMat(row + nMat.row, col, 0);
        for(unsigned i = 0; i < col; i++) {
            for(unsigned j = 0; j < row; j++) 
                newMat(j, i) = mat[j][i];
            for(unsigned j = row; j < row + nMat.get_row(); j++)
                newMat(j, i) = nMat(j - row, i);
        }
        (*this) = newMat;
    } 
    else {
        throw MatrixSizeException("Matrix is 2D. You need to reset axis");
    }
}

template<typename V>
void Matrix<V>::random(V s, V e, string distribution) {

    random_device rd;
    mt19937 gen(rd());

    if (distribution.compare("uniform") == 0) {
        uniform_real_distribution<V> dist(s, e);
        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j < col; j++)
                mat[i][j] = V(dist(gen));
        }
    }
    else if (distribution.compare("normal")) {
        // s is mean, e is std
        normal_distribution<double> dist(s, e);
        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j < col; j++)
                mat[i][j] = V(dist(gen));
        }
    }
}

template<>
void Matrix<int>::random(int s, int e, string distribution) {

    random_device rd;
    mt19937 gen(rd());

    if (distribution.compare("uniform") == 0) {
        uniform_int_distribution<int> dist(s, e);
        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j < col; j++)
                mat[i][j] = dist(gen);
        }
    }
    else if (distribution.compare("normal")) {
        // s is mean, e is std
        normal_distribution<double> dist(s, e);
        for (unsigned i = 0; i < row; i++) {
            for (unsigned j = 0; j < col; j++)
                mat[i][j] = int(dist(gen));
        }
    }
}

template<typename V>
inline void Matrix<V>::switching(const unsigned& s, const unsigned& e, int axis) {
    if (axis == 0) {
        vector<V> r = mat[s];
        mat[s] = mat[e];
        mat[e] = r;
    }
    else if (axis == 1) {
        for (unsigned i = 0; i < row; i++) {
            V c = mat[i][s];
            mat[i][s] = mat[i][e];
            mat[i][e] = c;
        }
    }
    else {
        throw MatrixSizeException("Matrix is 2D. You need to reset axis");
    }
}

// 4. Utils

template<typename U>
ostream& operator<<(ostream& output, const Matrix<U>& nMat) {
    for (unsigned i = 0; i < nMat.row; i++) {
        for (unsigned j = 0; j < nMat.col; j++) {
            output << nMat.mat[i][j] << " ";
        }
        output << endl;
    }
    return output;
}

template<typename U>
istream& operator>>(istream& input, Matrix<U>& nMat) {
    for (unsigned i = 0; i < nMat.row; i++) {
        vector<U> new_line;
        for (unsigned j = 0; j < nMat.col; j++) {
            cin >> nMat.mat[i][j];
        }
    }
    return input;
}


template<typename V>
template<typename U>
inline Matrix<U> Matrix<V>::typecast() const{
    Matrix<U> casted(row, col, 0);
    for (unsigned i = 0; i < row; i++) {
        for (unsigned j = 0; j < col; j++) {
            casted(i, j) = (U)(mat[i][j]);
        }
    }
    return casted;
}

template<typename V>
V& Matrix<V>::operator()(const unsigned& r, const unsigned& c) {
    return mat[r][c];
}

template<typename V>
const V& Matrix<V>::operator()(const unsigned& r, const unsigned& c) const {
    return mat[r][c];
}

template<typename V>
unsigned Matrix<V>::get_row() const {
    return row;
}

template<typename V>
unsigned Matrix<V>::get_col() const {
    return col;
}

template<typename V>
void Matrix<V>::rebuild_mat(unsigned r, unsigned c, V initial) {
    row = r;
    col = c;
    
    mat.resize(row);
    for (unsigned i = 0; i < mat.size(); i++) {
        mat[i].resize(col, initial);
    }
}

template<typename V>
void Matrix<V>::clear() {
    row = 0;
    col = 0;
    vector<vector<V>>().swap(mat);
}
 
#endif
