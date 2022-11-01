//mat.h
#ifndef mat_H
#define mat_H

#include <iostream>
using namespace std;


template<typename T>
class mat{

 public:
    mat(int NO_ROWS, int NO_COLS, T INIT);
    mat(const mat<T> &); // Not implemented (rule of three)
    mat<T>& operator=(const mat<T> & L); // Not implemented (rule of three)
    ~mat();

    bool operator==(const mat<T>& rhs);

    void print(ostream& ostr) const;
    void print_true(ostream& ostr) const;

    T* operator()(int rowNo) const;
    T& operator()(int rowNo, int colNo) const;

    // For boolean matrices:
    void compute_rref(); // rref = reduced row echelon form
    mat<T>& operator^=(const mat<T>& matrix_right);

    int get_no_rows() const;
    int get_no_cols() const;	
    int count_nonzero() const;

 private:
    int no_rows, no_cols;
    T** M;
    
    // For boolean matrices:
    int first_line_true(int rank, int j) const; // Returns the first i >= rank such that mat(i,j) = true
    void swap_rows(int i, int j);
};


template<typename T>
mat<T>::mat(int NO_ROWS, int NO_COLS, T INIT) {
    no_rows = NO_ROWS;
    no_cols = NO_COLS;

    M = new T*[no_rows];

    for (int i = 0; i < no_rows; i++) {
	M[i] = new T[no_cols];
	for (int j = 0; j < no_cols; j++) {
	    M[i][j] = INIT;
	}
    }

}

//Copy-constructor
template<typename T>
mat<T>::mat(const mat<T> & matrix) {
    no_rows = matrix.no_rows;
    no_cols = matrix.no_cols;

    M = new T*[no_rows];

    for (int i = 0; i < no_rows; i++) {
        M[i] = new T[no_cols];
        for (int j = 0; j < no_cols; j++) {
	    M[i][j] = matrix(i,j);
        }
    }
}

//Affectation
template<typename T>
mat<T>& mat<T>::operator=(const mat<T> & L) {
    if (this != &L){
	no_rows = L.no_rows;
	no_cols = L.no_cols;

	for (int i = 0; i < no_rows; i++) {
	    for (int j = 0; j < no_cols; j++) {
		(*this)(i,j) = L(i,j);
	    }
	}
    }
    return *this;
}

template<typename T>
mat<T>::~mat() {
    for (int i = 0; i < no_rows; i++) {
        delete [] M[i];
    }
    delete [] M;
}


// check for equality
template<typename T>
bool mat<T>::operator==(const mat<T>& matrix_right) {
    if (no_rows != matrix_right.no_rows || no_cols != matrix_right.no_cols) {
    	return false;
    } 

    for (int i = 0; i < no_rows; i++) {
        for (int j = 0; j < no_cols; j++) {
            if ((*this)(i,j) != matrix_right(i,j)) {
                return false;
            }
        }
    }

    return true;
}


// Prints the matrix beautifully
// template<typename T>
// void mat<T>::print() const {
//     std::cout << "Matrix: " << std::endl;
//     for (int i = 0; i < no_rows; i++) {
//         for (int j = 0; j < no_cols; j++) {
// 	    std::cout << "[" << (*this)(i,j) << "] ";
//         }
// 	std::cout << std::endl;
//     }
// }

template<typename T>
T* mat<T>::operator()(int rowNo) const {
#if DEBUG
    if (rowNo > no_rows) {
	throw std::invalid_argument("Out of bounds");
    }
#endif
    return this->M[rowNo];
}

template<typename T>
T& mat<T>::operator()(int rowNo, int colNo) const {
#if DEBUG
    if ((rowNo > no_rows) || (colNo > no_cols)) {
	cout << "Trying to access row: " << rowNo << " and col " << colNo << endl;
	throw std::invalid_argument("Out of bounds");
    }
#endif
    return this->M[rowNo][colNo];
}

template<typename T>
void mat<T>::swap_rows(int i, int j) {
    for (int a = 0; a < no_cols; a++) {
	int temp = (*this)(i,a);
	(*this)(i,a) = (*this)(j,a);
	(*this)(j,a) = temp;
    }
}


template<typename T>
void mat<T>::print_true(ostream& ostr) const {
    // ostr << "Matrix indices: " << std::endl;
    for (int i = 0; i < no_rows; i++) {
        for (int j = 0; j < no_cols; j++) {
            if ((*this)(i,j)) {
                ostr << "(" << i << "," << j << ") ";
            }
        }
    }
    // ostr << std::endl;
}

template<typename T>
void mat<T>::print(ostream& ostr) const {
    for (int i = 0; i < no_rows; i++) {
        for (int j = 0; j < no_cols; j++) {
            ostr << (*this)(i,j);
        }
        ostr << endl;
    }
}


template<typename T>
int mat<T>::get_no_rows() const{
    return no_rows;
}

template<typename T>
int mat<T>::get_no_cols() const{
    return no_cols;
}

template<typename T>
int mat<T>::count_nonzero() const{
    int sum = 0;
    for (int i = 0; i < no_rows; i++) {
        for (int j = 0; j < no_cols; j++) {
            if ((*this)(i,j)) {
                sum++;
            }
        }
    }
    return sum;
}

#endif
