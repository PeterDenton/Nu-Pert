#ifndef Matrix_H
#define Matrix_H

#include <cstddef>

template <class T>
class Matrix
{
	private:
		std::size_t rows, cols;
		T* buffer;
		std::size_t Get_index(std::size_t row, std::size_t col) const;

	public:
		Matrix(std::size_t rows, std::size_t cols);
		Matrix();
		Matrix(const Matrix&);
		~Matrix();

		T& operator () (std::size_t row, std::size_t col);
		T operator () (std::size_t row, std::size_t col) const;

		Matrix operator = (const Matrix other);
		Matrix operator += (const Matrix other);
		Matrix operator + (const Matrix other) const;
		Matrix operator * (const Matrix rhs) const;

		void Transpose();
		void Hermitian();

		std::size_t Get_rows() const;
		std::size_t Get_cols() const;
		T Get(std::size_t row, std::size_t col) const;
		void Set_all(T value);
};

Matrix<std::complex<double> > double2complex(Matrix<double> Matrix); // converts a double to a complex<double>

#endif
