#include <algorithm>
#include <cassert>
#include <cstddef>
#include <complex>

#include "Matrix.h"

template <class T>
Matrix<T>::Matrix(std::size_t rows, std::size_t cols)
: rows(rows), cols(cols), buffer(new T[rows * cols])
{
	Set_all(T());
}

template <class T>
Matrix<T>::Matrix()
: rows(0), cols(0)
{
	buffer = 0;
}

template <class T>
Matrix<T>::Matrix(const Matrix &other)
: rows(other.rows), cols(other.cols), buffer(new T[rows * cols])
{
	std::copy(other.buffer, other.buffer + rows * cols, buffer);
}

template <class T>
Matrix<T>::~Matrix()
{
	delete [] buffer;
}

template <class T>
T& Matrix<T>::operator () (std::size_t row, std::size_t col)
{
	assert(row >= 0 and row < rows);
	assert(col >= 0 and col < cols);
	return buffer[Get_index(row, col)];
}

template <class T>
T Matrix<T>::operator () (std::size_t row, std::size_t col) const
{
	assert(row >= 0 and row < rows);
	assert(col >= 0 and col < cols);
	return buffer[Get_index(row, col)];
}

template <class T>
Matrix<T> Matrix<T>::operator = (const Matrix other)
{
	rows = other.rows;
	cols = other.cols;
	buffer = new T[rows * cols];
	std::copy(other.buffer, other.buffer + rows * cols, buffer);
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator += (const Matrix other)
{
	assert(rows == other.rows);
	assert(cols == other.cols);

	for (std::size_t i = 0; i < rows * cols; i++)
		buffer[i] += other.buffer[i];
	return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator + (const Matrix other) const
{
	assert(rows == other.rows);
	assert(cols == other.cols);

	Matrix ret(cols, rows);
	for (std::size_t i = 0; i < rows * cols; i++)
		ret.buffer[i] = buffer[i] + other.buffer[i];
	return ret;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const Matrix rhs) const
{
	assert(cols == rhs.rows);

	Matrix ret(rows, rhs.cols);
	for (std::size_t i = 0; i < rows; i++)
	{
		for (std::size_t j = 0; j < rhs.cols; j++)
		{
			for (std::size_t k = 0; k < rows; k++)
			{
				ret(i, j) += this->Get(i, k) * rhs(k, j);
			}
		}
	}

	return ret;
}

template <class T>
void Matrix<T>::Transpose()
{
	T * tmp = new T[rows * cols];
	std::copy(buffer, buffer + rows * cols, tmp);

	for (std::size_t i = 0; i < rows; i++)
	{
		for (std::size_t j = 0; j < cols; j++)
			buffer[Get_index(i, j)] = tmp[Get_index(j, i)];
	}
	delete [] tmp;
}

template <class T>
void Matrix<T>::Hermitian()
{
	T * tmp = new T[rows * cols];
	std::copy(buffer, buffer + rows * cols, tmp);

	for (std::size_t i = 0; i < rows; i++)
	{
		for (std::size_t j = 0; j < cols; j++)
			buffer[Get_index(i, j)] = std::conj(tmp[Get_index(j, i)]);
	}
	delete [] tmp;
}

template <class T>
std::size_t Matrix<T>::Get_index(std::size_t row, std::size_t col) const
{
	return row * cols + col;
}

template <class T>
std::size_t Matrix<T>::Get_rows() const
{
	return rows;
}

template <class T>
std::size_t Matrix<T>::Get_cols() const
{
	return cols;
}

template <class T>
T Matrix<T>::Get(std::size_t row, std::size_t col) const
{
	assert(row >= 0 and row < rows);
	assert(col >= 0 and col < cols);
	return buffer[Get_index(row, col)];
}

template <class T>
void Matrix<T>::Set_all(T value)
{
	for (std::size_t i = 0; i < rows * cols; i++)
		buffer[i] = value;
}

Matrix<std::complex<double> > double2complex(const Matrix<double> old_matrix) // converts a double to a complex<double>
{
	Matrix<std::complex<double> > new_matrix(old_matrix.Get_rows(), old_matrix.Get_cols());

	for (std::size_t i = 0; i < old_matrix.Get_rows(); i++)
	{
		for (std::size_t j = 0; j < old_matrix.Get_cols(); j++)
		{
			new_matrix(i, j) = old_matrix(i, j);
		}
	}
	return new_matrix;
}

template class Matrix<double>;
template class Matrix<std::complex<double> >;
