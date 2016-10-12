#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <array>

template <class T, size_t ROW, size_t COL>
using Matrix = std::array<std::array<T, COL>, ROW>;
Matrix<float, 3, 4> mat;

#endif // MATRIX_H_INCLUDED
