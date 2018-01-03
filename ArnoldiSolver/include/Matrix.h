#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <vector>

//template <class T, size_t ROW, size_t COL>
//using Matrix = std::vector<std::vector<T, COL>, ROW>;

namespace MyMatrix{

template <class T>
using Matrix = std::vector<std::vector<T>>;

}

#endif // MATRIX_H_INCLUDED
