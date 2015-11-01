#pragma once
#include "math.hpp"

namespace libspgl1 {
namespace vector {

template<typename VectorType>
size_t n_elem(const VectorType& vector)
{
    return static_cast<size_t>(vector.size());
}

template<typename ElementType, typename VectorType>
ElementType get_element(const VectorType& vector, size_t index)
{
    return static_cast<ElementType>(vector(index));
}

template<typename ElementType, typename VectorType>
void set_element(VectorType& vector, const size_t index, ElementType value)
{
	vector(index) = value;
}

template<typename VectorType>
VectorType abs(const VectorType& x) {
    VectorType x_abs = x;
    auto n = libspgl1::vector::n_elem(x);
    for (size_t i=0; i<n; i++){
    	auto current_element = libspgl1::vector::get_element<double>(x, i);
    	libspgl1::vector::set_element(x_abs, i, std::abs(current_element));
    }
    return x_abs;
}

template<typename ElementType, typename VectorType>
ElementType dot(const VectorType& a, const VectorType& b){
	const size_t n_elems = n_elem(a);
	ElementType result{0};
	for (size_t i=0; i<n_elems; ++i)
		result += libspgl1::vector::get_element<ElementType>(a, i) * libspgl1::vector::get_element<ElementType>(b, i);
	return static_cast<ElementType>(result);
}

template<typename VectorType>
VectorType sign(const VectorType& x) {
    size_t n = libspgl1::vector::n_elem(x);
    VectorType sign_x = x;
    for (size_t i=0; i<n; i++){
    	auto current_element = libspgl1::vector::get_element<double>(x, i);
    	if (current_element < 0.0)
    		libspgl1::vector::set_element(sign_x, i, -1.0);
    	else
    		libspgl1::vector::set_element(sign_x, i, 1.0);
    }
    return sign_x;
}

template<typename VectorType>
VectorType elementwise_multiplication(const VectorType& x, const VectorType& y) {
    size_t n = libspgl1::vector::n_elem(x);
    VectorType z = x;
    for (size_t i=0; i<n; i++){
    	auto x_elem = libspgl1::vector::get_element<double>(x, i);
    	auto y_elem = libspgl1::vector::get_element<double>(y, i);
    	libspgl1::vector::set_element(z, i, x_elem * y_elem);
    }
    return z;
}

template<typename VectorType>
VectorType elementwise_division(const VectorType& x, const VectorType& y) {
    size_t n = libspgl1::vector::n_elem(x);
    VectorType z = x;
    for (size_t i=0; i<n; i++){
    	auto x_elem = libspgl1::vector::get_element<double>(x, i);
    	auto y_elem = libspgl1::vector::get_element<double>(y, i);
    	libspgl1::vector::set_element(z, i, x_elem / y_elem);
    }
    return z;
}

} // matrix
} // libspgl1

