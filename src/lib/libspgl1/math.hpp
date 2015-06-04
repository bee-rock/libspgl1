#pragma once
#include "vector_adapter.hpp"

namespace libspgl1 {
namespace math {

template<typename ElementType>
ElementType pow(const ElementType& x, const ElementType& p) {
    return std::pow(x, p);
}

template<typename ElementType, typename VectorType>
double norm(const VectorType& a, const double& p){
	double result{0};
	const size_t n_elems = libspgl1::vector::n_elem(a);
	for (size_t i=0; i<n_elems; ++i)
		result += libspgl1::math::pow<ElementType>(
				std::abs(libspgl1::vector::get_element<ElementType>(a, i)),p);
	return static_cast<ElementType>(libspgl1::math::pow<ElementType>(result, 1/p));
}

template<typename ElementType, typename VectorType>
ElementType NormL1_primal(const VectorType& a, const VectorType& b)
{
	size_t n_elems = libspgl1::vector::n_elem(a);
	ElementType result{0};
	for (size_t i=0; i<n_elems; ++i)
		result += std::abs(libspgl1::vector::get_element<ElementType>(a, i)) * libspgl1::vector::get_element<ElementType>(b, i);
	return result;
}

template<typename ElementType, typename VectorType>
double max(const VectorType& a){
	double current_max = a[0];
	const size_t n_elems = libspgl1::vector::n_elem(a);
	for (size_t i=0; i<n_elems; ++i){
		double current_elem = libspgl1::vector::get_element<ElementType>(a, i);
		if (current_elem > current_max){
			current_max = current_elem;
		}
	}
	return current_max;
}

template<typename ElementType, typename VectorType>
double min(const VectorType& a){
	double current_min = a[0];
	const size_t n_elems = libspgl1::vector::n_elem(a);
	for (size_t i=0; i<n_elems; ++i){
		double current_elem = libspgl1::vector::get_element<ElementType>(a, i);
		if (current_elem < current_min){
			current_min = current_elem;
		}
	}
}

} // math
} // libspgl1
