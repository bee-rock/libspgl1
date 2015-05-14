namespace libspgl1 {
namespace math {

template<typename ElementType, typename VectorType>
ElementType NormL1_primal(const VectorType& a, const VectorType& b)
{
	size_t n_elems = libspgl1::vector::n_elem(a);
	ElementType result{0};
	VectorType weighted_vector(n_elems, 1);
	for (size_t i=0; i<n_elems; ++i)
		result += abs(libspgl1::vector::get_element<ElementType>(a, i)) * libspgl1::vector::get_element<ElementType>(b, i);
	return result;
}

} // math
} // libspgl1
