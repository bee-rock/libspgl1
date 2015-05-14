namespace libspgl1 {
namespace vector {

template<typename VectorType>
size_t n_elem(const VectorType& vector){
	return vector.n_elem;
}

template<typename ElementType, typename VectorType>
ElementType get_element(const VectorType& vector, size_t index)
{
    return static_cast<ElementType>(vector(index));
}

template<typename ElementType, typename VectorType>
void set_element(VectorType& vector, size_t index, ElementType value)
{
	vector(index) = value;
}

template<typename ElementType, typename VectorType>
ElementType dot(const VectorType& a, const VectorType& b){
	size_t n_elems = n_elem(a);
	ElementType result{0};
	for (size_t i=0; i<n_elems; ++i)
		result += get_element<ElementType>(a, i) * get_element<ElementType>(b, i);
	return static_cast<ElementType>(result);
}

template<typename ElementType, typename VectorType>
ElementType norm(const VectorType& a, float p){
	ElementType result{0};
	size_t n_elems = n_elem(a);
	for (size_t i=0; i<n_elems; ++i)
		result += pow(abs(get_element<ElementType>(a, i)), p);
	return static_cast<ElementType>(pow(result,1/p));
}

} // matrix
} // libspgl1

