namespace libspgl1 {

	class ActiveVars {
		size_t nnzOld, nnzX, nnzG;//, nnzIdx, nnzDiff;

		ActiveVars(libspgl1::Parameters parameters) :
			nnzOld{0}, nnzX{0}, nnzG{0}//, nnzIdx{0}, nnzDiff{0}
		{}

		size_t compute_nnzX();
		size_t compute_nnzG();
		size_t compute_nnzIdx();
		size_t compute_nnzDiff();
	};

} // libspgl1
