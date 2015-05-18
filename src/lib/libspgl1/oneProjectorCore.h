#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>     /* provides DBL_EPSILON */
#include <sys/types.h>
#include <algorithm>

/* ----------------------------------------------------------------------- */
template<typename VectorType>
VectorType projectI(VectorType& c, const double tau, const size_t n)
/* ----------------------------------------------------------------------- */
{
   double
       b,          /* Current element of vector b */
       csb = 0,        /* Cumulative sum of b */
       alpha = 0,
       soft  = 0;  /* Soft thresholding value */
   c = libspgl1::vector::abs(c);
   auto c_bar = c;
//
//   /* The vector xPtr[] is initialized to bPtr[] prior to the function call */
//
//   /* Check if tau is essentially zero.  Exit with x = 0. */
   if (tau < DBL_EPSILON) {
       for (size_t i = 0; i < n; i++){
    	   libspgl1::vector::set_element(c_bar, i, 0);
       	   }
       return c_bar;
   }

   /* Check if ||b||_1 <= lambda.  Exit with x = b. */
   for (size_t i = 0; i < n; i++){
	   csb += libspgl1::vector::get_element<double>(c_bar, i);
   }
   if (csb <= tau){
       return c_bar;
   }
   std::make_heap(c_bar.begin(), c_bar.end());

   /* Initialise csb with -tau so we don't have to subtract this at every iteration */
   csb = -tau;
   /* Determine threshold value `soft' */
   for (size_t j = 0; j < n; j++)
   {
	  soft = alpha;
	  /* Get current maximum heap element         */
	  b = libspgl1::vector::get_element<double>(c_bar, 0);
      std::pop_heap(c_bar.begin(), c_bar.end());
      csb += b;          /* Update the cumulative sum of b           */

      /* Move heap to next maximum value */
      std::pop_heap(c_bar.begin(),c_bar.end());
      /* Compute the required step to satisfy the tau constraint */
      alpha  = csb / j;

      /* We are done as soon as the constraint can be satisfied    */
      /* without exceeding the current minimum value in `vector' b */
      if (alpha >= b){
          break;
      }
   }

   /* Set the solution by applying soft-thresholding with `soft' */
   for (size_t i = 0; i < n; i++)
   {  b = libspgl1::vector::get_element<double>(c, i);
      if (b <= soft){
           libspgl1::vector::set_element<double>(c_bar, i, 0);
      }
      else{
    	  libspgl1::vector::set_element<double>(c_bar, i, b - soft);
      }
   }

   return c_bar;
}
