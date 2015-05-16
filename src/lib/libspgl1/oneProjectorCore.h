#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>     /* provides DBL_EPSILON */
#include <sys/types.h>
#include "heap.h"

/* ----------------------------------------------------------------------- */
template<typename VectorType>
VectorType projectI(const VectorType& xPtr, const VectorType& bPtr, double tau, int n)
/* ----------------------------------------------------------------------- */
{  int
       i, j;
   double
       b,          /* Current element of vector b */
       csb,        /* Cumulative sum of b */
       alpha = 0,
       soft  = 0;  /* Soft thresholding value */
   VectorType x_return{libspgl1::vector::n_elem(xPtr)};
//
//   /* The vector xPtr[] is initialized to bPtr[] prior to the function call */
//
//   /* Check if tau is essentially zero.  Exit with x = 0. */
   if (tau < DBL_EPSILON) {
       for (size_t i = 0; i < n; i++)
    	   libspgl1::vector::set_element(x_return, i, 0);
       return x_return;
   }
//
//   /* Check if ||b||_1 <= lambda.  Exit with x = b. */
//   for (csb = 0, i = 0; i < n; i++) csb += bPtr[i];
//   if (csb <= tau)
//       return 0;
//
//   /* Set up the heap */
//   heap_build(n, xPtr);
//
//   /* Initialise csb with -tau so we don't have to subtract this at every iteration */
//   csb = -tau;
//
//   /* Determine threshold value `soft' */
//   for (i = n, j = 0; j < n; soft = alpha)
//   {
//      b = xPtr[0];       /* Get current maximum heap element         */
//      j ++;              /* Give compiler some room for optimization */
//      csb += b;          /* Update the cumulative sum of b           */
//
//      /* Move heap to next maximum value */
//      i = heap_del_max(i, xPtr);
//
//      /* Compute the required step to satisfy the tau constraint */
//      alpha  = csb / j;
//
//      /* We are done as soon as the constraint can be satisfied    */
//      /* without exceeding the current minimum value in `vector' b */
//      if (alpha >= b)
//          break;
//   }
//
//   /* Set the solution by applying soft-thresholding with `soft' */
//   for (i = 0; i < n; i++)
//   {  b = bPtr[i];
//      if (b <= soft)
//           xPtr[i] = 0;
//      else xPtr[i] = b - soft;
//   }
//
//   return j;
}


/* ----------------------------------------------------------------------- */
int projectD(double xPtr[], double bPtr[], double dPtr[], double dOrg[], double tau, int n)
/* ----------------------------------------------------------------------- */
{  int
       i, j;
   double
       csdb,        /* Cumulative sum of d.*b          */
       csd2,        /* Cumulative sum of d.^2          */
       b,           /* Current element of vector b     */
       d,           /* Current element of vector d     */
       bd,          /* Current element of vector b / d */
       alpha  = 0,
       soft   = 0;

   /* Check if tau is essentially zero.  Exit with x = 0. */
   if (tau < DBL_EPSILON)
   {   for (i = 0; i < n; i++) xPtr[i] = 0;
       return 0;
   }

   /* Preliminary check on trivial solution x = b (meanwhile, scale x) */
   for (csdb = 0, i = 0; i < n; i++)
   {  d = dPtr[i];
      b = xPtr[i];
      csdb += (d * b);
      xPtr[i] = b / d;
   }

   if (csdb <= tau)
   {
       /* Reset the entries of x to b */
      memcpy((void *)xPtr, (void *)bPtr, n * sizeof(double));
      return 0;
   }

   /* Set up the heap (we have to sort on b./d) */
   heap_build_2(n, xPtr, dPtr);

   /* Initialise csbd with -tau so we don't have to subtract this at every iteration */
   csdb = -tau;
   csd2 =  0;

   /* Determine the threshold level `soft' */
   for (i = n, j = 0; j < n; soft = alpha)
   {
      bd    = xPtr[0];        /* Get current maximum b / d                */
      j    ++;                /* Give compiler some room for optimization */
      d     = dPtr[0];        /* Get current value of d                   */
      d    *= d;              /* Compute d squared                        */
      csd2 += d;              /* Update the cumulative sum of d.*d        */
      csdb += bd * d;         /* Update the cumulative sum of d.*b        */

      /* Move heap to next maximum value */
      i = heap_del_max_2(i, xPtr, dPtr);

      /* Compute the required step to satisfy the lambda constraint */
      alpha  = csdb / csd2;

      /* We are done as soon as the constraint can be satisfied */
      /* without exceeding the current minimum value of b / d   */
      if (alpha >= bd) break;
   }

   /* Set the solution */
   for (i = 0; i < n; i++)
   {  b     = bPtr[i];
      alpha = dOrg[i] * soft; /* Use the original values of d here */
      if (b <= alpha)
           xPtr[i] = 0;
      else xPtr[i] = b - alpha;
   }

   return j;
}
