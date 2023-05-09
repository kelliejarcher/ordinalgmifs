/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_vector_sum.h

	Purpose:  Provide a (hopefully!) faster version of part of
				Dr. Kellie Archer's R code.
				This is to replace sum on a vector in R.

	Never forget, this is for COLUMN-MAJOR arrays!

 ***********************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_vector_sum(SEXP);
SEXP vector_sum(double*, int);

/////////////////////////////////////////////////////////////////////////
//  function definitions  //
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  do_vector_sum does the setup to call vector_sum
//
//  Inputs:
//    SEXP V, an R Vector to be summed
//	Outputs:
//    SEXP S, a double value containing the sum of the doubles
//      in input Vector V.
/////////////////////////////////////////////////////////////////////////
SEXP do_vector_sum(SEXP V) {
	if (!isVector(V)) {
		Rprintf("C code do_vector_sum:  Oops, please pass argument ");
		Rprintf("as an R vector:\n");
		Rprintf("\tdo_vector_sum(vector)\n");
		return R_NilValue;
	} 
	//  check to make sure something was actually passed in...
	if (!V) {
		Rprintf("Oops, can't use an empty vector in ");
		Rprintf("C code do_vector_sum...\n");
		return R_NilValue;
	}
	int protections=0, len;
	double* pV;

	PROTECT(V = AS_NUMERIC(V)); 
	protections++;
	pV = NUMERIC_POINTER(V);
	//  get the length of the vector
	len = length(V);

	SEXP S = vector_sum(pV, len);
	UNPROTECT(protections);
	return S;
} // do_vector_sum
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  vector_sum is used to sum the doubles stored in an R Vector.
//
//  Inputs:
//    double* V, containing the doubles to sum, and
//    int len, the length of V; the number of doubles to sum
//	Outputs:
//    SEXP S, a double value containing the sum of the doubles
//      in input Vector V.
/////////////////////////////////////////////////////////////////////////
SEXP vector_sum(double* V, int len) {
	if (!V) {
		Rprintf("C code vector_sum_rows:  Can't use NULL vector!\n");
		return R_NilValue;
	}
	//  keep track of any R/S types created
	int protections = 0;
	//  return a Vector of length 1 containing the sum:
	SEXP S;
	PROTECT(S = allocVector(REALSXP, 1));
	protections++;
	double* pS = NUMERIC_POINTER(S);
	//  uh-oh
	if (!pS) {
		Rprintf("C code vector_sum:  Couldn't allocate");
		Rprintf("vector to return!\n");
		return R_NilValue;
	}
	//  use some pointer math for speed:
	//  beginning & end markers of vector:
	double* pV = &*V;
	double* Vend = &*pV + len;
	//  step through items in vector:
	//double* pVi;
	//  use fastest memory in machine if available
	double sum = 0.0;
	for (; pV<Vend; pV++) {
	//for (pV=&*V; pV<Vend; pV++) {  -> 03/15/2023
		sum += *pV;
	}
	*pS = sum;
	UNPROTECT(protections);
	return S;
}  // vector_sum
/////////////////////////////////////////////////////////////////////////
