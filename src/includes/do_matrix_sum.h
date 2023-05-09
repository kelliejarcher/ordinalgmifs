/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_matrix_sum.h

	Purpose:  Provide a (hopefully!) faster version of part of 
				Dr. Kellie Archer's fn.cum function in her R code.
				This is to replace row sum on a matrix in R.

	Never forget, this is for COLUMN-MAJOR arrays!

 ***********************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_matrix_sum_rows(SEXP);
SEXP matrix_sum_rows(double*, int[2]);

/////////////////////////////////////////////////////////////////////////
//  function definitions  //
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  do_matrix_sum_rows does the setup to call matrix_sum_ptr.
//  This is just candy to make it easy for an R user to use this 
//    function.
//
//  Inputs:
//		A, an R matrix.
//  Outputs:
//		S, an R vector, containing in each element the sum of the
//			corresponding row in the R matrix.
//		Returns R_NilValue on failure.
/////////////////////////////////////////////////////////////////////////
SEXP do_matrix_sum_rows(SEXP A) {
	if (!isMatrix(A)) {
		Rprintf("do_matrix_sum_rows:  Oops, please pass argument ");
		Rprintf("as an R matrix:\n");
		Rprintf("\tdo_matrix_sum_rows(matrix)\n");
		return R_NilValue;
	} 
	//  check to make sure something was actually passed in...
	if (!A) {
		Rprintf("Oops, can't use an empty matrix in ");
		Rprintf("do_matrix_sum_rows...\n");
		return R_NilValue;
	}
	int protections=0, *dims;
	double* pA;

	PROTECT(A = AS_NUMERIC(A)); 
	protections++;
	pA = NUMERIC_POINTER(A);
	//  get the dimensions of the matrix
	dims = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!dims) {
		Rprintf("Oops, couldn't get the dimensions of the ");
		Rprintf("matrix in do_matrix_sum_rows(matrix)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	SEXP S = matrix_sum_rows(pA, dims);
	UNPROTECT(protections);
	return S;
} // do_matrix_sum_rows
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  matrix_sum_rows returns an array of doubles (size # rows in the
//    matrix) containing the sums of each row in the matrix.
//
//  Inputs:
//		A, a C array/matrix of type double
//		dims, an array of 2 ints containing the matrix dimensions.
//  Outputs:
//		S, an R vector, containing in each element the sum of the
//			corresponding row in the C matrix.
//		Returns R_NilValue on failure.
/////////////////////////////////////////////////////////////////////////
SEXP matrix_sum_rows(double* A, int dims[2]) {
	if (!A) {
		Rprintf("C code matrix_sum_rows:  Can't use NULL matrix!\n");
		return R_NilValue;
	}
	int protections = 0;
	//  return a vector/array of length # rows in A:
	SEXP S;
	PROTECT(S = allocVector(REALSXP, dims[0]));
	protections++;
	double* pS = NUMERIC_POINTER(S);
	if (!pS) {
		Rprintf("C code matrix_sum_rows:  Couldn't allocate");
		Rprintf("vector to return!\n");
		return R_NilValue;
	}
	int dims0 = dims[0], dims1 = dims[1];
	int disp  = dims0 * (dims1 - 1);
	//  beginning & end markers of matrix:
	double* pAR = &*A;
	double* pACend;
	double* ARend = &*pAR + dims0;
	//  step through items in rows in matrix:
	double* pAC;
	register double sum = 0.0;
	//  REMEMBER, R is COLUMN major!!!!!!!!!!!!!!!!!
	//  for each row in A:
	for (; pAR<ARend; pAR++, pS++) {
	//for (pAR; pAR<ARend; pAR++, pS++) {  -> 03/15/2023
		pACend = &*pAR + disp;
		//  for each column in Z:
		sum = 0.0;
		for (pAC=&*pAR; pAC<=pACend; pAC+=dims0) {
			sum += *pAC;
		}
		*pS = sum;
	}
	UNPROTECT(protections);
	return S;
}  // matrix_sum_rows
/////////////////////////////////////////////////////////////////////////
