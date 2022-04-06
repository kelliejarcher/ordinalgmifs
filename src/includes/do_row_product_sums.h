/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_row_product_sums.h

	Purpose:  Provide a (hopefully!) faster version of part  of
				Dr. Kellie Archer's fn.cum function in her R code.

 ***********************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_row_product_sums(SEXP, SEXP);
SEXP row_product_sums(double*, double*, int[2]);

/////////////////////////////////////////////////////////////////////////
//  function definitions follow: 
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/*
	do_row_product_sums accepts 2 R Matrices as arguments, sets up
	for a call to row_product_sums, and returns the R Vector returned
	by row_product_sums containing the sums of the rows of the result
	of multiplying the two matrices.

	On failure, R_NilValue is returned.
*/
/////////////////////////////////////////////////////////////////////////
SEXP do_row_product_sums(SEXP A, SEXP B) {
	int protections = 0;
	//  were R matrices passed in?
	if (!isMatrix(A)) {
		Rprintf("do_row_product_sums:  Oops, please pass 1st argument ");
		Rprintf("as an R matrix:  do_row_product_sums(matrix, matrix)\n");
		return R_NilValue;
	}
	PROTECT(A = AS_NUMERIC(A));
	protections++;

	if (!isMatrix(B)) {
		Rprintf("do_row_product_sums:  Oops, please pass 2nd argument ");
		Rprintf("as an R matrix:  do_row_product_sums(matrix, matrix)\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	PROTECT(B = AS_NUMERIC(B));
	protections++;

	//  check to make sure no NULL data passed in...
	if (!A || !B) {
		Rprintf("Oops, can't use an empty matrix in ");
		Rprintf("do_row_product_sums...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  get the dimensions of the matrices
	int *adims, *bdims;
	adims = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!adims) {
		Rprintf("Oops, couldn't get the dimensions of matrix A ");
		Rprintf("in do_row_product_sums(matrix A, matrix B)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	bdims = INTEGER_POINTER(getAttrib(B, R_DimSymbol));
	if (!bdims) {
		Rprintf("Oops, couldn't get the dimensions of matrix B ");
		Rprintf("in do_row_product_sums(matrix A, matrix B)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	if ((adims[0] != bdims[0]) ||
		(adims[1] != bdims[1])) {
		Rprintf("C code do_row_product_sums:  dimensions of ");
		Rprintf("both matrices must be the same...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  pointers to matrices:
	double* pA = NUMERIC_POINTER(A);
	double* pB = NUMERIC_POINTER(B);
	SEXP S = row_product_sums(pA, pB, adims);
	UNPROTECT(protections);
	return S;
}  // do_row_product_sums
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  return a double matrix containing the result of sum(A * B),
//	  over the ROWs of the result of matrix multiplication, A * B...
//  Returns R_NilValue upon failure.
/////////////////////////////////////////////////////////////////////////
SEXP row_product_sums(double* A, double* B, int dims[2]) {
	double* C = (double*)malloc(dims[0] * sizeof(double));
	if (!C) {
		Rprintf("row_product_sum:  Couldn't allocate C matrix...\n");
		return NULL;
	}
	int protections = 0;
	//  return an R vector:
	SEXP S;
	PROTECT(S = allocVector(REALSXP, dims[0]));
	protections++;
	double* pS = NUMERIC_POINTER(S);
	if (!pS) {
		Rprintf("Oops, couldn't allocate a vector to return from C ");
		Rprintf("code row_product_sums...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  beginning & end markers of matrices:
	double* pAR = &*A;
	double* ARend = &A[dims[0]];
	double* pBR = &*B;
	//  step through items in matrices:
	double* pAC;
	double* pACend;
	double* pBC;
	//  make access to this fastfastfast!
	register double sum;
	//  make it a bit easer using array
	//    dimensions
	int dims0=dims[0], dims1=dims[1];
	//  precalculate how much to step through
	//    array
	int disp = dims0 * (dims1 - 1);

	//  for each row in A & B:
	for (; pAR<ARend; pAR++, pBR++, pS++) {
		sum = 0.0;
		pACend = &*pAR + disp;
		//  for each column in A & B:
		for (pAC=&*pAR, pBC=&*pBR; pAC<=pACend; pAC+=dims0, pBC+=dims0) {
			sum += *pAC * *pBC;
		}
		*pS = sum;
	}
	UNPROTECT(protections);
	return S;
}  // row_product_sums
/////////////////////////////////////////////////////////////////////////
