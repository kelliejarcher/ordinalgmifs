/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_cumsum.h

	Purpose:  Provide a (hopefully!) faster version of part  of
				Dr. Kellie Archer's acat.fn function in her R code.

 ***********************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_row_cumsum(SEXP);
SEXP row_cumsum(double*, int[2]);
SEXP do_col_cumsum(SEXP);
SEXP col_cumsum(double*, int[2]);

/////////////////////////////////////////////////////////////////////////
//  function definitions follow: 
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/*
	do_row_cumsum accepts 1 R Matrix as argument, sets up
	  for a call to row_cumsum, and returns the R Matrix returned
	  by row_cumsum.
	Returns R_NilValue on failure.
*/
/////////////////////////////////////////////////////////////////////////
SEXP do_row_cumsum(SEXP A) {
	int protections = 0;
	//  were matrices passed in?
	if (!isMatrix(A)) {
		Rprintf("do_cumsum:  Oops, please pass 1st argument ");
		Rprintf("as an R matrix:  do_cumsum(matrix)\n");
		return R_NilValue;
    }
	PROTECT(A = AS_NUMERIC(A));
	protections++;

	//  check to make sure no NULL data passed in...
	if (!A) {
		Rprintf("Oops, can't use an empty matrix in do_cumsum...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  get the dimensions of the matrices
	int *adims;
	adims = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!adims) {
		Rprintf("Oops, couldn't get the dimensions of matrix A ");
		Rprintf("in do_row_cumsum(matrix A, matrix B)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  pointers to matrices:
	double* pA = NUMERIC_POINTER(A);
	SEXP S = row_cumsum(pA, adims);
	UNPROTECT(protections);
	return S;
}  // do_row_cumsum
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  Return a double matrix containing the result of R's 
//      apply(A, 1, cumsum),
//	  over the ROWs of A...
//	Returned matrix has reversed dimensions of A, so if A is MxN,
//	  return an NxM matrix...
//	Returns R_NilValue on failure.
/////////////////////////////////////////////////////////////////////////
SEXP row_cumsum(double* A, int dims[2]) {
	int protections = 0;
	int dimsa0=dims[0];
	int dimsb0=dims[1], dimsb1=dims[0];
	//  return an R matrix:
	SEXP B;
	PROTECT(B = allocMatrix(REALSXP, dimsb0, dimsb1));
	protections++;
	double* pB = NUMERIC_POINTER(B);
	if (!pB) {
		Rprintf("Oops, couldn't allocate a matrix to return from C ");
		Rprintf("code row_cumsum...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  using a lot of pointer math to speed things up; it
	//    helps me to think of things here in "blocks" like
	//    bricks in a wall, chunks of memory at a time used
	//    to step through an array until the last chunk is
	//    reached.
	//  beginning & end markers of matrices:
	double* pAR = &*A;
	double* ARend = &A[dimsa0];
	//  step through items in matrices:
	double* pARi;
	double* pBC = &*pB;
	double* pBCi;
	//  how big a memory chunk to step
	int dispb  = dimsb0;
	//  fast fast fast!  If you get it...
	register double sum;

	//  for each row in A:
	for (; pAR<ARend; pAR++, pBC+=dispb) {
		sum = 0.0;
		pARi = &*pAR;
		//  for each item in this row of A, and this
		//  column of B:
		for (pBCi=&*pBC; pBCi<(pBC+dispb); pBCi++) {
			sum += *pARi;
			*pBCi = sum;
			pARi += dimsa0;
		}
	}
	UNPROTECT(protections);
	return B;
}  // row_cumsum
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/*
	do_col_cumsum accepts 1 R Matrix as argument, sets up
	  for a call to col_cumsum, and returns the R Matrix
	  returned by col_cumsum.
	Returns R_NilValue on failure.
*/
/////////////////////////////////////////////////////////////////////////
SEXP do_col_cumsum(SEXP A) {
	int protections = 0;
	//  were matrices passed in?
	if (!isMatrix(A)) {
		Rprintf("do_cumsum:  Oops, please pass 1st argument ");
		Rprintf("as an R matrix:  do_cumsum(matrix)\n");
		return R_NilValue;
    }
	PROTECT(A = AS_NUMERIC(A));
	protections++;

	//  check to make sure no NULL data passed in...
	if (!A) {
		Rprintf("Oops, can't use an empty matrix in do_cumsum...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  get the dimensions of the matrices
	int *adims;
	adims = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!adims) {
		Rprintf("Oops, couldn't get the dimensions of matrix A ");
		Rprintf("in do_row_cumsum(matrix A, matrix B)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  pointers to matrices:
	double* pA = NUMERIC_POINTER(A);
	SEXP S = col_cumsum(pA, adims);
	UNPROTECT(protections);
	return S;
}  // do_row_cumsum
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  Return a double matrix containing the result of R's 
//      apply(A, 2, cumsum),
//	  over the COLUMNs of A...
//	Returned matrix has same dimensions of A, so if A is MxN,
//	  return an MxN matrix...
//  Returns R_NilValue on failure.
/////////////////////////////////////////////////////////////////////////
SEXP col_cumsum(double* A, int dims[2]) {
	int protections = 0;
	//  return an R matrix:
	SEXP B;
	PROTECT(B = allocMatrix(REALSXP, dims[0], dims[1]));
	protections++;
	double* pB = NUMERIC_POINTER(B);
	if (!pB) {
		Rprintf("Oops, couldn't allocate a matrix to return from C ");
		Rprintf("code col_cumsum...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  beginning & end markers of matrices:
	double* pAC = &*A;
	double* ACend = &A[dims[0] * dims[1] - dims[0]];
	//  step through items in matrices:
	double* pACi;
	double* pBC = &*pB;
	double* pBCi;
	//  easier to use dimensions of matrices
	int dimsa0 = dims[0];
	int dimsb0 = dims[0];
	//  size of memory chunks to step by
	int dispa  = dimsa0;
	int dispb  = dimsb0;
	register double sum;

	//  for each column in A & B:
	for (; pAC<=ACend; pAC+=dispa, pBC+=dispb) {
		sum = 0.0;
		pACi = &*pAC;
		pBCi = &*pBC;
		//  for each item in this column of A, and this
		//  column of B:
		for (; pACi<(pAC+dimsa0); pACi++, pBCi++) {
			sum += *pACi;
			*pBCi = sum;
		}
	}
	UNPROTECT(protections);
	return B;
}  // col_cumsum
/////////////////////////////////////////////////////////////////////////
