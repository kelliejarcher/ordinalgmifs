/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_matrix_mult.h

	Purpose:  Provide a (hopefully!) faster version of 
				matrix multiplication for Dr. Kellie
				Archer's R code...

 ***********************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_matrix_mult(SEXP, SEXP);
SEXP mult_matrix_ptr(double*, int[2], double*, int[2]);

/////////////////////////////////////////////////////////////////////////
//  function definitions follow: 
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/*
	do_matrix_mult accepts 2 R Matrices as arguments, sets up
	  for a call to matrix_mult_ptr, and returns the R Matrix
	  returned by matrix_mult_ptr.
	Returns R_NilValue upon failure.
*/
/////////////////////////////////////////////////////////////////////////
SEXP do_matrix_mult(SEXP A, SEXP B) {
	int protections = 0;
	//  were matrices passed in?
	if (!isMatrix(A)) {
		Rprintf("do_matrix_mult_sums:  Oops, please pass 1st argument ");
		Rprintf("as an R matrix:  do_matrix_mult_sums(matrix, matrix)\n");
		return R_NilValue;
	}
	PROTECT(A = AS_NUMERIC(A));
	protections++;

	if (!isMatrix(B)) {
		Rprintf("do_matrix_mult_sums:  Oops, please pass 2nd argument ");
		Rprintf("as an R matrix:  do_matrix_mult_sums(matrix, matrix)\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	PROTECT(B = AS_NUMERIC(B));
	protections++;

	//  check to make sure no NULL data passed in...
	if (!A || !B) {
		Rprintf("Oops, can't use an empty matrix in ");
		Rprintf("do_matrix_mult_sums...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  get the dimensions of the matrices
	int *adims, *bdims;
	adims = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!adims) {
		Rprintf("Oops, couldn't get the dimensions of matrix A ");
		Rprintf("in do_matrix_mult_sums(matrix A, matrix B)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	bdims = INTEGER_POINTER(getAttrib(B, R_DimSymbol));
	if (!bdims) {
		Rprintf("Oops, couldn't get the dimensions of matrix B ");
		Rprintf("in do_matrix_mult_sums(matrix A, matrix B)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  pointers to matrices:
	double* pA = NUMERIC_POINTER(A);
	double* pB = NUMERIC_POINTER(B);
	SEXP M = mult_matrix_ptr(pA, adims, pB, bdims);
	UNPROTECT(protections);
	return M;
}  // do_matrix_mult
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  return an R Matrix containing the result of (A %*% B),
//    R's matrix multiplication...
//  Returns R_NilValue upon failure.
/////////////////////////////////////////////////////////////////////////
SEXP mult_matrix_ptr(double* A, int dimsa[2], 
					 double* B, int dimsb[2]) {
	if ((0 >= dimsa[0]) || (0 >= dimsa[1]) ||
		(0 >= dimsb[0]) || (0 >= dimsb[1])) {
		Rprintf("C code matrix_mult_ptr:  Sorry, no ");
		Rprintf("dimensions <= 0 for matrices!\n");
		return R_NilValue;
	}   
	if (dimsa[1] != dimsb[0]) {
		Rprintf("C code matrix_mult_ptr:  # columns in ");
		Rprintf("left matrix must == # rows in right ");
		Rprintf("matrix...\n");
		return R_NilValue;
	}   
	if (!A || !B) {
		Rprintf("C code matrix_mult_ptr:  No NULL ");
		Rprintf("matrices!\n");
		return R_NilValue;
	}   

	//  make it a bit easier & more readable
	//  to handle array dimensions later on:
	int dimsa0 = dimsa[0]; //, dimsa1 = dimsa[1]; -> 05/03/2023
	int dimsb0 = dimsb[0], dimsb1 = dimsb[1];
	//  I use "displacements" as a means to
	//  step through the arrays in the proper
	//  manner according to the memory chunks
	//  that are needed to be stepped at a time.
	//  So I have to figure out, for instance,
	//  how many "double addresses" to step
	//  to advance 1 column, and so on...
	//int dispa  = dimsa0 * (dimsa1 - 1);  ->  05/03/2023
	int dispb  = dimsb0;
	int dispm  = dimsa0;
	/*  
	In order to multiply two matrices, the matrix on the left (A) must 
		have as many columns as the matrix on the right (B) has rows. 
		That way you can match up each pair while you're multiplying.
		The size of the final matrix is determined by the rows in the 
		left matrix (A) and the columns in the right (B).
	*/
	int protections = 0;
	SEXP M;
	PROTECT(M = allocMatrix(REALSXP, dimsa0, dimsb1));
	protections++;
	double* pM = NUMERIC_POINTER(M);
	if (!pM || !M) {
		Rprintf("Oops, couldn't allocate a matrix to return from C ");
		Rprintf("code mult_matrix_ptr...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	double* pMR = &*pM;
	double* pMCi;
	
	//  beginning & end markers of matrices:
	double* pAR = &*A;
	double* ARend = &*pAR + dimsa[0];
	double* pBR = &*B;
	//double* BRend = &*pBR + dimsb[0];  --> 03/15/2023
	double* pBCend = &B[dimsb[0] * dimsb[1] - dimsb[0]];
	//double* pACend; -> 05/01/2023
	//  incrementers for items in matrices:
	//double* pAC;  -> 03/15/2023
	double* pACi;
	double* pBC;
	double* pBCi;
	//  used a LOT, make it fastfastfast if we can
	register double sum;
	//  REMEMBER, R is COLUMN major!!!!!!!!!!!!!!!!!
	//  for each row in A:
	for (; pAR<ARend; pAR++, pMR++) {
	//for (pAR; pAR<ARend; pAR++, pMR++) {  -> 03/15/2023
		//pACend = &*pAR + dispa;  -> 05/01/2023
		pMCi = &*pMR;
		//  for each column in B:
		for (pBC=&*pBR; pBC<=pBCend; pBC+=dispb) {
			//  for each item in this row of A, and this
			//  column of B:
			pACi = &*pAR;
			sum = 0.0;
			for (pBCi=&*pBC; pBCi<(pBC+dispb); pBCi++) {
				sum += *pACi * *pBCi;
				pACi += dimsa0;
			}
			*pMCi = sum;
			pMCi += dispm;
		}
	}
	UNPROTECT(protections);
	return M;
}  // matrix_mult_ptr
/////////////////////////////////////////////////////////////////////////
