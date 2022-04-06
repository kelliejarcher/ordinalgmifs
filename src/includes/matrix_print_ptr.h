/***********************************************************************
	 Author:  Johnny Layne
	   File:  matrix_print_ptr.h

	Purpose:  Provide a clean, fast matrix printing function for
				Dr. Kellie Archer's R code.

	Note:  I'm not sure I'm really using this any more, but sure is
		nice to have around, so I want to hang onto it.
 ***********************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP matrix_print_ptr(SEXP);

/////////////////////////////////////////////////////////////////////////
//  function definitions follow: 
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
SEXP matrix_print_ptr(SEXP A) { 
	int protections = 0;
	//  was matrix passed in?
	if (!isMatrix(A)) {
		Rprintf("matrix_print_ptr:  Oops, please pass 1st argument ");
		Rprintf("as an R matrix:  matrix_print_ptr(matrix)\n");
		return R_NilValue;
    }
	PROTECT(A = AS_NUMERIC(A));
	protections++;
	if (!A) {
		Rprintf("matrix_print_ptr:  Oops, please pass 1st argument ");
		Rprintf("as an R matrix:  matrix_print_ptr(matrix)\n");
		UNPROTECT(protections);
		return R_NilValue;
    }

	//  get the dimensions of the matrix
	int *dimsa;
	dimsa = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!dimsa) {
		Rprintf("Oops, couldn't get the dimensions of matrix A ");
		Rprintf("in matrix_print_ptr(matrix)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	//  pointer to matrix
	double* pA = NUMERIC_POINTER(A);
	//  beginning & end markers of matrices:
	double* pAR = pA;
	double* ARend = &*pAR + dimsa[0];
	//  iterators for rows in matrices:
	double*  pAC;
	double*  pACend;
	int dims0 = dimsa[0];
	int dims1 = dimsa[1];
	int disp = dims0 * (dims1 - 1);
	//  REMEMBER, R is COLUMN major!!!!!!!!!!!!!!!!!
	//  for each row in A:
	for (pAR; pAR<ARend; pAR++) {
		pACend = &*pAR + disp;
		Rprintf("|");
		//  for each column in A:
		for (pAC=&*pAR; pAC<=pACend; pAC+=dimsa[0]) {
			if (pAC < pACend) {
				Rprintf("%.4f, ", *pAC);
			}
			else {
				Rprintf("%.4f|\n", *pAC);
			}
		}
	}
	UNPROTECT(protections);
	return R_NilValue;
}
/////////////////////////////////////////////////////////////////////////
