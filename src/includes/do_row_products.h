
/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_row_products.h

	Purpose:  Provide a (hopefully!) faster version of part of 
				Dr. Kellie Archer's fn.cum function in her R code.
				This is to replace row sum on a matrix in R.

	Never forget, this is for COLUMN-MAJOR arrays!

 ***********************************************************************/

#ifndef boolean
typedef enum boolean { true=1, false=0 } bool;
#endif

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_row_products(SEXP);
SEXP row_products(double*, int[2]);

/////////////////////////////////////////////////////////////////////////
//  function definitions  //
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  setup to call row_products
//	Inputs:
//		A, an R matrix of doubles.
//  Outputs:
//		An R vector returned from row_products, containing the
//			row products of A.
//		Returns R_NilValue on failure.
/////////////////////////////////////////////////////////////////////////
SEXP do_row_products(SEXP A) {
	if (!isMatrix(A)) {
		Rprintf("C code do_row_products:  Oops, please pass argument ");
		Rprintf("as an R matrix:\n");
		Rprintf("\tdo_row_products(matrix)\n");
		return R_NilValue;
	} 
	//  check to make sure something was actually passed in...
	if (!A) {
		Rprintf("Oops, can't use an empty matrix in ");
		Rprintf("C code do_row_products...\n");
		return R_NilValue;
	}
	int protections=0, 
		*dims;
	double* pA;

	PROTECT(A = AS_NUMERIC(A)); 
	protections++;
	pA = NUMERIC_POINTER(A);
	//  get the dimensions of the matrix
	dims = INTEGER_POINTER(getAttrib(A, R_DimSymbol));
	if (!dims) {
		Rprintf("Oops, couldn't get the dimensions of the ");
		Rprintf("matrix in C code do_row_products(matrix)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	SEXP S = row_products(pA, dims);
	UNPROTECT(protections);
	return S;
} // do_row_products
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  0's in rows are ignored!!!
//  row_products returns the results of calculating the products of
//    the items in a row of matrix.
//
//  Inputs:
//		double array A, the matrix passed in from R
//		2d int array dims, containing the # rows & # cols in A
//  Outputs:
//		An R vector containing the results of the row products 
//		  of A.
//		Returns R_NilValue on failure.
/////////////////////////////////////////////////////////////////////////
SEXP row_products(double* A, int dims[2]) {
	if (!A) {
		Rprintf("C code row_products:  Can't use NULL matrix!\n");
		return R_NilValue;
	}
	int protections = 0;
	//  return a vector/array of length # rows in A:
	SEXP S;
	PROTECT(S = allocVector(REALSXP, dims[0]));
	protections++;
	double* pS = NUMERIC_POINTER(S);
	if (!pS) {
		Rprintf("C code row_products:  Couldn't allocate");
		Rprintf("vector to return!\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	int dims0 = dims[0], dims1 = dims[1];
	int disp  = dims0 * (dims1 - 1);
	//  beginning & end markers of matrix:
	double* pAR = &*A;
	double* pACend;
	double* ARend = &*pAR + dims0;
	//double* pS = &*S;
	//  iterators for rows in matrix:
	double* pAC;
	register double p = 0.0;
	//  REMEMBER, R is COLUMN major!!!!!!!!!!!!!!!!!
	//  for each row in A:
	for (; pAR<ARend; pAR++, pS++) {
		pACend = &*pAR + disp;
		//  for each column in A (step across this row):
		pAC = &*pAR;
		//  p is used to ignore 0's in row products
		p = *pAC;
		for (pAC=pAC+dims0; pAC<=pACend; pAC+=dims0) {
			if (*pAC > 0.0) {
				if (p > 0.0) {
					p *= *pAC;
				}
				else {
					p = *pAC;
				}
			}
		}
		*pS = p;
		p = 0.0;
	}
	UNPROTECT(protections);
	return S;
}  // row_products
/////////////////////////////////////////////////////////////////////////
