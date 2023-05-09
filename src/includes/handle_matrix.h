/************************************************************************
	Author:  Johnny Layne
	  File:  handle_matrix.h

	Note:  I don't think that I actually use this for anything anymore,
		but as this code was pretty critical in testing & developing the
		other stuff, I think I'll keep it around.  It's a bunch of 
		utilities to print a C matrix, both to screen and a file, as
		well as play with column major vs. row major indexing.
 ************************************************************************/

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
void debug_print_Cmatrix(double*, int, int, int);
SEXP print_matrix(SEXP);
SEXP dum_print_matrix(SEXP, SEXP);
SEXP back_column(SEXP);
SEXP back_row(SEXP);
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  function definitions  //
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  debug_print_Cmatrix is just the usual version of stepping through
//    a C matrix.  It's purpose is only to help me develop a pointer-
//    math version (that's in another file now) that will perform the
//    same function, but hopefully be much faster.
//
//  Inputs:
//    double* mat, memory location of the first double value in the
//      matrix
//    int m, the # rows in the matrix
//    int n, the # columns in the matrix
//    int nrows, the # rows of the matrix to print, in case it's
//      desirable to print less rows than are in the matrix
//  Outputs:
//    the C matrix is printed to the screen, and the function
//    returns nothing.
/////////////////////////////////////////////////////////////////////////
void debug_print_Cmatrix(double* mat, int m, int n, int nrows) {
	int i, j, ndxA;
	i = j = ndxA = 0;

	//  step down the rows of the matrix
	for (i=0; i<m; i++) {
		//  step across the columns of the matrix
		for (j=0; j<n; j++) {
			//  get the proper index; R is column-major, C is row-major....
			if (i<nrows) {
				ndxA = j * m + i;
				Rprintf("%.4f\t", mat[ndxA]);
			}
		}
		Rprintf("\n");
	}
}  //  debug_print_Cmatrix
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  print_matrix(SEXP mat) is a bit of an improvement over the typical
//    means of stepping through and printing a C matrix.  Given a pointer
//    to the first double value in the matrix, and the dimensions of the
//    matrix, print_matrix does the indexing calculations itself to
//    figure out where to be in the matrix.  Still not as fast as using
//    pointer math to move through the array, but we're getting there!
//
//  Inputs:
//    SEXP mat, an R Matrix
//	Outputs:
//    The contents are printed to the screen, and the function returns 
//      R_NilValue.
//    On error though (non-Matrix passed in, memory allocation fails)
//      R_NilValue is returned.
/////////////////////////////////////////////////////////////////////////
SEXP print_matrix(SEXP mat) {
	if (!isMatrix(mat)) {
		Rprintf("Oops, please pass 1st argument as an R matrix:\n");
		Rprintf("\tprint_matrix(Matrix)\n");
		return R_NilValue;
	}

	int i, j, m, n, protections, *dims, ndxA;
	double *pA;
	i = j = protections = ndxA = 0;

	//  setup & check matrix for use
	PROTECT(mat = AS_NUMERIC(mat));
	protections++;
	//  check to make sure something was actually passed in...
	if (!mat) {
		Rprintf("Oops, can't use an empty matrix in print_matrix...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	//  get the dimensions of the matrix
	dims = INTEGER_POINTER(getAttrib(mat, R_DimSymbol));
	if (!dims) {
		Rprintf("Oops, couldn't get the dimensions of the matrix ");
		Rprintf("in print_matrix...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	m = dims[0];
	n = dims[1];

	//  pointer to matrix to allow moving through it easily
	pA = NUMERIC_POINTER(mat);
	//  check to make sure allocation worked...
	if (!pA) {
		Rprintf("Oops, couldn't allocate a matrix to work with in ");
		Rprintf("print_matrix...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  step down the rows of the matrix
	for (i=0; i<m; i++) {
		Rprintf("|");
		//  step across the columns of the matrix
		for (j=0; j<n; j++) {
			//  get the proper index; R is column-major, C is row-major....
			ndxA = j * m + i;
			Rprintf("%.4f\t", pA[ndxA]);
		}
		Rprintf("|\n");
	}

	UNPROTECT(protections);
	return R_NilValue;
}  //  print_matrix
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  pretty much the same thing as "print_matrix", above; lets the
//    user choose how many rows of the matrix to print though.
//  
//  Inputs:
//    SEXP mat, the R Matrix to print
//    SEXP nrows, the # of rows of mat to print
/////////////////////////////////////////////////////////////////////////
SEXP dum_print_matrix(SEXP mat, SEXP nrows) {
	if (!isMatrix(mat)) {
		Rprintf("Oops, please pass 1st argument as an R matrix:\n");
		Rprintf("\tdum_print_matrix(matrix)\n");
		return R_NilValue;
	}

	int i, j, m, n, nr, protections, *dims, ndxA;
	double *pA;
	i = j = protections = ndxA = 0;

	//  setup & check matrix for use
	PROTECT(mat = AS_NUMERIC(mat));
	protections++;
	//  check to make sure something was actually passed in...
	if (!mat) {
		Rprintf("Oops, can't use an empty matrix in dum_print_matrix...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  get the dimensions of the matrix
	dims = INTEGER_POINTER(getAttrib(mat, R_DimSymbol));
	if (!dims) {
		Rprintf("Oops, couldn't get the dimensions of the matrix ");
		Rprintf("in dum_print_matrix...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}
	m = dims[0];
	n = dims[1];

	//  get the limit on rows printed passed in by the user
	PROTECT(nrows = AS_INTEGER(nrows));
	protections++;
	nr = asInteger(nrows);

	//  pointer to matrix to allow moving through it easily
	pA = NUMERIC_POINTER(mat);
	//  check to make sure allocation worked...
	if (!pA) {
		Rprintf("Oops, couldn't allocate a matrix to work with in ");
		Rprintf("dum_print_matrix...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	//  step down the rows of the matrix
	for (i=0; i<m; i++) {
		//  step across the columns of the matrix
		for (j=0; j<n; j++) {
			//  get the proper index; R is column-major, C is row-major....
			if (i < nr) {
				ndxA = j * m + i;
				Rprintf("%.4f\t", pA[ndxA]);
			}
		}
		if (i < nr) {
			Rprintf("\n");
		}
	}

	UNPROTECT(protections);
	return R_NilValue;
}  //  dum_print_matrix
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  This is just a function to help me figure out how to do in C what
//    some R code is doing.
//  I needed to figure out how to use values in the PREVIOUS column for
//    some calculations...
/////////////////////////////////////////////////////////////////////////
SEXP back_column(SEXP mat) {
    if (!isMatrix(mat)) {
		Rprintf("back_column:  Oops, please pass argument ");
		Rprintf("as an R Matrix:  back_column(Matrix)\n");
		return R_NilValue;
    } 
    //  check to make sure something was actually passed in...
    if (!mat) {
        Rprintf("Oops, can't use an empty Matrix in back_column...\n");
        return R_NilValue;
    }
	int i=0, j=0, protections=0, *matdims, m, n;  //, ndx;
	PROTECT(mat = AS_NUMERIC(mat)); 
	protections++;
	matdims = INTEGER_POINTER(getAttrib(mat, R_DimSymbol));
    if (!matdims) {
        Rprintf("Oops, couldn't get the dimensions of the matrix ");
		Rprintf("in back_column...\n");
        UNPROTECT(protections);
        return R_NilValue;
    }
	m = matdims[0];
	n = matdims[1];
	Rprintf("matrix mat[%d][%d]...\n", m, n);
	//  step down the rows
	Rprintf("Printing matrix \"mat\", item INDEXES ");
	Rprintf("one column back marked by \"C\"...\n");
	for (i=0; i<m; i++) {
		Rprintf("i==%d) ", i);
		//  ...and across the columns
		for (j=0; j<n; j++) {
			if (0 < j) {
				//  j == C's column #,
				//  (j*m+i) == actual index in linear array,
				//  ((j-1)*m+i) == PREVIOUS column!
				Rprintf("(%d)[%d]C:%d ", j, (j*m+i), ((j-1)*m+i));
			}
			else {
				Rprintf("(%d)[%d]\t\t", j, (j*m+i));
			}
		}
		Rprintf("\n");
	}
	UNPROTECT(protections);
	return R_NilValue;
}  //  back_column
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  Just as with back_column, I wrote this as I needed to figure out
//    how to access some values stored in the previous row from the
//    current location in the matrix.
/////////////////////////////////////////////////////////////////////////
SEXP back_row(SEXP mat) {
    if (!isMatrix(mat)) {
		Rprintf("do_exp:  Oops, please pass 2nd argument ");
		Rprintf("as an R matrix:  back_row(matrix)\n");
		return R_NilValue;
    } 
    //  check to make sure something was actually passed in...
    if (!mat) {
        Rprintf("Oops, can't use an empty matrix in back_row...\n");
        return R_NilValue;
    }
	int i=0, j=0, protections=0, *matdims, m, n;
	PROTECT(mat = AS_NUMERIC(mat)); 
	protections++;
	matdims = INTEGER_POINTER(getAttrib(mat, R_DimSymbol));
    if (!matdims) {
        Rprintf("Oops, couldn't get the dimensions of the matrix ");
		Rprintf("in back_row...\n");
        UNPROTECT(protections);
        return R_NilValue;
    }
	m = matdims[0];
	n = matdims[1];
	Rprintf("matrix mat[%d][%d]...\n", m, n);
	//  step down the rows
	Rprintf("Printing matrix \"mat\", item INDEXES ");
	Rprintf("one row back marked by \"R\"...\n");
	for (i=0; i<m; i++) {
		//  ...and across the columns
		for (j=0; j<n; j++) {
			if (0 < i) {
				Rprintf("[%d]R:%d\t", (j*m+i), (j*m+(i-1)));
			}
			else {
				Rprintf("[%d]\t\t", (j*m+i));
			}
		}
		Rprintf("\n");
	}
	
	UNPROTECT(protections);
	return R_NilValue;
}  //  back_row
/////////////////////////////////////////////////////////////////////////
