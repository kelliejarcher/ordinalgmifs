/***********************************************************************
	 Author:  Johnny Layne
	   File:  do_exp.h

	Purpose:  Provide a (hopefully!) faster version of part of 
				Dr. Kellie Archer's fn.cum function in her R code.
 ***********************************************************************/

#ifndef boolean
typedef enum boolean { true=1, false=0 } bool;
#define boolean
#endif

/////////////////////////////////////////////////////////////////////////
//  function prototypes
/////////////////////////////////////////////////////////////////////////
SEXP do_exp(SEXP, SEXP, SEXP);
SEXP exp_ptr(int, int[2], double*, double*);

/////////////////////////////////////////////////////////////////////////
//  function definitions  //
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  setup to call exp_ptr
//  Inputs:
//		k, an integer containing # columns in pizmat (need k?)
//			I think that "k" is a leftover from the time when 
//			I was just getting started on this code...
//		Z, a matrix of doubles
//		pizmat, another matrix of doubles
//  Outputs:
//		pizmat, with values changed according to Kellie Archer's
//			original R code version.
//		R_NilValue upon failure.
/////////////////////////////////////////////////////////////////////////
SEXP do_exp(SEXP k, SEXP Z, SEXP pizmat) {
	if (!isMatrix(Z)) {
		Rprintf("do_exp:  Oops, please pass 2nd argument ");
		Rprintf("as an R matrix:\n");
		Rprintf("\tdo_exp(integer, matrix, matrix)\n");
		return R_NilValue;
	} 
	//  check to make sure something was actually passed in...
	if (!Z) {
		Rprintf("Oops, can't use an empty matrix in do_exp...\n");
		return R_NilValue;
	}
	//  do same checks for pizmat:
	if (!isMatrix(pizmat)) {
		Rprintf("do_exp:  Oops, please pass 3rd argument ");
		Rprintf("as an R matrix:\n");
		Rprintf("\tdo_exp(integer, matrix, matrix)\n");
		return R_NilValue;
	} 
	//  check to make sure something was actually passed in...
	if (!pizmat) {
		Rprintf("Oops, can't use an empty matrix in do_exp...\n");
		return R_NilValue;
	}

	int len=0, 
		protections=0, 
		*zdims;
	double* pZ;
	double* piz;
	
	PROTECT(k = AS_INTEGER(k));
	protections++;
	len = asInteger(k);

	PROTECT(Z = AS_NUMERIC(Z)); 
	protections++;
	pZ = NUMERIC_POINTER(Z);
	//  get the dimensions of the matrix
	zdims = INTEGER_POINTER(getAttrib(Z, R_DimSymbol));
	if (!zdims) {
		Rprintf("Oops, couldn't get the dimensions of the 1st ");
		Rprintf("matrix in do_exp(integer, matrix, matrix)...\n");
		UNPROTECT(protections);
		return R_NilValue;
	}

	PROTECT(pizmat = AS_NUMERIC(pizmat)); 
	protections++;
	piz = NUMERIC_POINTER(pizmat);
	exp_ptr(len, zdims, pZ, piz);
	UNPROTECT(protections);
	return pizmat;
} // do_exp
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
//  exp_ptr is C code to replace some R code in Kellie Archer's 
//    "CumultCommented.R" that does a lot of calls to "exp"
//    (I will leave it to the reader to look up that function)
//    inside of a big slow loop.  This loop was one of the big
//    bottlenecks (timewise) in this R code, so this C code 
//    provides a much faster implementation of her R code.
//
//  "k" is probably redundant, since I can just get the 
//    dimensions from R-to-C functions and dimension of
//    "p" is one column more than in Z.
//
//  Inputs:
//		k, # columns in p
//		dims, 2-element array of ints containing rows, cols in Z.
//		Z, a matrix of doubles containing values to be used to
//			fill p.
//		p, a matrix of doubles to be returned filled with the
//			results of using the exp function.
//  Outputs:
//		p, a matrix of doubles to be returned filled with the
//			results of using the exp function.
//		R_NilValue upon failure.
/////////////////////////////////////////////////////////////////////////
SEXP exp_ptr(int k, int dims[2], double* Z, double* p) {
	if (!Z || !p) {
		Rprintf("C code matrix_exp_ptr:  Can't use NULL matrix!\n");
		return R_NilValue;
	}

	int dims0 = dims[0], dims1 = dims[1];

	//   fastfastfast!  If compiler/system let you have it...
	register int j = 0;
	register double expZC, expZC_back;
	int disp  = dims0 * (dims1 - 1);
    //  using a lot of pointer math to speed things up; it
    //    helps me to think of things here in "blocks" like
    //    bricks in a wall, chunks of memory at a time used
    //    to step through an array until the last chunk is
    //    reached.
	//  beginning & end markers of matrices:
	double* pZR = &*Z;
	double* pZCend;
	double* ZRend = &*pZR + dims0;
	double* pPR = &*p;
	//  increment through rows in matrices:
	double* pZC;
	double* pZC_back;
	double* pPC;
	//  REMEMBER, R is COLUMN major!!!!!!!!!!!!!!!!!
	//  for each row in Z:
	for (; pZR<ZRend; pZR++, pPR++) {
		pZCend = &*pZR + disp;
		//  for each column in Z:
		j = 0;
		pPC = &*pPR;
		for (pZC=&*pZR; pZC<=pZCend+dims0; pZC+=dims0) {
			if (0 == j) {
				//  reduce calls to exp:
				expZC = exp(*pZC);
				*pPC =
					expZC / (1.0 + expZC);
			}
			else if (j == k - 1) {
				pZC_back = pZC - dims0;
				//  reduce calls to exp:
				expZC_back = exp(*pZC_back);
				*pPC =
					 1.0 - expZC_back /
					(1.0 + expZC_back);
			}
			else {
				pZC_back = pZC - dims0;
				//  reduce calls to exp:
				expZC = exp(*pZC);
				expZC_back = exp(*pZC_back);
				*pPC =
					expZC / (1.0 + expZC) -
					expZC_back / (1.0 + expZC_back);
			}
			j++;
			pPC += dims0;
		}
	}
	return R_NilValue;
}  // exp_ptr
/////////////////////////////////////////////////////////////////////////
