/*
 * $Id$
 * John Ashburner
 */

/*
 * This code is based on that of Philippe Thevenaz, which I took from:
 *	http://bigwww.epfl.ch/algorithms.html
 *
 * It has been substantially modified, so blame me (John Ashburner) if there
 * are any bugs. Many thanks to Philippe Thevenaz for advice with the code.
 *
 * See:
 *	M. Unser.
 *	"Splines: A Perfect Fit for Signal and Image Processing,"
 *	IEEE Signal Processing Magazine, 16(6):22-38 (1999)
 *
 *	P. Thevenaz and T. Blu and M. Unser.
 *	"Interpolation Revisited"
 *	IEEE Transactions on Medical Imaging 19(7):739-758 (2000).
 */


#include <math.h>
#include "mex.h"

/***************************************************************************************
Different degrees of B-splines
	x - position relative to origin
	returns value of basis function at x
*/

/*static double wt1(double x)
{
	x = fabs(x);
	return((x > 1.0) ? (0.0) : (1.0 - x));
}*/

static double wt2(double x)
{
	x = fabs(x);
	if (x < 0.5)
		return(0.75 - x*x);
	if (x < 1.5)
	{
		x = 1.5 - x;
		return(0.5*x*x);
	}
	return(0.0);
}

static double wt3(double x)
{
	x = fabs(x);
	if (x < 1.0)
		return(x*x*(x - 2.0)*0.5 + 2.0/3.0);
	if (x < 2.0)
	{
		x = 2.0 - x;
		return(x*x*x*(1.0/6.0));
	}
	return(0.0);
}

static double wt4(double x)
{
	x = fabs(x);
	if (x < 0.5)
	{
		x *= x;
		return(x*(x*0.25 - 0.625) + 115.0/192.0);
	}
	if (x < 1.5)
		return(x*(x*(x*(5.0/6.0 - x*(1.0/6.0)) - 1.25) + 5.0/24.0) + 55.0/96.0);
	if (x < 2.5)
	{
		x -= 2.5;
		x *= x;
		return(x*x*(1.0/24.0));
	}
	return(0.0);
}

static double wt5(double x)
{
	double y;
	x = fabs(x);
	if (x < 1.0)
	{
		y = x*x;
		return(y*(y*(0.25 - x*(1.0/12.0)) - 0.5) + 0.55);
	}
	if (x < 2.0)
		return(x*(x*(x*(x*(x*(1.0/24.0) - 0.375) + 1.25) - 1.75) + 0.625) + 0.425);
	if (x < 3.0)
	{
		y = 3.0 - x;
		x = y*y;
		return(y*x*x*(1.0/120.0));
	}
	return(0.0);
}

static double wt6(double x)
{
	x = fabs(x);
	if (x < 0.5)
	{
		x *= x;
		return(x*(x*(7.0/48.0 - x*(1.0/36.0)) - 77.0/192.0) + 5887.0/11520.0);
	}
	if (x < 1.5)
		return(x*(x*(x*(x*(x*(x*(1.0/48.0) - 7.0/48.0) + 0.328125)
			 - 35.0/288.0) - 91.0/256.0) - 7.0/768.0) + 7861.0/15360.0);
	if (x < 2.5)
		return(x*(x*(x*(x*(x*(7.0/60.0 - x*(1.0/120.0)) - 0.65625)
			+ 133.0/72.0) - 2.5703125) + 1267.0/960.0) + 1379.0/7680.0);
	if (x < 3.5)
	{
		x -= 3.5;
		x *= x*x;
		return(x*x*(1.0/720.0));
	}
	return(0.0);
}

static double wt7(double x)
{
	double y;

	x = fabs(x);
	if (x < 1.0)
	{
		y = x*x;
		return(y*(y*(y*(x*(1.0/144.0) - 1.0/36.0) + 1.0/9.0) - 1.0/3.0)
			+ 151.0/315.0);
	}
	if (x < 2.0)
		return(x*(x*(x*(x*(x*(x*(0.05 - x*(1.0/240.0)) - 7.0/30.0) + 0.5)
			- 7.0/18.0) - 0.1) - 7.0/90.0) + 103.0/210.0);
	if (x < 3.0)
		return(x*(x*(x*(x*(x*(x*(x*(1.0/720.0) - 1.0/36.0) + 7.0/30.0)
			- 19.0/18.0) + 49.0/18.0) - 23.0/6.0) + 217.0/90.0) - 139.0/630.0);
	if (x < 4.0)
	{
		y = 4.0 - x;
		x = y*y*y;
		return(x*x*y*(1.0/5040.0));
	}
	return(0.0);
}

/***************************************************************************************
Derivatives of different degrees of B-splines
	x - position relative to origin
	returns derivative of basis function at x
*/

static double dwt2(double x)
{
	int s;
	s = (x>0 ? 1 : -1);
	x = fabs(x);

	if (x < 0.5)
		return(-2*x*s);
	if (x < 1.5)
		return((x - 1.5)*s);
	return(0.0);
}

static double dwt3(double x)
{
	int s;
	s = (x>0 ? 1 : -1);
	x = fabs(x);


	if (x < 1.0)
		return(x*(1.5*x - 2.0)*s);
	if (x < 2.0)
	{
		x = x - 2.0;
		return(-0.5*x*x*s);
	}
	return(0.0);
}

static double dwt4(double x)
{
	int s;
	s = (x>0 ? 1 : -1);
	x = fabs(x);

	if (x < 0.5)
	{
		return((x*(x*x - 5.0/4.0))*s);
	}
	if (x < 1.5)
		return((x*(x*(x*(-2.0/3.0) + 2.5) - 5.0/2.0) + 5.0/24.0)*s);
	if (x < 2.5)
	{
		x = x*2.0 - 5.0;
		return((1.0/48.0)*x*x*x*s);
	}
	return(0.0);
}

static double dwt5(double x)
{
	int s;
	s = (x>0 ? 1 : -1);
	x = fabs(x);

	if (x < 1.0)
		return((x*(x*(x*(x*(-5.0/12.0) + 1.0)) - 1.0))*s);
	if (x < 2.0)
		return((x*(x*(x*(x*(5.0/24.0) - 1.5) + 3.75) - 3.5) + 0.625)*s);
	if (x < 3.0)
	{
		x -= 3.0;
		x *= x;
		return((-1.0/24.0)*x*x*s);
	}
	return(0.0);
}

static double dwt6(double x)
{
	double y;
	int s;
	s = (x>0 ? 1 : -1);
	x = fabs(x);

	if (x < 0.5)
	{
		y = x*x;
		return(x*((7.0/12.0)*y - (1.0/6.0)*y*y - (77.0/96.0))*s);
	}
	if (x < 1.5)
		return((x*(x*(x*(x*(x*0.125 - 35.0/48.0) + 1.3125) - 35.0/96.0)
			- 0.7109375) - 7.0/768.0)*s);
	if (x < 2.5)
		return((x*(x*(x*(x*(x*(-1.0/20.0) + 7.0/12.0) - 2.625) + 133.0/24.0)
			- 5.140625) + 1267.0/960.0)*s);
	if (x < 3.5)
	{
		x *= 2.0;
		x -= 7.0;
		y = x*x;
		return((1.0/3840.0)*y*y*x*s);
	}
	return(0.0);
}

static double dwt7(double x)
{
	double y;
	int s;
	s = (x>0 ? 1 : -1);
	x = fabs(x);

	if (x < 1.0)
	{
		y = x*x;
		return(x*(y*(y*(x*(7.0/144.0) - 1.0/6.0) + 4.0/9.0) - 2.0/3.0)*s);
	}
	if (x < 2.0)
		return((x*(x*(x*(x*(x*(x*(-7.0/240.0) + 3.0/10.0)
			- 7.0/6.0) + 2.0) - 7.0/6.0) - 1.0/5.0) - 7.0/90.0)*s);
	if (x < 3.0)
		return((x*(x*(x*(x*(x*(x*(7.0/720.0) - 1.0/6.0)
			+ 7.0/6.0) -38.0/9.0) + 49.0/6.0) - 23.0/3.0) + 217.0/90.0)*s);
	if (x < 4.0)
	{
		x -= 4;
		x *= x*x;
		x *= x;
		return((-1.0/720.0)*x*s);
	}
	return(0.0);
}

/***************************************************************************************
Generate B-spline basis functions
	d	- degree of spline
	x	- position relative to centre
	i	- pointer to first voxel position in convolution
	w	- vector of spline values

	Should really combine this function with wt2 to wt7 for most
	efficiency (as for case 0).

	Note that 0th degree B-spline returns nearest neighbour basis.
*/
static void weights(int d, double x, int *i, double w[])
{
	int k;

	*i = floor(x-(d-1)*0.5);
	x -= *i;

	switch (d){
	case 2:
		for(k=0; k<=2; k++) w[k] = wt2(x-k);
		break;
	case 3:
		for(k=0; k<=3; k++) w[k] = wt3(x-k);
		break;
	case 4:
		for(k=0; k<=4; k++) w[k] = wt4(x-k);
		break;
	case 5:
		for(k=0; k<=5; k++) w[k] = wt5(x-k);
		break;
	case 6:
		for(k=0; k<=6; k++) w[k] = wt6(x-k);
		break;
	case 7:
		for(k=0; k<=7; k++) w[k] = wt7(x-k);
		break;

	case 1:
		w[0] = 1.0-x;
		w[1] = x;
		break;
	case 0:
		w[0] = 1.0; /* Not correct at discontinuities */
		break;

	default:
		for(k=0; k<=7; k++) w[k] = wt7(x-k);
	}
}


/***************************************************************************************
Generate derivatives of B-spline basis functions
	d	- degree of spline
	x	- position relative to centre
	i	- pointer to first voxel position in convolution
	w	- vector of spline values

	Should really combine this function with dwt2 to dwt7 for most
	efficiency (as for case 0 and case 1).

	Note that 0th and 1st degree B-spline return derivatives of
	nearest neighbour and linear interpolation bases.
*/
static void dweights(int d, double x, int *i, double w[])
{
	int k;
	*i = floor(x-(d-1)*0.5);
	x -= *i;

	switch (d){
	case 2:
		for(k=0; k<=2; k++) w[k] = dwt2(x-k);
		break;
	case 3:
		for(k=0; k<=3; k++) w[k] = dwt3(x-k);
		break;
	case 4:
		for(k=0; k<=4; k++) w[k] = dwt4(x-k);
		break;
	case 5:
		for(k=0; k<=5; k++) w[k] = dwt5(x-k);
		break;
	case 6:
		for(k=0; k<=6; k++) w[k] = dwt6(x-k);
		break;
	case 7:
		for(k=0; k<=7; k++) w[k] = dwt7(x-k);
		break;

	case 1:
		w[0] = -1.0; /* Not correct at discontinuities */
		w[1] =  1.0; /* Not correct at discontinuities */
		break;
	case 0:
		w[0] = 0.0; /* Not correct at discontinuities */
		break;

	default:
		for(k=0; k<=7; k++) w[k] = dwt7(x-k);
	}
}


/***************************************************************************************
Work out what to do with positions outside the FOV
	i	- Co-ordinate (0<=i<m)
	m	- dimension
	returns reflected co-ordinate
*/
static int splins_mirror(int i, int m)
{
	int m2;
	i  = abs(i);
	if (i< m) return(i);
	if (m==1) return(0);
	m2 = (m-1)*2;
	i %= m2;
	return((i<m) ? i : m2-i);
}

/***************************************************************************************
Work out what to do with positions outside the FOV
        i       - Co-ordinate (0<=i<m)
        m       - dimension
        returns wrapped co-ordinate

        For MRI, it may be better to wrap the boundaries
        - especially in the read and phase encode directions.
*/
static int splins_wrap(int i, int m)
{
	if (i<0) return(m-1-((-i-1) % m));
	return(i % m);
}

/***************************************************************************************
Resample a point
	c	- Volume of B-spline coefficients
	m0,m1,m2	- dimensions of c
	x0,x1,x2	- co-ordinate to sample
	d	- degrees of splines used
	returns value of sampled point
*/
static double sample ( const double *coeffs, const mwSize *csize,
	double x0, double x1, double x2, int d[],
	int (*bnd[])())
{
	double w0[32], w1[32], w2[32]; /* B-spline weights */
	int    o0[32], o1[32], o2[32]; /* Offsets */
	int    i0,     i1,     i2;     /* Initial offsets */
	double d0,     d1,     d2;     /* Used by seperable convolution */
	int k;
	const double *cp;

	/* Generate seperable B-spline basis functions */
	weights(d[0], x0, &i0, w0);
	weights(d[1], x1, &i1, w1);
	weights(d[2], x2, &i2, w2);

	/* Create lookups of voxel locations - for coping with edges */
	for(k=0; k<=d[0]; k++) o0[k] = bnd[0](k+i0, csize [0]);
	for(k=0; k<=d[1]; k++) o1[k] = bnd[1](k+i1, csize [1])*csize [0];
	for(k=0; k<=d[2]; k++) o2[k] = bnd[2](k+i2, csize [2])*(csize [1]*csize [0]);

	/* Convolve coefficients with basis functions */
	d2 = 0.0;
	for(i2=0; i2<=d[2]; i2++)
	{
		d1 = 0.0;
		for(i1=0; i1<=d[1]; i1++)
		{
			cp = coeffs+o2[i2]+o1[i1];
			d0 = 0.0;
			for(i0=0; i0<=d[0]; i0++)
				d0 += cp[o0[i0]] * w0[i0];
			d1 += d0 * w1[i1];
		}
		d2 += d1 * w2[i2];
	}
	return(d2);
}


/***************************************************************************************
Resample a point and its gradients
	c	- Volume of B-spline coefficients
	m0,m1,m2	- dimensions of c
	x0,x1,x2	- co-ordinate to sample
	d	- degrees of splines used
	pg0,pg1,pg2	- gradients
	returns value of sampled point
*/
static double dsample ( const double *coeffs, const mwSize *csize,
	double x0, double x1, double x2,
	int d[], double *pg0, double *pg1, double *pg2,
	int (*bnd[])())
{
	double  w0[32],  w1[32],  w2[32]; /* B-spline weights */
	double dw0[32], dw1[32], dw2[32]; /* B-spline derivatives */
	int     o0[32],  o1[32],  o2[32]; /* Offsets */
	int     i0,      i1,      i2;     /* Initial offsets */
	double  d0,      d1,      d2;     /* Used by seperable convolution */
	double g00, g10,g11, g20,g21,g22; /* Used for generating gradients */
	int k;
	const double *cp;

	/* Generate seperable B-spline basis functions */
	weights(d[0], x0, &i0, w0);
	weights(d[1], x1, &i1, w1);
	weights(d[2], x2, &i2, w2);

	dweights(d[0], x0, &i0, dw0);
	dweights(d[1], x1, &i1, dw1);
	dweights(d[2], x2, &i2, dw2);

	/* Create lookups of voxel locations - for coping with edges */
	for(k=0; k<=d[0]; k++) o0[k] = bnd[0](k+i0, csize [0]);
	for(k=0; k<=d[1]; k++) o1[k] = bnd[1](k+i1, csize [1])*csize [0];
	for(k=0; k<=d[2]; k++) o2[k] = bnd[2](k+i2, csize [2])*(csize [0]*csize [1]);

	/* Convolve coefficients with basis functions */
	g20 = g21 = g22 = d2 = 0.0;
	for(i2=0; i2<=d[2]; i2++)
	{
		g10 = g11 = d1 = 0.0;
		for(i1=0; i1<=d[1]; i1++)
		{
			cp = coeffs+o2[i2]+o1[i1];
			g00 = d0  = 0.0;
			for(i0=0; i0<=d[0]; i0++)
			{
				d0  += cp[o0[i0]] *  w0[i0];
				g00 += cp[o0[i0]] * dw0[i0];
			}
			d1  += d0  *  w1[i1];
			g10 += g00 *  w1[i1];
			g11 += d0  * dw1[i1];
		}
		d2  += d1  *  w2[i2];
		g20 += g10 *  w2[i2];
		g21 += g11 *  w2[i2];
		g22 += d1  * dw2[i2];
	}
	*pg0 = g20;
	*pg1 = g21;
	*pg2 = g22;

	return(d2);
}


/***************************************************************************************
Loop through data and resample the points
	c	- Volume of B-spline coefficients
	m0,m1,m2	- dimensions of c
	n	- number of points to resample
	x0,x1,x2	- array of co-ordinate to sample
	d	- degree of spline used
	cond	- code determining boundaries to mask at
	bnd	- functions for dealing with edges
	f	- resampled data
*/
#define TINY 5e-2

static void fun ( const double *coeffs, const mwSize *csize,
	int n, double x0[], double x1[], double x2[], int d[],
	int cond, int (*bnd[])(), double f[])
{
	int j;
	double NaN = mxGetNaN();

	for(j=0; j<n; j++)
	{
		if (((cond&1) | (x0[j]>=1-TINY && x0[j]<=csize [0]+TINY)) &&
			((cond&2) | (x1[j]>=1-TINY && x1[j]<=csize [1]+TINY)) &&
			((cond&4) | (x2[j]>=1-TINY && x2[j]<=csize [2]+TINY)))
			f[j] = sample ( coeffs, csize, x0 [j] - 1, x1 [j] - 1, x2 [j] - 1, d, bnd );
		else
			f[j] = NaN;
	}
}


/***************************************************************************************
Loop through data and resample the points and their derivatives
	c	- Volume of B-spline coefficients
	m0,m1,m2	- dimensions of c
	n	- number of points to resample
	x0,x1,x2	- array of co-ordinate to sample
	d	- degrees of splines used
	cond	- code determining boundaries to mask at
	bnd	- functions for dealing with edges
	f	- resampled data
	df0, df1, df2	- gradients
*/
static void dfun ( const double *coeffs, const mwSize *csize,
	int n, double x0[], double x1[], double x2[],int d[],
	int cond, int (*bnd[])(),
	double f[], double df0[], double df1[], double df2[])
{
	int j;
	double NaN = mxGetNaN();

	for(j=0; j<n; j++)
	{
		if (((cond&1) | (x0[j]>=1-TINY && x0[j]<=csize [0]+TINY)) &&
			((cond&2) | (x1[j]>=1-TINY && x1[j]<=csize [1]+TINY)) &&
			((cond&4) | (x2[j]>=1-TINY && x2[j]<=csize [2]+TINY)))
			f[j] = dsample ( coeffs, csize, x0 [j] - 1,x1 [j] - 1,x2 [j] - 1, d,
				&df0 [j], &df1 [j], &df2 [j], bnd );
		else
			f[j] = NaN;
	}
}


/***************************************************************************************/
void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs [] ) {
    
	int k, dim, degree [3];
    double *coeffs, *vol;
	double *x0, *x1, *x2;
    double *df0, *df1, *df2;
    mwSize cdims, dims;
    const mwSize *csize, *tsize;
    mwSize size [3] = { 1, 1, 1 };
	int ( * splins [3] ) ();
	int cond = 0;

	/* Usage:
			f = function(c,x0,x1,x2,d)
				c - B-spline coefficients
				x0, x1, x2 - co-ordinates
				d	- B-spline degree
				f	- sampled function
	   or:
			[f,df0,df1,df2] = function(c,x0,x1,x2,d)
				c - B-spline coefficients
				x0, x1, x2 - co-ordinates
				d	- B-spline degree
				f	- sampled function
				df0, df1, df2	- sampled derivatives
	*/
	if (nrhs < 5 || nlhs>4)
		mexErrMsgTxt("Incorrect usage.");

	for(k=0; k<5; k++)
	{
		if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) ||
			mxIsSparse(prhs[k]) || !mxIsDouble(prhs[k]))
			mexErrMsgTxt("Input must be numeric, real, full and double precision.");
	}
    
    if ( mxGetNumberOfDimensions ( prhs [0] ) != 3 )
        mexErrMsgTxt ( "This function requires a 3-D volume as input." );
    
    if ( mxGetM ( prhs [4] ) * mxGetN ( prhs [4] ) != 1 && mxGetM ( prhs [4] ) * mxGetN ( prhs [4] ) != 3 )
		mexErrMsgTxt ( "The degree parameter must have 1 or 3 elements." );
    
    if ( nrhs > 5 )
        if ( mxGetM ( prhs [5] ) * mxGetN ( prhs [5] ) != 3 )
            mexErrMsgTxt ( "The warp parameter must have 3 elements." );
    
    
    /* Gets the degree in each dimmension. */
	for ( dim = 0; dim < 3; dim++ ) {
        if ( mxGetM ( prhs [4] ) * mxGetN ( prhs [4] ) == 3 )
            degree [dim] = (int) round ( mxGetPr ( prhs [4] ) [dim] );
        else
            degree [dim] = (int) round ( mxGetScalar ( prhs [4] ) );
        
		if ( degree [dim] < 0 || degree [dim] > 7 )
			mexErrMsgTxt ( "Invalid spline degree." );
	}
    
    /* Sets the spline function to use. */
	for ( dim = 0; dim < 3; dim++ )
        splins [dim] = splins_mirror;
    
    /* If warping in any dimension modifies the spline function. */
	if ( nrhs > 5)
		for ( dim = 0; dim < 3; dim++ )
			if ( mxGetPr ( prhs [5] ) [ dim + 3 ] ) {
				splins [dim] = splins_wrap;
				cond += 1<<k;
            }
    
    
    /* Gets the spline coefficent size. */
	cdims = mxGetNumberOfDimensions ( prhs [0] );
	csize = mxGetDimensions ( prhs [0] );
    
	/* Gets the output size. */
	dims = mxGetNumberOfDimensions ( prhs [1] );
	tsize = mxGetDimensions ( prhs [1] );
    
     /* Fix the output size for less than 3 dimensions. */
    for ( dim = 0; dim < dims; dim++ )
        size [dim] = tsize [dim];
    
    /* Makes sure that all the coordinate packs have the same dimensons. */
	if ( dims != mxGetNumberOfDimensions ( prhs [2] ) || dims != mxGetNumberOfDimensions ( prhs [3] ) )
		mexErrMsgTxt ( "Incompatible dimensions." );
	
	for ( dim = 0; dim < dims; dim++ )
		if ( size [dim] != mxGetDimensions ( prhs [2] ) [dim] || size [dim] != mxGetDimensions ( prhs [3] ) [dim] )
			mexErrMsgTxt("Incompatible dimensions.");
    
	/* Gets the input variables. */
	coeffs = mxGetPr ( prhs [0] );
	x0 = mxGetPr ( prhs [1] );
	x1 = mxGetPr ( prhs [2] );
	x2 = mxGetPr ( prhs [3] );

	/* Gets the output variable. */
	plhs [0] = mxCreateNumericArray ( dims, size, mxDOUBLE_CLASS, mxREAL );
	vol = mxGetPr ( plhs [0] );
    
	if ( nlhs <= 1 )
		fun ( coeffs, csize, size [0] * size [1] * size [2], x0, x1, x2, degree, cond, splins, vol );
	else
	{
		plhs [1] = mxCreateNumericArray ( dims, size, mxDOUBLE_CLASS, mxREAL );
		plhs [2] = mxCreateNumericArray ( dims, size, mxDOUBLE_CLASS, mxREAL );
		plhs [3] = mxCreateNumericArray ( dims, size, mxDOUBLE_CLASS, mxREAL );
		df0 = mxGetPr ( plhs [1] );
		df1 = mxGetPr ( plhs [2] );
		df2 = mxGetPr ( plhs [3] );
		dfun ( coeffs, csize, size [0] * size [1] * size [2], x0, x1, x2, degree, cond, splins, vol, df0, df1, df2 );
	}
}
