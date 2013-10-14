/*************************************
 File: misc_math.c
 Project: Light Curve Simulation (LC)
 Purpose: miscellaneous math functions
 Author: Eric Addison
 Original date: 25Apr13
 *************************************/

/* file containing some random mathy functions */

#include "LC.h"


/*------------------------------------------------------------	
 function: pnpoly
 purpose: Point in Polygon algorithm
 inputs: number of vertices nvert, pointer to coords of vertices *vertx and 		*verty,	test point coords testx and testy
 notes: stolen directly from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
 ------------------------------------------------------------*/
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

/*------------------------------------------------------------	
 function: npnpoly
 purpose: Point in Polygon algorithm
 inputs: number of vertices nvert, pointer to coords of vertices *vertx and 		*verty,	test point coords testx and testy
 outputs: returns 1 if NOT in the polygon, 0 if it is
 notes: stolen directly from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
 ------------------------------------------------------------*/
int npnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return !c;
}

/*------------------------------------------------------------	
 function: mat_array_mult
 purpose: multiply a matrix by an array of vectors
 inputs: the matrix M, array of vectors v, and number of vectors in the array n
 notes: overwrites v with result
 ------------------------------------------------------------*/
void mat_vect_array_mult(double M[][3], double **v, int n)
{
	double temp[3];

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<3; j++)
		{
			temp[j] = 0;
			for(int k=0; k<3; k++)
				temp[j] += M[j][k]*v[i][k];
		}
		v[i][0] = temp[0];
		v[i][1] = temp[1];
		v[i][2] = temp[2];
	}
}

/*------------------------------------------------------------	
 function: mat_array_mult2
 purpose: multiply a matrix by an array of vectors
 inputs: the matrix M, array of vectors v, and number of vectors in the array n
 notes: overwrites v with result, takes v[][3] instead of **v
 ------------------------------------------------------------*/
void mat_vect_array_mult2(double M[][3], double v[][3], int n)
{
	double temp[3];

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<3; j++)
		{
			temp[j] = 0;
			for(int k=0; k<3; k++)
				temp[j] += M[j][k]*v[i][k];
		}
		v[i][0] = temp[0];
		v[i][1] = temp[1];
		v[i][2] = temp[2];
	}
}


// just zero out some matrices
void zero_mats(double M1[][3], double M2[][3], double M3[][3], double M4[][3], double M5[][3])
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			M1[i][j] = M2[i][j] = M3[i][j] = M4[i][j] = M5[i][j] = 0;
		}
	}
}

// set up a rotation matrix
void build_rot_mat(double M[][3], double th, char axis)
{

	switch(axis)
	{

		case 'x':
			M[0][0] = 1;
			M[1][1] = cos(th);
			M[1][2] = -sin(th);
			M[2][1] = sin(th);
			M[2][2] = cos(th);
		break;

		case 'y':
			M[1][1] = 1;
			M[0][0] = cos(th);
			M[0][2] = sin(th);
			M[2][0] = -sin(th);
			M[2][2] = cos(th);
		break;

		case 'z':
			M[2][2] = 1;
			M[0][0] = cos(th);
			M[0][1] = -sin(th);
			M[1][0] = sin(th);
			M[1][1] = cos(th);
		break;

		default:
			fprintf(stderr,"Error: build_rot_mat(): invalid axis choice -- output set to identity.\n");
			M[0][0] = M[1][1] = M[2][2] = 1;
			M[0][1] = M[0][2] = M[1][0] = M[1][2] = M[2][0] = M[2][1] = 0;
	}
}


void mat_mult(double A[][3], double B[][3], double C[][3])
{
	
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			C[i][j] = 0;
			for(int k=0; k<3; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}			
}


int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
	/* Perform the carry for the later subtraction by updating y. */
	if (x->tv_usec < y->tv_usec) {
		int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000) {
		int nsec = (x->tv_usec - y->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}
	
	/* Compute the time remaining to wait.
	 tv_usec is certainly positive. */
	result->tv_sec = x->tv_sec - y->tv_sec;
	result->tv_usec = x->tv_usec - y->tv_usec;
	
	/* Return 1 if result is negative. */
	return x->tv_sec < y->tv_sec;
}

// coord = 0,1,2 for x,y,z
double Cmean(double **x, int n, int coord)
{
	double ave=0;
	for(int i=0; i<n; i++)
		ave +=	x[i][coord];
	return ave/n;
}
