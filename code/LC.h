/*************************************
 File: LC.h
 Project: Light Curve Simulation (LC)
 Purpose: Main program header file
 Author: Eric Addison
 Original date: 23Apr13
 *************************************/

/* This is the main header file for the LC simulatin program */

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>

/* GSL include files */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/* project include files */
#include "hullLC.h"
#include "cosmo.h"

/* Defined contants */
#define SIGN_BIT_MASK 0x80

/* data structures */


	// parameter struct to pass into the simulation
	typedef struct{
	
		/* Basic variables */
		double q;			// mass ratio
		double M;			// total mass (in units of solar masses)
		double Porb;		// orbital period (in seconds)
	
		/* location variables */
		double ia;			// inclination angle (degrees)
		double lan;			// longitude of the ascending node / GW polarization angle (degrees)
		double dL;			// Luminosity distance (LYs)
	
		/* temperature variables */	
		double Tpole;		// pole temperature of the secondary star, also used as temp scale factor
		double Tdisk;		// characteristic temp of the disc, based on BOB model
		double TWD;			// constant temp of the primary WD
		double THS;			// max temperature on the HS
	
		/* geometry variables */	
		double th0;			// initial binary phase (degrees)
		double hot_scale;	// size of HS radius relative to disc height -- 3 or 4 seems to work well for fitting AmCVN
		double hot_th;		// phase of the HS on the disc, measured from -x axis on the disc
		double Rd;			// scaled disc radius
	
		/* simulation parameters */
		double dt;			// time step (sampling period) of simulation, in seconds
		int fileOut:1;		// flag to signify whether to ouput to file or not
		char fileName[25];	// string to hold filename if outputting
		int verbose:1;		// flag to signify whether to output info to screen
	
		/* output variables */
		int Ns;				// number of samples in the output
	
	}LCparams;



/* Function Decls -- LC.c */
double **LCsim(LCparams *pars);


/* Function Decls -- LC_funcs.c */
void createSphere(double *center, double R, int n, hullLC *h); 
void sphere_splitter(hullLC *h0, hullLC *hT, hullLC *hB);
void createDisc(double *center, double Ri, double Ro, double H, int n, hullLC *h);
void disc_temp(double *T, hullLC *h, double R, double Tdisk, coordT *center);
void get_order(int *order, double cm[][3], double ia);
void build_objects(hullLC *hulls, double mu, double DiscCen[3], double rWD, double Rd, double Hd, double WDcen[3], double HScen[3], double rHS, int verb);
void assign_temps(double *T[], hullLC *hulls, double mu, double rWD, double Tdisk, double DiscCen[3], double HScen[3], double THS, double TWD, double Rd, double rHS, double expval, int num, int verb);
void split_WD_and_HS(hullLC *hulls, int *num, int verb);


/* Function Decls -- Roche.c */
double get_roche_cloud(double mu, int n, hullLC *h);
void roche_temp(double *T, hullLC *h, double mu);
double calc_L1(double mu, double tol);
double RL_VXp(double x, double mu);
double RL_VXpp(double x, double mu);
double RL_V(double x, double y, double mu, double r1, double r2);
void RL_common_r(double *r1, double *r2, double x, double y, double mu);
void RL_common_Vp(double *Vpx, double *Vpy, double x, double y, double mu,double r1, double r2);
double **RL_get_contour(double mu, int *n);


/* function Decls -- misc_math.c */
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
int npnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);
double Cmean(double **x, int n, int coord);
void mat_vect_array_mult(double M[][3], double **v, int n);
void mat_vect_array_mult2(double M[][3], double v[][3], int n);
void zero_mats(double M1[][3], double M2[][3], double M3[][3], double M4[][3], double M5[][3]);
void build_rot_mat(double M[][3], double th, char axis);
void mat_mult(double A[][3], double B[][3], double C[][3]);




/* Error Codes */
#define RL_CONT_ERROR 1		// error in function RL_get_contour
