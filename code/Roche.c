/*************************************
 File: Roche.c
 Project: Light Curve Simulation (LC)
 Purpose: LC Roche Lobe functions file
 Author: Eric Addison
 Original date: 23Apr23
 *************************************/

/* This file defines functions asociated with finding and building the Roche Lobe
 */

#include "LC.h"
int N;

/*------------------------------------------------------------	
	function: get_roche_cloud
	purpose: generate point cloud for Roche Lobe
	inputs: scaled mass mu, empty pointer *P, number of points n
------------------------------------------------------------*/
double get_roche_cloud(double mu, int n, hullLC *h)
{
	coordT **RLxy;		// array to hold RL contour points
	coordT *RLinterp;	// interpolated points along RL
	double th, dth, dx;	// theta variabes for making circles
	double lRL, L1;		// length of RL, L1 value
	int nc,np,n_cont;	// number of circles to make, number of points per circle, and number of points in the contour
	
	L1 = calc_L1(mu,1e-6);
	//printf("\nL1 = %.5f\n",L1);
	
	/* calculate top half of Roche Lobe */
	RLxy = RL_get_contour(mu,&n_cont);
	lRL = fabs(RLxy[0][0] - RLxy[0][n_cont-1]);		// length of RL
	
	/* interpolate to get the right number of points */

	// spline initialization
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n_cont);
	gsl_spline_init (spline, RLxy[0], RLxy[1], n_cont);
	
	// memory for RLinterp
	nc = np = floor(sqrt(n));
	RLinterp = malloc(nc*sizeof(coordT));
	// call spline
	dx = (lRL/(double)(nc));
	for(int i=0; i<nc; i++)
		RLinterp[i] = gsl_spline_eval (spline, L1-(i+1)*dx, acc);
	
	
	
	/* make circles traversing the RL surface */
	
	h->Pq = malloc((3*np*nc+6)*sizeof(coordT));
	
	dth = 2*PI/(double)np;
	// assign first point (no circle)
	h->Pq[0] = L1;
	h->Pq[1] = h->Pq[2] = 0;

	// assign circles
	for(int i=0; i<nc; i++)
	{
		th=0;
		for(int j=0; j<np; j++)
		{
			h->Pq[3+i*np*3 + 3*j + 0] = L1-(i+1)*dx;
			h->Pq[3+i*np*3 + 3*j + 1] = RLinterp[i]*cos(th);
			h->Pq[3+i*np*3 + 3*j + 2] = RLinterp[i]*sin(th);
			th += dth;
		}
		
	}

	// assign last point (no circle)
	h->Pq[6+nc*np*3-1] = 0;
	h->Pq[6+nc*np*3-2] = 0;
	h->Pq[6+nc*np*3-3] = RLxy[0][0]-0.01;
	
	h->nP = (unsigned)(np*nc+2);	
	
	/* free allocated memory */
	for(int i=0;i<2;i++)
		free(RLxy[i]);
	free(RLxy);
	free(RLinterp);
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	return L1;
}



/*------------------------------------------------------------	
 function: roche_temp
 purpose: calculates and assigns temperature for roche lobe points 
 inputs: result array *T, hullLC pointer of hull to assign temps, scaled m1 mass mu
 notes:	this assumes tpole set to Tpole=1, assuming all temps in the simulation
		are scaled by Tpole
 ------------------------------------------------------------*/
void roche_temp(double *T, hullLC *h, double mu) 
{
	
	int zi=0;
	double a = 1.0-mu, b = mu, q = (1.0-mu)/mu;
	double x,y,z;
	double c1, c2, dVdx, dVdy, dVdz;
	double g[h->nF], gg;
	
	
	/* calculate required partial dervatives for g */
	for(int i=0; i<h->nF;i++)
	{
		x = h->C[i][0];
		y = h->C[i][1];
		z = h->C[i][2];
		
		c1 = 1.0/pow( x*x - 2.0*x*a + a*a + y*y + z*z , 1.5 );
		c2 = q/pow(x*x + 2.0*x*b + b*b + y*y + z*z, 1.5 );
		
		dVdx = mu*(c1*(x-a)+c2*(x+b)) - x;
		dVdy = mu*y*(c1+c2) - y;
		dVdz = mu*z*(c1+c2);
		
		// now local g is
		g[i] = sqrt(dVdx*dVdx + dVdy*dVdy + dVdz*dVdz);
		
		T[i] = pow(g[i],0.25);
		
	}
	

	for(int i=0; i<h->nF;i++)
	{
		if(fabs(h->C[i][2]) > fabs(h->C[zi][2]))
			zi = i;
	}
	
	/* need to normalize by gpole^0.25 (g at max value of z) */
	gg = 1.0/pow(g[zi],0.25);
	for(int i=0; i<h->nF;i++)
		T[i] *= gg;
	
}



/*------------------------------------------------------------	
 function: calc_L1
 purpose: calculates the first Lagrange point L1 
 inputs: scaled mass mu, error tolerance tol
 ------------------------------------------------------------*/
double calc_L1(double mu, double tol)
{
    double temp = 999;
	double L1 = 0;

	/* Solve with Newton's method */
    while(fabs(L1 - temp) > tol)
    {
        temp = L1;
        L1 -= RL_VXp(L1,mu)/RL_VXpp(L1,mu);
    }
	
	return L1;
}


/*------------------------------------------------------------	
 function: RL_VXp
 purpose: calculates first x-derivative of two-body potential along the x-axis
			used in calculating L1 point
 inputs: position x, scaled mass mu, 
 ------------------------------------------------------------*/
double RL_VXp(double x, double mu)
{
    return -(1-mu)/(x+mu)/(x+mu) + mu/(x-1+mu)/(x-1+mu) + x;
}

/*------------------------------------------------------------	
 function: RL_VXpp
 purpose: calculates second x-derivative of two-body potential along the x-axis
 			used in calculating L1 point
 inputs: position x, scaled mass mu 
 ------------------------------------------------------------*/
double RL_VXpp(double x, double mu)
{
    return 2*(1-mu)/(x+mu)/(x+mu)/(x+mu) - 2*mu/(x-1+mu)/(x-1+mu)/(x-1+mu) + 1;
}


/*------------------------------------------------------------	
 function: RL_get_contour
 purpose: calculates Roche Lobe contour with 2nd order Runge-Kutta
 inputs: scaled mass mu 
 source: Astophysics with a PC by Paul Hellings
 ------------------------------------------------------------*/
double **RL_get_contour(double mu, int *n)
{
    double x,y,r2,r1,Vv,K,Vpx,Vpy,dx,dy,xh,yh,r2h,r1h,Vpxh,Vpyh,step;
    double L1 = calc_L1(mu,1e-6);
    double xi = L1;
    double yi = 0;
	double **pts;	// output array
	int i=0;
	
	if(mu >= 0.05)
		N = 200;
	else
	{
		N=0;
		fprintf(stderr,"Error: must use mu>0.05 in RL_get_contour (or adjust array size N)\n");
		exit(RL_CONT_ERROR);
	}

	/* malloc space for pts */
	pts = malloc(2*sizeof(double*));
	for(i=0;i<2;i++)
		pts[i] = calloc(N,sizeof(double));
	
		dx = 0.0;
		dy = 0.01;
		step = dy;
		x=xi-0.000001;
		y=yi+0.00001;
		//compute common values
		RL_common_r(&r1,&r2,x,y,mu);
		RL_common_Vp(&Vpx,&Vpy,x,y,mu,r1,r2);
		
		//compute potential V(x,y)
		K = RL_V(x,y,mu,r1,r2);

	i=0;
	//for(int i=1;i<200;i++)
	while(x<L1)
	{
            RL_common_r(&r1,&r2,x,y,mu);
            RL_common_Vp(&Vpx,&Vpy,x,y,mu,r1,r2);
            // if |dx| > |dy|, use equation (28), else use (29)
			if( fabs(dx) >= fabs(dy) )      //then use (28) to find dy
			{
				if(dx > 0)
					dx = step;
				else
					dx = -step;
				
				
				//half steps
				xh = x + 0.5*dx;
				yh = y - 0.5*dx*Vpx/Vpy;
				
				//full steps
				x += dx;
				RL_common_r(&r1h,&r2h,xh,yh,mu);
				RL_common_Vp(&Vpxh,&Vpyh,xh,yh,mu,r1h,r2h);
				dy = -dx*Vpxh/Vpyh;
				y += dy;
				
				//correction
                RL_common_r(&r1,&r2,x,y,mu);
                RL_common_Vp(&Vpx,&Vpy,x,y,mu,r1,r2);
                Vv = RL_V(x,y,mu,r1,r2);
                y -= (Vv - K)/Vpy;
			}
			else                            //so use (29) to find dx
			{
				if(dy > 0)
					dy = step;
				else
					dy = -step;
				
				//half steps
				yh = y + 0.5*dy;
				xh = x - 0.5*dy*Vpy/Vpx;
				
				//full steps
				y += dy;
				RL_common_r(&r1h,&r2h,xh,yh,mu);
				RL_common_Vp(&Vpxh,&Vpyh,xh,yh,mu,r1h,r2h);
				dx = -dy*Vpyh/Vpxh;
				x += dx;
								
				//correction
                RL_common_r(&r1,&r2,x,y,mu);
                RL_common_Vp(&Vpx,&Vpy,x,y,mu,r1,r2);
                Vv = RL_V(x,y,mu,r1,r2);
                x -= (Vv - K)/Vpx;
			}
		if(y<0)
		{
			pts[0][i] = x;
			pts[1][i] = -y;
			i++;
		}
    }
	
	pts[0][i] = L1;
	pts[1][i] = 0;	
	*n = i;
	
    return pts;
}

/*------------------------------------------------------------	
 function: RL_V
 purpose: calculates two-body potential V(x,y)
 used in computing Roche Lobe
 inputs: position x,y, scaled mass mu, distance from binary masses r1,r2 
 ------------------------------------------------------------*/
double RL_V(double x, double y, double mu, double r1, double r2)
{
    return (mu)/r1 + (1-mu)/r2 + (x*x + y*y)/2.0;
}


/*------------------------------------------------------------	
 function: RL_common_r
 purpose: calculates distance from binary masses r1,r2
 used in computing Roche Lobe
 inputs: pointers r1,r2, position x,y, scaled mass mu
 ------------------------------------------------------------*/
void RL_common_r(double *r1, double *r2, double x, double y, double mu)
{
    *r1 = sqrt( (x-1+mu)*(x-1+mu) + y*y );
    *r2 = sqrt( (x+mu)*(x+mu) + y*y );
}


/*------------------------------------------------------------	
 function: RL_common_VP
 purpose: calculates common potential terms
 used in computing Roche Lobe
 inputs: pointers Vpx, Vpy, position x,y, scaled mass mu, distances r1,r2
 ------------------------------------------------------------*/
void RL_common_Vp(double *Vpx, double *Vpy, double x, double y, double mu,double r1, double r2)
{
	*Vpx = -(x+mu)*(1-mu)/pow(r2,3.0) - mu*(x-1+mu)/pow(r1,3.0) + x;
	*Vpy = -y*(1-mu)/pow(r2,3.0) - mu*y/pow(r1,3.0) + y;
	return;
}
