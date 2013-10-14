/*
	Eric Addison
	Utah State University
	Summer 2009

	This file contains my cosmology routines
*/


#include "cosmo.h"


//UNIT CONVERSION ROUTINES

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: m_to_ly													**/
/** Purpose: converts meters to light years								**/
/**	Inputs: double meters												**/
/**	Output: double light years											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double m_to_ly(double x)
{
	return x/METERPERLY;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: m_to_au													**/
/** Purpose: converts meters to astronomical units						**/
/**	Inputs: double meters												**/
/**	Output: double au													**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double m_to_au(double x)
{
	return x/METERPERAU;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: ly_to_m													**/
/** Purpose: converts light years to meters								**/
/**	Inputs: double light years											**/
/**	Output: double meters												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double ly_to_m(double x)
{
	return x*METERPERLY;
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: m_to_par													**/
/** Purpose: converts meters to parsecs									**/
/**	Inputs: double meters												**/
/**	Output: double parsecs												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double m_to_par(double x)
{
	return x/METERPERPARSEC;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: par_to_m													**/
/** Purpose: converts parsecs to meters									**/
/**	Inputs: double parsecs												**/
/**	Output: double meters												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double par_to_m(double x)
{
	return x*METERPERPARSEC;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: par_to_ly													**/
/** Purpose: converts parsecs to light years							**/
/**	Inputs: double parsecs												**/
/**	Output: double light years											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double par_to_ly(double x)
{
	return x*LYPERPARSEC;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: ly_to_par													**/
/** Purpose: converts light years to parsecs							**/
/**	Inputs: double light years											**/
/**	Output: double parsecs												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double ly_to_par(double x)
{
	return x/LYPERPARSEC;
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: sec_to_yr													**/
/** Purpose: converts seconds to years									**/
/**	Inputs: double seconds												**/
/**	Output: double years												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double sec_to_yr(double x)
{
	return x/SECPERYEAR;
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: yr_to_sec													**/
/** Purpose: converts years to seconds									**/
/**	Inputs: double years												**/
/**	Output: double seconds												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double yr_to_sec(double x)
{
	return x*SECPERYEAR;
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

//Distance conversions


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: Rl_to_z													**/
/** Purpose: converts luminosity distance to Redshift					**/
/**	Inputs: double luminosity distance									**/
/**	Output: double Redshift												**/
/**																		**/
/*************************************************************************/
/*************************************************************************/
double Rl_to_z(double Rl)
{

	return bisect_z(0,100,f1,Rl,1e-10);

}


//needed for Rl_to_z bisection
double f1(double z1, double Rl)
{
	//printf("\nnonRl part of f1(%.1f) = %.3e\n2*(1+z1-sqrt(1+z1)) = %.3e\t\tH0 = %.3e",z1,(2*(1+z1-sqrt(1+z1))/(H0)),2*(1+z1-sqrt(1+z1)),H0);
	return Rl - (2*(1+z1-sqrt(1+z1))/(H0));
}



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: z_to_Rl													**/
/** Purpose: converts Redshift to luminosity distance					**/
/**	Inputs: double Redshift												**/
/**	Output: double luminosity distance									**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double z_to_Rl(double z)
{
	return (2*(1+z-sqrt(1+z))/(H0));
}

/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: E_to_a														**/
/** Purpose: computes semi-major axis from energy						**/
/**	Inputs: double energy, double masses								**/
/**	Output: double semi-major axis										**/
/**																		**/
/*************************************************************************/
/*************************************************************************/
double E_to_a(double E, double m1, double m2)
{

	if(E >=0)
		return -1.0;
	else
		return -G*(m1*m2)/(2*E);
}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: L_E_to_e													**/
/** Purpose: computes eccentricity from energy and angular momentum		**/
/**	Inputs: double energy, double ang. mom, double masses				**/
/**	Output: double eccentricity											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/
double L_E_to_e(double E, double l, double m1, double m2)
{

	double k = G*(m1*m2), m = m1*m2/(m1+m2);

	return sqrt(fabs(1+(2*E*l*l)/(m*k*k)));

}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

//math routines needed for previous routines



/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: bisect_t_c													**/
/** Purpose: finds the root of a single valued function bracketed by	**/
/**				values a and b to specified tolerance					**/
/**				this is specifically for finding t_c from separation	**/
/**					or from wave frequency								**/
/**	Inputs: bracket values a and b, function pointer f, tolerance		**/
/**	Output: x value of the root											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double bisect_t_c(double a, double b, double f(double x, double a0, double m1, double m2), double tol, double a0, double m1, double m2)
{

	double c1;		//c will hold the midpoint value

	int n;

	//first a little error checking, f(a)*f(b) < 0 is required
	if(((*f)(b,a0,m1,m2))*((*f)(a,a0,m1,m2)) >= 0)
	{
		printf("\nERROR: Invalid interval endpoints: f(a)*f(b) > 0\n");
		return 0;
	}



	//find the number of iterations needed

	n = ceil( (log(b - a) - log(tol)) / log(2) );

	for(int j=0;j<n;j++)
	{
		c1 = a + 0.5*(b-a);		//compute midpoint

		if( ((*f)(c1,a0,m1,m2))*((*f)(a,a0,m1,m2)) < 0 )	//then root lies between a and c
			b = c1;					//set new upper bound

		else if( ((*f)(c1,a0,m1,m2))*((*f)(b,a0,m1,m2)) < 0 )	//then root lies between b and c
			a = c1;						//set new lower bound
	}

	return c1;

}


/*************************************************************************/
/*************************************************************************/
/**																		**/
/** Funtion: bisect_z													**/
/** Purpose: finds the root of a single valued function bracketed by	**/
/**				values a and b to specified tolerance					**/
/**				this is specifically for finding z from R_l				**/
/**	Inputs: bracket values a and b, function pointer f, tolerance		**/
/**	Output: x value of the root											**/
/**																		**/
/*************************************************************************/
/*************************************************************************/

double bisect_z(double a, double b, double f(double z1, double Rl), double Rl, double tol)
{

	double c1;		//c will hold the midpoint value

	int n;

	//first a little error c1hec1king, f(a)*f(b) < 0 is required
	if(((*f)(b,Rl))*((*f)(a,Rl)) >= 0)
	{
		printf("\nERROR: Invalid interval endpoints: f(a)*f(b) > 0\n");
		return 0;
	}



	//find the number of iterations needed

	n = ceil( (log(b - a) - log(tol)) / log(2) );

	for(int j=0;j<n;j++)
	{

		c1 = a + 0.5*(b-a);						//c1ompute midpoint

		if( ((*f)(c1,Rl))*((*f)(a,Rl)) < 0 )	//then root lies between a and c1
			b = c1;								//set new upper bound

		else if( ((*f)(c1,Rl))*((*f)(b,Rl)) < 0 )	//then root lies between b and c1
			a = c1;								//set new lower bound
	}

	return c1;

}



