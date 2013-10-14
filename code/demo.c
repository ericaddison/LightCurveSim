/*************************************
 File: demo.c
 Project: Light Curve Simulation (LC)
 Purpose: demo file to show how to call the simulation
 Author: Eric Addison
 Original date: 23Apr13
 *************************************/

/* This file demonstrates how to setup and call LCsim */

#include "LC.h"


int main()
{
		
	double **F;							// pointer for Flux array returned from LC function
	
	LCparams params;					// Light Curve parameters struct
	
	params.q = 0.125/0.68;				// mass ratio
	params.M = 0.125+0.68;				// total mass
	params.Porb = 1028.73;				// orbital period
	
	params.ia = 43;						// inclination angle
	params.lan = 0;						// longitude of ascending node / GW polarization angle
	params.dL = 1976;					// Luminosity distance (LY)
	
	params.Tpole = 13800;				// Pole temp of secondary
	params.Tdisk = 94763;				// Char. disc temp
	params.TWD = 20537;					// WD temp
	params.THS = 120000;				// HS max temp
	
	params.th0 = 0;						// initial binary phase
	params.hot_scale = 4.1737;			// ratio of HS radius to disc height
	params.hot_th = 8.0;				// HS phase on disc
	params.Rd = 0.48;					// Scaled disc radius
	
	params.dt = 5;						// time step (sampling period) of simulation, in seconds
	params.fileOut = 1;					// YES output to file
	sprintf(params.fileName,"FluxOut.dat");	// output file name
	params.verbose = 1;					// YES print screen output
	
	printf("\nCalling LCsim()..."); fflush(stdout);
	F = LCsim(&params);					// calling LCsim, assigns flux output to pointer F
										// F is an array that holds the flux output
	printf("\nLCsim complete\n\n");
	
	free(F[0]);
	free(F[1]);
	free(F);
	
  return 0;
}
