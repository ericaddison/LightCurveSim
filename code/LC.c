/*************************************
 Project: Light Curve Simulation (LC)
 Purpose: Main code file
 Author: Eric Addison
 Original date: 23Apr13
 *************************************/

/* This file holds the main program files for the light curve simulation */

#include "LC.h"


double **LCsim(LCparams *pars)
{
  /* Take input parameters:
	Tpole -- Temp at pole of Roche Lobe
	Tdisk -- Characteristic temp used for disc temp profile
	TWD -- constant temp of primary
	THS -- max temp of HS
	Rdisk -- radius of accretion disc
	expval -- cooling rate of HS
	hot_scale -- radius of HS relative to disk height
	th0 -- initial phase value
	ia -- inclination
	lan -- longitude of ascending node
	q -- mass ratio
	M -- total mass
	Porb -- orbital period
	d -- Luminosity distance
  */
		
/* Variable decls */
	double Tpole, Tdisk, TWD, THS;		// temperature variables
	double hot_scale, hot_th, expval;	// HS geometry variables
	double Rd, Hd, DiscCen[3];			// Disk radius, height, and center 
	double rWD, WDcen[3];				// WD radius and center coords
	double rHS, HScen[3];				// HS radius and center coords
	double th0,ia,lan,d;				// binary extrinsic variables
	double q,M,Porb;					// binary intrinsic variables
	double mu, omega, Aa, a, b;			// derived variables
	int num=8;							// number of objects in simulation
	double cm[num][3];					// array to hold center of mass of objects
	hullLC hulls[num];					// array to hold RL points	(replace with hull array)
	double *T[num];						// 2D temperature array
	double *t;							// time array
	double Mth[3][3], Mi[3][3], Mlan[3][3], M0[3][3], Mf[3][3];	// rotation matrices 
	double iMi[3][3], iMlan[3][3], Mtemp[3][3];
	unsigned *Nseen[num];
	coordT *pts2, *pts2x, *pts2y;
	int ind, ind2;
	double nz;							// for use in main loop step 1
	unsigned char *bytes;				// pointer for use as a byte array
	double *F, **Fout, Ffactor, Fdenom;			// flux output
	//double movAve[3];
	//int movAveN = 3;

	zero_mats(Mth,Mi,Mlan,M0,Mf);		// set all values to zero

/* timing variables */
	struct timeval t1,t2,t3;

/* sim variables */
	int norbits, npts, order[num];
	double dt, dth, tmax; 

/* unpack parameters */
	q = pars->q;
	M = pars->M;
	Porb = pars->Porb;
	ia = pars->ia*PI/180.0;
	lan = pars->lan*PI/180.0;
	d = pars->dL;
	Tpole = pars->Tpole;
	Tdisk = pars->Tdisk;
	TWD = pars->TWD;
	THS = pars->THS;
	th0 = pars->th0*PI/180.0;
	hot_scale = pars->hot_scale;
	hot_th = pars->hot_th*PI/180.0;
	dt = pars->dt;
	Rd = pars->Rd;
	
/* Variable calulation and scaling */

	mu = 1.0/(1.0+q);							// scaled mass 1, mu = m1/M
	omega = 2*PI/Porb;							// angular orbital frequency
	Aa = pow(G*M*MSUN/(omega*omega),1.0/3.0);   // binary semi-major axis
	a = 1-mu; b = mu;							// scaled mass positions
	d *= METERPERLY;
	d /= Aa;									// scaled lum. distance
	rWD = R_WD/Aa;								// scaled WD radius
	rWD = 6378e3/Aa;
	expval = 2.9985;

	// scale temperatures
	Tdisk /= Tpole;
	TWD /= Tpole;
	THS /= Tpole;
	
	// factors used in flux calculation
	Ffactor = SB*pow(Tpole,4.0);
	Fdenom = (4.0*PI*d*d);

if(pars->verbose)	
	gettimeofday(&t1, NULL); 	

	
/* Object Data */
	DiscCen[0] = a; DiscCen[1] = 0; DiscCen[2] = 0; 
	Hd = rWD/50.0;
	WDcen[0] = a; WDcen[1] = 0; WDcen[2] = 0; 
	
/* HS rotation angle */
	// possibly require solving for HS phase from disc radius here
	// i.e. code up maple simulation here
	//hot_th = 8.0*PI/180.0;			// HS phase on disc
	HScen[0] = Rd*cos(hot_th+PI) + DiscCen[0];
	HScen[1] = Rd*sin(hot_th+PI) + DiscCen[1];
	HScen[2] = 0;
	rHS = hot_scale*Hd;
	printf("rHS = %g\nhot_scale = %g\n",rHS,hot_scale);
	
/* Object creation -- calls to qhull */

	build_objects(hulls,mu,DiscCen,rWD,Rd,Hd,WDcen,HScen,rHS,pars->verbose);
	

/* Splitting HS and WD objects */
	// since no rotation of the system has happened yet, split HS and WD by looking for centers with z coords
	// above or below zero
	
	split_WD_and_HS(hulls, &num, pars->verbose);


/* Temperature Profiles */
	
	assign_temps(T,hulls,mu,rWD,Tdisk,DiscCen,HScen,THS,TWD,Rd,rHS,expval,num,pars->verbose);

	
/* CM arrays used for object ordering */

// cm[0] and cm[1] contains cm of roche lobe and disc system
	cm[0][0] = -b;
	cm[1][0] = a;
	cm[0][1] = cm[0][2] = cm[1][1] = cm[1][2] = 0;

	for(int i=0; i<3; i++)
		for(int j=2; j<num; j++)
			cm[j][i] = Cmean(hulls[j].C,hulls[j].nF,i);



/* Sim paramter setup */

	norbits = 1;			// number of orbits to complete
	dt = 5*omega;			// sampling period
	dth = dt;				// angular step
	tmax = norbits*2*PI;	// max time
	npts = (int) (tmax/dt);
	
	//malloc more arrays
	t = malloc(npts*sizeof(double));
	F = malloc((num+1)*sizeof(double));
	Fout = malloc(2*sizeof(double*));
	Fout[0] = malloc(npts*sizeof(double));
	Fout[1] = malloc(npts*sizeof(double));
	for(int i=0; i<num; i++)
		Nseen[i] = calloc(hulls[i].nF,sizeof(unsigned));
	

	// build rotation matrices
	build_rot_mat(M0,th0,'z');
	build_rot_mat(Mi,ia,'x');
	build_rot_mat(Mlan,lan,'z');
	build_rot_mat(Mth,dth,'z');
	build_rot_mat(iMi,-ia,'x');		// inverse of Mi
	build_rot_mat(iMlan,-lan,'z');	// inverse of Mlan

	// initial rotation matrix
	mat_mult(Mlan,Mi,Mtemp);
	mat_mult(Mtemp,M0,Mf);
	
	// perform initial rotations
	for(int i=0; i<num; i++)
	{
		mat_vect_array_mult(Mf,hulls[i].C,hulls[i].nF);
		mat_vect_array_mult(Mf,hulls[i].N,hulls[i].nF);
		mat_vect_array_mult(Mf,hulls[i].P,hulls[i].nV);
	}
	mat_vect_array_mult2(Mf,cm,num);


	// build overall matrix to rotate stuff at each step
	// Mf = Mlan*Mi*Mth*inv(Mi)*inv(Mlan)
	mat_mult(Mlan,Mi,Mtemp);
	mat_mult(Mtemp,Mth,Mf);
	mat_mult(Mf,iMi,Mtemp);
	mat_mult(Mtemp,iMlan,Mf);

	// array for 2D convex hulls
	pts2 = malloc(2*hulls[1].nP*sizeof(coordT));
	pts2x = malloc(hulls[1].nP*sizeof(coordT));
	pts2y = malloc(hulls[1].nP*sizeof(coordT));


num=6;


//	FILE *fp1 = fopen("Ftest.dat","w");	

	
/* main simulation loop */	
if(pars->verbose)	
	printf("\nEntering main loop...\n");
// includes both qhull and PNPOLY calls
	for(int i=0; i<npts; i++)
	{

		/* for each object, determine which normal vectors have some component pointing toward observer (-z direction)
	   that is, the patches of area that the observer can see */
		for(int j=0; j<num; j++)
		{
			for(int k=0; k<hulls[j].nF; k++)	// loop over all facets in this object		
			{
				nz = hulls[j].N[k][2];	
				bytes = (unsigned char *)&nz;
				// this checks the sign bit of the number to see if it's negative or not -- bit set to 1 => negative
				Nseen[j][k] = (unsigned) ((bytes[7] & SIGN_BIT_MASK)>>7);
			}
		}
		

	/* for the second object and onward (ordered by location of cm z coordinate) see which of these normal
		vectors are blocked by closer objects */
		
		// decide order of projection based on cm
		get_order(order,cm,ia);

	/* find which are blocked */
	for(int l=0;l<num;l++)	// loop over first num-1 objects
	{	
		// build point array for 2D qhull
		ind = order[l];
		for(int j=0; j<hulls[ind].nV;j++)
		{
			pts2[2*j] = hulls[ind].P[j][0];
			pts2[2*j+1] = hulls[ind].P[j][1];
		}

		// call 2D convhull
			qh_new_qhull(2, hulls[ind].nV, pts2, 0, "qhull ", NULL, stderr);
			extract_qhull2D_P2(qh facet_list, NULL, pts2x,pts2y);

		// now check the rest of the objects in the order and see which are blocked
		for(int j=l+1; j<num; j++) 				// loop over higher objects
		{
			ind2 = order[j];
			for(int k=0; k<hulls[ind2].nF; k++)		// loop over visible facets
			{
				// multiply current Nseen value by 1 if outside of object
				// shadows, and 0 if inside
				Nseen[ind2][k] *= npnpoly(qh num_vertices,pts2x,pts2y,hulls[ind2].C[k][0],hulls[ind2].C[k][1]); 	
			}
		}
	}


/*		
		int sum=0;
		for(int i=0; i<hulls[1].nF;i++)
			sum += Nseen[1][i];
		printf("\nsum = %d\n",sum);
*/		
		// now facets that can be seen are listed in Nseen with a one, otherwise it's a zero

	/* sum up the flux */
	F[num] = 0;
	for(int j=0; j<num; j++)
	{
		F[j] = 0;
		for(int k=0; k<hulls[j].nF; k++)
		{
			F[j] -= Nseen[j][k] * hulls[j].A[k] * hulls[j].N[k][2] * 
				pow(T[j][k],4.0) / Fdenom;
		}
			if(!isnan(F[j]) && !isinf(F[j]))
			   F[num] += F[j]*Ffactor;		// total flux
	}

	// keep a moving average going for Disc
/*		
	if(i==0)
		for(int j=0; j<movAveN; j++)
			movAve[j]=F[1];
	F[num] -= F[1]*Ffactor;

	for(int j=0; j<movAveN-1; j++)
	{	
		movAve[j] = movAve[j+1];
	}
	
	movAve[movAveN-1] = F[1];

	F[1] = 0;

	for(int j=0; j<movAveN-1; j++)
		F[1] += movAve[j]/(double)movAveN;
	F[num] += F[1]*Ffactor;
*/
	// output flux update
		Fout[0][i] = i*dt/omega;		// time stamp
		Fout[1][i] = F[num];		// current flux value
	
//		for(int j=0; j<num; j++)
//			fprintf(fp1,"%g, ",F[j]*Ffactor);
//		fprintf(fp1,"%g\n",F[num]);
		
	/* rotate system */
	for(int i=0; i<num; i++)
	{
		mat_vect_array_mult(Mf,hulls[i].C,hulls[i].nF);
		mat_vect_array_mult(Mf,hulls[i].N,hulls[i].nF);
		mat_vect_array_mult(Mf,hulls[i].P,hulls[i].nV);
	}
	mat_vect_array_mult2(Mf,cm,num);

	if(pars->verbose)
	{
		printf("\b\b\b   \rProgress: %.1f%%",(double)i/(double)npts*100.0);
		fflush(stdout);
	}


}

	if(pars->verbose)
	{
		printf("\b\b\b   \rProgress: %.1f%%",100.0);
		fflush(stdout);
		printf("\n");
	}	
		
//	fclose(fp1);
	
/* Preparing output */
	
	// file output if requested
	if(pars->fileOut)		
	{
		FILE *fp = fopen(pars->fileName,"w");
		for(int j=0; j<npts; j++)
			fprintf(fp,"%.2f, %.10e\n ",Fout[0][j],Fout[1][j]);
			fclose(fp);
	}
	

/* Writing to Files */
	// calculate wall time
if(pars->verbose)
{	
	gettimeofday(&t2, NULL);
	timeval_subtract(&t3,&t2,&t1); 	
	printf("\nElapsed time = %d.%d s\n",t3.tv_sec,t3.tv_usec);
}

	//write hulls to file
/*
	printf("\nWriting hulls to files:\n");
	printf("\tWriting hull to file \"Roche.xxx\"...");
	write_hull_to_file(&(hulls[0]),"Roche");
	printf("Success!\n");		   
	printf("\tWriting hull to file \"disc.xxx\"...");
	write_hull_to_file(&(hulls[1]),"disc");
	printf("Success!\n");
	printf("\tWriting hull to file \"WDtop.xxx\"...");
	write_hull_to_file(&(hulls[2]),"WDtop");
	printf("Success!\n");	
	printf("\tWriting hull to file \"WDbot.xxx\"...");
	write_hull_to_file(&(hulls[3]),"WDbot");
	printf("Success!\n");	
	printf("\tWriting hull to file \"HStop.xxx\"...");
	write_hull_to_file(&(hulls[4]),"HStop");
	printf("Success!\n");
	printf("\tWriting hull to file \"HSbot.xxx\"...");
	write_hull_to_file(&(hulls[5]),"HSbot");
	printf("Success!\n");
*/
	
	
	
/* Free allocated memory */
	qh_freeqhull(!qh_ALL);	// free qhull memory
	for(int i=0; i<num; i++)
	{
		free(T[i]);
		free(Nseen[i]);
	}
	FreeHulls(hulls,num);
	free(t);
	free(pts2);
	free(pts2x);
	free(pts2y);
	free(F);

  return Fout;
}
