/*************************************
 File: LC_funcs.c
 Project: Light Curve Simulation functions (LC)
 Purpose: Generic Function File
 Author: Eric Addison
 Original date: 24Apr13
 *************************************/

/* This file holds functions crucial to the LC simulation */

#include "LC.h"



/*------------------------------------------------------------	
 function: build_objects
 purpose: call all necessary funcs to completely build all objects
 inputs: hulls struct array and all necessary parameters
 ------------------------------------------------------------*/
void build_objects(hullLC *hulls, double mu, double DiscCen[3], double rWD, double Rd, double Hd, double WDcen[3], double HScen[3], double rHS, int verb)
{
if(verb)
{
	printf("\nCreating Objects...\n\n");
	printf("\tCreating Roche Lobe...\n");
	printf("\t  calling get_roche_cloud()...");
}
	get_roche_cloud(mu,500,&(hulls[0]));
if(verb)
	printf("Success!\n");
	
if(verb)
	printf("\t  Calling Qhull for Roche Lobe points...");
	hulls[0].dim=3;
	call_qhull(&(hulls[0]));
if(verb)
	printf("Success!\n");
if(verb)
{
	printf("\n\tCreating Disc...\n");
	printf("\t  Calling createDisc2()...");
}
	createDisc(DiscCen, rWD, Rd, Hd, 5000, &(hulls[1]));
if(verb)
	printf("Success!\n");	
	
if(verb)
{
	printf("\n\tCreating WD sphere...\n");
	printf("\t  Calling createSphere()...");
}
	createSphere(WDcen, rWD, 500, &(hulls[2]));
if(verb)
	printf("Success!\n");
	
if(verb)	
	printf("\t  Calling Qhull for WD points...");
	hulls[2].dim=3;
	call_qhull(&(hulls[2]));
if(verb)
	printf("Success!\n");
	
if(verb)
{
	printf("\n\tCreating HS sphere...\n");
	printf("\t  Calling createSphere()...");
}
	if(rHS!=0.0)
		createSphere(HScen, rHS, 500, &(hulls[3]));
	else
	{
		printf("yo\n");
		createSphere(HScen, 1e-10, 10, &(hulls[3]));
	}
if(verb)
	printf("Success!\n");
	
if(verb)
	printf("\t  Calling Qhull for HS points...");
	hulls[3].dim=3;
	call_qhull(&(hulls[3]));
if(verb)
	printf("Success!\n");	
	
if(verb)	
	printf("\nFacet totals:\n\tRL -- %u\n\tDisc -- %u\n\tWD -- %u\n\tHS -- %u\n",hulls[0].nF,hulls[1].nF,hulls[2].nF,hulls[3].nF);
	
	
	
}


/*------------------------------------------------------------	
 function: assign_temps
 purpose: call all necessary funcs to assign temperature profiles
 inputs: hulls struct array and all necessary parameters
 ------------------------------------------------------------*/
void assign_temps(double *T[], hullLC *hulls, double mu, double rWD, double Tdisk, double DiscCen[3], double HScen[3], double THS, double TWD, double Rd, double rHS, double expval, int num, int verb)
{
	
	double Td, dc;
if(verb)	
	printf("\nAssigning Temperature Profiles\n");
	// malloc space for T arrays
	for(int i=0; i<num; i++)
		T[i] = malloc(hulls[i].nF*sizeof(double));
	
	// get Roche Lobe temp profile
if(verb)
	printf("\n\tAssigning Roche Lobe Temps...");	
	roche_temp(T[0],&(hulls[0]),mu);
if(verb)
	printf("Success!");
	
	// set Disc profile
if(verb)
	printf("\n\tAssigning Disc Temps...");	
	disc_temp(T[1],&(hulls[1]),rWD,Tdisk, DiscCen);
if(verb)
	printf("Success!");
	
	// set constant WD temp
if(verb)
	printf("\n\tAssigning WD Temps...");	
	for(int i=0; i<hulls[2].nF; i++)
	{
		T[2][i] = TWD;
		T[3][i] = TWD;
	}
if(verb)
	printf("Success!");
	
	// set exponential HS temp
if(verb)
	printf("\n\tAssigning HS Temps...");	
	Td = Tdisk * pow(rWD/(Rd-rHS),0.75) * pow(1.0 - sqrt(rWD/(Rd-rHS)), 0.25);	// matching temp on disc
	for(int i=0; i<hulls[4].nF; i++)
	{
		dc = -sqrt( pow(hulls[4].C[i][0],2.0) + pow(hulls[4].C[i][1],2.0) ) + sqrt( pow(HScen[0],2.0) + pow(HScen[1],2.0));	// distance from center of HS to point
		T[4][i] = (THS-Td)*exp(-expval*(dc/rHS + 1)) + Td;
	}
	for(int i=0; i<hulls[5].nF; i++)
	{
		dc = -sqrt( pow(hulls[5].C[i][0],2.0) + pow(hulls[5].C[i][1],2.0) ) + sqrt( pow(HScen[0],2.0) + pow(HScen[1],2.0));	// distance from center of HS to point
		T[5][i] = (THS-Td)*exp(-expval*(dc/rHS + 1)) + Td;
	}
if(verb)
	printf("Success!\n");
	/*
	 // file test
	 FILE *fp = fopen("temptest.dat","w");
	 for(int i=0;i<hulls[3].nF;i++)
	 fprintf(fp,"%.5f\n",T[3][i]);
	 fclose(fp);
	 */
	
}

void split_WD_and_HS(hullLC *hulls, int *num, int verb)
{

	if(verb)
		printf("\nSplitting sphere objects in half...");
	sphere_splitter(&(hulls[2]),&(hulls[4]),&(hulls[5]));
	sphere_splitter(&(hulls[3]),&(hulls[6]),&(hulls[7]));
	
	// free whole sphere hulls and switch hulls pointers around
	FreeHull(&(hulls[2]));
	FreeHull(&(hulls[3]));	
	hulls[2] = hulls[4];
	hulls[3] = hulls[5];
	hulls[4] = hulls[6];
	hulls[5] = hulls[7];
	*num=6;
	if(verb)
		printf("Success!\n");	
	
}


/*------------------------------------------------------------	
 function: createSphere
 purpose: generate point cloud for a sphere with ~n points
 inputs: center position *center, radius R, number of points n, hullLC struct result
 ------------------------------------------------------------*/
void createSphere(double *center, double R, int n, hullLC *h)
{

	double th, dth, dphi, z, s;
	int N = floor(sqrt(n)); // put the same number of latitudinal levels as points in the circles
	
	/* Allocate memory in h->Pq */
	h->Pq = malloc((3*N*N+6)*sizeof(coordT));
	
	/* make some circles */
	dth = PI/(double)(N+1);
	dphi = 2*PI/(double)N;


	for(int i=0; i<N; i++)
	{
		th = (i+1)*dth;
		z = R*cos(th)+center[2];
		s = sin(th);
		for(int j=0; j<N; j++)
		{
			h->Pq[3+i*N*3+j*3+0] = R*cos(j*dphi)*s+center[0]; 
			h->Pq[3+i*N*3+j*3+1] = R*sin(j*dphi)*s+center[1];
			h->Pq[3+i*N*3+j*3+2] = z;
		}
	}
	
	/* put in poles manually */
	
	h->Pq[0] = center[0]; h->Pq[1] = center[1]; h->Pq[2] = R + center[2];
	h->Pq[3+N*N*3] = center[0]; h->Pq[3+N*N*3+1] = center[1]; h->Pq[3+N*N*3+2] = -R + center[2];
	
	h->nP = (unsigned)(N*N+2);


	
}

/*------------------------------------------------------------	
 function: sphere_splitter
 purpose: split the hullLC of a sphere into two separate hemispheres
 inputs: three hullLC pointers: the original sphere h0, and the desinations
		  for the top and bottom hemispheres hT and hB
 ------------------------------------------------------------*/
void sphere_splitter(hullLC *h0, hullLC *hT, hullLC *hB)
{

	int iT=-1, iB=-1;
	/* malloc space for hT and hB */
	int Np = h0->nP;
	int Nf = h0->nF;
	
	double eps = 1e-8;
	
	hT->dim = hB->dim = 3;
	hT->nF = hB->nF = Nf;
	hT->nP = hB->nP = Np;
	hT->nV = hB->nV = Np;
	
	MallocHull(hT);
	MallocHull(hB);
	
	/* split up the points */
	for(int i=0; i<h0->nP; i++)
	{

		if(h0->P[i][2] > eps)
		{
			iT++;
			hT->P[iT][0] = h0->P[i][0];
			hT->P[iT][1] = h0->P[i][1];
			hT->P[iT][2] = h0->P[i][2];
		}
		else if (h0->P[i][2] < -eps)
		{
			iB++;
			hB->P[iB][0] = h0->P[i][0];
			hB->P[iB][1] = h0->P[i][1];
			hB->P[iB][2] = h0->P[i][2];
		}
		else
		{
			iT++;
			hT->P[iT][0] = h0->P[i][0];
			hT->P[iT][1] = h0->P[i][1];
			hT->P[iT][2] = h0->P[i][2];			
			iB++;
			hB->P[iB][0] = h0->P[i][0];
			hB->P[iB][1] = h0->P[i][1];
			hB->P[iB][2] = h0->P[i][2];
		}			
	}

	// reset point numbers
	hT->nP = hT->nV= iT;
	hB->nP = hB->nV = iB;	
	
	
	/* split up facet info */
	// not copying vertex V info -- that's really only needed for plotting pretty pictures
	// and we can really just plot spheres for pretty pictures...
	iT = iB = -1;
	for(int i=0; i<h0->nF; i++)
	{
		if(h0->C[i][2] > 0)
		{
			iT++;
			hT->C[iT][0] = h0->C[i][0];
			hT->C[iT][1] = h0->C[i][1];
			hT->C[iT][2] = h0->C[i][2];
			hT->N[iT][0] = h0->N[i][0];
			hT->N[iT][1] = h0->N[i][1];
			hT->N[iT][2] = h0->N[i][2];
			hT->A[iT] = h0->A[i];
		}
		else
		{
			iB++;
			hB->C[iB][0] = h0->C[i][0];
			hB->C[iB][1] = h0->C[i][1];
			hB->C[iB][2] = h0->C[i][2];
			hB->N[iB][0] = h0->N[i][0];
			hB->N[iB][1] = h0->N[i][1];
			hB->N[iB][2] = h0->N[i][2];
			hB->A[iB] = h0->A[i];
		}
	}	
	
	// reset facet numbers
	hT->nF = iT;
	hB->nF = iB;	
	
	
	
	
}



/*------------------------------------------------------------	
 function: disc_temp
 purpose: compute temp profile for acc disc (BOB ch 18, pg 663)
 inputs: temp array *T, hullLC *h, primary radius R, char. temp Tdisk, center point of disc
 ------------------------------------------------------------*/
void disc_temp(double *T, hullLC *h, double R, double Tdisk, coordT *center)
{

	double r;	// distance from center of disc to point
	double dx, dy, dz;
	
	for(int i=0; i<h->nF; i++)
	{
		dx = h->C[i][0] - center[0];
		dy = h->C[i][1] - center[1];
		dz = h->C[i][2] - center[2];
		
		r = sqrt( dx*dx + dy*dy + dz*dz );
		T[i] = Tdisk * pow(R/r,0.75) * pow(1.0 - sqrt(R/r), 0.25);
	}
	
	
	
}



/*------------------------------------------------------------	
 function: createDisc
 purpose: generate point cloud and complete hull for a Disc with ~n points
 inputs: center position *center, inner radius Ri, outer radius Ro,
			height h, number of points n, hullLC struct result
 ------------------------------------------------------------*/
void createDisc(double *center, double Ri, double Ro, double H, int n, hullLC *h)
{
	double cx = center[0],cy = center[1],cz = center[2];
	int i,j;							// loop iterators
	int Np = floor(sqrt(n));			// number of pointer per circle
	int N = (Np-1)/4;					// related to number of circles
	int Nc = 2*N-1;						// number of circles on top side
	int Ntop = Np*Nc;					// total number of points in the top side of hull
	int Nside = Np;						// total number of points along the side
	int Ftop = 4*Np*(N-1);				// number of facets on top side
	int Fside = 4*Np;					// number of facets along the side
	int Ntot = 2*Ntop+Nside;			// total number of points in the hull
	int Ftot = 2*Ftop+Fside;
	double dth = 2.0*PI/(double)Np;		// angular space between points on a circle
	double dr = (Ro-Ri)/(double)(N-1);	// radial space between circles
	double r, th, z, hT, l1, l2, dl, r1, c, s;	// various vars
	int fid, pid;

//printf("\nNp = %d\nN = %d\nNc = %d\nNtot = %d\nFtot = %d\n",Np,N,Nc,Ntot,Ftot);

	/* Allocate Memory */
	h->P = malloc(Ntot*sizeof(coordT*));
	for(i=0; i<Ntot; i++)
		h->P[i] = malloc(3*sizeof(coordT));

	h->V = malloc(Ftot*sizeof(unsigned*));
	h->N = malloc(Ftot*sizeof(coordT*));
	h->C = malloc(Ftot*sizeof(coordT*));
	h->A = malloc(Ftot*sizeof(double));
	for(i=0; i<Ftot; i++)
	{
		h->V[i] = malloc(3*sizeof(unsigned));
		h->N[i] = malloc(3*sizeof(coordT));
		h->C[i] = malloc(3*sizeof(coordT));
	}

	/* Create vertices for top and bottom */
	z = H/2.0;						// z value for top side
	for(i=0; i<N; i++)				// main circles
	{
		r = Ri+i*dr;				// current radius
		for(j=0; j<Np; j++)			// loop over angles in this circle
		{
			pid = 2*i*Np + j;
			th = j*dth;
			h->P[pid][0] = h->P[pid + Ntop][0] = r*cos(th)+cx;		// assign x value
			h->P[pid][1] = h->P[pid + Ntop][1] = r*sin(th)+cy;		// assign y value
			h->P[pid][2] = z+cz;						// assign z value
			h->P[pid + Ntop][2] = -z+cz;				// assign z value
		}
	}

	for(i=0; i<N-1; i++)			// in between circles for trap centers
	{	
		r1 = (Ri + i*dr);			// current inner radius
		l1 = 2.0*r1*sin(dth/2.0);		// length of short end of trap
		l2 = 2.0*(r1+dr)*sin(dth/2.0);	// lenght of long end of trap
		dl = l2-l1;
		hT = sqrt(dr*dr - dl*dl/4.0);
		r = r1+hT/2.0;//sqrt(r1*r1 - l1*l1/4.0) + l1/(l1+l2)*sqrt(dr*dr-dl*dl/4.0);
		for(j=0; j<Np; j++)						// loop over angles in this circle
		{
			pid = Np+2*i*Np + j; 
			th = (j+0.5)*dth;
			h->P[pid][0] = h->P[pid + Ntop][0] = r*cos(th)+cx;		// assign x value
			h->P[pid][1] = h->P[pid + Ntop][1] = r*sin(th)+cy;		// assign y value
			h->P[pid][2] = z+cz;						// assign z value
			h->P[pid + Ntop][2] = -z+cz;				// assign z value
		}
	}

/* create facets for top and bottom */
	for(i=0;i<Nc-1;i+=2)
	{
		r1 = (Ri + (i/2)*dr);				// current inner radius
		l1 = 2.0*r1*sin(dth/2.0);		// length of short end of trap
		l2 = 2.0*(r1+dr)*sin(dth/2.0);	// lenght of long end of trap
		dl = l2-l1;
		hT = sqrt(dr*dr - dl*dl/4.0);		
		for(j=0; j<Np; j++)
		{
			th = (j+0.5)*dth;
			fid = 4*Np*(i/2) + 4*j;		// current facet id base value
			c = cos(th); s = sin(th);
			// facet 1 -- inner triangle (top)
			h->V[fid+0][0] = i*Np + j;
			h->V[fid+0][1] = (i+1)*Np + j;
			h->V[fid+0][2] = i*Np + (j+1)%Np;
			h->N[fid+0][0] = 0;
			h->N[fid+0][1] = 0;
			h->N[fid+0][2] = 1;
			h->A[fid+0] = l1*hT/4.0;
			h->C[fid+0][0] = (r1+hT/4.0)*c+cx;
			h->C[fid+0][1] = (r1+hT/4.0)*s+cy;
			h->C[fid+0][2] = z+cz;
			// facet 1 -- inner triangle (bottom)
			h->V[Ftop+fid+0][0] = (i+Nc)*Np + j;
			h->V[Ftop+fid+0][1] = (i+Nc+1)*Np + j;
			h->V[Ftop+fid+0][2] = (i+Nc)*Np + (j+1)%Np;
			h->N[Ftop+fid+0][0] = 0;
			h->N[Ftop+fid+0][1] = 0;
			h->N[Ftop+fid+0][2] = -1;
			h->A[Ftop+fid+0] = l1*hT/4.0;
			h->C[Ftop+fid+0][0] = (r1+hT/4.0)*c+cx;
			h->C[Ftop+fid+0][1] = (r1+hT/4.0)*s+cy;
			h->C[Ftop+fid+0][2] = -z+cz;
			
			// facet 2 -- outer triangle
			h->V[fid+1][0] = (i+1)*Np + j;
			h->V[fid+1][1] = (i+2)*Np + (j+1)%Np;
			h->V[fid+1][2] = (i+2)*Np + j;
			h->N[fid+1][0] = 0;
			h->N[fid+1][1] = 0;
			h->N[fid+1][2] = 1;
			h->A[fid+1] = l2*hT/4.0;
			h->C[fid+1][0] = (r1+3.0*hT/4.0)*c+cx;
			h->C[fid+1][1] = (r1+3.0*hT/4.0)*s+cy;
			h->C[fid+1][2] = z+cz;
			// facet 2 -- outer triangle (bottom)
			h->V[Ftop+fid+1][0] = (i+Nc+1)*Np + j;
			h->V[Ftop+fid+1][1] = (i+Nc+2)*Np + (j+1)%Np;
			h->V[Ftop+fid+1][2] = (i+Nc+2)*Np + j;
			h->N[Ftop+fid+1][0] = 0;
			h->N[Ftop+fid+1][1] = 0;
			h->N[Ftop+fid+1][2] = -1;
			h->A[Ftop+fid+1] = l2*hT/4.0;
			h->C[Ftop+fid+1][0] = (r1+3.0*hT/4.0)*c+cx;
			h->C[Ftop+fid+1][1] = (r1+3.0*hT/4.0)*s+cy;
			h->C[Ftop+fid+1][2] = -z+cz;
			
			// facet 3 -- side 1 triangle
			h->V[fid+2][0] = (i+1)*Np + j;
			h->V[fid+2][1] = i*Np + (j+1)%Np;
			h->V[fid+2][2] = (i+2)*Np + (j+1)%Np;
			h->N[fid+2][0] = 0;
			h->N[fid+2][1] = 0;
			h->N[fid+2][2] = 1;
			h->A[fid+2] = (l1+l2)*hT/8.0;
			h->C[fid+2][0] = (r1+hT/2.0)*cos(th+dth/4.0)+cx;
			h->C[fid+2][1] = (r1+hT/2.0)*sin(th+dth/4.0)+cy;
			h->C[fid+2][2] = z+cz;
			// facet 3 -- side 1 triangle (bottom)
			h->V[Ftop+fid+2][0] = (i+Nc+1)*Np + j;
			h->V[Ftop+fid+2][1] = (i+Nc)*Np + (j+1)%Np;
			h->V[Ftop+fid+2][2] = (i+Nc+2)*Np + (j+1)%Np;
			h->N[Ftop+fid+2][0] = 0;
			h->N[Ftop+fid+2][1] = 0;
			h->N[Ftop+fid+2][2] = -1;
			h->A[Ftop+fid+2] = (l1+l2)*hT/8.0;
			h->C[Ftop+fid+2][0] = (r1+hT/2.0)*cos(th+dth/4.0)+cx;
			h->C[Ftop+fid+2][1] = (r1+hT/2.0)*sin(th+dth/4.0)+cy;
			h->C[Ftop+fid+2][2] = -z+cz;
			
			// facet 4 -- side 2 triangle
			h->V[fid+3][0] = i*Np + j;
			h->V[fid+3][1] = (i+1)*Np + j;
			h->V[fid+3][2] = (i+2)*Np + j;
			h->N[fid+3][0] = 0;
			h->N[fid+3][1] = 0;
			h->N[fid+3][2] = 1;
			h->A[fid+3] = (l1+l2)*hT/8.0;
			h->C[fid+3][0] = (r1+hT/2.0)*cos(th-dth/4.0)+cx;
			h->C[fid+3][1] = (r1+hT/2.0)*sin(th-dth/4.0)+cy;
			h->C[fid+3][2] = z+cz;
			// facet 4 -- side 2 triangle (bottom)
			h->V[Ftop+fid+3][0] = (i+Nc)*Np + j;
			h->V[Ftop+fid+3][1] = (i+Nc+1)*Np + j;
			h->V[Ftop+fid+3][2] = (i+Nc+2)*Np + j;
			h->N[Ftop+fid+3][0] = 0;
			h->N[Ftop+fid+3][1] = 0;
			h->N[Ftop+fid+3][2] = -1;
			h->A[Ftop+fid+3] = (l1+l2)*hT/8.0;
			h->C[Ftop+fid+3][0] = (r1+hT/2.0)*cos(th-dth/4.0)+cx;
			h->C[Ftop+fid+3][1] = (r1+hT/2.0)*sin(th-dth/4.0)+cy;
			h->C[Ftop+fid+3][2] = -z+cz;
		}
	}

	
	
/* points along the sides */
	for(i=0; i<1; i++)				// main circles
	{
		for(j=0; j<Np; j++)			// loop over angles in this circle
		{
			pid = 2*Ntop+i*Np+j;
			th = j*dth;
			h->P[pid][0] = Ro*cos(th)+cx;		// assign x value
			h->P[pid][1] = Ro*sin(th)+cy;		// assign y value
			h->P[pid][2] = 0+cz;						// assign z value
		}
	}
	
/* facets along the sides */	
	pid = Np*(Nc-1);	
	for(i=0;i<Np;i++)
	{
		th = i*dth;
		fid = 2*Ftop+4*i;
		// top triangle top row
		h->V[fid][0] = 2*Ntop + i;
		h->V[fid][1] = pid+(i+1)%Np;
		h->V[fid][2] = pid+i;
		h->N[fid][0] = cos(th+dth/4.0);
		h->N[fid][1] = sin(th+dth/4.0);
		h->N[fid][2] = 0;
		h->A[fid] = H/4.0*Ro*sin(dth/2.0);
		h->C[fid][0] = Ro*cos(th+dth/4.0)+cx;
		h->C[fid][1] = Ro*sin(th+dth/4.0)+cy;
		h->C[fid][2] = 3.0*H/8.0+cz;
		
		// bottom triangle top row
		h->V[fid+1][0] = 2*Ntop + i;
		h->V[fid+1][1] = pid+(i+1)%Np;
		h->V[fid+1][2] = 2*Ntop + (i+1)%Np;
		h->N[fid+1][0] = cos(th+3.0*dth/4.0);
		h->N[fid+1][1] = sin(th+3.0*dth/4.0);
		h->N[fid+1][2] = 0;
		h->A[fid+1] = H/4.0*Ro*sin(dth/2.0);
		h->C[fid+1][0] = Ro*cos(th+3.0*dth/4.0)+cx;
		h->C[fid+1][1] = Ro*sin(th+3.0*dth/4.0)+cy;
		h->C[fid+1][2] = 1.0*H/8.0+cz;
		
		// top triangle bottom row
		h->V[fid+2][0] = pid+Ntop+i;
		h->V[fid+2][1] = 2*Ntop + (i+1)%Np;
		h->V[fid+2][2] = 2*Ntop + i;
		h->N[fid+2][0] = cos(th+dth/4.0);
		h->N[fid+2][1] = sin(th+dth/4.0);
		h->N[fid+2][2] = 0;
		h->A[fid+2] = H/4.0*Ro*sin(dth/2.0);
		h->C[fid+2][0] = Ro*cos(th+dth/4.0)+cx;
		h->C[fid+2][1] = Ro*sin(th+dth/4.0)+cy;
		h->C[fid+2][2] = -1.0*H/8.0+cz;
		
		// bottom triangle bottom row
		h->V[fid+3][0] = pid+Ntop+i;
		h->V[fid+3][1] = pid+Ntop+(i+1)%Np;
		h->V[fid+3][2] = 2*Ntop + (i+1)%Np;
		h->N[fid+3][0] = cos(th+3.0*dth/4.0);
		h->N[fid+3][1] = sin(th+3.0*dth/4.0);
		h->N[fid+3][2] = 0;
		h->A[fid+3] = H/4.0*Ro*sin(dth/2.0);
		h->C[fid+3][0] = Ro*cos(th+3.0*dth/4.0)+cx;
		h->C[fid+3][1] = Ro*sin(th+3.0*dth/4.0)+cy;
		h->C[fid+3][2] = -3.0*H/8.0+cz;			
		
	}
	
h->nF = Ftot;
h->nP = Ntot;
h->nV = Ntot;
h->dim = 3;

/*
// write to file for testing
FILE *fp = fopen("test1.dat","w");
FILE *fp1 = fopen("test2.dat","w");
for(i=0;i<2*Nc+1;i++)
{
	for(j=0; j<Np; j++)
		fprintf(fp,"%.5f,%.5f,%.5f\n",h->P[i*Np+j][0],h->P[i*Np+j][1],h->P[i*Np+j][2]);
}

for(i=0;i<Ftot;i++)
	fprintf(fp1,"%u,%u,%u\n",h->V[i][0],h->V[i][1],h->V[i][2]);



fclose(fp);
fclose(fp1);
*/
}


/*------------------------------------------------------------	
 function: get_order
 purpose: determine order of object projection
 inputs: array of order values (to fill in), array of cm values
		 inclination angle ia
 ------------------------------------------------------------*/
void get_order(int *order, double cm[][3], double ia)
{
	int Oid;

	if(cm[0][2] < cm[1][2])		// if RL is lower than Disc system
	{
		order[0] = 0;			// project RL first
		Oid = 1;				// start projecting others next
	}
	else
	{
		order[5] = 0;			// project RL last
		Oid = 0;
	}

	if(ia <= PI/2.0 || ia > 3.0*PI/2.0)		// then we see more of the bottom
	{
		if(cm[3][2] < cm[5][2])		// if bottom half of WD is lower than HS
		{
			order[Oid] = 3;			// next in line is bottom of WD
			order[Oid+1] = 5;		// then bottom of HS
			order[Oid+2] = 1;		// then disc
			order[Oid+3] = 2;		// then top of WD
			order[Oid+4] = 4;		// lastly top of HS
		}
		else
		{
			order[Oid] = 5;			// next in line is bottom of HS
			order[Oid+1] = 3;		// then bottom of WD
			order[Oid+2] = 1;		// then disc
			order[Oid+3] = 4;		// then top of HS
			order[Oid+4] = 2;		// lastly top of WD
		}		
	}
	else							// if we see more of the top of the disc
	{
		if(cm[2][2] < cm[4][2])		// if top half of WD is lower than HS
		{
			order[Oid] = 2;			// next in line is top of WD
			order[Oid+1] = 4;		// then top of HS
			order[Oid+2] = 1;		// then disc
			order[Oid+3] = 3;		// then bottom of WD
			order[Oid+4] = 5;		// lastly bottom of HS
		}
		else
		{
			order[Oid] = 4;			// next in line is top of HS
			order[Oid+1] = 2;		// then top of WD
			order[Oid+2] = 1;		// then disc
			order[Oid+3] = 5;		// then bottom of HS
			order[Oid+4] = 3;		// lastly bottom of WD
		}	
	}

}


