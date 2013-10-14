/*************************************
	File: hullLC.h
	Project: Light Curve Simulation (LC)
	Purpose: Header file for convex hull calculations
	Author: Eric Addison
	Original date: 23Apr13
*************************************/

/* This file declares functions and data structures to interact with the
	qhull library. This requires that the qhull library is installed, and 
	appropriate linker flags are included in the comiler call.
	e.g. gcc -o bar -I/usr/local/include/libqhull -L/usr/local/lib foo.c -lqhullstatic -lm
*/

#ifndef HULLLC_H
#define HULLLC_H

/* qhull library inclusion */
#include <qhull_a.h>


/***************** data structures **********************/

/* struct to hold data relevent to one convex hull */
typedef struct {
	coordT *Pq;		// points list -- 1D array to feed to qhull (coordT is a real type)
	coordT **P;		// 2D array of the vertices from the hull
	coordT **C;		// center points of each facet
	coordT **N;		// normal vectors for each facet
	unsigned **V;	// array of vertices associated with each facet
	realT *A;		// area of each facet
	unsigned nF;	// number of facets in this hull
	unsigned nP;	// number of input points
	unsigned nV;	// number of vertices in the hull
	unsigned dim;	// dimension of the space
}hullLC;


/*************** hullLC.c functions decls ********************/

void call_qhull(hullLC *h);
void extract_qhull_CNA(hullLC *h);
void extract_qhull_VP(hullLC *h);
void extract_qhull2D_P(coordT *Px, coordT *Py);
void MallocHull(hullLC *h);
void FreeHulls(hullLC *h, int n);
void FreeHull(hullLC *h);
void write_hull_to_file(hullLC *h, const char* filenamebase);
void extract_qhull2D_P2(facetT *facetlist, setT *facets, double *px, double *py);


/*************** Error Codes ********************/
#define CALL_QHULL_ERR_INV_DIM 456
#define CALL_QHULL_ERR_INV_NP 457

#endif
