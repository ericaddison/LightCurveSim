/*************************************
	File: hullLC.c
	Project: Light Curve Simulation (LC)
	Purpose: Code file for convex hull calculations
	Author: Eric Addison
	Original date: 23Apr13
*************************************/

/* This file defines functions declared in hullLC.h. This requires that the qhull library is installed, and appropriate linker flags are included in the comiler call.
	e.g. gcc -o bar -I/usr/local/include/libqhull -L/usr/local/lib foo.c -lqhullstatic -lm
*/

#include <stdio.h>
#include "hullLC.h"

/* function to call qhull with triangulation option, allocate hull, and extract useful info all in one call */
void call_qhull(hullLC *h)
{
		char flags[25];	// string to hold qhull options
		sprintf(flags, "qhull FA Qt Qc");	// FA option computes areas of facets
									 // Qt option produces triangulation instead
									// 		of merging coplanar facets 
	/* Error checking -- make sure proper fields in h are filled in */
		if(h->dim <= 1 || h->dim > 9)
		{	
			fprintf(stderr,"Error: Invalid dimension in hullLC struct in call_qhull\n");
			exit(CALL_QHULL_ERR_INV_DIM);
		}
		if(h->nP <= 0 || h->nP > 1e6)
		{	
			fprintf(stderr,"Error: Invalid number of points h->nP in call_qhull\n");
			exit(CALL_QHULL_ERR_INV_NP);
		}
		
		qh_new_qhull(h->dim, h->nP, h->Pq, 0, flags, NULL, stderr);
		// the NULL there outputs the summary of the Qhull call to NULL
		// to see the summary, change NULL to stdout or stderr	
		h->nF = qh num_facets;	
		h->nV = qh num_vertices;
		MallocHull(h);
		extract_qhull_CNA(h);
		extract_qhull_VP(h);
}

/* function to extract centers, normals, and areas of all facets in facet_list of most recent call to Qhull */
void extract_qhull_CNA(hullLC *h)
{
	facetT *facet;	// iterator used in FORALLfacets
	pointT *center; // mem to hold center valuem
	int i=0;	// array counter

	FORALLfacets {
		h->A[i] = (double)facet->f.area;		// area of facet i
		center = qh_getcentrum(facet);	// compute center
		for(int j=0; j<3; j++) {	
			h->C[i][j] = center[j];		// center of facet i
			h->N[i][j] = facet->normal[j];	// normal of facet i
		}		
		i++;
	}
}

/* function to extract the points of the hull and the vertices of each facet and store in arrayis P and V */
void extract_qhull_VP(hullLC *h)
{
	// list the id values of each vertex
	unsigned i=0,j;
	vertexT *vertex;	// iterator to walk through vertex_list
	FORALLvertices{
		// write current point to h->P array
		for(int k=0;k<h->dim;k++)		
			h->P[i][k] = vertex->point[k];
		//printf("Vertex %d: id = %u, (%.1f,%.1f,%.1f)\n",i,vertex->id,vertex->point[0],vertex->point[1],vertex->point[2]);
		vertex->id = i;
		i++;
	}

	// find the vertices for each facet
	facetT *facet;		// required iterator for FORALLfacets
	vertexT **vertexp;	// required iterator for FOREACHvertex (also *vertex above)
	h->V = malloc(qh num_facets*sizeof(int*));
	for(i=0;i<qh num_facets;i++)
		h->V[i] = malloc(3*sizeof(unsigned));
	
	i=0;
	FORALLfacets{
		//printf("%d: Facet id = %u, vertices = (",i,facet->id);
		j=0;
		FOREACHvertex_(facet->vertices){
			//printf("%u, ",vertex->id);
			h->V[i][j++] = vertex->id;
		}
		i++;
		//printf("\b\b)\n");
	}
}


/* function to extract the points of the 2D hull only and store them in an array *Px and *Py */
void extract_qhull2D_P(coordT *Px, coordT *Py)
{
	// list the id values of each vertex
	unsigned i=0;
	vertexT *vertex;	// iterator to walk through vertex_list
//	FILE *fp = fopen("idTest.dat","w");
	FORALLvertices{
//		fprintf(fp,"%d\n",vertex->id);
		// write current point to h->P array		
			Px[vertex->id] = vertex->point[0];
			Py[vertex->id] = vertex->point[1];
		i++;
	}
//	fclose(fp);
}

/* function to allocate memory needed for one hullLC struct */
void MallocHull(hullLC *h)
{
	// allocate memory for C,N,A
	h->C = (double **)malloc(h->nF*sizeof(double*));
	h->N = (double **)malloc(h->nF*sizeof(double*));
	h->V = (unsigned **)malloc(h->nF*sizeof(unsigned*));
	h->P = (double **)malloc(h->nV*sizeof(double*));
	h->A = (double *)malloc(h->nF*sizeof(double));
	for(int i=0;i<h->nF;i++)
	{
		h->C[i] = (double *)malloc(h->dim*sizeof(double));
		h->N[i] = (double *)malloc(h->dim*sizeof(double));
		h->V[i] = (unsigned *)malloc(h->dim*sizeof(unsigned));
	}
	
	for(int i=0;i<h->nV;i++)
		h->P[i] = (double *)malloc(h->dim*sizeof(double));

}

/* function to free memory used by hullLC arrays */
void FreeHulls(hullLC *h, int n)
{

	for (int j=0; j<n; j++)
		FreeHull(&(h[j]));
}


void FreeHull(hullLC *h)
{
	for(int i=0;i<h->nF;i++)
	{
		free(h->C[i]);
		free(h->N[i]);
		free(h->V[i]);
	}
	for(int i=0;i<h->nV;i++)
		free(h->P[i]);
	free(h->C);
	free(h->N);
	free(h->P);
	free(h->V);
	free(h->A);	
}

/* function to write hull data out to a file */
void write_hull_to_file(hullLC *h, const char* filenamebase)
{
	// build filenames
	char facetFile[80],vertexFile[80],infoFile[80];
	sprintf(facetFile,"%s.facets.dat",filenamebase);
	sprintf(vertexFile,"%s.vertex.dat",filenamebase);
	sprintf(infoFile,"Hull.info");

	// open files
	FILE *fpF = fopen(facetFile,"w");
	FILE *fpV = fopen(vertexFile,"w");

	// write to facetFile
	for(int i=0; i<h->nF; i++)
	{
		//printf("\niii = %d",i);
		fprintf(fpF,"%d, ",i);
		for(int j=0; j<h->dim; j++)
			fprintf(fpF,"%u, ",h->V[i][j]);
		for(int j=0; j<h->dim; j++)
			fprintf(fpF,"%.10f, ",h->C[i][j]);
		for(int j=0; j<h->dim; j++)
			fprintf(fpF,"%.10f, ",h->N[i][j]);
		fprintf(fpF,"%.10f",h->A[i]);
		if(i != h->nF-1)
			fprintf(fpF,"\n");
	}



	// write to vertexFile
	for(int i=0; i<h->nV; i++)
	{
		//printf("\ni = %d",i);
		fprintf(fpV,"%d, ",i);
		for(int j=0; j<h->dim; j++)
			fprintf(fpV,"%.10f, ",h->P[i][j]);
		if(i != h->nP-1)
			fprintf(fpV,"\n");
	}
	
	fclose(fpF);
	fclose(fpV);
}

// stolen and modified from function qh_printextremes_2d() in qhull->io.c
// orders points in a 2D convex hull in CC order
// saves ordered points in arrays px and py

void extract_qhull2D_P2(facetT *facetlist, setT *facets, double *px, double *py) {
	int numfacets, numridges, totneighbors, numcoplanars, numsimplicial, numtricoplanars;
	setT *vertices;
	facetT *facet, *startfacet, *nextfacet;
	vertexT *vertexA, *vertexB;
	
	qh_countfacets(facetlist, facets, 0, &numfacets, &numsimplicial,
				   &totneighbors, &numridges, &numcoplanars, &numtricoplanars); /* marks qh visit_id */
	vertices= qh_facetvertices(facetlist, facets, 0);
	//printf("%d\n", qh_setsize(vertices));
	qh_settempfree(&vertices);
	if (!numfacets)
		return;
	facet= startfacet= facetlist ? facetlist : SETfirstt_(facets, facetT);
	qh vertex_visit++;
	qh visit_id++;
	int i=0;
	do {
		if (facet->toporient ^ qh_ORIENTclock) {
			vertexA= SETfirstt_(facet->vertices, vertexT);
			vertexB= SETsecondt_(facet->vertices, vertexT);
			nextfacet= SETfirstt_(facet->neighbors, facetT);
		}else {
			vertexA= SETsecondt_(facet->vertices, vertexT);
			vertexB= SETfirstt_(facet->vertices, vertexT);
			nextfacet= SETsecondt_(facet->neighbors, facetT);
		}
		if (facet->visitid == qh visit_id) {
			printf("Qhull internal error (qh_printextremes_2d): loop in facet list.  facet %d nextfacet %d\n",
					   facet->id, nextfacet->id);
			qh_errexit2 (qh_ERRqhull, facet, nextfacet);
		}
		if (facet->visitid) {
			if (vertexA->visitid != qh vertex_visit) {
				vertexA->visitid= qh vertex_visit;
				px[i] = vertexA->point[0];
				py[i] = vertexA->point[1];
				i++;
				//printf("%d\n", qh_pointid(vertexA->point));
			}
			if (vertexB->visitid != qh vertex_visit) {
				vertexB->visitid= qh vertex_visit;
				px[i] = vertexB->point[0];
				py[i] = vertexB->point[1];
				i++;				
				//printf("%d\n", qh_pointid(vertexB->point));
			}
		}
		facet->visitid= qh visit_id;
		facet= nextfacet;
	}while (facet && facet != startfacet);
} /* printextremes_2d */


