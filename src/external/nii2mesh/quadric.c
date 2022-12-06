/////////////////////////////////////////////
//
// Mesh Simplification Tutorial
//
// (C) by Sven Forstmann in 2014
//
// License : MIT
// http://opensource.org/licenses/MIT
//
//https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
//
// 5/2016: Chris Rorden created minimal version for OSX/Linux/Windows compile
// 1/2022: Chris Rorden ported from C++ to pure C

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _MSC_VER

#else
 #include <unistd.h>
#endif
#include "meshtypes.h"
#include "quadric.h"
#include <float.h> //FLT_EPSILON, DBL_EPSILON

typedef double TSymetricMatrix[10];
struct TRef{
	int tid,tvertex;
};
struct TVertex { vec3d p;int tstart,tcount;TSymetricMatrix q;int border;};
struct TTriangle{
	int v[3];
	double err[4];
	bool dirty, deleted;
	vec3d n;
};

static void symMat1(TSymetricMatrix ret, double c){
	for (int i = 0; i < 10; i++)
		ret[i] = c;
} // symMat()

static void symMat4(TSymetricMatrix ret, double a,double b,double c,double d){
	ret[0] = a*a; ret[1] = a*b; ret[2] = a*c; ret[3] = a*d;
	ret[4] = b*b; ret[5] = b*c; ret[6] = b*d;
	ret[7] = c*c; ret[8] = c*d;
	ret[9] = d*d;
}// symMat2()

static void symMat10(TSymetricMatrix ret, double m11, double m12, double m13, double m14, double m22, double m23, double m24, double m33, double m34, double m44){
	ret[0] = m11; ret[1] = m12; ret[2] = m13; ret[3] = m14;
	ret[4] = m22; ret[5] = m23; ret[6] = m24;
	ret[7] = m33; ret[8] = m34;
	ret[9] = m44;
} // symMat3()

static void symMatAdd(TSymetricMatrix ret, TSymetricMatrix n, TSymetricMatrix m) {
	symMat10(ret, n[0]+m[0], n[1]+m[1], n[2]+m[2], n[3]+m[3], n[4]+m[4],
	n[5]+m[5], n[6]+m[6], n[7]+m[7], n[8]+m[8], n[9]+m[9]);
} // symMatAdd()

static double symMatDet(TSymetricMatrix m, int a11, int a12, int a13, int a21, int a22, int a23, int a31, int a32, int a33) {
	return m[a11]*m[a22]*m[a33] + m[a13]*m[a21]*m[a32] + m[a12]*m[a23]*m[a31]
	- m[a13]*m[a22]*m[a31] - m[a11]*m[a23]*m[a32]- m[a12]*m[a21]*m[a33];
} // symMatDet()

static vec3d ptf(double x, double y, double z) {
	return (vec3d){.x = x, .y = y, .z = z};
}// ptf()

static vec3d vCross(vec3d v1, vec3d v2) { //cross-product
	return ptf(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x);
}

static vec3d vSum(vec3d a, vec3d b){ //add two vectors
	return ptf(a.x+b.x, a.y+b.y, a.z+b.z);
}

static vec3d vSubtract(vec3d a, vec3d b){
	return ptf(a.x-b.x, a.y-b.y, a.z-b.z);
}

static void vNormalize(vec3d *v){ //make vector unit length
	double len = sqrt( (v->x*v->x) + (v->y*v->y) + (v->z*v->z));
	if (len <= 0) len = 0.001;
	v->x = v->x / len;
	v->y = v->y / len;
	v->z = v->z / len;
}

static double vDot (vec3d a,vec3d b){ //dot product
	return a.x*b.x + a.y*b.y + a.z*b.z;
} // vDot()

static vec3d vMult(vec3d a, double v){ //multiply
	return ptf(a.x*v, a.y*v, a.z*v);
} // vMult()

static double vertex_error(TSymetricMatrix q, double x, double y, double z){
	return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y
	+ 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
} // vertex_error()

static double calculate_error(int id_v1, int id_v2, vec3d *p_result, struct TVertex vertices[]) {
	TSymetricMatrix q;
	symMatAdd(q, vertices[id_v1].q, vertices[id_v2].q);
	int border = vertices[id_v1].border + vertices[id_v2].border;
	double det = symMatDet(q, 0, 1, 2, 1, 4, 5, 2, 5, 7);
	if (( det != 0.0) && ( border == 0)) {
		// q_delta is invertible
		p_result->x = -1.0/det*(symMatDet(q,1, 2, 3, 4, 5, 6, 5, 7 , 8));  // vx = A41/det(q_delta)
		p_result->y =  1.0/det*(symMatDet(q,0, 2, 3, 1, 5, 6, 2, 7 , 8));  // vy = A42/det(q_delta)
		p_result->z = -1.0/det*(symMatDet(q,0, 1, 3, 1, 4, 6, 2, 5,  8));  // vz = A43/det(q_delta)
		return vertex_error(q, p_result->x, p_result->y, p_result->z);
		return 1.0;
	}
	// det = 0 -> try to find best result
	vec3d p1 = vertices[id_v1].p;
	vec3d p2 = vertices[id_v2].p;
	vec3d p3 = vMult(vSum(p1, p2), 0.5);
	double error1 = vertex_error(q, p1.x,p1.y,p1.z);
	double error2 = vertex_error(q, p2.x,p2.y,p2.z);
	double error3 = vertex_error(q, p3.x,p3.y,p3.z);
	double error = fmin(error1, fmin(error2, error3));
	if (error1 == error) *p_result = p1;
	if (error2 == error) *p_result = p2;
	if (error3 == error) *p_result = p3;
	return error;
}

#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
#define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

static void update_mesh(int iteration, struct TTriangle triangles[], struct TVertex vertices[], struct TRef refs[], int *nrefs, int* nTri, int* nVert) {
	if (iteration>0) { // compact triangles
		int dst = 0;
		for (int i = 0; i < *nTri; i++) {
			if(!triangles[i].deleted) {
				triangles[dst] = triangles[i];
				dst = dst + 1;
			}//if not deleted
		} //for each triangle
		*nTri = dst;
		//realloc(triangles, dst); //<- we never resize
	} //if iteration > 0
	loopi(0,*nVert) {
		vertices[i].tstart=0;
		vertices[i].tcount=0;
	}
	loopi(0,*nTri)
		loopj(0,3) vertices[triangles[i].v[j]].tcount++;
	int tstart=0;
	loopi(0,*nVert) {
		vertices[i].tstart=tstart;
		tstart+=vertices[i].tcount;
		vertices[i].tcount=0;
	}
	int maxRef = 0;
	loopi(0,*nTri) {
		struct TTriangle *t=&triangles[i];
		loopj(0,3) {
			struct TVertex* v=&vertices[t->v[j]];
			maxRef = MAX(v->tstart+v->tcount, maxRef);
		}
	}
	*nrefs = maxRef + 1; //0..maxRef, so maxRef+1 items
	loopi(0,*nTri) {
		struct TTriangle *t=&triangles[i];
		loopj(0,3) {
			struct TVertex* v=&vertices[t->v[j]];
			refs[v->tstart+v->tcount].tid=i;
			refs[v->tstart+v->tcount].tvertex=j;
			v->tcount++;
		}
	}
	if( iteration != 0 ) return;
	// Init Quadrics by Plane & Edge Errors
	//
	// required at the beginning ( iteration == 0 )
	// recomputing during the simplification is not required,
	// but mostly improves the result for closed meshes
	//
	// Identify boundary : vertices[].border=0,1
	//std::vector<int> vcount,vids;
	loopi(0,*nVert)
		vertices[i].border=0;
	int *vids = (int *)malloc(*nVert * sizeof(int));
	int *vcount = (int *)malloc(*nVert * sizeof(int));
	loopi(0,*nVert) {
		int nvcount = 0;
		struct TVertex* v=&vertices[i];
		loopj(0,v->tcount) {
			int k=refs[v->tstart+j].tid;
			struct TTriangle *t=&triangles[k];
			loopk(0,3) {
				int ofs=0,id=t->v[k];
				while(ofs < nvcount) {
					if(vids[ofs]==id)break;
					ofs++;
				}
				if(ofs == nvcount) {
					vcount[nvcount] = 1;
					vids[nvcount] = id;
					nvcount++;
				}
				else
					vcount[ofs]++;
			}
		}
		loopj(0,nvcount) if(vcount[j]==1)
			vertices[vids[j]].border=1;
	}
	free(vcount);
	free(vids);
	//initialize errors
	loopi(0,*nVert)
		symMat1(vertices[i].q, 0.0);
	loopi(0,*nTri) {
		struct TTriangle *t=&triangles[i];
		vec3d n,p[3];
		loopj(0,3) p[j]=vertices[t->v[j]].p;
		n = vCross(vSubtract(p[1],p[0]),vSubtract(p[2],p[0]));
		vNormalize(&n);
		t->n=n;
		loopj(0,3) {
			TSymetricMatrix q;
			symMat4(q, n.x,n.y,n.z,-vDot(n,p[0]));
			symMatAdd(vertices[t->v[j]].q, vertices[t->v[j]].q, q);
		}
	}
	loopi(0,*nTri) {
		// Calc Edge Error
		struct TTriangle *t=&triangles[i];
		vec3d p;
		loopj(0,3) t->err[j]=calculate_error(t->v[j],t->v[(j+1)%3],&p, vertices);
		t->err[3]=fmin(t->err[0],fmin(t->err[1],t->err[2]));
	}
}

static void compact_mesh(struct TTriangle triangles[], struct TVertex vertices[], int* nTri, int* nVert){
		int dst=0;
		loopi(0,*nVert)
			vertices[i].tcount=0;
		loopi(0,*nTri) {
			if(!triangles[i].deleted){
				struct TTriangle t=triangles[i];
				triangles[dst++]=t;
				loopj(0,3)vertices[t.v[j]].tcount=1;
			}
		}
		* nTri = dst;
		dst=0;
		loopi(0, *nVert) {
			if(vertices[i].tcount) {
				vertices[i].tstart=dst;
				vertices[dst].p=vertices[i].p;
				dst++;
			}
		}
		loopi(0,*nTri){
			struct TTriangle *t=&triangles[i];
			loopj(0,3)t->v[j]=vertices[t->v[j]].tstart;
		}
		* nVert = dst;
}

static void update_triangles(int i0, struct TVertex* v, bool *deleted, int* deleted_triangles, struct TTriangle triangles[], struct TRef refs[], struct TVertex vertices[], int * nrefs){
	vec3d p;
	loopk(0,v->tcount) {
		struct TRef r=refs[v->tstart+k];
		struct TTriangle *t=&triangles[r.tid];
		if(t->deleted)continue;
		if(deleted[k]) {
			t->deleted=1;
			*deleted_triangles = *deleted_triangles + 1;
			continue;
		}
		t->v[r.tvertex]=i0;
		t->dirty=1;
		t->err[0]=calculate_error(t->v[0],t->v[1],&p, vertices);
		t->err[1]=calculate_error(t->v[1],t->v[2],&p, vertices);
		t->err[2]=calculate_error(t->v[2],t->v[0],&p, vertices);
		t->err[3]=fmin(t->err[0],fmin(t->err[1],t->err[2]));
		refs[*nrefs] = r;
		*nrefs = *nrefs + 1;
	}
}

static bool flipped(vec3d p,int i0,int i1,struct TVertex v0, struct TVertex v1,bool *deleted, struct TTriangle triangles[], struct TRef refs[], struct TVertex vertices[]) {
	loopk(0,v0.tcount) {
		struct TTriangle *t=&triangles[refs[v0.tstart+k].tid];
		if(t->deleted)continue;
		int s=refs[v0.tstart+k].tvertex;
		int id1=t->v[(s+1)%3];
		int id2=t->v[(s+2)%3];
		if(id1==i1 || id2==i1) {// delete ?
			deleted[k]=1;
			continue;
		}
		vec3d d1 = vSubtract(vertices[id1].p, p);
		vNormalize(&d1);
		vec3d d2 = vSubtract(vertices[id2].p, p);
		vNormalize(&d2);
		if(fabs(vDot(d1, d2))>0.999) return true;
		vec3d n = vCross(d1,d2);
		vNormalize(&n);
		deleted[k]=0;
		if (vDot(n, t->n)<0.2) return true;
	}
	return false;
}

static void laplacian_smooth(vec3d *verts, vec3i *tris, int nvert, int ntri) {
	vec3d* sum = (vec3d*) malloc(nvert * sizeof(vec3d));
	memset(sum, 0, nvert * sizeof(vec3d));
	int* num = (int*) malloc(nvert * sizeof(int));
	memset(num, 0, nvert * sizeof(int));
	loopi(0, ntri) {
		//each point of a triangle has two neighbors:
		int p0 = tris[i].x;
		int p1 = tris[i].y;
		int p2 = tris[i].z;
		num[p0] += 2;
		sum[p0] = vSum(sum[p0], vSum(verts[p1], verts[p2]));
		num[p1] += 2;
		sum[p1] = vSum(sum[p1], vSum(verts[p0], verts[p2]));
		num[p2] += 2;
		sum[p2] = vSum(sum[p2], vSum(verts[p0], verts[p1]));
	}
	loopi(0, nvert) { //mean location of neighbors
		if (num[i] <= 0) continue;
		verts[i].x = sum[i].x / num[i];
		verts[i].y = sum[i].y / num[i];
		verts[i].z = sum[i].z / num[i];
	}
	free(sum);
	free(num);
}

void laplacian_smoothHC(vec3d *verts, vec3i *tris, int nvert, int ntri, double alpha, double beta, int iter, bool lockEdges) {
	// Laplacian smooth with Humphrey’s Classes to preserve volume: https://doi.org/10.1111/1467-8659.00334
	//  trimesh.smoothing.filter_humphrey(mesh, alpha=0.1, beta=0.5, iterations=10) https://trimsh.org/trimesh.smoothing.html
	double alpha1 = 1.0 - alpha;
	double beta1 = 1.0 - beta;
	vec3d* p = (vec3d*) malloc(nvert * sizeof(vec3d));
	vec3d* q = (vec3d*) malloc(nvert * sizeof(vec3d));
	vec3d* b = (vec3d*) malloc(nvert * sizeof(vec3d));
	memcpy(p, verts, nvert * sizeof(vec3d)); //dst,src,n
	loopj(0,iter) {
		memcpy(q, p, nvert * sizeof(vec3d)); //dst,src,n
		laplacian_smooth(p, tris, nvert, ntri);
		loopi(0,nvert)
			b[i] = vSubtract(p[i], vSum(vMult(verts[i], alpha), vMult(q[i], alpha1)));
		memcpy(q, b, nvert * sizeof(vec3d));
		laplacian_smooth(q, tris, nvert, ntri);
		loopi(0,nvert)
			p[i] = vSubtract(p[i], vSum(vMult(b[i], beta), vMult(q[i], beta1)));
	}
	free(q);
	free(b);
	if (!lockEdges) {
		loopi(0,nvert)
			verts[i] = p[i];
		free(p);
		return;
	}
	//find border edges:
	struct TVertex* vertices = (struct TVertex*) malloc(nvert * sizeof(struct TVertex));
	loopi(0, nvert)
		vertices[i].p = verts[i];
	struct TTriangle* triangles = (struct TTriangle*) malloc(ntri * sizeof(struct TTriangle));
	loopi(0, ntri) {
		triangles[i].deleted=0;
		triangles[i].v[0] = tris[i].x;
		triangles[i].v[1] = tris[i].y;
		triangles[i].v[2] = tris[i].z;
	}
	int nref = 0;
	struct TRef* refs = (struct TRef*) malloc(ntri * 6 * sizeof(struct TRef)); //overprovision initially *ntri * 3, allow room for growth
	int ntriOK = ntri;
	int vertex_count = nvert;
	update_mesh(0, triangles, vertices, refs, &nref, &ntriOK, &vertex_count);
	free(triangles);
	free(refs);
	loopi(0,nvert) {
		if (vertices[i].border) continue;
		verts[i] = p[i];
	}
	free(vertices);
	free(p);
}

void quadric_simplify_mesh(vec3d **vs, vec3i **ts, int* nvert, int *ntri, int target_count, double agressiveness, bool verbose, bool finishLossless) {
	// init: load vertices
	vec3d *verts = *vs;
	struct TVertex* vertices = (struct TVertex*) malloc(*nvert * sizeof(struct TVertex));
	loopi(0,*nvert)
		vertices[i].p = verts[i];
	free(*vs);
	//init: load triangle faces
	vec3i *tris = *ts;
	struct TTriangle* triangles = (struct TTriangle*) malloc(*ntri * sizeof(struct TTriangle));
	loopi(0,*ntri) {
		triangles[i].deleted=0;
		triangles[i].v[0] = tris[i].x;
		triangles[i].v[1] = tris[i].y;
		triangles[i].v[2] = tris[i].z;
	}
	free(*ts);
	int nref = 0;
	struct TRef* refs = (struct TRef*) malloc(*ntri * 9 * sizeof(struct TRef)); //overprovision initially *ntri * 3, allow room for growth
	//init other structures
	bool* deleted0 = (bool*) malloc(*ntri * 3 * sizeof(bool)); //overprovision so we never need to realloc
	bool* deleted1 = (bool*) malloc(*ntri * 3 * sizeof(bool)); //overprovision so we never need to realloc
	// main iteration loop
	int deleted_triangles=0;
	int vertex_count = *nvert;
	int triangle_count=*ntri;
	int ntriOK = triangle_count;
	int max_iter = 100;
	bool lossy = true;
	double threshold = DBL_EPSILON;
	if (target_count >= ntriOK) {
		lossy = false;
		max_iter = 1000;
	}
	int iterationStartCount = 0;
	for (int iteration = 0; iteration < max_iter; iteration ++) {
		if ((lossy) && ((triangle_count-deleted_triangles)<=target_count)) {
			if (!finishLossless) break;
			lossy = false;
			threshold = DBL_EPSILON;
			max_iter = 1000;
		}
		if (lossy) {
			//lossy: update mesh once in a while
			if(iteration%5==0)
				update_mesh(iteration, triangles, vertices, refs, &nref, &ntriOK, &vertex_count);
			threshold = 0.000000001*pow((double)(iteration+3.0),agressiveness);
		} else {
			if (iterationStartCount == (triangle_count-deleted_triangles)) break;
			//lossless: update mesh constantly
			update_mesh(iteration, triangles, vertices, refs, &nref, &ntriOK, &vertex_count);
		}
		iterationStartCount = triangle_count-deleted_triangles;
		// clear dirty flag
		loopi(0,ntriOK)
			triangles[i].dirty=0;
		//
		// All triangles with edges below the threshold will be removed
		//
		// The following numbers works well for most models.
		// If it does not, try to adjust the 3 parameters
		//
		// target number of triangles reached ? Then break
		if ((verbose) && (iteration%5==0))
			printf(" iteration %d - triangles %d threshold %g\n",iteration,triangle_count-deleted_triangles, threshold);
		// remove vertices & mark deleted triangles
		loopi(0,ntriOK) {
			struct TTriangle *t=&triangles[i];
			if(t->err[3]>threshold) continue;
			if(t->deleted) continue;
			if(t->dirty) continue;
			loopj(0,3)if(t->err[j]<threshold) {
				int i0=t->v[j];
				struct TVertex *v0 = &vertices[i0];
				int i1=t->v[(j+1)%3];
				struct TVertex *v1 = &vertices[i1];
				if(v0->border != v1->border)  continue;
				// Compute vertex to collapse to
				vec3d p;
				calculate_error(i0,i1,&p, vertices);
				if( flipped(p,i0,i1,*v0,*v1,deleted0, triangles, refs, vertices) ) continue;
				if( flipped(p,i1,i0,*v1,*v0,deleted1, triangles, refs, vertices) ) continue;
				// not flipped, so remove edge
				v0->p=p;
				symMatAdd(v0->q, v1->q, v0->q); //v0.q=v1.q+v0.q;
				int tstart=nref;
				update_triangles(i0,v0,deleted0,&deleted_triangles, triangles, refs, vertices, &nref);
				update_triangles(i0,v1,deleted1,&deleted_triangles, triangles, refs, vertices, &nref);
				int tcount=nref-tstart;
				if(tcount<=v0->tcount) {
					// save ram
						if(tcount)memcpy(&refs[v0->tstart],&refs[tstart],tcount*sizeof(struct TRef));
				}
				else
					// append
					v0->tstart=tstart;
				v0->tcount=tcount;
				break;
			}
			// done?
			//if(triangle_count-deleted_triangles<=target_count) threshold = DBL_EPSILON;
			if((lossy) && ((triangle_count-deleted_triangles)<=target_count)) break;
		} //for each triangle
	} //for each iteration
	free(refs);
	free(deleted0);
	free(deleted1);
	triangle_count = ntriOK;
	// clean up mesh
	compact_mesh(triangles, vertices, &triangle_count, &vertex_count);
	*ntri = triangle_count;
	*ts= (vec3i *) malloc(*ntri * sizeof(vec3i));
	tris = *ts;
	loopi(0,*ntri)
		tris[i] = (vec3i){ .x = triangles[i].v[0], .y = triangles[i].v[1], .z = triangles[i].v[2]};
	*nvert = vertex_count;
	*vs = (vec3d *) malloc(*nvert * sizeof(vec3d));
	verts = *vs;
	loopi(0,*nvert)
		verts[i] = (vec3d){ .x = vertices[i].p.x, .y = vertices[i].p.y, .z = vertices[i].p.z};
	free(triangles);
	free(vertices);
} //quadric_simplify_mesh()
