#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Q .809016994374947424102293 /* (1+sqrt(5))/4 */
#define P .309016994374947424102293
double gr[61*9]= {
	1,0,0,
	0,1,0,
	0,0,1,

	-.5,  P,  Q,
	 -P,  Q,-.5,
	 -Q,-.5, -P,

	-1, 0, 0,
	 0, 1, 0,
	 0, 0,-1
};
double vert[13*3];
double face[21*3]= { 1,1,1 };
double edge[31*3]= { 1,0,0 };
static int npts,nv,nf,ne;
icosagen()
{
	int i,j,k,n;
	n=3;
	for(i=1; i<n; i++) {
		if(n>60)break;
		for(j=1; j<i; j++) {
			mulmat(gr+i*9,gr+j*9,gr+n*9);
			for(k=0; k<n; k++)if(eq(gr+k*9,gr+n*9,9))break;
			if(k==n)n++;
		}
	}
	if(n!=60)fprintf(stderr,"icos botch\n");
	vert[0]=1./Q; vert[1]=2; vert[2]=0;
	fixpt(vert,12);
	fixpt(face,20);
	fixpt(edge,30);
}
mulmat(a,b,c)
double *a,*b,*c;
{
	int i,j,k;
	for(i=0; i<3; i++)for(j=0; j<3; j++) {
		c[i*3+j]=0;
		for(k=0; k<3; k++)c[i*3+j]+=a[i*3+k]*b[k*3+j];
	}
}
mulvec(a,b,c)
double *a,*b,*c;
{
	int i,j;
	for(i=0; i<3; i++) {
		c[i]=0;
		for(j=0; j<3; j++)c[i]+=a[i*3+j]*b[j];
	}
}
eq(a,b,nel)
double *a,*b;
{
	int i;
	for(i=0; i<nel; i++)if(a[i]<b[i]-1e-9||a[i]>b[i]+1e-9)return 0;
	return 1;
}
fixpt(a,n)
double *a;
{
	int i,j,k,m;
	double d;
	for(d=i=0; i<3; i++)d+=a[i]*a[i];
	for(i=0; i<3; i++)a[i]/=sqrt(d);
	for(m=k=1; k<60; k++) {
		if(m>n)break;
		mulvec(gr+k*9,a,a+m*3);
		for(j=0; j<m; j++)if(eq(a+m*3,a+j*3,3))break;
		if(j==m)m++;
	}
	if(m!=n) {
		fprintf(stderr,"fixpt botch\n");
		exit(1);
	}
}
icosnum(n)
{
	int i;
	npts=n;
	for(i=0; i<10000; i++) {
		for(nv=0; nv<=12; nv+=12)
		for(nf=0; nf<=20; nf+=20)
		for(ne=0; ne<=30; ne+=30)
			if(n==i*60+nv+nf+ne)return i;
	}
	fprintf(stderr,"bad icosahedral number\n");
	exit(1);
}
icosfill(a)
double *a;
{
	int i,j,k,n;
	n=icosnum(npts);
	icosagen();
	k=n*3;
	if(nv)for(i=0; i<12*3; i++)a[k++]=vert[i];
	if(nf)for(i=0; i<20*3; i++)a[k++]=face[i];
	if(ne)for(i=0; i<30*3; i++)a[k++]=edge[i];
	for(i=1; i<60; i++) {
		for(j=0; j<n; j++) {
			mulvec(gr+i*9,a+j*3,a+k);
			k+=3;
		}
	}
	return npts;
}
