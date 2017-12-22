#include <math.h>
#include <stdio.h>
#define SMALL (1e-10)
int last=0;
double dot(a,b,m)
double *a,*b;
{
	double sum;
	sum=0; while(--m>=0)sum+=*a++**b++;
	return sum;
}
double gram(a,v,n,m,rhs,normv)
double *a,*v,*normv;
double rhs;
{
	double d,dsq;
	int i,k,mpo;
	mpo=m+1;
	for(i=0; i<m; i++)normv[i]=v[i];
	normv[m]=rhs;
	for(k=0; k<n; k++) {
		d=dot(a+k*mpo,normv,m);
		for(i=0; i<mpo; i++)normv[i]-=d*a[k*mpo+i];
	}
	dsq=dot(normv,normv,m);
	if(dsq<SMALL)return 0;
	d=sqrt(dsq);
	for(i=0; i<mpo; i++)normv[i]/=d;
	return dsq;
}
randvec(v,m)
double *v;
{
	int i;
	double d;
	for(i=0; i<m; i++) {
		last=(last*101+55)&32767;
		v[i]=last-16384;
	}
	d=sqrt(dot(v,v,m));
	for(i=0; i<m; i++)v[i]/=d;
}
gsolve(a,n,m,ans)
double *a,*ans;
{
	int i,j,k,mpo;
	double d;
	mpo=m+1;
	for(k=n; k<m; k++) {
		do {
			randvec(a+k*mpo,m);
			d=gram(a,a+k*mpo,k,m,0.,a+k*mpo);
		} while(d<1./m);
	}
	for(i=0; i<m; i++) {
		ans[i]=0;
		for(j=0; j<n /*sic*/; j++)ans[i]+=a[j*mpo+i]*a[j*mpo+m];
	}
	return 0;
}
