#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Note: letters O and l missing */
char *table="0123456789ABCDEFGHIJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz+-*/";

int icos;
int tab[256];
bit()
{
	static int b,bitnum;
	int i;
	if(bitnum<=0) {
skip:
		i=getchar();
		if(i==EOF) {
			fprintf(stderr,"premature EOF\n");
			exit(1);
		}
		if(i=='\n'||i==' ')goto skip;
		b=tab[i];
		if(b<0) {
			fprintf(stderr,"illegal character '%c'",i);
			exit(1);
		}
		bitnum=6;
	}
	return (b&(1<<--bitnum))?1:0;
}
readint(fieldw)
{
	int i,j;
	for(i=j=0; i<fieldw; i++)j=j*2+bit();
	return j;
}
double readdouble(fieldw)
{
	int i;
	double d;
	fieldw--;
	d=0; for(i=0; i<fieldw; i++)d=d*2+bit();
	for(i=0; i<fieldw; i++)d/=2;
	if(bit())d=-d;
	return d;
}
textnum()
{
	int i,j;
	for(i=0; (j=getchar())!=EOF&&j>='0'&&j<='9'; )i=i*10+j-'0';
	if(j=='i') {
		icos=1;
		j=getchar();
	}
	if(j!=' ') {
		fprintf(stderr,"header number not followed by blank\n");
		exit(1);
	}
	return i;
}
readjump(vsize,vbreak)
{
	int i,j;
	if(vsize==0)return 1;
	else if(vsize==1)return bit()+1;
	i=readint(vsize-1);
	if(i<vbreak)j=i;
	else j=i*2+bit()-vbreak;
	j+=1;
	return j;
}
main()
{
	int i,j,k,m,n,nconst,njump,numberw,vsize,vbreak,jump,prev;
	double dsq,d,*a,*ans,*mat,dt,code;
	double gram();
	for(i=0; i<256; i++)tab[i]=-1;
	for(i=0; table[i]; i++)tab[table[i]]=i;
	tab['O']=tab['0']; tab['l']=tab['1']; /* O=0 l=1 */
	m=textnum(); n=textnum(); njump=textnum();
	numberw=textnum();
	a=(double*)malloc(((n*m)+(m)+((m+1)*m))*sizeof*a);
	if(a==0) {
		fprintf(stderr,"can't malloc\n");
		exit(1);
	}
	ans=a+n*m;
	mat=ans+m;
	if(njump)dt=readdouble(numberw);
	if(icos) {
		n=icosnum(n)+2;
		for(i=0; i<6; i++)a[i]=0;
		a[0]=a[5]=1;
	}
	for(k=icos?2:0; k<n; k++) {
		if(njump==0) {
			dsq=0;
			for(i=0; i<m-1; i++) {
				a[k*m+i]=readdouble(numberw);
				dsq+=a[k*m+i]*a[k*m+i];
			}
			if(dsq>1)dsq=1;
			d=sqrt(1.-dsq);
			a[k*m+i]=bit()?-d:d;
			continue;
		}
		prev=k;
		for(i=0; i<k&&i<m-1; i++) {
			j=prev-(k<m-1?k:m-1)+i+1;
			if(j>njump)j=njump;
			for(vsize=0; (1<<vsize)<j; vsize++);
			vbreak=(1<<vsize)-j;
			if(i==0&&bit()) {
				jump=1;
				code=dt;
			}
			else {
				jump=readjump(vsize,vbreak);
				code=bit()?dt:readdouble(numberw);
			}
			prev-=jump;
			dsq=gram(mat,a+m*prev,i,m,code,mat+i*(m+1));
			if(dsq==0) {
				fprintf(stderr,"problem with gram\n");
				exit(1);
			}
		}
		gsolve(mat,i,m,ans);
		dsq=0; for(i=0; i<m; i++)dsq+=ans[i]*ans[i];
		if(dsq>1)dsq=1; d=sqrt(1.-dsq);
		if(k&&bit())d=-d;
		for(i=0; i<m; i++)a[m*k+i]=ans[i]+d*mat[(m-1)*(m+1)+i];
	}
	if(icos)n=icosfill(a+=6);
	for(k=0; k<n; k++) {
		for(i=0; i<m; i++) {
			printf("%.12f\n",a[k*m+i]);
		}
	}
	exit(0);
}
