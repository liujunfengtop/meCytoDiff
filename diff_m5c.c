/* diff_m5c.c

   Junfeng Liu, 2017
   
   This is an program for the differential analysis for m5c level according the file m5c_A.tsv and m5c_B.tsv.

   gcc -o diff_m5c diff_m5c.c -lm
   
   diff_m5c <m5c_consitionA.tsv> <m5c_conditionB.tsv> <pvalue>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>


int Diffm5c (char *m5cAf, char *m5cBf, char *pvalue);
double lagam(double x);
double lbgam(double a, double x);
double lcerf(double x);
double ligas(double a, double d, double x);

int main (int argc, char* argv[])
{
   int l;
   
   char m5cAf[1000], m5cBf[1000], pvalue[1000];

   if(argc>3) {
   	strcpy(m5cAf, argv[1]);
   	strcpy(m5cBf, argv[2]);
   	strcpy(pvalue, argv[3]);
   } else {
   	return -1;
   }

   l=Diffm5c(m5cAf, m5cBf, pvalue);
     
   return (l);
}

int Diffm5c (char *m5cAf, char *m5cBf, char *pvalue)
{
	/* j for chracter loop about each row of the file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read;
	
	
	double es_A_m5c, es_A_m5c_mean, es_A_m5c_var, es_B_m5c, es_B_m5c_mean, es_B_m5c_var;
	double es_A_pvalue, es_B_pvalue;
	
	/*ls[] for each row of the file;
	  temp[] for each column of the row from the file;
	*/
	char ls[5000], temp[1000], target_id[1000];

	
	FILE *fm5cA, *fm5cB, *fout;
	
	if((fm5cA=fopen(m5cAf, "r"))==NULL) {
   		return -1;
   }
   
  if((fm5cB=fopen(m5cBf, "r"))==NULL) {
   		return -1;
   }
   
  
	if((fout=fopen("diff_out.tsv", "w"))==NULL) {
   		return -1;
   }
  
  fprintf(fout,"target_id\testimated_A_m5c\tA_mean\tA_variance\tA_pvalue\testimated_B_m5c\tB_mean\tB_variance\tB_pvalue\n");
  
  fgets(ls,5000,fm5cA);
  fgets(ls,5000,fm5cB);
  
	while(fgets(ls,5000,fm5cA) != NULL) {
		es_A_m5c=0.0;
		es_A_m5c_mean=0.0;
		es_A_m5c_var=0.0;
		es_B_m5c=0.0;
		es_B_m5c_mean=0.0;
		es_B_m5c_var=0.0;
		es_A_pvalue=0.0;
		es_B_pvalue=0.0;
		k=0;
		j=0;
		read=0;
		memset(temp,0,sizeof(temp));
		while(read<4) {
			temp[k]=ls[j];
			if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
				if(read==0) {
					memset(target_id,0,sizeof(target_id));
					for(l=0;l<strlen(temp)-1;l++) {
						target_id[l]=temp[l];
					}
				}
				if(read==1) {
					es_A_m5c=atof(temp);
				}
				if(read==2) {
					es_A_m5c_mean=atof(temp);
				}
				if(read==3) {
					es_A_m5c_var=atof(temp);
				}
				memset(temp,0,sizeof(temp));
				read=read+1;
				k=0;
			} else {
				k=k+1;
			}
			j=j+1;
		}
		fgets(ls,5000,fm5cB);
		k=0;
		j=0;
		read=0;
		memset(temp,0,sizeof(temp));
		while(read<4) {
			temp[k]=ls[j];
			if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
				if(read==1) {
					es_B_m5c=atof(temp);
				}
				if(read==2) {
					es_B_m5c_mean=atof(temp);
				}
				if(read==3) {
					es_B_m5c_var=atof(temp);
				}
				memset(temp,0,sizeof(temp));
				read=read+1;
				k=0;
			} else {
				k=k+1;
			}
			j=j+1;
		}
		es_A_pvalue=ligas(es_B_m5c_mean,sqrt(es_B_m5c_var),es_A_m5c);
		es_B_pvalue=ligas(es_A_m5c_mean,sqrt(es_A_m5c_var),es_B_m5c);
		if(((es_A_pvalue<atof(pvalue)/2)||(es_A_pvalue>(1-atof(pvalue)/2)))||((es_B_pvalue<atof(pvalue)/2)||(es_B_pvalue>(1-atof(pvalue)/2)))) {
			fprintf(fout,"%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",target_id,es_A_m5c,es_A_m5c_mean,es_A_m5c_var,es_A_pvalue,es_B_m5c,es_B_m5c_mean,es_B_m5c_var,es_B_pvalue);
		}
	}
	
	
	fclose(fm5cA);
	fclose(fm5cB);
	fclose(fout);
	return 0;
}

double lagam(double x) {
	int i;
	double y, t, s, u;
	static double a[11]={0.0000677106, -0.0003442342, 0.0015397681, -0.0024467480, 0.0109736958, -0.0002109075, 0.0742379071, 0.0815782188, 0.4118402518, 0.4227843370, 1.0};
	if(x<=0.0) {
		printf("err * * <=0!\n");
		return(-1.0);
	}
	y=x;
	if(y<=1.0) {
		t=1.0/(y*(y+1.0));
		y=y+2.0;
	} else if(y<=2.0) {
		t=1.0/y;
		y=y+1.0;
	} else if(y<=3.0) {
		t=1.0;
	} else {
		t=1.0;
		while (y>3.0) {
			y=y-1.0;
			t=t*y;
		}
	}
	s=a[0];
	u=y-2.0;
	for(i=1;i<=10;i++) {
		s=s*u+a[i];
	}
	s=s*t;
	return(s);	
}


double lbgam(double a, double x) {
	int n;
	double p, q, d, s, s1, p0, q0, p1, q1, qq;
	if((a<=0.0)||(x<0.0)) {
		if(a<=0.0) {
			printf("err * * a<=0!\n");
		}
		if(x<0.0) {
			printf("err * * x<0!\n");
			return(-1.0);
		}
	}
	if(x+1.0==1.0) {
		return(0.0);
	}
	if(x>1.0e+35) {
		return(1.0);
	}
	q=log(x);
	q=a*q;
	qq=exp(q);
	if(x<1.0+a) {
		p=a;
		d=1.0/a;
		s=d;
		for(n=1;n<=100;n++) {
			p=1.0+p;
			d=d*x/p;
			s=s+d;
			if(fabs(d)<fabs(s)*1.0e-07) {
				s=s*exp(-x)*qq/lagam(a);
				return(s);
			}
		}
	} else {
		s=1.0/x;
		p0=0.0;
		p1=1.0;
		q0=1.0;
		q1=x;
		for(n=1;n<=100;n++) {
			p0=p1+(n-a)*p0;
			q0=q1+(n-a)*q0;
			p=x*p0+n*p1;
			q=x*q0+n*q1;
			if(fabs(q)+1.0!=1.0) {
				s1=p/q;
				p1=p;
				q1=q;
				if(fabs((s1-s)/s1)<1.0e-07) {
					s=s1*exp(-x)*qq/lagam(a);
					return(1.0-s);
				}
				s=s1;
			}
			p1=p;
			q1=q;
		}
	}
	printf("a too large!\n");
	s=1.0-s*exp(-x)*qq/lagam(a);
	return(s);
}

double lcerf(double x) {
	double y;
	if(x>=0.0) {
		y=lbgam(0.5,x*x);
	} else {
		y=-lbgam(0.5,x*x);
	}
	return(y);
}

double ligas(double a, double d, double x) {
	double y;
	if(d<=0.0) {
		d=1.0e-10;
	}
	y=0.5+0.5*lcerf((x-a)/(sqrt(2.0)*d));
	return(y);
}