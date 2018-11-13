/* compute_m5c.c

   Junfeng Liu, 2017
   
   This is an program to compute the m5c level according the file summary_A.tsv and summary_methy_A.tsv.

   gcc -o compute_m5c compute_m5c.c -lm
   
   compute_m5c <summary_consition.tsv> <summary_methy_condition.tsv>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>


int Comm5c (char *scf, char *smcf);

int main (int argc, char* argv[])
{
   int l;
   
   char scf[1000], smcf[1000];

   if(argc>2) {
   	strcpy(scf, argv[1]);
   	strcpy(smcf, argv[2]);
   } else {
   	return -1;
   }

   l=Comm5c(scf, smcf);
     
   return (l);
}

int Comm5c (char *scf, char *smcf)
{
	/* j for chracter loop about each row of the file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read;
	
	
	double es_counts, es_mean, es_var, es_methy_counts, es_methy_mean, es_methy_var;
	double es_m5c, es_m5c_mean, es_m5c_var;
	
	/*ls[] for each row of the file;
	  temp[] for each column of the row from the file;
	*/
	char ls[5000], temp[1000], target_id[1000];

	
	FILE *fsc, *fsmc, *fout;
	
	if((fsc=fopen(scf, "r"))==NULL) {
   		return -1;
   }
   
  if((fsmc=fopen(smcf, "r"))==NULL) {
   		return -1;
   }
   
  
	if((fout=fopen("temp.tsv", "w"))==NULL) {
   		return -1;
   }
  
  fprintf(fout,"target_id\testimated_m5c\tmean\tvariance\n");
  
  fgets(ls,5000,fsc);
  fgets(ls,5000,fsmc);
  
	while(fgets(ls,5000,fsc) != NULL) {
		es_counts=0.0;
		es_mean=0.0;
		es_var=0.0;
		es_methy_counts=0.0;
		es_methy_mean=0.0;
		es_methy_var=0.0;
		es_m5c=0.0;
		es_m5c_mean=0.0;
		es_m5c_var=0.0;
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
					es_counts=atof(temp);
				}
				if(read==2) {
					es_mean=atof(temp);
				}
				if(read==3) {
					es_var=atof(temp);
				}
				memset(temp,0,sizeof(temp));
				read=read+1;
				k=0;
			} else {
				k=k+1;
			}
			j=j+1;
		}
		fgets(ls,5000,fsmc);
		k=0;
		j=0;
		read=0;
		memset(temp,0,sizeof(temp));
		while(read<4) {
			temp[k]=ls[j];
			if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
				if(read==1) {
					es_methy_counts=atof(temp);
				}
				if(read==2) {
					es_methy_mean=atof(temp);
				}
				if(read==3) {
					es_methy_var=atof(temp);
				}
				memset(temp,0,sizeof(temp));
				read=read+1;
				k=0;
			} else {
				k=k+1;
			}
			j=j+1;
		}
		if(es_counts>0.0&&es_mean>0.0) {
			es_m5c=es_methy_counts/es_counts;
			es_m5c_mean=(es_methy_mean/es_mean)+(es_methy_mean*es_var)/pow(es_mean,3);
			es_m5c_var=(es_methy_var/pow(es_mean,2))+(pow(es_methy_mean,2)*es_var)/pow(es_mean,4);
			fprintf(fout,"%s\t%f\t%f\t%f\n",target_id,es_m5c,es_m5c_mean,es_m5c_var);
		} else {
			fprintf(fout,"%s\t%s\t%s\t%s\n",target_id,"NA","NA","NA");
		}
	}
	
	
	fclose(fsc);
	fclose(fsmc);
	fclose(fout);
	return 0;
}

