/* sum_counts.c

   Junfeng Liu, 2017
   
   This is an program to compute the average estimated_counts according abundance_condition_order.tsv and the mean, variance of the estimated_counts according bootstrap_condition_order.tsv.

   gcc -o sum_counts sum_counts.c -lm
   
   sum_counts <abundance_consition_order.tsv> <bootstrap_condition_order.tsv> <number of sample> <number of bootstrap> <estimated_counts filterd value>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>


int SumCount (char *adf, char *bsf, char *sample_number, char *bs_number, char *filter_number);

int main (int argc, char* argv[])
{
   int l;
   
   char adf[1000], bsf[1000], sample_number[1000], bs_number[1000], filter_number[1000];

   if(argc>5) {
   	strcpy(adf, argv[1]);
   	strcpy(bsf, argv[2]);
   	strcpy(sample_number, argv[3]);
   	strcpy(bs_number, argv[4]);
   	strcpy(filter_number, argv[5]);
   } else {
   	return -1;
   }

   l=SumCount(adf, bsf, sample_number, bs_number, filter_number);
     
   return (l);
}

int SumCount (char *adf, char *bsf, char *sample_number, char *bs_number, char *filter_number)
{
	/* j for chracter loop about each row of the file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read, m1, m2;
	
	long n, m;
	
	double es_counts, es_mean, es_var;
	
	/*ls[] for each row of the file;
	  temp[] for each column of the row from the file;
	*/
	char ls[5000], temp[1000], target_id[1000];
	
	char **summ_ad, **summ_bs;
	
	
	
	FILE *fad, *fbs, *fout;
	
	if((fad=fopen(adf, "r"))==NULL) {
   		return -1;
   }
   
  if((fbs=fopen(bsf, "r"))==NULL) {
   		return -1;
   }
   
  
	if((fout=fopen("temp.tsv", "w"))==NULL) {
   		return -1;
   }
  
  fprintf(fout,"target_id\testimated_counts\tmean\tvariance\n");
  
	 n=0;
	 
	 while(fgets(ls,5000,fad) != NULL) {
	 	n=n+1;
	 }
	 
	 rewind(fad);
	 
	 n=n/atoi(sample_number);
	 
	 
	 
	 
	 summ_ad=(char**)malloc(atoi(sample_number)*sizeof(char*));
	 
	 summ_bs=(char**)malloc((atoi(sample_number)*atoi(bs_number))*sizeof(char*));
	 
	 
	 for(l=0;l<atoi(sample_number);l++) {
	 	summ_ad[l]=(char*)malloc(1000*sizeof(char));
	 }
	
	 for(l=0;l<(atoi(sample_number)*atoi(bs_number));l++) {
	 	summ_bs[l]=(char*)malloc(1000*sizeof(char));
	 }
	
   for(m=0;m<n;m++) {
   	es_counts=0.0;
   	es_mean=0.0;
   	es_var=0.0;
   	for(m1=0;m1<atoi(sample_number);m1++) {
   		fgets(ls,5000,fad);
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
   				if(read==3) {
   					for(l=0;l<strlen(temp)-1;l++) {
   						summ_ad[m1][l]=temp[l];
   					}
   				}
   				memset(temp,0,sizeof(temp));
   				read=read+1;
   				k=0;
   			} else {
   				k=k+1;
   			}
   			j=j+1;
   		}
   	}
   	for(m2=0;m2<(atoi(sample_number)*atoi(bs_number));m2++) {
   		fgets(ls,5000,fbs);
   		k=0;
   		j=0;
   		read=0;
   		memset(temp,0,sizeof(temp));
   		while(read<4) {
   			temp[k]=ls[j];
   			if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
   				if(read==3) {
   					for(l=0;l<strlen(temp)-1;l++) {
   						summ_bs[m2][l]=temp[l];
   					}
   				}
   				memset(temp,0,sizeof(temp));
   				read=read+1;
   				k=0;
   			} else {
   				k=k+1;
   			}
   			j=j+1;
   		}
   	}
   	for(l=0;l<atoi(sample_number);l++) {
   		es_counts=es_counts+atof(summ_ad[l]);
   		memset(summ_ad[l],0,sizeof(summ_ad[l]));
   	}
   	es_counts=es_counts/atoi(sample_number);
   	for(l=0;l<(atoi(sample_number)*atoi(bs_number));l++) {
   		es_mean=es_mean+atof(summ_bs[l]);
   	}
   	es_mean=es_mean/(atoi(sample_number)*atoi(bs_number));
   	for(l=0;l<(atoi(sample_number)*atoi(bs_number));l++) {
   		es_var=es_var+pow((atof(summ_bs[l])-es_mean),2);
   		memset(summ_bs[l],0,sizeof(summ_bs[l]));
   	}
   	es_var=es_var/((atoi(sample_number)*atoi(bs_number))-1);
   	if(es_counts>atoi(filter_number)) {
   		fprintf(fout,"%s\t%f\t%f\t%f\n",target_id,es_counts,es_mean,es_var);
   	}
   }
	

	free(summ_ad);
	free(summ_bs);
	
	fclose(fad);
	fclose(fbs);
	fclose(fout);
	return 0;
}

