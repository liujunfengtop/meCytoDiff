/* selmethy.c

   Junfeng Liu, 2017
   
   This is an program to select the methylated reads from RNA_Seq_peseudo data.

   gcc -o selmethy selmethy.c -lm
   
   selmethy <RNA_Seq_1.fq> <RNA_Seq_2.fq>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>


int SelMethy (char *fq1f, char *fq2f);

int main (int argc, char* argv[])
{
   int l;
   
   char fq1[1000], fq2[1000];

   if(argc>2) {
   	strcpy(fq1, argv[1]);
   	strcpy(fq2, argv[2]);
   } else {
   	return -1;
   }

   l=SelMethy(fq1,fq2);
     
   return (l);
}

int SelMethy (char *fq1f, char *fq2f)
{
	
	int l;
	
	
	/*ls[] for each row of the hits file;
	  temp[] for each column of the row;
	*/
	char ls[5000], ls1[5000];
	
	
	FILE *ffq1, *ffq2, *fout1, *fout2;
	
	if((ffq1=fopen(fq1f, "r"))==NULL) {
   		return -1;
  }
  
	if((ffq2=fopen(fq2f, "r"))==NULL) {
   		return -1;
  }
   
  
	if((fout1=fopen("methy_1.fq", "w"))==NULL) {
   		return -1;
  }
  
	if((fout2=fopen("methy_2.fq", "w"))==NULL) {
   		return -1;
  }
	
	 while(fgets(ls,5000,ffq1) != NULL) {
	 	fgets(ls1,5000,ffq2);
	 	if(strstr(ls,"methylated_liu")!=NULL) {
	 		fprintf(fout1,"%s",ls);
	 		fprintf(fout2,"%s",ls1);
	 		for(l=0;l<3;l++) {
	 			fgets(ls,5000,ffq1);
	 			fgets(ls1,5000,ffq2);
	 			fprintf(fout1,"%s",ls);
	 			fprintf(fout2,"%s",ls1);
	 		}
	 	}
	}
	fclose(ffq1);
	fclose(ffq2);
	fclose(fout1);
	fclose(fout2);
	return 0;
}

