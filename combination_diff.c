/* combination_diff.c

   Junfeng Liu, 2018
   
   This is an program to add site id to diff_out.tsv.

   gcc -o combination_diff combination_diff.c
   
   combination_diff <diff_out.tsv> <site_info.txt> <the number of site>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>


int ComDiff (char *diff, char *sitef, char *j);

int main (int argc, char* argv[])
{
   int l;
   
   char dif[1000], site[1000], j[1000];

   if(argc>3) {
   	strcpy(dif, argv[1]);
   	strcpy(site, argv[2]);
   	strcpy(j, argv[3]);
   } else {
   	return -1;
   }

   l=ComDiff(dif,site,j);
     
   return (l);
}

int ComDiff (char *diff, char *sitef, char *j)
{
	
	int k;
	long l;
		
	/*ls[] for each row of the hits file;
	  temp[] for each column of the row;
	*/
	char ls[5000], site_id[5000];
		
	FILE *fdif, *fsite, *fout;
	
	if((fdif=fopen(diff, "r"))==NULL) {
   		return -1;
  }  
	if((fsite=fopen(sitef, "r"))==NULL) {
   		return -1;
  }     
	if((fout=fopen("diff_out_single.tsv", "w"))==NULL) {
   		return -1;
  }
  
  memset(site_id,0,sizeof(site_id));
  for(l=0;l<atol(j);l++) {
  	fgets(ls,5000,fsite);
  	if(l==atol(j)-1) {
  		for(k=0;k<strlen(ls);k++) {
  			if((ls[k]=='\t')||(ls[k]=='\n')) {
  				break;
  			} else {
  				site_id[k]=ls[k];
  			}
  		}
  	}
  }
	
	fgets(ls,5000,fdif);
	fprintf(fout,"site_id\t%s",ls);
	
	while(fgets(ls,5000,fdif) != NULL) {
		fprintf(fout,"%s\t%s",site_id,ls);
	}

	fclose(fdif);
	fclose(fsite);
	fclose(fout);
	return 0;
}

