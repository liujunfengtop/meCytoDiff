/* m5c_filter.c

   Junfeng Liu, 2018
   
   This is an program to filter the m5c level file in order to drop the isofrom whose m5c level is more than 1.

   gcc -o m5c_filter m5c_filter.c
   
   m5c_filter <m5c level file>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>


/*The function is used to filter the isoform file in order to drop the isofrom whose FPKM is less than the assign number*/
int filter (char *m5cfile);

int main (int argc, char* argv[])
{
   int l;   
   char m5cfile[1000];

   if(argc>1) {
   	strcpy(m5cfile, argv[1]);
   } else {
   	return -1;
   }
   
   l=filter(m5cfile);
        
   return (l);
}

int filter (char *m5cfile)
{
	/* j for chracter loop about each row of the isoform file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
  */
	int k, j, l, read;	
	
	
	/*ls[] for each row of the isoform file;
	  temp[] for each column of the row from the isoform file;
	*/
	char ls[5000], temp[1000];
	
	FILE *fm5c, *fout;
	
	if((fm5c=fopen(m5cfile, "r"))==NULL) {
		return -1;
  }
	
	if((fout=fopen("filter_out", "w"))==NULL) {
		return -1;
  }

	/*skipping the first row*/
  for(l=0;l<1;l++) {
   	fgets(ls,5000,fm5c);
   	fprintf(fout,"%s",ls);
  }
   
   /*filter the isoform file*/
   while(fgets(ls,5000,fm5c) != NULL) {
   	k=0;
	 	j=0;
	 	read=0;
	 	memset(temp,0,sizeof(temp));
	 	while(read<2) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==1) {
	 				if(atof(temp)<=1) {
	 					fprintf(fout,"%s",ls);
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
	 
	fclose(fm5c);
	fclose(fout);
	return 0;
}
