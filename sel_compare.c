/* ls.c

   Junfeng Liu, 2015
   
   This is an program to complete the temporary task about operating txt file.

   gcc -o ls ls.c
   
   ls <txtfile>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>



FILE *fout, *fout1, *ftxt1f, *ftxt2f;

int main (int argc, char* argv[])
{
   char VerStr[32] = "Version 1.0, Aug 2015";
   char txt1f[128], txt2f[128];
   char ls1[2000], ls2[2000];
   char temp[100], temp1[100], temp2[100], temp3[100];
   char te;
   int  k, k1, j, j1, read, read1, l;

   if(argc>2) {
   	strcpy(txt1f, argv[1]);
   	strcpy(txt2f, argv[2]);
   } else {
   	printf("Please input two file to compare!\n");
   	return -1;
   }

   /*open the file*/
   if((fout=fopen("compare_out", "w"))==NULL) {
   		printf("Open %s file Error: %s\n","compare_out",strerror(errno));
   		return -1;
   }
   
   if((fout1=fopen("compare_out_1", "w"))==NULL) {
   		printf("Open %s file Error: %s\n","compare_out_1",strerror(errno));
   		return -1;
   }
   
   if((ftxt1f=fopen(txt1f,"r"))==NULL) {
   	  printf("Open %s file Error: %s\n",txt1f,strerror(errno));
   		return -1;
   }
   
   if((ftxt2f=fopen(txt2f,"r"))==NULL) {
   	  printf("Open %s file Error: %s\n",txt2f,strerror(errno));
   		return -1;
   }
   
   
     
   while(fgets(ls1,2000,ftxt1f) != NULL) {
   	  k=0;
   	  j=0;
   	  read=0;
   	  while(read<1) {
   	  	temp[k]=ls1[j];
   	  	/*read each column of each row*/
   	  	if(ls1[j]==' '||ls1[j]==','||ls1[j]=='\t'||ls1[j]=='\n') {
   	  		if(read==0) {
   	  			memset(temp2,0,sizeof(temp2));
   	  			for(l=0;l<k;l++) {
   	  				temp2[l]=temp[l];
   	  			}
   	  		}
   	  		te='0';
   	  		rewind(ftxt2f);
   	  		while(fgets(ls2,2000,ftxt2f) != NULL) {
   	  			k1=0;
   	  			j1=0;
   	  			read1=0;
   	  			while(read1<1) {
   	  				temp1[k1]=ls2[j1];
   	  				if(ls2[j1]==' '||ls2[j1]==','||ls2[j1]=='\t'||ls2[j1]=='\n') {
   	  					if(read1==0) {
   	  						memset(temp3,0,sizeof(temp3));
   	  						for(l=0;l<k1;l++) {
   	  							temp3[l]=temp1[l];
   	  						}
   	  						if(strcmp(temp2,temp3)==0) {
   	  							te='1';
   	  						}
   	  					}   	  					
   	  					memset(temp1,0,sizeof(temp1));
   	  					memset(temp3,0,sizeof(temp3));
   	  		      read1=read1+1;
   	  		      k1=0;
   	  				} else {
   	  					k1=k1+1;
   	  				}
   	  				j1=j1+1;
   	  			}
   	  			if(te=='1') {
   	  				fprintf(fout, "%s", ls1);
   	  				fprintf(fout1, "%s", ls2);
   	  				break;
   	  			}
   	  		}
   	  		memset(temp,0,sizeof(temp));
   	  		memset(temp2,0,sizeof(temp2));
   	  		read=read+1;
   	  		k=0;
   	  	} else {
   	  		k=k+1;
   	  	}
   	  	j=j+1;
   	  }
    }
 
   

   fclose(ftxt1f);
   fclose(ftxt2f);
   fclose(fout);
   fclose(fout1);
   return 0;
}

