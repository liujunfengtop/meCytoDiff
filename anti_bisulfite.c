/* anti_bisulfite.c

   Junfeng Liu, 2018
   
   This is an program to convert the bisulfited reads to the normal reads.

   gcc -o anti_bisulfite anti_bisulfite.c
   
   anti_bisulfite <ctlfile>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
/*include ctype.h because of implicit declaration of function on linux*/
#include<ctype.h>


struct CommonInfo {
   char *z[3], *spname[3], outf[128], outreadf[128], intransf[128], location[128], intxtf[128], inmultxtf[128], ratef[128], ctlf[128], fix_locusrate, flag[128], length[128], chrom_name[128], skipped_number[128];
   int model, ncode, cleandata, seed, npoints, ncatBeta, UseMedianBeta, getSE;
   int ndata, ngene, seqtype, ns, ls, posG[1+1], lgene[1], *pose, npatt, readpattern;
   int *Nij, nGtree;
   double *fpatt, kappa, alpha, rho, rgene[1], pi[4], piG[1][4];
   double *lnLmax, *locusrate;
   double *pDclass, *tau1beta, *bp0124[5], *wwprior[5];
}  com;




int GetOptions (char *ctlf);



FILE *fout, *foutread, *foutread1, *foutread2, *fintxt;

int main (int argc, char* argv[])
{

   /*i for the skipped rows from the txt file;
     j for chracter loop about each row of the txt file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
   */
   int i,j,read,k,l;
   
   
   /*ls[] for each row of the txt file;
     ls1[] for each row of the multxt file;
   */
   char ls[2000], ls1[2000];
   
   /*temp[] for each column of each row of the txt file;
     temp2[] for seq-name;
     temp3[] for transcript name;
     read_1[] for the first read on each row of the txt file;
     read_1_call[] for the first read methylation call on each row of the txt file;
     read_2[] for the second read on each row of the txt file;
     read_2_call[] for the second read methylation call on each row of the txt file;
   */
   char temp[2000],temp2[2000],temp3[2000],read_1[2000],read_1_call[2000],read_2[2000],read_2_call[2000];
   
   /*outlsfile1[] for recording outread1;
     outlsfile2[] for recording outread2;
   */
   char outlsfile1[128],outlsfile2[128];
   
   
   /*tag for indication whether the segment includes the methylation site, '1' means Yes, '0' means No;
   */
   char tag;
  
   

   /*read the control file*/
   strcpy(com.ctlf, "anti_bisulfite.ctl");
   if(argc>1) strcpy(com.ctlf, argv[1]);
   GetOptions (com.ctlf);
   
   
   /*open the file*/
   if((fout=fopen(com.outf, "w"))==NULL) {
   		printf("Open %s file Error: %s\n",com.outf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.outf,strerror(errno));
   		return -1;
   }
   

   if((fintxt=fopen(com.intxtf,"r"))==NULL) {
  	  printf("Open %s file Error: %s\n",com.intxtf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.intxtf,strerror(errno));
   		return -1;
   }
   
   

   if(strncmp(com.flag,"p",1)==0) {
    strcpy(outlsfile1,com.outreadf);
    if((foutread1=fopen(strcat(outlsfile1,"_1.fq"),"w"))==NULL) {
   	  printf("Open %s file Error: %s\n",outlsfile1,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",outlsfile1,strerror(errno));
   		return -1;
    }
    strcpy(outlsfile2,com.outreadf);
    if((foutread2=fopen(strcat(outlsfile2,"_2.fq"),"w"))==NULL) {
   	  printf("Open %s file Error: %s\n",outlsfile2,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",outlsfile2,strerror(errno));
   		return -1;
    }
   } else {
    if((foutread=fopen(strcat(com.outreadf,".fq"),"w"))==NULL) {
   	  printf("Open %s file Error: %s\n",com.outreadf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.outreadf,strerror(errno));
   		return -1;
    }
  }
   
     
   
   /*skipping the first com.skipped_number rows of the txt file*/
   for(i=0;i<atol(com.skipped_number);i++) {
   	fgets(ls,2000,fintxt);
   }
   
   /*read each row of intxt file */
   while(fgets(ls,2000,fintxt) != NULL) {
   	  k=0;
   	  j=0;
   	  read=0;
   	  tag='0';
   	  memset(temp,0,sizeof(temp));
   	  memset(temp2,0,sizeof(temp2));
   	  memset(temp3,0,sizeof(temp3));
      memset(read_1,0,sizeof(read_1));
      memset(read_2,0,sizeof(read_2));
      memset(read_1_call,0,sizeof(read_1_call));
      memset(read_2_call,0,sizeof(read_2_call));
   	  /*output the row of the read file(or the read1 file and the read2 file)*/
   	  while(read<15) {
   	  	temp[k]=ls[j];
   	  	/*read each column of each row from the txt file*/
   	  	if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
   	  		/*read the first column of the row and record the seq-name for temp2*/
   	  		if(read==0) {
   	  			temp2[0]='@';
   	  			for(l=0;l<k-2;l++) {
   	  				temp2[l+1]=temp[l];
   	  			}
   	  		}
   	  		if(read==2) {
   	  			for(l=0;l<k;l++) {
   	  				temp3[l]=temp[l];
   	  			}
   	  		}
   	  		
         /*read the sixth column of the row*/
         if(read==5) {
            for(l=0;l<k;l++) {
                read_1[l]=temp[l];
            }
         }
                        /*read the eigth column of the row and identify the tag*/
                        if(read==7) {
                        	for(l=0;l<k;l++) {
                        		read_1_call[l]=temp[l];
                        		if(temp[l]=='H'||temp[l]=='X'||temp[l]=='Z'||temp[l]=='U') {
                        			tag='1';
                        		}
                        	}
                        }
                        /*read the nineth column of the row*/
                        if(read==8) {
                          for(l=0;l<k;l++) {
                            read_2[l]=temp[l];
                          }
                        }
                        /*read the eleventh column of the row and identify the tag; 
                          output the seqquence name and the number of methylation site in one segment to the fsummary*/
                        if(read==10) {
                          for(l=0;l<k;l++) {
                        		read_2_call[l]=temp[l];
                        		if(temp[l]=='H'||temp[l]=='X'||temp[l]=='Z'||temp[l]=='U') {
                        			tag='1';
                        		}
                        	}
                        }                        
                        /*read the twelfth column of the row and output according the convert type*/
                        if(read==11) {
                          if((strncmp(temp,"CT",2)==0)) {
                            if(tag=='1') {
                            	for(l=0;l<strlen(temp2);l++) {
                            		if(l==0) {
                            			fprintf(foutread1,"@methylated_liu");
                            			fprintf(foutread2,"@methylated_liu");
                            		} else {
                            			fprintf(foutread1,"%c",temp2[l]);
                            			fprintf(foutread2,"%c",temp2[l]);
                            		}
                            		if(l==strlen(temp2)-1) {
                            			fprintf(foutread1,"/1\n");
                            			fprintf(foutread2,"/2\n");
                            		}
                            	}
                            } else {
                            	for(l=0;l<strlen(temp2);l++) {
                            		fprintf(foutread1,"%c",temp2[l]);
                            		fprintf(foutread2,"%c",temp2[l]);
                            		if(l==strlen(temp2)-1) {
                            			fprintf(foutread1,"/1\n");
                            			fprintf(foutread2,"/2\n");
                            		}
                            	}
                            }
                            for(l=0;l<strlen(read_1);l++) {
                              if(read_1_call[l]=='h'||read_1_call[l]=='x'||read_1_call[l]=='z'||read_1_call[l]=='u') {
                                fprintf(foutread1,"%c",'C');
                                if(l==strlen(read_1)-1) {
                                  fprintf(foutread1,"\n");
                                }
                              } else {
                                fprintf(foutread1,"%c",read_1[l]);
                                if(l==strlen(read_1)-1) {
                                  fprintf(foutread1,"\n");
                                }
                              }
                              if(read_2_call[l]=='h'||read_2_call[l]=='x'||read_2_call[l]=='z'||read_2_call[l]=='u') {
                                fprintf(foutread2,"%c",'G');
                                if(l==strlen(read_1)-1) {
                                  fprintf(foutread2,"\n");
                                }
                              } else {
                                fprintf(foutread2,"%c",read_2[l]);
                                if(l==strlen(read_1)-1) {
                                  fprintf(foutread2,"\n");
                                }
                              }
                            }
                          }
                          if((strncmp(temp,"GA",2)==0)) {
                            if(tag=='1') {
                            	for(l=0;l<strlen(temp2);l++) {
                            		if(l==0) {
                            			fprintf(foutread1,"@methylated_liu");
                            			fprintf(foutread2,"@methylated_liu");
                            		} else {
                            			fprintf(foutread1,"%c",temp2[l]);
                            			fprintf(foutread2,"%c",temp2[l]);
                            		}
                            		if(l==strlen(temp2)-1) {
                            			fprintf(foutread1,"/1\n");
                            			fprintf(foutread2,"/2\n");
                            		}
                            	}
                            } else {
                            	for(l=0;l<strlen(temp2);l++) {
                            		fprintf(foutread1,"%c",temp2[l]);
                            		fprintf(foutread2,"%c",temp2[l]);
                            		if(l==strlen(temp2)-1) {
                            			fprintf(foutread1,"/1\n");
                            			fprintf(foutread2,"/2\n");
                            		}
                            	}
                            }
                            for(l=0;l<strlen(read_1);l++) {
                              if(read_1_call[l]=='h'||read_1_call[l]=='x'||read_1_call[l]=='z'||read_1_call[l]=='u') {
                                fprintf(foutread1,"%c",'G');
                                if(l==strlen(read_1)-1) {
                                  fprintf(foutread1,"\n");
                                }
                              } else {
                                fprintf(foutread1,"%c",read_1[l]);
                                if(l==strlen(read_1)-1) {
                                  fprintf(foutread1,"\n");
                                }
                              }
                              if(read_2_call[l]=='h'||read_2_call[l]=='x'||read_2_call[l]=='z'||read_2_call[l]=='u') {
                                fprintf(foutread2,"%c",'C');
                                if(l==strlen(read_2)-1) {
                                  fprintf(foutread2,"\n");
                                }
                              } else {
                                fprintf(foutread2,"%c",read_2[l]);
                                if(l==strlen(read_2)-1) {
                                  fprintf(foutread2,"\n");
                                }
                              }
                            }
                          }
                        }
                        /*read the fourteenth column of the row and output the score*/
                        if(read==13) {
                          
                            fprintf(foutread1,"+\n");
                            for(l=0;l<k;l++) {
                              fprintf(foutread1,"%c",temp[l]);
                              if(l==k-1) {
                                fprintf(foutread1,"\n");
                              }
                            }
                          
                        }
                        /*read the fifteenth column of the row and output the score*/
                        if(read==14) {
                          
                            fprintf(foutread2,"+\n");
                            for(l=0;l<k;l++) {
                              fprintf(foutread2,"%c",temp[l]);
                              if(l==k-1) {
                                fprintf(foutread2,"\n");
                              }
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
 

   fprintf(fout,"Completing! Please check the file (%s)\n", com.outreadf);
   

   fclose(fout);
   fclose(fintxt);
   if(strncmp(com.flag,"p",1)==0) {
   	 fclose(foutread1);
   	 fclose(foutread2);
   } else {
  	 fclose(foutread);
   }
   return 0;
}


int GetOptions (char *ctlf)
{
   int iopt,i, nopt=6, lline=4096;
   char line[4096],*pline, opt[32], *comment="*#", *seqerrstr="0EF";
   char *optstr[] = {"outfile","intxtfile","outreadfile","flag","skipped_number","inmultxtfile"};
   double t=1;
   FILE  *fctl=fopen (ctlf, "r");

   if (fctl) {
   	
      for (;;) {
         if(fgets(line, lline, fctl) == NULL) break;
         if(line[0]=='/' && line[1]=='/') 
            break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (strchr(comment,line[i])) break;
         if (t==0) continue;
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "="))==NULL)
            continue;
         
         for (iopt=0; iopt<nopt; iopt++) {
            if (strncmp(opt, optstr[iopt], 11)==0)  {
               switch (iopt) {
                  case ( 0): sscanf(pline+1, "%s", com.outf);      break;
                  case ( 1): sscanf(pline+1, "%s", com.intxtf);    break;
                  case ( 2): sscanf(pline+1, "%s", com.outreadf);  break;
                  case ( 3): sscanf(pline+1, "%s", com.flag);      break;
                  case ( 4): sscanf(pline+1, "%s", com.skipped_number); break;
                  case ( 5): sscanf(pline+1, "%s", com.inmultxtf); break;
               }
               break;
            }
         }

         if (iopt==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
      }
      fclose(fctl);
   }

   return (0);
}
