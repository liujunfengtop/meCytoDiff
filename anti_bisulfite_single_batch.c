/* anti-bisulfite-single-batch.c

   Junfeng Liu, 2018
   
   This is an program to convert bisulfite reads which include assigned methylated site from the file methylation_genom-rate to normal reads.

   gcc -o anti_bisulfite_single_batch anti_bisulfite_single_batch.c -lm
   
   anti_bisulfite_single_batch <ctlfile> <methylation_genom-rate> <index>

*/

#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h> /*add head file for compile warning about isalnum function*/

struct CommonInfo {
   char *z[3], *spname[3], outf[128], outreadf[128], intransf[128], location[128], intxtf[128], ratef[128], ctlf[128], fix_locusrate, flag[128], length[128], chrom_name[128], skipped_number[128];
   int model, ncode, cleandata, seed, npoints, ncatBeta, UseMedianBeta, getSE;
   int ndata, ngene, seqtype, ns, ls, posG[1+1], lgene[1], *pose, npatt, readpattern;
   int *Nij, nGtree;
   double *fpatt, kappa, alpha, rho, rgene[1], pi[4], piG[1][4];
   double *lnLmax, *locusrate;
   double *pDclass, *tau1beta, *bp0124[5], *wwprior[5];
}  com;




int GetOptions (char *ctlf);

/*The function is used to convert the location in the chromsome to the location in the transcript*/
int GetNumber (char *transf);


FILE *fout, *foutread, *foutread1, *foutread2, *fintxt, *fsummary, *frate;

char **summ_chrom, **summ_location;

int main (int argc, char* argv[])
{   
   /*i for the skipped rows from the txt file;
     j for chracter loop about each row of the sam file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
     m for the number of mehtylation site in one read;
   */
   int i,j,read,k,l,m;
   
   /*n1 for the start position of the segment on each row of the txt file;
     n2 for the end position of the first read on each row of the txt file;
     n3 for the end position of the segment on each row of the txt file;
     n4 for the start position the second read on each row of the txt file;
     n5 is the number of the transcripts which include the assigned methylated site;
     location for the the location in the transcript according to the chromsome name which is from Bismark;
   */
   long n1,n2,n3,n4,n5,n6,location;
   
   /*ls[] for each row of the txt file;*/
   char ls[5000];
   
   /*temp[] for each column of each row of the txt file;
     temp2[] for seq-name;
     read_1[] for the first read on each row of the txt file;
     read_1_call[] for the first read methylation call on each row of the txt file;
     read_2[] for the second read on each row of the txt file;
     read_2_call[] for the second read methylation call on each row of the txt file;
   */
   char temp[2000],temp2[2000],read_1[2000],read_1_call[2000],read_2[2000],read_2_call[2000];
   
   /*outlsfile1[] for recording outread1;
     outlsfile2[] for recording outread2;
   */
   char outlsfile1[128],outlsfile2[128];
   
   /*pre_chromname for the chromosome name from the previous row*/
   char pre_chromname[128];
   
   /*tag for indication whether the position of the assigned methylation site lie in the segment, '1' means Yes, '0' means No;
     tag1 for indication whether the assigned methylation site lie in the segment, '1' means Yes, '0' means No;
     te for indication of the strand, '1' means the top strand, '0' means the bottom strand;
   */
   char tag, tag1, te, lsls;
  
   char ratef[1000], read_index[1000];

   /*read the control file and the other parameter*/
   if(argc>3) {
   	strcpy(com.ctlf, argv[1]);
   	strcpy(ratef, argv[2]);
   	strcpy(read_index, argv[3]);
   } else {
   	return -1;
   }
   
   GetOptions (com.ctlf);
   
   
   /*open the file*/
   if((fout=fopen(com.outf, "w"))==NULL) {
   		printf("Open %s file Error: %s\n",com.outf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.outf,strerror(errno));
   		return -1;
   }
   
   if((fsummary=fopen("methylation_summary_sam", "w"))==NULL) {
   		printf("Open %s file Error: %s\n","methylation_summary_sam",strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n","methylation_summary_sam",strerror(errno));
   		return -1;
   }

   if((fintxt=fopen(com.intxtf,"r"))==NULL) {
  	  printf("Open %s file Error: %s\n",com.intxtf,strerror(errno));
   		fprintf(fout,"Open %s file Error: %s\n",com.intxtf,strerror(errno));
   		return -1;
   }
   
   if((frate=fopen(ratef,"r"))==NULL) {
   	return -1;
   }
   
   /*reset the values of com.location&com.chrom_name according to the file methylation_genom-rate*/
   memset(com.location,0,sizeof(com.location));
   memset(com.chrom_name,0,sizeof(com.chrom_name));
   for(l=0;l<atol(read_index);l++) {
   	fgets(ls,5000,frate);
   	if(l==atol(read_index)-1) {
   		j=0;
   		for(k=0;k<strlen(ls);k++) {
   			if(ls[k]=='-') {
   				j=k;
   			} else {
   				if(j==0) {
   					com.chrom_name[k]=ls[k];
   				}
   				if(j>0&&k>j) {
   					if((ls[k]=='\t')||(ls[k]=='\n')) {
   						break;
   					} else {
   						com.location[k-j-1]=ls[k];
   					}
   				}
   			}
   		}
   	}
   }
   
   /*convert the position of the assigned methylation site in genome to in the transcript and return the number of the transcript which include the assigned methylation site*/
   n5=GetNumber(com.intransf);
   
   if(n5==0) {
   	 printf("No data include the assigned methylated site!\n");
   	 fprintf(fout,"No data include the assigned methylated site!\n");
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
   if(atol(com.skipped_number)>0) {
   	for(i=0;i<atol(com.skipped_number);i++) {
   		fgets(ls,5000,fintxt);
   	}
   }
   
   n6=0;
   /*read each row of intxt file */
   while(fgets(ls,5000,fintxt) != NULL) {
   	  k=0;
   	  j=0;
   	  read=0;
   	  m=0;
   	  tag1='0';
   	  memset(temp,0,sizeof(temp));
   	  memset(temp2,0,sizeof(temp2));
      memset(read_1,0,sizeof(read_1));
      memset(read_2,0,sizeof(read_2));
      memset(read_1_call,0,sizeof(read_1_call));
      memset(read_2_call,0,sizeof(read_2_call));
   	  /*output the row of the read file(or the read1 file and the read2 file) if the row includes the assigned methylated site*/
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
   	  		/*read the second column of the row and record the type of strand*/
   	  		if(read==1) {
   	  			te=temp[0];
   	  		}
   	  		/*read the third column of the row and record the chrosome name, the '-*' of the chrosome name is discarded; Moreover, getting the location*/
   	  		if(read==2) {
   	  			if(strlen(pre_chromname)==k&&strncmp(pre_chromname,temp,k)==0) {
   	  				if(location==0) {
   	  					break;
   	  				}
   	  			} else {
   	  				memset(pre_chromname,0,sizeof(pre_chromname));
   	  				for(l=0;l<k;l++) {
   	  				  pre_chromname[l]=temp[l];
   	  			  }
   	  			  for(l=0;l<n5;l++) {
   	  			  	if(strcmp(pre_chromname,summ_chrom[l])==0) {
   	  			  		location=atol(summ_location[l]);
   	  			  		break;
   	  			  	} else {
   	  			  		location=0;
   	  			  	}
   	  			  }
   	  			  if(location==0) {
   	  			  	break;
   	  			  } else {
   	  			  	n6=n6+1;
   	  			  }
   	  			}
   	  		}
   	  		/*read the fourth column of the row and record the start position and the end position of the first read*/
   	  		if(read==3) {
   	  			n1=atol(temp);
   	  			n2=n1+atol(com.length);
   	  		}
         /*read the fifth column of the row*/               
         if(read==4) {
            n3=atol(temp);
            n4=n3-atol(com.length);
            if((location>=n1&&location<n2)||(location>n4&&location<=n3)) {
              tag='1';
            } else {
              tag='0';
              break;
            }
         }
         /*read the sixth column of the row*/
         if(read==5) {
            for(l=0;l<k;l++) {
                read_1[l]=temp[l];
            }
         }
                        /*read the eigth column of the row and identify the tag1*/
                        if(read==7) {
                          if(te=='+') {
                            if((tag=='1')&&(location>=n1&&location<n2)) {
                              l=location-n1;
                              if(temp[l]=='H'||temp[l]=='X'||temp[l]=='Z'||temp[l]=='U') {
                                tag1='1';
                                m=m+1;
                              }
                              if(temp[l]=='h'||temp[l]=='x'||temp[l]=='z'||temp[l]=='u') {
                              	tag1='2';
                              }
                            }
                          }
                          if(te=='-') {
                            if((tag=='1')&&(location>n4&&location<=n3)) {
                              l=n3-location;
                              if(temp[l]=='H'||temp[l]=='X'||temp[l]=='Z'||temp[l]=='U') {
                                tag1='1';
                                m=m+1;
                              }
                              if(temp[l]=='h'||temp[l]=='x'||temp[l]=='z'||temp[l]=='u') {
                              	tag1='2';
                              }
                            }
                          }
                          for(l=0;l<k;l++) {
                          	read_1_call[l]=temp[l];
                          }
                        }
                        /*read the nineth column of the row*/
                        if(read==8) {
                          for(l=0;l<k;l++) {
                            read_2[l]=temp[l];
                          }
                        }
                        /*read the eleventh column of the row and identify the tag1*/
                        if(read==10) {
                          if(te=='-') {
                            if((tag=='1')&&(location>=n1&&location<n2)) {
                              l=location-n1;
                              if(temp[l]=='H'||temp[l]=='X'||temp[l]=='Z'||temp[l]=='U') {
                                tag1='1';
                                m=+1;
                              }
                              if(temp[l]=='h'||temp[l]=='x'||temp[l]=='z'||temp[l]=='u') {
                              	tag1='2';
                              }
                            }
                          }
                          if(te=='+') {
                            if((tag=='1')&&(location>n4&&location<=n3)) {
                              l=n3-location;
                              if(temp[l]=='H'||temp[l]=='X'||temp[l]=='Z'||temp[l]=='U') {
                                tag1='1';
                                m=m+1;
                              }
                              if(temp[l]=='h'||temp[l]=='x'||temp[l]=='z'||temp[l]=='u') {
                              	tag1='2';
                              }
                            }
                          }
                          for(l=0;l<k;l++) {
                          	read_2_call[l]=temp[l];
                          }
                          if(tag1=='1') {
                          	fprintf(fsummary,"methylated_liu");
                          	for(l=1;l<strlen(temp2);l++) {
                          		fprintf(fsummary,"%c",temp2[l]);
                          	}
                          	fprintf(fsummary,"\t%d\n",m);
                          } else {
                          	if(tag1=='2') {
                          		for(l=1;l<strlen(temp2);l++) {
                          			fprintf(fsummary,"%c",temp2[l]);
                          		}
                          		fprintf(fsummary,"\t%d\n",0);
                          	}
                          }
                        }
                        /*read the twelfth column of the row and output according the convert type*/
                        if(read==11) {
                          if((strncmp(temp,"CT",2)==0)&&(tag1!='0')) {
                            if(tag1=='1') {
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
                          if((strncmp(temp,"GA",2)==0)&&(tag1!='0')) {
                            if(tag1=='1') {
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
                        	if(tag1!='0') {
                            fprintf(foutread1,"+\n");
                            for(l=0;l<k;l++) {
                              fprintf(foutread1,"%c",temp[l]);
                              if(l==k-1) {
                                fprintf(foutread1,"\n");
                              }
                            }
                          }  
                        }
                        /*read the fifteenth column of the row and output the score*/
                        if(read==14) {
                        	if(tag1!='0') {
                            fprintf(foutread2,"+\n");
                            for(l=0;l<k;l++) {
                              fprintf(foutread2,"%c",temp[l]);
                              if(l==k-1) {
                                fprintf(foutread2,"\n");
                              }
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
   	  if(n6==n5+30&&location==0) {
   	  	break;
   	  }
    }
    

   free(summ_chrom);
   free(summ_location);

   fclose(fout);
   fclose(fintxt);
   fclose(fsummary);
   fclose(frate);
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
   int iopt,i, nopt=9, lline=4096;
   char line[4096],*pline, opt[32], *comment="*#", *seqerrstr="0EF";
   char *optstr[] = {"outfile","intxtfile","intransfile","outreadfile","location","flag","length","chrom_name","skipped_number"};
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
                  case ( 2): sscanf(pline+1, "%s", com.intransf);   break;
                  case ( 3): sscanf(pline+1, "%s", com.outreadf);  break;
                  case ( 4): sscanf(pline+1, "%s", com.location);  break;
                  case ( 5): sscanf(pline+1, "%s", com.flag);      break;
                  case ( 6): sscanf(pline+1, "%s", com.length);    break;
                  case ( 7): sscanf(pline+1, "%s", com.chrom_name); break;
                  case ( 8): sscanf(pline+1, "%s", com.skipped_number); break;
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

int GetNumber (char *transf)
{
	/* j for chracter loop about each row of the trans file;
     read for column loop about each row;
     k for character loop about each column;
     l for 'for loop';
   */
	int k, j, l, read, n;
	
	/*n1 for the start position of the exon;
	  n2 for the end position of the exon;
	  n3 for the location in the transcript;
	  n4 for the number of the transcript which include the assigned methylation site;
	*/
	long n1, n2, n3, n4;
	
	/*ls[] for each row of the trans file;
	  temp[] for each column of the row;
	  temp2[] for chromosome name and the transcript name;
	*/
	char ls[50000], temp[1000], temp2[1000];
	
	
	
	FILE *ftrans=fopen (transf, "r");
	FILE *ftrans_out=fopen(strcat(transf,"_out"), "w");
	
	n4=0;
	
	 /*read each row of the trans file*/
	 while(fgets(ls,50000,ftrans) != NULL) {
	 	k=0;
	 	j=0;
	 	read=0;
	 	n1=0;
	 	n2=0;
	 	n3=0;
	 	memset(temp,0,sizeof(temp));
   	memset(temp2,0,sizeof(temp2));
	 	while(read<1000) {
	 		temp[k]=ls[j];
	 		if(ls[j]==' '||ls[j]=='\t'||ls[j]=='\n') {
	 			if(read==0) {
	 				for(l=0;l<k;l++) {
	 					temp2[l]=temp[l];
	 				}
	 				if(strcmp(temp2,com.chrom_name)!=0) {
	 					break;
	 				}
	 			}
	 			if(read==1) {
	 				memset(temp2,0,sizeof(temp2));
	 				for(l=0;l<k;l++) {
	 					temp2[l]=temp[l];
	 				}
	 			}
	 			if((read>=2)&&(fmod(read,2)==0)) {
	 				n1=atol(temp);
	 			}
	 			if((read>=2)&&(fmod(read,2)==1)) {
	 				n2=atol(temp);
	 			}
	 			/*get the location in the transcript and output*/
	 			if((n2>=n1)&&(n1>0)) {
	 				if(atol(com.location)>=n1&&atol(com.location)<=n2) {
	 					n3=n3+atol(com.location)-n1+1;
	 					fprintf(ftrans_out,"%s\t%ld\n",temp2,n3);
	 					n4=n4+1;
	 					break;
	 				} else {
	 					if(ls[j]=='\n') {
	 						break;
	 					} else {
	 						n3=n3+n2-n1+1;
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
	fclose(ftrans);
	fclose(ftrans_out);
	
	ftrans_out=fopen(transf,"r");
	
	summ_chrom=(char**)malloc(n4*sizeof(char*));
	summ_location=(char**)malloc(n4*sizeof(char*));
	
	for(l=0;l<n4;l++) {
		summ_chrom[l]=(char*)malloc(1000*sizeof(char));
		summ_location[l]=(char*)malloc(1000*sizeof(char));
	}
	
	n=0;
	while(fgets(ls,50000,ftrans_out)!=NULL) {
		k=0;
		j=0;
		read=0;
		memset(temp,0,sizeof(temp));
		while(read<2) {
			temp[k]=ls[j];
			if(ls[j]=='\t'||ls[j]=='\n') {
				if(read==0) {
					memset(summ_chrom[n],0,sizeof(summ_chrom[n])); //by ljf 2017-05-09
					for(l=0;l<k;l++) {
						summ_chrom[n][l]=temp[l];
					}
				}
				if(read==1) {
					memset(summ_location[n],0,sizeof(summ_location[n])); //by ljf 2017-05-09
					for(l=0;l<k;l++) {
						summ_location[n][l]=temp[l];
					}
				}
				k=0;
				memset(temp,0,sizeof(temp));
				read=read+1;
			} else {
				k=k+1;
			}
			j=j+1;
		}
		n=n+1;
	}
	
	fclose(ftrans_out);
	
	return(n4);
}
