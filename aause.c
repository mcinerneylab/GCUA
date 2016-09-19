int  AAuse (  )/****************************AAuse********************/
{
FILE *outf;
char  outfilnam[36];
int   i = 0,
      r = 0,
      screen = FALSE,
      file = FALSE,
      cumul = FALSE,
      allall = FALSE,
      choice = TRUE;
/*#ifdef MAC
_fcreator = 'CRGR';
_ftype = 'TEXT';
#endif*/
      
    if(nbseq == 0){
       printf("\n\n\tThere are no sequences in memory\n\n");
      return 0;
     }
    for(i = 0; i < nbseq; i++){
      if(choice){
                 AAchoicemenu ( &screen, &file, &cumul, &choice, &allall ); /*Ask about the kind of output*/
		        }
     if(file == true && i == 0){
      if(cumul == FALSE){
          getstr("\n\nName of file to write Amino acid Usage for each gene?", outfilnam);
          }
      else
          getstr("\n\nName of file to write cumulative Amino acid Usage?", outfilnam);
	          if((outf = fopen (outfilnam,"w")) == NULL){                  /*Check to see if it can be opened for writing*/
	            printf("\n\nCannot open the file %s\n\n", outfilnam);
	             return 0;
                 }
               }
     if(file == true && (cumul == FALSE || allall == TRUE)){                  /*Write individual amino acid usage to a file*/
        if(i == 0){
          fprintf(outf, "Gene Name\tPhe\tLeu\tSer\tTyr\tCys\tTrp\tPro\tHis\tGln\tArg\tIle\tMet\tThr\tAsn\tLys\tVal\tAla\tAsp\tGlu\tGly\tUAA\tUAG\tUGA\n");
         }
          fprintf(outf, "%s", genname[i]);
        AAfileprint(i, outf);
       }
    if(screen && (cumul == FALSE || allall == TRUE)){                       /*Write individual amino acid usage to screen*/
         printf("\n\n");
         printf("%s", genname[i]);
        AAscreenprint(i);
     }
   }
 if(cumul){
  if(file){
      fprintf(outf, "Gene Name\tPhe\tLeu\tSer\tTyr\tCys\tTrp\tPro\tHis\tGln\tArg\tIle\tMet\tThr\tAsn\tLys\tVal\tAla\tAsp\tGlu\tGly\tUAA\tUAG\tUGA\n");
      fprintf(outf, "Cumulative");
      AAfileprint (MAXNUMGEN, outf);
      }
  if(screen){      
      printf("\n\n%s", "Cumulative Amino acid Usage");
      AAscreenprint (MAXNUMGEN);
      printf("\nPress [Return] to continue: ");
      getchar();
      }
     }
choice = TRUE;
if(file){
 fclose(outf);
 }
return 0;
}

void BaseComp (void)
{
FILE *outf;
char  outfilnam[36];
int   i = 0,
      j = 0,
      r = 0,
      screen = FALSE,
      file = FALSE,
      cumul = FALSE,
      allall = FALSE,
      choice = TRUE;
/*#ifdef MAC
_fcreator = 'CRGR';
_ftype = 'TEXT';
#endif*/
    if(nbseq == 0){
       printf("\n\n\tThere are no sequences in memory\n\n");
      return;
     }
    for(i = 0; i < nbseq; i++){
      if(choice){
                 basechoicemenu ( &screen, &file, &cumul, &choice, &allall ); /*Ask about the kind of output*/
		        }
      if(file && i == 0){
          getstr("\n\nName of file to write base composition data for each gene?", outfilnam);
	          if((outf = fopen (outfilnam,"w")) == NULL){                     /*Check to see if it can be opened for writing*/
	            printf("\n\nCannot open the file %s\n\n", outfilnam);
	             return;
                 }
                }
     if(file && (cumul == FALSE || allall == TRUE)){                          /*Write individual base composition data to a file*/
       if(i == 0){
        fprintf(outf, "Gene Name \tLength    \tAA      \tGC    \tGC1     \tGC2     \tGC3  \tGC3s  \tter\n");
        }
       fprintf(outf, "%-10s", genname[i]);
       basecompfileprint(i, outf);
       }
     if(screen && (cumul == FALSE || allall == TRUE)){                       /*Write individual base composition data to screen*/
       if(i == 0){
         printf("\n           Length      AA      GC     GC1     GC2     GC3    GC3s   ter\n");
         }
        printf("%-10s", genname[i]);     
        basecompscreenprint(i);
       }
    }
 if(cumul){
  if(file){
     fprintf(outf, "Gene Name \tLength    \tAA      \tGC    \tGC1     \tGC2     \tGC3  \tGC3s  \tter\n");
      basecompfileprint(MAXNUMGEN, outf);
      }
  if(screen){      
      printf("\n\n%s", "Cumulative base composition data\n");
      printf("\nLength      AA      GC     GC1     GC2     GC3    GC3s\n");
      basecompscreenprint(MAXNUMGEN);
      printf("\n\nPress [Return] to continue: ");
      getchar();
      }
     }
choice = TRUE;
if(file){
 fclose(outf);
 }
return;
}
