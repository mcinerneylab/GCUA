/****************************************************************>cutot<**/
void cutot ( void )
{
FILE *outf;
char outfilnam[36];
int   i = 0,
      j = 0,
      k = 0,
      r = 0,
      m = 0,
      n = 0, 
      ii = 0,
      screen = FALSE,
      file = FALSE,
      cumul = FALSE,
      allall = FALSE,
      choice = TRUE,
      ade = FALSE;
/*#ifdef MAC
_fcreator = 'R*CH';
_ftype = 'TEXT';
#endif*/
if(nbseq == 0){
   printf("\n\n\tThere are no sequences in memory\n\n");
   return;
   }
for(i = 0; i < nbseq; i++){
      if(choice){
                 choicemenu ( &screen, &file, &cumul, &choice, &allall, &ade ); /*Ask about the kind of output*/
		        }
      if(file && i == 0){
                getstr("\n\nName of file to write Codon Usage for each gene?", outfilnam);
	               if((outf = fopen (outfilnam,"w")) == NULL){  /*Check to see if it can be opened for writing*/
	                printf("\n\nCannot open the file %s\n\n", outfilnam);
	               return;
                   }
                }
     if(file && (cumul == FALSE || allall == TRUE) && ade == FALSE){ /*Write individual codon usage to a file*/
       fileprint(i, outf);
      }
     if(file && ade){
     ADEfileprint(i, outf);
     }
    if(screen && (cumul == FALSE || allall == TRUE)){ /*Write individual codon usage to screen*/
           printf("%s", genname[i]);
       screenprint(i);
       if(i == (nbseq-1)){
       printf("Press [Return] to continue: ");
       getchar();
       }
     }
  }
if(cumul){
      RSCUfunc (MAXNUMGEN);
if(file){
      fprintf(outf, "\n\n%s\n", "Cumulative Codon Usage");
      fileprint (MAXNUMGEN, outf);
      }
if(screen){      
      printf("\n\n%s", "Cumulative Codon Usage");
      screenprint (MAXNUMGEN);
      printf("Press [Return] to continue: ");
      getchar();
      printf("\n\n");
      }
  }
choice = TRUE;
if(file){
 fclose(outf);
 }
return;
}
