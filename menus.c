void main_menu( void )
{
char lin1[2];
for(;;)
 {
   printf("\n\n\t%s", "There "); 
  if (!loadedfile){
   printf("is no datafile in memory\n");
   }
  if (loadedfile && nbseq == 1){
   printf("is a single sequence in memory from \'%s\'\n", InputFilnam);
   }
  if (loadedfile && nbseq > 1){
   printf("are %d sequences in memory from \'%s\'\n", nbseq, InputFilnam);
   }
  if (loadedfile && nbseq == 0){
   printf("\n\n\tThe file exists but there are no sequences in it\n");
   }
printf("%s", BANNER);
printf("%s %s %s %s %s %s %s %s %s",
     "\n\n\n\t1.          Read in a fasta file             \n",
           "\t2.          Calculate codon usage            \n",
           "\t3.          Calculate amino acid usage       \n",
           "\t4.          Calculate base composition       \n",
           "\t5.          Multivariate Analysis            \n",
           "\t6.          Distances                        \n\n\n",
           "\tP.          (P)references                    \n",
           "\tH.          (H)elp                           \n\n",
           "\tQ.          (Q)uit program                   \n\n\n");
        
        getstr("enter your choice", lin1);

		switch(toupper(*lin1)) {
			case '1': fasta_in();
						break;
			case '2': cutot();
						break;
			case '3': AAuse();
			            break;
			case '4': BaseComp();
			            break;
			case '5': pca();
			            break;
			case '6': distances();
			            break;
			case 'P': prefsmenu();
			            break;
			case 'H': Help();
			            break;
			case 'Q': prog_exit();
						break;
			default: printf("%s","\n\n\tUnrecognized Command\n\n");
						break;

    }
  }
 return;
}


int choicemenu ( int *scr, int *fil, int *cumu, int *choi, int *allal, int *ade )
{
char lin1[2];

printf("%s", BANNER);
printf("%s %s %s %s %s %s %s %s %s %s",
     "\n\n\n\t1.      Print Codon Usage table for each gene to the screen\n",
           "\t2.      Print Codon Usage table for each gene to a file\n",
           "\t3.      Print Codon Usage table for each gene to the screen and a file\n\n",
           "\t4.      print cumulative Codon Usage table to screen \n",
           "\t5.      print cumulative Codon Usage table to a file\n",
           "\t6.      print cumulative Codon Usage table to screen and a file\n",
           "\t7.      print everything everywhere\n\n",
           "\t8.      print out ADE-formatted file\n\n",
           "\tQ.      Exit from program\n\n\n",
           "\tPress [Return] to exit this menu\n\n");
        
        getstr("enter your choice", lin1);

      		switch(toupper(*lin1)) {
			case '1': *scr = TRUE;
			          *fil = FALSE;
			          *cumu = FALSE;
			          *choi = FALSE;
						break;
			case '2': *scr = FALSE;
			          *fil = TRUE;
			          *cumu = FALSE;
			          *choi = FALSE;
						break;
		    case '3': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = FALSE;
		              *choi = FALSE;
		                break;
		    case '4': *scr = TRUE;
		              *fil = FALSE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '5': *scr = FALSE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '6': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '7': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *allal = TRUE;
		              *choi = FALSE;
		                break;
		   case '8': *scr = FALSE;
		             *fil = TRUE;
		             *cumu = FALSE;
		             *choi = FALSE;
		             *ade = TRUE;
		                break;
			case 'Q': prog_exit();
						break;
			default : *choi = FALSE;
			          break;
			}
	return 0;
}

int AAchoicemenu ( int *scr, int *fil, int *cumu, int *choi, int *allal )
{
char lin1[2];

printf("%s", BANNER);
printf("%s %s %s %s %s %s %s %s %s",
     "\n\n\n\t1.      Print Amino acid Usage table for each gene to the screen\n",
           "\t2.      Print Amino acid Usage table for each gene to a file\n",
           "\t3.      Print Amino acid Usage table for each gene to the screen and a file\n\n",
           "\t4.      print cumulative Amino acid Usage table to screen \n",
           "\t5.      print cumulative Amino acid Usage table to a file\n",
           "\t6.      print cumulative Amino acid Usage table to screen and a file\n",
           "\t7.      print everything everywhere\n\n\n",
           "\tQ.      Exit from program\n\n\n",
           "\tHit [RETURN] to exit this menu\n\n");
        
        getstr("enter your choice", lin1);

      		switch(toupper(*lin1)) {
			case '1': *scr = TRUE;
			          *fil = FALSE;
			          *cumu = FALSE;
			          *choi = FALSE;
						break;
			case '2': *scr = FALSE;
			          *fil = TRUE;
			          *cumu = FALSE;
			          *choi = FALSE;
						break;
		    case '3': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = FALSE;
		              *choi = FALSE;
		                break;
		    case '4': *scr = TRUE;
		              *fil = FALSE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '5': *scr = FALSE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '6': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '7': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *allal = TRUE;
		              *choi = FALSE;
		                break;
			case 'Q': prog_exit();
						break;
			default : *choi = FALSE;
			          break;
			}
	return 0;
}

int basechoicemenu ( int *scr, int *fil, int *cumu, int *choi, int *allal )
{
char lin1[2];

printf("%s", BANNER);
printf("%s %s %s %s %s %s %s %s %s",
     "\n\n\n\t1.      Print base composition data for each gene to the screen\n",
           "\t2.      Print base composition data for each gene to a file\n",
           "\t3.      Print base composition data for each gene to the screen and a file\n\n",
           "\t4.      print cumulative base composition data to screen \n",
           "\t5.      print cumulative base composition data to a file\n",
           "\t6.      print cumulative base composition data to screen and a file\n",
           "\t7.      print everything everywhere\n\n\n",
           "\tQ.      Exit from program\n\n\n",
           "\tHit [RETURN] to exit this menu\n\n");
        
        getstr("enter your choice", lin1);

      		switch(toupper(*lin1)) {
			case '1': *scr = TRUE;
			          *fil = FALSE;
			          *cumu = FALSE;
			          *choi = FALSE;
						break;
			case '2': *scr = FALSE;
			          *fil = TRUE;
			          *cumu = FALSE;
			          *choi = FALSE;
						break;
		    case '3': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = FALSE;
		              *choi = FALSE;
		                break;
		    case '4': *scr = TRUE;
		              *fil = FALSE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '5': *scr = FALSE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '6': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *choi = FALSE;
		                break;
		    case '7': *scr = TRUE;
		              *fil = TRUE;
		              *cumu = TRUE;
		              *allal = TRUE;
		              *choi = FALSE;
		                break;
			case 'Q': prog_exit();
						break;
			default : *choi = FALSE;
			          break;
			}
	return 0;
}

void RSCU_OR_AA(int *DNA, int *m ){
char lin1[2];

printf("%s", BANNER);
printf("\n\n\n\n%s %s %s",
       " 1.            Analyse RSCU values for codons [Default]\n",
       "2.            Analyse Amino Acid frequencies\n\n\n\n",
       "Press [Return] to go to main menu\n\n\n");
       
       getstr("enter your choice", lin1);
		switch(toupper(*lin1)) {
			case '1': *DNA = true;
			          *m = 59;
			            break;
            case '2': *DNA = false;
                       *m = 20;
                        break;
            default: *DNA = true; return;
            }
}

void COR_VAR_SSCP ( char *option ){
char lin1[2];

printf("%s", BANNER);
printf("\n\n\n\n%s %s %s %s %s",
       " 1.            Correlation Analysis            \n",
       "2.            Variance/Covariance Analysis     \n",
       "3.            SSCP Analysis                    \n",
       "4.            Correspondence Analysis [Default]\n\n\n\n",
       "Press [Return] to go to main menu\n\n\n");
       
      getstr("enter your choice", lin1);
		switch(toupper(*lin1)) {
			case '1': *option = 'r';
			              break;
			case '2': *option = 'v';
			              break;
			case '3': *option = 's';
			              break;
			case '4': *option = 'c';
			              break;
			default : *option = 'c';
			              break;
			}
}

void prefsmenu ( void ){
char lin1[2];

printf("%s", BANNER);
printf("\n\n\n\n%s %s %s",
       " 1.         Genetic code \n\n",
       "2.         Running preferences \n\n\n",
       "Press [Return] to go back to main menu\n\n\n");
       
      getstr("enter your choice", lin1);
		switch(toupper(*lin1)) {
		    case '1': gencodemenu();
		                break;
		    case '2': Errorsmenu();
		                break;
		    default : return;
		    }
	}

void gencodemenu (void){
char lin1[2];
if(BEEP) printf("\7");
printf("\n\n\tWARNING WARNING WARNING\n\n\tIn this version of the program, \n\tif you change the genetic code \n\tyou must re-load the datafile\n\n");
for(;;){
printf("%s", BANNER);
printf("\n\n1.          Universal genetic code");
if(code == 1) printf("    ON\n");
else printf("   OFF\n");
printf("2.          Mycoplasma/Spiroplasma code");
if(code == 7) printf("   ON\n\n\n");
else printf("   OFF\n\n\n");
printf("Press [Return] to go back to main menu\n\n\n");
       getstr("enter your choice", lin1);
		switch(toupper(*lin1)) {
              case '1' : code = 1;
                         continue;
              case '2' : code = 7;
                         continue;
              default  : break;
               }
               break;
            }
          return;
        }
           
void Errorsmenu ( void ){
char lin1[2];
for(;;){
printf("%s", BANNER);
printf("1.          Show a lot of the calculations");
if(calcs == 3) printf("   ON\n");
else printf("   OFF\n");
printf("2.          Show some calculations");
                    if(calcs == 2) printf("           ON\n");
                              else printf("           OFF\n");
printf("3.          Show only the results");
                   if(calcs == 1) printf("            ON\n");
                             else printf("            OFF\n");
printf("4.          Toggle System beep");
                      if(BEEP) printf("               ON\n");
                          else printf("               OFF\n");
printf("5.          Toggle Distance file format");
                 if(distance == 1) printf("      PAUP\n\n\n\n");
                                   else printf("      PHYLIP\n\n\n\n");
printf("Press [Return] to go back to main menu\n\n\n");

       getstr("enter your choice", lin1);
		switch(toupper(*lin1)) {
              case '1' : calcs = 3;
                           continue;
              case '2' : calcs = 2;
                           continue;
              case '3' : calcs = 1;
                           continue;
              case '4' : if(BEEP) BEEP = false;
                         else BEEP = true;
                           continue;
              case '5' : if(distance == 1) distance = 2;
                         else distance = 1;
                         continue;
              default  : break;
              }
              break;
             }
            return;
           }
