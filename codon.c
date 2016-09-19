/*
This program was written by James McInerney 1997.  
It is copyright © James McInerney.  You may not alter or re-distribute this
code or program without the express consent of the author.  This is as much
for your protection as mine.
*/


#include "codon.h"
#include "inout.c"
#include "math.c"
#include "menus.c"
#include "codonpca.c"
#include "cutot.c"
#include "aause.c"
#include "process.c"
#include "distance.c"

int main ( void )
{
 {
    printf("%s %s %s %s %s %s %s %s",
        " ****************************************************************\n",
        "           GCUA: General Codon Usage Analysis                   \n",
        "                        by                                      \n",
        "                   James McInerney,                             \n",
        "              http://bioinf.nuim.ie,                            \n",
        "       http://www.ncbi.nlm.nih.gov/pubmed/9632833               \n",
        "                    March 2015.                                 \n",
        "****************************************************************\n");
    printf(version);
    printf("\nThis version is formatted for a maximum of %d genes\neach with a maximum length of %d bases\n", MAXNUMGEN, MAXGENLEN);
    printf("\n\n\nMake sure the preferences are set before you start\n");
printf("\n\nPress [Return] to continue: ");
getchar();
printf("\n\n\n\n\n\n\n\n\n\n");
     }
 main_menu();
return 0;
}

/*Trims the genes down so that the length if the gene 'i' 
that is passes to it is a multiple of three
***/
int trimgenes ( int i) /************************************************Trimgenes**********/
{
        if(((length[i] + 3) % 3) != 0){
        length[i]--;
        }
        if(((length[i] + 3) % 3) != 0){
        length[i]--;
        }       
return length[i];
} 


void Help ( void )  /*******************************************Help********************/
{
FILE *helpf;
char c, lin[2];
int i = 1;
	               if((helpf = fopen ("codon.hlp","r")) == NULL){  /*Check to see if it can be opened for writing*/
	                printf("\n\nCannot open codon.hlp\n\n");
	               return;
                   }
                  printf("\n\n");
                 while(feof(helpf) == 0){
                 c = getc(helpf);
                  printf("%c", c);
                  if (c == '\n') i++;
                  if ((i % 20) == NULL){
                  i++;
                   getstr("\nPress [Return] to continue, X to eXit", lin);
                   
                   switch (toupper(*lin)){
                      case 'X': printf("\n");
                                fclose(helpf);
                                return;
                                 break;
                      default : printf("\n");
                                continue;
                  }
                 }
                }
               printf("\n\nPress [Return] to continue: ");
               getchar();
               printf("\n\n\n\n");
               fclose(helpf);
          return;
}

/* Lifted from Richard Thomas, who lifted it from CLUSTAL by uncle Desmo*/
void getstr(char *instr, char *outstr)
{
	fprintf(stdout, "%s: ", instr);
	gets(outstr);
}

void prog_exit (void)
{
  printf("\n\nQuitting %s \n", progname);
  exit (0);
}/* prog_exit() */
