/*********************************includes*****/

#include <stdio.h>
/*#include <unix.h>*/
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


/********************************defines*****/

#define MAXGENLEN 2000000     /*Max length of a gene*/
#define MAXNAMLEN 200        /*Max length of gene name*/
#define MAXNUMGEN 5000000      /*Max number of genes in a single analysis*/
#define BANNER "\n***************************************************************\n*************GCUA: General Codon Usage Analysis****************\n***************************************************************\n\n"
#define progname "GCUA: General Codon Usage Analysis\n"
#define version "\nVersion: 1.1.1b1   31.03.2015\n"
#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef true
#define true 1
#endif
#ifndef false
#define false 0
#endif
#ifdef MAC  
#undef MAC      /*Edit this for unix machines*/
#endif
/********************************prototypes*****/

void getstr              ( char *instr, char *outstr );
void prog_exit           ( void   );
void main_menu           ( void   );
void fasta_in            ( void   );
void cutot               ( void   );
int  sumcodons           ( int gennum );
int  RSCUfunc            ( int gennum );
int  choicemenu          ( int *screen, int *file, int *cumul, int *choice, int *allall, int *ade );
int  AAchoicemenu        ( int *screen, int *file, int *cumul, int *choice, int *allall );
int  AAuse               (        );
int  screenprint         ( int gennum );
int  fileprint           ( int gennum, FILE *outf );
int  AAfileprint         ( int i, FILE *outf);
int  AAscreenprint       ( int i );
void Help                ( void  );
int  trimgenes           ( int i );
int  FillCodonsArray     ( int i, int j);
void BaseComp            ( void );
int  basechoicemenu      ( int *screen, int *file, int *cumul, int *choice, int *allall );
void gencodemenu         ( void );
void prefsmenu           ( void );
void Errorsmenu          ( void );
int  basecompscreenprint ( int i );
int  basecompfileprint   ( int i, FILE *outfile );
int  ADEfileprint        ( int gennum, FILE *outf );
void pca                 ( void );
int  ca                  ( float **data, int n, int m, float **A);
void erhand              ( char *err_msg  );
float **matrix           ( int n, int m );
void free_vector         ( float *v, int n );
void free_matrix         ( float **mat, int n, int m );
void tred2               ( float **a, int n, float *d, float *e );
void tqli                ( float d[], float e[], int n, float **z );
int tql2                 ( float *d, float *e,int n, float **z );
void scpcol              ( float **data, int n, int m, float **symmat );
void covcol              ( float **data, int n, int m, float **symmat );
void corcol              ( float **data, int n, int m, float **symmat );
void RSCU_OR_AA          ( int *DNA, int *m );
void COR_VAR_SSCP        ( char *option);
int  process             ( int i );
int distances            ( void );
int enc                  (int gennum );
/********************************Global variables******/

static char genname[MAXNUMGEN][MAXNAMLEN + 1], 
            s[MAXGENLEN + 1],
            ter[MAXNUMGEN][3];

int         loadedfile = 0, 
            nbseq = 0,
            code = 1, /*1 = Universal, 7 = mycoplasma/Spiroplasma*/
            doneRSCU = false,
            donesumcodons = false,
            codonstotfill = false,
            BEEP = false,
            distance = 1, /*1 = PAUP, 2 = PHYLIP*/
            calcs = 2, /*how much calculations to show 1 = results only, 2 = some calculations, 3 = everything*/
            length[MAXNUMGEN +1], 
            codonstot[MAXNUMGEN + 1][65],
            GC[MAXNUMGEN + 1],
            GC3[MAXNUMGEN + 1],
            GC1[MAXNUMGEN + 1],
            GC2[MAXNUMGEN + 1],
            nonsyn[MAXNUMGEN + 1],
            stop[MAXNUMGEN + 1];
            
float       AminoArray[23][MAXNUMGEN + 1], /*Numbers of aminoacids per gene*/
            RSCU[MAXNUMGEN + 1][65], /*Relative synon. cod. use array*/
            Cmean[59],        /*Column means for correspondence analysis*/
            Rmean[MAXNUMGEN], /*Row means for correspondence analysis*/
            effco[MAXNUMGEN]; /*Effective number of codons*/

char        InputFilnam[36];
/*Order of the amino acids in the AminoArray

             Phe = 0, 
             Leu = 1, 
             Ser = 2, 
             Tyr = 3, 
          terUAA = 4, 
          terUAG = 5, 
             Cys = 6, 
          terUGA = 7,
             Trp = 8, 
             Pro = 9, 
             His = 10, 
             Gln = 11, 
             Arg = 12, 
             Ile = 13,       
             Met = 14, 
             Thr = 15, 
             Asn = 16,  
             Lys = 17, 
             Val = 18, 
             Ala = 19, 
             Asp = 20, 
             Glu = 21, 
             Gly = 22,*/



