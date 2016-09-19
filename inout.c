/****This function reads in fasta-formatted files 
and puts the names into an array called genname[]
and the sequences into an array called s[]*****/

void fasta_in ( void )   
{
FILE *fp;
char filnam[36];
char c;
int i = 0, 
    j = 0, 
    k = 0, 
    l = 0, 
    m = 0, 
    n = 0,
    i1 = 0,
    i2 = 0;

          getstr("\n\n\tName of the Fasta-formatted file", filnam);
	if((fp = fopen (filnam,"r")) == NULL){            /*Check to see if it is there*/
	  printf("\n\n\tCannot open the file %s\n\n", filnam);
	  return;
	  }
	  strcpy(InputFilnam, filnam );
	  if(calcs != 1) printf("Opened %s...\n", InputFilnam);
for(i = 0; i <= MAXNUMGEN; i++){                      /*initialise a number of arrays*/
GC[i] = 0; GC1[i] = 0; GC2[i] = 0; GC3[i] = 0; nonsyn[i] = 0; stop[i] = 0;
 for(j = 0; j < 23; j++){
  AminoArray[j][i] = 0;
  }
}
for(i = 0; i <= MAXNUMGEN; i++){                    /*Initialise RSCU and codonstot Arrays*/
 for(j = 0; j< 64; j++){
   RSCU[i][j] = 0;
   codonstot[i][j] = 0;
   }
  }
i = 0;
 do {
if((c = getc(fp)) == '>'){
     j = 0;
    do {                                             /*read in the sequence name**/
     if(c != ' '){
         genname[i][j] = c;
         }
        if(j < MAXNAMLEN){
         j++;
         }
        }while((c = getc(fp)) != '\n' && c != '\r' && feof(fp) == 0);
         genname[i][j] = '\0';
        }
 else{ /*If the first character in the file is not a '>' then we assume the format is wrong*/
     printf("\n\n\nThis does not seem to be a Fasta-formatted file\n\n");
     return;
     }
      j = 0;
    do {                                            /*now read in the sequence data*/
         if((s[j] = c) != '\n'){
         j++;
         }
        } while((c = getc(fp)) != '>' && feof(fp) == 0); /*the presence of a '>' or an EOF tells us the sequence is finished**/
         s[j] = '\0';                               /*Append a NULL terminator*/
        if(c == '>') ungetc(c, fp);                 /*put the '>' back to the stream*/
         length[l] = j;
         length[l] = trimgenes( l );                /*trim off any unwanted nucleotides**/
               for (i1 = (length[l] - 3), i2 = 0;  i2 < 3; i1++, i2++){
                ter[nbseq][i2] = s[i1];
                }                                   /*Fill the ter[] array with the termination codon */
     if(j > MAXGENLEN){
       printf("\n\n\nSequence \'%s\' is longer than the maximum allowable\n", genname[i]);
       printf("Increase the MAXGENLEN variable in codon.h and recompile the program\n\n");
       printf("\tWARNING: Terminating data input!\n\n");
       return;
       }
       process (l);
       sumcodons(l);
       RSCUfunc(l);
       /*effco[i] = enc(l);*/
       if(calcs == 3){ /*Print data to screen*/
       printf("%s\tlength:\t%d\n", genname[i], length[i]);
        for(j = 0; j < length[l]; j++){
        printf("%c", s[j]);
        if((j+1) % 60 == 0) printf("\n");
        }
        printf("\n");
       }
     if(nbseq < (MAXNUMGEN - 1)){                   /*If we haven't reached MAXNUMGEN -1 let's go again**/
           l++;
    nbseq = ++i;
     }
     else {                                         /*otherwise inform the user and return to main_menu*/
     printf("\n\n\nThere are too many sequences in the file\n\n");
     return;
     }
     if(nbseq % 10 == 0){                           /*display some signs of progress*/
      if(calcs != 1) printf("Loaded %d sequences\n", nbseq);
     }
     } while(feof(fp) == 0);                        /*As long as we are not at the EOF, go again**/
    fclose(fp);                                     /*Close the file**/
    loadedfile=TRUE;
if(calcs == 2){
    if(nbseq < 50){                                /*If there are less than 50 sequences, print their names to the screen**/
    for(i = 0; i < nbseq; i++){
      printf("%3d.\tlen: %-7d\t", i+1, length[i]);
      printf("%s\n", genname[i]);
      }
      printf("\n");
     }
    }
  if(BEEP) printf("\7");
    if(calcs != 1) printf("Successfully read %s!\n", InputFilnam);
return;
}

int fileprint( int gennum, FILE *outf )
{

      fprintf(outf,"\nAA Codon       N   RSCU  AA Codon       N   RSCU\n\n");  
      fprintf(outf, "Phe UUU %8d (%4.2f) Ser UCU %8d (%4.2f)\n", 
        codonstot[gennum][0], RSCU[gennum][0],  codonstot[gennum][4], RSCU[gennum][4]); 
      fprintf(outf, "    UUC %8d (%4.2f)     UCC %8d (%4.2f)\n", 
        codonstot[gennum][1], RSCU[gennum][1],  codonstot[gennum][5], RSCU[gennum][5]);
      fprintf(outf, "Leu UUA %8d (%4.2f)     UCA %8d (%4.2f) \n", 
        codonstot[gennum][2], RSCU[gennum][2],  codonstot[gennum][6], RSCU[gennum][6]);
      fprintf(outf, "    UUG %8d (%4.2f)     UCG %8d (%4.2f)\n", 
        codonstot[gennum][3], RSCU[gennum][3],  codonstot[gennum][7], RSCU[gennum][7]);
      fprintf(outf, "Tyr UAU %8d (%4.2f) Cys UGU %8d (%4.2f)\n",
        codonstot[gennum][8], RSCU[gennum][8], codonstot[gennum][12], RSCU[gennum][12] );
      fprintf(outf, "    UAC %8d (%4.2f)     UGC %8d (%4.2f)\n",
        codonstot[gennum][9], RSCU[gennum][9], codonstot[gennum][13], RSCU[gennum][13] );
if(code == 1){
     fprintf(outf,  "ter UAA %8d (%4.2f) ter UGA %8d (%4.2f)\n",
       codonstot[gennum][10], RSCU[gennum][10], codonstot[gennum][14], RSCU[gennum][14] );
       }
if(code == 7){
     fprintf(outf,  "ter UAA %8d (%4.2f) Trp UGA %8d (%4.2f)\n",
       codonstot[gennum][10], RSCU[gennum][10], codonstot[gennum][14], RSCU[gennum][14] );
      }
     fprintf(outf,  "ter UAG %8d (%4.2f) Trp UGG %8d (%4.2f)\n\n",
       codonstot[gennum][11], RSCU[gennum][11], codonstot[gennum][15], RSCU[gennum][15] );
    
      fprintf(outf, "Leu CUU %8d (%4.2f) Pro CCU %8d (%4.2f)\n", 
        codonstot[gennum][16], RSCU[gennum][16], codonstot[gennum][20], RSCU[gennum][20]); 
      fprintf(outf, "    CUC %8d (%4.2f)     CCC %8d (%4.2f)\n", 
        codonstot[gennum][17], RSCU[gennum][17], codonstot[gennum][21], RSCU[gennum][21]); 
      fprintf(outf, "    CUA %8d (%4.2f)     CCA %8d (%4.2f)\n", 
        codonstot[gennum][18], RSCU[gennum][18], codonstot[gennum][22], RSCU[gennum][22]);
      fprintf(outf, "    CUG %8d (%4.2f)     CCG %8d (%4.2f)\n", 
        codonstot[gennum][19], RSCU[gennum][19], codonstot[gennum][23], RSCU[gennum][23]);
      fprintf(outf, "His CAU %8d (%4.2f) Arg CGU %8d (%4.2f)\n",
        codonstot[gennum][24], RSCU[gennum][24], codonstot[gennum][28], RSCU[gennum][28] );
      fprintf(outf, "    CAC %8d (%4.2f)     CGC %8d (%4.2f)\n",
        codonstot[gennum][25], RSCU[gennum][25], codonstot[gennum][29], RSCU[gennum][29] );
      fprintf(outf, "Gln CAA %8d (%4.2f)     CGA %8d (%4.2f)\n",
        codonstot[gennum][26], RSCU[gennum][26], codonstot[gennum][30], RSCU[gennum][30] );
      fprintf(outf, "    CAG %8d (%4.2f)     CGG %8d (%4.2f)\n\n",
        codonstot[gennum][27], RSCU[gennum][27], codonstot[gennum][31], RSCU[gennum][31] );
      
      
      fprintf(outf, "Ile AUU %8d (%4.2f) Thr ACU %8d (%4.2f)\n", 
        codonstot[gennum][32], RSCU[gennum][32], codonstot[gennum][36], RSCU[gennum][36]); 
      fprintf(outf, "    AUC %8d (%4.2f)     ACC %8d (%4.2f)\n", 
        codonstot[gennum][33], RSCU[gennum][33], codonstot[gennum][37], RSCU[gennum][37]); 
      fprintf(outf, "    AUA %8d (%4.2f)     ACA %8d (%4.2f)\n", 
        codonstot[gennum][34], RSCU[gennum][34], codonstot[gennum][38], RSCU[gennum][38]); 
      fprintf(outf, "Met AUG %8d (%4.2f)     ACG %8d (%4.2f)\n", 
        codonstot[gennum][35], RSCU[gennum][35], codonstot[gennum][39], RSCU[gennum][39]); 
      fprintf(outf, "Asn AAU %8d (%4.2f) Ser AGU %8d (%4.2f)\n",
        codonstot[gennum][40], RSCU[gennum][40], codonstot[gennum][44], RSCU[gennum][44] );
      fprintf(outf, "    AAC %8d (%4.2f)     AGC %8d (%4.2f)\n",
         codonstot[gennum][41], RSCU[gennum][41], codonstot[gennum][45], RSCU[gennum][45] );
      fprintf(outf, "Lys AAA %8d (%4.2f) Arg AGA %8d (%4.2f)\n",
        codonstot[gennum][42], RSCU[gennum][42], codonstot[gennum][46], RSCU[gennum][46] );
      fprintf(outf, "    AAG %8d (%4.2f)     AGG %8d (%4.2f)\n\n",
        codonstot[gennum][43], RSCU[gennum][43], codonstot[gennum][47], RSCU[gennum][47] );
      
      fprintf(outf, "Val GUU %8d (%4.2f) Ala GCU %8d (%4.2f)\n", 
        codonstot[gennum][48], RSCU[gennum][48], codonstot[gennum][52], RSCU[gennum][52]); 
      fprintf(outf, "    GUC %8d (%4.2f)     GCC %8d (%4.2f)\n", 
        codonstot[gennum][49], RSCU[gennum][49], codonstot[gennum][53], RSCU[gennum][53]); 
      fprintf(outf, "    GUA %8d (%4.2f)     GCA %8d (%4.2f)\n", 
        codonstot[gennum][50], RSCU[gennum][50], codonstot[gennum][54], RSCU[gennum][54]);
      fprintf(outf, "    GUG %8d (%4.2f)     GCG %8d (%4.2f)\n", 
        codonstot[gennum][51], RSCU[gennum][51], codonstot[gennum][55], RSCU[gennum][55]); 
      fprintf(outf, "Asp GAU %8d (%4.2f) Gly GGU %8d (%4.2f)\n",
        codonstot[gennum][56], RSCU[gennum][56], codonstot[gennum][60], RSCU[gennum][60] );
      fprintf(outf, "    GAC %8d (%4.2f)     GGC %8d (%4.2f)\n",
        codonstot[gennum][57], RSCU[gennum][57], codonstot[gennum][61], RSCU[gennum][61] );
      fprintf(outf, "Glu GAA %8d (%4.2f)     GGA %8d (%4.2f)\n",
        codonstot[gennum][58], RSCU[gennum][58], codonstot[gennum][62], RSCU[gennum][62] );
      fprintf(outf, "    GAG %8d (%4.2f)     GGG %8d (%4.2f)\n\n\n",
        codonstot[gennum][59], RSCU[gennum][59], codonstot[gennum][63], RSCU[gennum][63] );
return 0;
}
/*******Screenprint*******/
int screenprint ( int gennum )
{
      printf( "\nAA Codon       N   RSCU  AA Codon       N   RSCU\n\n");  
      printf( "Phe UUU %8d (%4.2f) Ser UCU %8d (%4.2f)\n", 
        codonstot[gennum][0], RSCU[gennum][0],  codonstot[gennum][4], RSCU[gennum][4]); 
      printf( "    UUC %8d (%4.2f)     UCC %8d (%4.2f)\n", 
        codonstot[gennum][1], RSCU[gennum][1],  codonstot[gennum][5], RSCU[gennum][5]);
      printf( "Leu UUA %8d (%4.2f)     UCA %8d (%4.2f) \n", 
        codonstot[gennum][2], RSCU[gennum][2],  codonstot[gennum][6], RSCU[gennum][6]);
      printf( "    UUG %8d (%4.2f)     UCG %8d (%4.2f)\n", 
        codonstot[gennum][3], RSCU[gennum][3],  codonstot[gennum][7], RSCU[gennum][7]);
      printf( "Tyr UAU %8d (%4.2f) Cys UGU %8d (%4.2f)\n",
        codonstot[gennum][8], RSCU[gennum][8], codonstot[gennum][12], RSCU[gennum][12] );
      printf( "    UAC %8d (%4.2f)     UGC %8d (%4.2f)\n",
        codonstot[gennum][9], RSCU[gennum][9], codonstot[gennum][13], RSCU[gennum][13] );
if(code == 1){
     printf(  "ter UAA %8d (%4.2f) ter UGA %8d (%4.2f)\n",
       codonstot[gennum][10], RSCU[gennum][10], codonstot[gennum][14], RSCU[gennum][14] );
       }
else{ printf(  "ter UAA %8d (%4.2f) Trp UGA %8d (%4.2f)\n",
       codonstot[gennum][10], RSCU[gennum][10], codonstot[gennum][14], RSCU[gennum][14] );
}
     printf( "ter UAG %8d (%4.2f) Trp UGG %8d (%4.2f)\n\n",
       codonstot[gennum][11], RSCU[gennum][11], codonstot[gennum][15], RSCU[gennum][15] );
    
      printf( "Leu CUU %8d (%4.2f) Pro CCU %8d (%4.2f)\n", 
        codonstot[gennum][16], RSCU[gennum][16], codonstot[gennum][20], RSCU[gennum][20]); 
      printf( "    CUC %8d (%4.2f)     CCC %8d (%4.2f)\n", 
        codonstot[gennum][17], RSCU[gennum][17], codonstot[gennum][21], RSCU[gennum][21]); 
      printf( "    CUA %8d (%4.2f)     CCA %8d (%4.2f)\n", 
        codonstot[gennum][18], RSCU[gennum][18], codonstot[gennum][22], RSCU[gennum][22]);
      printf( "    CUG %8d (%4.2f)     CCG %8d (%4.2f)\n", 
        codonstot[gennum][19], RSCU[gennum][19], codonstot[gennum][23], RSCU[gennum][23]);
      printf( "His CAU %8d (%4.2f) Arg CGU %8d (%4.2f)\n",
        codonstot[gennum][24], RSCU[gennum][24], codonstot[gennum][28], RSCU[gennum][28] );
      printf( "    CAC %8d (%4.2f)     CGC %8d (%4.2f)\n",
        codonstot[gennum][25], RSCU[gennum][25], codonstot[gennum][29], RSCU[gennum][29] );
      printf( "Gln CAA %8d (%4.2f)     CGA %8d (%4.2f)\n",
        codonstot[gennum][26], RSCU[gennum][26], codonstot[gennum][30], RSCU[gennum][30] );
      printf( "    CAG %8d (%4.2f)     CGG %8d (%4.2f)\n\n",
        codonstot[gennum][27], RSCU[gennum][27], codonstot[gennum][31], RSCU[gennum][31] );
      
      
      printf( "Ile AUU %8d (%4.2f) Thr ACU %8d (%4.2f)\n", 
        codonstot[gennum][32], RSCU[gennum][32], codonstot[gennum][36], RSCU[gennum][36]); 
      printf( "    AUC %8d (%4.2f)     ACC %8d (%4.2f)\n", 
        codonstot[gennum][33], RSCU[gennum][33], codonstot[gennum][37], RSCU[gennum][37]); 
      printf( "    AUA %8d (%4.2f)     ACA %8d (%4.2f)\n", 
        codonstot[gennum][34], RSCU[gennum][34], codonstot[gennum][38], RSCU[gennum][38]); 
      printf( "Met AUG %8d (%4.2f)     ACG %8d (%4.2f)\n", 
        codonstot[gennum][35], RSCU[gennum][35], codonstot[gennum][39], RSCU[gennum][39]); 
      printf( "Asn AAU %8d (%4.2f) Ser AGU %8d (%4.2f)\n",
        codonstot[gennum][40], RSCU[gennum][40], codonstot[gennum][44], RSCU[gennum][44] );
      printf( "    AAC %8d (%4.2f)     AGC %8d (%4.2f)\n",
         codonstot[gennum][41], RSCU[gennum][41], codonstot[gennum][45], RSCU[gennum][45] );
      printf( "Lys AAA %8d (%4.2f) Arg AGA %8d (%4.2f)\n",
        codonstot[gennum][42], RSCU[gennum][42], codonstot[gennum][46], RSCU[gennum][46] );
      printf( "    AAG %8d (%4.2f)     AGG %8d (%4.2f)\n\n",
        codonstot[gennum][43], RSCU[gennum][43], codonstot[gennum][47], RSCU[gennum][47] );
      
      
      printf( "Val GUU %8d (%4.2f) Ala GCU %8d (%4.2f)\n", 
        codonstot[gennum][48], RSCU[gennum][48], codonstot[gennum][52], RSCU[gennum][52]); 
      printf( "    GUC %8d (%4.2f)     GCC %8d (%4.2f)\n", 
        codonstot[gennum][49], RSCU[gennum][49], codonstot[gennum][53], RSCU[gennum][53]); 
      printf( "    GUA %8d (%4.2f)     GCA %8d (%4.2f)\n", 
        codonstot[gennum][50], RSCU[gennum][50], codonstot[gennum][54], RSCU[gennum][54]);
      printf( "    GUG %8d (%4.2f)     GCG %8d (%4.2f)\n", 
        codonstot[gennum][51], RSCU[gennum][51], codonstot[gennum][55], RSCU[gennum][55]); 
      printf( "Asp GAU %8d (%4.2f) Gly GGU %8d (%4.2f)\n",
        codonstot[gennum][56], RSCU[gennum][56], codonstot[gennum][60], RSCU[gennum][60] );
      printf( "    GAC %8d (%4.2f)     GGC %8d (%4.2f)\n",
        codonstot[gennum][57], RSCU[gennum][57], codonstot[gennum][61], RSCU[gennum][61] );
      printf( "Glu GAA %8d (%4.2f)     GGA %8d (%4.2f)\n",
        codonstot[gennum][58], RSCU[gennum][58], codonstot[gennum][62], RSCU[gennum][62] );
      printf( "    GAG %8d (%4.2f)     GGG %8d (%4.2f)\n\n\n",
        codonstot[gennum][59], RSCU[gennum][59], codonstot[gennum][63], RSCU[gennum][63] );
return 0;
}

int AAfileprint (int i, FILE *outf)
{
float totAA = 0;
int  j;
 for(j = 0; j < 23; j ++){
 if(AminoArray[j][i] == -1){
  AminoArray[j][i] = 0;
  }
  totAA += AminoArray[j][i];
  }
     fprintf(outf, "\t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f \t%5.2f",
          (int)AminoArray[0][i]*100 / totAA,(int)AminoArray[1][i]*100 / totAA,(int)AminoArray[2][i]*100 / totAA,(int)AminoArray[3][i]*100 / totAA,
          (int)AminoArray[6][i]*100 / totAA,(int)AminoArray[8][i]*100 / totAA,(int)AminoArray[9][i]*100 / totAA,(int)AminoArray[10][i]*100 / totAA,
          (int)AminoArray[11][i]*100 / totAA,(int)AminoArray[12][i]*100 / totAA,(int)AminoArray[13][i]*100 / totAA,(int)AminoArray[14][i]*100 / totAA,
          (int)AminoArray[15][i]*100 / totAA,(int)AminoArray[16][i]*100 / totAA,(int)AminoArray[17][i]*100 / totAA,(int)AminoArray[18][i]*100 / totAA,
          (int)AminoArray[19][i]*100 / totAA,(int)AminoArray[20][i]*100 / totAA,(int)AminoArray[21][i]*100 / totAA,(int)AminoArray[22][i]*100 / totAA,
          (int)AminoArray[4][i]*100 / totAA,(int)AminoArray[5][i]*100 / totAA);
          if(code == 1) fprintf(outf, " \t%5.2f\n", (int)AminoArray[7][i]*100 / totAA);
          else fprintf(outf, "\n");
     return 0;
}

int AAscreenprint ( int i )
{
float totAA = 0;
int  j;
 for(j = 0; j < 23; j ++){
 if(AminoArray[j][i] == -1){
  AminoArray[j][i] = 0;
  }
  totAA += AminoArray[j][i];
  }
     printf("  %7.0f codons\n", totAA);
     printf(" AA         N   Percent\n");
     printf( " Phe: %7d\t%5.2f\n", (int)AminoArray[0][i], ((AminoArray[0][i] *100) / totAA));
     printf( " Leu: %7d\t%5.2f\n", (int)AminoArray[1][i], ((AminoArray[1][i] *100) / totAA));
     printf( " Ser: %7d\t%5.2f\n", (int)AminoArray[2][i], ((AminoArray[2][i] *100) / totAA));
     printf( " Tyr: %7d\t%5.2f\n", (int)AminoArray[3][i], ((AminoArray[3][i] *100) / totAA));
     printf( " Cys: %7d\t%5.2f\n", (int)AminoArray[6][i], ((AminoArray[6][i] *100) / totAA));
     printf( " Trp: %7d\t%5.2f\n", (int)AminoArray[8][i], ((AminoArray[8][i] *100) / totAA));
     printf( " Pro: %7d\t%5.2f\n", (int)AminoArray[9][i], ((AminoArray[9][i] *100) / totAA));
     printf( " His: %7d\t%5.2f\n", (int)AminoArray[10][i], ((AminoArray[10][i] *100) / totAA));
     printf( " Gln: %7d\t%5.2f\n", (int)AminoArray[11][i], ((AminoArray[11][i] *100) / totAA));
     printf( " Arg: %7d\t%5.2f\n", (int)AminoArray[12][i], ((AminoArray[12][i] *100) / totAA));
     printf( " Ile: %7d\t%5.2f\n", (int)AminoArray[13][i], ((AminoArray[13][i] *100) / totAA));
     printf( " Met: %7d\t%5.2f\n", (int)AminoArray[14][i], ((AminoArray[14][i] *100) / totAA));
     printf( " Thr: %7d\t%5.2f\n", (int)AminoArray[15][i], ((AminoArray[15][i] *100) / totAA));
     printf( " Asn: %7d\t%5.2f\n", (int)AminoArray[16][i], ((AminoArray[16][i] *100) / totAA));
     printf( " Lys: %7d\t%5.2f\n", (int)AminoArray[17][i], ((AminoArray[17][i] *100) / totAA));
     printf( " Val: %7d\t%5.2f\n", (int)AminoArray[18][i], ((AminoArray[18][i] *100) / totAA));
     printf( " Ala: %7d\t%5.2f\n", (int)AminoArray[19][i], ((AminoArray[19][i] *100) / totAA));
     printf( " Asp: %7d\t%5.2f\n", (int)AminoArray[20][i], ((AminoArray[20][i] *100) / totAA));
     printf( " Glu: %7d\t%5.2f\n", (int)AminoArray[21][i], ((AminoArray[21][i] *100) / totAA));
     printf( " Gly: %7d\t%5.2f\n", (int)AminoArray[22][i], ((AminoArray[22][i] *100) / totAA));
     printf( " UAA: %7d\t%5.2f\n", (int)AminoArray[4][i], ((AminoArray[4][i] *100) / totAA));
     printf( " UAG: %7d\t%5.2f\n", (int)AminoArray[5][i], ((AminoArray[5][i] *100) / totAA));
if(code == 1){
     printf( " UGA: %7d\t%5.2f\n", (int)AminoArray[7][i], ((AminoArray[7][i] *100) / totAA));
 }
     return 0;
}


int basecompscreenprint ( int i )
{
int j, k;

if(i == MAXNUMGEN){
length[MAXNUMGEN] = 0;
stop[MAXNUMGEN] = 0;
nonsyn[MAXNUMGEN] = 0;
 for(k = 0; k < nbseq; k++){
  length[MAXNUMGEN] += length[k];
  stop[MAXNUMGEN] += stop[k];
  nonsyn[MAXNUMGEN] += nonsyn[k];
  }
 }
     printf("%7d %7d %7.2f %7.2f %7.2f %7.2f %7.2f\t", 
        length[i] - (stop[i] * 3), 
        ((length[i] ) / 3) - stop[i], 
        (GC[i]   /  ((float)length[i] - (3 * stop[i]              )     )) * 100, 
        (GC1[i]  / (((float)length[i] - (3 * stop[i]              )) / 3)) * 100, 
        (GC2[i]  / (((float)length[i] - (3 * stop[i]              )) / 3)) * 100, 
        (GC3[i]  / (((float)length[i] - (3 * stop[i]              )) / 3)) * 100, 
        ((GC3[i] - (nonsyn[i] + stop[i])) / (((float)length[i] - (3 * (nonsyn[i] + stop[i]))) / 3)) * 100);
      if(i != MAXNUMGEN){
      if(stop[i] != 0){
       for(j = 0; j < 3; j++){
        printf("%c", ter[i][j]);
        }
        printf("\n");
       }
     else printf("\n");
     }
     return 0;
     }
int basecompfileprint ( int i, FILE *outf )
{
int j,k;
if(i == MAXNUMGEN){
length[MAXNUMGEN] = 0;
stop[MAXNUMGEN] = 0;
nonsyn[MAXNUMGEN] = 0;
 for(k = 0; k < nbseq; k++){
  length[MAXNUMGEN] += length[k];
  stop[MAXNUMGEN] += stop[k];
  nonsyn[MAXNUMGEN] += nonsyn[k];
  }
 }
     fprintf(outf,"\t%7d \t%7d \t%7.2f \t%7.2f \t%7.2f \t%7.2f \t%7.2f\t", 
        length[i] - (stop[i] * 3), 
        ((length[i] ) / 3) - stop[i], 
        (GC[i]   /  ((float)length[i] - (3 * stop[i]              )     )) * 100, 
        (GC1[i]  / (((float)length[i] - (3 * stop[i]              )) / 3)) * 100, 
        (GC2[i]  / (((float)length[i] - (3 * stop[i]              )) / 3)) * 100, 
        (GC3[i]  / (((float)length[i] - (3 * stop[i]              )) / 3)) * 100, 
        ((GC3[i] - (nonsyn[i] + stop[i])) / (((float)length[i] - (3 * (nonsyn[i] + stop[i]))) / 3)) * 100);
    if(i != MAXNUMGEN){
    if(stop[i] != 0){
       for(j = 0; j < 3; j++){
        fprintf(outf,"%c", ter[i][j]);
        }
        fprintf(outf, "\n");
       }
     else fprintf(outf, "\n");
     }
     return 0;
     }
     
int  ADEfileprint ( int gennum, FILE *outf )
{
       fprintf(outf, "%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n",
                    RSCU[gennum][0],RSCU[gennum][1],RSCU[gennum][2],RSCU[gennum][3],
                    RSCU[gennum][4],RSCU[gennum][5],RSCU[gennum][6],RSCU[gennum][7],
                    RSCU[gennum][8],RSCU[gennum][9],RSCU[gennum][12],RSCU[gennum][13],
                    RSCU[gennum][16],RSCU[gennum][17],RSCU[gennum][18],RSCU[gennum][19],
                    RSCU[gennum][20],RSCU[gennum][21],RSCU[gennum][22],RSCU[gennum][23],
                    RSCU[gennum][24],RSCU[gennum][25],RSCU[gennum][26],RSCU[gennum][27],
                    RSCU[gennum][28],RSCU[gennum][29],RSCU[gennum][30],RSCU[gennum][31],
                    RSCU[gennum][32],RSCU[gennum][33],RSCU[gennum][34],RSCU[gennum][36],
                    RSCU[gennum][37],RSCU[gennum][38],RSCU[gennum][39],RSCU[gennum][40],
                    RSCU[gennum][41],RSCU[gennum][42],RSCU[gennum][43],RSCU[gennum][44],
                    RSCU[gennum][45],RSCU[gennum][46],RSCU[gennum][47],RSCU[gennum][48],
                    RSCU[gennum][49],RSCU[gennum][50],RSCU[gennum][51],RSCU[gennum][52],
                    RSCU[gennum][53],RSCU[gennum][54],RSCU[gennum][55],RSCU[gennum][56],
                    RSCU[gennum][57],RSCU[gennum][58],RSCU[gennum][59],RSCU[gennum][60],
                    RSCU[gennum][61],RSCU[gennum][62],RSCU[gennum][63]);
return 0;
}
