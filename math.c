/******sumcodons calculates the number of times 
that a particular amino acid is used in a particular gene.  
If the amino acid is not used at all, then its value is set to -1
This is done so that it can be easily identified later on for RSCU calculations
an so on*******/

int sumcodons(int gennum)
{
int yy;
      if((AminoArray[0][gennum] = (codonstot[gennum][0] + codonstot[gennum][1])) == 0)                                                                                                    {AminoArray[0][gennum] = -1;}
      if((AminoArray[1][gennum] = (codonstot[gennum][2] + codonstot[gennum][3] + codonstot[gennum][16] + codonstot[gennum][17] + codonstot[gennum][18] + codonstot[gennum][19])) == 0)    {AminoArray[1][gennum] = -1;}
      if((AminoArray[2][gennum] = (codonstot[gennum][4] + codonstot[gennum][5] + codonstot[gennum][6] + codonstot[gennum][7] + codonstot[gennum][44] + codonstot[gennum][45])) == 0)      {AminoArray[2][gennum] = -1;}
      if((AminoArray[3][gennum] = (codonstot[gennum][8] + codonstot[gennum][9])) == 0)                                                                                                    {AminoArray[3][gennum] = -1;}
          AminoArray[4][gennum] = codonstot[gennum][10];
          AminoArray[5][gennum] = codonstot[gennum][11];
      if((AminoArray[6][gennum] = (codonstot[gennum][12] + codonstot[gennum][13])) == 0)                                                                                                  {AminoArray[6][gennum] = -1;}
          AminoArray[7][gennum] = codonstot[gennum][14];
 if(code == 1){      
          AminoArray[8][gennum] = codonstot[gennum][15];
         }
if(code == 7){
      if((AminoArray[8][gennum] = (codonstot[gennum][15] + codonstot[gennum][14])) == 0)                                                                                                  {AminoArray[8][gennum] = -1;}
         }
      if((AminoArray[9][gennum] = (codonstot[gennum][20] + codonstot[gennum][21] + codonstot[gennum][22] + codonstot[gennum][23])) == 0)                                                  {AminoArray[9][gennum] = -1;}
      if((AminoArray[10][gennum] = (codonstot[gennum][24] + codonstot[gennum][25])) == 0)                                                                                                 {AminoArray[10][gennum] = -1;}
      if((AminoArray[11][gennum] = (codonstot[gennum][26] + codonstot[gennum][27])) == 0)                                                                                                 {AminoArray[11][gennum] = -1;}
      if((AminoArray[12][gennum] = (codonstot[gennum][28] + codonstot[gennum][29] + codonstot[gennum][30] + codonstot[gennum][31] + codonstot[gennum][46] + codonstot[gennum][47])) == 0) {AminoArray[12][gennum] = -1;}
      if((AminoArray[13][gennum] = (codonstot[gennum][32] + codonstot[gennum][33] + codonstot[gennum][34])) == 0)                                                                         {AminoArray[13][gennum] = -1;}
      if((AminoArray[14][gennum] = codonstot[gennum][35]) == 0)                                                                                                                           {AminoArray[14][gennum] = -1;}
      if((AminoArray[15][gennum] = (codonstot[gennum][36] + codonstot[gennum][37] + codonstot[gennum][38] + codonstot[gennum][39])) == 0)                                                 {AminoArray[15][gennum] = -1;}
      if((AminoArray[16][gennum] = (codonstot[gennum][40] + codonstot[gennum][41])) == 0)                                                                                                 {AminoArray[16][gennum] = -1;}
      if((AminoArray[17][gennum] = (codonstot[gennum][42] + codonstot[gennum][43])) == 0)                                                                                                 {AminoArray[17][gennum] = -1;}
      if((AminoArray[18][gennum] = (codonstot[gennum][48] + codonstot[gennum][49] + codonstot[gennum][50] + codonstot[gennum][51])) == 0)                                                 {AminoArray[18][gennum] = -1;}
      if((AminoArray[19][gennum] = (codonstot[gennum][52] + codonstot[gennum][53] + codonstot[gennum][54] + codonstot[gennum][55])) == 0)                                                 {AminoArray[19][gennum] = -1;}
      if((AminoArray[20][gennum] = (codonstot[gennum][56] + codonstot[gennum][57])) == 0)                                                                                                 {AminoArray[20][gennum] = -1;}
      if((AminoArray[21][gennum] = (codonstot[gennum][58] + codonstot[gennum][59])) == 0)                                                                                                 {AminoArray[21][gennum] = -1;}
      if((AminoArray[22][gennum] = (codonstot[gennum][60] + codonstot[gennum][61] + codonstot[gennum][62] + codonstot[gennum][63])) == 0)                                                 {AminoArray[22][gennum] = -1;}
for(yy = 0; yy < 23; yy++){
 if(AminoArray[yy][gennum] != -1){
  AminoArray[yy][MAXNUMGEN] += AminoArray[yy][gennum];
   }
  }
  if(gennum == MAXNUMGEN){
  donesumcodons = true;
  }
return 0;
}
int RSCUfunc (int gennum)
{
int ii = 0;
      RSCU[gennum][0] = ((codonstot[gennum][0] * 2) / AminoArray[0][gennum]);
      RSCU[gennum][1] = ((codonstot[gennum][1] * 2) / AminoArray[0][gennum]);
      RSCU[gennum][2] = ((codonstot[gennum][2] * 6) / AminoArray[1][gennum]);
      RSCU[gennum][3] = ((codonstot[gennum][3] * 6) / AminoArray[1][gennum]);
      RSCU[gennum][4] = ((codonstot[gennum][4] * 6) / AminoArray[2][gennum]);
      RSCU[gennum][5] = ((codonstot[gennum][5] * 6) / AminoArray[2][gennum]);
      RSCU[gennum][6] = ((codonstot[gennum][6] * 6) / AminoArray[2][gennum]);
      RSCU[gennum][7] = ((codonstot[gennum][7] * 6) / AminoArray[2][gennum]);
      RSCU[gennum][8] = ((codonstot[gennum][8] * 2) / AminoArray[3][gennum]);
      RSCU[gennum][9] = ((codonstot[gennum][9] * 2) / AminoArray[3][gennum]);
      RSCU[gennum][10] = 0.;
      RSCU[gennum][11] = 0.;
      RSCU[gennum][12] = ((codonstot[gennum][12] * 2) / AminoArray[6][gennum]);
      RSCU[gennum][13] = ((codonstot[gennum][13] * 2) / AminoArray[6][gennum]);
if(code == 1){ 
      RSCU[gennum][14] = 0.;
      RSCU[gennum][15] = 1;
      }
if(code == 7){
      RSCU[gennum][14] = ((codonstot[gennum][14] * 2) / AminoArray[8][gennum]);
      RSCU[gennum][15] = ((codonstot[gennum][15] * 2) / AminoArray[8][gennum]);
      }
      RSCU[gennum][16] = ((codonstot[gennum][16] * 6) / AminoArray[1][gennum]);
      RSCU[gennum][17] = ((codonstot[gennum][17] * 6) / AminoArray[1][gennum]);
      RSCU[gennum][18] = ((codonstot[gennum][18] * 6) / AminoArray[1][gennum]);
      RSCU[gennum][19] = ((codonstot[gennum][19] * 6) / AminoArray[1][gennum]);
      RSCU[gennum][20] = ((codonstot[gennum][20] * 4) / AminoArray[9][gennum]);
      RSCU[gennum][21] = ((codonstot[gennum][21] * 4) / AminoArray[9][gennum]);
      RSCU[gennum][22] = ((codonstot[gennum][22] * 4) / AminoArray[9][gennum]);
      RSCU[gennum][23] = ((codonstot[gennum][23] * 4) / AminoArray[9][gennum]);
      RSCU[gennum][24] = ((codonstot[gennum][24] * 2) / AminoArray[10][gennum]);
      RSCU[gennum][25] = ((codonstot[gennum][25] * 2) / AminoArray[10][gennum]);
      RSCU[gennum][26] = ((codonstot[gennum][26] * 2) / AminoArray[11][gennum]);
      RSCU[gennum][27] = ((codonstot[gennum][27] * 2) / AminoArray[11][gennum]);
      RSCU[gennum][28] = ((codonstot[gennum][28] * 6) / AminoArray[12][gennum]);
      RSCU[gennum][29] = ((codonstot[gennum][29] * 6) / AminoArray[12][gennum]);
      RSCU[gennum][30] = ((codonstot[gennum][30] * 6) / AminoArray[12][gennum]);
      RSCU[gennum][31] = ((codonstot[gennum][31] * 6) / AminoArray[12][gennum]);
      RSCU[gennum][32] = ((codonstot[gennum][32] * 3) / AminoArray[13][gennum]);
      RSCU[gennum][33] = ((codonstot[gennum][33] * 3) / AminoArray[13][gennum]);
      RSCU[gennum][34] = ((codonstot[gennum][34] * 3) / AminoArray[13][gennum]);
      RSCU[gennum][35] = ((codonstot[gennum][35] * 1) / AminoArray[14][gennum]);
      RSCU[gennum][36] = ((codonstot[gennum][36] * 4) / AminoArray[15][gennum]);
      RSCU[gennum][37] = ((codonstot[gennum][37] * 4) / AminoArray[15][gennum]);
      RSCU[gennum][38] = ((codonstot[gennum][38] * 4) / AminoArray[15][gennum]);
      RSCU[gennum][39] = ((codonstot[gennum][39] * 4) / AminoArray[15][gennum]);
      RSCU[gennum][40] = ((codonstot[gennum][40] * 2) / AminoArray[16][gennum]);
      RSCU[gennum][41] = ((codonstot[gennum][41] * 2) / AminoArray[16][gennum]);
      RSCU[gennum][42] = ((codonstot[gennum][42] * 2) / AminoArray[17][gennum]);
      RSCU[gennum][43] = ((codonstot[gennum][43] * 2) / AminoArray[17][gennum]);
      RSCU[gennum][44] = ((codonstot[gennum][44] * 6) / AminoArray[2][gennum]);
      RSCU[gennum][45] = ((codonstot[gennum][45] * 6) / AminoArray[2][gennum]);
      RSCU[gennum][46] = ((codonstot[gennum][46] * 6) / AminoArray[12][gennum]);
      RSCU[gennum][47] = ((codonstot[gennum][47] * 6) / AminoArray[12][gennum]);
      RSCU[gennum][48] = ((codonstot[gennum][48] * 4) / AminoArray[18][gennum]);
      RSCU[gennum][49] = ((codonstot[gennum][49] * 4) / AminoArray[18][gennum]);
      RSCU[gennum][50] = ((codonstot[gennum][50] * 4) / AminoArray[18][gennum]);
      RSCU[gennum][51] = ((codonstot[gennum][51] * 4) / AminoArray[18][gennum]);
      RSCU[gennum][52] = ((codonstot[gennum][52] * 4) / AminoArray[19][gennum]);
      RSCU[gennum][53] = ((codonstot[gennum][53] * 4) / AminoArray[19][gennum]);
      RSCU[gennum][54] = ((codonstot[gennum][54] * 4) / AminoArray[19][gennum]);
      RSCU[gennum][55] = ((codonstot[gennum][55] * 4) / AminoArray[19][gennum]);
      RSCU[gennum][56] = ((codonstot[gennum][56] * 2) / AminoArray[20][gennum]);
      RSCU[gennum][57] = ((codonstot[gennum][57] * 2) / AminoArray[20][gennum]);
      RSCU[gennum][58] = ((codonstot[gennum][58] * 2) / AminoArray[21][gennum]);
      RSCU[gennum][59] = ((codonstot[gennum][59] * 2) / AminoArray[21][gennum]);
      RSCU[gennum][60] = ((codonstot[gennum][60] * 4) / AminoArray[22][gennum]);
      RSCU[gennum][61] = ((codonstot[gennum][61] * 4) / AminoArray[22][gennum]);
      RSCU[gennum][62] = ((codonstot[gennum][62] * 4) / AminoArray[22][gennum]);
      RSCU[gennum][63] = ((codonstot[gennum][63] * 4) / AminoArray[22][gennum]);
      for(ii = 0; ii < 64; ii++){
       if(RSCU[gennum][ii] <= 0.0){
        RSCU[gennum][ii] = 0.0;
        }
       }
       if(gennum == MAXNUMGEN){ /*let the rest of the program know that the RSCU values have been calculated*/
       doneRSCU = true;
       }
return 0;
}
int enc(int gennum) { /*Effective number of codons*/
    static float codat[128] = { 10.f,1.f,1.f,1.f,11.f,
	    1.f,12.f,1.f,10.f,2.f,1.f,2.f,11.f,2.f,12.f,2.f,2.f,1.f,1.f,3.f,
	    21.f,1.f,21.f,2.f,2.f,2.f,1.f,4.f,21.f,3.f,19.f,1.f,2.f,3.f,4.f,
	    1.f,13.f,1.f,3.f,1.f,2.f,4.f,4.f,2.f,13.f,2.f,3.f,2.f,2.f,5.f,4.f,
	    3.f,14.f,1.f,3.f,3.f,2.f,6.f,4.f,4.f,14.f,2.f,3.f,4.f,9.f,1.f,5.f,
	    1.f,15.f,1.f,1.f,5.f,9.f,2.f,5.f,2.f,15.f,2.f,1.f,6.f,9.f,3.f,5.f,
	    3.f,16.f,1.f,3.f,5.f,20.f,1.f,5.f,4.f,16.f,2.f,3.f,6.f,6.f,1.f,
	    7.f,1.f,17.f,1.f,8.f,1.f,6.f,2.f,7.f,2.f,17.f,2.f,8.f,2.f,6.f,3.f,
	    7.f,3.f,18.f,1.f,8.f,3.f,6.f,4.f,7.f,4.f,18.f,2.f,8.f,4.f };
/* System generated locals */
    int i__1;
    float r__1;
/* Local variables */
    float aveb[6];
    int noaa;
    float sumb, b[18];
    int j, k, l;
    float cu[126], sp2, saa[21];
    int ico[21], ijk;
    float effco, pka2;

/*calculates effective number of codons as per Wright (1990) */
    for (j = 1; j <= 64; ++j) {
	ico[(int) codat[(j << 1) - 2] - 1] = codat[(j << 1) - 1];
    }
    for (j = 1; j <= 64; ++j) {
	 cu[(int) codat[(j << 1) - 2] + (int) codat[(j << 1) - 1] * 21 - 22] = (float) codonstot[gennum][j-1];
    }
    ijk = 9;
    for (k = 1; k <= 18; ++k) {
	 saa[k - 1] = 0.f;
	 sp2 = 0.f;
	 i__1 = ico[k - 1];
	for (l = 1; l <= i__1; ++l) {
	    saa[k - 1] += cu[k + l * 21 - 22];
	}
	if (saa[k - 1] <= 1.f) {
	    b[k - 1] = 99.f;
	} else {
	    i__1 = ico[k - 1];
	    for (l = 1; l <= i__1; ++l) {
		 if (cu[k + l * 21 - 22] == 0.f) {
		    pka2 = 0.f;
		} 
		else { /*Computing 2nd power */
		    r__1 = cu[k + l * 21 - 22] / saa[k - 1];
		    pka2 = r__1 * r__1;
		}
		sp2 += pka2;
	    }
	    b[k - 1] = (saa[k - 1] * sp2 - 1.f) / (saa[k - 1] - 1.f);
	  }
     }
/*HEX CODERS - SER LEU ARG */
    noaa = 3;
    sumb = 0.f;
    for (k = 1; k <= 3; ++k) {
	 if (b[k - 1] == 99.f) {
	    --noaa;
	 } 
	else {
	    sumb += b[k - 1];
	 }
    }
    if (sumb == 0.f || noaa == 0) {
	 effco = 0.f;
	return effco;
    }
    aveb[5] = sumb / (float) noaa;
/*QUARTETS PRO THR VAL ALA GLY */
    noaa = 5;
    sumb = 0.f;
    for (k = 4; k <= 8; ++k) {
	if (b[k - 1] == 99.f) {
	    --noaa;
	} else {
	    sumb += b[k - 1];
	 }
    }
    if (sumb == 0.f || noaa == 0) {
	effco = 0.f;
	return effco;
    }
    aveb[3] = sumb / (float) noaa;
/*DUETS PHE TYR CYS HIS GLN ASN LYS ASP GLU */
    noaa = 9;
    sumb = 0.f;
    for (k = 10; k <= 18; ++k) {
	 if (b[k - 1] == 99.f) {
	    --noaa;
	} 
	else {
	    sumb += b[k - 1];
	 }
    }
    if (sumb == 0.f || noaa == 0) {
	effco = 0.f;
	return effco;
    }
    aveb[1] = sumb / (float) noaa;
/* TRIPLET ILE ONLY */
    if (b[8] == 99.f || b[8] == 0.f) {
	aveb[2] = (aveb[1] + aveb[3]) * .5f;
    } 
    else {
	aveb[2] = b[8];
    }
/*  SUMMING UP FOR EffECTIVE NUMBER OF CODONS */
    effco = 9.f / aveb[1] + 2.f + 1.f / aveb[2] + 5.f / aveb[3] + 3.f / aveb[5];
    if (effco > 61.f) {
	effco = 61.f;
    }
    /*printf("%f\n", effco);*/
    return (float) effco;
} /* enc */

