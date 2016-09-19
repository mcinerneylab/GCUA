/* ***************************ENC subroutine********************* */
int enc(integer *icu, real *effco)
{
    /* Initialized data */
    static float codat[128]	 = { 10.0,1.0,1.0,1.0,11.0,
	    1.0,12.0,1.0,10.0,2.0,1.0,2.0,11.0,2.0,12.0,2.0,2.0,1.0,1.0,3.0,
	    21.0,1.0,21.0,2.0,2.0,2.0,1.0,4.0,21.0,3.0,19.0,1.0,2.0,3.0,4.0,
	    1.0,13.0,1.0,3.0,1.0,2.0,4.0,4.0,2.0,13.0,2.0,3.0,2.0,2.0,5.0,4.0,
	    3.0,14.0,1.0,3.0,3.0,2.0,6.0,4.0,4.0,14.0,2.0,3.0,4.0,9.0,1.0,5.0,
	    1.0,15.0,1.0,1.0,5.0,9.0,2.0,5.0,2.0,15.0,2.0,1.0,6.0,9.0,3.0,5.0,
	    3.0,16.0,1.0,3.0,5.0,20.0,1.0,5.0,4.0,16.0,2.0,3.0,6.0,6.0,1.0,
	    7.0,1.0,17.0,1.0,8.0,1.0,6.0,2.0,7.0,2.0,17.0,2.0,8.0,2.0,6.0,3.0,
	    7.0,3.0,18.0,1.0,8.0,3.0,6.0,4.0,7.0,4.0,18.0,2.0,8.0,4.0 };

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static float aveb[6];
    static int noaa;
    static float sumb, b[18];
    static int j, k, l;
    static float cu[126], sp2, saa[21];
    static int ico[21], ijk;
    static float pka2;

/*  calculates effective number of codons as per Wright (1990) */

    for (j = 0; j < 64; j++) {
	 ico[(int) codat[j - 2]] = codat[j];
    }
    for (j = 0; j < 64; j++) {
	cu[(int) codat[j - 2] + (int) codat[j] * 21 - 22] = (float) icu[j];
    }
    ijk = 9;
    for (k = 0; k < 18; k++) {
	saa[k] = 0.0;
	sp2 = 0.0;
	for (l = 0; l < ico[k]; l++) {
	    saa[k] += cu[k + l * 21 - 22];
	}
	if (saa[k] <= 1.0) {
	    b[k] = 99.0;
	} 
	else {
	    for (l = 0; l < ico[k]; l++) {
		if (cu[k + l * 21 - 22] == 0.0) {
		    pka2 = 0.0;
		} 
		else {
/* Computing 2nd power */
		    r__1 = cu[k + l * 21 - 22] / saa[k - 1];
		    pka2 = r__1 * r__1;
		}
		sp2 += pka2;
	    }
	    b[k] = (saa[k] * sp2 - 1.0) / (saa[k] - 1.0);
	}
    }
/*  HEX CODERS - SER LEU ARG */
    noaa = 3;
    sumb = 0.0;
    for (k = 0; k < 3; k++) {
	if (b[k] == 99.0) {
	    --noaa;
	} 
	else {
	    sumb += b[k];
	}
    }
    if (sumb == 0.0 || noaa == 0) {
	*effco = 0.0;
	return;
    }
    aveb[5] = sumb / (real) noaa;
/*  QUARTETS PRO THR VAL ALA GLY */
    noaa = 5;
    sumb = 0.0;
    for (k = 3; k < 8; k++) {
	if (b[k] == 99.0) {
	    --noaa;
	} else {
	    sumb += b[k - 1];
	}
    }
    if (sumb == 0.0 || noaa == 0) {
	*effco = 0.0;
    return;
    }
    aveb[3] = sumb / (real) noaa;
/*  DUETS PHE TYR CYS HIS GLN ASN LYS ASP GLU */
    noaa = 9;
    sumb = 0.0;
    for (k = 9; k < 18; k++) {
	if (b[k] == 99.0) {
	    --noaa;
	} else {
	    sumb += b[k];
	}
    }
    if (sumb == 0.0 || noaa == 0) {
	*effco = 0.0;
	return;
    }
    aveb[1] = sumb / (real) noaa;
/* TRIPLET ILE ONLY */
    if (b[8] == 99.0 || b[8] == 0.0) {
	aveb[2] = (aveb[1] + aveb[3]) * 0.5f;
    } 
    else {
	aveb[2] = b[8];
    }
/*  SUMMING UP FOR EffECTIVE NUMBER OF CODONS */
    *effco = 9.0 / aveb[1] + 2.0 + 1.0 / aveb[2] + 5.0 / aveb[3] + 3.0 / aveb[5];
    if (*effco > 61.0) {
	*effco = 61.0;
    }
    return 0;
} /* enc */

