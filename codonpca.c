void pca ( void)
{
FILE *pcaout;
int   n = nbseq, 
      m = 59, 
      i, 
      j,
      l,
      j1, 
      j2, 
      k, 
      k2,
     ii = 0, 
     jj = 0, 
     kk = 0, 
     nn = 0, 
     mm = 0,
     DNA = true,
     screen = true,
     file = true;
float **data, **matrix(), **symmat, **symmat2, *vector(), *evals, *interm;
void free_matrix(), free_vector(), corcol(), covcol(), scpcol();
void tred2(), tqli();
float in_value;
char option = 'r',  lin1[36];
/*#ifdef MAC
_fcreator = 'CRGR';
_ftype = 'TEXT';
#endif*/
if(!loadedfile){
 printf("\n\tThere are no sequences in memory\n\n");
 return;
 }
if(code == 7){
 m = 61;
 }
RSCU_OR_AA( &DNA, &m); /*ask if the analysis is to be performed on DNA or AA data*/
COR_VAR_SSCP( &option ); /*Ask how to form matrix*/
       

    data = matrix(n, m);  /* Storage allocation for input data */
if(DNA){
   printf("No. of Genes: %d, no. of Codons: %d.\n", n, m);
    for (i = 1; i <= n; i++){
       for (j = 1, j2 = 0; j <= m, j2 < 64; j++, j2++){
       if(code == 1){
          while(j2 == 10 || j2 == 11 || j2 == 14 || j2 == 15 || j2 == 35){
            j2++;
            }
           }
        if(code == 7){
          while(j2 == 10 || j2 == 11 || j2 == 35){
            j2++;
            }
           }
            data[i][j] = RSCU[i-1][j2];
            }
           }
          }
else { /*if analysis is to be perfromed on Amino Acids*/
  printf("No. of Genes: %d, no. of Amino Acids %d.\n", n, m);
   for(i = 1; i <= n; i++){
    for(j = 1, j2 = 0; j <= m, j2 < m+3; j++, j2++){
     while(j2 == 4 || j2 == 5 || j2 == 7){
      j2++;
      }
     data[i][j] = AminoArray[j2][i-1];
     }
    }
   }
  symmat = matrix(m, m);  /* Allocation of correlation (etc.) matrix */

     switch(option) { /* Look at analysis option; branch in accordance with this. */
          case 'R':
          case 'r':
              printf("Analysis of correlations chosen.\n");
              corcol(data, n, m, symmat);
                           /*Output correlation matrix.*/
              if(calcs == 3){
                          for (i = 1; i <= m; i++) {
                           for (j = 1; j <= 8; j++)  {
                            printf("%7.4f", symmat[i][j]);  }
                            printf("\n");  }  
                           }                        
              break;
          case 'V':
          case 'v':
              printf("Analysis of variances-covariances chosen.\n");
              covcol(data, n, m, symmat);

                          /* Output variance-covariance matrix.*/
              if(calcs == 3){
                          for (i = 1; i <= m; i++) {
                          for (j = 1; j <= 8; j++)  {
                            printf("%7.1f", symmat[i][j]);  }
                            printf("\n");  }  
                           }                        
              break;
          case 'S':
          case 's':
              printf("Analysis of sums-of-squares-cross-products");
              printf(" matrix chosen.\n");
              scpcol(data, n, m, symmat);
                         /* Output SSCP matrix.*/
            if(calcs == 3){
                        for (i = 1; i <= m; i++) {
                          for (j = 1; j <= 8; j++)  {
                            printf("%7.1f", symmat[i][j]);  }
                            printf("\n");  
                           }
                          }
              break;              
         case 'C':
         case 'c':
            printf("Correspondence Analysis\n");
            ca(data, n, m, symmat/*, Cmean*/);                
                          /*Output CA matrix.*/
              if(calcs == 3){
                         for (i = 1; i <= m; i++) {
                          for (j = 1; j <= 8; j++)  {
                            printf("%7.4f", symmat[i][j]);  
                           }
                            printf("\n");  
                          }
                         }
               break;
          default:
              printf("Option: %c\n",option);
              printf("For option, please type R, V, or S\n");
              printf("(upper or lower case).\n");
              printf("Exiting to system.\n");
              exit(1);
              break;
          }

/*********************************************************************
    Eigen-reduction
**********************************************************************/

    /* Allocate storage for dummy and new vectors. */
    evals = vector(m);     /* Storage alloc. for vector of eigenvalues */
    interm = vector(m);    /* Storage alloc. for 'intermediate' vector */
    symmat2 = matrix(m, m);  /* Duplicate of correlation (etc.) matrix */
    for (i = 1; i <= m; i++) {
     for (j = 1; j <= m; j++) {
      symmat2[i][j] = symmat[i][j]; /* Needed below for col. projections */
      }
     }
    tred2(symmat, m, evals, interm);  /* Triangular decomposition */
    if(option != 'c'){
    tqli(evals, interm, m, symmat); /* Reduction of sym. trid. matrix */
    }
    else{
    tql2(evals, interm, m, symmat);
     for(j1 = 1; j1 < m; j1++){
       for(j2 = 1; j2 < 4; j2++){
        symmat[j1][j2] /= sqrt(Cmean[j1]);
       }
      }
     }
    /* evals now contains the eigenvalues,
       columns of symmat now contain the associated eigenvectors. */
    if(calcs == 3){
     printf("\nEigenvalues:\n");
     for (j = m; j >= 1; j--) {
       printf("%18.5f\n", evals[j]);
       }
     printf("\n(Eigenvalues should be strictly positive; limited\n");
     printf("precision machine arithmetic may affect this.\n");
     printf("Eigenvalues are often expressed as cumulative\n");
     printf("percentages, representing the 'percentage variance\n");
     printf("explained' by the associated axis or principal component.)\n");
     printf("\nEigenvectors:\n");
     printf("(First three; their definition in terms of original vbes.)\n");
     for (j = 1; j <= m; j++) {
       for (i = 1; i <= 3; i++)  {
          printf("%12.4f", symmat[j][m-i+1]);
         }
        printf("\n");
        }
       }
     /* Form projections of row-points on first three prin. components. */
     /* Store in 'data', overwriting original data. */
if(option != 'c'){
     for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
        interm[j] = data[i][j]; }   /* data[i][j] will be overwritten */
        for (k = 1; k <= 4; k++) {
          data[i][k] = 0.0;
          for (k2 = 1; k2 <= m; k2++) {
            data[i][k] += interm[k2] * symmat[k2][m-k+1]; }
        }
     }
    }
else{
for(k = 1; k <= n; k++){
 for(l = 1; l <= m; l++){
  interm[l] = data[k][l];
  }
 for(k2 = 1; k2 <= 4; k2++){
  data[k][k2] = 0.0;
 for(j = 1; j <= m; j++){
  data[k][k2] += (interm[j] * symmat[j][m-k2]);
  }
  if(Rmean[k] > 0.0) data[k][k2] /= Rmean[k];
  if(Rmean[k] == 0.0) data[k][k2] = 0.0;
  }
 }
}
     printf("\n\n\n\n\n\n\n\n\n\n");
     printf("%s", BANNER);
     printf("\t1.             Send the output to the screen\n");
     printf("\t2.             Send the output to a file\n");
     printf("\t3.             Send the output to both the screen and a file\n\n\n");
     printf("\tPress [Return] to go to main menu\n\n\n\n");
     
     
        getstr("enter your choice", lin1);

		switch(toupper(*lin1)) {
			case '1': screen = true;
			          file = false;
						break;
			case '2': screen = false;
			          file = true;
			            break;
			case '3': screen = true;
			          file = true;
			            break;
			default : return;
			}
if(file){
       getstr("Name of file to print multivariate analysis ", lin1);
	   if((pcaout = fopen(lin1, "w")) == NULL){  /*Check to see if it can be opened for writing*/
	   printf("\n\nCannot open %s\n\n", lin1);
	   return;
       }
     fprintf(pcaout, "Gene\tAx1\tAx2\tAx3\tAx4\n");
     for (i = 1; i <= n; i++) {
      fprintf(pcaout, "%-10s\t", genname[i-1]);
       for (j = 1; j <= 4; j++)  {
          fprintf(pcaout, "%12.4f\t", data[i][j]);
          }
         fprintf(pcaout, "\n");
        }
   fclose(pcaout);
  }
if(screen){
     printf("\nProjections of Genes on first 4 prin. comps.:\n");
     for (i = 1; i <= n; i++) {
      printf("%-10s\t", genname[i-1]);
       for (j = 1; j <= 4; j++)  {
          printf("%12.4f", data[i][j]);
          }
          printf("\n");  
          }
        }
     /* Form projections of col.-points on first three prin. components. */
     /* Store in 'symmat2', overwriting what was stored in this. */
if(option != 'c'){
     for (j = 1; j <= m; j++) {
      for (k = 1; k <= m; k++) {
        interm[k] = symmat2[j][k]; }  /*symmat2[j][k] will be overwritten*/
        for (i = 1; i <= 3; i++) {
          symmat2[j][i] = 0.0;
          for (k2 = 1; k2 <= m; k2++) {
            symmat2[j][i] += interm[k2] * symmat[k2][m-i+1]; }
          if (evals[m-i+1] > 0.0005)   /* Guard against zero eigenvalue */
             symmat2[j][i] /= sqrt(evals[m-i+1]);   /* Rescale */
          else
             symmat2[j][i] = 0.0;    /* Standard kludge */
        }
     }
    }
else{
for(j = 1; j <= m; j++){
 for(k = 1; k <= m; k++){
  interm[k] = symmat2[j][k];
  }
  for(i = 1; i <= 4; i++){
   symmat2[j][i] = 0.0;
   for(k2 = 1; k2 <= m; k2++){
    symmat2[j][i] += interm[k2] * symmat[k2][m-i] * sqrt(Cmean[k2]);
    }
   }
   if(evals[m-i] > 0.0 && Cmean[j] > 0.0){
    symmat2[j][i] /= sqrt(evals[m-i] * Cmean[j]);
    }
   if(evals[m-i] == 0.0 || Cmean[j] == 0.0){
    symmat2[j][i] = 0.0;
    }
   }
  }
if(screen){
     printf("\nProjections of Codons on first 4 prin. comps.:\n");
     for (j = 1; j <= m; j++) {
       for (k = 1; k <= 4; k++)  {
          printf("%12.4f", symmat2[j][k]);
         }
        printf("\n");
       }
}
    free_matrix(data, n, m);
    free_matrix(symmat, m, m);
    free_matrix(symmat2, m, m);
    free_vector(evals, m);
    free_vector(interm, m);
}

/**  Correlation matrix: creation  ***********************************/
void corcol( float **data, int n, int m, float **symmat )
/* Create m * m correlation matrix from given n * m data matrix. */
{
float eps = 0.005;
float x, *mean, *stddev, *vector();
int i, j, j1, j2;

/* Allocate storage for mean and std. dev. vectors */
mean = vector(m);
stddev = vector(m);

/* Determine mean of column vectors of input data matrix */
for (j = 1; j <= m; j++){
    mean[j] = 0.0;
    for (i = 1; i <= n; i++){
        mean[j] += data[i][j];
        }
    mean[j] /= (float)n;
    }
/* Determine standard deviations of column vectors of data matrix. */
for (j = 1; j <= m; j++){
    stddev[j] = 0.0;
    for (i = 1; i <= n; i++){
        stddev[j] += (   ( data[i][j] - mean[j] ) *
                         ( data[i][j] - mean[j] )  );
        }
        stddev[j] /= (float)n;
        stddev[j] = sqrt(stddev[j]);
        if (stddev[j] <= eps) stddev[j] = 1.0;
    }
/* Center and reduce the column vectors. */
for (i = 1; i <= n; i++){
    for (j = 1; j <= m; j++){
        data[i][j] -= mean[j];
        x = sqrt((float)n);
        x *= stddev[j];
        data[i][j] /= x;
        }
    }
/* Calculate the m * m correlation matrix. */
for (j1 = 1; j1 <= m-1; j1++){
    symmat[j1][j1] = 1.0;
    for (j2 = j1+1; j2 <= m; j2++){
        symmat[j1][j2] = 0.0;
        for (i = 1; i <= n; i++){
            symmat[j1][j2] += ( data[i][j1] * data[i][j2]);
            }
        symmat[j2][j1] = symmat[j1][j2];
        }
    }
    symmat[m][m] = 1.0;
return;
}

/**  Variance-covariance matrix: creation  *****************************/
void covcol( float **data, int n, int m, float **symmat)
/* Create m * m covariance matrix from given n * m data matrix. */
{
float *mean, *vector();
int i, j, j1, j2;

/* Allocate storage for mean vector */
mean = vector(m);

/* Determine mean of column vectors of input data matrix */
for (j = 1; j <= m; j++){
    mean[j] = 0.0;
    for (i = 1; i <= n; i++){
        mean[j] += data[i][j];
        }
    mean[j] /= (float)n;
    }

/* Center the column vectors. */
for (i = 1; i <= n; i++){
    for (j = 1; j <= m; j++){
        data[i][j] -= mean[j];
        }
       }

/* Calculate the m * m covariance matrix. */
for (j1 = 1; j1 <= m; j1++){
    for (j2 = j1; j2 <= m; j2++){
        symmat[j1][j2] = 0.0;
        for (i = 1; i <= n; i++){
            symmat[j1][j2] += data[i][j1] * data[i][j2];
            }
        symmat[j2][j1] = symmat[j1][j2];
        }
    }
return;
}

/**  Sums-of-squares-and-cross-products matrix: creation  **************/
void scpcol(float **data, int n, int m, float **symmat)
/* Create m * m sums-of-cross-products matrix from n * m data matrix. */
{
int i, j1, j2;

/* Calculate the m * m sums-of-squares-and-cross-products matrix. */
for (j1 = 1; j1 <= m; j1++)
    {
    for (j2 = j1; j2 <= m; j2++)
        {
        symmat[j1][j2] = 0.0;
        for (i = 1; i <= n; i++)
            {
            symmat[j1][j2] += data[i][j1] * data[i][j2];
            }
        symmat[j2][j1] = symmat[j1][j2];
        }
    }
return;
}

/* Calculate m * m Correspondence Analysis matrix ********************/
int ca(float **dmatrix, int n, int m, float **A)
{
int i, j, j1, j2;
float tot = 0.0;

for(i = 1; i <= n; i++){ /*form row sums and total**/
 Rmean[i] = 0.0;
 for(j = 1; j <= m; j++){
  tot += dmatrix[i][j];
  Rmean[i] += dmatrix[i][j];
  }
 }
 for(j = 1; j <= m; j++){ /*form column sums and then means*/
  Cmean[j] = 0.0;
   for(i = 1; i <= n; i++){
    Cmean[j] += dmatrix[i][j];
    }
    if(Cmean[j] == 0.0){
      printf("All values in column %d are 0.0 and must be removed\nEXITING\n", j);
    return 0;
    }
    Cmean[j] /= (float)tot;
   }  
   for(i = 1; i <= n; i++){ /*form row means and make dmatrix into frequencies*/
    if(Rmean[i] == 0.0){
      printf("All values in row %d are 0.0 and must be removed\nEXITING\n", i);
    return 0;
    }
   Rmean[i] /= (float)tot;
    for(j = 1; j <= m; j++){
 dmatrix[i][j] = dmatrix[i][j]/(float)tot;
    }
   }
for(j1 = 1; j1 <= m; j1++){
 for(j2 = 1; j2 <= m; j2++){
  A[j1][j2] = 0.0;
   for(i = 1; i <= n; i++){
    A[j1][j2] += dmatrix[i][j1] * dmatrix[i][j2] / (Rmean[i] * sqrt(Cmean[j1] * Cmean[j2]));
    }
   }
  } 
 return 0;
}
/**  Error handler  **************************************************/

void erhand(char *err_msg){
    printf("\tRun-time error:\n");
    printf("\t%s\n", err_msg);
    printf("Press [Return] to continue ");
    getchar();
    return;
    /*exit(1);*/
}

/**  Allocation of vector storage  ***********************************/

float *vector(int n)
/* Allocates a float vector with range [1..n]. */
{

    float *v;

    v = (float *) malloc ((unsigned) n*sizeof(float));
    if (!v) erhand("Allocation failure in vector().");
    return v-1;

}

/**  Allocation of float matrix storage  *****************************/

float **matrix( int n, int m)
/* Allocate a float matrix with range [1..n][1..m]. */
{
    int i;
    float **mat;

    /* Allocate pointers to rows. */
    mat = (float **) malloc((unsigned) (n)*sizeof(float*));
    if (!mat) erhand("Allocation failure 1 in matrix().");
    mat -= 1;

    /* Allocate rows and set pointers to them. */
    for (i = 1; i <= n; i++)
        {
        mat[i] = (float *) malloc((unsigned) (m)*sizeof(float));
        if (!mat[i]) erhand("Allocation failure 2 in matrix().");
        mat[i] -= 1;
        }

     /* Return pointer to array of pointers to rows. */
     return mat;

}

/**  Deallocate vector storage  *********************************/

void free_vector( float *v, int n )
/* Free a float vector allocated by vector(). */
{
   free((char*) (v+1));
}

/**  Deallocate float matrix storage  ***************************/

void free_matrix(float **mat, int n, int m)
/* Free a float matrix allocated by matrix(). */
{
   int i;

   for (i = n; i >= 1; i--)
       {
       free ((char*) (mat[i]+1));
       }
   free ((char*) (mat+1));
}

/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */

void tred2(float **a, int n, float *d, float *e)

/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
{
int l, k, j, i;
float scale, hh, h, g, f;
if(calcs == 3) printf("Reducing data matrix to tridiagonal form\n");
for (i = n; i >= 2; i--)
    {
    l = i - 1;
    h = scale = 0.0;
    if (l > 1)
       {
       for (k = 1; k <= l; k++)
           scale += fabs(a[i][k]);
       if (scale == 0.0)
          e[i] = a[i][l];
       else
          {
          for (k = 1; k <= l; k++)
              {
              a[i][k] /= scale;
              h += a[i][k] * a[i][k];
              }
          f = a[i][l];
          g = f>0 ? -sqrt(h) : sqrt(h);
          e[i] = scale * g;
          h -= f * g;
          a[i][l] = f - g;
          f = 0.0;
          for (j = 1; j <= l; j++)
              {
              a[j][i] = a[i][j]/h;
              g = 0.0;
              for (k = 1; k <= j; k++)
                  g += a[j][k] * a[i][k];
              for (k = j+1; k <= l; k++)
                  g += a[k][j] * a[i][k];
              e[j] = g / h;
              f += e[j] * a[i][j];
              }
          hh = f / (h + h);
          for (j = 1; j <= l; j++)
              {
              f = a[i][j];
              e[j] = g = e[j] - hh * f;
              for (k = 1; k <= j; k++)
                  a[j][k] -= (f * e[k] + g * a[i][k]);
              }
         }
    }
    else
        e[i] = a[i][l];
    d[i] = h;
    }
d[1] = 0.0;
e[1] = 0.0;
for (i = 1; i <= n; i++)
    {
    l = i - 1;
    if (d[i])
       {
       for (j = 1; j <= l; j++)
           {
           g = 0.0;
           for (k = 1; k <= l; k++)
               g += a[i][k] * a[k][j];
           for (k = 1; k <= l; k++)
               a[k][j] -= g * a[k][i];
           }
       }
       d[i] = a[i][i];
       a[i][i] = 1.0;
       for (j = 1; j <= l; j++)
           a[j][i] = a[i][j] = 0.0;
    }
}

/**  Tridiagonal QL algorithm -- Implicit  **********************/

void tqli(float d[], float e[], int n, float **z)
{
int m, l, iter, i, k;
float s, r, p, g, f, dd, c, b;
void erhand();

for (i = 2; i <= n; i++)
    e[i-1] = e[i];
    e[n] = 0.0;
for (l = 1; l <= n; l++)
    {
    iter = 0;
    do
      {
      for (m = l; m <= n-1; m++)
          {
          dd = fabs(d[m]) + fabs(d[m+1]);
          if (fabs(e[m]) + dd == dd) break;
          }
          if (m != l){
             if (iter++ == 30) {erhand("No convergence in TLQI."); return;}
             g = (d[l+1] - d[l]) / (2.0 * e[l]);
             r = sqrt((g * g) + 1.0);
             g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
             s = c = 1.0;
             p = 0.0;
             for (i = m-1; i >= l; i--){
                 f = s * e[i];
                 b = c * e[i];
                 if (fabs(f) >= fabs(g)){
                    c = g / f;
                    r = sqrt((c * c) + 1.0);
                    e[i+1] = f * r;
                    c *= (s = 1.0/r);
                    }
                 else {
                    s = f / g;
                    r = sqrt((s * s) + 1.0);
                    e[i+1] = g * r;
                    s *= (c = 1.0/r);
                    }
                 g = d[i+1] - p;
                 r = (d[i] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i+1] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++) {
                     f = z[k][i+1];
                     z[k][i+1] = s * z[k][i] + c * f;
                     z[k][i] = c * z[k][i] - s * f;
                     }
                 }
                 d[l] = d[l] - p;
                 e[l] = g;
                 e[m] = 0.0;
             }
          }  while (m != l);
      }
 }


int tql2( float *d, float *e,int n, float **z )
{
    /* Initialized data */
    float eps = 1e-12f;
    /* System generated locals */
    int i__1, i__2, i__3;
    float r__1, r__2;

    /* Local variables */
    float b, c, f, g, h;
    int i, iter, k, l, m, j;
    float p, r, s;
    int l1, ii, mml;

    /* Function Body */
if (n == 1) {
   return 0;
    }
   for (i = 2; i <= n; i++) {
	e[i - 1] = e[i];
    }
    f = 0.f;
    b = 0.f;
    e[n] = 0.f;
   for (l = 1; l <= n; l++) {
	iter = 0;
	h = eps * ((r__1 = d[l], abs(r__1)) + (r__2 = e[l], fabs(r__2)));
	if (b < h) {
	    b = h;
	}
do {  /*Look for small sub-diagonal element. */
	for (m = l; m <= n-1; m++) {
	    if ((r__1 = e[m], fabs(r__1)) <= b) {
		break;
	    }
	}
if(m!= l){
    	if (iter++ == 30) { erhand("No convergence after 30 iterations"); return 0;	}
    	if(calcs == 3) printf("Iteration %d\n", iter);
	l1 = l + 1; /*Form shift. */
	g = d[l];
	p = (d[l1] - g) / (e[l] * 2.f);
	r = sqrt(p * p + 1.f);
	d[l] = e[l] / (p + SIGN(r, p));
	h = g - d[l];
  for (i = l1; i <= n; i++) {
	    d[i] -= h;
	}
	f += h;
	p = d[m]; /*QL transformation. */
	c = 1.f;
	s = 0.f;
	mml = m - l;
	for (ii = 1; ii <= mml; ++ii) {
	    i = m - ii;
	    g = c * e[i];
	    h = c * p;
	    if (fabs(p) >= (r__1 = e[i], fabs(r__1))) {
	    c = e[i] / p;
	    r = sqrt(c * c + 1.f);
	    e[i + 1] = s * p * r;
	    s = c / r;
	    c = 1.f / r;
        }
    else{
	    c = p / e[i];
	    r = sqrt(c * c + 1.f);
	    e[i + 1] = s * e[i] * r;
	    s = 1.f / r;
	    c *= s;
        }
	    p = c * d[i] - s * g;
	    d[i + 1] = h + s * (c * g + s * d[i]);
	    for (k = 1; k <= n; k++) { /*Form vector. */
		h = z[k][i + 1];
		z[k][i + 1] = s * z[k][i] + c * h;
		z[k][i] = c * z[k][i] - s * h;
	    }
	}
	e[l] = s * p;
	d[l] = c * p;
	if ((r__1 = e[l], fabs(r__1)) <= b) {
	    break;
	}
   }
  }while (m != l);
	d[l] += f;
} 
     for (ii = 2; ii <= n; ii++) { /*Order eigenvectors and eigenvalues. */
	i = ii - 1;
	k = i;
	p = d[i];
	for (j = ii; j <= n; j++) {
	    if (d[j] < p) {
	    k = j;
	    p = d[j];
	   }
	 }
	if (k != i) {
	d[k] = d[i];
	d[i] = p;
	for (j = 1; j <= n; j++) {
	    p = z[j][i];
	    z[j][i] = z[j][k];
	    z[j][k] = p;
   }
  }
 }
 return 0;
} /* tql2 */


