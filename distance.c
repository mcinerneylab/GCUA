int distances( void )
{
int i, j, k;
FILE *outf;
float dist = 0.0;
char lin1[36];
/*#ifdef MAC
if(distance = 1){
_fcreator = 'PAUP';
_ftype = 'TEXT';
}
else {
_fcreator = 'R*CH';
_ftype = 'TEXT';
}
#endif*/
if(!loadedfile){
printf("Please load a file first.\n\n");
return 0;
}
if(distance == 1)printf("Name of (PAUP) file to write distances");
else printf("Name of (PHYLIP) file to write distances");
getstr("", lin1);
  if((outf = fopen(lin1, "w")) == NULL){
  printf("Cannot open %s for writing\n\n", lin1);
  return 0;
  }
if(calcs != 1)  printf("Opened %s...\n", lin1);
if(distance == 1){
  fprintf(outf, "#NEXUS\n");
  fprintf(outf, "[! Pairwise distances derived from RSCU values of the genes in \'%s\']\n", InputFilnam);
  fprintf(outf, "begin distances;\n");
  fprintf(outf, "dimensions\n");
  fprintf(outf, "ntax = %d;\n", nbseq);
  fprintf(outf, "matrix");
  }
else fprintf(outf, "%5d", nbseq);
for(i = 0; i < nbseq; i++){ /*Calculate pairwise distances between genes*/
 fprintf(outf, "\n%-10s", genname[i]);
 for(j = 0; j <= i; j++){
  for(k = 0; k < 64; k++){
  if( k !=  10 || k != 11 || k!= 14 || k!= 15 || k != 35 ){
   dist += fabs(RSCU[i][k] - RSCU[j][k]);
/*  printf("%f\t%f\t\tdist:%f\n", RSCU[i][k], RSCU[j][k], dist);*/
   }
  }
  if(dist == 0.0) dist = 0.00001;
  dist = dist / 59;
 fprintf(outf, "%8.5f", dist);
  dist = 0.0;
 }
}
if(distance == 1){
fprintf(outf, "\n;\n");
fprintf(outf, "END;\n");
}
fclose(outf);
if(calcs != 1) printf("Closed %s...\n\n", lin1);
return 0;
}
