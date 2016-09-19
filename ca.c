int ca(float data[][], int n, int m, float A[][])
{
int i, j, j1, j2;
float tot = 0.0;

for(i = 0; i < n; i++){ /*form row sums and total**/
 R[i] = 0.0;
 for(j = 0; j < m; j++){
  tot = (tot + data[i][j]);
  R[i] = R[i] + data[i][j];
  }
 }
 
 for(j = 0; j < m; j++){ /*form column sums and then means*/
  C[j] = 0.0;
   for(i = 0; i < n; i++){
    C[j] = C[j] + data[i][j];
    }
    if(C[j] !> 0.0){
    printf("All values in column %d are 0.0 and must be removed\nEXITING\n");
    return;
    }
    C[j] = C[j]/tot;
   }
   
   for(i = 0; i < n; i++){
    if(R[i] !> 0.0){
    printf("All values in row %d are 0.0 and must be removed\nEXITING\n");
    return;
    }
   R[i] = R[i]/tot;
    for(j = 0; j < m; j++){
     data[i][j] = data[i][j]/tot;
     }
    }
    
for(j1 = 0; j1 < m; j1++){
 for(j2 = 0; j2 < m; j2++){
  A[j1][j2] = 0.0;
   for(i = 0; i < n; i++){
    A[j1][j2] = A[j1][j2] + data[i][j1] * data[i][j2] / (R[i] * sqrt(C[j1] * C[j2]));
    }
    A[j2][j1] = A[j1][j2];
   }
   }
   return;
  }