int process (int i){
int goodcodon, 
    badcodon = 0,
    j,
    k;
    
GC[i] = 0; 
GC3[i] = 0, 
nonsyn[i] = 0, 
stop[i] = 0;
          for(j = 0; j < length[i]; j++){
                if((j + 3) % 3 == 0){
           goodcodon = TRUE;
           k = 0;
           switch(toupper(s[j])){
            case 'T':            break;
            case 'C': k += 16;   GC[i]++; GC1[i]++; GC[MAXNUMGEN]++; GC1[MAXNUMGEN]++; break;
            case 'A': k += 32;   break;
            case 'G': k += 48;   GC[i]++; GC1[i]++; GC[MAXNUMGEN]++; GC1[MAXNUMGEN]++; break;
            default : goodcodon = FALSE;   break;
            }
           }
       else if((j + 2) % 3 == 0){
           switch(toupper(s[j])){
            case 'T':            break; 
            case 'C': k += 4;    GC[i]++; GC2[i]++; GC[MAXNUMGEN]++; GC2[MAXNUMGEN]++; break; 
            case 'A': k += 8;    break; 
            case 'G': k += 12;   GC[i]++; GC2[i]++; GC[MAXNUMGEN]++; GC2[MAXNUMGEN]++; break; 
            default : goodcodon = FALSE;    break;
            }
           } 
         else if((j + 1) % 3 == 0){
           switch(toupper(s[j])){
            case 'T':            break;
            case 'C': k += 1;    GC[i]++; GC3[i]++; GC[MAXNUMGEN]++; GC3[MAXNUMGEN]++; break;
            case 'A': k += 2;    break;
            case 'G': k += 3;    GC[i]++; GC3[i]++; GC[MAXNUMGEN]++; GC3[MAXNUMGEN]++; break;
            default : goodcodon = FALSE;    break;
            }
      if(goodcodon){
       codonstot[i][k]++;
       codonstot[MAXNUMGEN][k]++;
if(code == 1){
      if(k == 14 || k == 10 || k == 11){
       stop[i]++;
       stop[MAXNUMGEN]++;
       }
      if(k == 15 || k == 35){
        nonsyn[i]++;
        nonsyn[MAXNUMGEN]++;
        }
     }
if(code == 7){
      if( k == 10 || k == 11){
       stop[i]++;
       stop[MAXNUMGEN]++;
       }
       if(k == 35){
        nonsyn[i]++;
        nonsyn[MAXNUMGEN]++;
        }
      }
     }
      else badcodon += 1;
    }
   }
if(badcodon > 0){
             printf("\n%s has %4d bad codon", genname[i], badcodon); /*Tell about bad codons*/
if(badcodon > 1){
             printf("s\n");
              }
else printf("\n");
           badcodon = 0;
           }
}
