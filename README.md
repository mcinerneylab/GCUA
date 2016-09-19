# GCUA
GCUA: General Codon Usage Analysis

I have made this code available as open source.  It is provided as-is and needs to be compiled by the user.

compilation is easiest using the GCC compiler by typing:

gcc codon.c -o GCUA

The memory requirements are not particularly onerous right now and should be fine on most modern computers (~3Gb).
The program can handle "genes" of length 2 million bp and it can handle 5 million genes in a single analysis.

If this memory requirement is too big, then you can either change this line in codon.h:

&#35;define MAXGENLEN 2000000     /*Max length of a gene*/

- if you dont have genes of length 2 million bp

or this line

&#35;define MAXNUMGEN 5000000      /*Max number of genes in a single analysis*/

 - if you have fewer than 5 million genes.


and recompile the program.
