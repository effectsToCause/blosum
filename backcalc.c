/*=======================================================================
    BACKCALC

    Reads an sij matrix assumed to be in bit units and back-calculates
    qij using Blosum 62 marginal pi:
	qij = pipje^(Lsij), L = ln2
-------------------------------------------------------------------------
  6/20/94  J. Henikoff 
=========================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define NO 0
#define YES 1
#define AAS 20
#define MAXLINE 240

void read_file();
double entropy();

/*---- Global scoring matrix , order is :
  A R N D C Q E G H I L K M F P S T W Y V B Z J   -----------*/

char Alphabet[23]={'A','R','N','D','C','Q','E','G','H','I','L','K',
		   'M','F','P','S','T','W','Y','V','B','Z','J'};

/*--------Blosum50 62% clustering marginal AA frequencies--------------*/
double Blosum[20] = {0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054,
	  0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057,
	  0.051, 0.013, 0.032, 0.073};
/*=======================================================================*/
void main(argc, argv)
int argc;
char *argv[];
{
   char filename1[80];
   FILE *fin;
   double qij[AAS][AAS], sij[AAS][AAS];

   if (argc > 1) strcpy(filename1, argv[1]);
   else
   {
      printf("\nEnter name of file containing bit scores (sij): " );
      gets(filename1);
   }
   if ((fin = fopen(filename1, "rt")) == NULL)
   {
      printf("\nCould not open %s", filename1);
      exit(-1);
   }
   else
      printf("\nReading %s ...\n", filename1);
   read_file(fin, sij);
   fclose(fin);

   printf("\nH = % 8.4f\n", entropy(qij, sij, Blosum) );

   exit(0);
}
/*======================================================================*/
void read_file(fin, scores)
FILE *fin;
double scores[AAS][AAS];
{
   char line[MAXLINE], *ptr;
   int alpha[AAS], nrows, ncols, row, col, i;

/*----------Read file until first non-blank line --------------*/
/* Skip comments at beginning of file - 1st char = #, > or ;   */
   line[0] = '\0';
   while ((strlen(line) < 1 || line[0]=='#' || line[0]=='>' || line[0]==';')
          && fgets(line, sizeof(line), fin) != NULL)
	    ;
/*------See if the first line has characters on it ------------*/
   for (col=0; col < AAS; col++) alpha[col] = -1;
   if (strstr(line, "A") != NULL)	/* This line has characters */
   {
      row = 0;	/* # of alphabetic characters on the line */
      for (i=0; i<strlen(line); i++)
      {
	 col = -1;
	 if (line[i] == 'A') col = 0;
	 if (line[i] == 'R') col = 1;
	 if (line[i] == 'N') col = 2;
	 if (line[i] == 'D') col = 3;
	 if (line[i] == 'C') col = 4;
	 if (line[i] == 'Q') col = 5;
	 if (line[i] == 'E') col = 6;
	 if (line[i] == 'G') col = 7;
	 if (line[i] == 'H') col = 8;
	 if (line[i] == 'I') col = 9;
	 if (line[i] == 'L') col = 10;
	 if (line[i] == 'K') col = 11;
	 if (line[i] == 'M') col = 12;
	 if (line[i] == 'F') col = 13;
	 if (line[i] == 'P') col = 14;
	 if (line[i] == 'S') col = 15;
	 if (line[i] == 'T') col = 16;
	 if (line[i] == 'W') col = 17;
	 if (line[i] == 'Y') col = 18;
	 if (line[i] == 'V') col = 19;
	 if (col >= 0)
	 {
	    alpha[row] = col;
	    row++;
	 }
	 else if (isalpha(line[i])) row++;
      }
   }
/*-------Get the data values now ------------*/
   for (row=0; row<AAS; row++)
     for (col=0; col<AAS; col++)
	scores[row][col] = -999;		/* Null value */
   nrows = 0;
   line[0] = '\0';
   while (fgets(line, sizeof(line), fin) != NULL)
   {
      if (strlen(line) > 1 && nrows < AAS)
      {
	 if (alpha[nrows] >= 0 && alpha[nrows] < AAS)
	 {
	    row = alpha[nrows]; ncols = 0;
	    ptr = strtok(line, " ,\n");
	    while (ptr != NULL)
	    {
	       if (strspn(ptr, ".+-0123456789") == strlen(ptr))
	       {
		  col = alpha[ncols];
		  if (col >= 0 && col < AAS)
		     scores[row][col] = (double) atof(ptr);
		  ncols++;
	       }
	       ptr = strtok(NULL, " ,\n");
	    }
	 }
	 nrows++;
      }
   }

/*-------If some entries are still missing, assume symmetry ---------*/
   for (row=0; row<AAS; row++)
     for (col=0; col<AAS; col++)
	if (scores[row][col] == -999) scores[row][col] = scores[col][row];

}  /* end of read_file */
/*========================================================================*/
double entropy(qij, sij, pi)
double qij[AAS][AAS], sij[AAS][AAS], pi[AAS];
{
   int row, col;
   double entropy, lambda, dtemp;

   lambda = log(2.0);
   entropy = 0.0;
   for (row = 0; row < AAS; row++)
      for (col = 0; col < AAS; col++)
      {
	   dtemp = exp(lambda * sij[row][col]);
	   qij[row][col] = pi[row] * pi[col] * dtemp;
	   entropy += qij[row][col] * sij[row][col];
      }

   return (entropy);
}  /* end of entropy */
