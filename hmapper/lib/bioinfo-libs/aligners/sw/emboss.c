#include "emboss.h"


//-------------------------------------------------------------------------

int emboss_cnt = 0;
extern double emboss_matrix_t, emboss_tracking_t;

//-------------------------------------------------------------------------

float sw_emboss(char* seq_a, char* seq_b, float gapopen, float gapextend, 
		char* m, char* n, int* start1, int* start2,
		double *matrix_time, double *tracking_time, double *total_time) {
  int len_a = strlen(seq_a);
  int len_b = strlen(seq_b);
  int total = (len_a + 1) * (len_b + 1);
  
  float score;
  float* path = (float*) calloc(total, sizeof(float));
  int* compass = (int*) calloc(total, sizeof(int));

  double partial_t, partial1_t;
  partial1_t = sw_tic();
  partial_t = sw_tic();
  score = AlignPathCalcSW(seq_a, seq_b, len_a, len_b, gapopen, gapextend, path, compass);
  *matrix_time += sw_toc(partial_t);

  partial_t = sw_tic();
  AlignWalkSWMatrix(path, compass, gapopen, gapextend, seq_a, seq_b,
  		    (char*) m, (char *) n, len_a, len_b, start1, start2);
  *tracking_time += sw_toc(partial_t);
  
  *total_time += sw_toc(partial1_t);

  free(path);
  free(compass);

  return score;
}

/* @func embAlignPathCalcSW ***************************************************
**
** Create path matrix for Smith-Waterman
** Nucleotides or proteins as needed.
**
** @param [r] a [const char *] first sequence
** @param [r] b [const char *] second sequence
** @param [r] lena [ajint] length of first sequence
** @param [r] lenb [ajint] length of second sequence
** @param [r] gapopen [float] gap opening penalty
** @param [r] gapextend [float] gap extension penalty
** @param [w] path [float *] path matrix
** @param [w] compass [ajint *] Path direction pointer array
**
** @return [float] Maximum score
** @@
** Optimised to keep a maximum value to avoid looping down or left
** to find the maximum. (il 29/07/99)
******************************************************************************/

float AlignPathCalcSW(const char *a, const char *b, int lena, int lenb,
                      float gapopen, float gapextend, float* path, int* compass) {
  float ret;
  long xpos;
  long ypos;
  long i;
  long j;

  double match;
  double mscore;
  double result;
  double fnew;
  double* maxa;

  double bx;
  char compasschar;

  ret = -FLT_MAX;

  /* Create stores for the maximum values in a row or column */
  maxa = (double*) calloc(lena, sizeof(double));

  /* First initialise the first column and row */
  for (i = 0; i < lena; ++i) {
    result = (a[i]==b[0] ? 5.0 : -4.0);

    fnew = i==0 ? 0. : path[(i-1)*lenb] -(compass[(i-1)*lenb]==DOWN ? gapextend : gapopen);

    if (result > fnew && result > 0) {
      path[i*lenb] = (float) result;
      compass[i*lenb] = DIAGONAL;
    } else if (fnew > 0) {
      path[i*lenb] = (float) fnew;
      compass[i*lenb] = DOWN;
    } else {
      path[i*lenb] = 0.;
      compass[i*lenb] = ZERO;
    }

    maxa[i] = i==0 ? path[i*lenb]-gapopen : path[i*lenb] - (compass[(i-1)*lenb]==DOWN ? gapextend : gapopen);
  }

  for (j = 0; j < lenb; ++j) {
    result = (a[0]==b[j] ? 5.0 : -4.0); 
    fnew = j==0 ? 0. : path[j-1] -(compass[j-1]==LEFT ? gapextend : gapopen);

    if (result > fnew && result > 0) {
      path[j] = (float) result;
      compass[j] = DIAGONAL;
    } else if (fnew >0) {
      path[j] = (float) fnew;
      compass[j] = LEFT;
    } else {
      path[j] = 0.;
      compass[j] = ZERO;
    }
  }

  /* xpos and ypos are the diagonal steps so start at 1 */
  xpos = 1;
  float aux;

  while (xpos!=lenb) {
    ypos  = 1;
    bx = path[xpos]-gapopen-gapextend;

    while (ypos < lena) {
      /* get match for current xpos/ypos */
      match = (a[ypos]==b[xpos] ? 5.0 : -4.0);

      /* Get diag score */
      mscore = path[(ypos-1)*lenb+xpos-1] + match;
      aux = mscore;

      /* Set compass to diagonal value 0 */
      compass[ypos*lenb+xpos] = DIAGONAL;
      path[ypos*lenb+xpos] = (float) mscore;


      /* Now parade back along X axis */
      maxa[ypos] -= gapextend;
      fnew=path[(ypos)*lenb+xpos-1];
      fnew-=gapopen;


      if (fnew > maxa[ypos]) {
        maxa[ypos] = fnew;
      }

      if (maxa[ypos] > mscore) {
        mscore = maxa[ypos];
        path[ypos*lenb+xpos] = (float) mscore;

        /* Score comes from left */
        compass[ypos*lenb+xpos] = LEFT;
      }

      /* And then bimble down Y axis */
      bx -= gapextend;
      fnew = path[(ypos-1)*lenb+xpos];
      fnew-=gapopen;

      if (fnew > bx) {
        bx = fnew;
      }

      if (bx > mscore) {
        mscore = bx;
        path[ypos*lenb+xpos] = (float) mscore;

        /* Score comes from bottom */
        compass[ypos*lenb+xpos] = DOWN;
      }



      if (mscore > ret) {
        ret = (float) mscore;
      }

      result = path[ypos*lenb+xpos];

      if (result < 0.) {
        path[ypos*lenb+xpos] = 0.;
        compass[ypos*lenb+xpos] = ZERO;
      }

      ypos++;
    }

    ++xpos;
  }

  free(maxa);
  return ret;
}

//-------------------------------------------------------------------------

void AlignWalkSWMatrix(const float* path, const int* compass,
		       float gapopen, float gapextend, const char*  a, const char* b,
		       char* m, char* n, int lena, int lenb, int *start1, int *start2) {
  long i;
  long j;
  long k;
  long gapcnt;
  double pmax;
  double score;
  double bimble;

  long ix;
  long iy;
  
  long xpos = 0;
  long ypos = 0;
  const char *p;
  const char *q;
  
  int ic;
  double errbounds = (double) 0.01;

  /* Get maximum path score and save position */
  pmax = -FLT_MAX;
  k = (long)lena*(long)lenb-1;

  for (i = lena-1; i>=0; --i) {
    for (j=lenb-1; j>=0; --j) {
      if ((path[k--] > pmax) || E_FPEQ(path[k+1],pmax,U_FEPS)) {
        pmax = path[k+1];
	      xpos = j;
	      ypos = i;
      }
    }
  }

  p = a;
  q = b;

  while (xpos >= 0 && ypos >= 0) {
    /* diagonal */
    if (compass[ypos*lenb+xpos] == DIAGONAL) {
      strncat(m, &p[ypos--], 1);
      strncat(n, &q[xpos--], 1);

      if (ypos >= 0 && xpos>=0 && path[(ypos)*lenb+xpos] <= 0.) {
        break;
      }

      continue;
    }
    /* Left, gap(s) in vertical */
    else if (compass[ypos*lenb+xpos]==LEFT) {
      score  = path[ypos*lenb+xpos];
      gapcnt = 0;
      ix     = xpos-1;

      while (1) {
        bimble = path[ypos*lenb+ix]-gapopen-(gapcnt*gapextend);

        if(!ix || fabs((double)score-(double)bimble) < errbounds) {
          break;
        }

        --ix;
        ++gapcnt;
      }

      if (bimble <= 0.0) {
        break;
      }

      for (ic = 0; ic <= gapcnt; ++ic) {
        strcat(m, "-");
        strncat(n, &q[xpos--], 1);
      }

      continue;
    } else if (compass[ypos*lenb+xpos] == DOWN) {
      /* Down, gap(s) in horizontal */
      score  = path[ypos*lenb+xpos];
      gapcnt = 0;
      iy = ypos-1;

      while (1) {
        bimble = path[iy*lenb+xpos]-gapopen-(gapcnt*gapextend);

        if (!iy || fabs((double)score-(double)bimble) < errbounds) {
          break;
        }

        --iy;

        if (iy < 0) {
        printf("SW: Error walking down");
        exit(-1);
        }

        ++gapcnt;
      }

      if (bimble <= 0.0) {
        break;
      }

      for (ic = 0; ic <= gapcnt; ++ic) {
        strncat(m, &p[ypos--], 1);
        strcat(n, "-");
      }

      continue;
    } else {
      printf("Walk Error in SW");
      exit(-1);
    }
  }

  *start1 = (int) (ypos + 1); /* Potential lossy cast */
  *start2 = (int) (xpos + 1); /* Potential lossy cast */


  revstr(m);
  revstr(n);
  return;
}

//------------------------------------------------------------------------

void revstr(char* str) {
  int i;
  int len = strlen(str);
  char cpstr[len+1];

  for(i=0; i < len ; i++) {
    cpstr[i] = str[len-i-1];
  }

  cpstr[i] = '\0';
  strcpy(str, cpstr);
}
