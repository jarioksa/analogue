/* Implement Gower distance with Podani amendment to ordered factors. */

/* Standard R headers */

#include <R.h>
#include <Rmath.h>
#include <math.h>

/* The variable type (vtype) must be one of the following
   integers. Type PRESENCE is the same as "asymmetric binary" and
   implements a Jaccard-like disismilarity. Ordered factors are
   handled either as their internal codes (METRICINTERNAL) like in
   cluster::daisy(), as their ranks (METRICORDINAL) or as ranks with
   Podani's tie correction. Normal quantitative variables are
   METRIC. TIE is a special case: the input matrix must be amended
   with columns which give number of tied values for each observation
   of each ORDINAL variable.
*/

#define BINARY 1
#define PRESENCE 2
#define NOMINAL 3
#define ORDINAL 4
#define METRICORDINAL 5
#define METRICINTERNAL 6
#define METRIC 7
#define TIE 8

/* Function adds two arguments to the normal call. Both of these are
   vectors of lenght nc:

   vtype: variable type as #define'd above

   scale: variable to scale each column to equal weight in Gower
   distance. For most vtype's this is 1, but for METRIC, METRICORDINAL
   and METRICINTERNAL it is the range, and for ORDINAL and TIE it is
   the Podani denominator. These must be calculated in the calling
   program.

   NB. This does not have weights but adding these is trivial (or they
   can be incorporated in scale).

*/

double xx_gowerpodani(double *x, int nr, int nc, int i1, int i2,
		      int *vtype, double *scale)
{
    double dist;
    int j, count;

    dist = 0.0;
    count = 0;

    for (j=0; j<nc; j++ ) {
	if (R_FINITE(x[i1]) && R_FINITE(x[i2])) {
	    switch(vtype[j]) {
	    case METRICORDINAL:
	    case METRICINTERNAL:
	    case METRIC:
	    case ORDINAL:
		dist += fabs(x[i1]-x[i2])/scale[j];
		count++;
		break;
	    case BINARY:
	    case PRESENCE:
		if (x[i1] > 0 || x[i2] > 0) {
		    dist += fabs(x[i1]-x[i2])/scale[j];
		    count++;
		}
		break;
	    case NOMINAL:
		dist += (x[i1] == x[i2]) ? 0.0 : 1.0 ;
		count++;
		break;
	    case TIE:
		dist -= ((x[i1] + x[i2]) / 2.0 - 1.0) / scale[j];
		break;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    return dist / (double) count;
}

/* A note on implementation:
   
   Handling of PRESENCE only differs in treating double zeros: these
   do not update count and hence missing species do not reduce
   dissimilarity (i.e., do not increase count in return
   dist/count). PRESENCE does not mean that data should be 0/1: they
   can quite well be quantitative and hence we divide with scale. The
   point is exactly in treating double-zeros.

   Most other variables are treated identically as Manhattan: the only
   difference comes from scale. only NOMINAL variables are handled
   with matching coefficient. If NOMINAL variables are replaced with
   their contrast matrix, they could be handled as PRESENCE.

   The handling of Podani tie correction is more opaque. It very much
   ties around defining scale and making scale correction. For other
   variables the scale is plain range, but in Podani correction the
   scale is range - (Timax-1)/2 - (Timin-1)/2 where Timax and Timin
   are number of items tied to max and min of the variable. This is a
   constant to each variable and can be supplied as scale. The
   equation 2b of Podani 1999 can be rearranged as 

   (Manhattan - Tiecorrection)/scale = 
            Manhattan/scale - Tiecorrection/scale 
   
   where

   Tiecorrection = (Tik-1)/2 + (Tjk-1)/2, 

   where T is the number of tied values for rows i and j in variable
   k. That is the column we must add to the input data. 

 */
