/******************************
 * common.h
 * utility functions
 *****************************/

#ifndef COMMON
#define COMMON

#include <stdio.h>
#include <cstring>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>

using namespace std;

#define MAX_TP_DIST -0.0091
#define MIN_TP_DIST -0.0033

const double EulerConstant = exp(1.0);

/* STRUCTS */
/****************************************************/
/* structure of 3D points, denoting its coordinates */
typedef struct points{
  double x;
  double y;
  double z;
}POINT_S;

/***************************************/
/* interactions between two loci/beads */
typedef struct interactions{
  int loc1; //row
  int loc2; //column
  double num_freq; //the vaule of row loc1 and column loc2
	double native_dist; //distance
  //INTERACT_S *next;
}INTERACT_S;

/* End STRUCTS */

/* FUNCTIONS */

void chartoint(char *, vector<int> &);

/* print error message */
void print_error(const char *);

/* check if a file exist */
bool file_exist(const char *);

/* given a file, get the number of lines */
int get_num_lines(char* );

/* given a file including a matrix, get the number of nonzero values in the matrix */
int get_num_nonzero(char* );

/* randomize a point in the unit sphere */
void rand_point(POINT_S* );

/* get a new point nearby a given point */
void get_nearby_point(POINT_S* , POINT_S* );

/* given two points in the 3D space, calculate their Euclidean distance */
double cal_euc_distance(POINT_S* , POINT_S* );

/* given a distance D of two points, calculate the native distance F */
double cal_nat_distance(double , double , double );

void get_distance_matrix(char* , vector< vector<double> > &, double);

#endif

