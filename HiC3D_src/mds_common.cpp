/******************************************
 * common.cpp
 * implementation of some common functions
 *****************************************/

#include "mds_common.h"

/* FUNCTIONS */
void chartoint(char *p, vector<int> &tmp) {
    char *tok;
    tok = strtok(p, ",");
    while (tok != NULL) {
        int index = atoi(tok);
        tmp.push_back(index);
        tok = strtok(NULL, ",");
    }
}

void print_error(const char *info_error) {
    printf("%s\n", info_error);
    exit(1);
}

bool file_exist(const char *file_name) {
    std::ifstream myfile(file_name);
    return myfile.good();
}

static double cal_distance(POINT_S* p) {
    return sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
}

static bool valid_point(POINT_S* p) {
    if (cal_distance(p) < 1.) {
      return true;
    }
    return false;
}

void rand_point(POINT_S* p) {
    do {
        p->x = (2.*drand48()) - 1;
        p->y = (2.*drand48()) - 1;
        p->z = (2.*drand48()) - 1;
    } while (!valid_point(p));
}

void get_nearby_point(POINT_S* pre, POINT_S* post) {
    static POINT_S r;
    do {
        rand_point(&r);
        r.x *= MAX_TP_DIST; 
        r.y *= MAX_TP_DIST;
        r.z *= MAX_TP_DIST;
        post->x = pre->x + r.x;
        post->y = pre->y + r.y;
        post->z = pre->z + r.z;
    } while (!valid_point(post));
}

double cal_euc_distance(POINT_S* p1, POINT_S* p2) {
    double tmp1 = p1->x - p2->x;
    double tmp2 = p1->y - p2->y;
    double tmp3 = p1->z - p2->z;
    return sqrt(tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3);
}

double cal_nat_distance(double d, double alpha, double beta) {
    if (d<=0) {
        print_error("Error, there is a frequency with zero value."); 
    }
    double tmp1 = (double)1/d;
    double tmp2 = (double)1/alpha;
    return beta*pow(tmp1, tmp2);
}

int get_num_lines(char* file_name) {
    ifstream myfile;
    int lines = 0;
    string line;
    myfile.open(file_name);
    if (!myfile.is_open()) {
        print_error("Error, when open the file to compute the number of lines.");
    }
    while (getline(myfile, line)) {
        lines++;
    }
    myfile.close();

    return lines;
}

int get_num_nonzero(char* file_name) { //off-diagonal
    ifstream myfile;
    int nonzero = 0;
    string line;
    double temp;
    myfile.open(file_name);
    if( !myfile.is_open()) {
        print_error("Error, when open the file to compute the number of lines.");
    }
    int rows = 0;
    while (getline(myfile, line)) {
        rows++;
        int cols = 0;
        istringstream istr(line);
        while (istr >> temp) {
            cols++;
            if (temp != 0 && rows != cols) {
                nonzero++;
            }
        }
        line.clear();
        istr.clear();
    }
    myfile.close();

    return nonzero;
}

void get_distance_matrix(char* distance_matrix_file, vector< vector<double> > &mat_distance, double beta) {
    ifstream myfile;
    string line;
    double temp;
    myfile.open(distance_matrix_file);
    if (!myfile.is_open()) {
        print_error("Error, when opening the interaction file.");
    }
    while (getline(myfile, line)) {
        istringstream istr(line);
        vector<double> tmp;
        while (istr >> temp) {
            temp *= beta;
            tmp.push_back(temp);
        }
        line.clear();
        istr.clear();
        mat_distance.push_back(tmp);
    }
    myfile.close();
}

// End
