#include "unistd.h"
#include "IpIpoptApplication.hpp"
#include "mds_nlp.hpp"

#define MAX_FILE_NAME 1000

using namespace Ipopt;

char distance_matrix_file[MAX_FILE_NAME]; //-d
char output_file[MAX_FILE_NAME]; //-o

double threshold = 1.0;
double beta = 1.0;

void parse_input_arguments(int, char **);

int main(int argc, char* argv[])
{
    parse_input_arguments(argc, argv);

    /*****************************/
    /* read distance matrix file */ 
    vector< vector<double> > mat_distance; 
    get_distance_matrix(distance_matrix_file, mat_distance, beta); 

    int num_nonzero = get_num_nonzero(distance_matrix_file);
    num_nonzero /= 2;
    INTERACT_S* interacts = new INTERACT_S[num_nonzero];
    
    int num_interacts = 0;
    for (unsigned int i=0; i < mat_distance.size(); i++) {
        for (unsigned int j=(i+1); j < mat_distance.at(i).size(); j++) {
            if (mat_distance.at(i).at(j)) {
                interacts[num_interacts].loc1 = i;
                interacts[num_interacts].loc2 = j;
                interacts[num_interacts++].native_dist = mat_distance.at(i).at(j);
            }
        }
    } 
		cout << "Start to run IPOPT..." << endl;
    /************************/
    /* start to run the job */
    SmartPtr<TNLP> mynlp = new MDS_NLP(mat_distance, interacts, num_nonzero, threshold, output_file);

    SmartPtr<IpoptApplication> app = new IpoptApplication();

    //Change some options
    app->Options()->SetNumericValue("tol", 1e-1);
    app->Options()->SetNumericValue("acceptable_tol", 1e-1);
    app->Options()->SetNumericValue("constr_viol_tol", 1e10);
    app->Options()->SetNumericValue("mu_init", 0.0001);
    app->Options()->SetNumericValue("dual_inf_tol", 10000);
    app->Options()->SetNumericValue("compl_inf_tol", 1000000);
    app->Options()->SetIntegerValue("acceptable_iter", 0);
    app->Options()->SetIntegerValue("max_iter", 100000);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    /* initialize the IpoptApplication and process the options */
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        printf("\n* Error during initialization!\n");
        return (int)status;
    }
    status = app->OptimizeTNLP(mynlp);
    if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
        printf("\n* Done! the problem solved!\n");
    } else {
        printf("\n* Error! the problem failed!\n");
        exit(1);
    }
  
    delete [] interacts;
	
    return (int)status;
}/* End main() */

void parse_input_arguments(int argc, char *argv[]) {
    int curr;
    opterr = 0;
    int num_argv = 0;
    while ((curr = getopt(argc, argv, "b:d:o:t:h")) != -1) {
        switch (curr) {
            case 'd': strcpy(distance_matrix_file, optarg); num_argv++; break;
            case 'o': strcpy(output_file, optarg); num_argv++; break;
            case 't': threshold = atof(optarg); break;
            case 'b': beta = atof(optarg); break;
            case '?': 
            case 'h':
            default:
                opterr = 1;
        }
    }
}

/* End */
