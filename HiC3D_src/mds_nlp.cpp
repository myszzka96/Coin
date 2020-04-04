/************************************
 * mds_nlp.cpp
 ************************************/
#include "mds_nlp.hpp"
#include <cassert>
#include <iostream>

using namespace Ipopt;

// constructor
MDS_NLP::MDS_NLP()
{}

// self-defined constructor
MDS_NLP::MDS_NLP(vector< vector<double> > &mat_distance, INTERACT_S* interactions, int num_nonzero, double max_th, char* output_file){

    this->interactions = interactions;
    this->num_interactions = num_nonzero;

    this->num_beads = (int)mat_distance.size();
    this->num_variables = 3*num_beads;

    //std::cout << "The number of beads is " << num_beads << std::endl;
    //std::cout << "The number of variables is " << num_variables << std::endl;
    this->max_th = max_th;
    this->output_file = output_file;

    //this->starting_point_file = "./mds.input.starting.points.txt";

    this->num_constraints_container = num_beads;
    this->num_constraints = this->num_constraints_container;
    //std::cout << "Constraints:" << std::endl;
    //std::cout << "The number of container constraints is " << num_constraints_container << std::endl; 
    //std::cout << "The total number of constraints is " << num_constraints << std::endl << std::endl;
}

// destructor
MDS_NLP::~MDS_NLP()
{}

// return the size of the problem
bool MDS_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    /* the number of variables for this problem */ 
    n = this->num_variables;

    /* the number of constraints */
    m = this->num_constraints;

    /* the number of nonzero elements in the Jacobian of the constraints */
    nnz_jac_g = 3*this->num_constraints_container;

    /* the Hessian of the Lagrangian function */
    // nnz_h_lag = ;

    /* if a quasi-Newton option is chosen to approximate the second derivatives, the above variable is not required */
  
    /* use the C style indexing (0-based) */  
    index_style = TNLP::C_STYLE;
  
    return true;
}

// return the variable bounds
bool MDS_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u, 
                            Index m, Number* g_l, Number* g_u)
{
    /* here, the n and m we gave Ipopt in get_nlp_info are passed back to us. */
    /* If desired, we could assert to make sure they are what we think they are. */
    assert(n == this->num_variables);
    assert(m == this->num_constraints);
 
    /* the variables have lower bound of -d and upper bound d */
    for (Index i=0; i < num_variables; i++) {
        x_l[i] = -100000;
				x_u[i] = 100000;
    }

    /* the functions of constraints lower/upper bounds */
    /* g_l[0..num_constraints], g_u[0..num_constraints] */
    for (Index i=0; i < this->num_constraints_container; i++) {
        g_l[i] = -2e19;
        g_u[i] = this->max_th*this->max_th;
    }
  
    /* the second constraint g2 */
  
    /* the third constraint g3 */

    return true;
}

bool MDS_NLP::get_starting_point(Index n, bool init_x, Number* x, 
                                   bool init_z, Number* z_L, Number* z_U, 
                                   Index m, bool init_lambda, Number* lambda)
{
    /* here, we assume we only have starting values for x, if you code your own 
     * NLP, you can provide starting values for the dual variables, 
     * if you wish to use a warmstart option */
  
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    assert(n == this->num_variables);

    /* first, set every x[i] to zero */
    for (Index i=0; i < num_variables; i++) {
        x[i] = 0;
    }
    
    POINT_S pre, post;
    for (Index i=0; i<1; i++) {
        rand_point(&pre);
        x[i] = pre.x;
        x[i+1] = pre.y;
        x[i+2] = pre.z;
    }
    for (Index i=3; i < num_variables; i+=3) {
        get_nearby_point(&pre, &post);
        x[i] = post.x;
        x[i+1] = post.y;
        x[i+2] = post.z;
        pre = post;
    }

    return true;
}

/*************************************************************/
/* return the value of the objective function at the point x */
bool MDS_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    assert(n == this->num_variables);
    obj_value = 0;
  
    for (Index i=0; i < this->num_interactions; i++) {
        int i1 = 3*this->interactions[i].loc1;
        int i2 = 3*this->interactions[i].loc2;
        double nat_dist = this->interactions[i].native_dist;
	      double w = 1.;
	      double obj_single = sqrt((x[i1]-x[i2])*(x[i1]-x[i2])+
                             (x[i1+1]-x[i2+1])*(x[i1+1]-x[i2+1])+
                             (x[i1+2]-x[i2+2])*(x[i1+2]-x[i2+2])) - nat_dist;
	      w = 1./(nat_dist*nat_dist);
	      obj_value += w * obj_single * obj_single;
     }
  
    return true;
}

/****************************************************************/
/* return the gradient of the objective function at the point x */
bool MDS_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	assert(n == this->num_variables);
  
	for (Index i=0; i < this->num_variables; i++) {
		grad_f[i] = 0;
	}
	for (Index i=0; i < this->num_interactions; i++) {
			int i1 = 3*this->interactions[i].loc1;
			int i2 = 3*this->interactions[i].loc2;
			double nat_dist = this->interactions[i].native_dist;
			double w = 1.;
			double tp_dist = sqrt((x[i1]-x[i2])*(x[i1]-x[i2])+
                          	(x[i1+1]-x[i2+1])*(x[i1+1]-x[i2+1])+
                          	(x[i1+2]-x[i2+2])*(x[i1+2]-x[i2+2]));
			w = 1./(nat_dist*nat_dist);
			if (tp_dist != 0) {
				grad_f[i1] += 2 * w * (tp_dist - nat_dist) * (x[i1] - x[i2]) / tp_dist;
				grad_f[i1+1] += 2 * w * (tp_dist - nat_dist) * (x[i1+1] - x[i2+1]) / tp_dist;
				grad_f[i1+2] += 2 * w * (tp_dist - nat_dist) * (x[i1+2] - x[i2+2]) / tp_dist;

				grad_f[i2] -= 2 * w * (tp_dist - nat_dist) * (x[i1] - x[i2]) / tp_dist;
				grad_f[i2+1] -= 2 * w * (tp_dist - nat_dist) * (x[i1+1] - x[i2+1]) / tp_dist;
				grad_f[i2+2] -= 2 * w * (tp_dist - nat_dist) * (x[i1+2] - x[i2+2]) / tp_dist;
			}
	}

	return true;
}

/**************************************************************/
/* return the value of the constraint function at the point x */
bool MDS_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	assert(n == this->num_variables);
	assert(m == this->num_constraints);
  
	/* the first constraint function */
	for (Index i=0; i < 3*this->num_constraints_container; i+=3) {
		g[i/3] = x[i]*x[i] + x[i+1]*x[i+1] + x[i+2]*x[i+2];
	}
	//int index_now = this->num_constraints_container;
  
	/* the second constraint function */

	return true;
}

/**********************************************************************************/
bool MDS_NLP::eval_jac_g(Index n, const Number* x, bool new_x, Index m, 
                           Index nele_jac, Index* iRow, Index* jCol, Number* values)
{
  
   int index = 0;
	/* if */
	/* return the structure of the Jacobian of constraints g(x) */
	/* the first constraint function g[0]=x[i]*x[i]+x[i+1]*x[i+1]+x[i+2]*x[i+2] */
	/* else */
	/* return the values of the Jacobian of the constraints g(x) */
	/* values for the Jacobian of g[0] */

	if (values == NULL) {
		for (Index i=0; i < 3*num_constraints_container; i+=3) {
			iRow[index] = i/3; jCol[index++] = i;
			iRow[index] = i/3; jCol[index++] = i+1;
			iRow[index] = i/3; jCol[index++] = i+2;
		}
	} else {
		for (Index i=0; i < 3*num_constraints_container; i+=3) {
			values[index++] = 2*x[i];
			values[index++] = 2*x[i+1];
			values[index++] = 2*x[i+2];
		}
	}

	return true;
}

/*******************************************************************************/
bool MDS_NLP::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, 
                       Index m, const Number* lambda, bool new_lambda, 
                       Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
	/* the structure or values of the Hessian */
  
	return false;
	/* denoting that a quasi-Newton is choosen to approximate the second derivatives, and the Hessian of Lagrangian is not required */
}

/*******************************************************************/
/* this method is called by IPOPT after the algorithm has finished */
void MDS_NLP::finalize_solution(SolverReturn status, Index n, const Number* x, 
       const Number* z_L, const Number* z_U, Index m, const Number* g, 
       const Number* lambda, Number obj_value, const IpoptData* ip_data, 
       IpoptCalculatedQuantities* ip_cq)
{
	// here is where we would store the solution to variables, or write to a file, etc
	// so we could use the solution.

	// For this example, we write the solution to the "output_file" file
	// std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
	if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
		ofstream myfile;
		myfile.open(output_file);  
		for (Index i=0; i < n; i++) { 
			myfile << x[i] << " ";
			if ((i+1) % 3 == 0){ 
				myfile << std::endl; 
			}
		}
		myfile.close();
	}
}


