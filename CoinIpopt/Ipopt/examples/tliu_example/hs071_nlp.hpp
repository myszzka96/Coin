
#ifndef __HS071_NLP_HPP__
#define __HS071_NLP_HPP__

#include "IpTNLP.hpp"

using namespace Ipopt;

class HS071_NLP : public TNLP
{
public:
  /* default constructor */
  HS071_NLP();

  /* default destructor */
  virtual ~HS071_NLP();

  // overloaded from TNLP
 
  /* give Ipopt the information about the size of the problem */ 
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, 
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /* give Ipopt the value of the bounds on the variables and constraints */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, 
                           Index m, Number* g_l, Number* g_u); 

  /* give Ipopt the starting point before it begins iterating */
  virtual bool get_starting_point(Index n, bool init_x, Number* x, 
                                  bool init_z, Number* z_L, Number* z_U, 
                                  Index m, bool init_lambda, Number* lambda);

  /* return the value of the objective function at the point x */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /* return the gradient of the objective function at the point x */ 
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /* return the constraint function at the point x */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /* return either the sparsity structure of the Jacobian of the constraints, 
   * or the values for the Jacobian of the constraints at the point x. */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x, 
                          Index m, Index nele_jac, Index* iRow, Index* jCol, 
                          Number* values);

  /* return either the sparsity structure of the Hessian of the Lagrangian,
   * or the values of the Hessian of the Lagrangian for the given values for x, obj_factor, and lambda1/2 */
  virtual bool eval_h(Index n, const Number* x, bool new_x, 
                      Number obj_factor, Index m, const Number* lambda, 
                      bool new_lambda, Index nele_hess, Index* iRow, 
                      Index* jCol, Number* values);

  /* This method is called by Ipopt after the algorithm has finished.
   * (successfully or even with more errors) */
  virtual void finalize_solution(SolverReturn status, 
                                 Index n, const Number* x, const Number* z_L, 
                                 const Number* z_U, Index m, const Number* g, 
                                 const Number* lambda, Number obj_value, 
                                 const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   r  you want (for all but the most simple classes), we usually
r: expected unqualified-id before 'using'
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  //@{
  //  HS071_NLP();
  HS071_NLP(const HS071_NLP&);
  HS071_NLP& operator=(const HS071_NLP&);
  //@}

};

#endif

