
#include "IpIpoptApplication.hpp"
#include "hs071_nlp.hpp"

#include <iostream>

using namespace Ipopt;

int main(int argc, char* argv[])
{
  // create a new instance of your NLP
  // use a SmartPtr, not raw
  SmartPtr<TNLP> mynlp = new HS071_NLP();  
  // create a new instance of IpoptApplication
  // use a SmartPtr, nor raw
  // we are using the factory, since this allows us to compile this example
  // with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->RethrowNonIpoptException(true);

  // change some options
  // note: only for this example
  app->Options()->SetNumericValue("tol", 1e-7);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");

  // initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if(status != Solve_Succeeded){
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int)status;
  }

  // ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);
 
  if(status == Solve_Succeeded){
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
  }else{
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
  }
 
  /* as the SmartPtrs go out of scope, the reference count will be decremented
   * and the objects will automatically be deleted. */  

  return (int) status;
}
