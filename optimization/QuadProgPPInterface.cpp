#include "QuadProgPPInterface.h"
#include <iostream>
using namespace Optimization;
using namespace std;

#if HAVE_QUADPROGPP

#include "QuadProg++.hh"
typedef double RowVector [MATRIX_DIM];

LinearProgram::Result QuadProgPPInterface::Solve(const QuadraticProgram& qp,Vector& xopt)
{
  if(qp.Pobj.m > MATRIX_DIM) {
    cerr<<"QuadProgPPInterface: QPs must not exceed "<<MATRIX_DIM<<" dimensions"<<endl;
    return LinearProgram::Error;
  }

  int n=qp.Pobj.m;
  Matrix Aeq,Aineq;
  Vector beq,bineq;
  qp.GetSimpleForm(Aeq,beq,Aineq,bineq);
  int p=Aeq.m;
  int m=Aineq.m;
  if(p > MATRIX_DIM) {
    cerr<<"QuadProgPPInterface: QPs must not exceed "<<MATRIX_DIM<<" dimensions"<<endl;
    return LinearProgram::Error;
  }
  if(m > MATRIX_DIM) {
    cerr<<"QuadProgPPInterface: QPs must not exceed "<<MATRIX_DIM<<" dimensions"<<endl;
    return LinearProgram::Error;
  }
  RowVector* G = new RowVector[n];
  double* g0 = new double[n];
  RowVector* CE = new RowVector[n];
  double* ce0 = new double[p];
  RowVector* CI = new RowVector[n];
  double* ci0 = new double[m];
  double* x = new double[n];

  //copy CE and CI, negative for ce0, ci0
  for(int i=0;i<n;i++) {
    g0[i] = qp.qobj(i);
    for(int j=0;j<n;j++)
      G[i][j] = qp.Pobj(i,j);
  }
  for(int i=0;i<p;i++) {
    ce0[i] = -beq(i);
    for(int j=0;j<n;j++)
      CE[j][i] = Aeq(i,j);
  }
  for(int i=0;i<m;i++) {
    ci0[i] = -bineq(i);
    for(int j=0;j<n;j++)
      CI[j][i] = Aineq(i,j);
  }

  printf("Calling solve_quadprog...\n");
  double res=solve_quadprog(G,g0,n, 
			    CE,ce0,p, 
			    CI,ci0,m,
			    x);
  printf("Done with solve_quadprog\n");

  delete [] G;
  delete [] g0;
  delete [] CE;
  delete [] ce0;
  delete [] CI;
  delete [] ci0;

  if(!IsInf(std::numeric_limits<double>::infinity())) {
    cout<<"Warning, numeric limits infinity is not infinity..."<<endl;
    getchar();
  }

  if(IsInf(res) || res== std::numeric_limits<double>::infinity()) {
    delete [] x;
    return LinearProgram::Infeasible;
  }
  else {
    xopt.resize(n);
    for(int i=0;i<n;i++)
      xopt(i) = x[i];
    return LinearProgram::Feasible;
  }
}

#else

LinearProgram::Result QuadProgPPInterface::Solve(const QuadraticProgram& qp,Vector& xopt)
{
  cout<<"QuadProgPPInterface: not defined"<<endl;
  return LinearProgram::Error;
}

#endif // HAVE_QUADPROGPP
