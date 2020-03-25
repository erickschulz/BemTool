#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "bemtool/tools.hpp"
#include <fstream>

using namespace bemtool;


int main(int argc, char* argv[]){

  Real kappa = 0.5; // What does this mean?
  Real kappa2 = kappa*kappa;

  // Loading the nodes into the geometry object? This is how to read .msh files
  Geometry node("mesh/circle.msh");
  // Mesh1D object mesh? loading the node into the mesh?
  Mesh1D mesh; mesh.Load(node,1);
  // Orienting the mesh
  Orienting(mesh);
  // Number of elements in the mesh
  int nb_elt = NbElt(mesh);
  std::cout << "nb_elt:\t" << nb_elt << std::endl;

  // Defining the operator type (Laplacian double layer 2d)
  // P1XP1 most probably means trial and test spaces as seen in a few lines
  typedef LA_DL_2D_P1xP1 OperatorType;

  // Defining the BIO based on the operator type
  BIOp<OperatorType> V(mesh,mesh,kappa);
  // Dof object. Templated on the BEM space.
  Dof<P1_1D> dof(mesh);
  // Degrees of freedom. The BEM space info must have been used by now to get this value.
  int nb_dof = NbDof(dof);
  std::cout << "nb_dof:\t" << nb_dof << std::endl;
  // Initializing the Galerkin matrix A for the BIOp? A is most probably the Galerkin matrix
  EigenDense  A(nb_dof,nb_dof); Clear(A);

  const int n = 1;
  // Vector of complex numbers. Size = degrees of freedom. What does it represent?
  std::vector<Cplx> En(nb_dof);
  // Looping over the number of elements to fill En vector
  for(int j=0; j<nb_elt; j++){
    const N2&         jdof = dof[j];
    // Point in 2D
    const array<2,R3> xdof = dof(j);
    for(int k=0; k<2; k++){
      En[jdof[k]] = pow( xdof[k][0]+iu*xdof[k][1], n);}
  }

  /*
  progress bar("Assemblage",nb_elt);
  for(int j=0; j<nb_elt; j++){
    bar++;
    for(int k=0; k<nb_elt; k++){
      A(dof[j],dof[k]) += V(j,k);
    }
  }
  */

  // Assembly of the Galerkin matrix?
  progress bar("Assemblage",nb_dof);
  for(int j=0; j<nb_dof; j++){
    bar++;
    for(int k=0; k<nb_dof; k++){
      A(j,k) += V(dof.ToElt(j),dof.ToElt(k)); // dof to elt represents global to local map?
    }
  }
  bar.end();

  // What is this sum for? It represents the 'solution' as it is used to determine the error.
  Cplx sum =0.;
  for(int j=0; j<nb_dof; j++){
    for(int k=0; k<nb_dof; k++){
      sum+= A(j,k)*conj(En[j])*En[k]; // conj(En)^T A En
    }
  }

  //  Cplx refsol = RefSol<OperatorType>::Compute(n,1.,kappa);
  Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,1.,kappa);
  std::cout << "erreur relative:\t";
  //std::cout << abs(sum) << " %" << std::endl;
  std::cout << 100*sqrt( abs(sum - refsol)/abs(refsol) )<< " %" << std::endl;

}
