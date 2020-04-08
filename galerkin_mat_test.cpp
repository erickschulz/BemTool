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

  Real kappa = 0; // Is this parameter necessary for the Laplace equation?
  Real kappa2 = kappa*kappa;

  // Loading the nodes into the geometry object? This is how to read .msh files
  Geometry node("mesh/circle05.msh");
  // Mesh1D object mesh? loading the node into the mesh?
  Mesh1D mesh; mesh.Load(node,1);
  // Orienting the mesh
  Orienting(mesh);
  // Number of elements in the mesh
  int nb_elt = NbElt(mesh);
  std::cout << "nb_elt:\t" << nb_elt << std::endl;

  // Defining the operators, this is the right way! (test,trial)
  typedef LA_SL_2D_P0xP0 Voperator;
  typedef LA_DL_2D_P0xP1 Koperator;

  using Mass_2D_P0xP1 = LocalMatrix<P0_1D,P1_1D>;

  //typedef CST_1D_P0xP1 Moperator; // CST is for double integrals!

  // Defining the BIOp objects based on the operators defined above
  BIOp<Voperator> VV(mesh,mesh,kappa);
  BIOp<Koperator> KK(mesh,mesh,kappa);
  Mass_2D_P0xP1 MM(mesh);

  // Dof objects for P1 and P0 (provided as template arguments)
  Dof<P1_1D> dof1(mesh);
  Dof<P0_1D> dof0(mesh);

  // Checking degrees of freedom
  int nb_dof1 = NbDof(dof1);
  int nb_dof0 = NbDof(dof0);
  std::cout << "dof1 : " << nb_dof1 << std::endl;
  std::cout << "dof0 : " << nb_dof0 << std::endl;

  EigenDense V(nb_dof0,nb_dof0); Clear(V);
  EigenDense K(nb_dof0,nb_dof1); Clear(K); // test, trial
  EigenDense M(nb_dof0,nb_dof1); Clear(M);

  // Non-panel oriented assembly ---------->

  // Initializing the Galerkin matrices with the right sizes. Matrix type = Complex?
  /*
  // Assembly of the Galerkin matrices
  progress bar("Assemblage V",nb_dof0);
  for(int j=0; j<nb_dof0; j++){
    bar++;
    for(int k=0; k<nb_dof0; k++){
      V(j,k) += VV(dof0.ToElt(j),dof0.ToElt(k)); // dof to elt represents global to local map?
    }
  }
  bar.end();

  progress bar1("Assemblage K",nb_dof0);
  for(int j=0; j<nb_dof0; j++){
    bar++;
    for(int k=0; k<nb_dof1; k++){
        K(j,k) += KK(dof0.ToElt(j),dof1.ToElt(k));
  }
}
  bar1.end();

  progress bar2("Assemblage M",nb_dof0);
  for(int j=0; j<nb_dof0; j++){
    bar++;
    for(int k=0; k<nb_dof1; k++){
      M(j,k) += MM(dof0.ToElt(j),dof1.ToElt(k));
    }
  }
  bar2.end();*/

  // Panel-oriented assembly -------->
  progress bar("Assemblage V",nb_elt);
  for(int j=0; j<nb_elt; j++){
    bar++;
    for(int k=0; k<nb_elt; k++){
      V(dof0[j],dof0[k]) += VV(j,k);
    }
  }

  progress bar1("Assemblage K",nb_elt);
  for(int j=0; j<nb_elt; j++){
    bar1++;
    for(int k=0; k<nb_elt; k++){
      K(dof0[j],dof1[k]) += KK(j,k);
    }
  }

  progress bar2("Assemblage M",nb_elt);
  for(int j=0; j<nb_elt; j++){
    bar2++;
      M(dof0[j],dof1[j]) += MM(j);
  }

  // Vector encoding the Dirichlet boundary condition
  Eigen::VectorXd g = Eigen::VectorXd::Constant(nb_dof1,5);

  // Galerkin matrices of double type.
  Eigen::MatrixXd V1(nb_dof0,nb_dof0);
  Eigen::MatrixXd K1(nb_dof0,nb_dof1);
  Eigen::MatrixXd M1(nb_dof0,nb_dof1);

  // Extracting the real parts of the Galerkin matrices to solve the linear system
  for (int i = 0 ; i < nb_dof0 ; ++i) {
    for (int j = 0 ; j < nb_dof0 ; ++j) {
      V1(i,j) = V(i,j).real();
    }
  }

  for (int i = 0 ; i < nb_dof0 ; ++i) {
    for (int j = 0 ; j < nb_dof1 ; ++j) {
      K1(i,j) = K(i,j).real();
      M1(i,j) = M(i,j).real();
    }
  }

  // Solving the linear system given by direct first kind formulation
  Eigen::VectorXd Neumann_trace = V1.lu().solve((0.5*M1+K1)*g);
  std::cout << "Neumann trace evaluated \n" << Neumann_trace.transpose() << std::endl;
  std::cout << "RHS (expected to be 0 vector): \n" << ((0.5*M1+K1)*g).transpose() << std::endl;

}
