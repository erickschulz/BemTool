//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_SUBOPERATOR_HPP
#define BEMTOOL_SUBOPERATOR_HPP

#include <vector>
#include <map>
#include <iostream>
#include "operator.hpp"

namespace bemtool {

  class BlockMat{

  private:
    int           nr,nc;
    std::vector<Cplx> v;

  public:
    BlockMat(const int& nr0 = 0, const int& nc0 = 0): nr(nr0), nc(nc0), v(nr*nc,0.){}
    Cplx& operator()(const int& j, const int& k){return v[j*nr+k];}
    const Cplx& operator()(const int& j, const int& k) const {return v[j*nr+k];}
    template <typename r_t> void operator=(const r_t& r){
      for(int j=0; j<nr; j++){for(int k=0; k<nc; k++){v[j*nr+k] = r(j,k);}}}
    void Clear(){for(int j=0;j<v.size();j++){v[j]=0.;}}
    void Resize(const int& newnr, const int& newnc){
      Clear(); nr = newnr; nc = newnc; v.resize(nr*nc);}
    friend const int& NbRow(const BlockMat& m){return m.nr;}
    friend const int& NbCol(const BlockMat& m){return m.nc;}
    friend std::ostream& operator<<(std::ostream& o, const BlockMat& m){
      for(int j=0; j<NbRow(m);j++){for(int k=0; k<NbCol(m);k++){
	  o << m(j,k) << "\t";} o << "\n";}
      return o;}

  };


  template <int d> class NDofLoc{

    typedef NDofLoc<d> this_t;
    static const int dim = d;

  private:
    int v[d];

  public:
    NDofLoc<d>(){
      for(int j=0;j<d;j++){v[j]=-1;}}
    NDofLoc<d>(const this_t& I){
      for(int j=0;j<d;j++){v[j]=I[j];}}
    void operator=(const this_t& I){
      for(int j=0;j<d;j++){v[j]=I[j];}}
    int& operator[](const int& j){return v[j];}
    const int& operator[](const int& j) const {return v[j];}
    friend std::ostream& operator<<(std::ostream& o, const this_t& I){
      for(int j=0; j<this_t::dim; j++){ o<<I[j]<<"\t";} return o;}
  };


  template<typename BIOpType>
  class SubBIOp{

    typedef typename BIOpType::KernelTypeTrait  KernelTrait;
    typedef typename KernelTrait::ShapeFctX     PhiX;
    typedef typename KernelTrait::ShapeFctY     PhiY;
    typedef typename KernelTrait::MatType       MatType;
    typedef Dof<PhiX>                           DofX;
    typedef Dof<PhiY>                           DofY;

    static const int nb_dof_loc_x = KernelTrait::nb_dof_x;
    static const int nb_dof_loc_y = KernelTrait::nb_dof_y;
    typedef NDofLoc<nb_dof_loc_x> Nlocx;
    typedef NDofLoc<nb_dof_loc_y> Nlocy;
    typedef typename std::map<int,Nlocx>::iterator ItTypeX;
    typedef typename std::map<int,Nlocy>::iterator ItTypeY;


  private:
    BIOpType      biop;
    const DofX&    dofx;
    const DofY&    dofy;

    BlockMat             block_mat;
    MatType              elt_mat;
    std::map<int,Nlocx>  Ix;
    std::map<int,Nlocy>  Iy;


  public:
    SubBIOp<BIOpType>(BIOpType& biop0,
		      const DofX&     dofx0,
		      const DofY&     dofy0): biop(biop0), dofx(dofx0), dofy(dofy0){}

    const BlockMat& operator()(const std::vector<int>& jjx, const std::vector<int>& jjy){
      block_mat.Resize(jjx.size(),jjy.size());
      elt_mat = 0.; Ix.clear(); Iy.clear();

      for(int k=0; k<jjx.size(); k++){
	const std::vector<N2>& jj = dofx.ToElt(jjx[k]);
	for(int l=0; l<jj.size(); l++){
	  const N2& j = jj[l]; Ix[j[0]][j[1]] = k;
	}
      }

      for(int k=0; k<jjy.size(); k++){
	const std::vector<N2>& jj = dofy.ToElt(jjy[k]);
	for(int l=0; l<jj.size(); l++){
	  const N2& j = jj[l]; Iy[j[0]][j[1]] = k;
	}
      }

      for(ItTypeX itx = Ix.begin(); itx!=Ix.end(); itx++){
	const int&   jx = itx->first;
	const Nlocx& nx = itx->second;

	for(ItTypeY ity = Iy.begin(); ity!=Iy.end(); ity++){
	  const int&   jy = ity->first;
	  const Nlocy& ny = ity->second;

	  elt_mat = biop(jx,jy);
	  for(int kx=0; kx<nb_dof_loc_x; kx++){ if(nx[kx]!=-1){
	      for(int ky=0; ky<nb_dof_loc_y; ky++){ if(ny[ky]!=-1){
		  block_mat(nx[kx],ny[ky]) += elt_mat(kx,ky);
		}}
	    }}

	}
      }

      return block_mat;
    }

    template<typename Matrix>
    void compute_block(const std::vector<int>& jjx, const std::vector<int>& jjy, Matrix& mat){
      elt_mat = 0.; Ix.clear(); Iy.clear();

      for(int k=0; k<jjx.size(); k++){
    const std::vector<N2>& jj = dofx.ToElt(jjx[k]);
    for(int l=0; l<jj.size(); l++){
      const N2& j = jj[l]; Ix[j[0]][j[1]] = k;
    }
      }

      for(int k=0; k<jjy.size(); k++){
    const std::vector<N2>& jj = dofy.ToElt(jjy[k]);
    for(int l=0; l<jj.size(); l++){
      const N2& j = jj[l]; Iy[j[0]][j[1]] = k;
    }
      }

      for(ItTypeX itx = Ix.begin(); itx!=Ix.end(); itx++){
    const int&   jx = itx->first;
    const Nlocx& nx = itx->second;

    for(ItTypeY ity = Iy.begin(); ity!=Iy.end(); ity++){
      const int&   jy = ity->first;
      const Nlocy& ny = ity->second;

      elt_mat = biop(jx,jy);
      for(int kx=0; kx<nb_dof_loc_x; kx++){ if(nx[kx]!=-1){
          for(int ky=0; ky<nb_dof_loc_y; ky++){ if(ny[ky]!=-1){
          mat(nx[kx],ny[ky]) += elt_mat(kx,ky);
        }}
        }}

    }
      }

    }

    template<typename Matrix>
    void compute_neumann_block(const std::vector<int>& jjx, const std::vector<int>& jjy, Matrix& mat){
      elt_mat = 0.; Ix.clear(); Iy.clear();

      for(int k=0; k<jjx.size(); k++){
	const std::vector<N2>& jj = dofx.ToElt(jjx[k]);
	for(int l=0; l<jj.size(); l++){
	  const N2& j = jj[l]; Ix[j[0]][j[1]] = k;
	}
      }

      for(int k=0; k<jjy.size(); k++){
	const std::vector<N2>& jj = dofy.ToElt(jjy[k]);
	for(int l=0; l<jj.size(); l++){
	  const N2& j = jj[l]; Iy[j[0]][j[1]] = k;
	}
      }

      for(ItTypeX itx = Ix.begin(); itx!=Ix.end(); itx++){
	const int&   jx = itx->first;
	const Nlocx& nx = itx->second;

	for(ItTypeY ity = Iy.begin(); ity!=Iy.end(); ity++){
	  const int&   jy = ity->first;
	  const Nlocy& ny = ity->second;

      bool testx=1;
      bool testy=1;
      for(int kx=0; kx<nb_dof_loc_x; kx++){
          if(nx[kx]==-1){
                testx=0;
          }
      }
      for(int ky=0; ky<nb_dof_loc_y; ky++){
          if(ny[ky]==-1){
                testy=0;
          }
      }
      if (testx && testy){
	  elt_mat = biop(jx,jy);
	  for(int kx=0; kx<nb_dof_loc_x; kx++){
	      for(int ky=0; ky<nb_dof_loc_y; ky++){
		  mat(nx[kx],ny[ky]) += elt_mat(kx,ky);
		}
	    }
    }
	}
      }

    }


  };


}
#endif