#ifndef ALGO_H
#define ALGO_H

#include <vector>
#include <queue>
#include "../calculus/calculus.hpp"
#include "mesh.hpp"

template <typename m_t>
class Adjacency{
  
public:
  static const int dim =       m_t::dim;
  typedef Adjacency<m_t>       this_t;
  typedef Elt<dim>             e_t;
  typedef Elt<dim-1>           f_t;
  typedef array<dim+1,f_t>     FaceArray;
  typedef array<dim+1,int>     NumArray;
  
private:
  const m_t&             mesh;
  const vector<e_t>&     elt;
  vector<f_t>            face;
  vector<NumArray>       neig;
  vector<NumArray>       back;
  
 public:
  // Constructeur
  Adjacency(const m_t&);
  
  // Acces aux donnees
  const NumArray&
  operator[](const int& j){
    return neig[j];}
  
  const NumArray&
  operator()(const int& j){
    return back[j];}
  
  friend const vector<f_t>&
  FacesOf(const this_t& adj){
    return adj.face;}
  
};




template <typename m_t>
Adjacency<m_t>::Adjacency(const m_t&  m):
  mesh(m), elt(EltOf(m)) {
  
  int nbnode   = NbNode(GeometryOf(m));
  int end      = -1; 
  
  vector<int>  first; 
  vector<int>  next;
  vector<N2>   num;  
  
  first.resize(nbnode,end);
  neig.resize(NbElt(m));
  back.resize(NbElt(m));  
  
  for(int j=0; j<NbElt(m); j++){
    bool exist;
    FaceArray aux = FacesOf(m[j]);    

    for(int k=0; k<dim+1; k++){
      
      f_t f = aux[k]; Order(f);
      exist = false;
      int& begin = first[Key(f)];
      
      //==================================//
      for(int q = begin; q!=end; q=next[q]){
	if(f==face[q]){
	  exist=true;
	  neig[j][k] = num[q][0];
	  back[j][k] = num[q][1];
	  neig[ num[q][0] ][ num[q][1] ] = j;
	  back[ num[q][0] ][ num[q][1] ] = k;
	  break;
	}
      }
      
      //==================================//
      if(!exist){
	face.push_back(f);
	next.push_back(begin);
	begin = face.size()-1;
	//N2 J; J[0] = j; J[1] = k;
	//num.push_back(J);
	num.push_back(N2_(j,k));
      }
      //==================================//      
      
    }
  }
}

typedef Adjacency<Mesh1D> Adjacency1D;
typedef Adjacency<Mesh2D> Adjacency2D;
typedef Adjacency<Mesh3D> Adjacency3D;



template <typename m_t>
class Connected{
  
 public:
  static const int dim =    m_t::dim;  
  typedef Connected<m_t>    this_t;
  typedef Elt<dim>          e_t;
  typedef Elt<dim-1>        f_t;

 private:
  // tableau  avec les no. des elements
  // de chaque composante.
  // num[j][k] est le no du k ieme elt du
  // la composante no. j
  vector< vector<int> >  num;

  // nbre de composantes
  int                    nbc; 

  // reference vers le maillage 
  const m_t&             mesh;
  
 public:
  Connected(const m_t&);  
  const vector<int>& operator[](const int& j)   {return num[j];}  
  friend const int&  NbOf      (const this_t& c){return c.nbc;}
  
};


//===============================//
//     Breadth First Search      //
//===============================//
template <typename m_t>
Connected<m_t>::Connected(const m_t& m): mesh(m) {  
  
  nbc = 1;
  num.resize(nbc);
  int nbelt  = NbElt(mesh);
  Adjacency<m_t> adj(mesh);
  int nb_visited = 0;
  queue<int> visit;
  vector<bool> visited(nbelt,false);
  
  // Initialisation de l'algo
  int j0 = 0;
  visit.push(j0);
  visited[j0]=true;
  num[nbc-1].push_back(j0);
    
  // Lancement de l'algo
  while(nb_visited<nbelt){
    
    // Reinitialisation dans le cas
    // de plusieurs composantes connexes
    if(visit.empty()){
      nbc++; num.resize(nbc);
      j0=0; while(visited[j0]){j0++;}
      visit.push(j0);
      visited[j0]=true;
      num[nbc-1].push_back(j0);      
    }
    else{
      j0 = visit.front();
      visit.pop();
    }
    nb_visited++;
    
    // Boucle sur les voisins de
    // l'element courant
    for(int k0=0; k0<dim+1; k0++){
      const int& j1 = adj[j0][k0];
      
      // Si voisin pas deja visite:
      // propagation de l'algo
      if(!visited[j1]){
	visited[j1]=true;	
	visit.push(j1);
	num[nbc-1].push_back(j1);      	
      }
    }
  }
  
  
}




#endif