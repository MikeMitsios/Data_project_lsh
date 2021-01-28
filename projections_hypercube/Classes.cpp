#include <iostream>
#include "Classes.h"
#include "Functions.h"
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>


Dianisma::Dianisma(int L,int K){
  //a=vector<vector<vector<int> > >(L,vector<vector<int> >(K ,vector<int>(0,0)));
  h=vector<vector<int> >(L,vector<int>(0,0));
}

void Dianisma::set_id(string id){
  this->id=id;
}
void Dianisma::set_gxi(int x,int i){
  this->gx[i]=x;
}
void Dianisma::set_coordi(int i,double val){
  this->coords[i]=val;
}
void Dianisma::set_hi(int i,int j,int val){
  this->h[i][j]=val;
}
void Dianisma::set_cati(int x,int i){
  this->categ[i]=x;
}
void Dianisma::set_mother(Curves *x){
  this->mother=x;
}
void Dianisma::set_fx(string x){
  this->fx=x;
}


void Dianisma::push_coord(double coord){
  this->coords.push_back(coord);
}
void Dianisma::push_h(int i,int h){
  this->h[i].push_back(h);
}
void Dianisma::push_gx(int gx){
  this->gx.push_back(gx);
}
void Dianisma::push_cat(int cat){
  this->categ.push_back(cat);
}



vector<int> Dianisma::get_gx(){
  return this->gx;
}
string Dianisma::get_id(){
  return this->id;
}
vector<double> Dianisma::get_coords(){
  return this->coords;
}
vector<vector<int> > Dianisma::get_h(){
  return this->h;
}
vector<int> Dianisma::get_cat(){
  return this->categ;
}
Curves * Dianisma::get_mother(){
  return this->mother;
}
string Dianisma::get_fx(){
  return this->fx;
}



double Dianisma::get_coordi(int i){
  return this->coords[i];
}
int Dianisma::get_hi(int i,int j){
  return this->h[i][j];
}
int Dianisma::get_gxi(int i){
  return this->gx[i];
}
int Dianisma::get_cati(int i){
  return this->categ[i];
}

void Dianisma::erase_coordi(int i){
  this->coords.erase(this->coords.begin() + i);
}



Alfa::Alfa(int L,int K){
  a=vector<vector<vector<int> > >(L,vector<vector<int> >(K ,vector<int>(0,0)));
}

void Alfa::set_ai(int i,int j,int z,int val){
  this->a[i][j][z]=val;
}
void Alfa::push_a(int i,int j,int a){
  this->a[i][j].push_back(a);
}
vector< vector<vector<int> > > Alfa::get_a(){
  return this->a;
}
int Alfa::get_ai(int i,int j,int z){
  return this->a[i][j][z];
}



Cell::Cell(){
  this->dianismata=vector<Dianisma*>(0,0);
}
void Cell::push_dian(Dianisma *v){
  this->dianismata.push_back(v);
}
Dianisma Cell::get_diani(int i){
  return *(this->dianismata[i]);
}
vector<Dianisma*> Cell::get_dianismata(){
  return this->dianismata;
}


Hash_table::Hash_table(){

}
Hash_table::Hash_table(int hash_categories){
  for(int i=0;i<hash_categories;i++){
    Cell c1;
    this->cells.push_back(c1);
  }

}
void Hash_table::push_cell(int i,Dianisma *v){
  this->cells[i].push_dian(v);
}
Dianisma Hash_table::get_celli_dianj(int i,int j){
  return this->cells[i].get_diani(j);
}
vector<Dianisma*> Hash_table::get_cell(int i){
  return this->cells[i].get_dianismata();
}




Curves::Curves(int m,string id){
  this->m=m;
  this->id=id;
  this->coords=vector<vector<double> >(m,vector<double>(0,0));

}
void Curves::init_dian(int maxM){
  for(int i=0;i<maxM;i++){
    vector<Dianisma> d;
    this->dianismata.push_back(d);
  }
}
void Curves::set_h_dianisi(int i,int j,int  L,int  K,Alfa a1,long long m,int M){
  h_calc(&(this->dianismata[i][j]), L, K, a1, m, M);
}
void Curves::set_gx_dianisi(int i,int j,int  L,int  K){
  gx_calc(&(this->dianismata[i][j]), L, K);
}
void Curves::set_hashing_dianisi(gx_map_pointer gxmap,int i,int j,int  L,int  hash_categories){
  hashing(gxmap,&(this->dianismata[i][j]), L, hash_categories);
}
void Curves::push_cell_dianis(Hash_table *Hash,int i,int j){
  Hash->push_cell(this->dianismata[i][j].get_cat()[0],&(this->dianismata[i][j]));
}
void Curves::push_coordsi(int i,double x,double y){
  this->coords[i].push_back(x);
  this->coords[i].push_back(y);
}
void Curves::push_boxi_dianis(int i,Dianisma d){
  this->dianismata[i].push_back(d);
}

string Curves::get_id(){
  return this->id;
}
int Curves::get_m(){
  return this->m;
}
double Curves::get_coordi_dianis(int i,int j,int z){
  return this->dianismata[i][j].get_coordi(z);
}
vector<vector<double> > Curves::get_coords(){
  return this->coords;
}
Dianisma Curves::get_diani(int i,int j){
  return this->dianismata[i][j];
}
vector<double> Curves::get_coordsi(int i){
  return this->coords[i];
}


void traversals::set_paths(int m,int n){
  this->paths=vector<vector<point>>(1,vector<point>(0));
  find_paths(m,n,&(this->paths));
  this->size=this->paths.size();
}
vector<vector<point>> traversals::get_paths(){
  return this->paths;
}

int traversals::get_x(int i,int j){
  return this->paths[i][j].x;
}
int traversals::get_y(int i,int j){
  return this->paths[i][j].y;
}
int traversals::get_size(){
  return this->size;
}
int traversals::get_path_size(int i){
  return this->paths[i].size();
}


Hyper::Hyper(int L,int K,int d,int W,int hash_categories){
  random_device rd2; //Will be used to obtain a seed for the random number engine
  mt19937 gen2(rd2()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<double> dis2(0.0,W);

  Hash_table h1(hash_categories);
  Hash=h1;

  s=vector<vector<vector<double> > >( L, vector<vector<double> >(K, vector<double>(0,0)));
  for( int i=0;i<L;i++){
    for(int j=0;j<K;j++){
      for(int z=0;z<d;z++){
        s[i][j].push_back(dis2(gen2));
      }

    }
  }

}
