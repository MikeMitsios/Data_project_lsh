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




Curves::Curves(int L_grid,int m,string id){
  this->m=m;
  this->id=id;
  this->coords=vector<vector<double> >(m,vector<double>(0,0));
}
void Curves::set_h_dianisi(int i,int  L,int  K,Alfa a1,long long m,int M){
  h_calc(&(this->dianismata[i]), L, K, a1, m, M);
}
void Curves::set_gx_dianisi(int i,int  L,int  K){
  gx_calc(&(this->dianismata[i]), L, K);
}
void Curves::set_hashing_dianisi(gx_map_pointer gxmap,int i,int  L,int  hash_categories){
  hashing(gxmap,&(this->dianismata[i]), L, hash_categories);
}
void Curves::push_cell_dianis(Hash_table *Hash,int j){
  Hash->push_cell(this->dianismata[j].get_cat()[0],&(this->dianismata[j]));
}
void Curves::push_coordsi(int i,double x,double y){
  this->coords[i].push_back(x);
  this->coords[i].push_back(y);
}
void Curves::push_dianis(Dianisma x){
  this->dianismata.push_back(x);
}
void Curves::push_coord_dianis(int i,double coord){
  this->dianismata[i].push_coord(coord);
}
string Curves::get_id(){
  return this->id;
}
int Curves::get_m(){
  return this->m;
}
double Curves::get_coordi_dianis(int i,int j){
  return this->dianismata[i].get_coordi(j);
}
vector<vector<double> > Curves::get_coords(){
  return this->coords;
}
vector<double> Curves::get_coordsi(int i){
  return this->coords[i];
}
Dianisma Curves::get_diani(int i){
  return this->dianismata[i];
}

void Curves::erase_diani_coordj(int i,int j){
  this->dianismata[i].erase_coordi(j);
}
