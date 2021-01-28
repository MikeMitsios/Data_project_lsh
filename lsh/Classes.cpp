#include <iostream>
#include "Classes.h"
#include "Functions.h"
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>


Dianisma::Dianisma(int L,int K){
  h=vector<vector<int> >(L,vector<int>(0,0));
}

void Dianisma::set_id(string id){
  this->id=id;
}
void Dianisma::set_gxi(int x,int i){
  this->gx[i]=x;
}
void Dianisma::set_coordi(int i,int val){
  this->coords[i]=val;
}
void Dianisma::set_hi(int i,int j,int val){
  this->h[i][j]=val;
}
void Dianisma::set_cati(int x,int i){
  this->categ[i]=x;
}


void Dianisma::push_coord(int coord){
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
vector<int> Dianisma::get_coords(){
  return this->coords;
}
vector<vector<int> > Dianisma::get_h(){
  return this->h;
}
vector<int> Dianisma::get_cat(){
  return this->categ;
}



int Dianisma::get_coordi(int i){
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
