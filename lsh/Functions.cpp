#include <iostream>
#include <random>
#include <cstring>
#include <fstream>
#include <bits/stdc++.h>
#include "Classes.h"
#include "Functions.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <unordered_map>


void flags(int *L,int *K,string *input_file,string *query_file,string *output_file,int argc,char **argv){
  string flags[5];

  flags[0] = string("-d");
  flags[1] = string("-q");
  flags[2] = string("-k");
  flags[3] = string("-L");
  flags[4] = string("-o");

  if(argc<=4){
		printf("WRONG INPUT\n");
		exit(1);
	}

  for(int j=1;j<argc;j=j+2){
    for(int u=0;u<5;u++){
      if(flags[u]==string(argv[j])){

        if(flags[u]=="-d"){
          *input_file=string(argv[j+1]);
        }
        else if(flags[u]=="-q"){
          *query_file=string(argv[j+1]);
        }
        else if(flags[u]=="-k"){
          *K=stoi(argv[j+1]);
        }
        else if(flags[u]=="-L"){
          *L=stoi(argv[j+1]);
        }
        else if(flags[u]=="-o"){
          *output_file=string(argv[j+1]);
        }
      }
    }
  }
}

int a(int x,double s,int W){
  return floor((x-s)/(double)W);
}


int h(int a,int m,int M){
  return modulo_(modulo_(a,M)*m,M);
}

int manhatan_dist(Dianisma *v1,Dianisma *v2,int dim){
  int manh_dist=0;
  for(int i=0;i<dim;i++){
    manh_dist=manh_dist+abs(v1->get_coordi(i)-v2->get_coordi(i));
  }
  return manh_dist;
}

void h_calc(Dianisma *v,int L,int K,Alfa a,long long m,int M){
  int j,z,d;

  for (j = 0; j < L; j++){
    for(int z=0; z<K; z++){
      int total_h=0;
      int Dim=a.get_a()[j][z].size();
      for(int d=0;d<Dim;d++){
        total_h=total_h+h(a.get_ai(j,z,d),MOD_power(m,Dim-1-d,M),M);
      }
      total_h=modulo_(total_h, M);
      v->push_h(j,total_h);
    }
  }
}


void gx_calc(Dianisma *v,int L,int K){
  for (int j = 0; j < L; j++){
    string gx;
    unsigned int gx_int=0;
    for(int z=0; z<K; z++){
      gx=gx+ to_string(v->get_hi(j,z));
      gx_int=gx_int|v->get_hi(j,z)<<(32/K)*(K-z-1);
    }
    v->push_gx(gx_int);
  }
}


int modulo_(int a, int b) {
  int m = a % b;
  if (m < 0) {
    m = (b < 0) ? m - b : m + b;
  }
  return m;
}


int MOD_power(long long x,  int y, int p) {
    int res = 1;      // Initialize result

    x = x % p;  // Update x if it is more than or
                // equal to p

    while (y > 0){
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res*x) % p;
        // y must be even now
        y = y>>1; // y = y/2
        x = (x*x) % p;
    }
    return res;
}


int Big_mod(long long int x, int size){
    int pos;
  	long long int y, mult;
  	y=x/(long long int)size;
  	mult=y*size;
  	pos=x-mult;

  	return pos;

}


void hashing(Dianisma *v,int L,int hash_categories){
  for(int i=0;i<L;i++){
    int g=v->get_gxi(i);

    int category=0;

    category=modulo_(g,hash_categories);
    v->push_cat(category);

  }
}




int get_W(string path){

  vector<Dianisma> All_a;
  vector<int> distances;

  string data;
  ifstream myfile;
  myfile.open (path);
  int d=0,Data_num=0,i,j;


  while (getline(myfile,data)){
    d=0;
    istringstream ss(data);
    string id;
    ss >> id;
    Dianisma vec(0,0);
    vec.set_id(id);
    while (ss){
      string coor;
      ss >> coor;
      if(coor!=""){

        vec.push_coord(stoi(coor));
        d++;
      }

    }
    All_a.push_back(vec);
    Data_num++;
  }

  for(int i=0;i<Data_num;i++){
    int min=manhatan_dist(&All_a[i],&All_a[0],d);

    for(int j=1;j<Data_num;j++){
      //cout<<"iID= "<<All_a.at(i).get_id()<<"jID= "<<All_a.at(j).get_id()<<endl;
      if(All_a[i].get_id()!=All_a[j].get_id()){
        int manh=manhatan_dist(&All_a[i],&All_a[j],d);
        //cout<<"MANHATAN= "<<manh<<endl;
        if(min>manh){
          min=manh;
        }
      }
    }
    distances.push_back(min);
  }

  int total_sum=accumulate(distances.begin(), distances.end(), 0)/Data_num;

  myfile.close();

  return total_sum;

}
