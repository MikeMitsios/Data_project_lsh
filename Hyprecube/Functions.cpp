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


void flags(int *thresh,int *probes,int *k,string *input_file,string *query_file,string *output_file,int argc,char **argv){
  string flags[6];

  flags[0] = string("-d");
  flags[1] = string("-q");
  flags[2] = string("-k");
  flags[3] = string("-o");
  flags[4] = string("-M");
  flags[5] = string("-probes");

  if(argc<=4){
		printf("WRONG INPUT\n");
		exit(1);
	}

  for(int j=1;j<argc;j=j+2){
    for(int u=0;u<6;u++){
      if(flags[u]==string(argv[j])){

        if(flags[u]=="-d"){
          *input_file=string(argv[j+1]);
        }
        else if(flags[u]=="-q"){
          *query_file=string(argv[j+1]);
        }
        else if(flags[u]=="-k"){
          *k=stoi(argv[j+1]);
        }
        else if(flags[u]=="-o"){
          *output_file=string(argv[j+1]);
        }
        else if(flags[u]=="-M"){
          *thresh=stoi(argv[j+1]);
        }
        else if(flags[u]=="-probes"){
          *probes=stoi(argv[j+1]);
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
    unsigned int gx_int=0;
    for(int z=0; z<K; z++){
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


void hashing(gx_map_pointer gxmap,Dianisma *v,int L,int hash_categories){


  random_device rd2;
  mt19937 e2(rd2());
  uniform_real_distribution<> dist(0, 1);

  int rand;
  string fx;

  for(int i=0;i<L;i++){

    int g=v->get_gxi(i);

    if(gxmap->find(g)==gxmap->end()){
      rand=round(dist(e2));
      gxmap->insert({g,rand});
    }
    else{
      rand=gxmap->at(g);
    }
    fx=fx+ to_string(rand);
  }
    v->set_fx(fx);
    unsigned int category;
    category=stoi(fx,nullptr,2);
    v->push_cat(category);

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
      if(All_a[i].get_id()!=All_a[j].get_id()){
        int manh=manhatan_dist(&All_a[i],&All_a[j],d);
        if(min>manh){
          min=manh;
        }
      }
    }
    distances.push_back(min);
  }

  int total_sum=accumulate(distances.begin(), distances.end(), 0)/Data_num;
std::cout << "TELOS" << '\n';

  myfile.close();

  return total_sum;

}
