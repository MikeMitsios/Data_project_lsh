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
#include <limits>


void flags(double *e,int *L,int *K,string *input_file,string *query_file,string *output_file,int argc,char **argv){
  string flags[6];

  flags[0] = string("-d");
  flags[1] = string("-q");
  flags[2] = string("-k_vec");
  flags[3] = string("-L_vec");
  flags[4] = string("-o");
  flags[5] = string("-e");

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
        else if(flags[u]=="-k_vec"){
          *K=stoi(argv[j+1]);
        }
        else if(flags[u]=="-L_vec"){
          *L=stoi(argv[j+1]);
        }
        else if(flags[u]=="-o"){
          *output_file=string(argv[j+1]);
        }
        else if(flags[u]=="-e"){
          *e=stod(argv[j+1]);
        }
      }
    }
  }
}

int a(double x,double s,int W){
  return floor((x-s)/(double)W);
}

void removeCharsFromString( string &str, char* charsToRemove ) {
   for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
      str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   }
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

double Euclidean_dist(vector<double> x,vector<double> y){
  double e_dist=0.0;
  e_dist=pow(x[0]-y[0],2)+pow(x[1]-y[1],2);
  e_dist=sqrt(e_dist);
  return e_dist;
}

double minimum(double x1,double x2, double x3){
  double min=x1;

  if(min>x2){
    min=x2;
  }
  if(min>x3){
    min=x3;
  }
  return min;
}


double DTW(Curves& curv1,Curves& curv2){

  int m1=curv1.get_m(),m2=curv2.get_m();
  double DTW_array[m1][m2];
  DTW_array[0][0]=Euclidean_dist(curv1.get_coordsi(0),curv2.get_coordsi(0));
  int i,j;
  for(i=1;i<m2;i++){
    DTW_array[0][i]=DTW_array[0][i-1]+Euclidean_dist(curv1.get_coordsi(0),curv2.get_coordsi(i));
  }
  for(i=1;i<m1;i++){
    DTW_array[i][0]=DTW_array[i-1][0]+Euclidean_dist(curv1.get_coordsi(i),curv2.get_coordsi(0));
  }
  for(i=1;i<m1;i++){
    for(j=1;j<m2;j++){
      DTW_array[i][j]=minimum(DTW_array[i][j-1],DTW_array[i-1][j-1],DTW_array[i-1][j])+Euclidean_dist(curv1.get_coordsi(i),curv2.get_coordsi(j));
    }
  }

  return DTW_array[m1-1][m2-1];
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
    unsigned int g=v->get_gxi(i);

    int category=0;

    category=modulo_(g,hash_categories);

    v->push_cat(category);

  }
}

int big_count=0;

void find_all_paths(int *mat, int i, int j, int m, int n, point *path, int pi,vector<vector<point>> *p)
{
    double lg,tlg;
    if(m%n==0){
      if(m>n){
        lg=m/(double)n;
        if(j+1!=ceil((i+1)/(double)lg)){
          (*p)[big_count].clear();
          return;
        }

      }
      else if(n>m){
        lg=n/(double)m;
        if(i+1!=ceil((j+1)/(double)lg)){
          (*p)[big_count].clear();
          return;
        }
      }
      else{
        lg=n/(double)m;
        tlg=(j+1)/(double)(i+1);
        if(lg-tlg!=0){
          (*p)[big_count].clear();
          return;
        }
      }
    }
    else{
      if(m>n){
        lg=m/(double)n;
        if(j+1!=round((i+1)/(double)lg)){
          if(j+1!=ceil((i+1)/(double)lg) ){
            (*p)[big_count].clear();
            return;
          }

        }

      }
      else if(n>m){
        lg=n/(double)m;
        if(i+1!=round((j+1)/(double)lg)){
          if(i+1!=ceil((j+1)/(double)lg) ){
            (*p)[big_count].clear();
            return;
          }
        }
      }
      else{
        lg=n/(double)m;
        tlg=(j+1)/(double)(i+1);
        if(lg-tlg!=0){
          (*p)[big_count].clear();
          return;
        }
      }
    }

    // if(abs(i-j)>=2){
    //   (*p)[big_count].clear();
    //   return;
    // }
    // Reached the bottom of the matrix so we are left with
    // only option to move right
    if (i == m - 1)
    {
        for (int k = j; k < n; k++){
          path[pi + k - j].x = i; //= *((mat + k*n) + j);
          path[pi + k - j].y = k;
        }
        for (int l = 0; l < pi + n - j; l++){
            (*p)[big_count].push_back(path[l]);
        }
        std::vector<point> v;
        (*p).push_back(v);
        big_count++;
        return;
    }

    // Reached the right corner of the matrix we are left with
    // only the downward movement.
    if (j == n - 1)
    {
        for (int k = i; k < m; k++){
            path[pi + k - i].x = k; //= *((mat + k*n) + j);
            path[pi + k - i].y = j;
        }
        for (int l = 0; l < pi + m - i; l++){
            (*p)[big_count].push_back(path[l]);

        }
        std::vector<point> v;
        (*p).push_back(v);
        big_count++;
        return;
    }

    // Add the current cell to the path being generated
    path[pi].x = i;
    path[pi].y = j;
    // Print all the paths that are possible after moving down
    find_all_paths(mat, i+1, j, m, n, path, pi + 1,p);

    // Print all the paths that are possible after moving right
    find_all_paths(mat, i, j+1, m, n, path, pi + 1,p);

    find_all_paths(mat, i+1, j+1, m, n, path, pi + 1,p);

}

// The main function that prints all paths from top left to bottom right
// in a matrix 'mat' of size mXn
void find_paths(int m, int n,vector<vector<point>> *p)
{
    int mat[m][n];
    for(int i=0;i<m;i++){
      for(int j=0;j<n;j++){
        mat[i][j] = 0;
      }
    }

    point *path = new point[m+n];
    big_count=0;
    find_all_paths(*mat, 0, 0, m, n, path, 0, p);
}




int get_W(string path){

  string data;
  ifstream myfile;
  myfile.open (path);

  int Data_num=0;

  vector<Curves> All_c;
  char chrs[4]="(),";

  int counter2=0;
  while (getline(myfile,data)) {
    // if(counter2==100){
    //   break;
    // }

    istringstream ss(data);
    string id,length;
    ss >> id;
    ss >> length;
    int mlength=stoi(length);
    if(mlength>10){
      continue;
    }
    counter2++;
    Curves cur1(mlength,id);
    int d=0;    //min=0,
    while (ss){
        string coor1;
        ss >> coor1;
        removeCharsFromString(coor1,chrs);
        string coor2;
        ss >> coor2;
        removeCharsFromString(coor2,chrs);

        if(coor1!=""){
          double x,y;
          x=stod(coor1);
          y=stod(coor2);
          cur1.push_coordsi(d,x,y);
          d++;
        }
    }
    Data_num++;
    All_c.push_back(cur1);
  }
  double total_sum=0;
  for(int i=0;i<Data_num-1;i++){
    total_sum=total_sum+DTW(All_c[i],All_c[i+1]);
  }

  myfile.close();

  return (total_sum*4)/(double)Data_num;

}
