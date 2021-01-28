#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include "Classes.h"
#include "Functions.h"
#include <vector>
#include <math.h>
#include <ctime>
#include <unordered_map>


using namespace std;


int main(int argc,char *argv[])
{

  int k=3;
  int K=4,L=k,W=4000;//4656
  int thresh=10,probes=2;

  string input_file;
  string query_file;
  string output_file;
  flags(&thresh,&probes,&k , &input_file, &query_file, &output_file, argc,argv);
  L=k;
  long long m=pow(2,32)-5,M=pow(2,32/K);
  random_device rd; //Will be used to obtain a seed for the random number engine
  mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<double> dis(0.0,W);


  //W=get_W("./siftsmall/input_small_id");
  //cout<<"W= "<<W<<endl;

  string data;
  ifstream myfile;
  myfile.open (input_file);

  int Data_num=0;

  vector<Dianisma> All_v;
  int flag=0;

  vector<vector<vector<double> > > s = vector<vector<vector<double> > >( L, vector<vector<double> >(K, vector<double>(0,0)));

  cout<<"erxetai"<<endl;
  while (getline(myfile,data)) {

    istringstream ss(data);
    Dianisma vec1(L,K);
    Alfa a1(L,K);
    string id;
    ss >> id;
    vec1.set_id(id);
    int d=0;
    while (ss){
        string coor;
        ss >> coor;

        if(coor!=""){
          int a_val;
          for( int i=0;i<L;i++){
            for(int j=0;j<K;j++){
              if(flag==0){
                s[i][j].push_back(dis(gen));
              }
              a_val=a(stoi(coor),s[i][j][d],W);
              a1.push_a(i,j,a_val);
            }
          }

          vec1.push_coord(stoi(coor));
          d++;
        }
    }

    flag=1;
    h_calc(&vec1, L, K, a1, m, M);
    Data_num++;

    All_v.push_back(vec1);
  }

  myfile.close();

  int i,j,z;


  int hash_categories=pow(2,L);
  unordered_map<unsigned int,int> gxmap;
  unordered_map<unsigned int, int>:: iterator itr;



  Hash_table Hash(hash_categories);

  for (i = 0; i < All_v.size(); i++){
    gx_calc( &All_v[i], L, K);
    hashing(&gxmap,&All_v[i],L,hash_categories);

    Hash.push_cell(All_v[i].get_cati(0),&All_v[i]);
  }



  ofstream write_f;
  write_f.open(output_file);



  string path=query_file;

  cout<<"Start query"<<endl;
  while(path.compare("exit")!=0){

    myfile.open (path);
    int Query_num=0;
    vector<Dianisma> Query_v;

    //get the radius
    getline(myfile,data);
    istringstream ssR(data);
    string sR;
    ssR>>sR;
    ssR>>sR;
    double R;
    R=stod(sR);

    int dim;
    while (getline(myfile,data)) {
      istringstream ss(data);
      Dianisma vec1(L,K);
      Alfa a1(L,K);
      string id;
      ss >> id;
      vec1.set_id(id);
      dim=0;
      while (ss){
          string coor;
          ss >> coor;
          if(coor!=""){
            int a_val;
            for(int i=0;i<L;i++){
              for(int j=0;j<K;j++){
                a_val=a(stoi(coor),s[i][j][dim],W);
                a1.push_a(i,j,a_val);
              }
            }
            vec1.push_coord(stoi(coor));
            dim++;
          }
      }

      h_calc(&vec1, L, K, a1, m, M);
      Query_num++;
      Query_v.push_back(vec1);
    }

    double total_sum=0.0;

    for (i = 0; i < Query_v.size(); i++){
      int lsh_dist;
      int real_dist;

      vector<Dianisma> neighboor;
      vector<int> distances;

      gx_calc( &Query_v[i], L, K);
      hashing(&gxmap,&Query_v[i],L,hash_categories);

      clock_t begin_lsh = clock();
      int same_bucket=0,diff_bucket=0;
      string fx=Query_v[i].get_fx();
      int cat=Query_v[i].get_cati(0);
      //cout<<"CATEGORY: "<<cat<<endl;
      int cau=0;
      for(int z=0;z<Hash.get_cell(cat).size();z++){
          neighboor.push_back(Hash.get_celli_dianj(cat,z));
          distances.push_back(manhatan_dist(&(Query_v[i]),&(neighboor.back()),dim));
          cau++;
          if(cau==thresh)
            break;

      }
      //cout<<"--"<<Hash.get_cell(cat).size()<<"/"<<cau<<endl;
      vector<int> numbers;

      for(int j=0; j<L; j++)
          numbers.push_back(j);

      std::shuffle(numbers.begin(), numbers.end(),default_random_engine());

      if(cau!=thresh){
        int cau2=0;
        for(int j=0;j<L;j++){
          int count1=0;
          unsigned int nei=cat^(1u<<(numbers[j]));
          for(int z=0;z<Hash.get_cell(nei).size();z++){
              neighboor.push_back(Hash.get_celli_dianj(nei,z));
              distances.push_back(manhatan_dist(&(Query_v[i]),&(neighboor.back()),dim));
              cau++;
              count1++;
              if(cau==thresh)
                break;
          }
          //cout<<j<<"--"<<Hash.get_cell(nei).size()<<"/"<<count1<<endl;
          if(count1!=0)
            cau2++;
          if(cau2==probes)
              break;
          if(cau==thresh)
            break;

        }
      }

      clock_t end_lsh = clock();
      if(!distances.empty()){
        vector<int>::iterator min_dis=min_element(distances.begin(), distances.end());
        int min_pos=distance(distances.begin(), min_dis);

        lsh_dist=distances[min_pos];
        cout<<"Q_Item: "<<Query_v[i].get_id()<<endl;
        write_f<<"Q_Item: "<<Query_v[i].get_id()<<endl;
        cout<<"Hypercube_Neighbor: "<<neighboor[min_pos].get_id()<<endl;
        write_f<<"Hypercube_Neighbor: "<<neighboor[min_pos].get_id()<<endl;
        cout<<"Hypercube_Distance: "<<distances[min_pos]<<endl;
        write_f<<"Hypercube_Distance: "<<distances[min_pos]<<endl;
        cout<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
        write_f<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;

      }
      else{
        cout<<"Q_Item: "<<Query_v[i].get_id()<<endl;
        write_f<<"Q_Item: "<<Query_v[i].get_id()<<endl;
        cout<<"Hypercube_Neighbor: ERROR_NO_ITEM_COULD_BE_FOUND"<<endl;
        write_f<<"Hypercube_Neighbor: ERROR_NO_ITEM_COULD_BE_FOUND"<<endl;
        cout<<"Hypercube_Distance: -"<<endl;
        write_f<<"Hypercube_Distance: -"<<endl;
        cout<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
        write_f<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
      }

      clock_t begin_tr = clock();
      int min_d=manhatan_dist(&(Query_v[i]),&(All_v[0]),dim);
      string min_id=All_v[0].get_id();
      for(int j=1;j<All_v.size();j++){
        if(min_d>manhatan_dist(&(Query_v[i]),&(All_v[j]),dim)){
          min_d=manhatan_dist(&(Query_v[i]),&(All_v[j]),dim);
          min_id=All_v[j].get_id();
        }
      }
      clock_t end_tr = clock();

      real_dist=min_d;

      cout<<"TRUE_Time: "<<double(end_tr - begin_tr)/CLOCKS_PER_SEC<<endl;
      write_f<<"TRUE_Time: "<<double(end_tr - begin_tr)/CLOCKS_PER_SEC<<endl;
      cout<<"True_min_dist: "<<min_d<<endl;
      write_f<<"True_min_dist: "<<min_d<<endl;
      cout<<"True_Nei: "<<min_id<<endl;
      write_f<<"True_Nei: "<<min_id<<endl;
      cout<<"R-near neighbors:"<<endl;
      write_f<<"R-near neighbors:"<<endl;
      for(int j=0;j<distances.size();j++){
        if(distances[j]<R){
          cout<<"Item: "<<neighboor[j].get_id()<<"-Distance: "<<distances[j]<<endl;
        }
      }

      total_sum=total_sum+(lsh_dist/(double)real_dist);
      distances.clear();
      neighboor.clear();
    }
    myfile.close();
    cout<<"MO= "<<total_sum/(double)Query_v.size()<<endl;
    cout << "Enter the query name or exit if you want to leave" << '\n';
    cin>>path;
  }


  write_f.close();


  return 0;
}
