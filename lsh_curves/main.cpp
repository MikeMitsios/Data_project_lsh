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


using namespace std;

int main(int argc,char *argv[])
{

  int K=4,L=1,W=50;//4656
  double delta=0.0016;
  int L_grid=4;
  int thresh=100*L;
  long long m=pow(2,32)-5,M=pow(2,32/K);

  string input_file;
  string query_file;
  string output_file;
  flags( &L_grid, &K, &input_file, &query_file, &output_file, argc,argv);


  //W=get_W("./siftsmall/input_small_id");
  //cout<<"W= "<<W<<endl;

  string data;
  ifstream myfile;
  myfile.open (input_file);

  int Data_num=0;

  vector<Curves> All_curvs;
  int flag;
  char chrs[4]="(),";
  int maxM=-1;

  int max_coord=0;

  vector<vector<double> > t=vector<vector<double> >(L_grid,vector<double>(0,0));

  random_device rd; //Will be used to obtain a seed for the random number engine
  mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<double> dis(0.0,2.0);

  for(int i=0;i<L_grid;i++){
    t[i].push_back(dis(gen));
    t[i].push_back(dis(gen));
  }

  cout<<"erxetai"<<endl;

  int counter2=0;
  while (getline(myfile,data)) {
    flag=0;

    counter2++;

    istringstream ss(data);

    //Alfa a1(L,K);
    string id,length;
    ss >> id;
    ss >> length;
    int mlength=stoi(length);


    if(maxM<mlength){
      maxM=mlength;
    }
    Curves cur1(L_grid,mlength,id);
    int d=0;    //min=0,
    vector<Dianisma> dianismata;
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
          if(x>max_coord)
            max_coord=x;
          if(y>max_coord)
            max_coord=y;
          //printf("%.15f,%.15f\n",x,y );
          cur1.push_coordsi(d,x,y);
          for(int i=0;i<L_grid;i++){
            double a1,a2;
            a1=(x-t[i][0])/(double)delta;
            a2=(y-t[i][1])/(double)delta;
            double x1=round(a1)*delta+t[i][0];
            double y1=round(a2)*delta+t[i][1];
            if(flag==0){
              Dianisma vec1(1,K);
              vec1.set_id(to_string(Data_num));
              vec1.set_mother(&cur1);
              vec1.push_coord(x1);
              vec1.push_coord(y1);
              cur1.push_dianis(vec1);

            }
            else{
              cur1.push_coord_dianis(i,x1);
              cur1.push_coord_dianis(i,y1);
            }
          }
          flag=1;
          d++;
        }
    }

    Data_num++;

    All_curvs.push_back(cur1);
  }

  int max_v=-1;
  for(int i=0; i<All_curvs.size();i++){
    for(int j=0;j<L_grid;j++){
        int count=0;
      for(int z=0;z<All_curvs[i].get_diani(j).get_coords().size()-2;z=z+2){
        if(All_curvs[i].get_coordi_dianis(j,z)==All_curvs[i].get_coordi_dianis(j,z+2) && All_curvs[i].get_coordi_dianis(j,z+1)==All_curvs[i].get_coordi_dianis(j,z+3)){
          count++;
          All_curvs[i].erase_diani_coordj(j,z+2);
          All_curvs[i].erase_diani_coordj(j,z+2);
          z=z-2;
        }

      }
      int vsize=All_curvs[i].get_diani(j).get_coords().size();
      if(max_v<vsize){
        max_v=vsize;
      }
    }
  }
  max_coord=max_coord*2;

  myfile.close();

  for(int i=0;i<All_curvs.size();i++){
    for(int j=0;j<L_grid;j++){
      int size=All_curvs[i].get_diani(j).get_coords().size();
      size=max_v-size;
      for(int z=0;z<size;z++){
        All_curvs[i].push_coord_dianis(j,max_coord);
      }
    }
  }


  random_device rd2; //Will be used to obtain a seed for the random number engine
  mt19937 gen2(rd2()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<double> dis2(0.0,W);

  vector<vector<vector<double> > > s = vector<vector<vector<double> > >( L_grid, vector<vector<double> >(K, vector<double>(0,0)));

  for(int i=0;i<L_grid;i++){
    for(int j=0;j<K;j++){
      for(int z=0;z<max_v;z++){
        s[i][j].push_back(dis2(gen2));
      }
    }
  }

  for(int i=0;i<All_curvs.size();i++){
    for(int j=0;j<L_grid;j++){
      Alfa a1(L,K);
      for(int z=0;z<max_v;z++){
        int a_val;
        for(int l=0;l<K;l++){
          a_val=a(All_curvs[i].get_coordi_dianis(j,z),s[j][l][z],W);
          a1.push_a(0,l,a_val);
        }
      }
      All_curvs[i].set_h_dianisi(j, L, K, a1, m, M);
    }
  }


  int hash_categories=Data_num/8;

  vector<Hash_table> Hash;

  for(int i=0;i<L_grid;i++){
    Hash_table h1(hash_categories);
    Hash.push_back(h1);
  }

  for(int i=0;i<All_curvs.size();i++){
    for(int j=0;j<L_grid;j++){
      All_curvs[i].set_gx_dianisi(j,L,K);
      All_curvs[i].set_hashing_dianisi(j,L,hash_categories);

      All_curvs[i].push_cell_dianis(&(Hash[j]),j);

    }
  }

  ofstream outfile;
  outfile.open(output_file);

  string path=query_file;

  while(path.compare("exit")!=0){

    double maxdiff=-1;
    double Average_time=0;

    myfile.open (path);
    int Query_num=0;
    vector<Curves> Query_curvs;
    counter2=0;

    while (getline(myfile,data)) {
      flag=0;
      // if(counter2==5){
      //   break;
      // }
      // counter2++;

      istringstream ss(data);

      string id,length;
      ss >> id;
      ss >> length;
      int mlength=stoi(length);

      Curves cur1(L_grid,mlength,id);
      int d=0;    //min=0,
      vector<Dianisma> dianismata;
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
            for(int i=0;i<L_grid;i++){
              double a1,a2;
              a1=(x-t[i][0])/(double)delta;
              a2=(y-t[i][1])/(double)delta;
              double x1=round(a1)*delta+t[i][0];
              double y1=round(a2)*delta+t[i][1];
              if(flag==0){
                Dianisma vec1(1,K);
                vec1.set_id(id);
                vec1.push_coord(x1);
                vec1.push_coord(y1);
                cur1.push_dianis(vec1);
              }
              else{
                cur1.push_coord_dianis(i,x1);
                cur1.push_coord_dianis(i,y1);
              }
            }
            flag=1;
            d++;
          }
      }
      Query_num++;
      Query_curvs.push_back(cur1);
    }

    for(int i=0; i<Query_curvs.size();i++){
      for(int j=0;j<L_grid;j++){
          int count=0;
          // cout<<All_curvs[i].get_diani(j).get_coords().size()<<endl;
        for(int z=0;z<Query_curvs[i].get_diani(j).get_coords().size()-2;z=z+2){
          if(Query_curvs[i].get_coordi_dianis(j,z)==Query_curvs[i].get_coordi_dianis(j,z+2) && Query_curvs[i].get_coordi_dianis(j,z+1)==Query_curvs[i].get_coordi_dianis(j,z+3)){
            count++;
            Query_curvs[i].erase_diani_coordj(j,z+2);
            Query_curvs[i].erase_diani_coordj(j,z+2);
            z=z-2;
          }
        }
      }
    }

    for(int i=0;i<Query_curvs.size();i++){
      for(int j=0;j<L_grid;j++){
        int size=Query_curvs[i].get_diani(j).get_coords().size();
        size=max_v-size;
        for(int z=0;z<size;z++){
          Query_curvs[i].push_coord_dianis(j,max_coord);
        }
      }
    }


    for(int i=0;i<Query_curvs.size();i++){
      for(int j=0;j<L_grid;j++){
        Alfa a1(L,K);
        for(int z=0;z<max_v;z++){
          int a_val;
          for(int l=0;l<K;l++){
            a_val=a(Query_curvs[i].get_coordi_dianis(j,z),s[j][l][z],W);
            a1.push_a(0,l,a_val);
          }
        }
        Query_curvs[i].set_h_dianisi(j, L, K, a1, m, M);
      }
    }

    double total_sum=0.0;

    for(int i=0;i<Query_curvs.size();i++){
      double lsh_dist;
      double real_dist;
      vector<Dianisma> neighboor;
      vector<double> distances;
      clock_t begin_lsh = clock();
      for(int j=0;j<L_grid;j++){
        Query_curvs[i].set_gx_dianisi(j,L,K);
        Query_curvs[i].set_hashing_dianisi(j,L,hash_categories);
        int gx=Query_curvs[i].get_diani(j).get_gxi(0);
        int cat=Query_curvs[i].get_diani(j).get_cati(0);
        //cout<<"CATEGORY= "<<cat<<endl;
        int cau=0;
        for(int z=0;z<Hash[j].get_cell(cat).size();z++){
          if(gx==Hash[j].get_celli_dianj(cat,z).get_gxi(0)){
            neighboor.push_back(Hash[j].get_celli_dianj(cat,z));
            double apostasi=DTW(All_curvs[stoi(Hash[j].get_celli_dianj(cat,z).get_id())],Query_curvs[i]);
            distances.push_back(apostasi);
            cau++;
            if(cau==thresh){
              break;
            }
          }
        }
        //cout<<"--"<<Hash[j].get_cell(cat).size()<<"/"<<cau<<endl;
      }

      clock_t end_lsh = clock();
      int min_pos;
      if(!distances.empty()){
        vector<double>::iterator min_dis=min_element(distances.begin(), distances.end());
        min_pos=distance(distances.begin(), min_dis);

        lsh_dist=distances[min_pos];
        cout<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        outfile<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        cout<<"LSH_Neighbor: "<<All_curvs[stoi(neighboor[min_pos].get_id())].get_id()<<endl;
        outfile<<"LSH_Neighbor: "<<All_curvs[stoi(neighboor[min_pos].get_id())].get_id()<<endl;
        // cout<<"LSH_Distance: "<<distances[min_pos]<<endl;
        printf("LSH_Distance: %.15f\n",distances[min_pos] );
        outfile<<"LSH_Distance: "<<distances[min_pos]<<endl;
        cout<<"LSH_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
        outfile<<"LSH_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;

      }
      else{
        cout<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        outfile<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        cout<<"LSH_Neighbor: ERROR_NO_ITEM_COULD_BE_FOUND"<<endl;
        outfile<<"LSH_Neighbor: ERROR_NO_ITEM_COULD_BE_FOUND"<<endl;
        cout<<"LSH_Distance: -"<<endl;
        outfile<<"LSH_Distance: -"<<endl;
        cout<<"LSH_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
        outfile<<"LSH_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
      }

      clock_t begin_tr = clock();
      double min_d=DTW(All_curvs[0],Query_curvs[i]);
      string min_id=All_curvs[0].get_id();
      for(int j=1;j<All_curvs.size();j++){
        double dist=DTW(All_curvs[j],Query_curvs[i]);
        if(min_d>dist){
          min_d=dist;
          min_id=All_curvs[j].get_id();
        }
      }
      clock_t end_tr = clock();

      real_dist=min_d;

      cout<<"TRUE_Time: "<<double(end_tr - begin_tr)/CLOCKS_PER_SEC<<endl;
      outfile<<"TRUE_Time: "<<double(end_tr - begin_tr)/CLOCKS_PER_SEC<<endl;
      cout<<"True_min_dist: "<<min_d<<endl;
      outfile<<"True_min_dist: "<<min_d<<endl;
      cout<<"True_Nei: "<<min_id<<endl;
      outfile<<"True_Nei: "<<min_id<<endl;

      if(!distances.empty() && min_d!=0){
        double temp1=distances[min_pos]/(double)min_d;
        if(maxdiff<temp1){
          maxdiff=temp1;
        }
      }
      Average_time=Average_time+double(end_lsh - begin_lsh)/CLOCKS_PER_SEC;

      if(real_dist!=0)
        total_sum=total_sum+(lsh_dist/(double)real_dist);

      distances.clear();
      neighboor.clear();
    }
    cout<<"Average Time= "<<Average_time<<"/"<<Query_num<<Average_time/(double)Query_num<<endl;
    cout<<"MAX LSH/REAL= "<<maxdiff<<endl;
    cout<<"MO= "<<total_sum/(double)Query_num<<endl;
    myfile.close();
    cout << "Enter the query name or exit if you want to leave" << '\n';
    cin>>path;
  }
  outfile.close();


  return 0;
}
