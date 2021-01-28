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

  int k=3;
  int K=4,L=k,W=2;
  int thresh=10,probes=2;
  double e=0.5;
  long long m=pow(2,32)-5,M=pow(2,32/K);



  string input_file;
  string query_file;
  string output_file;
  flags(&e,&thresh,&probes,&k , &input_file, &query_file, &output_file, argc,argv);
  L=k;

  cout<<"k="<<L<<" thresh="<<thresh<<" probes="<<probes<<" e="<<e<<endl;

  int dim=2;
  int k1=(-dim*log2(e))/(double)pow(e,2);

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,1.0);

  double G[k1][dim];
  for(int i=0;i<k1;i++){
    for(int j=0;j<dim;j++){
      G[i][j]=abs(distribution(generator));
    }
  }

  // W=get_W(input_file);
  // cout<<"W= "<<W<<endl;

  string data;
  ifstream myfile;
  myfile.open (input_file);

  int Data_num=0;

  vector<Curves> All_curvs;
  int flag;
  char chrs[4]="(),";
  int maxM=-1;

  int max_coord=0;


  cout<<"erxetai"<<endl;

  int counter2=0;
  while (getline(myfile,data)) {
    flag=0;
    istringstream ss(data);

    //Alfa a1(L,K);
    string id,length;
    ss >> id;
    ss >> length;
    int mlength=stoi(length);

    if(mlength>10){
      continue;
    }
    counter2++;
    if(maxM<mlength){
      maxM=mlength;
    }
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
          if(x>max_coord)
            max_coord=x;
          if(y>max_coord)
            max_coord=y;
          //printf("%.15f,%.15f\n",x,y );
          cur1.push_coordsi(d,x,y);
          d++;
        }
    }

    Data_num++;

    All_curvs.push_back(cur1);
  }

  traversals Biga[maxM][maxM];

  for(int i=0;i<maxM;i++){
    for(int j=0;j<maxM;j++){
      Biga[i][j].set_paths(i+1,j+1);
    }
  }

  int Count_b[maxM];
  for(int j=0;j<maxM;j++){
    Count_b[j]=0;
  }

  for(int i=0;i<All_curvs.size();i++){
    All_curvs[i].init_dian(maxM);
    int size=All_curvs[i].get_m()-1;
    Count_b[size]++;
    for(int j=0;j<maxM;j++){
      int path_s=Biga[size][j].get_size();
      for(int a=0;a<path_s-1;a++){
        Dianisma d(L,K);
        d.set_id(to_string(i));
        for(int b=0;b<Biga[size][j].get_path_size(a);b++){

          double result;
          for(int z=0;z<k1;z++){
            int x1=Biga[size][j].get_x(a,b);
            result=G[z][0]*All_curvs[i].get_coordsi(x1)[0]+G[z][1]*All_curvs[i].get_coordsi(x1)[1];
            d.push_coord(result);
          }
        }
        All_curvs[i].push_boxi_dianis(j,d);
      }
    }
  }

  vector<Hyper> All_hyper[maxM][maxM];
  int hash_categories=pow(2,L);

  for(int i=0;i<maxM;i++){
    for(int j=0;j<maxM;j++){
      int path_s=Biga[i][j].get_size();

      for(int z=0;z<path_s-1;z++){
        int d1=Biga[i][j].get_path_size(z);
        Hyper l(L,K,d1,W,hash_categories);
        All_hyper[i][j].push_back(l);
      }
    }
  }

  for(int i=0;i<All_curvs.size();i++){
    int size=All_curvs[i].get_m()-1;
    for(int j=0;j<maxM;j++){
      int path_s=Biga[size][j].get_size()-1;
      for(int a0=0;a0<path_s;a0++){
        Alfa a1(L,K);
        for(int z=0;z<All_curvs[i].get_diani(j,a0).get_coords().size();z++){
          int a_val;
          for(int o=0;o<L;o++){
            for(int l=0;l<K;l++){
              a_val=a(All_curvs[i].get_coordi_dianis(j,a0,z),All_hyper[size][j][a0].s[o][l][z],W);
              a1.push_a(o,l,a_val);
            }
          }

        }
        All_curvs[i].set_h_dianisi(j,a0, L, K, a1, m, M);

        All_curvs[i].set_gx_dianisi(j,a0,L,K);
        All_curvs[i].set_hashing_dianisi(&(All_hyper[size][j][a0].gxmap),j,a0,L,hash_categories);
        All_curvs[i].push_cell_dianis(&(All_hyper[size][j][a0].Hash),j,a0);
      }
    }
  }

  myfile.close();

  ofstream outfile;
  outfile.open(output_file);


  string path=query_file;

  while(path.compare("exit")!=0){
    double maxdiff=-1;
    double Average_time=0;

    myfile.open(path);

    int Query_num=0;
    vector<Curves> Query_curvs;
    counter2=0;

    while (getline(myfile,data)) {
      flag=0;

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

      Query_num++;
      Query_curvs.push_back(cur1);
    }

    myfile.close();

    for(int i=0;i<Query_curvs.size();i++){
      Query_curvs[i].init_dian(maxM);
      int size=Query_curvs[i].get_m()-1;
      for(int j=0;j<maxM;j++){
        int path_s=Biga[j][size].get_size();
        for(int a=0;a<path_s-1;a++){
          Dianisma d(L,K);
          d.set_id(to_string(i));
          for(int b=0;b<Biga[j][size].get_path_size(a);b++){

            double result;
            for(int z=0;z<k1;z++){
              int x1=Biga[j][size].get_y(a,b);
              result=G[z][0]*Query_curvs[i].get_coordsi(x1)[0]+G[z][1]*Query_curvs[i].get_coordsi(x1)[1];
              d.push_coord(result);
            }
          }
          Query_curvs[i].push_boxi_dianis(j,d);
        }
      }
    }


    for(int i=0;i<Query_curvs.size();i++){
      int size=Query_curvs[i].get_m()-1;
      for(int j=0;j<maxM;j++){
        int path_s=Biga[j][size].get_size()-1;
        for(int a0=0;a0<path_s;a0++){
          Alfa a1(L,K);
          for(int z=0;z<Query_curvs[i].get_diani(j,a0).get_coords().size();z++){
            int a_val;
            for(int o=0;o<L;o++){
              for(int l=0;l<K;l++){
                a_val=a(Query_curvs[i].get_coordi_dianis(j,a0,z),All_hyper[j][size][a0].s[o][l][z],W);
                a1.push_a(o,l,a_val);
              }
            }

          }
          Query_curvs[i].set_h_dianisi(j,a0, L, K, a1, m, M);
        }
      }
    }

    double total_sum=0.0;

    for(int i=0;i<Query_curvs.size();i++){
      double lsh_dist;
      double real_dist;
      vector<Dianisma> neighboor;
      vector<double> distances;

      int size=Query_curvs[i].get_m()-1;

      clock_t begin_lsh = clock();
      for(int o=0;o<maxM;o++){
        if(abs(o-size)>=4)
          continue;
        int path_s=Biga[o][size].get_size()-1;

        for(int a0=0;a0<path_s;a0++){

          Query_curvs[i].set_gx_dianisi(o,a0,L,K);
          Query_curvs[i].set_hashing_dianisi(&(All_hyper[o][size][a0].gxmap),o,a0,L,hash_categories);

          string fx=Query_curvs[i].get_diani(o,a0).get_fx();
          int cat=Query_curvs[i].get_diani(o,a0).get_cati(0);
          // cout<<"CATEGORY= "<<cat<<endl;
          int cau=0;
          for(int z=0;z<All_hyper[o][size][a0].Hash.get_cell(cat).size();z++){
              neighboor.push_back(All_hyper[o][size][a0].Hash.get_celli_dianj(cat,z));
              double apostasi=DTW(All_curvs[stoi(All_hyper[o][size][a0].Hash.get_celli_dianj(cat,z).get_id())],Query_curvs[i]);
              distances.push_back(apostasi);
              cau++;
              if(cau==thresh)
                break;
          }
          // cout<<"--"<<All_hyper[o][size][a0].Hash.get_cell(cat).size()<<"/"<<cau<<endl;


          vector<int> numbers;

          for(int y=0; y<L; y++)
              numbers.push_back(y);

          std::shuffle(numbers.begin(), numbers.end(),default_random_engine());

          if(cau!=thresh){
            int cau2=0;
            //int cau=0;
            for(int y=0;y<L;y++){
              int count1=0;
              unsigned int nei=cat^(1u<<(numbers[y]));
              for(int z=0;z<All_hyper[o][size][a0].Hash.get_cell(nei).size();z++){
                  neighboor.push_back(All_hyper[o][size][a0].Hash.get_celli_dianj(nei,z));
                  double apostasi=DTW(All_curvs[stoi(All_hyper[o][size][a0].Hash.get_celli_dianj(nei,z).get_id())],Query_curvs[i]);
                  distances.push_back(apostasi);
                  cau++;
                  count1++;
                  if(cau==thresh)
                    break;
              }
              // cout<<"--"<<All_hyper[o][size][a0].Hash.get_cell(nei).size()<<"/"<<count1<<endl;
              if(count1!=0)
                cau2++;
              if(cau2==probes)
                  break;
              if(cau==thresh)
                break;
            }
          }
        }
      }
      clock_t end_lsh = clock();
      int min_pos;
      if(!distances.empty()){
        vector<double>::iterator min_dis=min_element(distances.begin(), distances.end());
        min_pos=distance(distances.begin(), min_dis);

        lsh_dist=distances[min_pos];
        cout<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        outfile<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        cout<<"Hypercube_Neighbor: "<<All_curvs[stoi(neighboor[min_pos].get_id())].get_id()<<endl;
        outfile<<"Hypercube_Neighbor: "<<All_curvs[stoi(neighboor[min_pos].get_id())].get_id()<<endl;
        // cout<<"LSH_Distance: "<<distances[min_pos]<<endl;
        printf("Hypercube_Distance: %.15f\n",distances[min_pos] );
        outfile<<"Hypercube_Distance: "<<distances[min_pos]<<endl;
        cout<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
        outfile<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;

      }
      else{
        cout<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        outfile<<"Q_Item: "<<Query_curvs[i].get_id()<<endl;
        cout<<"Hypercube_Neighbor: ERROR_NO_ITEM_COULD_BE_FOUND"<<endl;
        outfile<<"Hypercube_Neighbor: ERROR_NO_ITEM_COULD_BE_FOUND"<<endl;
        cout<<"Hypercube_Distance: -"<<endl;
        outfile<<"Hypercube_Distance: -"<<endl;
        cout<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
        outfile<<"Hypercube_Time: "<<double(end_lsh - begin_lsh)/CLOCKS_PER_SEC<<endl;
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
      printf("TRUE_Distance: %.15f\n",min_d );
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
    cout<<"Average Time= "<<Average_time/(double)Query_num<<endl;
    cout<<"MAX LSH/REAL= "<<maxdiff<<endl;
    myfile.close();
    cout<<"MO= "<<total_sum/(double)Query_num<<endl;
    cout << "Enter the query name or exit if you want to leave" << '\n';
    cin>>path;
  }
  outfile.close();

  return 0;
}
