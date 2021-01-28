#include <vector>
#include <iostream>


using namespace std;


class Cell;
class Dianisma;
class Curves;


class Hash_table{
private:
  vector<Cell> cells;
public:
  Hash_table(int hash_categories);
  void push_cell(int i,Dianisma *v);
  Dianisma get_celli_dianj(int i,int j);
  vector<Dianisma*> get_cell(int i);
};

class Cell{
private:
  vector<Dianisma*> dianismata;
public:
  Cell();
  void push_dian(Dianisma *v);
  Dianisma get_diani(int i);
  vector<Dianisma*> get_dianismata();
};

class Dianisma{

private:
  vector<double> coords;
  string id;
  vector<vector<int> >h;
  vector<int> gx;
  vector<int> categ;
  Curves *mother;

public:
  Dianisma(int L,int K);

  void set_id(string id);
  void set_coordi(int i,double val);
  void set_gxi(int x,int i);
  void set_hi(int i,int j,int val);
  void set_cati(int x,int i);
  void set_mother(Curves *x);

  void push_coord(double coord);
  void push_h(int i,int h);
  void push_gx(int gx);
  void push_cat(int cat);

  string get_id();
  vector<int> get_gx();
  vector<double> get_coords();
  vector<vector<int> > get_h();
  vector<int> get_cat();
  Curves * get_mother();

  int get_hi(int i,int j);
  double get_coordi(int i);
  int get_gxi(int i);
  int get_cati(int i);

  void erase_coordi(int i);

};

class Alfa{
private:
  vector<vector<vector<int> > >a;

public:
  Alfa(int L,int K);
  void set_ai(int i,int j,int z,int val);
  void push_a(int i,int j,int a);
  vector< vector<vector<int> > > get_a();
  int get_ai(int i,int j,int z);
};

class Curves{
private:
  vector<Dianisma> dianismata;
  string id;
  int m;
  vector<vector<double> >coords;

public:
  Curves(int L_grid,int m,string id);
  void set_h_dianisi(int i,int  L,int  K,Alfa a1,long long m,int M);
  void set_gx_dianisi(int i,int  L,int  K);
  void set_hashing_dianisi(int i,int  L,int  hash_categories);
  void push_cell_dianis(Hash_table *Hash,int j);
  void push_coordsi(int i,double x,double y);
  void push_dianis(Dianisma x);
  void push_coord_dianis(int i,double coord);
  int get_m();
  string get_id();
  double get_coordi_dianis(int i,int j);
  vector<vector<double> > get_coords();
  vector<double> get_coordsi(int i);
  void erase_diani_coordj(int i,int j);
  Dianisma get_diani(int i);
};
