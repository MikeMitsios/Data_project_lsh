#include <vector>
#include <iostream>


using namespace std;


class Cell;
class Dianisma;


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
  vector<int> coords;
  string id;
  vector<vector<int> >h;
  vector<int> gx;
  vector<int> categ;

public:
  Dianisma(int L,int K);

  void set_id(string id);
  void set_coordi(int i,int val);
  void set_gxi(int x,int i);
  void set_hi(int i,int j,int val);
  void set_cati(int x,int i);

  void push_coord(int coord);
  void push_h(int i,int h);
  void push_gx(int gx);
  void push_cat(int cat);

  string get_id();
  vector<int> get_gx();
  vector<int> get_coords();
  vector<vector<int> > get_h();
  vector<int> get_cat();

  int get_hi(int i,int j);
  int get_coordi(int i);
  int get_gxi(int i);
  int get_cati(int i);

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
