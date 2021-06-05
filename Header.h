#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
using namespace std; 
static int m = 3, t = 2, b = 0; //3,2,0
static int dmin = 2 * t + 1;
static int n = pow(2, m) - 1;
static int k = n - dmin + 1;

vector<bool> *makegf(int n);
static vector<bool>*gf = makegf(n);

void print_vector(vector<int>v);
void print_bector(vector<bool>v);
vector<bool> operator ^ (vector<bool> a, vector<bool> b);
int in_GF_range(int power);
vector<bool> A(int num);
int L(vector<bool> v);
vector<bool> operator / (vector<bool> a, vector<bool> b);
vector<bool> operator * (vector<bool> a, vector<bool> b);
vector<int> cut(vector<int>a);
vector<int> operator + (vector<int> a, vector<int> b);
vector<int> operator * (vector<int> a, vector<int> b);
vector<int> operator / (vector<int> a, vector<int> b);
vector<int> operator % (vector<int> a, vector<int> b);
vector <int> shift_n_zeros_to(vector<int>a, int n);
vector<int> return_OM();
vector<int> return_syndroms(vector<int> code);
int Fx(vector<int> f, int x);
vector<int> return_Z(vector<int>S, vector<int> C);
vector<int> BMA(vector<int> S);
vector<int> PGZ(vector<int> S);
vector<int> CHEN(vector<int> C);
vector<int> FORNEY(vector<int>j, vector<int> Zx);
vector<bool> det(vector<int> *mtx);
vector<bool> pow_gf(int a, int power);