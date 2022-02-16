#include "../cudd-complex-big/cplusplus/cuddObj.hh"
#include <iostream>
#include <algorithm>
#include <numbers>
#include <mpfr.h>
#include <cmath>
#include <ctime>
#include <ratio>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <typeinfo>
#include <fstream>
using namespace std::chrono;
#define RND_TYPE MPFR_RNDN

void adjustPrecision(Cudd& mgr, int n, int num){
	if (n >= num){
		mpfr_set_default_prec(300 + (n - num)*100);
		CUDD_VALUE_TYPE epsilon;
		mpfr_init_set_si(epsilon.real, -1 * (200 + (n-num)*100) , RND_TYPE);
		mpfr_exp10(epsilon.real, epsilon.real, RND_TYPE);
		mpfr_init_set_si(epsilon.imag, 0, RND_TYPE);
		mgr.SetEpsilon(epsilon);
	}
}

ADD identity_matrix(ADD x, ADD y){
  return (~x * ~y) + (x * y);
}

ADD exchange_matrix(ADD x, ADD y){
  return (~x * y) + (x * ~y);
}

ADD hadamard_matrix(Cudd& mgr, ADD x, ADD y){
	CUDD_VALUE_TYPE v;
	mpfr_init_set_d(v.real, 1.0/sqrt(2), RND_TYPE);
	mpfr_init_set_d(v.imag, 0.0, RND_TYPE);
  	ADD ret = (~x + x * (~y - y)) * mgr.constant(v);
  	mpfr_clear(v.real);mpfr_clear(v.imag);
  	return ret;
}
ADD CNOT_matrix(ADD x, ADD y, ADD w, ADD z){
  return (~x * ~y * ~z * ~w) + (~x * w * ~y * z) + (x * ~w * y * z) + (x * w * y * ~z);
}

ADD Voc14And23_matrix(ADD x, ADD y, ADD w, ADD z){
  return (~x * ~y * ~z * ~w) + (~x * w * y * ~z) + (x * ~w * ~y * z) + (x * w * y * z);
}

ADD identity_n(Cudd &mgr, unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 1)
    return identity_matrix(x_vars[start], y_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return identity_n(mgr, start, mid, x_vars, y_vars) * identity_n(mgr, mid, end, x_vars, y_vars);
}
ADD hadamard_n(Cudd &mgr,unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 1)
    return hadamard_matrix(mgr, x_vars[start], y_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return hadamard_n(mgr, start, mid, x_vars, y_vars) * hadamard_n(mgr, mid, end, x_vars, y_vars);
}
ADD C_n(Cudd &mgr,unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 2)
    return CNOT_matrix(x_vars[start], y_vars[start], x_vars[start+1], y_vars[start+1]);
  unsigned int mid = (end - start)/2 + start;
  return C_n(mgr, start, mid, x_vars, y_vars) * C_n(mgr, mid, end, x_vars, y_vars);
}

ADD exchange_n(Cudd &mgr, unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 1)
    return exchange_matrix(x_vars[start], y_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return exchange_n(mgr, start, mid, x_vars, y_vars) * exchange_n(mgr, mid, end, x_vars, y_vars);
}

ADD Voc14And23_n(Cudd &mgr, unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 2)
    return Voc14And23_matrix(x_vars[start], y_vars[start], x_vars[start+1], y_vars[start+1]);
  unsigned int mid = (end - start)/2 + start;
  return Voc14And23_n(mgr, start, mid, x_vars, y_vars) * Voc14And23_n(mgr, mid, end, x_vars, y_vars);
}

std::string getBits(unsigned int n, int len){
  std::string s = "";
  if (n == 0){
    for (int i = 0; i < len; i++)
      s += '0';
  }
  while (n != 0){
    int tmp = n % 2;
    s += (tmp == 0) ? '0' : '1';
    n = n/2;
  }
  while (s.length() < (unsigned int)len)
    s += '0';
  std::reverse(s.begin(), s.end());
  return s;
}


ADD D_n(Cudd& mgr, unsigned int start, unsigned int end, unsigned int size, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){
	unsigned int N = end - start;
	CUDD_VALUE_TYPE w;
	mpfr_init(w.real);
	mpfr_init(w.imag);
	mpfr_const_pi(w.real, RND_TYPE);
	mpfr_const_pi(w.imag, RND_TYPE);
	mpfr_mul_d(w.real, w.real, 2.0/pow(2, size), RND_TYPE);
	mpfr_mul_d(w.imag, w.imag, 2.0/pow(2, size), RND_TYPE);

	ADD ans = mgr.addZero();
	unsigned int counter = pow(2, N);
	for (unsigned int i = 0; i < counter; i++){
		ADD tmp_ans = mgr.addOne();
		unsigned int tmp_i = i;
		unsigned int count = 0;

		while (count < N){
			int j = tmp_i % 2;
			tmp_i /= 2;

			if (j == 1){
				tmp_ans = tmp_ans * x_vars[start + N - 1 - count] * y_vars[start + N - 1 - count];
			}else{
				tmp_ans = tmp_ans * ~x_vars[start + N - 1 - count] * ~y_vars[start + N - 1 - count];
			}
			count++;
		}
		CUDD_VALUE_TYPE val;
		mpfr_init_set(val.real, w.real, RND_TYPE);
		mpfr_init_set(val.imag, w.imag, RND_TYPE);
		mpfr_mul_si(val.real, val.real, (i) % ((unsigned int)pow(2, size)), RND_TYPE);
		mpfr_mul_si(val.imag, val.imag, (i) % ((unsigned int)pow(2, size)), RND_TYPE);
		mpfr_cos(val.real, val.real, RND_TYPE);
		mpfr_sin(val.imag, val.imag, RND_TYPE);
		ans = ans + tmp_ans * mgr.constant(val);
		mpfr_clear(val.real); mpfr_clear(val.imag);
	}

	mpfr_clear(w.real); mpfr_clear(w.imag);

	return ans;
}


ADD MkCyclicKMatrix(Cudd& mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, unsigned int N, unsigned int start){
	ADD ans = mgr.addZero();

	for (unsigned int i = 0; i < pow(2, N)/2; i++){
		ADD tmp_ans = mgr.addOne();

		unsigned int tmp_i = i;
		unsigned int val_i = 2*i;
		unsigned int count = 0;
		while (count < N){
			int tmp_j = tmp_i % 2;
			tmp_i /= 2;

			int val_j = val_i % 2;
			val_i /= 2;

			if (tmp_j == 1){
				tmp_ans = tmp_ans * x_vars[N-1-count + start];
			}else{
				tmp_ans = tmp_ans * (~x_vars[N-1-count + start]);
			}

			if (val_j == 1){
				tmp_ans = tmp_ans * y_vars[N-1-count + start];
			}else{
				tmp_ans = tmp_ans * (~y_vars[N-1-count + start]);
			}

			count++;

		}
		ans = ans + tmp_ans;
	}

	for (unsigned int i = 0; i < pow(2, N)/2; i++){
		ADD tmp_ans = mgr.addOne();

		unsigned int tmp_i = i;
		unsigned int val_i = 2*i + 1;
		unsigned int count = 0;
		while (count < N){
			int tmp_j = tmp_i % 2;
			tmp_i /= 2;

			int val_j = val_i % 2;
			val_i /= 2;

			if (tmp_j == 1){
				tmp_ans = tmp_ans * x_vars[N-1-count + start];
			}else{
				tmp_ans = tmp_ans * (~x_vars[N-1-count + start]);
			}

			if (val_j == 1){
				tmp_ans = tmp_ans * y_vars[N-1-count + start];
			}else{
				tmp_ans = tmp_ans * (~y_vars[N-1-count + start]);
			}

			count++;

		}
		ans = ans + tmp_ans;
	}

	return ans;
}



ADD Fourier(Cudd& mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, 
	std::vector<ADD>& w_vars, std::vector<ADD>& z_vars, unsigned int N, 
	unsigned int start){
	if (N == 1){
		return hadamard_matrix(mgr, x_vars[start], y_vars[start]);
	}

	ADD K = MkCyclicKMatrix(mgr, z_vars, y_vars, N, start);
	ADD F_recurse = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, N-1, start+1);

	ADD F = (~x_vars[start]*~y_vars[start] + x_vars[start]*y_vars[start]) * F_recurse;
	F = F.SwapVariables(x_vars, w_vars);
	F = F.SwapVariables(y_vars, z_vars);

	ADD Id = identity_n(mgr, start + 1, N, x_vars, y_vars);
	ADD D = D_n(mgr, start + 1, N, N, x_vars, y_vars);
	ADD ID = (~y_vars[start]*Id) + y_vars[start] * (~x_vars[start] - x_vars[start]) * D;
	ID = ID.SwapVariables(y_vars, w_vars);
	ADD IDF = ID.MatrixMultiply(F, w_vars);
	ADD ans = IDF.MatrixMultiply(K, z_vars);
	CUDD_VALUE_TYPE val;
	mpfr_init_set_d(val.real, 1.0/sqrt(2), RND_TYPE);
	mpfr_init_set_d(val.imag, 1.0/sqrt(2), RND_TYPE);
	ADD constant = mgr.constant(val);
	mpfr_clear(val.real); mpfr_clear(val.imag);
	ans = ans * constant;
	return ans;
}



unsigned int Fourier(Cudd& mgr, int n){
	std::vector<ADD> x_vars, w_vars, z_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(4*i));
    w_vars.push_back(mgr.addVar(4*i+1));
    z_vars.push_back(mgr.addVar(4*i+2));
    y_vars.push_back(mgr.addVar(4*i+3));
  }

  unsigned int N = (unsigned int)pow(2,n);
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD ans = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, N, 0);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
  ans.print(4*N,2);
  return ans.nodeCount();
}

int main (int argc, char** argv)
{
	if (argc < 4)
      return 0;
	Cudd mgr(0,0);
	if (strcmp(argv[3], "enable") == 0)
		mgr.AutodynEnable();
	unsigned int nodeCount = 0;

	srand(time(NULL));

	if (strcmp(argv[1], "fourier") == 0)
      nodeCount = Fourier(mgr, atoi(argv[2])); 

	return 0;
}