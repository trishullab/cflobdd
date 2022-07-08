#include "../cudd-fourier/cplusplus/cuddObj.hh"
#include <iostream>
#include <algorithm>
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
#include <random>
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
	// CUDD_VALUE_TYPE v;
	// mpfr_init_set_d(v.real, 1.0/sqrt(2), RND_TYPE);
	// mpfr_init_set_d(v.imag, 0.0, RND_TYPE);
 //  	ADD ret = (~x + x * (~y - y)) * mgr.constant(v);
 //  	mpfr_clear(v.real);mpfr_clear(v.imag);
	CUDD_VALUE_TYPE val;
	val.val = 0; val.base = 2;
	CUDD_VALUE_TYPE val2;
	val2.val = 1; val2.base = 2;
	val.is_complex_assigned = 0;
	val2.is_complex_assigned = 0;
	ADD ret = (~x + x * ~y) * mgr.constant(val);
	ret = ret + (x * y) * mgr.constant(val2);
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
	// CUDD_VALUE_TYPE w;
	// mpfr_init(w.real);
	// mpfr_init(w.imag);
	// mpfr_const_pi(w.real, RND_TYPE);
	// mpfr_const_pi(w.imag, RND_TYPE);
	// mpfr_mul_d(w.real, w.real, 2.0/pow(2, size), RND_TYPE);
	// mpfr_mul_d(w.imag, w.imag, 2.0/pow(2, size), RND_TYPE);

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
		// mpfr_init_set(val.real, w.real, RND_TYPE);
		// mpfr_init_set(val.imag, w.imag, RND_TYPE);
		// mpfr_mul_si(val.real, val.real, (i) % ((unsigned int)pow(2, size)), RND_TYPE);
		// mpfr_mul_si(val.imag, val.imag, (i) % ((unsigned int)pow(2, size)), RND_TYPE);
		// mpfr_cos(val.real, val.real, RND_TYPE);
		// mpfr_sin(val.imag, val.imag, RND_TYPE);
		val.val = (i) % ((unsigned int)pow(2, size));
		val.base = ((unsigned int)pow(2, size));
		val.is_complex_assigned = 0;
		ans = ans + tmp_ans * mgr.constant(val);
		// mpfr_clear(val.real); mpfr_clear(val.imag);
	}

	// mpfr_clear(w.real); mpfr_clear(w.imag);

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

		unsigned int tmp_i = i + pow(2, N)/2;
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
	std::cout << "start: " << N << std::endl;
	if (N == 1){
		return hadamard_matrix(mgr, x_vars[start], y_vars[start]);
	}

	ADD K = MkCyclicKMatrix(mgr, z_vars, y_vars, N, start);
	// K.print(4*N,2);
	ADD F_recurse = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, N-1, start+1);

	ADD F = (~x_vars[start]*~y_vars[start] + x_vars[start]*y_vars[start]) * F_recurse;
	F = F.SwapVariables(x_vars, w_vars);
	F = F.SwapVariables(y_vars, z_vars);
	ADD Id = identity_n(mgr, start + 1, x_vars.size(), x_vars, y_vars);
	ADD D = D_n(mgr, start + 1, x_vars.size(), N, x_vars, y_vars);
	// if (N == 3){
	// 	D.print(4*N,2);
	// }
	CUDD_VALUE_TYPE val;
	val.val = pow(2, N)/2; val.base = pow(2, N);val.is_complex_assigned = 0;
	ADD ID = (~y_vars[start]*Id) + (y_vars[start] * (~x_vars[start] + mgr.constant(val) * x_vars[start]) * D);
	ID = ID.SwapVariables(y_vars, w_vars);
	ADD IDF = ID.MatrixMultiply(F, w_vars);
	ADD ans = IDF.MatrixMultiply(K, z_vars);
	// CUDD_VALUE_TYPE val;
	// mpfr_init_set_d(val.real, 1.0/sqrt(2), RND_TYPE);
	// mpfr_init_set_d(val.imag, 0.0, RND_TYPE);
	// ADD constant = mgr.constant(val);
	// ans = ans * constant;
	// mpfr_clear(val.real); mpfr_clear(val.imag);
	ans = ans.ConvertToBase(pow(2, N));
	std::cout << "end: " << N << std::endl;
	return ans;
}



unsigned int Fourier(Cudd& mgr, int n){
	std::vector<ADD> x_vars, w_vars, z_vars, y_vars;
	unsigned int N = pow(2, n);
  for (unsigned int i = 0; i < N; i++){
    x_vars.push_back(mgr.addVar(4*i));
    y_vars.push_back(mgr.addVar(4*i+1));
    w_vars.push_back(mgr.addVar(4*i+2));
    z_vars.push_back(mgr.addVar(4*i+3));
  }

  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD ans = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, N, 0);
  ans = ans.ConvertToComplex();
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
  // ans.print(4*N,2);
  return ans.nodeCount();
}

std::string getBitString(unsigned int i, int n){
	std::string s(n, '0');
	int count = n - 1;
	while (i > 0){
		int k = i % 2;
		i = i / 2;
		s[count] = (k == 0) ? '0' : '1';
		count--;
	}
	return s;
}

ADD getVector(Cudd&mgr, std::string index_s, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){
	ADD F = mgr.addOne();
	F = F.SetToComplex();
  for (int i = 0; i < index_s.length(); i++){
  	if (i % 2 == 0)
  	{
  		if (index_s[i] == '0')
  		{
  			F *= ~x_vars[i/2];
  		}else{
  			F *= x_vars[i/2];
  		}
  	}else{
  		if (index_s[i] == '0')
  		{
  			F *= ~y_vars[i/2];
  		}else{
  			F *= y_vars[i/2];
  		}
  	}
  }
  return F;
}

ADD ShorsMainAlgo(Cudd&mgr, int l, ADD U, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars,
	std::vector<ADD>& w_vars, std::vector<ADD>& z_vars){
	std::cout << "Algo start.." << std::endl;
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD C = CNOT_matrix(x_vars[0], y_vars[0], w_vars[0], z_vars[0]);
  C = C.SetToComplex();
  unsigned int N = l;
  for (unsigned int i = 1; i < N; i++){
  	ADD tmp = CNOT_matrix(x_vars[i], y_vars[i], w_vars[i], z_vars[i]);
  	tmp = tmp.SetToComplex();
    C *= tmp;
  }
  C = C.SetToComplex();
  
  ADD ans = ~z_vars[0];
  ans = ans.SetToComplex();
  for (unsigned int i = 1; i < N; i++){
  	ADD tmp = ~z_vars[i];
  	tmp = tmp.SetToComplex();
    ans *= tmp;
  }
  ans = ans.SetToComplex();
  std::vector<ADD> mult_array;
  for (unsigned int i = 0; i < N; i++){
  	mult_array.push_back(y_vars[i]);
  	mult_array.push_back(z_vars[i]);
  }
  std::vector<ADD> swap_array;
  for (unsigned int i = 0; i < N; i++){
  	swap_array.push_back(x_vars[i]);
  	swap_array.push_back(w_vars[i]);
  }
  // std::cout << "Step 1 : " << ans.nodeCount() << std::endl;
  ans = C.MatrixMultiply(ans, mult_array);
  ans = ans.SwapVariables(swap_array, mult_array);
  ans.SetToComplex();
  CUDD_VALUE_TYPE val;
  val.is_complex_assigned = 1;
  mpfr_init_set_si(val.real, N, RND_TYPE);
  mpfr_exp2(val.real, val.real, RND_TYPE);
  mpfr_d_div(val.real, 1.0, val.real, RND_TYPE);
  mpfr_init_set_si(val.imag, 1, RND_TYPE);
  ADD Q = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, N, 0);
  Q = Q.ConvertToComplex();
  ADD V = mgr.constant(val);
  V = V.SetToComplex();
  // Q = Q * V;
  ADD QU = Q * U;
  // U.print(4*N,2);
  // std::cout << "Step 2 : " << ans.nodeCount() << std::endl;
  ans = QU.MatrixMultiply(ans, swap_array);
  ans.print(4*N, 2);
  ans = ans.SwapVariables(swap_array, mult_array);
  std::cout << "QU done" << std::endl; 
  // std::cout << "Step 3 : " << ans.nodeCount() << std::endl;
  ans = ans.SquareTerminalValues();
  // std::cout << "Step 4 : " << ans.nodeCount() << std::endl;
  ans = ans.UpdatePathInfo(2, 2*N);
  // ans.PrintPathInfo();
  std::cout << "Path updated" << std::endl;
  high_resolution_clock::time_point mid = high_resolution_clock::now();
  unsigned int iter = 1;
  while (iter <= 2 * N)
	{
		std::string s = "";
		s = ans.SamplePath(2*N, 2, "simons");
		// std::cout << "iter: " << iter << std::endl;
		s = s.substr(0, N);

		iter++;
	}

  high_resolution_clock::time_point end = high_resolution_clock::now();
  return ans;
  return Q;
}

unsigned int ShorsAlgo(Cudd& mgr){
	int n = 4;
	int M = 15;
	
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> dist(1, M - 1);
	auto start = high_resolution_clock::now();
	unsigned int a = 1;
	while (a == 1 || a % 2 == 0){
		a = dist(rng);
	}
	a = 7;
	std::cout << "a: " << a << std::endl;
	int l = 2 * n; // M^2 <= Q < 2M^2 and Q = 2^l
	int level = ceil(log2(l));
	unsigned int N = l;
  std::vector<ADD> x_vars, y_vars, w_vars, z_vars; // x & y belong to input and w & z belong to output
  for (unsigned int i = 0; i < N; i++){
    x_vars.push_back(mgr.addVar(4*i));
    y_vars.push_back(mgr.addVar(4*i+1));
    w_vars.push_back(mgr.addVar(4*i+2));
    z_vars.push_back(mgr.addVar(4*i+3));
  }

	// use w and z vars for creating F
	std::string s = getBitString(1, 2*l);
	ADD F = mgr.addZero();
	F = F.SetToComplex();
	// F.print(4*N,2);
	F = F + getVector(mgr, s, w_vars, z_vars);
	unsigned int f_i = 1;
	for (unsigned int i = 1; i < pow(2, l); i++){
		std::string index_s(2 * l, '0');
		f_i = (f_i * a) % M;
		std::string i_s = getBitString(i, l);
		std::string f_i_s = getBitString(f_i, l);
		int k = 0;
		for (unsigned int j = 0; j < index_s.length(); j+=2){
			index_s[j] = i_s[k];
			index_s[j + 1] = f_i_s[k];
			k++;
		}
		F = F + getVector(mgr, index_s, w_vars, z_vars);
	}
	F = F.SetToComplex();
	std::cout << "l: " << l << std::endl;
	// F.print(4*N,2);
	auto ans = ShorsMainAlgo(mgr, l, F, x_vars, y_vars, w_vars, z_vars);
	auto end = high_resolution_clock::now();
	return 0;
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
  if (strcmp(argv[1], "shors") == 0)
      nodeCount = ShorsAlgo(mgr); 

	return 0;
}
