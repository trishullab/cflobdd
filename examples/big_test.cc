#include "../cudd-big/cplusplus/cuddObj.hh"
#include <iostream>
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
using namespace std::chrono;
#define RND_TYPE MPFR_RNDN

void adjustPrecision(Cudd& mgr, int n, int num){
	if (n >= num){
		mpfr_set_default_prec(300 + (n - num)*100);
		CUDD_VALUE_TYPE epsilon;
		mpfr_init_set_si(epsilon.t_val, -1 * (200 + (n-num)*100) , RND_TYPE);
		mpfr_exp10(epsilon.t_val, epsilon.t_val, RND_TYPE);
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
	mpfr_init_set_d(v.t_val, 1.0/sqrt(2), RND_TYPE);
  	ADD ret = (~x + x * (~y - y)) * mgr.constant(v);
  	mpfr_clear(v.t_val);
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



ADD PermutationMatrix(Cudd& mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, unsigned int N){
	if (N <= 8){

		unsigned vals_size = pow(2, N);
		ADD ans = mgr.addZero();
		for (unsigned int i = 0; i < pow(2, N); i++){
			ADD tmp_ans = mgr.addOne();
			unsigned int y_val = rand() % vals_size;

			std::string x_string = getBits(i, N);
			std::string y_string = getBits(y_val, N);

			for (unsigned int j = 0; j < N; j++){
				if (x_string[j] == '1'){
					tmp_ans = tmp_ans * x_vars[j];
				}else{
					tmp_ans = tmp_ans * (~x_vars[j]);
				}
			}

			for (unsigned int j = 0; j < N; j++){
				if (y_string[j] == '1'){
					tmp_ans = tmp_ans * y_vars[j];
				}else{
					tmp_ans = tmp_ans * (~y_vars[j]);
				}
			}

			ans += tmp_ans;

		}

		return ans;

	}
	int rand_val = rand() % 3;
	if (rand_val == 0){
		return Voc14And23_n(mgr, 0, N, x_vars, y_vars);
	}else if (rand_val == 1){
		return C_n(mgr, 0, N, x_vars, y_vars);
	}else{
		return exchange_n(mgr, 0, N, x_vars, y_vars);
	}
}


unsigned int exclusive_or(Cudd &mgr, int n){
  high_resolution_clock::time_point start = high_resolution_clock::now();
    BDD f = mgr.bddVar();
    for (unsigned int i = 1; i < pow(2, n); i++){
      	f = f ^ mgr.bddVar();
    }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " << f.nodeCount() << " time: " << time_taken.count() << std::endl;
    return f.nodeCount();
}

unsigned int hi_sum(Cudd &mgr, int n){
  std::vector<ADD> x_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(2*i));
    y_vars.push_back(mgr.addVar(2*i + 1));
  }
  ADD H = hadamard_n(mgr, 0, pow(2, n), x_vars, y_vars);
  ADD I = identity_n(mgr, 0, pow(2, n), x_vars, y_vars);
  std::cout << "Matrices created .. " << std::endl;
    high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD ans = H + I;
  for (unsigned int i = 1; i < 1024; i++){
    ans += H + I;
  }
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> time_taken = duration_cast<duration<double>>(end - start);
    std::cout << time_taken.count() << " " << H.nodeCount() << " " <<  ans.nodeCount() << std::endl;
  return ans.nodeCount();
}

std::string compute_xor(std::string inp, std::string s){
  std::string out_string = "";
  for (unsigned int i = 0; i < s.length(); i++){
    if (s[i] == '1')
      out_string += (inp[i] == '1') ? '0' : '1';
    else
      out_string += inp[i];
  }
  return out_string;
}

ADD SMatrix(Cudd& mgr, std::string s, std::vector<ADD>& w_vars){
	int index = -1;
	for (unsigned int i = 0; i < s.length(); i++){
		if (s[i] == '1'){
			index = i;
			break;
		}
	}
	ADD ans;
	if (index == -1){
		ans = mgr.addOne();
	}
	else{
		ans = w_vars[index] - (~w_vars[index]);
	}
	return ans;
}

ADD Create2To1Func(Cudd& mgr, ADD F1, ADD F2, std::string s1, std::string s2, unsigned int n, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars, std::vector<ADD>& x_vars){
	ADD S1 = SMatrix(mgr, s1, w_vars);
	std::vector<ADD> new_w, new_z, new_x;
	new_w.insert(new_w.end(), w_vars.begin() + n/2, w_vars.begin() + n);
	new_z.insert(new_z.end(), z_vars.begin() + n/2, z_vars.begin() + n);
	new_x.insert(new_x.end(), x_vars.begin() + n/2, x_vars.begin() + n);
	ADD S2 = SMatrix(mgr, s2, new_w);
	ADD F1Prime = F1 * S1;
	ADD F2Prime = F2 * S2;
	ADD F1Permute = F1Prime;
	ADD F2Permute = F2Prime;

	if (!(s1.find('1') == std::string::npos || s2.find('1') == std::string::npos)){
		std::vector<ADD> first_half_z, first_half_x;
		first_half_z.insert(first_half_z.end(), z_vars.begin(), z_vars.begin() + n/2);
		first_half_x.insert(first_half_x.end(), x_vars.begin(), x_vars.begin() + n/2);
		ADD P1 = PermutationMatrix(mgr, x_vars, z_vars, n/2);
		ADD P2 = PermutationMatrix(mgr, new_x, new_z, n/2);

		F1Permute = F1Permute.MatrixMultiply(P1, first_half_z);
		F2Permute = F2Permute.MatrixMultiply(P2, new_z);

		F1Permute = F1Permute.SwapVariables(first_half_x, first_half_z);
		F2Permute = F2Permute.SwapVariables(new_z, new_x);
	}

	ADD t1 = F1Prime * F2Prime;
	t1 = t1.SimonsRemoveMinusOne();
	ADD t2 = F1Permute * F2Permute;
	t2 = t2.SimonsRemoveOne();
	ADD ans = t1 + t2;
	return ans;
}

ADD CreateFnMatrix(Cudd &mgr, std::string s, unsigned int n, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){

	if (n <= 2){
	  unsigned long long int numVals = pow(2,n);
	  ADD ans = mgr.addZero();
	  std::unordered_map<std::string, std::string> map_values;
	  for (unsigned long long int i = 0; i < pow(2,n); i++){
	  	unsigned long long int val = rand() % numVals; 
	    std::string inp_string = getBits(i, n); 
	    std::string xor_string = compute_xor(inp_string, s);
	    std::string val_string;
	    if (map_values.find(xor_string) != map_values.end())
	      val_string = map_values[xor_string]; 
	    else{
	      val_string = getBits(val, n);
	      map_values[inp_string] = val_string;
	    }
	    ADD tmp_ans = mgr.addOne();
	    for (unsigned int j = 0; j < n; j++){
	      if (inp_string[j] == '1')
			tmp_ans *= w_vars[j];
	      else
			tmp_ans *= ~w_vars[j];
	    }
	    for (unsigned int j = 0; j < n; j++){
	      if (val_string[j] == '1')
			tmp_ans *= z_vars[j];
	      else
			tmp_ans *= ~z_vars[j];
	    }
	    ans += tmp_ans;
	  }
	  return ans;
	}
	else
	{
		ADD F1 = CreateFnMatrix(mgr, s.substr(0, s.length()/2), n/2, w_vars, z_vars, x_vars, y_vars);
		std::vector<ADD> new_x, new_y, new_w, new_z;
		new_x.insert(new_x.end(), x_vars.begin() + n/2, x_vars.begin() + n);
		new_y.insert(new_y.end(), y_vars.begin() + n/2, y_vars.begin() + n);
		new_w.insert(new_w.end(), w_vars.begin() + n/2, w_vars.begin() + n);
		new_z.insert(new_z.end(), z_vars.begin() + n/2, z_vars.begin() + n);
		ADD F2 = CreateFnMatrix(mgr, s.substr(s.length()/2), n/2, new_w, new_z, new_x, new_y);
		ADD ans = Create2To1Func(mgr, F1, F2, s.substr(0, s.length()/2), s.substr(s.length()/2), n, w_vars, z_vars, x_vars);
		return ans;
	}
}

unsigned int simons(Cudd &mgr, int n){

	adjustPrecision(mgr, n, 6);
	unsigned int N = pow(2, n);
  std::vector<ADD> x_vars, y_vars, w_vars, z_vars; // x & y belong to input and w & z belong to output
  for (unsigned int i = 0; i < N; i++){
    x_vars.push_back(mgr.addVar(4*i));
    y_vars.push_back(mgr.addVar(4*i+1));
    w_vars.push_back(mgr.addVar(4*i+2));
    z_vars.push_back(mgr.addVar(4*i+3));
  }
  srand(2);
  std::string s(N, '0');
  for (unsigned int i = 0; i < N; i++){
  	s[i] = (rand() % 2) ? '1' : '0';
  }
  // s = "1110";
  ADD U = CreateFnMatrix(mgr, s, N, w_vars, z_vars, x_vars, y_vars);
  // U = U.SwapVariables(z_vars, w_vars);
  //ADD C = C_n(mgr, 0, pow(2, n), x_vars, y_vars, w_vars, z_vars);
  std::cout << "Algo start.." << std::endl;
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD C = CNOT_matrix(x_vars[0], y_vars[0], w_vars[0], z_vars[0]);
  for (unsigned int i = 1; i < N; i++){
    C *= CNOT_matrix(x_vars[i], y_vars[i], w_vars[i], z_vars[i]);
  }
  
  ADD ans = ~z_vars[0];
  for (unsigned int i = 1; i < N; i++){
    ans *= ~z_vars[i];
  }
  std::vector<ADD> mult_array;
  // mult_array.insert(mult_array.end(), y_vars.begin(), y_vars.end());
  // mult_array.insert(mult_array.end(), z_vars.begin(), z_vars.end());
  for (unsigned int i = 0; i < N; i++){
  	mult_array.push_back(y_vars[i]);
  	mult_array.push_back(z_vars[i]);
  }
  std::vector<ADD> swap_array;
  // swap_array.insert(swap_array.end(), x_vars.begin(), x_vars.end());
  // swap_array.insert(swap_array.end(), w_vars.begin(), w_vars.end());
  for (unsigned int i = 0; i < N; i++){
  	swap_array.push_back(x_vars[i]);
  	swap_array.push_back(w_vars[i]);
  }
  ans = C.MatrixMultiply(ans, mult_array);
  // ans = ans.SwapVariables(swap_array, mult_array);
  CUDD_VALUE_TYPE val;
  mpfr_init_set_si(val.t_val, N/2, RND_TYPE);
  mpfr_exp2(val.t_val, val.t_val, RND_TYPE);
  mpfr_d_div(val.t_val, 1.0, val.t_val, RND_TYPE);
  ADD H = hadamard_n(mgr, 0, N, x_vars, y_vars) * mgr.constant(val);
  mpfr_clear(val.t_val);
  ADD tmp = identity_n(mgr, 0, N, x_vars, y_vars) * U;
  ADD HU = H * U;
  std::cout << ans.nodeCount() << " " << HU.nodeCount() << std::endl;
  ans = HU.MatrixMultiply(ans, swap_array);
  ans = ans.SwapVariables(swap_array, mult_array);
  ans = ans.SquareTerminalValues();
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "string s: " << s << " U nodeCount: " << U.nodeCount() << std::endl;
  std::cout << "nodeCount: " << ans.nodeCount() << " time: " << time_taken.count() << std::endl;
  return ans.nodeCount();
}




ADD MkU_w(Cudd &mgr, std::string s, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){
  int n = x_vars.size();
  ADD I = identity_n(mgr, 0, n, x_vars, y_vars);
  CUDD_VALUE_TYPE val;
  mpfr_init_set_si(val.t_val, -2, RND_TYPE);
  ADD F = mgr.constant(val);
  for (int i = 0; i < n; i++){
    if (s[i] == '0'){
      F *= ~x_vars[i] * ~y_vars[i];
    }else{
    	// std::cout << "go HEad" << std::endl;
      F *= x_vars[i] * y_vars[i];
    }
  }
  return I + F;
}

ADD MkU_s(Cudd &mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){//, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars){
  int n = x_vars.size();
  ADD I = identity_n(mgr, 0, n, x_vars, y_vars);
  CUDD_VALUE_TYPE val;
  mpfr_init_set_si(val.t_val, n, RND_TYPE);
  mpfr_exp2(val.t_val, val.t_val, RND_TYPE);
  mpfr_d_div(val.t_val, 2.0, val.t_val, RND_TYPE);
  // mpfr_printf("val: %.128Rf\n", val.t_val);
  ADD X = (mgr.constant(val) - I);// * mgr.constant(1.0/pow(2,n));
  // X.print(4*n, 2);
  return X;
}

unsigned int grover(Cudd &mgr, int n){

	adjustPrecision(mgr, n, 8);
	
  std::string s;
  for (int i = 0; i < pow(2, n); i++){
    s += (rand()%2 == 0) ? '0' :'1';
  }
  std::vector<ADD> x_vars, w_vars, z_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(4*i));
    w_vars.push_back(mgr.addVar(4*i+1));
    z_vars.push_back(mgr.addVar(4*i+2));
    y_vars.push_back(mgr.addVar(4*i+3));
  }
  unsigned int N = (unsigned int)pow(2,n);
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD U_w = MkU_w(mgr,s, w_vars, z_vars);
  ADD U_s = MkU_s(mgr, x_vars, w_vars);//, z_vars, y_vars);//* mgr.constant(1.0/pow(2,n));
  ADD M = U_s.MatrixMultiply(U_w, w_vars); 
  std::cout << "M leaf counts: " << M.CountLeaves() << std::endl;
  M = M.SwapVariables(w_vars, z_vars);
  double const pi = 4 * std::atan(1);
  unsigned long long int iters = (unsigned long long int)floor(pi * 0.25 * pow(2, N/2)); 
  CUDD_VALUE_TYPE val;
  mpfr_init_set_si(val.t_val, N/2, RND_TYPE);
  mpfr_exp2(val.t_val, val.t_val, RND_TYPE);
  mpfr_d_div(val.t_val, 1.0, val.t_val, RND_TYPE);
  ADD ans = mgr.constant(val);
  std::cout << "M count: " << M.nodeCount() << " ans count: " << ans.nodeCount() << std::endl;
  std::cout << "iters count: " << iters << std::endl;
  for (unsigned long long int i = 0; i < iters; i++){
    ans = M.MatrixMultiply(ans, w_vars);// * mgr.constant(1.0/pow(2,n));
    // std::cout << "i : " << i << " ans count: " << ans.nodeCount() << " M leaves: " << M.CountLeaves() << " ans leaves: " << ans.CountLeaves() << std::endl;
    ans = ans.SwapVariables(x_vars, w_vars);
  }
  ans = ans.SwapVariables(x_vars, w_vars);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "string: " << s << std::endl;
  std::cout << "nodeCount: " <<  ans.nodeCount() << " leaf count: " << ans.CountLeaves() << " time_taken: " << time_taken.count() << std::endl;
  // ans.PrintMinterm();
  // ans.print(4*n,2);
  return ans.nodeCount();
}

unsigned int GHZ(Cudd &mgr, unsigned int n){

	adjustPrecision(mgr, n, 8);
  unsigned int N = pow(2,n) + 1;
  std::vector<ADD> x_vars, w_vars, z_vars, y_vars;
  for (unsigned int i = 0; i < N; i++){
    x_vars.push_back(mgr.addVar(3*i));
    y_vars.push_back(mgr.addVar(3*i+1));
    z_vars.push_back(mgr.addVar(3*i+2));
    /*
    y_vars.push_back(mgr.addVar(4*i+3));
    */
  }
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD I = identity_n(mgr, 0, N-1, x_vars, y_vars);
  ADD I_shift = identity_n(mgr, 0, N-1, y_vars, z_vars);
  ADD U = I * (CNOT_matrix(x_vars[0], y_vars[0], x_vars[N-1], y_vars[N-1]));
  for (unsigned int i = 1; i < N-1; i++){
      ADD tmp = I_shift * (CNOT_matrix(y_vars[i], z_vars[i], y_vars[N-1], z_vars[N-1]));
      U = U.MatrixMultiply(tmp, y_vars);
      U = U.SwapVariables(y_vars, z_vars);
  }
  ADD ans =  y_vars[N-1];
  for (unsigned int i = 0; i < N-1; i++)
    ans = ans * ~y_vars[i];
  ADD H = hadamard_n(mgr, 0, N-1, x_vars, y_vars);
  ADD I_0 = identity_matrix(x_vars[N-1], y_vars[N-1]);
  H = H * I_0; //* mgr.constant(1.0/pow(2,N-1));
  ans = H.MatrixMultiply(ans, y_vars);
  ans = ans.SwapVariables(x_vars, y_vars);
  ans = U.MatrixMultiply(ans, y_vars);
  ans = ans.SwapVariables(x_vars, y_vars);
  H = hadamard_n(mgr, 0, N-1, x_vars, y_vars);
  ADD H_0 = hadamard_matrix(mgr, x_vars[N-1], y_vars[N-1]);
  H = H * H_0;
  ans = H.MatrixMultiply(ans, y_vars);
  // ans.print(3*N,2);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
  return ans.nodeCount();
}

unsigned int BV(Cudd &mgr, int n){

	adjustPrecision(mgr, n, 8);
  unsigned int N = (unsigned int)pow(2, n);
  std::vector<ADD> x_vars, z_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n)+1; i++){
    x_vars.push_back(mgr.addVar(3*i));
    y_vars.push_back(mgr.addVar(3*i+1));
    z_vars.push_back(mgr.addVar(3*i+2));
  }
  std::string s;
  for (int i = 0; i < pow(2, n); i++){
    s += (rand()%2 == 0) ? '0' :'1';
  }
  long long int index = -1;
  for (unsigned int i = 0; i < s.length(); i++){
    if (s[i] == '1'){
      index = i;
      break;
    }
  }
  std::cout << "s created... " << s << std::endl;
  ADD I = identity_n(mgr, 0, N, x_vars, y_vars);
  std::cout << "I created..." << std::endl;
  ADD I_shift = identity_n(mgr, 0, N, y_vars, z_vars);
  std::cout << "I_shift created..." << std::endl;
  ADD U;
  if (index == -1)
    U = I;
  else
    U = I * (CNOT_matrix(x_vars[index], y_vars[index], x_vars[N], y_vars[N]));
  for (unsigned int i = index+1; i < s.length() && index != -1; i++){
    if (s[i] == '1'){
      ADD tmp = I_shift * (CNOT_matrix(y_vars[i], z_vars[i], y_vars[N], z_vars[N]));
      U = U.MatrixMultiply(tmp, y_vars);
      U = U.SwapVariables(y_vars, z_vars);
    }
  }

  std::cout << "Starting BV..." << std::endl;
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD ans =  y_vars[N];
  for (unsigned int i = 0; i < N; i++)
    ans = ans * ~y_vars[i];
  ADD H = hadamard_n(mgr, 0, N, x_vars, y_vars);
  ADD I_0 = identity_matrix(x_vars[N], y_vars[N]);
  CUDD_VALUE_TYPE value;
  mpfr_init_set_si(value.t_val, N, RND_TYPE);
  mpfr_exp2(value.t_val, value.t_val, RND_TYPE);
  mpfr_d_div(value.t_val, 1.0, value.t_val, RND_TYPE);
  H = H * I_0; //* mgr.constant(value);
  ans = H.MatrixMultiply(ans, y_vars);
  //ans = ans * mgr.constant(1.0/pow(2,N/2));
  ans = ans.SwapVariables(x_vars, y_vars);
  ans = U.MatrixMultiply(ans, y_vars);
  ans = ans.SwapVariables(x_vars, y_vars);
  ans = H.MatrixMultiply(ans, y_vars);
  //ans = ans * mgr.constant(1.0/pow(2,N/2));

  //ans.print(4*N,2);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "string: " << s << std::endl;
  std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
  // ans.print(3*n,2);
  return ans.nodeCount();
}

ADD CreateBalancedFn(Cudd& mgr, unsigned int N, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, std::vector<ADD>& z_vars){
	// if (N <= 2){
	// 	unsigned int count = 0;
	// 	CUDD_VALUE_TYPE zero;
	// 	mpfr_init_set_si(zero.t_val, 0, RND_TYPE);
	// 	ADD ans = mgr.constant(zero);
	// 	mpfr_clear(zero.t_val);
	// 	for (unsigned int i = 0; i < pow(2, N-1); i++){
	// 		int f_x = (rand() % 2);
	// 		if (count >= pow(2, N)/2)
	// 			f_x = 0;
	// 		if (f_x == 1)
	// 			count++;
	// 		std::string x_string = getBits(i, N-1);
	// 		CUDD_VALUE_TYPE one;
	// 		mpfr_init_set_si(one.t_val, 1, RND_TYPE);
	// 		ADD tmp_ans = mgr.constant(one);
	// 		mpfr_clear(one.t_val);
	// 		for (unsigned int j = 0; j < N-1; j++){
	// 			if (x_string[j] == '1')
	// 				tmp_ans *= x_vars[j] * y_vars[j];
	// 			else
	// 				tmp_ans *= ~x_vars[j] * ~y_vars[j];
	// 		}

	// 		if (f_x == 0){
	// 			ans += (tmp_ans * ~x_vars[N-1] * ~y_vars[N-1]) + (tmp_ans * x_vars[N-1] * y_vars[N-1]);
	// 		}
	// 		else{
	// 			ans += (tmp_ans * ~x_vars[N-1] * y_vars[N-1]) + (tmp_ans * x_vars[N-1] * ~y_vars[N-1]);				
	// 		}

	// 	}
	// 	return ans;
	// }
	// else{
	// 	ADD F = CreateBalancedFn(mgr, N/2, x_vars, y_vars, z_vars);
	// 	F.print(3*N, 2);
	// 	std::vector<ADD> new_x_vec(x_vars.begin() + N/2, x_vars.begin() + N);
	// 	std::vector<ADD> new_y_vec(y_vars.begin() + N/2, y_vars.begin() + N);
	// 	ADD P = PermutationMatrix(mgr, new_x_vec, new_y_vec, N/2);
	// 	F = F * P;
	// 	std::cout << "F created" << std::endl;
	// 	ADD permute_matrix = PermutationMatrix(mgr, z_vars, y_vars, N);
	// 	permute_matrix = permute_matrix.SwapVariables(y_vars, z_vars);
	// 	std::cout << "P created" << std::endl;
	// 	F = F.MatrixMultiply(permute_matrix, y_vars);
	// 	F = F.SwapVariables(y_vars, z_vars);
	// 	F.print(3 * N, 2);
	// 	return F;
	// }

	std::string s(N, '0');
	for (unsigned int i = 0; i < N; i++){
		s[i] = rand() % 2 ? '1' : '0';
	}


	ADD ans, V;
	if (s[0] == '1'){
		ans = exchange_matrix(z_vars[0], y_vars[0]);
		V = exchange_matrix(z_vars[0], x_vars[0]);
	}else{
		ans = identity_matrix(z_vars[0], y_vars[0]);
		V = identity_matrix(z_vars[0], x_vars[0]);
	}

	for (unsigned int i = 1; i < N; i++){
		ADD tmp1, tmp2;
		if (s[i] == '1'){
			tmp1 = exchange_matrix(z_vars[i], y_vars[i]);
			tmp2 = exchange_matrix(z_vars[i], x_vars[i]);
		}else{
			tmp1 = identity_matrix(z_vars[i], y_vars[i]);
			tmp2 = identity_matrix(z_vars[i], x_vars[i]);
		}
		ans *= tmp1;
		V *= tmp2;

	}

	ans *= identity_matrix(z_vars[N], y_vars[N]);
	V *= identity_matrix(z_vars[N], x_vars[N]);

	ADD I = identity_n(mgr, 0, N, x_vars, y_vars);
  	ADD I_shift = identity_n(mgr, 0, N, y_vars, z_vars);

	ADD C = I * CNOT_matrix(x_vars[0], y_vars[0], x_vars[N], y_vars[N]);
	for (unsigned int i = 1; i < N; i++){
		ADD tmp = I_shift * CNOT_matrix(y_vars[i], z_vars[i], y_vars[N], z_vars[N]);
		C = C.MatrixMultiply(tmp, y_vars);
		C = C.SwapVariables(y_vars, z_vars);
	}
	
	C = C.SwapVariables(y_vars, z_vars);
	ans = C.MatrixMultiply(ans, z_vars);
	ans = V.MatrixMultiply(ans, x_vars);
	ans = ans.SwapVariables(x_vars, z_vars);

	return ans;

}

ADD CreateFMatrix_DJ(Cudd& mgr, unsigned int N, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, std::vector<ADD>& z_vars){
	// srand(time(NULL));
	srand(2);
	bool isConstant_F = (rand() % 2);
	if (isConstant_F){
		int value = (rand() % 2);
		CUDD_VALUE_TYPE val;
		mpfr_init_set_si(val.t_val, value, RND_TYPE);
		ADD F = mgr.constant(val);
		mpfr_clear(val.t_val);
		return F;
	} else{
		// srand(time(NULL));
		return CreateBalancedFn(mgr, N, x_vars, y_vars, z_vars);
	}
}

unsigned int DJ(Cudd &mgr, int n){

	adjustPrecision(mgr, n, 6);

	unsigned int N = (unsigned int)pow(2, n);
  	std::vector<ADD> x_vars, y_vars, z_vars;
	for (unsigned int i = 0; i < pow(2, n)+1; i++){
	  x_vars.push_back(mgr.addVar(3*i));
	  y_vars.push_back(mgr.addVar(3*i+1));
	  z_vars.push_back(mgr.addVar(3*i+2));
	}

	ADD F = CreateFMatrix_DJ(mgr, N, x_vars, y_vars, z_vars);
	std::cout << "F nodeCount: " <<  F.nodeCount() << std::endl;
	high_resolution_clock::time_point start = high_resolution_clock::now();
	ADD ans = ~y_vars[N] - y_vars[N];
	ans = F.MatrixMultiply(ans, y_vars);
	ans = ans.SwapVariables(x_vars, y_vars);
	ADD H = hadamard_n(mgr, 0, N, x_vars, y_vars);
	ADD HI = H * identity_matrix(x_vars[N], y_vars[N]);
	ans = HI.MatrixMultiply(ans, y_vars);
  	high_resolution_clock::time_point end = high_resolution_clock::now();
  	bool isOutCorrect = false;
  	if (F.IsOne() || F.IsZero()){
  		isOutCorrect = ans.IsZero();
  	}
  	else{
  		isOutCorrect = !ans.IsZero();
  	}
  	duration<double> time_taken = duration_cast<duration<double>>(end - start);
  	std::cout << "Is Output correct: " << isOutCorrect << std::endl;
  	std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
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

	if (strcmp(argv[1], "xor") == 0)
      nodeCount = exclusive_or(mgr, atoi(argv[2])); 
    else if (strcmp(argv[1], "hi_sum") == 0)
      nodeCount = hi_sum(mgr, atoi(argv[2]));
  	else if (strcmp(argv[1], "simons") == 0)
      nodeCount = simons(mgr, atoi(argv[2]));
	else if (strcmp(argv[1], "grover") == 0)
      nodeCount = grover(mgr, atoi(argv[2]));
  	else if (strcmp(argv[1], "GHZ") == 0)
      nodeCount = GHZ(mgr, atoi(argv[2]));
  	else if (strcmp(argv[1], "BV") == 0)
      nodeCount = BV(mgr, atoi(argv[2]));
  	else if (strcmp(argv[1], "DJ") == 0)
      nodeCount = DJ(mgr, atoi(argv[2]));


	return 0;
}