#include "../cudd-3.0.0/cplusplus/cuddObj.hh"
#include <iostream>
#include <cmath>
#include <ctime>
#include <ratio>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <typeinfo>
#include <complex>
using namespace std::chrono;

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

ADD identity_matrix(ADD x, ADD y){
  return (~x * ~y) + (x * y);
}

ADD hadamard_matrix(ADD x, ADD y){
  return ~x + x * (~y - y);
}
ADD CNOT_matrix(ADD x, ADD y, ADD w, ADD z){
  return (~x * ~y * ~z * ~w) + (~x * w * ~y * z) + (x * ~w * y * z) + (x * w * y * ~z);
}
ADD addition_relation(ADD x, ADD y, ADD z){
  return (~z * ~(x + y)) + (z * (x + y));
}
ADD identity_n(Cudd &mgr, unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 1)
    return identity_matrix(x_vars[start], y_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return identity_n(mgr, start, mid, x_vars, y_vars) * identity_n(mgr, mid, end, x_vars, y_vars);
}

ADD hadamard_n(Cudd &mgr,unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars) {
  if (start == end - 1)
    return hadamard_matrix(x_vars[start], y_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return hadamard_n(mgr, start, mid, x_vars, y_vars) * hadamard_n(mgr, mid, end, x_vars, y_vars);
}
ADD C_n(Cudd &mgr,unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars) {
  if (start == end - 1)
    return CNOT_matrix(x_vars[start], y_vars[start], w_vars[start], z_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return C_n(mgr, start, mid, x_vars, y_vars, w_vars, z_vars) * C_n(mgr, mid, end, x_vars, y_vars, w_vars, z_vars);
}

ADD addition_n(Cudd &mgr,unsigned int start, unsigned int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, std::vector<ADD>& z_vars) {
  if (start == end - 1)
    return addition_relation(x_vars[start], y_vars[start], z_vars[start]);
  unsigned int mid = (end - start)/2 + start;
  return addition_n(mgr, start, mid, x_vars, y_vars, z_vars) * 
  addition_n(mgr, mid, end, x_vars, y_vars, z_vars);
}

unsigned int addition(Cudd &mgr, int n){
  std::vector<BDD> x_vars, y_vars, z_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.bddVar(3*i));
    y_vars.push_back(mgr.bddVar(3*i + 1));
    z_vars.push_back(mgr.bddVar(3*i + 2));
  }
  high_resolution_clock::time_point start = high_resolution_clock::now();
  //ADD A = addition_n(mgr, 0, pow(2, n), x_vars, y_vars, z_vars);
  BDD A = mgr.Dxyeqz(x_vars, y_vars, z_vars);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " << A.nodeCount() << " time: " << time_taken.count() << std::endl;
  return A.nodeCount();
}



unsigned int hi_sum(Cudd &mgr, int n){

/*
  ADD x = mgr.addVar(0);
  ADD y = mgr.addVar(1);
  ADD I = identity_matrix(x, y);
  ADD H = hadamard_matrix(x, y);

  for (unsigned int i = 1; i < pow(2, n); i++){
    ADD var1 = mgr.addVar(2 * i);
    ADD var2 = mgr.addVar(2 * i + 1);
    I *= identity_matrix(var1, var2);
    H *= hadamard_matrix(var1, var2);
  }
  */
  std::vector<ADD> x_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(2*i));
    y_vars.push_back(mgr.addVar(2*i + 1));
  }
  unsigned int N = (unsigned int)pow(2, n);
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD H = mgr.Walsh(x_vars, y_vars);
  ADD I = mgr.Xeqy(x_vars, y_vars);
  std::cout << "Matrices created .. " << std::endl;
  ADD ans = H + I;
  for (unsigned int i = 1; i < 1024; i++){
    ans += H + I;
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "time_taken: " << time_taken.count() << " H nodeCount: " << H.nodeCount() << " ans nodeCount: " <<  ans.nodeCount() << std::endl;
  return ans.nodeCount();
}

unsigned int matmult(Cudd &mgr, int n){

  std::vector<ADD> x_vars, y_vars, z_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(3*i));
    y_vars.push_back(mgr.addVar(3*i + 1));
    z_vars.push_back(mgr.addVar(3*i + 2));
  }
  unsigned int N = (unsigned int)pow(2, n);
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD H = mgr.Walsh(x_vars, y_vars);
  ADD I = mgr.Xeqy(x_vars, y_vars);
  ADD X = mgr.Xneqy(x_vars, y_vars);
  ADD ans = mgr.addZero();
  for (unsigned int i = 1; i < 4; i++){
	  if (i % 3 == 0){
		 ADD tmpH = H.SwapVariables(y_vars, z_vars); 
		 ADD tmpI = I.SwapVariables(x_vars, z_vars); 
		 ans = ans + tmpH.MatrixMultiply(tmpI, z_vars);
	  }
	  if (i % 3 == 1){
		 ADD tmpX = X.SwapVariables(y_vars, z_vars); 
		 ADD tmpH = H.SwapVariables(x_vars, z_vars); 
		 ans = ans + tmpX.MatrixMultiply(tmpH, z_vars);
	  }
	  if (i % 3 == 2){
		 ADD tmpI = I.SwapVariables(y_vars, z_vars); 
		 ADD tmpX = X.SwapVariables(x_vars, z_vars); 
		 ans = ans + tmpI.MatrixMultiply(tmpX, z_vars);
	  }
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "time_taken: " << time_taken.count() << " H nodeCount: " << H.nodeCount() << " ans nodeCount: " <<  ans.nodeCount() << std::endl;
  return ans.nodeCount();

}

unsigned int mult(Cudd &mgr, int n){

  high_resolution_clock::time_point start = high_resolution_clock::now();
  std::vector<ADD> x_vars, y_vars, z_vars;
  for (unsigned int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(3*i));
    y_vars.push_back(mgr.addVar(3*i + 1));
    z_vars.push_back(mgr.addVar(3*i + 2));
  }

  ADD H = hadamard_n(mgr, 0, pow(2,n), x_vars, y_vars);
  ADD I = identity_n(mgr, 0, pow(2,n), y_vars, z_vars);
  ADD ans = H.MatrixMultiply(I, y_vars);
  ans = ans.SwapVariables(y_vars, z_vars);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  std::cout << (H == ans) << std::endl;
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << time_taken.count() << " " << H.nodeCount() << " " <<  ans.nodeCount() << std::endl;
  return ans.nodeCount();
}

/*
ADD Func2To1Matrix(std::string s, std::vector<int> v, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars, int w_start, int w_end, int z_start, int z_end){
  unsigned int numOfBits = s.length();
  unsigned int sizeOfMatrix = pow(2, numOfBits);

  std::unordered_map<std::string, unsigned int> halfMatrixValuesMap;
  std::vector<ADD> baseMatrices = formMatrixWithOneNonZeroValueUtil();
  int v_index = 0;
  unsigned int random_value;
  if (v.size() != 0)
    random_value = v[v_index++];
  else
    random_value = rand() % sizeOfMatrix;

  halfMatrixValuesMap.insert(std::make_pair(to_string(compute_xor(0, s)), random_value));
  ADD Func2To1 = formMatrixWithOneZeroValue(baseMatrices, 0, random_value, sizeOfMatrix);

  for (unsigned int i = 1; i < sizeOfMatrix; i++){
    std::string strIndex = to_string(i);
    if (halfMatrixValuesMap.find(strIndex) == halfMatrixValuesMap.end()){
      unsigned int random_value;
      if (v.size() != 0)
	random_value = v[v_index++];
      else
	random_value = rand() % sizeOfMatrix;
      halfMatrixValuesMap.insert(std::make_pair(to_string(compute_xor(i, s)), random_value));
      Func2To1 += formMatrixWithOneNonZeroValue(baseMatrices, i, random_value, sizeOfMatrix);
    }else{
      Func2To1 += formMatrixWithOneNonZeroValue(baseMatrices, i, halfMatrixValuesMap[strIndex], sizeOfMatrix);
    }
  }
  return Func2To1;
}


ADD Func2To1Matrix_DivideAndConquer(std::string s, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars, int w_start, int w_end, int z_start, int z_end){

  if (s.length() == 2){
    if (s == "00"){
      std:vector<int> v;
      for (int i = 0; i < 4; i++)
	v.push_back(rand() % 4);
      return Func2To1Matrix(s,v, w_vars, z_vars, w_start, w_end, z_start, z_end);
    }else{
      std:vector<int> v;
      for (int i = 0; i < 2; i++)
	v.push_back(rand() % 4);
      return Func2To1Matrix(s,v);
    }
  }
  int w_mid = (w_end - w_start)/2 + w_start;
  int z_mid = (z_end - z_start)/2 + z_start;
  ADD F1 = Func2To1Matrix_DivideAndConquer(s.substr(0, s.length()/2), w_vars, z_vars, w_start, w_mid, z_start, z_mid);
  ADD F2 = Func2To1Matrix_DivideAndConquer(s.substr(s.length()/2), w_vars, z_vars, w_mid, w_end, z_mid, z_end);

  ADD ans = Create2To1Func(F1, s.substr(0, s.length()/2), F2, s.substr(s.length()/2));
  return ans;

  return (~w_vars[0] * ~w_vars[1] * ~z_vars[0] * ~z_vars[1])
  + (~w_vars[0] * w_vars[1] * ~z_vars[0] * z_vars[1])
  + (w_vars[0] * ~w_vars[1] * ~z_vars[0] * ~z_vars[1])
  + (w_vars[0] * w_vars[1] * ~z_vars[0] * z_vars[1]);
  
}
*/

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

ADD CreateFnMatrix(Cudd &mgr, std::string s, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars){
  int n = w_vars.size();
  unsigned long long int numVals = pow(2,n);
  ADD ans = mgr.constant(0);
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
    ADD tmp_ans = mgr.constant(1);
    for (int j = 0; j < n; j++){
      if (inp_string[j] == '1')
	tmp_ans *= w_vars[j];
      else
	tmp_ans *= ~w_vars[j];
    }
    for (int j = 0; j < n; j++){
      if (val_string[j] == '1')
	tmp_ans *= z_vars[j];
      else
	tmp_ans *= ~z_vars[j];
    }
    ans += tmp_ans;
  }

  return ans;
}

unsigned int simons(Cudd &mgr, int n){
  std::vector<ADD> x_vars, y_vars, w_vars, z_vars; // x & y belong to input and w & z belong to output
  for (int i = 0; i < n; i++){
    x_vars.push_back(mgr.addVar(4*i));
    y_vars.push_back(mgr.addVar(4*i+1));
    w_vars.push_back(mgr.addVar(4*i+2));
    z_vars.push_back(mgr.addVar(4*i+3));
  }
  std::string s = "1001100110011001";
  ADD U = CreateFnMatrix(mgr, s, w_vars, z_vars);
  U = U.SwapVariables(z_vars, w_vars);
  //ADD C = C_n(mgr, 0, pow(2, n), x_vars, y_vars, w_vars, z_vars);
  std::cout << "Algo start.." << std::endl;
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD C = CNOT_matrix(x_vars[0], y_vars[0], w_vars[0], z_vars[0]);
  for (int i = 1; i < n; i++){
    C *= CNOT_matrix(x_vars[i], y_vars[i], w_vars[i], z_vars[i]);
  }
  ADD e_x = (y_vars[0] + ~y_vars[0]);
  for (int i = 1; i < n; i++){
    e_x *= (y_vars[i] + ~y_vars[i]);
  }
  ADD e_0 = ~z_vars[0];
  for (int i = 1; i < n; i++){
    e_0 *= ~z_vars[i];
  }
  ADD tmp = e_x * e_0;
  std::vector<ADD> mult_array;
  mult_array.insert(mult_array.end(), y_vars.begin(), y_vars.end());
  mult_array.insert(mult_array.end(), z_vars.begin(), z_vars.end());
  std::vector<ADD> swap_array;
  swap_array.insert(swap_array.end(), x_vars.begin(), x_vars.end());
  swap_array.insert(swap_array.end(), w_vars.begin(), w_vars.end());
  tmp = C.MatrixMultiply(tmp, mult_array);
  tmp = tmp.SwapVariables(swap_array, mult_array);
  ADD H = hadamard_n(mgr, 0, n, x_vars, y_vars);
  ADD tmp1 = H * U;
  tmp1 = tmp1.MatrixMultiply(tmp, mult_array);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " << tmp1.nodeCount() << " time: " << time_taken.count() << std::endl;
  return tmp1.nodeCount();
}

ADD MkCyclicKMatrix(Cudd &mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){
  int n = x_vars.size();
  ADD ans = mgr.constant(0);
  for (int i = 0; i < pow(2, n)/2; i++){
    std::string inp_string = getBits(i, n); 
    std::string val_string = getBits(2*i, n);
    ADD tmp_ans = mgr.constant(1);
    for (int j = 0; j < n; j++){
      if (inp_string[j] == '1')
	     tmp_ans *= x_vars[j];
      else
	     tmp_ans *= ~x_vars[j];
    }
    for (int j = 0; j < n; j++){
      if (val_string[j] == '1')
	     tmp_ans *= y_vars[j];
      else
	     tmp_ans *= ~y_vars[j];
    }
    ans += tmp_ans;
  }

  for (int i = 0; i < pow(2, n)/2; i++){
    std::string inp_string = getBits(i + pow(2, n)/2, n); 
    std::string val_string = getBits(2*i+1, n);
    ADD tmp_ans = mgr.constant(1);
    for (int j = 0; j < n; j++){
      if (inp_string[j] == '1')
	     tmp_ans *= x_vars[j];
      else
	     tmp_ans *= ~x_vars[j];
    }
    for (int j = 0; j < n; j++){
      if (val_string[j] == '1')
	     tmp_ans *= y_vars[j];
      else
	     tmp_ans *= ~y_vars[j];
    }
    ans += tmp_ans;
  }
  return ans;
}

ADD D_n(Cudd &mgr, int start, int end, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, std::string comp) {
  int n = end - start;
  double const pi = 4 * std::atan(1);
  //double temp = static_cast<double>(n) * static_cast<double>(n);
  double temp = pow(2, n); 
  double temp2 = 1.0/temp;
  double w = static_cast<double>(-2.0) * pi * temp2;
  ADD ans = mgr.constant(0);
  for (int i = 0; i < pow(2, n); i++){
    std::string inp_string = getBits(i, n); 
    ADD tmp_ans = mgr.constant(1);
    for (int j = 0; j < n; j++){
      if (inp_string[j] == '1')
        tmp_ans *= x_vars[start + j];
      else
        tmp_ans *= ~x_vars[start + j];
    }
    for (int j = 0; j < n; j++){
      if (inp_string[j] == '1')
        tmp_ans *= y_vars[start + j];
      else
        tmp_ans *= ~y_vars[start + j];
    }
    double val = w * i;
    if (comp == "real")
      tmp_ans = tmp_ans * mgr.constant(cos(val));
    else if (comp == "img")
      tmp_ans = tmp_ans * mgr.constant(sin(val));
    else{
      //std::complex<double> w_pow = std::polar(1.0, val);
      double w_pow = val;
      tmp_ans = tmp_ans * mgr.constant(w_pow); 
    }
    ans += tmp_ans;
  }

  return ans;

}

std::pair<ADD,ADD> Fourier(Cudd &mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars, int start){
  int n = x_vars.size();
  if (n == start + 1){
    return std::make_pair(hadamard_matrix(x_vars[start], y_vars[start]) * mgr.constant(1.0/sqrt(2)), mgr.constant(0));
  }
  ADD K = MkCyclicKMatrix(mgr, z_vars, y_vars); 
  std::pair<ADD,ADD> F_recurse = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, start + 1);
  if (true){
    ADD F_real = (~x_vars[start] * ~y_vars[start] + x_vars[start] * y_vars[start]) * F_recurse.first;
    ADD F_img = (~x_vars[start] * ~y_vars[start] + x_vars[start] * y_vars[start]) * F_recurse.second;
    F_real = F_real.SwapVariables(x_vars, w_vars);
    F_img = F_img.SwapVariables(x_vars, w_vars);
    F_real = F_real.SwapVariables(y_vars, z_vars);
    F_img = F_img.SwapVariables(y_vars, z_vars);
    ADD I = identity_n(mgr, 1, n,  x_vars, y_vars);
    ADD D_real = D_n(mgr, 1, n, x_vars, y_vars, "real");
    ADD D_img = D_n(mgr, 1, n, x_vars, y_vars, "img");
    ADD ID_real = ~y_vars[start] * I + y_vars[start] * (~x_vars[start] - x_vars[start]) * D_real; 
    ADD ID_img = y_vars[start] * (~x_vars[start] - x_vars[start]) * D_img; 
    ID_real = ID_real.SwapVariables(y_vars, w_vars);
    ID_img = ID_img.SwapVariables(y_vars, w_vars);
    ADD real = ID_real.MatrixMultiply(F_real, w_vars) - ID_img.MatrixMultiply(F_img, w_vars);
    ADD img = ID_real.MatrixMultiply(F_img, w_vars) - ID_img.MatrixMultiply(F_real, w_vars);
    real = real.MatrixMultiply(K, z_vars) * mgr.constant(1.0/sqrt(2));
    img = img.MatrixMultiply(K, z_vars) * mgr.constant(1.0/sqrt(2));
    return std::make_pair(real, img);
  }else{
    ADD F = (~x_vars[start] * ~y_vars[start] + x_vars[start] * y_vars[start]) * F_recurse.first;
    F = F.SwapVariables(x_vars, w_vars);
    F = F.SwapVariables(y_vars, z_vars);
    ADD I = identity_n(mgr, 1, n,  x_vars, y_vars);
    ADD D = D_n(mgr, 1, n, x_vars, y_vars, "both");
    ADD ID = ~y_vars[start] * I + y_vars[start] * (~x_vars[start] - x_vars[start]) * D; 
    ID = ID.SwapVariables(y_vars, w_vars);
    ADD IDF = ID.MatrixMultiply(F, w_vars);
    ADD ans = IDF.MatrixMultiply(K, z_vars) * mgr.constant(1.0/sqrt(2));
    return std::make_pair(ans, mgr.constant(0));
  }
}

unsigned int MkFourierTransform(Cudd &mgr, int n){
  std::vector<ADD> x_vars, y_vars, w_vars, z_vars;
  for (int i = 0; i < pow(2, n); i++){
    x_vars.push_back(mgr.addVar(4*i));
    y_vars.push_back(mgr.addVar(4*i+1));
    w_vars.push_back(mgr.addVar(4*i+2));
    z_vars.push_back(mgr.addVar(4*i+3));
  }
  high_resolution_clock::time_point start = high_resolution_clock::now();
  std::pair<ADD,ADD> ans = Fourier(mgr, x_vars, y_vars, w_vars, z_vars, 0); 
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << ans.first.nodeCount() << " " << ans.second.nodeCount() << std::endl;
  std::cout << "time taken: " << time_taken.count() << std::endl;
  return ans.first.nodeCount() + ans.second.nodeCount();
}

ADD MkU_w(Cudd &mgr, std::string s, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){
  int n = x_vars.size();
  ADD I = identity_n(mgr, 0, n, x_vars, y_vars);
  ADD F = mgr.constant(-2);
  for (int i = 0; i < n; i++){
    if (s[i] == '0'){
      F *= ~x_vars[i] * ~y_vars[i];
    }else{
      F *= x_vars[i] * y_vars[i];
    }
  }
  return I + F;
}

ADD MkU_s(Cudd &mgr, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){//, std::vector<ADD>& w_vars, std::vector<ADD>& z_vars){
  int n = x_vars.size();
  ADD I = identity_n(mgr, 0, n, x_vars, y_vars);
  /*
  ADD tmp = mgr.constant(2);
  for (int i = 0; i < n; i++){
    tmp *= ~w_vars[i] * ~z_vars[i];
  }
  I =  tmp - mgr.constant(1)*I;
  ADD H1 = hadamard_n(mgr, 0, n, x_vars, w_vars) * mgr.constant(1.0/sqrt(pow(2,n)));
  ADD H2 = hadamard_n(mgr, 0, n, z_vars, y_vars) * mgr.constant(1.0/sqrt(pow(2,n)));
  ADD ans = I.MatrixMultiply(H2, z_vars);
  ans = H1.MatrixMultiply(ans, w_vars);
  return ans;
  */
  return (mgr.constant(2.0/pow(2,n)) - I);// * mgr.constant(1.0/pow(2,n));
}

std::pair<ADD, ADD> MultiplyRec(Cudd& mgr, ADD M, unsigned long long int iters, std::vector<ADD>& x_vars, 
  std::vector<ADD>& y_vars, std::vector<ADD>& z_vars, 
  std::unordered_map<unsigned long long int, ADD>& memo,
  unsigned int n,
  unsigned int* level){
  auto it = memo.find(iters);
  if (it != memo.end())
    return std::make_pair(it->second, it->second);
  if (iters == 1){
    if ((*level) >= n)
      memo[iters] = M;
    if ((*level) < n){
      *level = *level + 1;
      return std::make_pair(mgr.constant(1.0/sqrt(2)) * M, M);
    }
    return std::make_pair(M, M);
  }
  unsigned long long int half_iters = iters/2;
  auto X = MultiplyRec(mgr, M, half_iters, x_vars, y_vars, z_vars, memo, n, level);
  auto Y = MultiplyRec(mgr, M, iters - half_iters, x_vars, y_vars, z_vars, memo, n, level);
  ADD M1 = X.first;
  ADD M2 = Y.first;
  M1 = M1.SwapVariables(y_vars, z_vars);
  M2 = M2.SwapVariables(z_vars, x_vars);
  ADD ans = M1.MatrixMultiply(M2, z_vars);

  ADD M1P = X.second;
  ADD M2P = Y.second;
  M1P = M1P.SwapVariables(y_vars, z_vars);
  M2P = M2P.SwapVariables(z_vars, x_vars);
  ADD ansP = M1P.MatrixMultiply(M2P, z_vars);

  if ((*level) >= n)
    memo[iters] = ansP;
  return std::make_pair(ans, ansP);

}

unsigned int grover(Cudd &mgr, int n){
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
  M.print(4*n,2);
  std::cout << "M leaf counts: " << M.CountLeaves() << std::endl;
  M = M.SwapVariables(w_vars, z_vars);
  double const pi = 4 * std::atan(1);
  unsigned long long int iters = (unsigned long long int)floor(pi * 0.25 * pow(2, N/2)); 
  ADD ans = mgr.constant(1.0);
  std::cout << "M count: " << M.nodeCount() << std::endl;
  std::cout << "iters count: " << iters << std::endl;
  std::unordered_map<unsigned long long int, ADD> memo;
  unsigned int ulevel = 0;
  if (N > iters){
    ans = ans * mgr.constant(1.0/pow(2,(N-iters)/2));
  }
  auto M_rec = MultiplyRec(mgr, M, iters, x_vars, w_vars, y_vars, memo, N, &ulevel);
  ADD M_actual = M_rec.first;
  // M_actual.print(4*N, 2);
  // M_rec.second.print(4*N, 2);
  // ans.print(4*N, 2);
  ans = M_actual.MatrixMultiply(ans, w_vars);
  // for (unsigned long long int i = 0; i < iters; i++){
  //   ans = M.MatrixMultiply(ans, w_vars);// * mgr.constant(1.0/pow(2,n));
  //   // ans.print(4*n, 2);
  //   std::cout << "i : " << i << " ans count: " << ans.nodeCount() << " M leaves: " << M.CountLeaves() << " ans leaves: " << ans.CountLeaves() << std::endl;
  //   //std::cout << "i : " << i << std::endl;
  //   ans = ans.SwapVariables(x_vars, w_vars);
  // }
  // // ans.print(4*n, 2);
  // ans = ans.SwapVariables(x_vars, w_vars);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "string: " << s << std::endl;
  std::cout << "nodeCount: " <<  ans.nodeCount() << " leaf count: " << ans.CountLeaves() << " time_taken: " << time_taken.count() << std::endl;
  ans.print(4*n,2);
  return ans.nodeCount();
}

unsigned int GHZ(Cudd &mgr, unsigned int n){
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
  ADD H_0 = hadamard_matrix(x_vars[N-1], y_vars[N-1]);
  H = H * H_0;
  ans = H.MatrixMultiply(ans, y_vars);
  ans.print(4*N,2);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> time_taken = duration_cast<duration<double>>(end - start);
  std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
  return ans.nodeCount();
}

unsigned int BV(Cudd &mgr, int n){
  unsigned int N = (unsigned int)pow(2, n);
  std::vector<ADD> x_vars, w_vars, z_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n)+1; i++){
    x_vars.push_back(mgr.addVar(3*i));
    y_vars.push_back(mgr.addVar(3*i+1));
    z_vars.push_back(mgr.addVar(3*i+2));
    /*
    y_vars.push_back(mgr.addVar(4*i+3));
    */
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
  ADD I = identity_n(mgr, 0, N, x_vars, y_vars);
  ADD I_shift = identity_n(mgr, 0, N, y_vars, z_vars);
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
  U.print(4*n,2);
  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD ans =  y_vars[N];
  for (unsigned int i = 0; i < N; i++)
    ans = ans * ~y_vars[i];
  ADD H = hadamard_n(mgr, 0, N, x_vars, y_vars);
  ADD I_0 = identity_matrix(x_vars[N], y_vars[N]);
  H = H * I_0 ;//* mgr.constant(1.0/pow(2,N));
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
  return ans.nodeCount();
}


ADD CreateFMatrix_DJ(Cudd& mgr, unsigned int N, std::vector<ADD>& x_vars, std::vector<ADD>& y_vars){
  srand (2);
  bool isConstant_F = (rand() % 2);
  if (isConstant_F){
    int value = (rand() % 2);
    ADD F = mgr.constant(value);
    return F;
  } else{
    unsigned int count = 0;
    ADD ans = mgr.constant(0);
    for (unsigned int i = 0; i < pow(2, N); i++){
      int f_x = (rand() % 2);
      if (count >= pow(2, N)/2)
        f_x = 0;
      if (f_x == 1)
        count++;
      std::cout << i << " " << f_x << std::endl;
      std::string x_string = getBits(i, N);
      ADD tmp_ans = mgr.constant(1);
      for (unsigned int j = 0; j < N; j++){
        if (x_string[j] == '1')
          tmp_ans *= x_vars[j] * y_vars[j];
        else
          tmp_ans *= ~x_vars[j] * ~y_vars[j];
      }

      if (f_x == 0){
        ans += (tmp_ans * ~x_vars[N] * ~y_vars[N]) + (tmp_ans * x_vars[N] * y_vars[N]);
      }
      else{
        ans += (tmp_ans * ~x_vars[N] * y_vars[N]) + (tmp_ans * x_vars[N] * ~y_vars[N]);       
      }

    }
    return ans;
  }
}

unsigned int DJ(Cudd &mgr, int n){
  unsigned int N = (unsigned int)pow(2, n);
    std::vector<ADD> x_vars, y_vars;
  for (unsigned int i = 0; i < pow(2, n)+1; i++){
    x_vars.push_back(mgr.addVar(2*i));
    y_vars.push_back(mgr.addVar(2*i+1));
  }
  ADD F = CreateFMatrix_DJ(mgr, N, x_vars, y_vars);
  F.print(2*N, 2);

  high_resolution_clock::time_point start = high_resolution_clock::now();
  ADD ans = ~y_vars[N] - y_vars[N];
  ans = F.MatrixMultiply(ans, y_vars);
  ans = ans.SwapVariables(x_vars, y_vars);
  ADD H = hadamard_n(mgr, 0, N, x_vars, y_vars);
  ADD HI = H * identity_matrix(x_vars[N], y_vars[N]);
  ans = HI.MatrixMultiply(ans, y_vars);
    high_resolution_clock::time_point end = high_resolution_clock::now();
    ans.print(2*n,2);
    duration<double> time_taken = duration_cast<duration<double>>(end - start);
    std::cout << "nodeCount: " <<  ans.nodeCount() << " time_taken: " << time_taken.count() << std::endl;
    return ans.nodeCount();
}


int main (int argc, char** argv)
{
    if (argc < 2)
      return 0;
    Cudd mgr(0,0);
    //mgr.AutodynEnable();
    unsigned int nodeCount = 0;
    if (strcmp(argv[1], "xor") == 0)
      nodeCount = exclusive_or(mgr, atoi(argv[2])); 
    else if (strcmp(argv[1], "hi_sum") == 0)
      nodeCount = hi_sum(mgr, atoi(argv[2])); 
    else if (strcmp(argv[1], "matmult") == 0)
      nodeCount = matmult(mgr, atoi(argv[2])); 
    else if (strcmp(argv[1], "mult") == 0)
      nodeCount = mult(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "simons") == 0)
      nodeCount = simons(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "fourier") == 0)
      nodeCount = MkFourierTransform(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "grover") == 0)
      nodeCount = grover(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "BV") == 0)
      nodeCount = BV(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "GHZ") == 0)
      nodeCount = GHZ(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "DJ") == 0)
      nodeCount = DJ(mgr, atoi(argv[2]));
    else if (strcmp(argv[1], "add") == 0)
      nodeCount = addition(mgr, atoi(argv[2]));
    return 0;
}
