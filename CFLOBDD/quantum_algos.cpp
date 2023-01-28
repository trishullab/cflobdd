
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <ctime>
#include <time.h>
#include <chrono>
#include "cflobdd_int.h"
#include "cflobdd_node.h"
#include "cflobdd_top_node_t.h"
#include "cflobdd_top_node_int.h"
#include "matrix1234_node.h"
#include "matrix1234_double_top_node.h"
#include "matrix1234_double.h"
#include "vector_double.h"
#include "matrix1234_float_boost_top_node.h"
#include "matrix1234_float_boost.h"
#include "vector_float_boost.h"
#include "vector_complex_float_boost.h"
#include "matrix1234_complex_float_boost.h"
#include "matrix1234_fourier.h"
#include "Solver/uwr/matrix/HowellMatrix.h"
#include "Solver/uwr/matrix/ModularSquareMatrix.h"

using namespace std::chrono;

namespace CFL_OBDD {

	namespace QuantumAlgos {

		void QuantumAlgosInitializer()
		{
			return;
		}
	
		CFLOBDD_DOUBLE DeutschJozsaAlgo(unsigned int n, CFLOBDD_DOUBLE F)
		{
			unsigned int N = 1 << n;
			unsigned int level = ceil(log2(n));
			CFLOBDD_DOUBLE e0n = VectorDouble::MkBasisVector(level,0);
			CFLOBDD_DOUBLE e1 = VectorDouble::MkBasisVector(0, 1);
			if (e0n.root->level > e1.root->level)
			{
				e1 = VectorDouble::VectorPadWithZeros(e1, e0n.root->level);
			}
			
			unsigned int numtemp = n;
			unsigned int Walshlevel = floor(log2(numtemp)) + 1;
			CFLOBDD_DOUBLE H = Matrix1234Double::MkWalshInterleaved(Walshlevel);
			numtemp = numtemp - (1 << (Walshlevel-1));
			while (numtemp > 0)
			{
				Walshlevel = floor(log2(numtemp)) + 1;
				CFLOBDD_DOUBLE Htemp = Matrix1234Double::MkWalshInterleaved(Walshlevel);
				Htemp = Matrix1234Double::MatrixPadWithZeros(Htemp, H.root->level);
				//Htemp = Matrix1234Double::PromoteInterleavedTo12(Htemp);
				//H = Matrix1234Double::PromoteInterleavedTo12(H);
				H = Matrix1234Double::KroneckerProduct2Vocs(H, Htemp);
				numtemp = numtemp - (1 << (Walshlevel - 1));
			}
			
			CFLOBDD_DOUBLE Htemp = Matrix1234Double::MkWalshInterleaved(Walshlevel);
			if (Htemp.root->level != H.root->level)
			{
				Htemp = Matrix1234Double::MatrixPadWithZeros(Htemp, H.root->level);
			}
			Htemp = Matrix1234Double::PromoteInterleavedTo12(Htemp);
			CFLOBDD_DOUBLE HH = Matrix1234Double::PromoteInterleavedTo12(H);
			HH = Matrix1234Double::KroneckerProduct(HH, Htemp);

			CFLOBDD_DOUBLE I = Matrix1234Double::MkIdRelationInterleaved(1);
			if (I.root->level != H.root->level)
			{
				I = Matrix1234Double::MatrixPadWithZeros(I, H.root->level);
			}
			I = Matrix1234Double::PromoteInterleavedTo12(I);
			CFLOBDD_DOUBLE HI = Matrix1234Double::PromoteInterleavedTo12(H);
			HI = Matrix1234Double::KroneckerProduct(HI, I);

			assert(HI.root->level == HH.root->level);
			assert(HI.root->level == F.root->level);

			//if (F.root->level != HH.root->level)
			//{
			//	F = Matrix1234Int::MatrixPadWithZeros(F, HH.root->level);
			//}

			HI = Matrix1234Double::PromoteInterleavedTo12(HI);
			F = Matrix1234Double::PromoteInterleavedTo12(F);
			HH = Matrix1234Double::PromoteInterleavedTo12(HH);
			
			CFLOBDD_DOUBLE HIF = Matrix1234Double::MatrixMultiplyV4(HI, F);
			
			CFLOBDD_DOUBLE v = Matrix1234Double::MatrixMultiplyV4(HIF, HH);
			e0n = VectorDouble::MkVectorWithVoc12(e0n);
			e1 = VectorDouble::MkVectorWithVoc12(e1);
			CFLOBDD_DOUBLE vectorTensor = VectorDouble::KroneckerProduct(e0n,e1);
			vectorTensor = VectorDouble::VectorToMatrixInterleaved(vectorTensor);
			vectorTensor = Matrix1234Double::PromoteInterleavedTo12(vectorTensor);
			CFLOBDD_DOUBLE ans = Matrix1234Double::MatrixMultiplyV4(v, vectorTensor);
			ans = Matrix1234Double::Demote12ToInterleaved(v);
			ans = VectorDouble::MatrixToVector(v);
			return ans;
		}

		CFLOBDD_DOUBLE KroneckerPower(CFLOBDD_DOUBLE(*f)(unsigned int), unsigned int n, unsigned int basicLevel)
		{
			unsigned int numtemp = n;
			unsigned int Walshlevel = floor(log2(numtemp)) + 1;
			CFLOBDD_DOUBLE H = f(Walshlevel + basicLevel-1);
			numtemp = numtemp - (1 << (Walshlevel - 1));
			while (numtemp > 0)
			{
				Walshlevel = floor(log2(numtemp)) + 1;
				CFLOBDD_DOUBLE Htemp = f(Walshlevel + basicLevel-1);
				Htemp = Matrix1234Double::MatrixPadWithZeros(Htemp, H.root->level);
				Htemp = Matrix1234Double::PromoteInterleavedTo12(Htemp);
				H = Matrix1234Double::PromoteInterleavedTo12(H);
				H = Matrix1234Double::KroneckerProduct(H, Htemp);
				numtemp = numtemp - (1 << (Walshlevel - 1));
			}

			return H;
		}

		CFLOBDD_DOUBLE CNOTk(unsigned int k, unsigned int n)
		{
			CFLOBDD_DOUBLE Ia = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, 2 * n - 1 - k, 1);
			CFLOBDD_DOUBLE A = VectorDouble::MkBasisVector(Ia.root->level - 1, 0);
			A = VectorDouble::VectorToMatrixInterleaved(A);
			A = Matrix1234Double::MatrixPadWithZeros(A, Ia.root->level);
			A = Matrix1234Double::PromoteInterleavedTo12(A);
			Ia = Matrix1234Double::PromoteInterleavedTo12(Ia);
			CFLOBDD_DOUBLE AI = Matrix1234Double::KroneckerProduct(A, Ia);

			if (k != 0)
			{
				CFLOBDD_DOUBLE temp = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, k, 1);
				temp = Matrix1234Double::PromoteInterleavedTo12(temp);
				AI = Matrix1234Double::PromoteInterleavedTo12(AI);
				AI = Matrix1234Double::KroneckerProduct(temp, AI);
			}

			CFLOBDD_DOUBLE Ib = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, n - 1, 1);
			CFLOBDD_DOUBLE Ibk = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, n - 1 - k, 1);
			CFLOBDD_DOUBLE B = VectorDouble::MkBasisVector(0, 1);
			B = VectorDouble::VectorToMatrixInterleaved(B);
			B = Matrix1234Double::ReverseColumns(B);
			B = Matrix1234Double::MatrixPadWithZeros(B, Ib.root->level);
			B = Matrix1234Double::PromoteInterleavedTo12(B);
			CFLOBDD_DOUBLE X = Matrix1234Double::MkNegationMatrixInterleaved(1);
			X = Matrix1234Double::MatrixPadWithZeros(X, Ibk.root->level);
			X = Matrix1234Double::PromoteInterleavedTo12(X);
			Ib = Matrix1234Double::PromoteInterleavedTo12(Ib);
			CFLOBDD_DOUBLE Btemp = Matrix1234Double::KroneckerProduct(B, Ib);
			Ibk = Matrix1234Double::PromoteInterleavedTo12(Ibk);
			CFLOBDD_DOUBLE Xtemp = Matrix1234Double::KroneckerProduct(X, Ibk);
			Btemp = Matrix1234Double::PromoteInterleavedTo12(Btemp);
			Xtemp = Matrix1234Double::PromoteInterleavedTo12(Xtemp);
			CFLOBDD_DOUBLE BI = Matrix1234Double::KroneckerProduct(Btemp, Xtemp);

			if (k != 0)
			{
				CFLOBDD_DOUBLE temp = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, k, 1);
				temp = Matrix1234Double::PromoteInterleavedTo12(temp);
				BI = Matrix1234Double::PromoteInterleavedTo12(BI);
				BI = Matrix1234Double::KroneckerProduct(temp, BI);
			}

			if (AI.root->level < BI.root->level){
				AI = Matrix1234Double::MatrixPadWithZeros(AI, BI.root->level);
			}
			else if (AI.root->level > BI.root->level){
				BI = Matrix1234Double::MatrixPadWithZeros(BI, AI.root->level);
			}
			CFLOBDD_DOUBLE C = AI + BI;
			return C;
		}

		CFLOBDD_DOUBLE CreateFuncMatrix(CFLOBDD_DOUBLE F, int n)
		{
			CFLOBDD_DOUBLE C = KroneckerPower(Matrix1234Double::MkCNOTInterleaved, n, 2);
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_DOUBLE I = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, n, 1);
			F = Matrix1234Double::PromoteInterleavedTo12(F);
			I = Matrix1234Double::PromoteInterleavedTo12(I);
			CFLOBDD_DOUBLE FI = Matrix1234Double::KroneckerProduct(F, I);

			F = Matrix1234Double::Demote12ToInterleaved(F);
			CFLOBDD_DOUBLE Ft = Matrix1234Double::MatrixTranspose(F);
			Ft = Matrix1234Double::PromoteInterleavedTo12(Ft);
			CFLOBDD_DOUBLE FtI = Matrix1234Double::KroneckerProduct(Ft, I);
			if (FtI.root->level < C.root->level)
				FtI = Matrix1234Double::MatrixPadWithZeros(FtI, C.root->level);
			
			FtI = Matrix1234Double::PromoteInterleavedTo12(FtI);
			C = Matrix1234Double::PromoteInterleavedTo12(C);
			CFLOBDD_DOUBLE temp = Matrix1234Double::MatrixMultiplyV4(FtI, C);
			std::cout << "FtI * C done .." << std::endl;
			if (FI.root->level < C.root->level - 1)
				FI = Matrix1234Double::MatrixPadWithZeros(FI, C.root->level - 1);
			FI = Matrix1234Double::PromoteInterleavedTo12(FI);
			CFLOBDD_DOUBLE ans = Matrix1234Double::MatrixMultiplyV4(temp, FI);
			std::cout << "temp * FI done.." << std::endl;
			ans = Matrix1234Double::Demote12ToInterleaved(ans);
			return ans;
		}

		void addEquations(HowellMatrix::HowellMatrix<BitVector::BV1>* h, std::string s)
		{
			auto zero = BitVector::BV1::zero();
			auto one = BitVector::BV1::one();

			HowellMatrix::SparseTArray<BitVector::BV1> *array = new HowellMatrix::SparseTArray<BitVector::BV1>(s.length() + 1);
			array->set(0, zero);
			for (unsigned i = 0; i < s.length(); i++)
			{
				if (s[i] == '0')
					array->set(i + 1, zero);
				else if (s[i] == '1')
					array->set(i + 1, one);
				else
					abort();
			}
			h->insertRow(array, true);
		}

		std::vector<std::string> SimonsAlgo(int n, CFLOBDD_DOUBLE F)
		{
			CFLOBDD_DOUBLE UFunc = CreateFuncMatrix(F, n);
			//std::cout << "UFunc created ..." << std::endl;
			unsigned int level = ceil(log2(n));
			CFLOBDD_DOUBLE SumEx = VectorDouble::NoDistinctionNode(level);
			CFLOBDD_DOUBLE e0 = VectorDouble::MkBasisVector(level, 0);
			/*
			SumEx = VectorDouble::MkVectorWithVoc12(SumEx);
			e0 = VectorDouble::MkVectorWithVoc12(e0);
			*/
			SumEx = VectorDouble::VectorToMatrixInterleaved(SumEx);
			e0 = VectorDouble::VectorToMatrixInterleaved(e0);
			SumEx = Matrix1234Double::PromoteInterleavedTo12(SumEx);
			e0 = Matrix1234Double::PromoteInterleavedTo12(e0);
			CFLOBDD_DOUBLE vectorTensor = Matrix1234Double::KroneckerProduct(SumEx, e0);

			//CFLOBDD_DOUBLE vectorTensor = VectorDouble::KroneckerProduct(SumEx, e0);
			//vectorTensor = VectorDouble::VectorToMatrixInterleaved(vectorTensor);
			//vectorTensor = (0.5 * pow(2, -n / 2)) * vectorTensor;
			vectorTensor = (1.0 * pow(2, -n / 4)) * vectorTensor;

			if (UFunc.root->level > vectorTensor.root->level){
				vectorTensor = Matrix1234Double::MatrixPadWithZeros(vectorTensor, UFunc.root->level);
			}
			vectorTensor = Matrix1234Double::PromoteInterleavedTo12(vectorTensor);
			UFunc = Matrix1234Double::PromoteInterleavedTo12(UFunc);
			CFLOBDD_DOUBLE a = Matrix1234Double::MatrixMultiplyV4(UFunc, vectorTensor);
			a = Matrix1234Double::Demote12ToInterleaved(a);
			//a = Matrix1234Double::NormalizeOutputTo1(a);
			/*a = 0.5 * a;
			a = (1.0 / pow(2, n/2)) * a;*/
			//std::cout << "a matrix created.." << std::endl;
			CFLOBDD_DOUBLE H = KroneckerPower(Matrix1234Double::MkWalshInterleaved, n, 1);
			CFLOBDD_DOUBLE I = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, n, 1);

			H = Matrix1234Double::PromoteInterleavedTo12(H);
			I = Matrix1234Double::PromoteInterleavedTo12(I);
			CFLOBDD_DOUBLE HI = Matrix1234Double::KroneckerProduct(H, I);
			
			if (HI.root->level < a.root->level)
				HI = Matrix1234Double::MatrixPadWithZeros(HI, a.root->level);
			//std::cout << "HI matrix created.." << std::endl;

			HI = Matrix1234Double::PromoteInterleavedTo12(HI);
			a = Matrix1234Double::PromoteInterleavedTo12(a);
			CFLOBDD_DOUBLE b = Matrix1234Double::MatrixMultiplyV4(HI, a);
			b = Matrix1234Double::Demote12ToInterleaved(b);
			//std::cout << "b matrix created.." << std::endl;
			b = VectorDouble::MatrixToVector(b);
			b = VectorDouble::VectorWithAmplitude(b);
			//std::cout << "b vector created.." << std::endl;
			CFLOBDD_DOUBLE ans = b;

			/*CFLOBDD_DOUBLE column1Matrix = VectorDouble::MkColumn1Matrix(ans.root->level);
			ans = Matrix1234Double::PromoteInterleavedTo12(ans);
			column1Matrix = Matrix1234Double::PromoteInterleavedTo12(column1Matrix);
			ans = Matrix1234Double::MatrixMultiplyV4(ans, column1Matrix);
			ans = Matrix1234Double::Demote12ToInterleaved(ans);
			ans = VectorDouble::MatrixToVector(ans);
			ans = Matrix1234Double::NormalizeOutputTo1(ans);
			ans = (2.0 / pow(2, n)) * ans;*/
			std::cout << ans.root->rootConnection.returnMapHandle.Size() << " ";
			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n+1,false);
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
#ifdef PATH_COUNTING_ENABLED
				s = VectorDouble::Sampling(ans);
#endif
				s = s.substr(0, n);
				addEquations(howellMatrix, s);
				iter++;
				//std::cout << iter << std::endl;
			}
			//std::cout << iter << " ";
			
			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));
			
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return ans;
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return answers;
		}

		CFLOBDD_FLOAT_BOOST KroneckerPower(CFLOBDD_FLOAT_BOOST(*f)(unsigned int), unsigned int n, unsigned int basicLevel)
		{
			unsigned int numtemp = n;
			unsigned int Walshlevel = floor(log2(numtemp)) + 1;
			CFLOBDD_FLOAT_BOOST H = f(Walshlevel + basicLevel - 1);
			numtemp = numtemp - (1 << (Walshlevel - 1));
			while (numtemp > 0)
			{
				Walshlevel = floor(log2(numtemp)) + 1;
				CFLOBDD_FLOAT_BOOST Htemp = f(Walshlevel + basicLevel - 1);
				Htemp = Matrix1234FloatBoost::MatrixPadWithZeros(Htemp, H.root->level);
				Htemp = Matrix1234FloatBoost::PromoteInterleavedTo12(Htemp);
				H = Matrix1234FloatBoost::PromoteInterleavedTo12(H);
				H = Matrix1234FloatBoost::KroneckerProduct(H, Htemp);
				numtemp = numtemp - (1 << (Walshlevel - 1));
			}

			return H;
		}

		CFLOBDD_COMPLEX_BIG KroneckerPower(CFLOBDD_COMPLEX_BIG(*f)(unsigned int), unsigned int n, unsigned int basicLevel)
		{
			unsigned int numtemp = n;
			unsigned int Walshlevel = floor(log2(numtemp)) + 1;
			CFLOBDD_COMPLEX_BIG H = f(Walshlevel + basicLevel - 1);
			numtemp = numtemp - (1 << (Walshlevel - 1));
			/*while (numtemp > 0)
			{
				Walshlevel = floor(log2(numtemp)) + 1;
				CFLOBDD_COMPLEX_BIG Htemp = f(Walshlevel + basicLevel - 1);
				Htemp = Matrix1234ComplexFloatBoost::MatrixPadWithZeros(Htemp, H.root->level);
				Htemp = Matrix1234ComplexFloatBoost::PromoteInterleavedTo12(Htemp);
				H = Matrix1234ComplexFloatBoost::PromoteInterleavedTo12(H);
				H = Matrix1234ComplexFloatBoost::KroneckerProduct(H, Htemp);
				numtemp = numtemp - (1 << (Walshlevel - 1));
			}*/

			return H;
		}

		CFLOBDD_FLOAT_BOOST CreateFuncMatrix(CFLOBDD_FLOAT_BOOST F, int n)
		{
			CFLOBDD_FLOAT_BOOST C = KroneckerPower(Matrix1234FloatBoost::MkCNOTInterleaved, n, 2);
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			std::cout << "I matrix created..." << std::endl;
			F = Matrix1234FloatBoost::PromoteInterleavedTo12(F);
			I = Matrix1234FloatBoost::PromoteInterleavedTo12(I);
			CFLOBDD_FLOAT_BOOST FI = Matrix1234FloatBoost::KroneckerProduct(F, I);
			FI = Matrix1234FloatBoost::NormalizeOutputTo1(FI);

			F = Matrix1234FloatBoost::Demote12ToInterleaved(F);
			CFLOBDD_FLOAT_BOOST Ft = Matrix1234FloatBoost::MatrixTranspose(F);
			Ft = Matrix1234FloatBoost::PromoteInterleavedTo12(Ft);
			CFLOBDD_FLOAT_BOOST FtI = Matrix1234FloatBoost::KroneckerProduct(Ft, I);
			FtI = Matrix1234FloatBoost::NormalizeOutputTo1(FtI);
			if (FtI.root->level < C.root->level)
				FtI = Matrix1234FloatBoost::MatrixPadWithZeros(FtI, C.root->level);

			FtI = Matrix1234FloatBoost::PromoteInterleavedTo12(FtI);
			C = Matrix1234FloatBoost::PromoteInterleavedTo12(C);
			CFLOBDD_FLOAT_BOOST temp = Matrix1234FloatBoost::MatrixMultiplyV4(FtI, C);
			temp = Matrix1234FloatBoost::NormalizeOutputTo1(temp);
			std::cout << "FtI * C done .." << std::endl;
			if (FI.root->level < C.root->level - 1)
				FI = Matrix1234FloatBoost::MatrixPadWithZeros(FI, C.root->level - 1);
			FI = Matrix1234FloatBoost::PromoteInterleavedTo12(FI);
			CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4(temp, FI);
			ans = Matrix1234FloatBoost::NormalizeOutputTo1(ans);
			std::cout << "temp * FI done.." << std::endl;
			ans = Matrix1234FloatBoost::Demote12ToInterleaved(ans);
			return ans;
		}

		/*
		std::vector<std::string> SimonsAlgo(int n, CFLOBDD_FLOAT_BOOST F)
		{
			auto m1 = high_resolution_clock::now();
			
			CFLOBDD_FLOAT_BOOST UFunc = CreateFuncMatrix(F, n);
			std::cout << "UFunc created ..." << std::endl;
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);

			SumEx = VectorFloatBoost::VectorToMatrixInterleaved(SumEx);
			e0 = VectorFloatBoost::VectorToMatrixInterleaved(e0);
			SumEx = Matrix1234FloatBoost::PromoteInterleavedTo12(SumEx);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo12(e0);
			CFLOBDD_FLOAT_BOOST vectorTensor = Matrix1234FloatBoost::KroneckerProduct(SumEx, e0);

			BIG_FLOAT coeff = 0.5 * mp::pow(BIG_FLOAT(2), -n / 2);
			vectorTensor = coeff * vectorTensor;

			if (UFunc.root->level > vectorTensor.root->level){
				vectorTensor = Matrix1234FloatBoost::MatrixPadWithZeros(vectorTensor, UFunc.root->level);
			}
			vectorTensor = Matrix1234FloatBoost::PromoteInterleavedTo12(vectorTensor);
			UFunc = Matrix1234FloatBoost::PromoteInterleavedTo12(UFunc);
			CFLOBDD_FLOAT_BOOST a = Matrix1234FloatBoost::MatrixMultiplyV4(UFunc, vectorTensor);
			a = Matrix1234FloatBoost::Demote12ToInterleaved(a);

			printMemory();
			std::cout << CFLOBDDNodeHandle::canonicalNodeTable->Size() << std::endl;
			std::cout << CFLOBDDReturnMapHandle::canonicalReturnMapBodySet->Size() << std::endl;
			UFunc = VectorFloatBoost::NoDistinctionNode(UFunc.root->level);
			SumEx = VectorFloatBoost::NoDistinctionNode(SumEx.root->level);
			e0 = VectorFloatBoost::NoDistinctionNode(e0.root->level);
			vectorTensor = VectorFloatBoost::NoDistinctionNode(vectorTensor.root->level);
			DisposeOfPairProductCache();
			InitPairProductCache();
			CFLOBDDNodeHandle::DisposeOfReduceCache();
			CFLOBDDNodeHandle::InitReduceCache();
			std::cout << CFLOBDDNodeHandle::canonicalNodeTable->Size() << std::endl;
			std::cout << CFLOBDDReturnMapHandle::canonicalReturnMapBodySet->Size() << std::endl;
			printMemory();
			
			std::cout << "a matrix created.." << std::endl;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);

			H = Matrix1234FloatBoost::PromoteInterleavedTo12(H);
			I = Matrix1234FloatBoost::PromoteInterleavedTo12(I);
			CFLOBDD_FLOAT_BOOST HI = Matrix1234FloatBoost::KroneckerProduct(H, I);

			if (HI.root->level < a.root->level)
				HI = Matrix1234FloatBoost::MatrixPadWithZeros(HI, a.root->level);
			std::cout << "HI matrix created.." << std::endl;
			
			HI = Matrix1234FloatBoost::PromoteInterleavedTo12(HI);
			a = Matrix1234FloatBoost::PromoteInterleavedTo12(a);
			
			CFLOBDD_FLOAT_BOOST b = Matrix1234FloatBoost::MatrixMultiplyV4(HI, a);
			b = Matrix1234FloatBoost::Demote12ToInterleaved(b);
			std::cout << "b matrix created.." << std::endl;
			std::cout << b.IsValid() << std::endl;
			//b = VectorFloatBoost::MatrixToVector(b);
			b = VectorFloatBoost::VectorWithAmplitude(b);
			//std::cout << "b vector created.." << std::endl;
			CFLOBDD_FLOAT_BOOST ans = b;
			unsigned int nodeCount;
			ans.CountNodes(nodeCount);
			std::cout << nodeCount << " ";
			std::cout << ans.root->rootConnection.returnMapHandle.Size() << " ";
			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			auto m2 = high_resolution_clock::now();
			auto duration = duration_cast<seconds>(m2 - m1);
			std::cout << duration.count() << " ";
			std::cout << std::endl;
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = VectorFloatBoost::Sampling(ans);
				s = s.substr(0, n);
				//std::cout << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
				//std::cout << iter << std::endl;
			}
			auto m3 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m3 - m2);
			std::cout << duration.count() << " ";

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m4 - m3);
			std::cout << duration.count() << " ";

			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return ans;
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return answers;
		}
		*/
		
		
		std::vector<std::string> SimonsAlgo(int n, CFLOBDD_FLOAT_BOOST F)
		{
			auto m1 = high_resolution_clock::now();
			CFLOBDD_FLOAT_BOOST UFunc = CreateFuncMatrix(F, n);
			DisposeOfPairProductCache();
			InitPairProductCache();
			CFLOBDDNodeHandle::DisposeOfReduceCache();
			CFLOBDDNodeHandle::InitReduceCache();
			std::cout << CFLOBDDNodeHandle::canonicalNodeTable->Size() << std::endl;
			std::cout << "UFunc returnMapSize: " << UFunc.root->rootConnection.returnMapHandle.Size() << std::endl;
			std::cout << "UFunc created ..." << std::endl;
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level, 1);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);

			SumEx = VectorFloatBoost::NoDistinctionNode(level+1, 1);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo13(e0);
			//SumEx = VectorFloatBoost::VectorToMatrixInterleaved(SumEx);
			//e0 = VectorFloatBoost::VectorToMatrixInterleaved(e0);
			SumEx = Matrix1234FloatBoost::PromoteInterleavedTo12(SumEx);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo12(e0);
			CFLOBDD_FLOAT_BOOST vectorTensor = Matrix1234FloatBoost::KroneckerProduct(SumEx, e0);
			
			std::cout << "Is vectorTensor valid: " << vectorTensor.IsValid() << std::endl;
			//BIG_FLOAT coeff = 0.5 * mp::pow(BIG_FLOAT(2), -n / 2);
			//vectorTensor = coeff * vectorTensor;

			if (UFunc.root->level > vectorTensor.root->level)
				vectorTensor = Matrix1234FloatBoost::MatrixPadWithZeros(vectorTensor, UFunc.root->level);

			vectorTensor = Matrix1234FloatBoost::PromoteInterleavedTo12(vectorTensor);
			UFunc = Matrix1234FloatBoost::PromoteInterleavedTo12(UFunc);
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4(UFunc, vectorTensor);
			vectorTensor = Matrix1234FloatBoost::Demote12ToInterleaved(vectorTensor);
			std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle.Size() << std::endl;
			UFunc = VectorFloatBoost::NoDistinctionNode(UFunc.root->level, 1);
			SumEx = VectorFloatBoost::NoDistinctionNode(SumEx.root->level, 1);
			e0 = VectorFloatBoost::NoDistinctionNode(e0.root->level, 1);
			vectorTensor = Matrix1234FloatBoost::NormalizeOutputTo1(vectorTensor);
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			H = Matrix1234FloatBoost::PromoteInterleavedTo12(H);
			I = Matrix1234FloatBoost::PromoteInterleavedTo12(I);
			CFLOBDD_FLOAT_BOOST HI = Matrix1234FloatBoost::KroneckerProduct(H, I);
			std::cout << "HI matrix created.." << std::endl;

			CFLOBDD_FLOAT_BOOST a = HI * vectorTensor;
			std::cout << "a returnMapSize: " << a.root->rootConnection.returnMapHandle.Size() << std::endl;
			CFLOBDD_FLOAT_BOOST allOnes = VectorFloatBoost::NoDistinctionNode(a.root->level + 1, 1);
			std::cout << "a matrix created.." << std::endl;
			a = Matrix1234FloatBoost::PromoteInterleavedTo12(a);
			a = Matrix1234FloatBoost::NormalizeOutputTo1(a);
			DisposeOfPairProductCache();
			InitPairProductCache();
			CFLOBDDNodeHandle::DisposeOfReduceCache();
			CFLOBDDNodeHandle::InitReduceCache();
			HI = VectorFloatBoost::NoDistinctionNode(HI.root->level, 1);
			vectorTensor = VectorFloatBoost::NoDistinctionNode(vectorTensor.root->level, 1);
			CFLOBDD_FLOAT_BOOST b = Matrix1234FloatBoost::MatrixMultiplyV4(allOnes, a);
			b = Matrix1234FloatBoost::Demote12ToInterleaved(b);
			//CFLOBDD_FLOAT_BOOST b = Matrix1234FloatBoost::AddMatrixRows(a);
			std::cout << "b matrix created.." << std::endl;
			std::cout << b.IsValid() << std::endl;
			//b = VectorFloatBoost::MatrixToVector(b);
			b = VectorFloatBoost::VectorWithAmplitude(b);
			b = Matrix1234FloatBoost::NormalizeOutputTo1(b);
			//std::cout << "b vector created.." << std::endl;
			CFLOBDD_FLOAT_BOOST ans = b;
			unsigned int numNodesOfAns = 0;
			ans.CountNodes(numNodesOfAns);
			std::cout << "numNodesOfAns: " << numNodesOfAns << std::endl;
			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			auto m2 = high_resolution_clock::now();
			auto duration = duration_cast<seconds>(m2 - m1);
			std::cout << duration.count() << " ";
			std::cout << std::endl;
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
#ifdef PATH_COUNTING_ENABLED
				s = VectorFloatBoost::SamplingV2(ans);
#endif
				//std::cout << s << std::endl;
				s = s.substr(0, n);
				std::cout << "iter: " << iter << " s: " << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
				//std::cout << iter << std::endl;
			}
			auto m3 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m3 - m2);
			std::cout << duration.count() << " ";

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m4 - m3);
			std::cout << duration.count() << " ";
			std::cout << matrix_s << std::endl;
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return ans;
			}
			
			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return answers;
		}


		std::vector<std::string> SimonsAlgoV2(int n, CFLOBDD_FLOAT_BOOST F)
		{
			CFLOBDD_FLOAT_BOOST C = KroneckerPower(Matrix1234FloatBoost::MkCNOTInterleaved, n, 2);
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			std::cout << "I matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			std::cout << "H matrix created..." << std::endl;
			
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST ans = VectorFloatBoost::NoDistinctionNode(level + 1, 1);
			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level + 1, 1);

			F = Matrix1234FloatBoost::PromoteInterleavedTo12(F);
			SumEx = Matrix1234FloatBoost::PromoteInterleavedTo12(SumEx);

			ans = Matrix1234FloatBoost::MatrixMultiplyV4(F, SumEx);

			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo12(e0);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo12(e0);
			ans = Matrix1234FloatBoost::KroneckerProduct(ans, e0);
			e0 = VectorFloatBoost::NoDistinctionNode(level, 1);
			SumEx = VectorFloatBoost::NoDistinctionNode(level, 1);

			unsigned int numNodes = 0;
			ans.CountNodes(numNodes);
			std::cout << "ans created with numNodes: " << numNodes << std::endl;

			F = Matrix1234FloatBoost::Demote12ToInterleaved(F);
			CFLOBDD_FLOAT_BOOST Ft = Matrix1234FloatBoost::MatrixTranspose(F);
			F = VectorFloatBoost::NoDistinctionNode(level, 1);
			Ft = Matrix1234FloatBoost::PromoteInterleavedTo12(Ft);
			H = Matrix1234FloatBoost::PromoteInterleavedTo12(H);
			CFLOBDD_FLOAT_BOOST tmp = VectorFloatBoost::NoDistinctionNode(level + 1, 1);
			tmp = Matrix1234FloatBoost::MatrixMultiplyV4(H, Ft);
			I = Matrix1234FloatBoost::PromoteInterleavedTo12(I);
			tmp = Matrix1234FloatBoost::KroneckerProduct(tmp, I);
			H = VectorFloatBoost::NoDistinctionNode(level, 1);
			I = VectorFloatBoost::NoDistinctionNode(level, 1);
			Ft = VectorFloatBoost::NoDistinctionNode(level, 1);

			tmp.CountNodes(numNodes);
			std::cout << "tmp created with numNodes: " << numNodes << std::endl;

			C = Matrix1234FloatBoost::PromoteInterleavedTo12(C);
			ans = Matrix1234FloatBoost::PromoteInterleavedTo12(ans);
			ans = Matrix1234FloatBoost::MatrixMultiplyV4(C, ans);

			ans.CountNodes(numNodes);
			std::cout << "C * ans created with numNodes: " << numNodes << std::endl;

			//ans = Matrix1234FloatBoost::Demote12ToInterleaved(ans);
			tmp = Matrix1234FloatBoost::PromoteInterleavedTo12(tmp);
			ans = Matrix1234FloatBoost::KroneckerProduct(tmp, ans);
			CFLOBDDInternalNode* x = (CFLOBDDInternalNode *)ans.root->rootConnection.entryPointHandle->handleContents;
			std::cout << "numBConnections: " << x->numBConnections << std::endl;

			tmp = Matrix1234FloatBoost::PromoteInterleavedTo12(tmp);
			ans = Matrix1234FloatBoost::MatrixMultiplyV4(tmp, ans);
			tmp = VectorFloatBoost::NoDistinctionNode(level, 1);
			ans.CountNodes(numNodes);
			std::cout << "tmp * C * ans created with numNodes: " << numNodes << std::endl;

			ans = Matrix1234FloatBoost::Demote12ToInterleaved(ans);
			ans = VectorFloatBoost::VectorWithAmplitude(ans);
			ans.CountNodes(numNodes);
			std::cout << "ans vector amp with numNodes: " << numNodes << std::endl;

			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
#ifdef PATH_COUNTING_ENABLED
				s = VectorFloatBoost::Sampling(ans, false);
#endif
				//std::cout << s << std::endl;
				s = s.substr(0, n);
				std::cout << "iter: " << iter << " s: " << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
				//std::cout << iter << std::endl;
			}

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			std::cout << matrix_s << std::endl;
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ansV;
				ansV.push_back(s);
				return ansV;
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return answers;
		}

		std::vector<std::string> SimonsAlgoV3(int n, CFLOBDD_FLOAT_BOOST F)
		{
			auto m1 = high_resolution_clock::now();
			CFLOBDD_FLOAT_BOOST C = KroneckerPower(Matrix1234FloatBoost::MkCNOTInterleaved, n, 2);
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			std::cout << "H matrix created..." << std::endl;
			F = Matrix1234FloatBoost::PromoteInterleavedTo12(F);
			H = Matrix1234FloatBoost::PromoteInterleavedTo12(H);
			CFLOBDD_FLOAT_BOOST HF = Matrix1234FloatBoost::KroneckerProduct(H, F);
			std::cout << "HF matrix created..." << std::endl;
			CFLOBDDNodeHandle::DisposeOfReduceCache();
			CFLOBDDNodeHandle::InitReduceCache();
			std::cout << CFLOBDDNodeHandle::canonicalNodeTable->Size() << std::endl;
			HF = Matrix1234FloatBoost::PromoteInterleavedTo12(HF);
			C = Matrix1234FloatBoost::PromoteInterleavedTo12(C);
			//CFLOBDD_FLOAT_BOOST temp = Matrix1234FloatBoost::MatrixMultiplyV4(HF, C);
			CFLOBDDNodeHandle::DisposeOfReduceCache();
			CFLOBDDNodeHandle::InitReduceCache();
			std::cout << CFLOBDDNodeHandle::canonicalNodeTable->Size() << std::endl;

			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level, 1);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);
			SumEx = VectorFloatBoost::VectorToMatrixInterleaved(SumEx);
			SumEx = Matrix1234FloatBoost::PromoteInterleavedTo12(SumEx);
			e0 = VectorFloatBoost::VectorToMatrixInterleaved(e0);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo12(e0);
			CFLOBDD_FLOAT_BOOST vectorTensor = Matrix1234FloatBoost::KroneckerProduct(SumEx, e0);
			/*Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(vectorTensor, std::cout);
			CFLOBDD_FLOAT_BOOST x = Matrix1234FloatBoost::Demote12ToInterleaved(temp);
			Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(x, std::cout);*/

			std::cout << "Is vectorTensor valid: " << vectorTensor.IsValid() << std::endl;
			
			if (C.root->level-1 > vectorTensor.root->level)
				vectorTensor = Matrix1234FloatBoost::MatrixPadWithZeros(vectorTensor, C.root->level);
			
			vectorTensor = Matrix1234FloatBoost::PromoteInterleavedTo12(vectorTensor);
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4(C, vectorTensor);
			std::cout << "C * vectorTensor created..." << std::endl;
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4(HF, vectorTensor);
			std::cout << "HF * vectorTensor created..." << std::endl;
			vectorTensor = Matrix1234FloatBoost::Demote12ToInterleaved(vectorTensor);
			
			std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle.Size() << std::endl;
			F = VectorFloatBoost::NoDistinctionNode(F.root->level, 1);
			//temp = VectorFloatBoost::NoDistinctionNode(temp.root->level);
			SumEx = VectorFloatBoost::NoDistinctionNode(SumEx.root->level, 1);
			e0 = VectorFloatBoost::NoDistinctionNode(e0.root->level, 1);
			//Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(vectorTensor, std::cout);
			BIG_FLOAT coeff = mp::pow(BIG_FLOAT(2), -n);
			vectorTensor = coeff * vectorTensor;
			vectorTensor = VectorFloatBoost::VectorWithAmplitude(vectorTensor);

			//std::cout << "b vector created.." << std::endl;
			CFLOBDD_FLOAT_BOOST ans = vectorTensor;
			//Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			unsigned int numNodesOfAns = 0;
			ans.CountNodes(numNodesOfAns);
			std::cout << "numNodesOfAns: " << numNodesOfAns << std::endl;
			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			auto m2 = high_resolution_clock::now();
			auto duration = duration_cast<seconds>(m2 - m1);
			std::cout << duration.count() << " ";
			std::cout << std::endl;
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
#ifdef PATH_COUNTING_ENABLED
				s = VectorFloatBoost::Sampling(ans, false);
#endif
				//std::cout << s << std::endl;
				s = s.substr(0, n);
				std::cout << "iter: " << iter << " s: " << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
				//std::cout << iter << std::endl;
			}
			auto m3 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m3 - m2);
			std::cout << duration.count() << " ";

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m4 - m3);
			std::cout << duration.count() << " ";
			std::cout << matrix_s << std::endl;
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return ans;
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return answers;
		}

		std::pair<CFLOBDD_FLOAT_BOOST, std::vector<std::string>> SimonsAlgoV4(int n, CFLOBDD_FLOAT_BOOST F)
		{
			
			auto m1 = high_resolution_clock::now();
			CFLOBDD_FLOAT_BOOST C = KroneckerPower(Matrix1234FloatBoost::MkCNOTInterleaved, n, 2);
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			std::cout << "H matrix created..." << std::endl;
			F = Matrix1234FloatBoost::PromoteInterleavedTo12(F);
			H = Matrix1234FloatBoost::PromoteInterleavedTo12(H);
			CFLOBDD_FLOAT_BOOST HF = Matrix1234FloatBoost::KroneckerProduct(H, F);
			std::cout << "HF matrix created..." << std::endl;
			std::cout << HF.root->rootConnection.returnMapHandle << std::endl;

			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level, 1);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);
			SumEx = VectorFloatBoost::VectorToMatrixInterleaved(SumEx);
			SumEx = Matrix1234FloatBoost::PromoteInterleavedTo12(SumEx);
			e0 = VectorFloatBoost::VectorToMatrixInterleaved(e0);
			e0 = Matrix1234FloatBoost::PromoteInterleavedTo12(e0);
			CFLOBDD_FLOAT_BOOST vectorTensor = Matrix1234FloatBoost::KroneckerProduct(SumEx, e0);

			std::cout << "Is vectorTensor valid: " << vectorTensor.IsValid() << std::endl;

			if (C.root->level - 1 > vectorTensor.root->level)
				vectorTensor = Matrix1234FloatBoost::MatrixPadWithZeros(vectorTensor, C.root->level);
			// unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			// vectorTensor.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			// std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(C, vectorTensor);
			std::cout << "C * vectorTensor created... " << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			// vectorTensor.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			// std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(HF, vectorTensor);
			std::cout << "HF * vectorTensor created..." << std::endl;
			// vectorTensor.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			// std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			
			std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle.Size() << std::endl;
			BIG_FLOAT coeff = mp::pow(BIG_FLOAT(2), -n);
			vectorTensor = coeff * vectorTensor;
			//std::cout << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			vectorTensor = VectorFloatBoost::VectorWithAmplitude(vectorTensor);
			// vectorTensor.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			// std::cout << "Step 4 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			
			CFLOBDD_FLOAT_BOOST ans = vectorTensor;
			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			auto m2 = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(m2 - m1);
			std::cout << std::endl << "quantum_time: " << duration.count() << " ";
			std::cout << std::endl;
			//return std::make_pair(ans, std::vector<std::string>());
			vectorTensor.CountPaths();
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
//#ifdef PATH_COUNTING_ENABLED
				s = VectorFloatBoost::Sampling(ans, false);
//#endif
				//std::cout << s << std::endl;
				s = s.substr(0, n);
				//std::cout << "iter: " << iter << " s: " << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
			}
			auto m3 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m3 - m2);
			std::cout << "classical time: " << duration.count() << std::endl;

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m4 - m3);
			std::cout << " duration: " << duration.count() << " ";
			std::cout << matrix_s << std::endl;
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return std::make_pair(vectorTensor, ans);
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return std::make_pair(vectorTensor, answers);
			
		}

		std::pair<CFLOBDD_FLOAT_BOOST, std::vector<std::string>> SimonsAlgoV4New(int n, std::string s)
		{
			auto m1 = high_resolution_clock::now();
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			std::cout << "H matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST HI = Matrix1234FloatBoost::KroneckerProduct2Vocs(H, I);
			std::cout << "HF matrix created..." << std::endl;

			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level, 1);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);
			SumEx = VectorFloatBoost::VectorToMatrixInterleaved(SumEx);
			e0 = VectorFloatBoost::VectorToMatrixInterleaved(e0);
			CFLOBDD_FLOAT_BOOST vectorTensor = Matrix1234FloatBoost::KroneckerProduct2Vocs(SumEx, e0);


			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(Matrix1234FloatBoost::MkCNOT(level + 2, 2 * n, 0, n), vectorTensor);
			for (unsigned int i = 1; i < n; i++) {
				std::cout << i << std::endl;
				CFLOBDD_FLOAT_BOOST tmp = Matrix1234FloatBoost::MkCNOT(level + 2, 2 * n, i, n + i);
				vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(tmp, vectorTensor);
			}
			std::cout << "C * vectorTensor created... " << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			CFLOBDD_FLOAT_BOOST t = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, 2 * n, 1);
			int k = 0, m = 0;
			for (int i = n - 1; i >= 0; i--) {
				if (s[i] == '1') {
					m = n;
					for (int j = n - 1; j >= 0; j--) {
						if (s[j] == '1') {
							std::cout << i << " " << j << " " << k << " " << m << std::endl;
							CFLOBDD_FLOAT_BOOST tmp = Matrix1234FloatBoost::MkCNOT(level + 2, 2 * n, k, m);
							t = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(tmp, t);
						}
						m += 1;
					}
					break;
				}
				k += 1;
			}
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(t, vectorTensor);
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(HI, vectorTensor);
			std::cout << "HF * vectorTensor created..." << std::endl;

			//std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			BIG_FLOAT coeff = mp::pow(BIG_FLOAT(2), -n);
			vectorTensor = coeff * vectorTensor;
			//std::cout << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			vectorTensor = VectorFloatBoost::VectorWithAmplitude(vectorTensor);
			//std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			CFLOBDD_FLOAT_BOOST ans = vectorTensor;
			HowellMatrix::HowellMatrix<BitVector::BV1>* howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			auto m2 = high_resolution_clock::now();
			auto duration = duration_cast<seconds>(m2 - m1);
			std::cout << std::endl << " quantum_time: " << duration.count() << " ";
			std::cout << std::endl;
			//return std::make_pair(ans, std::vector<std::string>());
			ans.CountPaths();
			//std::cout << "path counts: " << ans.root->rootConnection.entryPointHandle.handleContents->numPathsToExit[0] << std::endl;
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
				s = VectorFloatBoost::Sampling(ans, true);
				s = s.substr(0, n);
				//std::cout << "iter: " << iter << " s: " << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
			}
			auto m3 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m3 - m2);
			std::cout << "classical time: " << duration.count() << std::endl;

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m4 - m3);
			std::cout << " duration: " << duration.count() << " ";
			std::cout << matrix_s << std::endl;
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return std::make_pair(vectorTensor, ans);
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return std::make_pair(vectorTensor, answers);

		}

		std::vector<std::string> SimonsAlgoV4_Voc2(int n, CFLOBDD_FLOAT_BOOST F)
		{

			auto m1 = high_resolution_clock::now();
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST C = Matrix1234FloatBoost::MkCNOT(level + 2, 2 * n, 0, n);
			for (unsigned int i = 1; i < n; i++){
				std::cout << i << std::endl;
				C = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(C, Matrix1234FloatBoost::MkCNOT(level + 2, 2 * n, i, n + i));
				/*unsigned int c_nodes;
				C.CountNodes(c_nodes);
				std::cout << "C nodeCount: " << c_nodes << " leaves: " << C.root->rootConnection.returnMapHandle.Size() << std::endl;*/
			}
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			std::cout << "H matrix created..." << std::endl;
			CFLOBDD_FLOAT_BOOST HF = Matrix1234FloatBoost::KroneckerProduct2Vocs(H, F);
			std::cout << "HF matrix created..." << std::endl;
			std::cout << HF.root->rootConnection.returnMapHandle << std::endl;

			CFLOBDD_FLOAT_BOOST SumEx = VectorFloatBoost::NoDistinctionNode(level, 1);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level, 0);
			SumEx = VectorFloatBoost::VectorToMatrixInterleaved(SumEx);
			e0 = VectorFloatBoost::VectorToMatrixInterleaved(e0);
			CFLOBDD_FLOAT_BOOST vectorTensor = Matrix1234FloatBoost::KroneckerProduct2Vocs(SumEx, e0);

			std::cout << "Is vectorTensor valid: " << vectorTensor.IsValid() << std::endl;

			if (C.root->level - 1 > vectorTensor.root->level)
				vectorTensor = Matrix1234FloatBoost::MatrixPadWithZeros(vectorTensor, C.root->level);

			//vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4(Matrix1234FloatBoost::PromoteInterleavedTo12(C), Matrix1234FloatBoost::PromoteInterleavedTo12(vectorTensor));
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(C, vectorTensor);
			std::cout << "C * vectorTensor created..." << std::endl;
			vectorTensor = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(HF, vectorTensor);
			std::cout << "HF * vectorTensor created..." << std::endl;

			std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle.Size() << std::endl;
			F = VectorFloatBoost::NoDistinctionNode(F.root->level, 1);
			//temp = VectorFloatBoost::NoDistinctionNode(temp.root->level);
			SumEx = VectorFloatBoost::NoDistinctionNode(SumEx.root->level, 1);
			e0 = VectorFloatBoost::NoDistinctionNode(e0.root->level, 1);
			//Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(vectorTensor, std::cout);
			BIG_FLOAT coeff = mp::pow(BIG_FLOAT(2), -n);
			vectorTensor = coeff * vectorTensor;
			vectorTensor = VectorFloatBoost::VectorWithAmplitude(vectorTensor);

			//std::cout << "b vector created.." << std::endl;
			CFLOBDD_FLOAT_BOOST ans = vectorTensor;
			//Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			unsigned int numNodesOfAns = 0;
			ans.CountNodes(numNodesOfAns);
			std::cout << "numNodesOfAns: " << numNodesOfAns << std::endl;
			HowellMatrix::HowellMatrix<BitVector::BV1> *howellMatrix = new HowellMatrix::HowellMatrix<BitVector::BV1>(n + 1, false);
			auto m2 = high_resolution_clock::now();
			auto duration = duration_cast<seconds>(m2 - m1);
			std::cout << duration.count() << " ";
			std::cout << std::endl;
			unsigned int iter = 1;
			//std::cout << "Sampling start.." << std::endl;
			while (iter <= 2 * n)
			{
				std::string s = "";
#ifdef PATH_COUNTING_ENABLED
				s = VectorFloatBoost::Sampling(ans, true);
#endif
				//std::cout << s << std::endl;
				s = s.substr(0, n);
				std::cout << "iter: " << iter << " s: " << s << std::endl;
				addEquations(howellMatrix, s);
				iter++;
				//std::cout << iter << std::endl;
			}
			auto m3 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m3 - m2);
			std::cout << duration.count() << " ";

			HowellMatrix::ModularSquareMatrix<BitVector::BV1> modMatrix = howellMatrix->getSquareMatrix();
			auto sqmatrix_s = HowellMatrix::ModularSquareMatrix<BitVector::BV1>::dualize(modMatrix);
			HowellMatrix::HowellMatrix<BitVector::BV1> matrix_s = HowellMatrix::HowellMatrix<BitVector::BV1>(*(sqmatrix_s.get_ptr()));

			auto m4 = high_resolution_clock::now();
			duration = duration_cast<seconds>(m4 - m3);
			std::cout << duration.count() << " ";
			std::cout << matrix_s << std::endl;
			if (matrix_s.size() == 1)
			{
				auto vector_s = matrix_s.getSparseRow(0).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}
				std::vector<std::string> ans;
				ans.push_back(s);
				return ans;
			}

			std::vector<std::string> answers;
			for (unsigned int j = 1; j < matrix_s.num_rows(); j++)
			{
				auto vector_s = matrix_s.getSparseRow(j).get_ptr();
				std::string s;

				for (unsigned int i = 1; i < vector_s->GetLength(); i++)
				{
					s += (((*vector_s)[i] == BitVector::BV1::zero()) ? "0" : "1");
				}

				answers.push_back(s);
			}

			return answers;

		}

		unsigned int getIndexFromBitString(std::string s){
			unsigned int index = 0;
			for (int i = 0; i < s.length(); i++){
				index = index * 2 + ((s[i] == '0') ? 0 : 1);
				index = index * 2 + ((s[i] == '0') ? 0 : 1);
			}
			return index;
		}

		std::string getIndexStringFromBitString(std::string s){
			std::string ans(s.length() * 2, '0');
			unsigned int j = 0;
			for (int i = 0; i < s.length(); i++){
				ans[j] = s[i];
				ans[j + 1] = s[i];
				j+=2;
			}
			return ans;
		}
		
		CFLOBDD_FLOAT_BOOST MkU_w(std::string s, int n){
			std::string index = getIndexStringFromBitString(s);
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST F = VectorFloatBoost::MkBasisVector(level + 1, index);
			// std::cout << index << std::endl;
			// Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(F, std::cout);
			F = -2 * F;
			return I + F;
		}

		CFLOBDD_FLOAT_BOOST MkU_s(int n){
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			BIG_FLOAT val = BIG_FLOAT(1.0) / mp::pow(BIG_FLOAT(2), n - 1);
			CFLOBDD_FLOAT_BOOST ans = val * VectorFloatBoost::NoDistinctionNode(level+1, 1) + (-1 *I);
			return ans;
		}
	
		std::string GroversAlgo(int n, std::string s){
			CFLOBDD_FLOAT_BOOST U_w = MkU_w(s, n);
			std::cout << "U_w created" << std::endl;
			CFLOBDD_FLOAT_BOOST U_s = MkU_s(n);
			std::cout << "U_s created" << std::endl;
			CFLOBDD_FLOAT_BOOST M = Matrix1234FloatBoost::MatrixMultiplyV4(Matrix1234FloatBoost::PromoteInterleavedTo12(U_s), Matrix1234FloatBoost::PromoteInterleavedTo12(U_w));
			std::cout << "M created" << std::endl;
			double const pi = 4 * std::atan(1);
			unsigned long long int iters = (unsigned long long int)floor(pi * 0.25 * pow(2, n / 2));
			unsigned int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST ans = (1.0 / pow(2, n / 2)) * VectorFloatBoost::NoDistinctionNode(level+1, 1);
			ans = Matrix1234FloatBoost::PromoteInterleavedTo12(ans);
			unsigned int count;
			M.CountNodes(count);
			std::cout << "Iter start, num iters: " << iters << " nodeCount: " << count << std::endl;
			for (unsigned long long int i = 0; i < iters; i++){
				ans = Matrix1234FloatBoost::MatrixMultiplyV4(M, ans);
				ans.CountNodes(count);
				std::cout << "iter : " << i << " nodeCount: " << count << " return size: " << ans.root->rootConnection.returnMapHandle.Size() << std::endl;
			}
			std::cout << "Iter end" << std::endl;
			ans = Matrix1234FloatBoost::Demote12ToInterleaved(ans);
			unsigned int nodeCount;
			ans.CountNodes(nodeCount);
			std::cout << "NodeCount: " << nodeCount << std::endl;
			std::cout << ans.root->rootConnection.returnMapHandle.Size() << std::endl;
			//Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			return s;
		}

		std::pair<CFLOBDD_FLOAT_BOOST, CFLOBDD_FLOAT_BOOST>
			MultiplyRec(CFLOBDD_FLOAT_BOOST M, mp::cpp_int iters,
			boost::unordered_map<mp::cpp_int, CFLOBDD_FLOAT_BOOST>& memo, unsigned int n, unsigned int* level){
			//std::cout << iters << std::endl;
			auto it = memo.find(iters);
			if (it != memo.end())
				return std::make_pair(it->second, it->second);
			if (iters == 1){
				BIG_FLOAT val = BIG_FLOAT(1.0);// / mp::pow(BIG_FLOAT(2), n / 4);
				CFLOBDD_FLOAT_BOOST ans = M;
				if ((*level) >= n) {
					memo[iters] = ans;
					//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
					//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
					//std::cout << "Step " << iters << " : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				}
				if ((*level) < n){
					*level = *level + 1;
					return std::make_pair((BIG_FLOAT(1.0) / sqrt(2)) * ans, ans);
				}
				return std::make_pair(ans, ans);
			}
			mp::cpp_int half_iters = iters / 2;
			auto F1 = MultiplyRec(M, half_iters, memo, n, level);
			auto F2 = MultiplyRec(M, iters - half_iters, memo, n, level);
			CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F1.first, F2.first);
			CFLOBDD_FLOAT_BOOST true_ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F1.second, F2.second);
			if ((*level) >= n) {
				memo[iters] = true_ans;
				//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
				//true_ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
				//std::cout << "Step " << iters << " : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
				//std::cout << "iter: " << iters << " return: " << true_ans.root->rootConnection.returnMapHandle << std::endl;
			}
			/*if ((*level) < n) {
				*level = *level + 1;
				return std::make_pair((BIG_FLOAT(1.0) / sqrt(2)) * ans, true_ans);
			}*/
			return std::make_pair(ans, true_ans);
		}

		std::pair<std::string, CFLOBDD_FLOAT_BOOST> GroversAlgoWithV4(int n, std::string s){
			CFLOBDD_FLOAT_BOOST U_w = MkU_w(s, n);
			std::cout << "U_w created" << std::endl;
			CFLOBDD_FLOAT_BOOST U_s = MkU_s(n);
			std::cout << "U_s created" << std::endl;
			// Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(U_w, std::cout);
			// Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(U_s, std::cout);
			CFLOBDD_FLOAT_BOOST M = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(U_s, U_w);
			// Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(M, std::cout);
			std::cout << "M created" << std::endl;
			std::cout << "M leaves count: " << M.root->rootConnection.returnMapHandle.Size() << std::endl;
			std::cout << "M leaves: " << M.root->rootConnection.returnMapHandle << std::endl;
			
			double const pi = 4 * std::atan(1);
			BIG_FLOAT tmp_val(mp::pow(BIG_FLOAT(2.0), n / 2) * boost::math::constants::pi<BIG_FLOAT>() * 0.25);
			mp::cpp_int iters = mp::floor(tmp_val).convert_to<mp::cpp_int>();
			unsigned int level = ceil(log2(n));
			BIG_FLOAT val = BIG_FLOAT(1.0);
			if (n - iters > 0)
				val = val / mp::pow(BIG_FLOAT(2), (BIG_FLOAT(n) - BIG_FLOAT(iters)) / 2);
			CFLOBDD_FLOAT_BOOST ans = val * VectorFloatBoost::NoDistinctionNode(level + 1, 1);
			std::cout << "Iter start, num iters: " << iters << std::endl;
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			boost::unordered_map<mp::cpp_int, CFLOBDD_FLOAT_BOOST> memo;
			unsigned int ulevel = 0;
			auto M_rec = MultiplyRec(M, iters, memo, n, &ulevel);
			CFLOBDD_FLOAT_BOOST M_actual = M_rec.first;
			// Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(M_actual, std::cout);
			std::cout << ulevel << std::endl;
			std::cout << M_actual.root->rootConnection.returnMapHandle << std::endl;
			ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(M_actual, ans);
			// Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step i : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;

			std::cout << "Iter end" << std::endl;
			ans = VectorFloatBoost::VectorWithAmplitude(ans);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step j : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans.CountPaths();
			std::cout << ans.root->rootConnection.returnMapHandle << std::endl;
			std::cout << ans.root->rootConnection.entryPointHandle->handleContents->numPathsToExit[1] << std::endl;
			// std::cout << ans << std::endl;
			std::string ans_s = "";
			int count = 0;
			// while (s != ans_s){
				ans_s = VectorFloatBoost::Sampling(ans, true, "Grovers").substr(0, n);
			// 	count++;
			// }
			// std::cout << "count: " << count << std::endl;
			return std::make_pair(ans_s, ans);
		}


		CFLOBDD_DOUBLE MkU_w_double(std::string s, int n){
			unsigned int index = getIndexFromBitString(s);
			unsigned int level = ceil(log2(n));
			CFLOBDD_DOUBLE I = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, n, 1);
			CFLOBDD_DOUBLE F = VectorDouble::MkBasisVector(level + 1, index);
			F = -2 * F;
			return I + F;
		}

		CFLOBDD_DOUBLE MkU_s_double(int n){
			unsigned int level = ceil(log2(n));
			CFLOBDD_DOUBLE I = KroneckerPower(Matrix1234Double::MkIdRelationInterleaved, n, 1);
			CFLOBDD_DOUBLE ans = (2.0 / pow(2, n)) * VectorDouble::NoDistinctionNode(level + 1) + (-1 * I);
			return ans;
		}

		std::string GroversAlgoWithV4_double(int n, std::string s){
			CFLOBDD_DOUBLE U_w = MkU_w_double(s, n);
			std::cout << "U_w created" << std::endl;
			CFLOBDD_DOUBLE U_s = MkU_s_double(n);
			std::cout << "U_s created" << std::endl;
			CFLOBDD_DOUBLE M = Matrix1234Double::MatrixMultiplyV4(U_s, U_w);
			std::cout << "M created" << std::endl;
			double const pi = 4 * std::atan(1);
			unsigned long long int iters = (unsigned long long int)floor(pi * 0.25 * pow(2, n / 2));
			unsigned int level = ceil(log2(n));
			CFLOBDD_DOUBLE ans = (1.0 / pow(2, n / 2)) * VectorDouble::NoDistinctionNode(level + 1);
			unsigned int count;
			M.CountNodes(count);
			std::cout << " M node count: " << count << " leaf count: " << M.root->rootConnection.returnMapHandle.Size() << std::endl;
			std::cout << M.root->rootConnection.returnMapHandle << std::endl;
			//std::cout << "Iter start, num iters: " << iters << " nodeCount: " << count << std::endl;
			/*for (unsigned long long int i = 0; i < iters; i++){
				ans = Matrix1234Double::MatrixMultiplyV4(M, ans);
				//ans.CountNodes(count);
				//std::cout << "iter : " << i << " nodeCount: " << count << " return size: " << ans.root->rootConnection.returnMapHandle.Size() << std::endl;
			}
			std::cout << "Iter end" << std::endl;
			//unsigned int nodeCount;
			//ans.CountNodes(nodeCount);
			//std::cout << "NodeCount: " << nodeCount << std::endl;
			std::cout << ans.root->rootConnection.returnMapHandle.Size() << std::endl;*/
			//Matrix1234FloatBoost::MatrixPrintRowMajorInterleaved(ans, std::cout);
			return s;
		}

		std::pair<std::string, CFLOBDD_FLOAT_BOOST> DeutschJozsaAlgo(unsigned int n, CFLOBDD_FLOAT_BOOST F)
		{
			int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST X = VectorFloatBoost::NoDistinctionNode(level, 1);
			std::string last_one(n, '0');
			last_one[0] = '1';
			CFLOBDD_FLOAT_BOOST Y = VectorFloatBoost::MkBasisVector(level, last_one);
			Y = VectorFloatBoost::VectorToMatrixInterleaved(Y);
			X = VectorFloatBoost::VectorToMatrixInterleaved(X);
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST HI = Matrix1234FloatBoost::KroneckerProduct2Vocs(H, I);
			Y = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(H, Y);
			CFLOBDD_FLOAT_BOOST ans = Matrix1234FloatBoost::KroneckerProduct2Vocs(X, Y);
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F, ans);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(HI, ans);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans = VectorFloatBoost::VectorWithAmplitude(ans);
			//ans.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 4 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			ans.CountPaths();
			std::string ans_s = VectorFloatBoost::Sampling(ans, true).substr(0, n);
			return std::make_pair(ans_s, ans);
		}

		std::pair<std::string, CFLOBDD_FLOAT_BOOST> BV(long long int n, CFLOBDD_FLOAT_BOOST F){
			int level = ceil(log2(n));
			CFLOBDD_FLOAT_BOOST e = VectorFloatBoost::NoDistinctionNode(level, 1);
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::MkBasisVector(level + 1, 0);
			e = VectorFloatBoost::VectorToMatrixInterleaved(e);
			CFLOBDD_FLOAT_BOOST stateV = Matrix1234FloatBoost::KroneckerProduct2Vocs(e, e0);
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 0 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F, stateV);
			BIG_FLOAT coeff = mp::pow(BIG_FLOAT(2), -n / 2);
			stateV = coeff * stateV;
			CFLOBDD_FLOAT_BOOST H = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST I = KroneckerPower(Matrix1234FloatBoost::MkIdRelationInterleaved, n, 1);
			CFLOBDD_FLOAT_BOOST HI = Matrix1234FloatBoost::KroneckerProduct2Vocs(H, I);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(HI, stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
		    //std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = VectorFloatBoost::VectorWithAmplitude(stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV.CountPaths();
			std::string ans_s = "";
			while (ans_s.find('1') == std::string::npos){
				ans_s = VectorFloatBoost::Sampling(stateV, true).substr(0, n);
			}
			return std::make_pair(ans_s, stateV);
		}

		std::pair<std::string, CFLOBDD_FLOAT_BOOST> GHZ(unsigned long long int n){
			int level = ceil(log2(n)) + 2;
			CFLOBDD_FLOAT_BOOST F = Matrix1234FloatBoost::MkCNOT(level, n, 0, n);
			std::cout << "Starting loop" << std::endl;
			for (unsigned int i = 1; i < n; i++){
				CFLOBDD_FLOAT_BOOST tmp = Matrix1234FloatBoost::MkCNOT(level, n, i, n);
				F = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F, tmp);
			}
			CFLOBDD_FLOAT_BOOST e0 = VectorFloatBoost::NoDistinctionNode(level - 1, 1);
			//CFLOBDD_FLOAT_BOOST e1 = VectorFloatBoost::MkBasisVector(level - 1, pow(2, 4 * n) - 1);
			std::string last_one(2 * n, '0');
			last_one[2 * n - 1] = '1';
			CFLOBDD_FLOAT_BOOST e1 = VectorFloatBoost::MkBasisVector(level - 1, last_one);
			//std::cout << "e1, etmp: " << (e1 == etmp) << std::endl;
			CFLOBDD_FLOAT_BOOST stateV = Matrix1234FloatBoost::KroneckerProduct2Vocs(e0, e1);
			//unsigned int nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount;
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 1 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			CFLOBDD_FLOAT_BOOST H2 = KroneckerPower(Matrix1234FloatBoost::MkWalshInterleaved, 2 * n, 1);
			stateV = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(F, stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 2 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = Matrix1234FloatBoost::MatrixMultiplyV4WithInfo(H2, stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 3 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV = VectorFloatBoost::VectorWithAmplitude(stateV);
			//stateV.CountNodesAndEdges(nodeCount, edgeCount, returnEdgeCount, returnEdgeObjCount);
			//std::cout << "Step 4 : " << nodeCount << " " << edgeCount << " " << (nodeCount + edgeCount) << std::endl;
			stateV.CountPaths();
			std::string ans_s = VectorFloatBoost::Sampling(stateV, true).substr(0, n + 1);
			return std::make_pair(ans_s, stateV);
		}

		std::pair<CFLOBDD_COMPLEX_BIG, std::vector<std::string>> ShorsAlgo(int n, CFLOBDD_COMPLEX_BIG F)
		{
			int level = ceil(log2(n));
			CFLOBDD_COMPLEX_BIG C = KroneckerPower(Matrix1234ComplexFloatBoost::MkCNOTInterleaved, n, 2);
			std::cout << "C matrix created..." << std::endl;
			CFLOBDD_FOURIER Q = Matrix1234Fourier::MkFourierMatrixInterleavedV4WithInfo(level + 1);
			CFLOBDD_COMPLEX_BIG QP = Matrix1234ComplexFloatBoost::ConvertToComplex(Q);
			std::cout << "QP matrix created..." << std::endl;
			F = Matrix1234ComplexFloatBoost::PromoteInterleavedTo12(F);
			QP = Matrix1234ComplexFloatBoost::PromoteInterleavedTo12(QP);
			CFLOBDD_COMPLEX_BIG QF = Matrix1234ComplexFloatBoost::KroneckerProduct(QP, F);
			std::cout << "HF matrix created..." << std::endl;

			CFLOBDD_COMPLEX_BIG SumEx = VectorComplexFloatBoost::NoDistinctionNode(level, BIG_COMPLEX_FLOAT(1));
			CFLOBDD_COMPLEX_BIG e0 = VectorComplexFloatBoost::MkBasisVector(level, 0);
			SumEx = VectorComplexFloatBoost::VectorToMatrixInterleaved(SumEx);
			SumEx = Matrix1234ComplexFloatBoost::PromoteInterleavedTo12(SumEx);
			e0 = VectorComplexFloatBoost::VectorToMatrixInterleaved(e0);
			e0 = Matrix1234ComplexFloatBoost::PromoteInterleavedTo12(e0);
			CFLOBDD_COMPLEX_BIG vectorTensor = Matrix1234ComplexFloatBoost::KroneckerProduct(SumEx, e0);

			std::cout << "Is vectorTensor valid: " << vectorTensor.IsValid() << std::endl;

			vectorTensor = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, vectorTensor);
			std::cout << "C * vectorTensor created... " << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			vectorTensor = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(QF, vectorTensor);
			std::cout << "HF * vectorTensor created..." << std::endl;

			std::cout << "vectorTensor returnMapSize: " << vectorTensor.root->rootConnection.returnMapHandle.Size() << std::endl;
			BIG_COMPLEX_FLOAT coeff = mp::pow(BIG_COMPLEX_FLOAT(2), -n);
			vectorTensor = coeff * vectorTensor;
			//std::cout << vectorTensor.root->rootConnection.returnMapHandle << std::endl;
			vectorTensor = VectorComplexFloatBoost::VectorWithAmplitude(vectorTensor);

			CFLOBDD_COMPLEX_BIG ans = vectorTensor;
			vectorTensor.CountPaths();
			for (int k = 0; k < 10; k++){
				std::string s = VectorComplexFloatBoost::Sampling(ans, false).substr(0, n);
				std::cout << s << std::endl;
			}
			std::vector<std::string> answers{ "" };
			return std::make_pair(vectorTensor, answers);

		}

		CFLOBDD_COMPLEX_BIG Hadamard(unsigned int n, unsigned int i)
		{
			if (n == 1)
			{
				return Matrix1234ComplexFloatBoost::MkWalshInterleaved(1);
			}
			else {
				if (i < n/2)
				{
					CFLOBDD_COMPLEX_BIG T = KroneckerPower(Matrix1234ComplexFloatBoost::MkIdRelationInterleaved, n/2, 1);
					CFLOBDD_COMPLEX_BIG H = Hadamard(n/2, i);
					return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(H, T);
				}
				else
				{
					CFLOBDD_COMPLEX_BIG T = KroneckerPower(Matrix1234ComplexFloatBoost::MkIdRelationInterleaved, n/2, 1);
					return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(T, Hadamard(n/2, i - n/2)); 
				}
			}
		}

		CFLOBDD_COMPLEX_BIG QFT(long long int n, std::string s)
		{
			unsigned int level = ceil(log2(n));
			
			// std::cout << "s: " << s << std::endl;
			// std::reverse(s.begin(), s.end());
			CFLOBDD_COMPLEX_BIG stateV = VectorComplexFloatBoost::MkBasisVector(level, s);
			stateV = VectorComplexFloatBoost::VectorToMatrixInterleaved(stateV);
			std::cout << "start" << std::endl;

			for (long long int i = 0; i < n/2; i++)
			{
				CFLOBDD_COMPLEX_BIG SwapM = Matrix1234ComplexFloatBoost::MkSwapGate(level+1, i, n-i-1);
				// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(SwapM, std::cout);
				stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(SwapM, stateV);
				// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(stateV, std::cout);
			}

			std::cout << "loop start" << std::endl;

			// CFLOBDD_COMPLEX_BIG CP = Matrix1234ComplexFloatBoost::MkCPGate(level+1, 2, 3, 0.5);
			// // std::cout << *(CP.root->rootConnection.entryPointHandle) << std::endl;
			// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(CP, std::cout);
			// CP.CountPaths();
			// std::cout << CP.root->rootConnection.entryPointHandle->handleContents->numPathsToExit[0] << std::endl;
			// std::cout << CP.root->rootConnection.entryPointHandle->handleContents->numPathsToExit[1] << std::endl;

			for (long long int i = n-1; i >= 0; i--)
			{
				CFLOBDD_COMPLEX_BIG H = Hadamard(n, i);
				// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(H, std::cout);	
				stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, stateV);
				// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(stateV, std::cout);
				// std::cout << "hey" << std::endl;
				// CFLOBDD_COMPLEX_BIG tmp = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level+1);
				for (long int j = 0; j < i; j++)
				{
					double theta = std::pow(2, j - i);
					//std::cout << j << " " << i << std::endl;
					CFLOBDD_COMPLEX_BIG CP = Matrix1234ComplexFloatBoost::MkCPGate(level+1, j, i, theta);
					// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(CP, std::cout);
					// tmp = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CP, tmp);
					stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CP, stateV);
					// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(stateV, std::cout);
					// std::cout << "over " << theta << std::endl;
				}
				// stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(tmp, stateV);
			}
			std::cout << "done" << std::endl;

			// Matrix1234ComplexFloatBoost::MatrixPrintRowMajorInterleaved(stateV, std::cout);

			return stateV;
		}

		CFLOBDD_COMPLEX_BIG ShorsAlgoNew(int a, int N)
		{
			unsigned int level = ceil(log2(N));
			CFLOBDD_COMPLEX_BIG allOnes = VectorComplexFloatBoost::NoDistinctionNode(level+1,1);
			std::string s(N, '0');
			s[3] = '1';
			CFLOBDD_COMPLEX_BIG e = VectorComplexFloatBoost::MkBasisVector(level, s);
			e = VectorComplexFloatBoost::VectorToMatrixInterleaved(e);
			CFLOBDD_COMPLEX_BIG stateV = Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(allOnes, e);

			for (int q = N-1; q >= 0; q--)
			{
				// std::cout << "q: " << q << std::endl;
				unsigned int power = std::pow(2, N-1-q);
				for (unsigned int i = 0; i < power; i++)
				{
					if (a == 2 || a == 13)
					{
						CFLOBDD_COMPLEX_BIG CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N, N+1);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
						CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N+1, N+2);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
						CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N+2, N+3);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
					}
					if (a == 7 || a == 8)
					{
						CFLOBDD_COMPLEX_BIG CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N+2, N+3);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
						CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N+1, N+2);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
						CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N, N+1);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
					}
					if (a == 4 || a == 11)
					{
						CFLOBDD_COMPLEX_BIG CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N+1, N+3);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
						CSWAP = Matrix1234ComplexFloatBoost::MkCSwapGate(level+2, q, N, N+2);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CSWAP, stateV);
					}
					if (a == 7 || a == 11 || a == 13)
					{
						CFLOBDD_COMPLEX_BIG X = Matrix1234ComplexFloatBoost::MkExchangeInterleaved(level);
						CFLOBDD_COMPLEX_BIG Id = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level);
						CFLOBDD_COMPLEX_BIG XI = Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(X, Id);
						CFLOBDD_COMPLEX_BIG I = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level+1);
						CFLOBDD_COMPLEX_BIG IXI = Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(I, XI);
						stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(IXI, stateV);
					}
				}
			}
			CFLOBDD_COMPLEX_BIG I = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level+1);
			for (long long int i = 0; i < N; i++)
			{
				for (long int j = 0; j < i; j++)
				{
					double theta = -1*std::pow(2, j - i);
					// std::cout << j << " " << i << std::endl;
					CFLOBDD_COMPLEX_BIG CP = Matrix1234ComplexFloatBoost::MkCPGate(level+1, j, i, theta);
					CFLOBDD_COMPLEX_BIG CPI = Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(CP, I);
					stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(CPI, stateV);
				}
				CFLOBDD_COMPLEX_BIG H = Hadamard(2*N, i);
				stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, stateV);
			}
			for (long long int i = 0; i < N/2; i++)
			{
				CFLOBDD_COMPLEX_BIG SwapM = Matrix1234ComplexFloatBoost::MkSwapGate(level+1, i, N-i-1);
				CFLOBDD_COMPLEX_BIG SwapMI = Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(SwapM, I);
				stateV = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(SwapMI, stateV);
			}

			return stateV;
		}
	}
}

