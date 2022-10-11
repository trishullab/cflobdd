#include <stdio.h>
#include "../cudd-complex/cudd/cudd.h"
#include "../cudd-complex/util/util.h"
#include <complex.h>
#include <math.h>
#include <time.h>

DdNode* identity_matrix(DdManager* manager, DdNode* x, DdNode* y){
    DdNode* x_neg = Cudd_addCmpl(manager, x);
    DdNode* y_neg = Cudd_addCmpl(manager, y);

    Cudd_Ref(x_neg);
    Cudd_Ref(y_neg);

    DdNode* f1 = Cudd_addApply(manager, Cudd_addTimes, x, y);
    Cudd_Ref(f1);

    DdNode* f2 = Cudd_addApply(manager, Cudd_addTimes, x_neg, y_neg);
    Cudd_Ref(f2);

    DdNode* F = Cudd_addApply(manager, Cudd_addPlus, f1, f2);
    Cudd_Ref(F);

    Cudd_RecursiveDeref(manager, f1);
    Cudd_RecursiveDeref(manager, f2);
    Cudd_RecursiveDeref(manager, x_neg);
    Cudd_RecursiveDeref(manager, y_neg);

    return F;
}

DdNode* hadamard_matrix(DdManager* manager, DdNode* x, DdNode* y){
    DdNode* x_neg = Cudd_addCmpl(manager, x);
    DdNode* y_neg = Cudd_addCmpl(manager, y);

    Cudd_Ref(x_neg);
    Cudd_Ref(y_neg);

    DdNode* F = Cudd_addApply(manager, Cudd_addMinus, y_neg, y);
    Cudd_Ref(F);

    DdNode* tmp = Cudd_addApply(manager, Cudd_addTimes, F, x);
    Cudd_RecursiveDeref(manager, F);
    Cudd_Ref(tmp);
    F = tmp;

    tmp = Cudd_addApply(manager, Cudd_addPlus, x_neg, F);
    Cudd_RecursiveDeref(manager, F);
    Cudd_Ref(tmp);
    F = tmp;

    Cudd_RecursiveDeref(manager, x_neg);
    Cudd_RecursiveDeref(manager, y_neg);

    return F;
}


DdNode* identity_n(DdManager* manager, unsigned int start, unsigned int end, DdNode** x, DdNode** y){
    if (start == end - 1)
        return identity_matrix(manager, x[start], y[start]);
    unsigned int mid = (end - start)/2 + start;
    DdNode* F1 = identity_n(manager, start, mid, x, y);
    DdNode* F2 = identity_n(manager, mid, end, x, y);
    Cudd_Ref(F1);
    Cudd_Ref(F2);
    DdNode* F = Cudd_addApply(manager, Cudd_addTimes, F1, F2);
    Cudd_Ref(F);
    Cudd_RecursiveDeref(manager, F1);
    Cudd_RecursiveDeref(manager, F2);
    return F;
}

DdNode* D_n(DdManager* manager, unsigned int start, unsigned int end, unsigned int size, DdNode** x, DdNode** y){
    int N = end - start;
    long int nodeCount = 0;
    double const pi = 4 * atan(1);
    double temp = pow(2, size); 
    double temp2 = 1.0/temp;
    double w = 2.0 * pi * temp2;
    DdNode* ans = Cudd_addConst(manager, 0);
    Cudd_Ref(ans);
    unsigned int counter = pow(2, N);
    for (unsigned int i = 0; i < counter; ++i)
    {
        DdNode* tmp_ans = Cudd_addConst(manager, 1);
        Cudd_Ref(tmp_ans);
        unsigned int tmp_i = i;
        unsigned int count = 0;
        while (count < N){
            int j = tmp_i % 2;
            tmp_i = tmp_i/2;

            if (j == 1){
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, x[start + N-1-count]);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
            }
            else{
                DdNode* x_neg = Cudd_addCmpl(manager, x[start + N-1-count]);
                Cudd_Ref(x_neg);
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, x_neg);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
                Cudd_RecursiveDeref(manager, x_neg);
            }

            if (j == 1){
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, y[start + N-1-count]);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
            }
            else{
                DdNode* y_neg = Cudd_addCmpl(manager, y[start + N-1-count]);
                Cudd_Ref(y_neg);
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, y_neg);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
                Cudd_RecursiveDeref(manager, y_neg);
            }

            count++;
        }
        double val = w * ((i) % ((unsigned int)pow(2, size)));
        double complex w_pow = cos(val) + I * sin(val);
        DdNode* w_pow_node = Cudd_addConst(manager, w_pow);
        Cudd_Ref(w_pow_node);
        DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, w_pow_node);
        Cudd_Ref(node);
        Cudd_RecursiveDeref(manager, tmp_ans);
        tmp_ans = node;

        node = Cudd_addApply(manager, Cudd_addPlus, ans, tmp_ans);
        Cudd_Ref(node);
        Cudd_RecursiveDeref(manager, ans);
        ans = node;
        Cudd_RecursiveDeref(manager, tmp_ans);
        // Cudd_RecursiveDeref(manager, node);
        Cudd_RecursiveDeref(manager, w_pow_node);
        nodeCount = Cudd_ReadNodeCount(manager);
        // printf("i: %d NodeCount: %ld\n", i, nodeCount);
    }
    Cudd_Deref(ans);
    return ans;
}


DdNode* MkCyclicMatrix(DdManager* manager, DdNode** x, DdNode** y, unsigned int N){
    DdNode* ans = Cudd_addConst(manager, 0);
    Cudd_Ref(ans);
    for (unsigned int i = 0; i < pow(2, N)/2; ++i)
    {
        DdNode* tmp_ans = Cudd_addConst(manager, 1);
        Cudd_Ref(tmp_ans);

        unsigned int tmp_i = i;
        unsigned int val_i = 2*i;
        unsigned int count = 0;
        while (count < N){
            int tmp_j = tmp_i % 2;
            tmp_i = tmp_i/2;

            int val_j = val_i % 2;
            val_i = val_i/2;

            if (tmp_j == 1){
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, x[N-1-count]);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
            }else{
                DdNode* x_neg = Cudd_addCmpl(manager, x[N-1-count]);
                Cudd_Ref(x_neg);
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, x_neg);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
                Cudd_RecursiveDeref(manager, x_neg);
            }

            if (val_j == 1){
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, y[N-1-count]);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
            }else{
                DdNode* y_neg = Cudd_addCmpl(manager, y[N-1-count]);
                Cudd_Ref(y_neg);
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, y_neg);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
                Cudd_RecursiveDeref(manager, y_neg);
            }
            count++;
        }

        DdNode* node = Cudd_addApply(manager, Cudd_addPlus, ans, tmp_ans);
        Cudd_Ref(node);
        Cudd_RecursiveDeref(manager, ans);
        ans = node;
        Cudd_RecursiveDeref(manager, tmp_ans);
    }

    for (unsigned int i = 0; i < pow(2, N)/2; ++i)
    {
        DdNode* tmp_ans = Cudd_addConst(manager, 1);
        Cudd_Ref(tmp_ans);

        unsigned int tmp_i = i + pow(2, N)/2;
        unsigned int val_i = 2*i+1;
        unsigned int count = 0;
        while (count < N){
            int tmp_j = tmp_i % 2;
            tmp_i = tmp_i/2;

            int val_j = val_i % 2;
            val_i = val_i/2;

            if (tmp_j == 1){
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, x[N-1-count]);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
            }else{
                DdNode* x_neg = Cudd_addCmpl(manager, x[N-1-count]);
                Cudd_Ref(x_neg);
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, x_neg);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
                Cudd_RecursiveDeref(manager, x_neg);
            }

            if (val_j == 1){
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, y[N-1-count]);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
            }else{
                DdNode* y_neg = Cudd_addCmpl(manager, y[N-1-count]);
                Cudd_Ref(y_neg);
                DdNode* node = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, y_neg);
                Cudd_Ref(node);
                Cudd_RecursiveDeref(manager, tmp_ans);
                tmp_ans = node;
                Cudd_RecursiveDeref(manager, y_neg);
            }
            count++;
        }

        DdNode* node = Cudd_addApply(manager, Cudd_addPlus, ans, tmp_ans);
        Cudd_Ref(node);
        Cudd_RecursiveDeref(manager, ans);
        ans = node;
        Cudd_RecursiveDeref(manager, tmp_ans);
    }

    return ans;
}

DdNode* Fourier(DdManager* manager, DdNode** x, DdNode** y, DdNode** w, DdNode** z, unsigned int N){
    if (N == 1){
        DdNode* ans = hadamard_matrix(manager, x[0], y[0]);
        Cudd_Ref(ans);
        DdNode* constant = Cudd_addConst(manager, 1.0/sqrt(2));
        Cudd_Ref(constant);
        DdNode* tmp_ans = Cudd_addApply(manager, Cudd_addTimes, ans, constant);
        Cudd_Ref(tmp_ans);
        Cudd_RecursiveDeref(manager, ans);
        Cudd_RecursiveDeref(manager, constant);
        return tmp_ans;

    }

    printf("start %d\n", N);

    DdNode* K = MkCyclicMatrix(manager, z, y, N);
    Cudd_Ref(K);
    
    DdNode* F_recurse = Fourier(manager, x+1, y+1, w+1, z+1, N-1);
    Cudd_Ref(F_recurse);
    // printf("F_recurse\n");
    // Cudd_PrintDebug(manager, F_recurse, 2*N, 2);
    

    DdNode* x_neg = Cudd_addCmpl(manager, x[0]);
    DdNode* y_neg = Cudd_addCmpl(manager, y[0]);
    Cudd_Ref(x_neg);
    Cudd_Ref(y_neg);
    DdNode* F1 = Cudd_addApply(manager, Cudd_addTimes, x_neg, y_neg);
    DdNode* F2 = Cudd_addApply(manager, Cudd_addTimes, x[0], y[0]);
    Cudd_Ref(F1);
    Cudd_Ref(F2);
    DdNode* F = Cudd_addApply(manager, Cudd_addPlus, F1, F2);
    Cudd_RecursiveDeref(manager, F2);
    Cudd_RecursiveDeref(manager, F1);
    Cudd_Ref(F);
    
    DdNode* tmp = Cudd_addApply(manager, Cudd_addTimes, F, F_recurse);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, F_recurse);
    Cudd_RecursiveDeref(manager, F);
    F = tmp;


    tmp = Cudd_addSwapVariables(manager, F, x, w, N);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, F);
    F = tmp;


    tmp = Cudd_addSwapVariables(manager, F, y, z, N);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, F);
    F = tmp;
    // Cudd_Ref(F);
    

    DdNode* Id = identity_n(manager, 1, N, x, y);
    Cudd_Ref(Id);
    

    DdNode* D = D_n(manager, 1, N, N, x, y);
    Cudd_Ref(D);


    DdNode* tmp1 = Cudd_addApply(manager, Cudd_addTimes, y_neg, Id);
    Cudd_Ref(tmp1);


    DdNode* tmp2 = Cudd_addApply(manager, Cudd_addMinus, x_neg, x[0]);
    Cudd_Ref(tmp2);

    tmp = Cudd_addApply(manager, Cudd_addTimes, y[0], tmp2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, tmp2);
    tmp2 = tmp;
    
    

    tmp = Cudd_addApply(manager, Cudd_addTimes, tmp2, D);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, D);
    Cudd_RecursiveDeref(manager, tmp2);
    tmp2 = tmp;
    


    DdNode* ID = Cudd_addApply(manager, Cudd_addPlus, tmp1, tmp2);
    Cudd_Ref(ID);
    Cudd_RecursiveDeref(manager, tmp1);
    Cudd_RecursiveDeref(manager, tmp2);
    

    tmp = Cudd_addSwapVariables(manager, ID, y, w, N);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, ID);
    ID = tmp;
    
    // printf("here\n");
    DdNode* IDF = Cudd_addMatrixMultiply(manager, ID, F, w, N);
    Cudd_Ref(IDF);
    // printf("here\n");
    // printf("IDF\n");
    // Cudd_PrintDebug(manager, IDF, 2*N, 2);
    // printf("K\n");
    // Cudd_PrintDebug(manager, K, 2*N, 2);

    Cudd_RecursiveDeref(manager, ID);
    Cudd_RecursiveDeref(manager, F);
    DdNode* tmp_ans = Cudd_addMatrixMultiply(manager, IDF, K, z, N);
    Cudd_Ref(tmp_ans);
    Cudd_RecursiveDeref(manager, IDF);
    Cudd_RecursiveDeref(manager, K);
    // printf("tmp_ans\n");
    // Cudd_PrintDebug(manager, tmp_ans, 2*N, 2);
    
    DdNode* constant = Cudd_addConst(manager, 1.0/sqrt(2));
    Cudd_Ref(constant);
    DdNode* ans = Cudd_addApply(manager, Cudd_addTimes, tmp_ans, constant);
    Cudd_Ref(ans);
    // printf("ans\n");
    // Cudd_PrintDebug(manager, ans, 2*N, 2);

    Cudd_RecursiveDeref(manager, constant);
    Cudd_RecursiveDeref(manager, tmp_ans);
    Cudd_RecursiveDeref(manager, y_neg);
    Cudd_RecursiveDeref(manager, x_neg);
    printf("end %d \n", N);

    return ans;
}

long int MkFourierTransform(DdManager* manager, int n){
    DdNode **x, **y, **w, **z;
    unsigned int N = pow(2, n);
    x = (DdNode **)malloc(sizeof(DdNode*) * N);
    y = (DdNode **)malloc(sizeof(DdNode*) * N);
    w = (DdNode **)malloc(sizeof(DdNode*) * N);
    z = (DdNode **)malloc(sizeof(DdNode*) * N);

    for (unsigned int i = 0; i < N; i++){
        x[i] = Cudd_addIthVar(manager, 4*i);
        y[i] = Cudd_addIthVar(manager, 4*i+1);
        w[i] = Cudd_addIthVar(manager, 4*i+2);
        z[i] = Cudd_addIthVar(manager, 4*i+3);

        Cudd_Ref(x[i]);
        Cudd_Ref(y[i]);
        Cudd_Ref(z[i]);
        Cudd_Ref(w[i]);
    }
    clock_t t;
    t = clock();
    DdNode* ans = Fourier(manager, x, y, w, z, N);
    // DdNode* ans = D_n(manager, 0, N, N, x, y);
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

    Cudd_Ref(ans);
    for (unsigned int i = 0; i < N; i++){
        Cudd_RecursiveDeref(manager, x[i]);
        Cudd_RecursiveDeref(manager, y[i]);
        Cudd_RecursiveDeref(manager, z[i]);
        Cudd_RecursiveDeref(manager, w[i]);
    }
    // Cudd_PrintDebug(manager, ans, 2*N, 2);
    long int nodeCount = Cudd_ReadNodeCount(manager);
    long int leafCount = Cudd_CountLeaves(ans);
    printf("NodeCount: %ld LeafCount: %ld time_taken: %f sec\n", nodeCount, leafCount, time_taken);
    return nodeCount;

}


int main (int argc, char** argv)
{
    if (argc < 2)
        return 1;
    unsigned int nvars = atoi(argv[2]);
    DdManager *manager;  /* pointer to DD manager */
    manager = Cudd_Init(nvars,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0);
    long int nodeCount = 0;

    if (strcmp(argv[1], "fourier") == 0)
      nodeCount = MkFourierTransform(manager, atoi(argv[2]));

    Cudd_Quit(manager);
    return 0;
}
