#ifndef _QUANTUM_VERIFIER
#define _QUANTUM_VERIFIER

#include "matrix1234_complex_float_boost.h"

using namespace CFL_OBDD;

class QuantumVerifier {
    public:
        // Constructor
        QuantumVerifier(unsigned int numQubits);
        // Constructor
        QuantumVerifier();
        // Destructor
        ~QuantumVerifier();
        // set qubit count;
        void setNumQubits(unsigned int numQubits);
        // I = [[1 0] [0 1]]
        // For no-op
        void ApplyIdentityGate(unsigned int index);
        // H = [[1 1] [1 -1]]
        // Also called Walsh Gate
        void ApplyHadamardGate(unsigned int index);
        // X = [[0 1] [1 0]]
        // Also called Pauli-X or bit-flip
        void ApplyNOTGate(unsigned int index);
        // Y = [[0 -i] [i 0]]
        void ApplyPauliYGate(unsigned int index);
        // Z = [[1 0] [0 -1]]
        // Also called phase-flip
        void ApplyPauliZGate(unsigned int index);
        // S = [[1 0] [0 i]]
        // Also called √(Z) Gate
        void ApplySGate(unsigned int index);
        // CNOT = [[1 0 0 0] [0 1 0 0] [0 0 0 1] [0 0 1 0]]
        void ApplyCNOTGate(long int controller, long int controlled);
        // Ph = e^{i phase} I
        void ApplyGlobalPhase(double phase);
        // SWAP = [[1 0 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1]]
        void ApplySwapGate(long int index1, long int index2);
        // iSWAP = [[1 0 0 0] [0 0 i 0] [0 i 0 0] [0 0 0 1]]
        void ApplyiSwapGate(long int index1, long int index2);
        // CZ = [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 -1]]
        // Also called Controlled Phase Flip Gate
        void ApplyCZGate(long int controller, long int controlled);
        // CPhase = [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 e^(i * π * θ)]]
        void ApplyCPGate(long int controller, long int controlled, double theta);
        // P = [[1 0] [0 e^{i * π * θ}]]
        void ApplyPhaseShiftGate(unsigned int index, double theta);
        // T = P(1/4) = [[1 0] [0 e^{i * π * 1/4}]]
        void ApplyTGate(unsigned int index);
        // CS = [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 e^(i * π * 1/2)]]
        void ApplyCSGate(long int controller, long int controlled);
        // CCNOT or Toffoli gate
        void ApplyCCNOTGate(long int controller1, long int controller2, long int controlled);
        // CSWAP or Fredkin gate
        void ApplyCSwapGate(long int controller, long int index1, long int index2); 
    private:
        CFLOBDD_COMPLEX_BIG stateVector;
        unsigned int numQubits;
};

#endif