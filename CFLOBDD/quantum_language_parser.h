#ifndef _QUANTUM_LANGUAGE_PARSER
#define _QUANTUM_LANGUAGE_PARSER

#include "quantum_verifier.h"
#include <fstream>
#include <unordered_map>

enum INST_NAME
{
    Identity,
    Hadamard,
    NOT,
    PauliY,
    PauliZ,
    S,
    CNOT,
    GlobalPhase,
    Swap,
    iSwap,
    CZ,
    CP,
    T,
    CS,
    CCNOT,
    CSwap
};

std::unordered_map<std::string, INST_NAME> inst_name_to_enum = {
    {"ID" , INST_NAME::Identity},
    {"H" , INST_NAME::Hadamard},
    {"NOT" , INST_NAME::NOT},
    {"PauliY" , INST_NAME::PauliY},
    {"PauliZ" , INST_NAME::PauliZ},
    {"S" , INST_NAME::S},
    {"CNOT" , INST_NAME::CNOT},
    {"GP" , INST_NAME::GlobalPhase},
    {"SWAP" , INST_NAME::Swap},
    {"iSWAP" , INST_NAME::iSwap},
    {"CZ" , INST_NAME::CZ},
    {"CP" , INST_NAME::CP},
    {"T" , INST_NAME::T},
    {"CS" , INST_NAME::CS},
    {"CCNOT" , INST_NAME::CCNOT},
    {"CSWAP" , INST_NAME::CSwap}
};

struct Instruction
{
    INST_NAME inst_name;
    std::vector<std::string> indices;
};

class QuantumLanguageParser
{
    public:
        QuantumLanguageParser(std::string filename);
        QuantumLanguageParser();
        ~QuantumLanguageParser();
        void ReadFromFile(std::string filename);
        void Execute();
    private:
        QuantumVerifier qv;
        std::vector<Instruction> insts;
        std::ifstream inp_file;
};

#endif