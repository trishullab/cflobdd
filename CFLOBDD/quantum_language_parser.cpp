#include "quantum_language_parser.h"
#include <cstring>
#include <iostream>
#include <fstream>

/* Tokenizing a string */
std::vector<std::string> tokenizer(const std::string& line)  {
    std::vector<std::string> tokens;
    std::stringstream   mySstream( line );
    std::string         temp;

    while( std::getline( mySstream, temp, ' ' ) ) {
        tokens.push_back( temp );
    }

    return tokens;
}

QuantumLanguageParser::QuantumLanguageParser() {}

QuantumLanguageParser::QuantumLanguageParser(std::string filename)
{
    ReadFromFile(filename);
}

QuantumLanguageParser::~QuantumLanguageParser()
{
    inp_file.close();
}

void QuantumLanguageParser::ReadFromFile(std::string filename)
{
    inp_file.open(filename);
    std::string inst;
    bool first_line = true;
    long int numQubits;
    while (std::getline (inp_file, inst)) {
        // Output the text from the file

        std::vector<std::string> tokens = tokenizer(inst);

        if (tokens[0] == "n")
        {
            if (first_line != true)
            {
                std::cout << "Error: number of qubits should be the first line" << std::endl;
                abort();
            }
            first_line = false;
            numQubits = std::stol(tokens[1]);
        }
        else
        {
            Instruction instruct;
            instruct.inst_name = inst_name_to_enum[tokens[0]];
            for (unsigned int i = 1; i < tokens.size(); i++)
            {
                instruct.indices.push_back(tokens[i]);
            }
            insts.push_back(instruct);
        }
    }
    qv = QuantumVerifier(numQubits);
}

void QuantumLanguageParser::Execute()
{
    for (unsigned int i = 0; i < insts.size(); i++)
    {
        switch (insts[i].inst_name)
        {
            case Identity:
                qv.ApplyIdentityGate(std::stoi(insts[i].indices[0]));
                break;
            case Hadamard:
                qv.ApplyHadamardGate(std::stoi(insts[i].indices[0]));
                break;
            case NOT:
                qv.ApplyNOTGate(std::stoi(insts[i].indices[0]));
                break;
            case PauliY:
                qv.ApplyPauliYGate(std::stoi(insts[i].indices[0]));
                break;
            case PauliZ:
                qv.ApplyPauliZGate(std::stoi(insts[i].indices[0]));
                break;
            case S:
                qv.ApplySGate(std::stoi(insts[i].indices[0]));
                break;
            case CNOT:
                qv.ApplyCNOTGate(std::stol(insts[i].indices[0]), stol(insts[i].indices[1]));
                break;
            case GlobalPhase:
                qv.ApplyGlobalPhase(std::stod(insts[i].indices[0]));
                break;
            case Swap:
                qv.ApplySwapGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]));
                break;
            case iSwap:
                qv.ApplyiSwapGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]));
                break;
            case CZ:
                qv.ApplyCZGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]));
                break;
            case CP:
                qv.ApplyCPGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]), std::stod(insts[i].indices[2]));
                break;
            case T:
                qv.ApplyTGate(std::stoi(insts[i].indices[0]));
                break;
            case CS:
                qv.ApplyCSGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]));
                break;
            case CCNOT:
                qv.ApplyCCNOTGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]), std::stol(insts[i].indices[2]));
                break;
            case CSwap:
                qv.ApplyCSwapGate(std::stol(insts[i].indices[0]), std::stol(insts[i].indices[1]), std::stol(insts[i].indices[2]));
                break; 
            default:
                break;
        }
    }
}