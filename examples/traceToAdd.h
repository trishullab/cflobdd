/**
  @file traceToAdd.h

  @brief Build an ADD directly from a token stream (trace).

  The input of length N is the left-justified yield of a decision tree
  of height h = ceil(log_2(N)), padded on the right with -1.

  The token stream is provided by a Tokenizer object.  The default
  tokenizer (CharTokenizer) treats each byte as a separate token with
  value 0..255.  A flex-based tokenizer can be substituted by providing
  a different Tokenizer implementation.
*/

#ifndef TRACE_TO_ADD_H_
#define TRACE_TO_ADD_H_

#include "../cudd-3.0.0/cplusplus/cuddObj.hh"
#include <string>

// =========================================================================
// Tokenizer interface
// =========================================================================

/// Abstract tokenizer: produces a stream of unsigned int token values.
class Tokenizer {
public:
    virtual ~Tokenizer() {}

    /// Open the input source.  Returns true on success.
    virtual bool open(const std::string &filename) = 0;

    /// Read the next token.  Returns true if a token was read,
    /// false at end of input.  The token value is stored in *value.
    virtual bool next(unsigned int *value) = 0;

    /// Return the total number of tokens in the input,
    /// or 0 if the count is not known in advance.
    /// If known, the tokenizer should be reset to the beginning
    /// after counting (ready for next() calls).
    virtual unsigned long count() = 0;

    /// Close the input source.
    virtual void close() = 0;
};

// =========================================================================
// Result type
// =========================================================================

struct TraceADDResult {
    ADD add;
    unsigned int numVars;        // number of Boolean variables used
    unsigned long traceLength;   // number of tokens read
};

// =========================================================================
// Public interface
// =========================================================================

/// Build an ADD from a token stream using the binary-counter merge algorithm.
/// The tokenizer must be opened before calling this function.
/// The ADD uses ceil(log2(N)) variables, with variable 0 as MSB.
TraceADDResult traceToADD(Cudd &mgr, Tokenizer &tok);

/// Build an ADD from a file using the default character tokenizer.
/// Convenience wrapper around traceToADD.
TraceADDResult traceFileToADD(Cudd &mgr, const std::string &filename);

#endif // TRACE_TO_ADD_H_
