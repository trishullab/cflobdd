#ifndef INTPAIR_GUARD
#define INTPAIR_GUARD

//
//    Copyright (c) 1999 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <iostream>
#include <fstream>
#include <cstdint>

// Packed representation: first in upper 32 bits, second in lower 32 bits.
// Single uint64_t enables one-instruction comparison and fast hashing.
class intpair {
 public:
  intpair() : packed(0) {}
  intpair(const int i1, const int i2)
    : packed(((uint64_t)(unsigned int)i1 << 32) | (unsigned int)i2) {}
  intpair operator! () {
    return intpair(!First(), !Second());
  }
  bool operator!= (const intpair& p) const { return packed != p.packed; }
  friend bool operator==(const intpair& lhs, const intpair& rhs) { return lhs.packed == rhs.packed; }
  int First() const { return (int)(packed >> 32); }
  int Second() const { return (int)(packed & 0xFFFFFFFF); }
  uint64_t getPacked() const { return packed; }
  struct intpair_hash {
    size_t operator()(const intpair& p) const {
      // fmix64 finalizer from MurmurHash3
      uint64_t h = p.packed;
      h ^= h >> 33;
      h *= 0xff51afd7ed558ccdULL;
      h ^= h >> 33;
      h *= 0xc4ceb9fe1a85ec53ULL;
      h ^= h >> 33;
      return (size_t)h;
    }
  };
 private:
  uint64_t packed;
};

typedef uint64_t packed_intpair;

struct packed_intpair_hash {
  size_t operator()(packed_intpair p) const {
    // fmix64 finalizer from MurmurHash3
    p ^= p >> 33;
    p *= 0xff51afd7ed558ccdULL;
    p ^= p >> 33;
    p *= 0xc4ceb9fe1a85ec53ULL;
    p ^= p >> 33;
    return (size_t)p;
  }
};

std::ostream& operator<< (std::ostream & out, const intpair &p);

#endif
