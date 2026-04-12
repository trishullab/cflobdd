// gen_luby: write a Luby sequence out to 2^n to a file named lubyM.txt
// where M = 2^n.
// Uses letters a..z for 2^0 through 2^25.
// Usage: gen_luby <n>   (0 <= n <= 25)
//
// The Luby sequence for parameter n is defined recursively:
//   luby(0) = "a"
//   luby(k) = luby(k-1) luby(k-1) char('a'+k)

#include <cstdlib>
#include <cstdio>
#include <string>

static void writeLuby(FILE *f, int k) {
    if (k == 0) {
        fputc('a', f);
        return;
    }
    writeLuby(f, k - 1);
    writeLuby(f, k - 1);
    fputc('a' + k, f);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: gen_luby <n>  (0 <= n <= 25)\n");
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 0 || n > 25) {
        fprintf(stderr, "Error: n must be between 0 and 25\n");
        return 1;
    }
    long long m = 1LL << n;
    std::string filename = "luby" + std::to_string(m) + ".txt";
    FILE *f = fopen(filename.c_str(), "wb");
    if (!f) {
        fprintf(stderr, "Cannot open %s for writing\n", filename.c_str());
        return 1;
    }
    writeLuby(f, n);
    fclose(f);
    fprintf(stderr, "Wrote %s\n", filename.c_str());
    return 0;
}
