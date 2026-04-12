// gen_repeat: write n copies of 'a' to a file (no spaces, no trailing newline)
// Usage: gen_repeat <n> <filename>

#include <cstdlib>
#include <cstdio>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: gen_repeat <n> <filename>\n");
        return 1;
    }
    int n = atoi(argv[1]);
    FILE *f = fopen(argv[2], "wb");
    if (!f) {
        fprintf(stderr, "Cannot open %s for writing\n", argv[2]);
        return 1;
    }
    for (int i = 0; i < n; i++)
        fputc('a', f);
    fclose(f);
    return 0;
}
