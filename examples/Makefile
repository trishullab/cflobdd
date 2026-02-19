.PHONY : all 
DIR = /home/ubuntu/Research/cudd-3.0.0
C_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-complex
B_DIR = /home/ubuntu/Research/cudd-big
C_B_DIR = /home/ubuntu/Research/cudd-complex
F_DIR = /home/ubuntu/CUDD/cudd-fourier
S_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-addmin
REL_DIR = ../cudd-3.0.0
REL_C_DIR = ../cudd-complex
REL_B_DIR = ../cudd-big
REL_C_B_DIR = ../cudd-complex
REL_F_DIR = ../cudd-fourier
REL_S_DIR = ../cudd-addmin


tmp:
	g++ -O3 -DHAVE_CONFIG_H -I./$(DIR)  -I./$(DIR)/cudd -I./$(DIR)/mtr -I./$(DIR)/epd -I./$(DIR)/st   -Wall -Wextra  -Wl,--stack,4096 -std=c++0x -g -O3 -MT test.o -MD -MP -MF test.Tpo -c -o test.o `test -f 'test.cc' || echo './'`test.cc
	/bin/bash ./$(REL_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o test test.o $(DIR)/cplusplus/libobj.la $(DIR)/cudd/libcudd.la

multmod:
	g++ -O3 -DHAVE_CONFIG_H -I./$(REL_DIR)  -I./$(REL_DIR)/cudd -I./$(REL_DIR)/mtr -I./$(REL_DIR)/epd -I./$(REL_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT addMultiplication.o -MD -MP -MF addMultiplication.Tpo -c -o addMultiplication.o `test -f 'addMultiplication.cc' || echo './'`addMultiplication.cc
	/bin/bash ./$(REL_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o addMultiplication addMultiplication.o $(REL_DIR)/cplusplus/libobj.la $(REL_DIR)/cudd/libcudd.la

# Direct build target (bypasses libtool - works on Windows with mixed shell environments)
CUDD_SRCS = $(wildcard $(REL_DIR)/cudd/cudd*.c) $(REL_DIR)/mtr/mtrBasic.c $(REL_DIR)/mtr/mtrGroup.c \
            $(REL_DIR)/st/st.c $(REL_DIR)/epd/epd.c \
            $(REL_DIR)/util/cpu_stats.c $(REL_DIR)/util/cpu_time.c $(REL_DIR)/util/cstringstream.c \
            $(REL_DIR)/util/datalimit.c $(REL_DIR)/util/pathsearch.c $(REL_DIR)/util/pipefork.c \
            $(REL_DIR)/util/prtime.c $(REL_DIR)/util/safe_mem.c $(REL_DIR)/util/strsav.c \
            $(REL_DIR)/util/texpand.c $(REL_DIR)/util/ucbqsort.c
CUDD_OBJS = $(patsubst $(REL_DIR)/%.c, cudd_build/%.o, $(CUDD_SRCS))
CUDDOBJ_SRC = $(REL_DIR)/cplusplus/cuddObj.cc

CUDD_CFLAGS = -DHAVE_CONFIG_H -I$(REL_DIR) -I$(REL_DIR)/cudd -I$(REL_DIR)/mtr -I$(REL_DIR)/epd -I$(REL_DIR)/st -I$(REL_DIR)/util -Wall -Wextra -g -O3
CUDD_CXXFLAGS = $(CUDD_CFLAGS) -std=c++14

cudd_build:
	mkdir -p cudd_build/cudd cudd_build/mtr cudd_build/st cudd_build/epd cudd_build/util cudd_build/cplusplus

cudd_build/%.o: $(REL_DIR)/%.c | cudd_build
	gcc $(CUDD_CFLAGS) -c -o $@ $<

cudd_build/cplusplus/cuddObj.o: $(CUDDOBJ_SRC) | cudd_build
	g++ $(CUDD_CXXFLAGS) -c -o $@ $<

cudd_build/libcudd.a: $(CUDD_OBJS) cudd_build/cplusplus/cuddObj.o
	ar rcs $@ $^

multmod-direct: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication addMultiplication.o cudd_build/libcudd.a

# Boost include path for multiprecision types
BOOST_INC = -I../../boost_1_81_0

# Build targets for different bit-widths
multmod-32: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=32 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication32 addMultiplication.o cudd_build/libcudd.a

multmod-64: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=64 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication64 addMultiplication.o cudd_build/libcudd.a

multmod-128: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=128 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication128 addMultiplication.o cudd_build/libcudd.a

multmod-256: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=256 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication256 addMultiplication.o cudd_build/libcudd.a

multmod-512: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=512 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication512 addMultiplication.o cudd_build/libcudd.a

multmod-1024: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=1024 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication1024 addMultiplication.o cudd_build/libcudd.a

multmod-2048: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=2048 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication2048 addMultiplication.o cudd_build/libcudd.a

multmod-4096: cudd_build/libcudd.a
	g++ $(CUDD_CXXFLAGS) $(BOOST_INC) -DNUM_BITS=4096 -c -o addMultiplication.o addMultiplication.cc
	g++ $(CUDD_CXXFLAGS) -o addMultiplication4096 addMultiplication.o cudd_build/libcudd.a

MULT = addMultiplication32.exe
spec-all-test:
	./$(MULT) spec-all
	./$(MULT) spec-all
	./$(MULT) spec-all
	./$(MULT) spec-all

shiftadd-all-test:
	./$(MULT) shiftadd-all
	./$(MULT) shiftadd-all
	./$(MULT) shiftadd-all
	./$(MULT) shiftadd-all

karatsuba-all-test:
	./$(MULT) karatsuba-all
	./$(MULT) karatsuba-all
	./$(MULT) karatsuba-all
	./$(MULT) karatsuba-all

mult-tests:
	make spec-all-test MULT=addMultiplication32.exe
	make shiftadd-all-test MULT=addMultiplication32.exe
	make karatsuba-all-test MULT=addMultiplication32.exe
	make spec-all-test MULT=addMultiplication64.exe
	make shiftadd-all-test MULT=addMultiplication64.exe
	make spec-all-test MULT=addMultiplication128.exe
	make shiftadd-all-test MULT=addMultiplication128.exe
	make spec-all-test MULT=addMultiplication256.exe
	make shiftadd-all-test MULT=addMultiplication256.exe

complex:
	gcc -DHAVE_CONFIG_H -I./$(C_DIR)  -I./$(C_DIR)/cudd -I./$(C_DIR)/mtr -I./$(C_DIR)/epd -I./$(C_DIR)/st -I./$(C_DIR)/util  -MT complex_test.o -MD -MP -MF complex_test.Tpo -c -o complex_test.o `test -f 'complex_test.c' || echo './'`complex_test.c
	/bin/sh ./$(REL_C_DIR)/libtool  --tag=CXX   --mode=link gcc  -Wall -Wextra -g -O3   -o complex_test complex_test.o $(C_DIR)/cudd/libcudd.la 

big:
	g++ -O3 -DHAVE_CONFIG_H -I./$(B_DIR)  -I./$(B_DIR)/cudd -I./$(B_DIR)/mtr -I./$(B_DIR)/epd -I./$(B_DIR)/st -I./Solver2/assert -I./Solver2/bit_vector -I./Solver2/matrix -I./Solver2/parsing  -Wall -Wextra -std=c++2a -g -O3 -MT big_test.o -MD -MP -MF big_test.Tpo -c -o big_test.o `test -f 'big_test.cc' || echo './'`big_test.cc
	/bin/bash ./$(REL_B_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++2a -g -O3 -Wdeprecated-copy -lgmp -lmpfr -lgmpxx  -o big_test big_test.o Solver2/bit_vector/value.o Solver2/bit_vector/bit_vector_1.o Solver2/parsing/parsing.o $(B_DIR)/cplusplus/libobj.la $(B_DIR)/cudd/libcudd.la 
complex_big:
	g++ -O3 -DHAVE_CONFIG_H -I./$(C_B_DIR)  -I./$(C_B_DIR)/cudd -I./$(C_B_DIR)/mtr -I./$(C_B_DIR)/epd -I./$(C_B_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT complex_big_test.o -MD -MP -MF complex_big_test.Tpo -c -o complex_big_test.o `test -f 'complex_big_test.cc' || echo './'`complex_big_test.cc
	/bin/sh ./$(REL_C_B_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o complex_big_test complex_big_test.o $(C_B_DIR)/cplusplus/libobj.la $(C_B_DIR)/cudd/libcudd.la 

fourier:
	g++ -O3 -DHAVE_CONFIG_H -I./$(F_DIR)  -I./$(F_DIR)/cudd -I./$(F_DIR)/mtr -I./$(F_DIR)/epd -I./$(F_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT fourier_test.o -MD -MP -MF fourier_test.Tpo -c -o fourier_test.o `test -f 'fourier_test.cc' || echo './'`fourier_test.cc
	/bin/bash ./$(REL_F_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3 -o fourier_test fourier_test.o $(F_DIR)/cplusplus/libobj.la $(F_DIR)/cudd/libcudd.la -lgmp -lmpfr  
shortest_path:
	g++ -g -DHAVE_CONFIG_H -I./$(S_DIR)  -I./$(S_DIR)/cudd -I./$(S_DIR)/mtr -I./$(S_DIR)/epd -I./$(S_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT shortest_path.o -MD -MP -MF shortest_path_test.Tpo -c -o shortest_path_test.o `test -f 'shortest_path_test.cc' || echo './'`shortest_path_test.cc
	/bin/sh ./$(REL_S_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o shortest_path_test shortest_path_test.o $(S_DIR)/cplusplus/libobj.la $(S_DIR)/cudd/libcudd.la 


dep:
	cd ../cudd-3.0.0/ && make all && cd ../examples
	cd ../cudd-complex/ && make all && cd ../examples 
	cd ../cudd-big/ && make all && cd ../examples 
	cd ../cudd-complex-big/ && make all && cd ../examples 
	cd ../cudd-fourier/ && make all && cd ../examples 
	cd ../cudd-addMin/ && make all && cd ../examples 
all:  dep tmp complex big complex_big fourier shortest_path
clean: 
	rm *.o *.Tpo test complex_test big_test complex_big_test fourier_test shortest_path
