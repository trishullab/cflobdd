.PHONY : all 
DIR = /users/meghana/desktop/meghana/research/cudd/cudd-3.0.0
C_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-complex
B_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-big
C_B_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-complex-big
F_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-fourier
REL_DIR = ../cudd-3.0.0
REL_C_DIR = ../cudd-complex
REL_B_DIR = ../cudd-big
REL_C_B_DIR = ../cudd-complex-big
REL_F_DIR = ../cudd-fourier


tmp:
	g++ -DHAVE_CONFIG_H -I./$(DIR)  -I./$(DIR)/cudd -I./$(DIR)/mtr -I./$(DIR)/epd -I./$(DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT test.o -MD -MP -MF test.Tpo -c -o test.o `test -f 'test.cc' || echo './'`test.cc
	/bin/sh ./$(REL_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o test test.o $(DIR)/cplusplus/libobj.la $(DIR)/cudd/libcudd.la 

complex:
	gcc -DHAVE_CONFIG_H -I./$(C_DIR)  -I./$(C_DIR)/cudd -I./$(C_DIR)/mtr -I./$(C_DIR)/epd -I./$(C_DIR)/st -I./$(C_DIR)/util  -MT complex_test.o -MD -MP -MF complex_test.Tpo -c -o complex_test.o `test -f 'complex_test.c' || echo './'`complex_test.c
	/bin/sh ./$(REL_C_DIR)/libtool  --tag=CXX   --mode=link gcc  -Wall -Wextra -g -O3   -o complex_test complex_test.o $(C_DIR)/cudd/libcudd.la 

big:
	g++ -DHAVE_CONFIG_H -I./$(B_DIR)  -I./$(B_DIR)/cudd -I./$(B_DIR)/mtr -I./$(B_DIR)/epd -I./$(B_DIR)/st -I./Solver2/assert -I./Solver2/bit_vector -I./Solver2/matrix -I./Solver2/parsing  -Wall -Wextra -std=c++20 -g -O3 -MT  big_test.o -MD -MP -MF big_test.Tpo -c -o big_test.o `test -f 'big_test.cc' || echo './'`big_test.cc
	/bin/sh ./$(REL_B_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++20 -g -O3   -o big_test big_test.o Solver2/bit_vector/value.o Solver2/bit_vector/bit_vector_1.o Solver2/parsing/parsing.o $(B_DIR)/cplusplus/libobj.la $(B_DIR)/cudd/libcudd.la 
complex_big:
	g++ -DHAVE_CONFIG_H -I./$(C_B_DIR)  -I./$(C_B_DIR)/cudd -I./$(C_B_DIR)/mtr -I./$(C_B_DIR)/epd -I./$(C_B_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT complex_big_test.o -MD -MP -MF complex_big_test.Tpo -c -o complex_big_test.o `test -f 'complex_big_test.cc' || echo './'`complex_big_test.cc
	/bin/sh ./$(REL_C_B_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o complex_big_test complex_big_test.o $(C_B_DIR)/cplusplus/libobj.la $(C_B_DIR)/cudd/libcudd.la 

fourier:
	g++ -g -DHAVE_CONFIG_H -I./$(F_DIR)  -I./$(F_DIR)/cudd -I./$(F_DIR)/mtr -I./$(F_DIR)/epd -I./$(F_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT fourier_test.o -MD -MP -MF fourier_test.Tpo -c -o fourier_test.o `test -f 'fourier_test.cc' || echo './'`fourier_test.cc
	/bin/sh ./$(REL_F_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o fourier_test fourier_test.o $(F_DIR)/cplusplus/libobj.la $(F_DIR)/cudd/libcudd.la 


dep:
	cd ../cudd-3.0.0/ && make all && cd ../examples
	cd ../cudd-complex/ && make all && cd ../examples 
	cd ../cudd-big/ && make all && cd ../examples 
	cd ../cudd-complex-big/ && make all && cd ../examples 
	cd ../cudd-fourier/ && make all && cd ../examples 
all:  dep tmp complex big complex_big fourier
clean: 
	rm *.o *.Tpo test complex_test big_test complex_big_test fourier_test
