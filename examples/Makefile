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
