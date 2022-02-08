.PHONY : all 
DIR = /users/meghana/desktop/meghana/research/cudd/cudd-3.0.0
C_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-complex
B_DIR = /users/meghana/desktop/meghana/research/cudd/cudd-big
REL_DIR = ../cudd-3.0.0
REL_C_DIR = ../cudd-complex
REL_B_DIR = ../cudd-big


tmp:
	g++ -DHAVE_CONFIG_H -I./$(DIR)  -I./$(DIR)/cudd -I./$(DIR)/mtr -I./$(DIR)/epd -I./$(DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT test.o -MD -MP -MF test.Tpo -c -o test.o `test -f 'test.cc' || echo './'`test.cc
	/bin/sh ./$(REL_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o test test.o $(DIR)/cplusplus/libobj.la $(DIR)/cudd/libcudd.la 

complex:
	gcc -DHAVE_CONFIG_H -I./$(C_DIR)  -I./$(C_DIR)/cudd -I./$(C_DIR)/mtr -I./$(C_DIR)/epd -I./$(C_DIR)/st -I./$(C_DIR)/util  -MT complex_test.o -MD -MP -MF complex_test.Tpo -c -o complex_test.o `test -f 'complex_test.c' || echo './'`complex_test.c
	/bin/sh ./$(REL_C_DIR)/libtool  --tag=CXX   --mode=link gcc  -Wall -Wextra -g -O3   -o complex_test complex_test.o $(C_DIR)/cudd/libcudd.la 

big:
	g++ -DHAVE_CONFIG_H -I./$(B_DIR)  -I./$(B_DIR)/cudd -I./$(B_DIR)/mtr -I./$(B_DIR)/epd -I./$(B_DIR)/st   -Wall -Wextra -std=c++0x -g -O3 -MT big_test.o -MD -MP -MF big_test.Tpo -c -o big_test.o `test -f 'big_test.cc' || echo './'`big_test.cc
	/bin/sh ./$(REL_B_DIR)/libtool  --tag=CXX   --mode=link g++  -Wall -Wextra -std=c++0x -g -O3   -o big_test big_test.o $(B_DIR)/cplusplus/libobj.la $(B_DIR)/cudd/libcudd.la 

dep:
	cd ../cudd-3.0.0/ && make all && cd ../examples
	cd ../cudd-complex/ && make all && cd ../examples 
	cd ../cudd-big/ && make all && cd ../examples 
all:  dep tmp complex big
clean: 
	rm *.o *.Tpo test complex_test big_test
