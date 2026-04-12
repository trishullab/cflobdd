CXX = g++
CXXFLAGS = -O3 -std=c++11

all: gen_repeat.exe gen_luby.exe

gen_repeat.exe: gen_repeat.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

gen_luby.exe: gen_luby.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

test-luby: gen_luby.exe
	rm -f testLuby.txt
	for n in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25; do \
		m=$$((1 << $$n)); \
		echo "n=$$n (luby$$m)"; \
		./gen_luby.exe $$n; \
		echo "=== n=$$n, m=$$m ===" >> testLuby.txt; \
		../testSequiturToADD.exe luby$$m.txt >> testLuby.txt 2>&1; \
		rm -f luby$$m.txt; \
	done
	@echo "Results in testLuby.txt"

test-luby-trace: gen_luby.exe
	rm -f testLubyTrace.txt
	for n in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25; do \
		m=$$((1 << $$n)); \
		echo "n=$$n (luby$$m)"; \
		./gen_luby.exe $$n; \
		echo "=== n=$$n, m=$$m ===" >> testLubyTrace.txt; \
		../testTraceToADD.exe luby$$m.txt >> testLubyTrace.txt 2>&1; \
		rm -f luby$$m.txt; \
	done
	@echo "Results in testLubyTrace.txt"

clean:
	rm -f gen_repeat.exe gen_luby.exe testLuby.txt testLubyTrace.txt

.PHONY: all clean test-luby test-luby-trace
