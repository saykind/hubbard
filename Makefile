CC		= clang
CXX		= clang++
CPPFLAGS= -Wall -lm -fopenmp -O2
LDFLAGS	= -lm 
SRC		= main.cpp sqrmtx.cpp qmc.cpp

all: build run

build: $(SRC)
	$(CXX) $(CPPFLAGS) $(SRC) -o a.out

run: 
	for beta in 4 2 1 .5 .25 .1 .05 .025 .01 .005 ; do \
		for u in .25 .5 1 2 4 6 8 10 12; do \
			./a.out u=$$u beta=$$beta ; \
		done ; \
	done
clean:
	rm *.out *.eps
