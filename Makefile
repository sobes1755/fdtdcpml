.POSIX:

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++26 -O3 -I./include/
# CXXFLAGS = -Wall -Wextra -std=c++26 -O3 -I./include/ -fsanitize=address
LDLIBS = -L. -lfdtdcpml -ltbb

AR = ar
ARFLAGS = rcs

SOLVER_OBJ = fdtdcpml.o fdtd.o cpml.o
SOLVER_LIB = libfdtdcpml.a

all: oamwave

$(SOLVER_LIB): $(SOLVER_OBJ)
	$(AR) $(ARFLAGS) $(SOLVER_LIB) $(SOLVER_OBJ)

oamwave: oamwave.cpp $(SOLVER_LIB)
	$(CXX) oamwave.cpp -o oamwave $(CXXFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $< -c -o $@ $(CXXFLAGS)

clean:
	-rm -f *.o oamwave $(SOLVER_LIB)

run:
	./oamwave 0.0 1.0 256 0.0 1.0 256 0.0 1.0 256 0.0 1E-6 128
