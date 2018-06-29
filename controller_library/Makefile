CC=gcc
CXX=g++

BACKEND_OBJS=casadi_interface.o
BACKEND_OBJS+=sqp.o
BACKEND_OBJS+=regularize_hessian.o

AUTOGEN_SRCS=$(wildcard codegen/*.c)
IFACE_OBJS=$(AUTOGEN_SRCS:.c=.o)
IFACE_OBJS+=ocp.o
IFACE_OBJS+=transverse_ocp.o
IFACE_OBJS+=tangential_ocp.o
IFACE_OBJS+=duo_nmpc.o

INC=-I$(ACADOS)/include
INC+=-I$(ACADOS)/include/blasfeo/include
INC+=-I$(ACADOS)/include/hpipm/include
CFLAGS=-O2 -fPIC
CXXFLAGS=-O2 -fPIC -std=c++11

LDFLAGS=-L${ACADOS}/lib
LDFLAGS+=-rdynamic -lacados_c -lacore -lhpipm -lblasfeo -lqpOASES_e -lm -ldl

INSTALLDIR=../simulation/

all: libsqp_interface.so _sqp_interface.so

install: _sqp_interface.so
	cp -f _sqp_interface.so $(INSTALLDIR)
	cp -f sqp_interface.py $(INSTALLDIR)

libsqp_interface.so: $(BACKEND_OBJS) $(IFACE_OBJS)
	$(CXX) -o $@ -shared $^ $(CXXFLAGS) $(INC) $(LDFLAGS)

_sqp_interface.so: $(BACKEND_OBJS) $(IFACE_OBJS) sqp_interface_python.o
	$(CXX) -o $@ -shared $^ $(CXXFLAGS) $(INC) $(LDFLAGS)

sqp_interface_python.o: sqp_interface.i $(DEPS)
	swig -c++ -python -o sqp_interface_python.cxx sqp_interface.i
	$(CXX) -c -o $@ sqp_interface_python.cxx $(CXXFLAGS) $(INC) -I$(PYTHONINC)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC)

%.o: %.cxx $(DEPS)
	$(CC) -c -o $@ $< $(CXXFLAGS) $(INC)

clean:
	rm -rf *.o codegen/*.o *.so sqp_interface __pycache__ sqp_interface.py *_python.cxx *_python.h
