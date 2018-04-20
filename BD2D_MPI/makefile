# choose compiler: icc or g++
CXX = icc
CXXFLAGS += -O2 -std=c++11 -I${NCDIR}/include -I${CMDLINE_DIR}
LDFLAGS += -L${NCDIR}/lib -L${H5DIR}/lib -L${CURLDIR} -Wl,-rpath=${H5DIR}/lib,-rpath=${CURLDIR}
# lib path
LIBS = -lnetcdf -lhdf5_hl -lhdf5 -lcurl -lm
# source
SRC = $(wildcard *.cpp)
# source object
SRC_OBJ = $(SRC:.cpp=.o)
# source binary
SRC_BIN =main.out

$(SRC_BIN): $(SRC_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm *o *.out
