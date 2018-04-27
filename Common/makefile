CXX = icc
CXXFLAGS += -O2 -std=c++11

libcomn.a: comn.o
	ar rcs libcomn.a comn.o

comn.o: comn.cpp comn.h
	$(CXX) $(CXXFLAGS) -c comn.cpp

clean:
	rm *o
