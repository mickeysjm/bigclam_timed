cc = g++
std = c++11
exe = bigclam.exe
objects = .\utils.o .\graph.o .\bigclam.o .\model.o

$(exe) : $(objects)
			$(cc) -O2 -fopenmp -std=$(std) -o $(exe) $(objects)

.\utils.o : .\utils.cpp common.h utils.h
	$(cc) -std=$(std) -O2 -fopenmp -c .\utils.cpp -o .\utils.o
.\graph.o : .\graph.cpp graph.h
	$(cc) -std=$(std) -O2 -fopenmp -c .\graph.cpp -o .\graph.o
.\bigclam.o : .\bigclam.cpp graph.h model.h utils.h
	$(cc) -std=$(std) -O2 -fopenmp -c .\bigclam.cpp -o .\bigclam.o
.\model.o : .\model.cpp graph.h utils.h common.h model.h
	$(cc) -std=$(std) -O2 -fopenmp -c .\model.cpp -o .\model.o

clean :
		del $(objects) $(exe)
