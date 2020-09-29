all: func

func: main.o task.o
	g++ -O3 main.o task.o  -pthread -o matrix

main.o: main.cpp 1.hpp
	g++ -c -O3  main.cpp -pthread

task.o: task.cpp 1.hpp
	g++ -c -O3  task.cpp -pthread

clean:
	rm -rf *.o func 
