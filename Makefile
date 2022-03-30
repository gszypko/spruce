EXEC=run

gnu:
	g++ -I./source -I./source/user-interface -I./source/mhd -I./source/modules -I./source/solar -I./source/ucnp\
		./main.cpp ./source/modules/*.cpp ./source/mhd/*.cpp  ./source/solar/*.cpp ./source/ucnp/*.cpp ./source/user-interface/*.cpp\
		-fopenmp -std=c++17 -o ./$(EXEC)

clean:
	rm mhdtoy
	
