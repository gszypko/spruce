.PHONY: clang intel gnu clean unspec

unspec:
	@echo 'Please specify which compiler to use, from {clang, intel, gnu}.'

clang:
	clang++ -I./source -I./source/user-interface -I./source/mhd -I./source/solar -I./source/ucnp\
		./source/mhd/*.cpp ./source/solar/*.cpp ./source/ucnp/*.cpp ./source/user-interface/*.cpp ./main.cpp\
		-Xpreprocessor -fopenmp -lomp -std=c++17 -o ./run

intel:
	icc -I./source -I./source/user-interface -I./source/mhd -I./source/solar -I./source/ucnp\
		./main.cpp ./source/mhd/*.cpp ./source/solar/*.cpp ./source/ucnp/*.cpp ./source/user-interface/*.cpp\
		-fopenmp -std=c++17 -o ./run

gnu:
	g++ -I./source -I./source/user-interface -I./source/mhd -I./source/solar -I./source/ucnp\
		./main.cpp ./source/mhd/*.cpp ./source/solar/*.cpp ./source/ucnp/*.cpp ./source/user-interface/*.cpp\
		-fopenmp -std=c++17 -o ./run

clean:
	rm mhdtoy
	
