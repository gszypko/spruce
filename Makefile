.PHONY: clang intel gnu clean unspec

unspec:
	@echo 'Please specify which compiler to use, from {clang, intel, gnu}.'

clang:
	clang++ -I./source/mhd -I./source/solar ./main.cpp ./source/mhd/*.cpp ./source/solar/*.cpp -Xpreprocessor -fopenmp -lomp -std=c++17 -o ./main

intel:
	icc -I./source/mhd -I./source/solar ./main.cpp ./source/mhd/*.cpp ./source/solar/*.cpp -fopenmp -std=c++17 -o ./main

gnu:
	g++ -I./source/mhd -I./source/solar ./main.cpp ./source/mhd/*.cpp ./source/solar/*.cpp -fopenmp -std=c++17 -o mhdtoy

clean:
	rm mhdtoy
	
