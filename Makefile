.PHONY: clang intel gnu clean unspec

unspec:
	@echo 'Please specify which compiler to use, from {clang, intel, gnu}.'

clang:
	clang++ -I ./source/mhd ./main.cpp ./source/mhd/*.cpp -Xpreprocessor -fopenmp -lomp -std=c++17 -o ./main

intel:
	icc -I ./source/mhd ./main.cpp ./source/mhd/*.cpp -fopenmp -std=c++17 -o ./main

gnu:
	g++ -I ./source/mhd ./main.cpp ./source/mhd/*.cpp -fopenmp -std=c++17 -o mhdtoy

clean:
	rm mhdtoy
	
