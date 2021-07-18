.PHONY: clang intel gnu clean unspec

unspec:
	@echo 'Please specify which compiler to use, from {clang, intel, gnu}.'

clang:
	clang++ -Xpreprocessor -fopenmp -lomp -std=c++11 grid.cpp utils.cpp derivs.cpp plasmadomain.cpp fileio.cpp evolution.cpp main.cpp -o mhdtoy

intel:
	icc grid.cpp utils.cpp derivs.cpp plasmadomain.cpp main.cpp -std=c++11 -fopenmp -o mhdtoy

gnu:
	g++ -fopenmp -std=c++11 grid.cpp utils.cpp derivs.cpp plasmadomain.cpp main.cpp -o mhdtoy

clean:
	rm mhdtoy
	
