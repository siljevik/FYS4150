all: analytical MCMC_spin main 

analytical: analytical.cpp
	g++ -o analytical analytical.cpp

MCMC_spin: MCMC_spin.cpp
	g++ -o MCMC_spin MCMC_spin.cpp

main: main.cpp
	g++ -o main main.cpp -std=c++11

clean:
	rm analytical MCMC_spin main
