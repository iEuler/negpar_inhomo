This project implements the negative particle method for collisional plasma described in the following work.

B. Yan, R. Caflisch, A Monte Carlo method with negative particles for Coulomb collisions, J. Comput. Phys., 298 (2015), pp. 711-740.

B. Yan, A hybrid method with deviational particles for spatial inhomogeneous plasma, J.Comput. Phys., 309 (2016), pp. 18-36.

To run the program under Linux, 

g++ *.cpp -o out -lfftw3 -std=c++17 -fopenmp

./out
