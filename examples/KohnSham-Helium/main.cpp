/*
 * solve ground state Hellium following the Kohn-Sham model, with LDA given by NIST database and Hartree atomic unit
 * use initial state \psi = \exp(-2r)/\|\exp(-2r)\|_2
 * consider energy error for break criterion
 * use multipole expansion when calculating Hartree potential
 * use LOBPCG method for linear eigenvalue problem, with preconditioner from free electron
 */

#include "SCFIteration.h"

int main(int argc, char *argv[])
{
    if (argc != 2){
	std::cout << "Usage: "
		  << argv[0]
		  << " filename_evolveInfo"
		  << std::endl;
	return 1;
    }
    try{
	SCFIteration scf_iteration(argv[1]);
	scf_iteration.run();
    }
    catch(std::exception& e) {
	std::cerr << "Exception caughted:" << std::endl
		  << e.what ()
		  << std::endl;
    }
    catch(...) {
	std::cerr << "Exception caughted:" << std::endl
		  << "unknown exception caughted."
		  << std::endl;
    }
    
    return 0;
}
