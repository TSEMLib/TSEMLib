/*
 * solve Poisson problem on [-1, 1]^2 with Dirichlet boundary condtions
 * use MovingMesh3D in AFEPack for r-adaptive
 * for gradient g, adopt 1/\sqrt{1+|g|^2} as indicator
 */

#include "Poisson.h"

int main(int argc, char* argv[])
{
    if (argc != 4){
	std::cout << "Usage: "
		  << argv[0]
		  << " filename_mesh polynomial_order move_step_number"
		  << std::endl;
	return 1;
    }
    try{
	Poisson poisson(argv[1], atoi(argv[2]), atoi(argv[3]));
	poisson.run();
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
