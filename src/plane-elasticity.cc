#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>     
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>


#include "driver.hh"

int main(int argc, char** argv){
	Dune::MPIHelper::instance(argc, argv);

	// read input file
	Dune::ParameterTree input_data;
	std::string filename (std::string(argv[0])+".input");

	if (argc > 1){
		filename = argv[1];
	}
	try{
		Dune::ParameterTreeParser::readINITree (filename, input_data);
	}
	catch (...){
		std::cerr << "The configuration file \"" << filename << "\" "
					 "could not be read. Exiting..." << std::endl;
		std::exit(1);
	}

	int   level   =  input_data.get<int>("level");  // refine level
	double E      =  input_data.get<double>("E");   // Young Modulus
	double nu     =  input_data.get<double>("nu");  // Poisson Ratio
	double g_vert =  input_data.get<double>("g_vert");
	double rho    =  input_data.get<double>("rho");  // Mass Density
	std::string name = input_data.get<std::string>("output"); 

	//computing Lame constants lambda and mui
	double numerator = E * nu;
	double denominator = (1 + nu) * (1 - 2*nu);
	double lambda = numerator / denominator;
	
	denominator = 2 * (1 + nu);
	double mu = E / denominator;

	g_vert *= (mu + lambda);

	constexpr int dim = 2;  // grid dimension
	using GridType = Dune::YaspGrid<dim>;
	Dune::FieldVector<GridType::ctype,dim> L(2.0); // 

	std::array<int,dim> s = {20, 20};
	GridType grid(L, s);
	if(level > 0){
		grid.globalRefine(level);
	}

	auto gv = grid.leafGridView();
	driver(gv, E, nu, g_vert, rho, name);

	return 0;
}
