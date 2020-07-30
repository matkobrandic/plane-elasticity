#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/uggrid.hh>     
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/grid/io/file/gmshreader.hh>



#include "driver.hh"

int main(int argc, char** argv){
	Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

	const int dim = 2;
	using GridType = Dune::UGGrid<dim>;
	using GridView = GridType::LeafGridView;
	bool verbosity = true;
	bool insertBoundarySegments = false;  // Bez toga Dune::GmshReader zna podbaciti (barem u 3D)
	// mashgrid filename
	std::string file("domain.msh");
	// Read mashgrid from msh file
	GridType* pgrid = Dune::GmshReader<GridType>::read(file, verbosity, insertBoundarySegments);
	GridView gv = pgrid->leafGridView();
	//In case of parallel computing -> distribute grid 
	pgrid->loadBalance();

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

	// int   level   =  input_data.get<int>("level");  // refine level
	double rho    =  input_data.get<double>("rho");  // Mass Density
	std::string name = input_data.get<std::string>("output");
	double E = input_data.get<double>("E");	// Young Modulus
	double nu = input_data.get<double>("nu");	// Poisson Ratio
	
	double mu = E / (2 * (1 + nu));
	double lambda = E * nu / ((1 + nu) * (1 - 2*nu));

	double g = 0.02 * (mu + lambda); 
	std::cout << "---------- Entering Driver routine ----------" << std::endl;
	driver(gv, mu, lambda, g, rho, name);

	return 0;
}