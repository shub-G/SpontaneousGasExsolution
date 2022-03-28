// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#define _USE_MATH_DEFINES
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<exception>
#include<chrono>

/*********************************************************************/
/* DEFAULT SETTINGS:
 * --
 */
#define USE_YASP
#define SEQUENTIAL
//#define FRAUENHOFER_MODEL

/* PROBLEM DESCRIPTION
 * 1D scenarios
 * 
 * 
 */

/*********************************************************************/

#include"include_dune.hh"
#include"case1_Test0_1D/include_problem.hh"

/*********************************************************************/

int main(int argc, char** argv)
{
	  try{
	    // Maybe initialize MPI
	    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
	    if(helper.size()==1){
	    std::cout << "This is case 1 of project Pockmark2P2C." << std::endl;
	    }
	    if(Dune::MPIHelper::isFake){
	      std::cout<< "This is a sequential program." << std::endl;
	    }
	    else {
	    	if(helper.size()==1){
	    		std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()<<" processes!"<<std::endl;
	    	}
	    }

		/**************************************************************************************************/
		// INPUTS
	    if (argc!=2)
	    {
	    	if(helper.rank()==0){
	    		std::cout << "usage: ./case1 <input_file.ini> " << std::endl;
	    	}
	        return 1;
	    }

	    // PARAMETER TREE

	    char input[40];
	    sscanf(argv[1],"%39s", input);
	    std::string input_file = "/sfs/fs1/work-geomar5/smomw325/dune_inputs/dune_2_8/Pockmark2P2C/case1_Test0_1D/";
	    input_file += input;
	    std::cout<< "input file: " << input_file << std::endl ;

	    Dune::ParameterTree ptree;
	    Dune::ParameterTreeParser ptreeparser;
	    ptreeparser.readINITree(input_file,ptree);
	    ptreeparser.readOptions(argc,argv,ptree);

		/**************************************************************************************************/
		// MESH
	    MeshParameters<Dune::ParameterTree> mesh(ptree);
	    const int dim = mesh.dimension;

#ifdef USE_YASP
		/*************
		 *  YASP
		 *************/
	    Dune::FieldVector<double,dim> L;
        L[0] = mesh.get_Xmax();
        L[1] = mesh.get_Zmax();
        std::array<int,dim> N;
        N[0] = mesh.get_nX();
        N[1] = mesh.get_nZ();

		std::bitset<dim> periodic(false);
		int overlap=1;
        typedef Dune::YaspGrid<dim> Grid;
        std::shared_ptr<Grid> grid = std::shared_ptr<Grid>(new Grid(L,N,periodic,overlap,helper.getCommunicator()));
        typedef Grid::LeafGridView GV;
        GV gv=grid->leafGridView();
        grid->loadBalance();

#elif defined(USE_UG)
		/*************
		 *  UG
		 *************/
		typedef std::vector<int> GmshIndexMap;
		GmshIndexMap boundary_index_map;
		GmshIndexMap element_index_map;
		typedef Dune::UGGrid<dim> GridType;
		GridType grid_type;
		const std::string grid_file = ptree.get("grid.ug.name") ;
		Dune::GmshReader<GridType> gmshreader;
		Dune::shared_ptr<GridType> grid(gmshreader.read(grid_file,boundary_index_map, element_index_map,true,false));

		typedef GridType::LeafGridView GV;
		GV gv = grid->leafGridView();
        grid->loadBalance();
#endif

		/**************************************************************************************************/
		// DRIVER
		driver( gv,ptree,helper);
		/**************************************************************************************************/

	  }
	  catch (Dune::Exception &e){
	    std::cerr << "Dune reported error: " << e << std::endl;
	  }
	  catch (...){
	    std::cerr << "Unknown exception thrown!" << std::endl;
	  }
}
