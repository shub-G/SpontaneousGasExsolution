template<typename GV,typename Properties>
class ProblemBoundaryConditions
{
private :
	const GV& gv ;
	const Properties& property;
	const static int dim = GV::dimension;

public :

	// ! construct from gridview
	ProblemBoundaryConditions (const GV& gv_,const Properties& property_)
	: gv ( gv_ ),
	  property(property_)
	{}
	
	/* NOTE:
	 * dirichlet: BC_w -> Pw, BC_g -> Sg, 
	 * neumann: total mass flux (convective+diffusive) BC_w -> f_w, BC_g -> f_g
	 * */
	
	/* boundary types */
	template<typename I> std::vector< int >
	type( I& intersection,
		  const Dune::FieldVector<double,dim-1>& xlocal,
		  double time/*s*/,
		  double dt/*s*/ ) const {
		
		auto xglobal = intersection.geometry().global(xlocal);
		
		std::vector< int > bct(Indices::numOfBCs);
		if( property.mesh.isTop(xglobal) ){
			bct[Indices::BC_w] = Indices::dirichlet;
			bct[Indices::BC_g] = Indices::dirichlet;
		}else{
			bct[Indices::BC_w] = Indices::neumann;
			bct[Indices::BC_g] = Indices::neumann;
		}
		
		return bct;
	}

	/* boundary values */
	template<typename I> std::vector< double >
	value ( I& intersection,
			const Dune::FieldVector<double,dim-1>& xlocal,
			double time/*s*/,
			double dt/*s*/ ) const {

		auto xglobal = intersection.geometry().global(xlocal);

		std::vector< double > bcv(Indices::numOfBCs,0.);
		if( property.mesh.isTop(xglobal) ){
			bcv[Indices::BC_w] = property.parameter.Pw_SF(xglobal,(time+dt));
			bcv[Indices::BC_g] = 0. ;
		}

		return bcv;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView () { return gv ; }
};
