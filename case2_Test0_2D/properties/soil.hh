template<typename GV, typename PTree>
class Soil
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  Parameters<PTree> parameter;

	  CharacteristicValues characteristicValue;
	  const static int dim = GV::dimension;

public:

  //! construct from grid view
  Soil (const GV& gv_, const PTree& ptree_)
  : gv( gv_ ), ptree(ptree_), parameter(ptree_)
  {}


	double SedimentPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();
		auto isRandom = parameter.layer_property_distribution_flag();

		double por = 0.;

		const typename GV::IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(element);
		srand (cell_number);
		if( x[dim-1]>0. and x[dim-1]<z_L[1] ){
			por = prop_L[0][0];
			if( isRandom[0] ){
				double max = prop_L[0][0] * 1.;
				double min = prop_L[0][0] * 0.5;
				por = ((double) rand() / (RAND_MAX+1.0)) * (max-min) + min;
			}
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1] and x[dim-1]<z_L[i] ){
				por = prop_L[i][0];
				if( isRandom[i] ){
					double max = prop_L[i][0] * 1.;
					double min = prop_L[i][0] * 0.5;
					por = ((double) rand() / (RAND_MAX+1.0)) * (max-min) + min;
				}
			}
		}

		return por;
	}

	double SedimentCompressibilityFactor
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double dP ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();

		double compressibility = 0.;
		double denom = 500.;

		double dP0plus=0.;//1200.
		double dP0minus= -dP0plus;

		if( x[dim-1]>0. and x[dim-1]<z_L[1] ){
			auto C_posdP = prop_L[0][7];
			auto C_negdP = prop_L[0][8];
//			compressibility = C_negdP + (C_posdP-C_negdP) * (1./(1.+exp(-dP/denom)));
			if(dP>dP0plus){
				compressibility = C_posdP;
//				std::cout << x[0] << "," << x[1] << '\t' << dP << '\t' << compressibility << std::endl;
			}else if(dP<dP0minus){
				compressibility = C_negdP;
//				std::cout << x[0] << "," << x[1] << '\t' << dP << '\t' << compressibility << std::endl;
			}
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1] and x[dim-1]<z_L[i] ){
				auto C_posdP = prop_L[i][7];
				auto C_negdP = prop_L[i][8];
//				compressibility = C_negdP + (C_posdP-C_negdP) * (1./(1.+exp(-dP/denom)));
				if(dP>dP0plus){
					compressibility = C_posdP;
//					std::cout<< x[0] << "," << x[1] << '\t' << dP << '\t' << compressibility << std::endl;
				}else if(dP<dP0minus){
					compressibility = C_negdP;
//					std::cout<< x[0] << "," << x[1] << '\t' << dP << '\t' << compressibility << std::endl;
				}
			}
		}

		return compressibility;
	}

	double SedimentPermeability
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();
		auto isRandom = parameter.layer_property_distribution_flag();

		double K = 0.; /*m^2*/

		const typename GV::IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(element);
		srand (cell_number);
		if( x[dim-1]>0. and x[dim-1]<z_L[1] ){
			K = prop_L[0][1];
			if( isRandom[0] ){
//				double Kmax = prop_L[0][1] * 1.;
//				double Kmin = prop_L[0][1] * 0.1;
//				K = ((double) rand() / (RAND_MAX+1.0)) * (Kmax-Kmin) + Kmin;
				auto por = SedimentPorosity(element,xlocal);
				auto por0= prop_L[0][0];
				K *= pow( por/por0 , 3. );
			}
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1] and x[dim-1]<z_L[i] ){
				K = prop_L[i][1];
				if( isRandom[i] ){
//					double Kmax = prop_L[i][1] * 1.;
//					double Kmin = prop_L[i][1] * 0.1;
//					K = ((double) rand() / (RAND_MAX+1.0)) * (Kmax-Kmin) + Kmin;
					auto por = SedimentPorosity(element,xlocal);
					auto por0= prop_L[i][0];
					K *= pow( por/por0 , 3. );
				}
			}
		}

		return K/characteristicValue.permeability_c; /*ndim*/
	}

	// vector coefficient
	Dune::FieldVector<double,dim>
	SedimentPermeabilityVector
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		double K_xx = SedimentPermeability(element,xlocal);
		double K_yy = K_xx;
		Dune::FieldVector<double,dim> PermeabilityVector;

		PermeabilityVector[0] = K_xx ;
		PermeabilityVector[1] = K_yy ;

		return PermeabilityVector; /*ndim*/
	}

	double Density() const {
		/* unit -> kg/m^3 */
		double rho = 2600.0;
		return rho/characteristicValue.density_c; /*ndim*/
	}

	double Tortuosity( double porosity ) const {
		return porosity * porosity ;/*ndim*/
	}

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
