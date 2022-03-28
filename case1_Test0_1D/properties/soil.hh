template<typename GV, typename PTree>
class Soil
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  Parameters<PTree> parameter;

	  CharacteristicValues characteristicValue;
	  double xc;
	  double Pc;
	  double Kc;
	  double rhoc;
	  const static int dim = GV::dimension;

public:

  //! construct from grid view
  Soil (const GV& gv_, const PTree& ptree_)
  : gv( gv_ ), ptree(ptree_), parameter(ptree_)
  {
	  xc=characteristicValue.x_c;
	  Pc=characteristicValue.P_c;
	  Kc=characteristicValue.permeability_c;
	  rhoc=characteristicValue.density_c;
  }


	double SedimentPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.z_layers();//in m
		auto prop_L = parameter.layer_properties();

		double por = 0.;

		if( x[dim-1]>0. and x[dim-1]<z_L[1]/xc ){
			por = prop_L[0][0];
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1]/xc and x[dim-1]<z_L[i]/xc ){
				por = prop_L[i][0];
			}
		}

		return por;
	}

	double SedimentCompressibilityFactor
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double dP /*Pa*/ ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();

		double compressibility = 0.;

		if( x[dim-1]>0. and x[dim-1]<z_L[1]/xc ){
			if(dP>0)
				compressibility = prop_L[0][7];
			else
				compressibility = prop_L[0][8];
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1]/xc and x[dim-1]<z_L[i]/xc ){
				if(dP>0)
					compressibility = prop_L[i][7];
				else
					compressibility = prop_L[i][8];
			}
		}

		return compressibility*Pc;
	}

	double SedimentPermeability
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();

		double K = 0.; /*m^2*/

		if( x[dim-1]>0. and x[dim-1]<z_L[1]/xc ){
			K = prop_L[0][1];
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1]/xc and x[dim-1]<z_L[i]/xc ){
				K = prop_L[i][1];
			}
		}

		return K/Kc; /*ndim*/
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
		return rho/rhoc; /*ndim*/
	}

	double Tortuosity( double porosity ) const {
		return porosity * porosity ;/*ndim*/
	}

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
