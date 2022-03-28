template<typename GV,typename PTree>
class HydraulicProperties
{
private:
	  const GV& gv;
	  const PTree& ptree;

	  const static int dim = GV::dimension;

	  const static int numOfParams  = 5;
	  const static int id_Pentry 	= 0;
	  const static int id_lambda 	= 1;
	  const static int id_Swr 		= 2;
	  const static int id_Sgr 		= 3;
	  const static int id_beta		= 4;

	  Parameters<PTree> parameter;
	  Soil<GV,PTree> soil;
	  CharacteristicValues Xc;
	  double x_c;
	  double P_c;

public:

  //! construct from grid view
  HydraulicProperties (const GV& gv_,const PTree& ptree_)
: gv( gv_ ),
  ptree(ptree_),
  parameter(ptree_),
  soil(gv_,ptree_)
{
	  x_c = Xc.x_c;
	  P_c = Xc.P_c;
}

	/* PARAMETERS FOR THE HYDRAULIC PROPERTY CONSTITUTIVE LAW (BROOKS-COREY) */
	std::vector<double>
	BrooksCoreyParameters( const typename GV::Traits::template Codim<0>::Entity& element,
	   	   	 	 	 	   const Dune::FieldVector<double,dim>& xlocal ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		std::vector<double> BCParams (numOfParams,0.);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();

		if( x[dim-1]>0. and x[dim-1]<z_L[1]/x_c ){
			BCParams[id_Pentry] = prop_L[0][2] ; /*Pa*/
			BCParams[id_lambda] = prop_L[0][3] ;
			BCParams[id_Swr] 	= prop_L[0][4] 	 ;
			BCParams[id_Sgr] 	= prop_L[0][5] 	 ;
			BCParams[id_beta] 	= prop_L[0][6] 	 ;
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1]/x_c and x[dim-1]<z_L[i]/x_c ){
				BCParams[id_Pentry] = prop_L[i][2] ; /*Pa*/
				BCParams[id_lambda] = prop_L[i][3] ;
				BCParams[id_Swr] 	= prop_L[i][4] 	 ;
				BCParams[id_Sgr] 	= prop_L[i][5] 	 ;
				BCParams[id_beta] 	= prop_L[i][6] 	 ;
			}
		}

		return BCParams;
	}

	/* EFFECTIVE SATURATION */

	double EffectiveSw( double Sw,
						double Swr,
						double Sgr ) const {


		double Sw_max = 1. - Sgr;
		double Swe = std::max(( Sw - Swr )/( Sw_max - Swr ),0.);
		
		return Swe;
	}

	double dSwe_dSw( double Sw,
					 double Swr,
					 double Sgr) const {

		double dSwe =  1./( 1. - Sgr - Swr );

		return dSwe;
	}


	/* SUCTION/CAPILLARY PRESSURE */

	double CapillaryPressure(  const typename GV::Traits::template Codim<0>::Entity& element,
  	   	   	 	 	 	 	   const Dune::FieldVector<double,dim>& xlocal ,
							   double Sw,
							   double porosity) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);

		double Pentry/*Pa*/	= BCParams[id_Pentry];
		double lambda 		= BCParams[id_lambda];
		double Sgr 			= BCParams[id_Sgr];
		double Swr 			= BCParams[id_Swr];
		double beta 		= BCParams[id_beta];

		double eta = (1/lambda);
		double Swe = EffectiveSw( Sw,Swr,Sgr );
		double Pc /*Pa*/ = std::min( Pentry * pow( Swe, -eta ) , 1.e8);

		double porosity_0 = soil.SedimentPorosity( element,xlocal );
		double SF_por = PcSF( porosity, porosity_0, beta );

		Pc *=SF_por; /*Pa*/
		return Pc/P_c; /*ndim*/
	}

	/* Pc SCALING FACTORS */

	double PcSF( double phi,
				 double phi_0,
				 double beta) const {

		double term = ( phi_0/phi ) * ( ( 1-phi )/( 1. - phi_0 ));
		double pcsf =pow( term, beta );

		return pcsf ;
	}

	/* RELATIVE PERMEABILITIES */

	double krw( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	const Dune::FieldVector<double,dim>& xlocal ,
				double Sw ) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);
		double lambda 	= BCParams[id_lambda];
		double Sgr 		= BCParams[id_Sgr];
		double Swr 		= BCParams[id_Swr];

		double Swe = EffectiveSw(Sw,Swr,Sgr);

		double kr = std::pow(Swe, (2.0/lambda + 3.0) );
		if( Swe>1.){
			kr=1.;
		}
		if( Swe<0.){
			kr=0.;
		}

		return kr ;
	}

	double krg( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	const Dune::FieldVector<double,dim>& xlocal ,
				double Sw ) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);
		double lambda 	= BCParams[id_lambda];
		double Sgr 		= BCParams[id_Sgr];
		double Swr 		= BCParams[id_Swr];

		double Swe = EffectiveSw(Sw,Swr,Sgr);

		double kr = std::pow(1.0-Swe, 2.0) * ( 1.0 - std::pow(Swe, (2.0/lambda + 1.0) ) );

		if( Swe>1.){
			kr=0.;
		}
		if( Swe<0.){
			kr=1.;
		}
		return kr;
	}

	/* PERMEABILITY SCALING FACTORS */

	double PermeabilityScalingFactor( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	 	 	 	 	 	  const Dune::FieldVector<double,dim>& xlocal ,
									  double porosity ) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);
		double Pentry 	= BCParams[id_Pentry];
		double lambda 	= BCParams[id_lambda];
		double Sgr 		= BCParams[id_Sgr];
		double Swr 		= BCParams[id_Swr];
		double beta 	= BCParams[id_beta];

		double porosity_0 = soil.SedimentPorosity( element,xlocal );
		double SF_por = KSF( porosity, porosity_0, beta );

		double SF =SF_por;
		return SF;
	}

	double KSF( double phi,
				double phi_0,
				double beta ) const {
		// POWER LAW MODEL PROPOSED BY CIVAN (2001)
		// Read " kinetic simulation of methane hydrate formation and issociation in porous media " by Xuefei Sun and Kishore Mohanty
		// phi_0 = basePorosity_initial (and NOT sediemntPorosity!! )
		// phi is basePorosity at current time

		double term1 = phi / phi_0 ;
		double term2 = ( 1-phi_0 ) / ( 1. - phi ) ;

		double ksf = term1 * std::pow( ( term1 * term2 ) , 2*beta );

		return ksf;
	}

	/*SCALED POROSITY */

	double ScaledPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double porosity0, double dPw /*(Pw-Pw0)*/ /*Pa*/ ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		double por = porosity0;
		por += soil.SedimentCompressibilityFactor(element,xlocal,dPw)/P_c * dPw;

		return por;
	}

//	double ScaledPorosity
//	(const typename GV::Traits::template Codim<0>::Entity& element,
//	 const Dune::FieldVector<double,dim>& xlocal,
//	 double porosity0, double dPw_dt /*(Pw-Pw0)/dt*/ /*Pa/s*/ ) const {
//
//		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);
//
//		double por = porosity0;
//		por += soil.SedimentCompressibilityFactor(element,xlocal,dPw_dt)/P_c * dPw_dt;
//
//		return por;
//	}


  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
