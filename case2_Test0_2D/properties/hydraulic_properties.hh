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
	  CharacteristicValues characteristicValue;

public:

  //! construct from grid view
  HydraulicProperties (const GV& gv_,const PTree& ptree_)
: gv( gv_ ),
  ptree(ptree_),
  parameter(ptree_),
  soil(gv_,ptree_)
{}

	/* PARAMETERS FOR THE HYDRAULIC PROPERTY CONSTITUTIVE LAW (BROOKS-COREY) */
	std::vector<double>
	BrooksCoreyParameters( const typename GV::Traits::template Codim<0>::Entity& element,
	   	   	 	 	 	   const Dune::FieldVector<double,dim>& xlocal ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		std::vector<double> BCParams (numOfParams,0.);

		auto z_L = parameter.z_layers();
		auto prop_L = parameter.layer_properties();
		auto isRandom = parameter.layer_property_distribution_flag();

		double K = 0.; /*m^2*/

		const typename GV::IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(element);
		srand (cell_number);
		if( x[dim-1]>0. and x[dim-1]<z_L[1] ){
			BCParams[id_Pentry] = prop_L[0][2] ; /*Pa*/
			BCParams[id_lambda] = prop_L[0][3] ;
			BCParams[id_Swr] 	= prop_L[0][4] 	 ;
			BCParams[id_Sgr] 	= prop_L[0][5] 	 ;
			BCParams[id_beta] 	= prop_L[0][6] 	 ;
			if( isRandom[0] ){
//				double Pemax = prop_L[0][2] * 1.;
//				double Pemin = prop_L[0][2] * 0.01;
//				BCParams[id_Pentry] = ((double) rand() / (RAND_MAX+1.0)) * (Pemax-Pemin) + Pemin;
				auto por = soil.SedimentPorosity(element,xlocal);
				auto por0= prop_L[0][0];
				BCParams[id_Pentry] *= pow( por0/por , 2. );
			}
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i-1] and x[dim-1]<z_L[i] ){
				BCParams[id_Pentry] = prop_L[i][2] ; /*Pa*/
				BCParams[id_lambda] = prop_L[i][3] ;
				BCParams[id_Swr] 	= prop_L[i][4] 	 ;
				BCParams[id_Sgr] 	= prop_L[i][5] 	 ;
				BCParams[id_beta] 	= prop_L[i][6] 	 ;
				if( isRandom[i] ){
//					double Pemax = prop_L[i][2] * 1.;
//					double Pemin = prop_L[i][2] * 0.01;
//					BCParams[id_Pentry] = ((double) rand() / (RAND_MAX+1.0)) * (Pemax-Pemin) + Pemin;
					auto por = soil.SedimentPorosity(element,xlocal);
					auto por0= prop_L[i][0];
					BCParams[id_Pentry] *= pow( por0/por , 2. );
				}
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
		double Pc = 0.; /*Pa*/
		double a = 0.05 ;

		Pc = std::min( Pentry * pow( Swe, -eta ) , 1.e8);

//		if( Swe > a ){
//			Pc = Pentry * pow( Swe, -eta );
//		}
//		else if ( Swe <= a ){
//			double Pc_a /*Pa*/  = Pentry * pow( a, -eta );
//			double dPc_a /*Pa*/ = dPc_dSwe( a,Pentry,lambda ) ;
//			Pc/*Pa*/ = Pc_a/*Pa*/ + dPc_a/*Pa*/ * ( Swe - a );
//		}
//		else {
//			std::cout<< " ERROR in " << __FILE__
//					 << " function: SuctionPressure( element,xlocal,Sw,porosity )"
//					 << "  , Swe = " << Swe
//					 << "  , Sw  = " << Sw
//					 << "  , por = " << porosity
//					 << "  , Pc  = " << Pc << std::endl;
////			exit(0);
//		}
//
//		if( Pc < -1e-3 ){
//			std::cout<< " Pc is -ve " << std::endl;
//			std::cout<< " ERROR in " << __FILE__
//					 << " function: SuctionPressure( element,xlocal,Sw,porosity )"
//					 << "  , Swe = " << Swe
//					 << "  , Sw  = " << Sw
//					 << "  , por = " << porosity
//					 << "  , Pc  = " << Pc << std::endl;
////			exit(0);
//		}

		double porosity_0 = soil.SedimentPorosity( element,xlocal );
		double SF_por = PcSF( porosity, porosity_0, beta );

		Pc *=SF_por; /*Pa*/
		return Pc/characteristicValue.P_c; /*ndim*/
	}

	double dPc_dSwe( double Swe,
					 double Pentry, /*Pa*/
					 double lambda ) const {

		double eta = (1/lambda);
		double dPc = 0.; /*Pa*/
		double a = 0.05 ;

		if( Swe > a ){
			dPc/*Pa*/ = Pentry * (-1./lambda) * std::pow( Swe , -(1./lambda) - 1. );
		}
		else if ( Swe <= a ){
			double dPc_a  = Pentry * (-1./lambda) * std::pow( a , -(1./lambda) - 1. ) ;
			double ddPc_a = Pentry * (-1./lambda) * (-1./lambda-1.) * std::pow( a , -(1./lambda) - 2. );
			dPc/*Pa*/ = dPc_a + ddPc_a * ( Swe - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::dPc_dSwe( Swe, Pentry, lambda ) "
					 << "  , Swe = " << Swe
					 << "  , dPc  = " << dPc << std::endl;
//			exit(0);
		}

		return dPc; /*Pa*/
	}

	/* Pc SCALING FACTORS */

	double PcSF( double phi,
				 double phi_0,
				 double beta) const {

		double pcsf = 0.;

		double term = ( phi_0/phi ) * ( ( 1-phi )/( 1. - phi_0 ));
		pcsf =pow( term, beta );

//		double a = 0.05 ;
//		if ( phi > a ){
//			double term = ( phi_0/phi ) * ( ( 1-phi )/( 1. - phi_0 ));
//			pcsf =pow( term, beta );
//		}
//		else if( phi <= a ){
//			double term_a = ( phi_0/a ) * ( ( 1-a )/( 1. - phi_0 ));
//			double pcsf_a = std::pow( term_a,beta );
//			double dpcsf_a = 0.;
//			double C = std::pow( phi_0/( 1. - phi_0 ) , beta );
//			dpcsf_a -= beta * C * std::pow( 1.-a , beta-1. ) * std::pow( a , -beta-1. );
//			pcsf = pcsf_a + dpcsf_a * ( phi - a );
//		}
//		else {
//			std::cout<< " ERROR in " << __FILE__
//					 << " function: PcSF( phi, phi_0, beta )" << std::endl;
////			exit(0);
//		}

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

		double ksf=0.;

//		double a = 0.95 ;
//		if( phi < a ){
//			ksf = term1 * std::pow( ( term1 * term2 ) , 2*beta );
//		}
//		else if ( phi >= a ){
//			double C = std::pow( 1-phi_0 , 2*beta ) / std::pow( phi_0 , 2*beta +1 );
//			double ksf_a = C * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta );
//			double dksf_a =    C * ( 2*beta +1 ) * std::pow( a, 2*beta   ) * std::pow( 1.-a , -2*beta      )
//							 + C * ( -2*beta   ) * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta - 1. );
//			ksf = ksf_a + dksf_a * ( phi - a );
//		}
//		else {
//			std::cout<< " ERROR in " << __FILE__
//					 << " function: KSF2( phi, phi_0, beta )" << std::endl;
////			exit(0);
//		}

		ksf = term1 * std::pow( ( term1 * term2 ) , 2*beta );

		return ksf;
	}

	/*SCALED POROSITY */

	double ScaledPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal,
	 double porosity0, double dPw_dt /*(Sg*Pw-Sg0*Pw0)/dt*/ ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		double por = porosity0 * exp( soil.SedimentCompressibilityFactor(element,xlocal,dPw_dt) * dPw_dt );

//		double por = std::min(porosity0 + soil.SedimentCompressibilityFactor(element,xlocal,dPw) * dPw, 0.8 );
//		if(por<0.01){
//			por = 0.01;
//		}

		return por;
	}


  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
