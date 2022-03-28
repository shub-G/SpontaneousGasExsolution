/*********************************************************
 * EVALUATE OUTPUT VARIABLES
 *********************************************************/

template< class GV,
		  class Params,
		  class Evaluation_Pw,
		  class Evaluation_Sg,
		  class Evaluation_XCH4,
		  class Evaluation_YH2O,
		  class Evaluation_por0,
		  class Evaluation_Pw0,
		  class Evaluation_Sg0,
		  class GFS_PP, typename U_pp,
		  class GFS0,	typename U0>
class PostProcess{
private:
	const GV&  gv;
	const Params& param;
	Evaluation_Pw    *evaluation_Pw;
	Evaluation_Sg    *evaluation_Sg;
	Evaluation_XCH4  *evaluation_XCH4;
	Evaluation_YH2O  *evaluation_YH2O;
	Evaluation_por0  *evaluation_por_old;
	Evaluation_Pw0   *evaluation_pw_old;
	Evaluation_Sg0   *evaluation_sg_old;
	GFS_PP gfs_pp;
	U_pp *u_pp;
	GFS0 gfs0;
	U0 *u0_por_up;
	U0 *u0_pw_up;
	U0 *u0_sg_up;
	double *time;
	double *dt;

	double Xc_K ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_D;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;

//    typedef typename GV::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	PostProcess(	const GV& gv_,
					const Params& param_,
					Evaluation_Pw 	*evaluation_Pw_,
					Evaluation_Sg 	*evaluation_Sg_,
					Evaluation_XCH4 *evaluation_XCH4_,
					Evaluation_YH2O *evaluation_YH2O_,
					Evaluation_por0 *evaluation_por_old_,
					Evaluation_Pw0 	*evaluation_pw_old_,
					Evaluation_Sg0 	*evaluation_sg_old_,
					GFS_PP	 gfs_pp_,
					U_pp	*u_pp_,
					GFS0	 gfs0_,
					U0		*u0_por_up_,
					U0		*u0_pw_up_,
					U0		*u0_sg_up_,
					double *time_,
					double *dt_)
	: gv(gv_),
	  param(param_),
	  evaluation_Pw(evaluation_Pw_),
	  evaluation_Sg(evaluation_Sg_),
	  evaluation_XCH4(evaluation_XCH4_),
	  evaluation_YH2O(evaluation_YH2O_),
	  evaluation_por_old(evaluation_por_old_),
	  evaluation_pw_old(evaluation_pw_old_),
	  evaluation_sg_old(evaluation_sg_old_),
	  gfs_pp(gfs_pp_),
	  u_pp(u_pp_),
	  gfs0(gfs0_),
	  u0_por_up(u0_por_up_),
	  u0_pw_up(u0_pw_up_),
	  u0_sg_up(u0_sg_up_),
	  time(time_),
	  dt(dt_)
	{
		  Xc_K 		= param.characteristicValue.permeability_c;
		  Xc_mu 	= param.characteristicValue.viscosity_c;
		  Xc_rho 	= param.characteristicValue.density_c;
		  Xc_D		= param.characteristicValue.dispersivity_c;
		  Xc_P 		= param.characteristicValue.P_c;
		  Xc_T 		= param.characteristicValue.T_c;
		  Xc_t 		= param.characteristicValue.t_c;
		  Xc_x 		= param.characteristicValue.x_c;
	  }

	virtual ~PostProcess()
	{}

	void evaluate(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS_PP > LFS_PP;
		LFS_PP lfs_pp(gfs_pp);
		typedef typename LFS_PP::template Child<Indices::SVId_Pg>::Type LFS_PP_Pg;
		const LFS_PP_Pg& lfs_pp_Pg = lfs_pp.template child<Indices::SVId_Pg>();
		typedef typename LFS_PP::template Child<Indices::SVId_Pw>::Type LFS_PP_Pw;
		const LFS_PP_Pw& lfs_pp_Pw = lfs_pp.template child<Indices::SVId_Pw>();
		typedef typename LFS_PP::template Child<Indices::SVId_Pc>::Type LFS_PP_Pc;
		const LFS_PP_Pc& lfs_pp_Pc = lfs_pp.template child<Indices::SVId_Pc>();

		typedef typename LFS_PP::template Child<Indices::SVId_Sg>::Type LFS_PP_Sg;
		const LFS_PP_Sg& lfs_pp_Sg = lfs_pp.template child<Indices::SVId_Sg>();
		typedef typename LFS_PP::template Child<Indices::SVId_Sw>::Type LFS_PP_Sw;
		const LFS_PP_Sw& lfs_pp_Sw = lfs_pp.template child<Indices::SVId_Sw>();

		typedef typename LFS_PP::template Child<Indices::SVId_XCH4>::Type LFS_PP_XCH4;
		const LFS_PP_XCH4& lfs_pp_XCH4 = lfs_pp.template child<Indices::SVId_XCH4>();
		typedef typename LFS_PP::template Child<Indices::SVId_XH2O>::Type LFS_PP_XH2O;
		const LFS_PP_XH2O& lfs_pp_XH2O = lfs_pp.template child<Indices::SVId_XH2O>();
		typedef typename LFS_PP::template Child<Indices::SVId_YCH4>::Type LFS_PP_YCH4;
		const LFS_PP_YCH4& lfs_pp_YCH4 = lfs_pp.template child<Indices::SVId_YCH4>();
		typedef typename LFS_PP::template Child<Indices::SVId_YH2O>::Type LFS_PP_YH2O;
		const LFS_PP_YH2O& lfs_pp_YH2O = lfs_pp.template child<Indices::SVId_YH2O>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_rhow>::Type LFS_PP_rhow;
		const LFS_PP_rhow& lfs_pp_rhow = lfs_pp.template child<Indices::SVId_rhow>();
		typedef typename LFS_PP::template Child<Indices::SVId_rhog>::Type LFS_PP_rhog;
		const LFS_PP_rhog& lfs_pp_rhog = lfs_pp.template child<Indices::SVId_rhog>();

		typedef typename LFS_PP::template Child<Indices::SVId_K>::Type LFS_PP_K;
		const LFS_PP_K& lfs_pp_K = lfs_pp.template child<Indices::SVId_K>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_krw>::Type LFS_PP_krw;
		const LFS_PP_krw& lfs_pp_krw = lfs_pp.template child<Indices::SVId_krw>();
		typedef typename LFS_PP::template Child<Indices::SVId_krg>::Type LFS_PP_krg;
		const LFS_PP_krg& lfs_pp_krg = lfs_pp.template child<Indices::SVId_krg>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_muw>::Type LFS_PP_muw;
		const LFS_PP_muw& lfs_pp_muw = lfs_pp.template child<Indices::SVId_muw>();
		typedef typename LFS_PP::template Child<Indices::SVId_mug>::Type LFS_PP_mug;
		const LFS_PP_mug& lfs_pp_mug = lfs_pp.template child<Indices::SVId_mug>();

		typedef typename LFS_PP::template Child<Indices::SVId_zCH4>::Type LFS_PP_zCH4;
		const LFS_PP_zCH4& lfs_pp_zCH4 = lfs_pp.template child<Indices::SVId_zCH4>();

		typedef typename LFS_PP::template Child<Indices::SVId_por>::Type LFS_PP_por;
		const LFS_PP_por& lfs_pp_por = lfs_pp.template child<Indices::SVId_por>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_DH2O>::Type LFS_PP_DH2O;
		const LFS_PP_DH2O& lfs_pp_DH2O = lfs_pp.template child<Indices::SVId_DH2O>();
		typedef typename LFS_PP::template Child<Indices::SVId_DCH4>::Type LFS_PP_DCH4;
		const LFS_PP_DCH4& lfs_pp_DCH4 = lfs_pp.template child<Indices::SVId_DCH4>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_Pwsat>::Type LFS_PP_Pwsat;
		const LFS_PP_Pwsat& lfs_pp_Pwsat = lfs_pp.template child<Indices::SVId_Pwsat>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_HCH4>::Type LFS_PP_HCH4;
		const LFS_PP_HCH4& lfs_pp_HCH4 = lfs_pp.template child<Indices::SVId_HCH4>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_tau>::Type LFS_PP_tau;
		const LFS_PP_tau& lfs_pp_tau = lfs_pp.template child<Indices::SVId_tau>();
		

		typedef Dune::PDELab::LFSIndexCache<LFS_PP> LFSCache_PP;
		LFSCache_PP lfs_cache_pp(lfs_pp);
		typedef typename U_pp::template LocalView<LFSCache_PP> VectorView_PP;
		VectorView_PP u_pp_view( (*u_pp) );

		typedef typename LFS_PP::template Child<0>::Type::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_volume = geo.volume();

	        RF Pw=0.;
	        evaluation_Pw->evalFunction(cell,cell_center_local,&Pw);
	        RF Sg=0.;
	        evaluation_Sg->evalFunction(cell,cell_center_local,&Sg);
	        RF XCH4=0.;
	        evaluation_XCH4->evalFunction(cell,cell_center_local,&XCH4);
	        RF YH2O=0.;
	        evaluation_YH2O->evalFunction(cell,cell_center_local,&YH2O);
	        RF Pw0=0.;
	        evaluation_pw_old->evalFunction(cell,cell_center_local,&Pw0);
	        RF Sg0=0.;
	        evaluation_sg_old->evalFunction(cell,cell_center_local,&Sg0);
	        RF por0=0.;
	        evaluation_por_old->evalFunction(cell,cell_center_local,&por0);


	        lfs_pp.bind(*self);
			lfs_cache_pp.update();
			u_pp_view.bind(lfs_cache_pp);
	        std::vector<double> ul_pp( lfs_pp.size() );
	        for(int i = 0. ; i < lfs_pp.size() ; i++){
	        	ul_pp[i] = 0.;
	        }

	        RF Pg0 = Pw0 + param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg0, por0 )/Xc_P;
	        RF Pg_tmp = Pw + param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg, por0 )/Xc_P;
//	        double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg0,0.)*Pg0 ) /(*dt);
//	        double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg0,0.)*Pg0 );
	        double dppg = ( std::max(Sg,0.)*Pg_tmp+(1.-Sg)*Pw ) - ( std::max(Sg0,0.)*Pg0+(1.-Sg0)*Pw0 );
			RF porosity = param.hydraulicProperty.ScaledPorosity(cell, cell_center_local, por0,dppg*Xc_P );
			RF K = param.soil.SedimentPermeability( cell,cell_center_local )
				 * param.hydraulicProperty.PermeabilityScalingFactor( cell,cell_center_local, porosity );

			RF S = param.parameter.ReferenceSalinity();
			RF Xc = param.parameter.ReferenceSaltConcentration();
			RF T = param.parameter.ReferenceTemperature()/Xc_T; /*ndim*/
			RF Sw = 1.-Sg;
			RF Pc = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw, porosity )/Xc_P;
			RF Pg = Pw + Pc;
			
			RF zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );
			
			RF XH2O = param.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,S);
			RF YCH4 = param.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,S,zCH4);
			
			RF rhog = param.gas.Density( T*Xc_T, Pg*Xc_P, zCH4 );
			RF rhow = param.water.Density( T*Xc_T, Pw*Xc_P, S );
			
			RF krg = param.hydraulicProperty.krg( cell, cell_center_local, Sw );
			RF krw = param.hydraulicProperty.krw( cell, cell_center_local, Sw );
			
			RF mug = param.gas.DynamicViscosity( T*Xc_T, Pg*Xc_P );
			RF muw = param.water.DynamicViscosity( T*Xc_T, Pw*Xc_P, S );
			
			RF tau = param.soil.Tortuosity( porosity );
			RF DH2O_g	= tau * porosity * Sg * param.mixture.DiffCoeffH2OInGas( T*Xc_T,Pg*Xc_P );
			RF DCH4_w	= tau * porosity * Sw * param.mixture.DiffCoeffCH4InLiquid( T*Xc_T,Pw*Xc_P );

			RF Pwsat = param.water.SaturatedVaporPressure( T*Xc_T,S );
			RF Hch4 = param.gas.SolubilityCoefficient(T*Xc_T,S);

			ul_pp[lfs_pp_Pg.localIndex(0)]	  = Pg*Xc_P ;
	        ul_pp[lfs_pp_Pw.localIndex(0)] 	  = Pw*Xc_P ;
	        ul_pp[lfs_pp_Pc.localIndex(0)]	  = Pc*Xc_P ;
	        ul_pp[lfs_pp_Sg.localIndex(0)]	  = Sg ;
	        ul_pp[lfs_pp_Sw.localIndex(0)] 	  = Sw ;
	        ul_pp[lfs_pp_XCH4.localIndex(0)]  = XCH4 ;
	        ul_pp[lfs_pp_XH2O.localIndex(0)]  = XH2O ;
	        ul_pp[lfs_pp_YCH4.localIndex(0)]  = YCH4 ;
	        ul_pp[lfs_pp_YH2O.localIndex(0)]  = YH2O ;
	        ul_pp[lfs_pp_rhow.localIndex(0)]  = rhow*Xc_rho ;
	        ul_pp[lfs_pp_rhog.localIndex(0)]  = rhog*Xc_rho ;
	        ul_pp[lfs_pp_K.localIndex(0)] 	  = K*Xc_K ;
	        ul_pp[lfs_pp_krw.localIndex(0)]   = krw ;
	        ul_pp[lfs_pp_krg.localIndex(0)]   = krg ;
	        ul_pp[lfs_pp_muw.localIndex(0)]   = muw*Xc_mu ;
	        ul_pp[lfs_pp_mug.localIndex(0)]   = mug*Xc_mu ;
	        ul_pp[lfs_pp_zCH4.localIndex(0)]  = zCH4 ;
	        ul_pp[lfs_pp_por.localIndex(0)]   = porosity ;
	        ul_pp[lfs_pp_DH2O.localIndex(0)]  = DH2O_g*Xc_D ;
	        ul_pp[lfs_pp_DCH4.localIndex(0)]  = DCH4_w*Xc_D ;
	        ul_pp[lfs_pp_Pwsat.localIndex(0)] = Pwsat*Xc_P ;
	        ul_pp[lfs_pp_HCH4.localIndex(0)]  = Hch4*Xc_P ;
	        ul_pp[lfs_pp_tau.localIndex(0)]   = tau ;

			u_pp_view.write( ul_pp );
			u_pp_view.commit();
			u_pp_view.unbind();

		}//END:iterate over each volume

	}


	void update_porosity(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS0 > LFS0;
		LFS0 lfs0(gfs0);

		typedef Dune::PDELab::LFSIndexCache<LFS0> LFSCache0;
		LFSCache0 lfs_cache0(lfs0);
		typedef typename U0::template LocalView<LFSCache0> VectorView0;
		VectorView0 u0_por_view((*u0_por_up));

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_volume = geo.volume();

	        double Pw=0.;
	        evaluation_Pw->evalFunction(cell,cell_center_local,&Pw);
	        double Pw0=0.;
	        evaluation_pw_old->evalFunction(cell,cell_center_local,&Pw0);
	        double Sg=0.;
	        evaluation_Sg->evalFunction(cell,cell_center_local,&Sg);
	        double Sg0=0.;
	        evaluation_sg_old->evalFunction(cell,cell_center_local,&Sg0);
	        double por0=0.;
	        evaluation_por_old->evalFunction(cell,cell_center_local,&por0);

	        lfs0.bind(*self);
			lfs_cache0.update();
			u0_por_view.bind(lfs_cache0);
	        std::vector<double> ul0_por_up( lfs0.size(),0. );

	        double Pg0 = Pw0 + param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg0, por0 )/Xc_P;
	        double Pg_tmp = Pw + param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg, por0 )/Xc_P;
//	        double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg0,0.)*Pg0 ) /(*dt);
//	        double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg0,0.)*Pg0 );
	        double dppg = ( std::max(Sg,0.)*Pg_tmp+(1.-Sg)*Pw ) - ( std::max(Sg0,0.)*Pg0+(1.-Sg0)*Pw0 );
	        ul0_por_up[lfs0.localIndex(0)] = param.hydraulicProperty.ScaledPorosity(cell, cell_center_local, por0, dppg*Xc_P );

			u0_por_view.write( ul0_por_up );
			u0_por_view.commit();
			u0_por_view.unbind();

		}//END:iterate over each volume
	}


	void update_pressure(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS0 > LFS0;
		LFS0 lfs0(gfs0);

		typedef Dune::PDELab::LFSIndexCache<LFS0> LFSCache0;
		LFSCache0 lfs_cache0(lfs0);
		typedef typename U0::template LocalView<LFSCache0> VectorView0;
		VectorView0 u0_pw_view((*u0_pw_up));

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_volume = geo.volume();

	        double Pw=0.;
	        evaluation_Pw->evalFunction(cell,cell_center_local,&Pw);

	        lfs0.bind(*self);
			lfs_cache0.update();
			u0_pw_view.bind(lfs_cache0);
	        std::vector<double> ul0_pw_up( lfs0.size(),0. );
	        ul0_pw_up[lfs0.localIndex(0)] = Pw;

			u0_pw_view.write( ul0_pw_up );
			u0_pw_view.commit();
			u0_pw_view.unbind();

		}//END:iterate over each volume
	}


	void update_saturation(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS0 > LFS0;
		LFS0 lfs0(gfs0);

		typedef Dune::PDELab::LFSIndexCache<LFS0> LFSCache0;
		LFSCache0 lfs_cache0(lfs0);
		typedef typename U0::template LocalView<LFSCache0> VectorView0;
		VectorView0 u0_sg_view((*u0_sg_up));

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_volume = geo.volume();

	        double Sg=0.;
	        evaluation_Sg->evalFunction(cell,cell_center_local,&Sg);

	        lfs0.bind(*self);
			lfs_cache0.update();
			u0_sg_view.bind(lfs_cache0);
	        std::vector<double> ul0_sg_up( lfs0.size(),0. );
	        ul0_sg_up[lfs0.localIndex(0)] = Sg;

			u0_sg_view.write( ul0_sg_up );
			u0_sg_view.commit();
			u0_sg_view.unbind();

		}//END:iterate over each volume
	}


};
