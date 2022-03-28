template< typename GV, typename PTree>
void driver( const GV& gv, 						// GridView
			 const PTree& ptree, 		// Input Parameters
			 Dune::MPIHelper& helper){


	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
		typedef typename GV::Grid::ctype Coord;
		const int dim = GV::dimension;
		typedef double RF;

	//	MATERIAL PROPERTIES, NUMERICAL AND TEST PARAMETERS, CONSTITUTIVE RELATIONSHIPS
		typedef Properties<GV,PTree> Properties;
		Properties property(gv,ptree);

		std::string pathName = "/sfs/fs1/work-geomar5/smomw325/dune_outputs/dune_2_8/Pockmark2P2C/case1_Test0_1D/";
		auto pathExt = ptree.get("output.path_name",(std::string)"test0");
		pathName += pathExt;
		pathName += "/";
		auto fileName = ptree.get("output.file_name",(std::string)"test");

		/*Non-dimensionalize time prams*/
		RF dt_ACTUAL  = ptree.get("time.dt_initial",(double)0.001);
		dt_ACTUAL 	 *= 1./property.characteristicValue.t_c;
		RF dt_INIT = ptree.get("before_storm.dt_initial",(double)0.001);
		dt_INIT   *= 1./property.characteristicValue.t_c;
		RF t_END  = ptree.get("time.time_end",(double)86400.);
		t_END 	 *= 1./property.characteristicValue.t_c;
		RF t_INIT  = ptree.get("water_column.wave_start_time",(double)0.);
		t_INIT 	 *= 1./property.characteristicValue.t_c;
		// output time interval
		RF t_OP_ACTUAL = ptree.get("output.time_interval",(double)1.);
		t_OP_ACTUAL   *= 1./property.characteristicValue.t_c;
		RF t_OP_INIT   = ptree.get("before_storm.output_time_interval",(double)1.);
		t_OP_INIT 	  *= 1./property.characteristicValue.t_c;
		//adaptive time control
		RF dt_min = ptree.get("adaptive_time_control.dt_min",(double)0.001);
		dt_min   *= 1./property.characteristicValue.t_c;
		RF dt_max = ptree.get("adaptive_time_control.dt_max",(double)1.);
		dt_max   *= 1./property.characteristicValue.t_c;
		int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)4);
		int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)2);

		RF dt=dt_INIT;
		RF t_OP=t_OP_INIT;
		RF time = 0.0;
		RF dtstart = dt;
		RF time_op = time;
		RF clock_time_elapsed = 0.;


	//	COMPOSITE GFS FOR PRIMARY VARIABLES
		typedef typename GV::Grid::ctype Coord;
		typedef Dune::PDELab::P0ParallelConstraints CON0;
		using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;					// default block size: 1
		typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord,RF,0,dim,Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM0;
		FEM0 fem0;
		typedef Dune::PDELab::GridFunctionSpace<GV, FEM0, CON0, VBE0> GFS0;
		GFS0 gfs0(gv, fem0);
		// gfs for composite system: Pw , Sg , XCH4 , YH2O
		using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFS0,
													  Indices::numOfPVs,
													  VBE,
													  Dune::PDELab::EntityBlockedOrderingTag > GFS;
		GFS gfs(gfs0);
		typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
		CC cc;
		cc.clear();

	//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Pw	>> > SUBGFS_Pw;
		SUBGFS_Pw	subgfs_Pw(gfs);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sg	>> > SUBGFS_Sg;
		SUBGFS_Sg	subgfs_Sg(gfs);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_XCH4>> > SUBGFS_XCH4;
		SUBGFS_XCH4	subgfs_XCH4(gfs);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_YH2O>> > SUBGFS_YH2O;
		SUBGFS_YH2O	subgfs_YH2O(gfs);

	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
		using U = Dune::PDELab::Backend::Vector<GFS,double>;
		U u_old(gfs,0.);
		U u_new(gfs,0.);

	//  EVALUATION FUNCTIONS FOR PRIMARY VARIABLES
		Dune::PDELab::Evaluation<SUBGFS_Pw	,U> evaluation_Pw(	 subgfs_Pw	,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_Sg	,U> evaluation_Sg(	 subgfs_Sg	,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_XCH4,U> evaluation_XCH4( subgfs_XCH4,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_YH2O,U> evaluation_YH2O( subgfs_YH2O,&u_new );

	//	INITIAL CONDITIONS
		//	Make function for initial values
		typedef Pw_Initial<GV,Properties> ICV_Pw;
		ICV_Pw Pw_initial(gv,property);
		typedef Sg_Initial<GV,Properties> ICV_Sg;
		ICV_Sg Sg_initial(gv,property);
		typedef XCH4_Initial<GV,Properties> ICV_XCH4;
		ICV_XCH4 XCH4_initial(gv,property);
		typedef YH2O_Initial<GV,Properties> ICV_YH2O;
		ICV_YH2O YH2O_initial(gv,property);
		typedef Dune::PDELab::CompositeGridFunction< ICV_Pw,
													 ICV_Sg,
													 ICV_XCH4,
													 ICV_YH2O > InitialValues;
		InitialValues icv( Pw_initial, Sg_initial, XCH4_initial, YH2O_initial );

		// 	Initialize the solution at t=0 (uold) with the given initial values
		Dune::PDELab::interpolate( icv, gfs, u_old );
		u_new = u_old;

	//	BOUNDARY CONDITIONS
		typedef ProblemBoundaryConditions<GV,Properties> BoundaryConditions ;
		BoundaryConditions bc( gv,property ) ;

	//  OLD-POROSFTY AND OLD-PEFF VARIABLES FOR EVALUATION OF NEW-POROSITY
		using U0 = Dune::PDELab::Backend::Vector<GFS0,double>;
		U0 u0_por(gfs0,0.);
		U0 u0_pw(gfs0,0.);
		U0 u0_sg(gfs0,0.);
		Dune::PDELab::Evaluation<GFS0,U0> evaluation_por_old( gfs0	,&u0_por );
		Dune::PDELab::Evaluation<GFS0,U0> evaluation_pw_old( gfs0	,&u0_pw  );
		Dune::PDELab::Evaluation<GFS0,U0> evaluation_sg_old( gfs0	,&u0_sg  );
		typedef Por_Initial<GV,Properties> POR_INIT;
		POR_INIT initial_porosity(gv,property);
		Dune::PDELab::interpolate( initial_porosity,  gfs0, u0_por );
		typedef Pw_Initial<GV,Properties> PW_INIT;
		PW_INIT initial_pressure(gv,property);
		Dune::PDELab::interpolate( initial_pressure,  gfs0, u0_pw  );
		typedef Sg_Initial<GV,Properties> SG_INIT;
		SG_INIT initial_saturation(gv,property);
		Dune::PDELab::interpolate( initial_saturation,gfs0, u0_sg  );

	//	MAKE INSTATIONARY GRID OPERATOR SPACE

		typedef LocalOperator< GV,
							   Properties,
							   BoundaryConditions,
							   Dune::PDELab::Evaluation<GFS0,U0>,
							   Dune::PDELab::Evaluation<GFS0,U0>,
							   Dune::PDELab::Evaluation<GFS0,U0> > LOP;	// spatial part
		LOP lop( gv,
				property,
				bc,
				&evaluation_por_old,
				&evaluation_pw_old,
				&evaluation_sg_old,
				&time, &dt );

		typedef TimeOperator< GV, Properties,
							  Dune::PDELab::Evaluation<GFS0,U0>,
							  Dune::PDELab::Evaluation<GFS0,U0>,
							  Dune::PDELab::Evaluation<GFS0,U0> > TLOP; // temporal part
		TLOP tlop( gv, property, &evaluation_por_old, &evaluation_pw_old, &evaluation_sg_old, &time, &dt );

		typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
		MBE mbe(6);

		typedef Dune::PDELab::GridOperator< GFS, GFS, LOP, MBE, RF, RF, RF, CC, CC> GOLOP;
		GOLOP golop(gfs, cc, gfs, cc, lop, mbe);

		// How well did we estimate the number of entries per matrix row?
		// => print Jacobian pattern statistics
		typename GOLOP::Traits::Jacobian jac(golop);
		if(helper.rank()==0){
			std::cout << jac.patternStatistics() << std::endl;
		}

		typedef Dune::PDELab::GridOperator< GFS, GFS, TLOP, MBE, RF, RF, RF, CC, CC > GOTLOP;
		GOTLOP gotlop(gfs, cc, gfs, cc, tlop, mbe);

		typedef Dune::PDELab::OneStepGridOperator< GOLOP, GOTLOP > IGO;
		IGO igo( golop, gotlop );

	// SELECT A LINEAR SOLVER BACKEND
		typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
		LS ls(0);
//		typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
//		LS ls(gfs,100,1,true,true);
//		typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS; //works
//		LS ls(gfs,500,1,false,true);

	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
		typedef Dune::PDELab::Newton< IGO, LS, U > PDESOLVER;
		PDESOLVER pdesolver( igo, ls );
	// 	select control parameters for non-linear PDE-solver
		pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::noLineSearch);
		pdesolver.setReassembleThreshold(0.0);
		pdesolver.setVerbosityLevel(2);
		pdesolver.setReduction(1e-6);
		pdesolver.setMinLinearReduction(1e-6);
		pdesolver.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
		pdesolver.setForceIteration(true);
		pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-6));


	//	SELECT TIME-STEPPER
		Dune::PDELab::ImplicitEulerParameter<RF> method;

		Dune::PDELab::OneStepMethod< RF, IGO, PDESOLVER, U, U > osm( method, igo, pdesolver );
		osm.setVerbosityLevel(2);


	//	POST-PROCESS
		typedef Dune::PDELab::GridFunctionSpace<GV, FEM0, CON0, VBE0> GFS_PP0;
		GFS_PP0 gfs_pp0(gv, fem0);
		typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none > VBE_PP;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFS_PP0,
													  Indices::numOfSVs,
													  VBE_PP,
													  Dune::PDELab::LexicographicOrderingTag > GFS_PP;
		GFS_PP gfs_pp(gfs_pp0);
		typedef typename GFS_PP::template ConstraintsContainer<RF>::Type CC_PP;
		CC_PP cc_pp;
		cc_pp.clear();

		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pg	>> > SUBPP_Pg;
		SUBPP_Pg 	subpp_Pg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pw	>> > SUBPP_Pw;
		SUBPP_Pw 	subpp_Pw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pc	>> > SUBPP_Pc;
		SUBPP_Pc 	subpp_Pc(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sg	>> > SUBPP_Sg;
		SUBPP_Sg 	subpp_Sg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sw	>> > SUBPP_Sw;
		SUBPP_Sw 	subpp_Sw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XCH4	>> > SUBPP_XCH4;
		SUBPP_XCH4 	subpp_XCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XH2O >> > SUBPP_XH2O;
		SUBPP_XH2O 	subpp_XH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_YCH4 >> > SUBPP_YCH4;
		SUBPP_YCH4 	subpp_YCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_YH2O >> > SUBPP_YH2O;
		SUBPP_YH2O 	subpp_YH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_rhow	>> > SUBPP_rhow;
		SUBPP_rhow 	subpp_rhow(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_rhog	>> > SUBPP_rhog;
		SUBPP_rhog 	subpp_rhog(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_K	>> > SUBPP_K;
		SUBPP_K 	subpp_K(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_krw	>> > SUBPP_krw;
		SUBPP_krw 	subpp_krw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_krg	>> > SUBPP_krg;
		SUBPP_krg 	subpp_krg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_muw	>> > SUBPP_muw;
		SUBPP_muw 	subpp_muw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_mug	>> > SUBPP_mug;
		SUBPP_mug 	subpp_mug(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_zCH4 >> > SUBPP_zCH4;
		SUBPP_zCH4 	subpp_zCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_por	 >> > SUBPP_por;
		SUBPP_por 	subpp_por(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_DH2O	 >> > SUBPP_DH2O;
		SUBPP_DH2O 	subpp_DH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_DCH4	 >> > SUBPP_DCH4;
		SUBPP_DCH4 	subpp_DCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pwsat >> > SUBPP_Pwsat;
		SUBPP_Pwsat 	subpp_Pwsat(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_HCH4	 >> > SUBPP_HCH4;
		SUBPP_HCH4 	subpp_HCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_tau	 >> > SUBPP_tau;
		SUBPP_tau 	subpp_tau(gfs_pp);

		using U_PP = Dune::PDELab::Backend::Vector<GFS_PP,double>;
		U_PP u_pp(gfs_pp, 0.0);

		PostProcess< GV, Properties,
					 Dune::PDELab::Evaluation<SUBGFS_Pw  ,U>,
					 Dune::PDELab::Evaluation<SUBGFS_Sg	 ,U>,
					 Dune::PDELab::Evaluation<SUBGFS_XCH4,U>,
					 Dune::PDELab::Evaluation<SUBGFS_YH2O,U>,
					 Dune::PDELab::Evaluation<GFS0,U0>,
					 Dune::PDELab::Evaluation<GFS0,U0>,
					 Dune::PDELab::Evaluation<GFS0,U0>,
					 GFS_PP, U_PP,
					 GFS0, U0 > postprocess(gv, property,
											&evaluation_Pw,
											&evaluation_Sg,
											&evaluation_XCH4,
											&evaluation_YH2O,
											&evaluation_por_old,
											&evaluation_pw_old,
											&evaluation_sg_old,
											gfs_pp, &u_pp,
											gfs0, &u0_por, &u0_pw, &u0_sg,
											&time, &dt);
		postprocess.evaluate();

		// QOI
		int csvcount=0;
		std::string qoi_file = pathName;
		qoi_file +=fileName;
		QoI<GV,Properties,GFS_PP,U_PP> qoi( qoi_file, gv, property,
											gfs_pp, &u_pp,
											&time, &dt, &csvcount);
		qoi.generate_csv_output();

		// prepare VTK writer and write first file

		// Make a grid function out of it
		//	GRAPHICS FOR INITIAL GUESS
		// secondary variables
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pg, U_PP > DGF_Pg;
		DGF_Pg dgf_Pg( subpp_Pg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pw, U_PP > DGF_Pw;
		DGF_Pw dgf_Pw( subpp_Pw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pc, U_PP > DGF_Pc;
		DGF_Pc dgf_Pc( subpp_Pc, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sg, U_PP > DGF_Sg;
		DGF_Sg dgf_Sg( subpp_Sg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sw, U_PP > DGF_Sw;
		DGF_Sw dgf_Sw( subpp_Sw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XCH4, U_PP > DGF_XCH4;
		DGF_XCH4 dgf_XCH4( subpp_XCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XH2O, U_PP > DGF_XH2O;
		DGF_XH2O dgf_XH2O( subpp_XH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_YCH4, U_PP > DGF_YCH4;
		DGF_YCH4 dgf_YCH4( subpp_YCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_YH2O, U_PP > DGF_YH2O;
		DGF_YH2O dgf_YH2O( subpp_YH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_rhow, U_PP > DGF_rhow;
		DGF_rhow dgf_rhow( subpp_rhow, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_rhog, U_PP > DGF_rhog;
		DGF_rhog dgf_rhog( subpp_rhog, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_K, U_PP > DGF_K;
		DGF_K dgf_K( subpp_K, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_krw, U_PP > DGF_krw;
		DGF_krw dgf_krw( subpp_krw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_krg, U_PP > DGF_krg;
		DGF_krg dgf_krg( subpp_krg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_muw, U_PP > DGF_muw;
		DGF_muw dgf_muw( subpp_muw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_mug, U_PP > DGF_mug;
		DGF_mug dgf_mug( subpp_mug, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_zCH4, U_PP > DGF_zCH4;
		DGF_zCH4 dgf_zCH4( subpp_zCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_por, U_PP > DGF_por;
		DGF_por dgf_por( subpp_por, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_DH2O, U_PP > DGF_DH2O;
		DGF_DH2O dgf_DH2O( subpp_DH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_DCH4, U_PP > DGF_DCH4;
		DGF_DCH4 dgf_DCH4( subpp_DCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pwsat, U_PP > DGF_Pwsat;
		DGF_Pwsat dgf_Pwsat( subpp_Pwsat, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_HCH4, U_PP > DGF_HCH4;
		DGF_HCH4 dgf_HCH4( subpp_HCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_tau, U_PP > DGF_tau;
		DGF_tau dgf_tau( subpp_tau, u_pp );


		int subsampling = 1;
		using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
		VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
		using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
		VTKSEQUENCEWRITER vtkSequenceWriter( std::make_shared<VTKWRITER>(vtkwriter),fileName,pathName,"");

		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pg    > >( dgf_Pg    , "Pg"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pw    > >( dgf_Pw    , "Pw"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pc    > >( dgf_Pc    , "Pc"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Sg    > >( dgf_Sg    , "Sg"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Sw    > >( dgf_Sw    , "Sw"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_XCH4  > >( dgf_XCH4  , "XCH4" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_XH2O  > >( dgf_XH2O  , "XH2O" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_YCH4  > >( dgf_YCH4  , "YCH4" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_YH2O  > >( dgf_YH2O  , "YH2O" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_rhow  > >( dgf_rhow  , "rhow" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_rhog  > >( dgf_rhog  , "rhog" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_K     > >( dgf_K     , "K"    ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_krw   > >( dgf_krw   , "krw"  ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_krg   > >( dgf_krg   , "krg"  ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_muw   > >( dgf_muw   , "muw"  ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_mug   > >( dgf_mug   , "mug"  ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_zCH4  > >( dgf_zCH4  , "z"    ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_por   > >( dgf_por   , "porosity"));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_DH2O  > >( dgf_DH2O  , "DH2O" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_DCH4  > >( dgf_DCH4  , "DCH4" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pwsat > >( dgf_Pwsat , "Pwsat"));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_HCH4  > >( dgf_HCH4  , "HCH4" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_tau   > >( dgf_tau   , "tau"  ));

		vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

		/***********************************************/
		//	BEGIN TIME LOOP
		/***********************************************/
		int opcount = 1;
		bool flip_opcount_flag = false;
		double timecount = time;
		double dtLast = dtstart;
		int dtFlag = 0;
		bool exceptionCaught = false;
		int newton_iterations = 0;

		while( time < t_END - 1e-8/CharacteristicValues::t_c ){

			if( exceptionCaught==false ){
				dt = std::max(dt,dt_min);
			}

			if(helper.rank()==0){
			std::cout<< "_____________________________________________________" <<std::endl;
			std::cout<< " current opcount = " << opcount - 1 << std::endl;
			}

			clock_t start = clock();
			try{
				if(helper.rank()==0){
				std::cout<<"****************************" << std::endl;
				std::cout<<"  CALLING osm.apply() !"	  << std::endl;
				std::cout<<"****************************" << std::endl;
				}

				osm.apply( time, dt, u_old, u_new );

				newton_iterations = osm.getPDESolver().result().iterations;

				exceptionCaught = false;

			}catch ( Dune::Exception &e ) {
				exceptionCaught = true;
				if( dt*property.characteristicValue.t_c > 1e-8 ){

					if(helper.rank()==0){
					std::cout << "Catched Error, Dune reported error: " << e << std::endl;
					}

					u_new = u_old;

					newton_iterations = 0;

					dt *= 0.5;
					continue;
				}
				else
				{
					if(helper.rank()==0){
						std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
					}
					exit(0);
				}
			}
			clock_t end = clock();
			double clock_time_this_step = (double) (end-start) / CLOCKS_PER_SEC;
			clock_time_elapsed += clock_time_this_step;

			if(helper.rank()==0){
			std::cout<<"DONE"<<std::endl;
			std::cout<<"_____________________________________________________"<<std::endl;
			}


			/*********************************************************************************************
			 * OUTPUT
			 *********************************************************************************************/
			/* At each time step: **Statistics**
			 * t_new,
			 * dt,
			 * fixed point iterations,
			 * newton iterations per fixed point iteration,
			 * total newton terations
			 */

			if(helper.rank()==0){
			std::string statistics_file = pathName;
			statistics_file +=fileName;
			statistics_file +="_statistics";
			statistics_file += ".txt";
			property.ReportStatistics( 	statistics_file,
										time*CharacteristicValues::t_c,
										dt*CharacteristicValues::t_c,
										newton_iterations,
										clock_time_elapsed );
			}

			/* At t_OP
			 *
			 */
			if( 	( (time<t_INIT)  and ( time+dt > t_OP*opcount - dt_min ) and ( time+dt <= t_OP*opcount ) )
				or  ( (time>=t_INIT) and ( time+dt > t_INIT + t_OP*opcount - dt_min ) and ( time+dt <= t_INIT + t_OP*opcount ) ) )
			{
				// POST PROCESS FOR NEW OUTPUT
				postprocess.evaluate();

				// WRITE CSV OUTPUT
				if(time>t_INIT){
					csvcount +=1;
					qoi.generate_csv_output();
				}

//				// WRITE OUTPUT
				vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

				if(helper.rank()==0){
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< " OUTPUT WRITTEN " << opcount << std::endl;
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< std::flush;
				}

				timecount = time;
				opcount += 1;
			}

	//		PREPARE FOR NEXT TIME INTEGRATION
	//		1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
			u_old = u_new;
			postprocess.update_porosity();
			postprocess.update_pressure();
			postprocess.update_saturation();

	//		2. ADVANCE TIME:
			time += dt;

			//

					if( !flip_opcount_flag and time>=t_INIT ){
						opcount = 1;
						flip_opcount_flag = true;
					}

					if( time>=t_INIT ){
						t_OP = t_OP_ACTUAL;
						dtstart = dt_ACTUAL;
					}

					if(helper.rank()==0){
					std::cout<<" "<< std::endl;
					std::cout<< " time = " << time*property.characteristicValue.t_c ;
					std::cout<< std::flush;
					}

					if( ptree.get("adaptive_time_control.flag", (bool)false) ){
						if(newton_iterations>maxAllowableIterations){
							dt=std::max(dt*0.9,dt_min);
						}
						else if(newton_iterations<=minAllowableIterations){
							dt=std::min(dt*1.1,dt_max);
						}
					}
					else{
						dt = dtstart;
					}

					if(time>=t_INIT){
						if(helper.rank()==0){
						std::cout << " , time+dt = " << (time + dt)*property.characteristicValue.t_c
								  << " , opNext = "  << (t_INIT + t_OP * opcount) * property.characteristicValue.t_c ;
						std::cout<< std::flush;
						}
						if( time + dt  >= t_INIT + t_OP * opcount){
							dtLast = dt;
							dt = t_INIT + t_OP * opcount - time ;

							if(helper.rank()==0){
							std::cout<< " . Because timeNext > opNext , dt set to : " << dt*property.characteristicValue.t_c << std::endl;
							std::cout<< std::flush;
							}

							dtFlag = 0;
						}
						dtFlag += 1;

					}else{
						if(helper.rank()==0){
						std::cout << " , time+dt = " << (time + dt)*property.characteristicValue.t_c
								  << " , opNext = "  << t_OP * opcount * property.characteristicValue.t_c ;
						std::cout<< std::flush;
						}
						if( time + dt  >= t_OP * opcount){
							dtLast = dt;
							dt = t_OP * opcount - time ;

							if(helper.rank()==0){
							std::cout<< " . Because timeNext > opNext , dt set to : " << dt*property.characteristicValue.t_c << std::endl;
							std::cout<< std::flush;
							}

							dtFlag = 0;
						}
						dtFlag += 1;
					}

					if(helper.rank()==0){
					std::cout<< " , dt  : " << dt*property.characteristicValue.t_c << std::endl;
					std::cout<<" "<< std::endl;
					std::cout << " READY FOR NEXT ITERATION. " << std::endl;
					std::cout<< std::flush;
					}
				}

		}
