template <	class GV, class Params, class BC, class EVALPor, class EVALPw, class EVALSg >
class LocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperator<GV,Params,BC,EVALPor,EVALPw,EVALSg> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperator<GV,Params,BC,EVALPor,EVALPw,EVALSg> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperator<GV,Params,BC,EVALPor,EVALPw,EVALSg> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperator<GV,Params,BC,EVALPor,EVALPw,EVALSg> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperator<GV,Params,BC,EVALPor,EVALPw,EVALSg> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperator<GV,Params,BC,EVALPor,EVALPw,EVALSg> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
    const GV& gv;
	const Params&	  param;
	const BC&	 	  bc;
	EVALPor			 *eval_porosity_old;
	EVALPw			 *eval_pw_old;
	EVALSg			 *eval_sg_old;
	double 			 *time;
	double 			 *dt;
	constexpr static double eps 	= 1.0e-6;
	constexpr static double eps_ap	= 0.;//1.e-9;
	constexpr static double pi 		= atan(1.)*4 ;
	double Xc_conv_m ;
	double Xc_source_m ;
	double Xc_diff_m ;
	double Xc_grav ;
	double Xc_vs ;
	double Xc_K ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_D;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;


public:
	  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor stores parameters
	  LocalOperator( const GV& 		 gv_,
					 const Params&	 param_,
					 const BC& 	 	 bc_	,
					 EVALPor		*eval_porosity_old_,
					 EVALPw			*eval_pw_old_,
					 EVALSg			*eval_sg_old_,
					 double			*time_	,
					 double 		*dt_	)
	  :  gv(gv_),
		 param( param_ ),
		 bc( bc_ ),
		 eval_porosity_old(eval_porosity_old_),
		 eval_pw_old(eval_pw_old_),
		 eval_sg_old(eval_sg_old_),
		 time( time_ ),
		 dt( dt_ )
	  {
		  Xc_conv_m 	= param.characteristicValue.X_convective_mass;
		  Xc_source_m 	= param.characteristicValue.X_source_mass;
		  Xc_diff_m 	= param.characteristicValue.X_diffusive_mass;
		  Xc_grav 		= param.characteristicValue.X_gravity;
		  Xc_vs			= param.characteristicValue.X_solidvelocity;

		  Xc_K 		= param.characteristicValue.permeability_c;
		  Xc_mu 	= param.characteristicValue.viscosity_c;
		  Xc_rho 	= param.characteristicValue.density_c;
		  Xc_D		= param.characteristicValue.dispersivity_c;
		  Xc_P 		= param.characteristicValue.P_c;
		  Xc_T 		= param.characteristicValue.T_c;
		  Xc_t 		= param.characteristicValue.t_c;
		  Xc_x 		= param.characteristicValue.x_c;
	  }

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
	        typedef typename LFSU::template Child< Indices::PVId_Pw	 >::Type LFS_Pw ;
	        const LFS_Pw&	lfs_Pw	= lfsu.template child< Indices::PVId_Pw  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_Sg	 >::Type LFS_Sg;
	        const LFS_Sg& 	lfs_Sg	= lfsu.template child< Indices::PVId_Sg  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_XCH4>::Type LFS_XCH4;
	        const LFS_XCH4&	lfs_XCH4= lfsu.template child< Indices::PVId_XCH4>() ;
	        typedef typename LFSU::template Child< Indices::PVId_YH2O>::Type LFS_YH2O;
	        const LFS_YH2O&	lfs_YH2O= lfsu.template child< Indices::PVId_YH2O>() ;

			// Reference to cell
	        const auto& cell = eg.entity();
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);

	        // get geometry
	        auto geo = eg.geometry();

			// dimension
			const auto dim = geo.mydimension;

	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_center_global = geo.center();
	        auto cell_volume = geo.volume();

			// compute PVs at local center
			double Pw   = x(lfs_Pw  ,0);
			double Sg   = x(lfs_Sg  ,0);
			double XCH4 = x(lfs_XCH4,0);
			double YH2O = x(lfs_YH2O,0);

			// reference salinity and temperature
			auto S = param.parameter.ReferenceSalinity();
			auto Xc = param.parameter.ReferenceSaltConcentration();
			auto T /*ndim*/ = param.parameter.ReferenceTemperature()/Xc_T;

	        // compute the old porosity and old partial Pg (ppg) at self and neighbour elements
	        double porosity_old=0.;
	        eval_porosity_old->evalFunction(cell,cell_center_local,&porosity_old);
	        double Pw_old=0.;
	        eval_pw_old->evalFunction(cell,cell_center_local,&Pw_old);
	        double Sg_old=0.;
	        eval_sg_old->evalFunction(cell,cell_center_local,&Sg_old);
	        double Sw_old=1.-Sg_old;
	        double Pc_old=param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw_old, porosity_old );
	        double Pg_old=Pw_old+Pc_old;
			double Pg_tmp=Pw+param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg, porosity_old );
//			double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg_old,0.)*Pg_old ) /(*dt);
//			dppg *= (Xc_P/Xc_t);
//			double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg_old,0.)*Pg_old );
			double dppg = ( std::max(Sg,0.)*Pg_tmp+(1.-Sg)*Pw ) - ( std::max(Sg_old,0.)*Pg_old+(1.-Sg_old)*Pw_old );
			dppg *= Xc_P;

			// secondary variables
			auto porosity = param.hydraulicProperty.ScaledPorosity(cell, cell_center_local, porosity_old, dppg );
			auto permeability = param.soil.SedimentPermeability( cell,cell_center_local )
							  * param.hydraulicProperty.PermeabilityScalingFactor( cell,cell_center_local, porosity );
			auto Sw = 1.-Sg;
			auto Pc = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw, porosity );
			auto Pg = Pw + Pc;
			auto zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );
			
			// compute source terms
			auto q_CH4  = param.source.GasGenerationRate(T*Xc_T,Pg*Xc_P);
			auto q_H2O  = param.source.WaterGenerationRate(T*Xc_T,Pw*Xc_P);

			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp -= Xc_source_m * q_CH4 ;
//			std::cout<< "alpha_vol 0: " << tmp << std::endl;
			r.accumulate(lfs_Sg, 0, +tmp*cell_volume);

			// H2O-component-wise mass-balance
			tmp = 0.;
			tmp -= Xc_source_m * q_H2O ;
//			std::cout<< "alpha_vol 1: " << tmp << std::endl;
			r.accumulate(lfs_Pw, 0, +tmp*cell_volume);

			// NCP -> water phase
			tmp = 0.;
			auto XH2O_alg = param.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,S);
			if( ( Sw - ( 1. - XCH4 - XH2O_alg - Xc ) ) > eps_ap ){//active set.
				tmp += 1. - XCH4 - XH2O_alg - Xc;//Active => phase is present => summation condition holds
//				std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
			}else{
				tmp += Sw; // inactive set. Inactive => phase is absent => Sw=0
			}
			r.accumulate(lfs_XCH4 , 0., +tmp*cell_volume);

			// NCP -> gas phase
			tmp = 0.;
			auto YCH4_alg = param.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,S,zCH4);
			if( ( Sg - ( 1. - YCH4_alg - YH2O ) ) > eps_ap ){//active set.
				tmp +=  1. - YCH4_alg - YH2O ;//Active => phase is present => summation condition holds
//				std::cout<< "alpha_vol YH2O: " << tmp << std::endl;
			}else{
				tmp += Sg;// inactive set. Inactive => phase is absent => Sg=0
			}
			r.accumulate(lfs_YH2O , 0., +tmp*cell_volume);

		  
	  }
	  
	  
	  // skeleton integral depending on test and ansatz functions
	  // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig,
			  	  	  	   const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
						   const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
						   R& r_s, R& r_n) const
	  {
			typedef typename LFSU::template Child<Indices::PVId_Pw>::Type LFS_Pw_s;
			const LFS_Pw_s& lfs_Pw_s = lfsu_s.template child<Indices::PVId_Pw>();
			typedef typename LFSU::template Child<Indices::PVId_Sg>::Type LFS_Sg_s;
			const LFS_Sg_s& lfs_Sg_s = lfsu_s.template child<Indices::PVId_Sg>();
			typedef typename LFSU::template Child<Indices::PVId_XCH4>::Type LFS_XCH4_s;
			const LFS_XCH4_s& lfs_XCH4_s = lfsu_s.template child<Indices::PVId_XCH4>();
			typedef typename LFSU::template Child<Indices::PVId_YH2O>::Type LFS_YH2O_s;
			const LFS_YH2O_s& lfs_YH2O_s = lfsu_s.template child<Indices::PVId_YH2O>();

			typedef typename LFSU::template Child<Indices::PVId_Pw>::Type LFS_Pw_n;
			const LFS_Pw_n& lfs_Pw_n = lfsu_n.template child<Indices::PVId_Pw>();
			typedef typename LFSU::template Child<Indices::PVId_Sg>::Type LFS_Sg_n;
			const LFS_Sg_n& lfs_Sg_n = lfsu_n.template child<Indices::PVId_Sg>();
			typedef typename LFSU::template Child<Indices::PVId_XCH4>::Type LFS_XCH4_n;
			const LFS_XCH4_n& lfs_XCH4_n = lfsu_n.template child<Indices::PVId_XCH4>();
			typedef typename LFSU::template Child<Indices::PVId_YH2O>::Type LFS_YH2O_n;
			const LFS_YH2O_n& lfs_YH2O_n = lfsu_n.template child<Indices::PVId_YH2O>();

	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
	        const auto& cell_outside = ig.outside();
			const IndexSet &indexSet = gv.indexSet();
			int inside_cell_number  = indexSet.index(cell_inside);
			int outside_cell_number = indexSet.index(cell_outside);

	        // get geometries
	        auto geo = ig.geometry();
	        const auto dim = geo.mydimension;
	        auto geo_inside  = cell_inside.geometry();
	        auto geo_outside = cell_outside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto ref_el_outside = referenceElement(geo_outside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto outside_cell_center_local 	= ref_el_outside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();
	        auto outside_cell_center_global = geo_outside.center();

	        // distance of cell centers
	        auto d = outside_cell_center_global;
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);

			// compute primary vars at local self and neighbour centers
	        double Pw_s   = x_s(lfs_Pw_s,  0);
			double Sg_s   = x_s(lfs_Sg_s,  0);
			double XCH4_s = x_s(lfs_XCH4_s,0);
			double YH2O_s = x_s(lfs_YH2O_s,0);

	        double Pw_n   = x_n(lfs_Pw_n,  0);
			double Sg_n   = x_n(lfs_Sg_n,  0);
			double XCH4_n = x_n(lfs_XCH4_n,0);
			double YH2O_n = x_s(lfs_YH2O_n,0);

			// reference salinity and temp.
			auto S_s = param.parameter.ReferenceSalinity();
			auto Xc_s = param.parameter.ReferenceSaltConcentration();
			auto T_s /*ndim*/ = param.parameter.ReferenceTemperature()/Xc_T;
			auto S_n = S_s;
			auto Xc_n = Xc_s;
			auto T_n = T_s;

	        // compute the old porosity and partial Pg (ppg) at self and neighbour elements
	        double porosity_old_s=0.;
	        eval_porosity_old->evalFunction(cell_inside,inside_cell_center_local,&porosity_old_s);
	        double porosity_old_n=0.;
	        eval_porosity_old->evalFunction(cell_outside,outside_cell_center_local,&porosity_old_n);
	        double Pw_old_s=0.;
	        eval_pw_old->evalFunction(cell_inside,inside_cell_center_local,&Pw_old_s);
	        double Pw_old_n=0.;
	        eval_pw_old->evalFunction(cell_outside,outside_cell_center_local,&Pw_old_n);
	        double Sg_old_s=0.;
	        eval_sg_old->evalFunction(cell_inside,inside_cell_center_local,&Sg_old_s);
	        double Sg_old_n=0.;
	        eval_sg_old->evalFunction(cell_outside,outside_cell_center_local,&Sg_old_n);
	        double Sw_old_s=1.-Sg_old_s;
	        double Sw_old_n=1.-Sg_old_n;
			auto Pc_old_s = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_old_s, porosity_old_s );
			auto Pc_old_n = param.hydraulicProperty.CapillaryPressure( cell_outside,outside_cell_center_local, Sw_old_n, porosity_old_n );
			auto Pg_old_s = Pw_old_s + Pc_old_s;
			auto Pg_old_n = Pw_old_n + Pc_old_n;
			double Pg_tmp_s = Pw_s+param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, 1.-Sg_s, porosity_old_s );
			double Pg_tmp_n = Pw_n+param.hydraulicProperty.CapillaryPressure( cell_outside,outside_cell_center_local, 1.-Sg_n, porosity_old_n );
//			double dppg_s = ( std::max(Sg_s,0.)*Pg_tmp_s-std::max(Sg_old_s,0.)*Pg_old_s ) /(*dt);
//			dppg_s *= (Xc_P/Xc_t);
//			double dppg_n = ( std::max(Sg_n,0.)*Pg_tmp_n-std::max(Sg_old_n,0.)*Pg_old_n ) /(*dt);
//			dppg_n *= (Xc_P/Xc_t);
//			double dppg_s = ( std::max(Sg_s,0.)*Pg_tmp_s-std::max(Sg_old_s,0.)*Pg_old_s );
			double dppg_s = ( std::max(Sg_s,0.)*Pg_tmp_s+(1.-Sg_s)*Pw_s ) - ( std::max(Sg_old_s,0.)*Pg_old_s+(1.-Sg_old_s)*Pw_old_s );
			dppg_s *= Xc_P;
//			double dppg_n = ( std::max(Sg_n,0.)*Pg_tmp_n-std::max(Sg_old_n,0.)*Pg_old_n );
			double dppg_n = ( std::max(Sg_n,0.)*Pg_tmp_n+(1.-Sg_n)*Pw_n ) - ( std::max(Sg_old_n,0.)*Pg_old_n+(1.-Sg_old_n)*Pw_old_n );
			dppg_n *= Xc_P;
			
			// compute properties at local cell center and neighbour cell center

			// Hydraulic properties
			// 1. porosity
			auto porosity_s = param.hydraulicProperty.ScaledPorosity(cell_inside, inside_cell_center_local,  porosity_old_s, dppg_s );
			auto porosity_n = param.hydraulicProperty.ScaledPorosity(cell_outside,outside_cell_center_local, porosity_old_n, dppg_n );
			// 2. abs. permeability
			auto permeability_s =  param.soil.SedimentPermeabilityVector( cell_inside,inside_cell_center_local ) * ig.unitOuterNormal(face_center_local) ;
			permeability_s *= param.hydraulicProperty.PermeabilityScalingFactor( cell_inside,inside_cell_center_local, porosity_s );
			permeability_s = abs(permeability_s);
			auto permeability_n =  param.soil.SedimentPermeabilityVector( cell_inside,inside_cell_center_local ) * ig.unitOuterNormal(face_center_local);
			permeability_n *= param.hydraulicProperty.PermeabilityScalingFactor( cell_outside,outside_cell_center_local, porosity_s );
			permeability_n = abs(permeability_n);
			// 3. tortuosity
			auto tau_s = param.soil.Tortuosity( porosity_s );
			auto tau_n = param.soil.Tortuosity( porosity_n );

			// 4. secondary variables
			auto Sw_s = 1.-Sg_s;
			auto Sw_n = 1.-Sg_n;
			auto Pc_s = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_s, porosity_s );
			auto Pc_n = param.hydraulicProperty.CapillaryPressure( cell_outside,outside_cell_center_local, Sw_n, porosity_n );
			auto Pg_s = Pw_s + Pc_s;
			auto Pg_n = Pw_n + Pc_n;
			auto Peff_s = Pg_s*Sg_s + Pw_s*Sw_s;
			auto Peff_n = Pg_n*Sg_n + Pw_n*Sw_n;
			auto zCH4_s = param.eos.EvaluateCompressibilityFactor( T_s*Xc_T,Pg_s*Xc_P );
			auto zCH4_n = param.eos.EvaluateCompressibilityFactor( T_n*Xc_T,Pg_n*Xc_P );

			// 5. relative permeability
			// self
			auto kr_g_s = param.hydraulicProperty.krg( cell_inside, inside_cell_center_local, Sw_s );
			auto kr_w_s = param.hydraulicProperty.krw( cell_inside, inside_cell_center_local, Sw_s );
			// neighbour
			auto kr_g_n = param.hydraulicProperty.krg( cell_outside, outside_cell_center_local, Sw_n );
			auto kr_w_n = param.hydraulicProperty.krw( cell_outside, outside_cell_center_local, Sw_n );

			// 6. dynamic viscosity
			// self
			auto mu_g_s = param.gas.DynamicViscosity( T_s*Xc_T, Pg_s*Xc_P );
			auto mu_w_s = param.water.DynamicViscosity( T_s*Xc_T, Pw_s*Xc_P, S_s );
			// neighbour
			auto mu_g_n = param.gas.DynamicViscosity( T_n*Xc_T, Pg_n*Xc_P );
			auto mu_w_n = param.water.DynamicViscosity( T_n*Xc_T, Pw_n*Xc_P, S_n );

			// 7. density
			// self
			auto rho_mass_g_s = param.gas.Density( T_s*Xc_T, Pg_s*Xc_P, zCH4_s );
			auto rho_mass_w_s = param.water.Density( T_s*Xc_T, Pw_s*Xc_P, S_s );
			// neighbour
			auto rho_mass_g_n = param.gas.Density( T_n*Xc_T, Pg_n*Xc_P, zCH4_n );
			auto rho_mass_w_n = param.water.Density( T_n*Xc_T, Pw_n*Xc_P, S_n );

			// 8. Diffusion coefficient
			// self
			auto DH2O_g_s	= tau_s * porosity_s * Sg_s * param.mixture.DiffCoeffH2OInGas( T_s*Xc_T,Pg_s*Xc_P );
			auto DCH4_w_s	= tau_s * porosity_s * Sw_s * param.mixture.DiffCoeffCH4InLiquid( T_s*Xc_T,Pw_s*Xc_P );
			// neighbour
			auto DH2O_g_n	= tau_n * porosity_n * Sg_n * param.mixture.DiffCoeffH2OInGas( T_n*Xc_T,Pg_n*Xc_P );
			auto DCH4_w_n	= tau_n * porosity_n * Sw_n * param.mixture.DiffCoeffCH4InLiquid( T_n*Xc_T,Pw_n*Xc_P );
			// interface values
			auto DH2O_g_int = 2.*DH2O_g_s*DH2O_g_n/( DH2O_g_s+DH2O_g_n + 1.e-15 );
			auto DCH4_w_int = 2.*DCH4_w_s*DCH4_w_n/( DCH4_w_s+DCH4_w_n + 1.e-15 );


//			//sediment velocity
			auto normalvelocity_s = param.parameter.SedimentVelocity( (*time)*Xc_t, (*dt)*Xc_t )
								  * ig.unitOuterNormal(face_center_local) ;
			normalvelocity_s *= Xc_vs;

			// upwinding wrt gas-phase velocity
			auto gravity = param.parameter.g() * ig.unitOuterNormal(face_center_local) ;
			gravity *= Xc_grav;
			auto normalpotential_g = (Pg_n - Pg_s)/distance
								   + ( 0.5*(rho_mass_g_s + rho_mass_g_n) ) * gravity ;
			double omegaup_g_s = 0., omegaup_g_n = 0.;
			if( normalpotential_g>0.){
				omegaup_g_s = 0.;
				omegaup_g_n = 1.;
			}else{
				omegaup_g_s = 1.;
				omegaup_g_n = 0.;
			}
			auto normalvelocity_g = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
									* ( omegaup_g_s * kr_g_s/mu_g_s + omegaup_g_n * kr_g_n/mu_g_n )
									* normalpotential_g;
			normalvelocity_g += porosity_s*Sg_s * normalvelocity_s;

			//upwinding wrt water-phase velocity
			auto normalpotential_w = (Pw_n - Pw_s)/distance
								   + ( 0.5*(rho_mass_w_s + rho_mass_w_n) ) * gravity ;
			double omegaup_w_s = 0., omegaup_w_n = 0.;
			if( normalpotential_w>0.){
				omegaup_w_s = 0.;
				omegaup_w_n = 1.;
			}else{
				omegaup_w_s = 1.;
				omegaup_w_n = 0.;
			}
			auto normalvelocity_w = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
									* ( omegaup_w_s * kr_w_s/mu_w_s + omegaup_w_n * kr_w_n/mu_w_n )
									* normalpotential_w;
			normalvelocity_w += porosity_s*Sw_s * normalvelocity_s;


			//convective flux terms
			auto convectiveflux_CH4 = ( omegaup_g_s*rho_mass_g_s*(1.-YH2O_s) + omegaup_g_n*rho_mass_g_n*(1.-YH2O_n) ) * normalvelocity_g
									+ ( omegaup_w_s*rho_mass_w_s*XCH4_s + omegaup_w_n*rho_mass_w_n*XCH4_n ) * normalvelocity_w;

			auto convectiveflux_H2O = ( omegaup_g_s*rho_mass_g_s*YH2O_s + omegaup_g_n*rho_mass_g_n*YH2O_n ) * normalvelocity_g
									+ ( omegaup_w_s*rho_mass_w_s*(1.-XCH4_s) + omegaup_w_n*rho_mass_w_n*(1.-XCH4_n) ) * normalvelocity_w;

			//diffusive flux terms
			double j_CH4_w  = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * DCH4_w_int * ( XCH4_n - XCH4_s )/distance;
			double j_H2O_w  = - j_CH4_w;
			double j_H2O_g  = - 0.5*( rho_mass_g_s + rho_mass_g_n ) * DH2O_g_int * ( YH2O_n - YH2O_s )/distance;
			double j_CH4_g  = - j_H2O_g;

			auto diffusiveflux_CH4 = j_CH4_w + j_CH4_g;

			auto diffusiveflux_H2O = j_H2O_w + j_H2O_g;
			

			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_CH4 + Xc_diff_m * diffusiveflux_CH4 ;
			r_s.accumulate(lfs_Sg_s, 0, +tmp*face_volume);
			r_n.accumulate(lfs_Sg_n, 0, -tmp*face_volume);

			// H2O-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_H2O + Xc_diff_m * diffusiveflux_H2O ;
			r_s.accumulate(lfs_Pw_s, 0, +tmp*face_volume);
			r_n.accumulate(lfs_Pw_n, 0, -tmp*face_volume);
			
	  }
	  
	  
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary ( const IG& ig,
			  	  	  	  	const LFSU& lfsu, const X& x, const LFSV& lfsv,
							R& r ) const
	  {
			typedef typename LFSU::template Child<Indices::PVId_Pw>::Type LFS_Pw_s;
			const LFS_Pw_s& lfs_Pw_s = lfsu.template child<Indices::PVId_Pw>();
			typedef typename LFSU::template Child<Indices::PVId_Sg>::Type LFS_Sg_s;
			const LFS_Sg_s& lfs_Sg_s = lfsu.template child<Indices::PVId_Sg>();
			typedef typename LFSU::template Child<Indices::PVId_XCH4>::Type LFS_XCH4_s;
			const LFS_XCH4_s& lfs_XCH4_s = lfsu.template child<Indices::PVId_XCH4>();
			typedef typename LFSU::template Child<Indices::PVId_YH2O>::Type LFS_YH2O_s;
			const LFS_YH2O_s& lfs_YH2O_s = lfsu.template child<Indices::PVId_YH2O>();

	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
			const IndexSet &indexSet = gv.indexSet();
			int inside_cell_number  = indexSet.index(cell_inside);

	        // get geometries
	        auto geo = ig.geometry();
	        const auto dim = geo.mydimension;
	        auto geo_inside = cell_inside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);

	        // distance of cell centers
	        auto d = geo.global(face_center_local);
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

			// compute primary vars at local self and neighbour centers
	        double Pw_s   = x(lfs_Pw_s,  0);
			double Sg_s   = x(lfs_Sg_s,  0);
			double XCH4_s = x(lfs_XCH4_s,0);
			double YH2O_s = x(lfs_YH2O_s,0);

			double Pw_n   = Pw_s;
			double Sg_n   = Sg_s;
			double XCH4_n = XCH4_s;
			double YH2O_n = YH2O_s;
			
			// reference salinity and temp.
			auto S_s = param.parameter.ReferenceSalinity();
			auto Xc_s = param.parameter.ReferenceSaltConcentration();
			auto T_s /*ndim*/ = param.parameter.ReferenceTemperature()/Xc_T;
			auto S_n = S_s;
			auto Xc_n = Xc_s;
			auto T_n = T_s;

	        // compute the old porosity and partial Pd (ppg) at self and neighbour elements
	        double porosity_old_s=0.;
	        eval_porosity_old->evalFunction(cell_inside,inside_cell_center_local,&porosity_old_s);
	        double porosity_old_n=porosity_old_s;
	        double Pw_old_s=0.;
	        eval_pw_old->evalFunction(cell_inside,inside_cell_center_local,&Pw_old_s);
	        double Pw_old_n=Pw_old_s;
	        double Sg_old_s=0.;
	        eval_sg_old->evalFunction(cell_inside,inside_cell_center_local,&Sg_old_s);
	        double Sg_old_n=Sg_old_s;
			auto Pg_old_s = Pw_old_s + param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, 1.-Sg_old_s, porosity_old_s );
			double Pg_tmp_s=Pw_s + param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, 1.-Sg_s, porosity_old_s );
//			double dppg_s = ( std::max(Sg_s,0.)*Pg_s-std::max(Sg_old_s,0.)*Pg_old_s ) /(*dt);
//			dppg_s *= (Xc_P/Xc_t);
//			double dppg_s = ( std::max(Sg_s,0.)*Pg_tmp_s-std::max(Sg_old_s,0.)*Pg_old_s );
			double dppg_s = ( std::max(Sg_s,0.)*Pg_tmp_s+(1.-Sg_s)*Pw_s ) - ( std::max(Sg_old_s,0.)*Pg_old_s+(1.-Sg_old_s)*Pw_old_s );
			dppg_s *= Xc_P;


			double newtime = *time + *dt;
			newtime *= Xc_t;
			// evaluate boundary condition types for {Pw,Sg} or {Fw,Fg} 
			auto bctype = bc.type( ig , face_center_local, (*time)*Xc_t, (*dt)*Xc_t ) ;

			// evaluate boundary condition values for {Pw,Sg} or {Fw,Fg} 
			auto bcvalue = bc.value( ig , face_center_local, (*time)*Xc_t, (*dt)*Xc_t ) ;


			/*TODO: documentation for boundary conditions*/
			double totalflux_CH4 = 0.;
			double totalflux_H2O = 0.;
			
			if( (bctype[Indices::BC_w] == Indices::dirichlet) and (bctype[Indices::BC_g] == Indices::dirichlet) ) {
				Pw_n = bcvalue[Indices::BC_w];
				Sg_n = bcvalue[Indices::BC_g];
				
				// compute properties at local cell center and neighbour cell center
				// Hydraulic properties for self and neighbour
				// 1. porosity
				auto porosity_s = param.hydraulicProperty.ScaledPorosity(cell_inside, inside_cell_center_local,  porosity_old_s, dppg_s );
				auto porosity_n = porosity_s;
				// 2. abs. permeability
				auto permeability_s =  param.soil.SedimentPermeabilityVector( cell_inside,inside_cell_center_local ) * ig.unitOuterNormal(face_center_local);
				permeability_s *= param.hydraulicProperty.PermeabilityScalingFactor( cell_inside,inside_cell_center_local, porosity_s );
				permeability_s = abs(permeability_s);
				auto permeability_n = permeability_s;
				// 3. tortuosity
				auto tau_s = param.soil.Tortuosity( porosity_s );
				auto tau_n = param.soil.Tortuosity( porosity_n );

				// 4. secondary variables
				// at self
				auto Sw_s = 1.-Sg_s;
				auto Pc_s = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_s, porosity_s );
				auto Pg_s = Pw_s + Pc_s;
				auto zCH4_s = param.eos.EvaluateCompressibilityFactor( T_s*Xc_T,Pg_s*Xc_P );
				// at neighbour
				auto Sw_n = 1.-Sg_n;
				auto Pc_n = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_n, porosity_n );
				auto Pg_n = Pw_n + Pc_n;
				auto zCH4_n = param.eos.EvaluateCompressibilityFactor( T_n*Xc_T,Pg_n*Xc_P );

				// 5. relative permeability
				// self
				auto kr_g_s = param.hydraulicProperty.krg( cell_inside, inside_cell_center_local, Sw_s );
				auto kr_w_s = param.hydraulicProperty.krw( cell_inside, inside_cell_center_local, Sw_s );
				// neighbour
				auto kr_g_n = param.hydraulicProperty.krg( cell_inside, inside_cell_center_local, Sw_n );
				auto kr_w_n = param.hydraulicProperty.krw( cell_inside, inside_cell_center_local, Sw_n );

				// 6. dynamic viscosity
				// self
				auto mu_g_s = param.gas.DynamicViscosity( T_s*Xc_T, Pg_s*Xc_P );
				auto mu_w_s = param.water.DynamicViscosity( T_s*Xc_T, Pw_s*Xc_P, S_s );
				// neighbour
				auto mu_g_n = param.gas.DynamicViscosity( T_n*Xc_T, Pg_n*Xc_P );
				auto mu_w_n = param.water.DynamicViscosity( T_n*Xc_T, Pw_n*Xc_P, S_n );

				// 7. density
				// self
				auto rho_mass_g_s = param.gas.Density( T_s*Xc_T, Pg_s*Xc_P, zCH4_s );
				auto rho_mass_w_s = param.water.Density( T_s*Xc_T, Pw_s*Xc_P, S_s );
				// neighbour
				auto rho_mass_g_n = param.gas.Density( T_n*Xc_T, Pg_n*Xc_P, zCH4_n );
				auto rho_mass_w_n = param.water.Density( T_n*Xc_T, Pw_n*Xc_P, S_n );

				// 8. Diffusion coefficient
				// self
				auto DH2O_g_s	= tau_s * porosity_s * Sg_s * param.mixture.DiffCoeffH2OInGas( T_s*Xc_T,Pg_s*Xc_P );
				auto DCH4_w_s	= tau_s * porosity_s * Sw_s * param.mixture.DiffCoeffCH4InLiquid( T_s*Xc_T,Pw_s*Xc_P );
				// neighbour
				auto DH2O_g_n	= tau_n * porosity_n * Sg_n * param.mixture.DiffCoeffH2OInGas( T_n*Xc_T,Pg_n*Xc_P );
				auto DCH4_w_n	= tau_n * porosity_n * Sw_n * param.mixture.DiffCoeffCH4InLiquid( T_n*Xc_T,Pw_n*Xc_P );
				// interface values
				auto DH2O_g_int = 2.*DH2O_g_s*DH2O_g_n/( DH2O_g_s+DH2O_g_n + 1.e-15 );
				auto DCH4_w_int = 2.*DCH4_w_s*DCH4_w_n/( DCH4_w_s+DCH4_w_n + 1.e-15 );

	//			//sediment velocity
				auto normalvelocity_s = 0.;//
				normalvelocity_s *= Xc_vs;
	//					param.parameter.SedimentVelocity( (*time), (*dt) )
	//								  * ig.unitOuterNormal(face_center_local) ;

				// upwinding wrt gas-phase velocity
				auto gravity = param.parameter.g() * ig.unitOuterNormal(face_center_local) ;
				gravity *= Xc_grav;
				auto normalpotential_g = (Pg_n - Pg_s)/distance
									   + ( 0.5*(rho_mass_g_s + rho_mass_g_n) ) * gravity ;
				double omegaup_g_s = 0., omegaup_g_n = 0.;
				if( normalpotential_g>0.){
					omegaup_g_s = 0.;
					omegaup_g_n = 1.;
				}else{
					omegaup_g_s = 1.;
					omegaup_g_n = 0.;
				}

				//upwinding wrt water-phase velocity
				auto normalpotential_w = (Pw_n - Pw_s)/distance
									   + ( 0.5*(rho_mass_w_s + rho_mass_w_n) ) * gravity ;
				double omegaup_w_s = 0., omegaup_w_n = 0.;
				if( normalpotential_w>0.){
					omegaup_w_s = 0.;
					omegaup_w_n = 1.;
				}else{
					omegaup_w_s = 1.;
					omegaup_w_n = 0.;
				}

				double normalvelocity_g = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
									    * ( omegaup_g_s * kr_g_s/mu_g_s + omegaup_g_n * kr_g_n/mu_g_n )
									    * normalpotential_g;
				double normalvelocity_w = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
									    * ( omegaup_w_s * kr_w_s/mu_w_s + omegaup_w_n * kr_w_n/mu_w_n )
									    * normalpotential_w;
				normalvelocity_g += porosity_s*Sg_s * normalvelocity_s;
				normalvelocity_w += porosity_s*Sw_s * normalvelocity_s;
				
				
				/* FLUX TERMS */

				//convective flux terms
				auto convectiveflux_CH4 = ( omegaup_g_s*rho_mass_g_s*(1.-YH2O_s) + omegaup_g_n*rho_mass_g_n*(1.-YH2O_n) ) * normalvelocity_g
										+ ( omegaup_w_s*rho_mass_w_s*XCH4_s + omegaup_w_n*rho_mass_w_n*XCH4_n ) * normalvelocity_w;

				auto convectiveflux_H2O = ( omegaup_g_s*rho_mass_g_s*YH2O_s + omegaup_g_n*rho_mass_g_n*YH2O_n ) * normalvelocity_g
										+ ( omegaup_w_s*rho_mass_w_s*(1.-XCH4_s) + omegaup_w_n*rho_mass_w_n*(1.-XCH4_n) ) * normalvelocity_w;

				//diffusive flux terms
				double j_CH4_w  = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * DCH4_w_int * ( XCH4_n - XCH4_s )/distance;
				double j_H2O_w  = - j_CH4_w;
				double j_H2O_g  = - 0.5*( rho_mass_g_s + rho_mass_g_n ) * DH2O_g_int * ( YH2O_n - YH2O_s )/distance;
				double j_CH4_g  = - j_H2O_g;

				auto diffusiveflux_CH4 = j_CH4_w + j_CH4_g;

				auto diffusiveflux_H2O = j_H2O_w + j_H2O_g;
				
				//total fluxes
				totalflux_CH4 /*ndim*/ = Xc_conv_m * convectiveflux_CH4 + Xc_diff_m * diffusiveflux_CH4;
				totalflux_H2O /*ndim*/ = Xc_conv_m * convectiveflux_H2O + Xc_diff_m * diffusiveflux_H2O;
				
			}else if( (bctype[Indices::BC_w] == Indices::neumann) and (bctype[Indices::BC_g] == Indices::neumann) ){
				totalflux_H2O /*ndim*/ = bcvalue[Indices::BC_w];
				totalflux_CH4 /*ndim*/ = bcvalue[Indices::BC_g];
			}else{
				std::cout<< __FILE__ << '\t' << __LINE__ << '\n' << "Invalid Boundary Conditions." << std::endl;
				exit(0);
			}


			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp = totalflux_CH4 ;
			r.accumulate(lfs_Sg_s, 0, +tmp*face_volume);

			// H2O-component-wise mass-balance
			tmp = totalflux_H2O ;
			r.accumulate(lfs_Pw_s, 0, +tmp*face_volume);

	  }
	
};
