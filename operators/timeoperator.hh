template < class GV, class Params, class EVALPor, class EVALPw, class EVALSg >
class TimeOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume< TimeOperator<GV,Params,EVALPor,EVALPw,EVALSg> >,
    public Dune::PDELab::NumericalJacobianVolume	 < TimeOperator<GV,Params,EVALPor,EVALPw,EVALSg> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV& gv;
	const Params&	  param;
	EVALPor	*eval_porosity_old;
	EVALPw  *eval_pw_old;
	EVALSg  *eval_sg_old;
	double 	*time;
	double 	*dt;
	constexpr static double eps 	= 1.0e-6;
	constexpr static double eps_ap	= 0.;//1.e-9;
	double Xc_K ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_D;
	double Xc_P;
	double Xc_T;

public:
	  // pattern assembly flags
	  enum { doPatternVolume = true };

	  // residual assembly flags
	  enum { doAlphaVolume = true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor remembers parameters
	  TimeOperator( const GV& gv_, const Params& param_, EVALPor *eval_porosity_old_, EVALPw *eval_pw_old_, EVALSg *eval_sg_old_,double	*time_,double *dt_ )
	  :  gv(gv_), param(param_), eval_porosity_old(eval_porosity_old_), eval_pw_old(eval_pw_old_), eval_sg_old(eval_sg_old_), time( time_ ),dt( dt_ )
	  {
		  Xc_K 		= param.characteristicValue.permeability_c;
		  Xc_mu 	= param.characteristicValue.viscosity_c;
		  Xc_rho 	= param.characteristicValue.density_c;
		  Xc_D 		= param.characteristicValue.dispersivity_c;
		  Xc_P 		= param.characteristicValue.P_c;
		  Xc_T 		= param.characteristicValue.T_c;
	  }

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const {
		  
	        typedef typename LFSU::template Child< Indices::PVId_Pw	  >::Type LFS_Pw ;
	        const LFS_Pw&	lfs_Pw	= lfsu.template child< Indices::PVId_Pw  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_Sg   >::Type LFS_Sg;
	        const LFS_Sg&	lfs_Sg	= lfsu.template child< Indices::PVId_Sg  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_XCH4 >::Type LFS_XCH4;
	        const LFS_XCH4& lfs_XCH4= lfsu.template child< Indices::PVId_XCH4>() ;
	        typedef typename LFSU::template Child< Indices::PVId_YH2O >::Type LFS_YH2O;
	        const LFS_YH2O& lfs_YH2O= lfsu.template child< Indices::PVId_YH2O>() ;

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
	        auto cell_volume = geo.volume();

			// compute PVs at local center
			double Pw   = x(lfs_Pw  ,0);
			double Sg   = x(lfs_Sg  ,0);
			double XCH4 = x(lfs_XCH4,0);
			double YH2O = x(lfs_YH2O,0);

			// reference salinity and temp.
			auto S = param.parameter.ReferenceSalinity();
			auto Xc = param.parameter.ReferenceSaltConcentration();
			auto T /*ndim*/ = param.parameter.ReferenceTemperature()/Xc_T;
			
	        // compute the old porosity and partial Pg (ppg) at self and neighbour elements
	        double porosity_old=0.;
	        eval_porosity_old->evalFunction(cell,cell_center_local,&porosity_old);
	        double Pw_old=0.;
	        eval_pw_old->evalFunction(cell,cell_center_local,&Pw_old);
	        double Sg_old=0.;
	        eval_sg_old->evalFunction(cell,cell_center_local,&Sg_old);
	        double Pg_old=Pw_old+param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg_old, porosity_old );
			double Pg_tmp=Pw+param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, 1.-Sg, porosity_old );
//			double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg_old,0.)*Pg_old ) /(*dt) * (Xc_P/Xc_t);
//			double dppg = ( std::max(Sg,0.)*Pg_tmp-std::max(Sg_old,0.)*Pg_old ) * Xc_P;
			double dppg = ( std::max(Sg,0.)*Pg_tmp+(1.-Sg)*Pw ) - ( std::max(Sg_old,0.)*Pg_old+(1.-Sg_old)*Pw_old );
			dppg *= Xc_P;

			// secondary variables
			auto porosity = param.hydraulicProperty.ScaledPorosity(cell, cell_center_local, porosity_old,dppg );
			auto permeability = param.soil.SedimentPermeability( cell,cell_center_local )
							  * param.hydraulicProperty.PermeabilityScalingFactor( cell,cell_center_local, porosity );
			auto Sw = 1.-Sg;
			auto Pc = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw, porosity );
			auto Pg = Pw + Pc;
			auto zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );

			// compute properties at local center
			// density
			auto rho_mass_g = param.gas.Density( T*Xc_T, Pg*Xc_P, zCH4 );
			auto rho_mass_w = param.water.Density( T*Xc_T, Pw*Xc_P, S );

			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp = 0.;
			tmp += porosity * rho_mass_g * (1.-YH2O) * Sg ;
			tmp += porosity * rho_mass_w * XCH4 * Sw ;
//			std::cout<< "alpha-time 0: " << tmp << std::endl;
			r.accumulate(lfs_Sg, 0, +tmp*cell_volume);

			// H2O-component-wise mass-balance
			tmp = 0.;
			tmp += porosity * rho_mass_g * YH2O * Sg ;
			tmp += porosity * rho_mass_w * (1.-XCH4) * Sw ;
			r.accumulate(lfs_Pw, 0, +tmp*cell_volume);

	  }
	
};
