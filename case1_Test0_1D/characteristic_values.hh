class CharacteristicValues{

public:

	constexpr static double x_c = 1.;
	constexpr static double t_c = 1.;
	constexpr static double P_c = 1.;
	constexpr static double T_c = 1.;//dummy
	constexpr static double permeability_c = 1.;
	constexpr static double density_c = 1.;
	constexpr static double viscosity_c = 1.;
	constexpr static double dispersivity_c = 1.;

	constexpr static double X_source_mass 		= t_c/density_c;
	constexpr static double X_convective_mass 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );
	constexpr static double X_diffusive_mass 	= ( dispersivity_c * t_c ) / ( x_c*x_c ) ;
	constexpr static double X_gravity 			= density_c*x_c/P_c;
	constexpr static double X_solidvelocity 	= viscosity_c*x_c/(P_c*permeability_c);

};
