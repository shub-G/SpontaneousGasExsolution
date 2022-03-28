class Water
{
private:
	CharacteristicValues characteristicValue;

public:

	double CriticalTemperature( ) const {
		return 647.096 ; /* [K] */
	}

	double CriticalPressure( ) const {
		return 22.064 * 1.0e6 ; /* [Pa] */
	}

	double MolarMass( ) const {
		/* unit -> kg/mol */
		return 18.0/1000;
	}

	double Density( double T/*K*/, double Pw/*Pa*/, double S ) const {

		double rho;
		/* rho: unit -> kg/m^3 */

		/*
		 * averages values & expansion coefficients: ρ0=1027 kg/m^3,  T0=10°C,  S_0=35 g/kg
		 * Thermal expansion: \alpha_T=0.15 kg/(m^3 °C)
		 * Salinity contraction: \alpha_S=0.78 kg/(m^3 g/kg)
		 * Pressure compressibility: \alpha_P=0.0045 kg/(m^3 dbar)
		 * UNESCO EOS-80 : Equation of state for seawater
		 * We use a linear EOS (web ref:http://www.ccpo.odu.edu/~atkinson/OEAS405/Chapter2_Stratified_Ocean/Lec_04_DensityEOS.pdf)
		 */

		double rho_0 = 1027.0;
		double T_0 = 10.;
		double S_0 = 0.035;
		double alpha_T = -0.15;
		double alpha_S = 0.78*1e3;
		double alpha_P = 0.0045;

		rho = rho_0
			+ (   alpha_P*(Pw*1.e-4)
				+ alpha_T*((T-273.15)-T_0)
				+ alpha_S*(S-S_0)
			  );

		return rho/characteristicValue.density_c; /*ndim*/
	}

	double DynamicViscosity( double T/*K*/, double Pw/*Pa*/, double S ) const {
		double mu;
		/* mu: unit -> Pa.s */

#ifdef FRAUENHOFER_MODEL
		mu = 0.5 * 1.0e-3 ;
#else
		// REFERENCE:
		double mu_0 = 0.001792 ; // kg/m/s
		double a = - 1.94 ;
		double b = - 4.80 ;
		double c =  6.74 ;
		double T0 = 273.15 ; // K
		double Tr = T0/T;
		mu = mu_0 * exp( a + b * Tr + c * Tr*Tr );
#endif

		return mu/characteristicValue.viscosity_c;/*ndim*/
	}


	double SaturatedVaporPressure( double T/*K*/, double S ) const {

		double psat;   /* [Pa] */

		// REF: SUGAR TOOLBOX

		double Pc = CriticalPressure(); // in Pa
		double Tc = CriticalTemperature(); // in K
		double Tr = T/Tc;

		double c1 = -7.85951783;
		double c2 = 1.84408259;
		double c3 = -11.7866497;
		double c4 = 22.6807411;
		double c5 = -15.9618719;
		double c6 = 1.80122502;

		double lnppc = 1./Tr * (  c1*(1-Tr)
								+ c2*pow((1-Tr),1.5)
								+ c3*pow((1-Tr),3)
								+ c4*pow((1-Tr),3.5)
								+ c5*pow((1-Tr),4)
								+ c6*pow((1-Tr),7.5) );

		psat = Pc * exp(lnppc);   /* [Pa] */

		return psat/characteristicValue.P_c; /*ndim*/
	}

};