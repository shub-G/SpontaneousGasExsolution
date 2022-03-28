class Methane
{
private:
	CharacteristicValues characteristicValue;
	
public:

	/* http://en.wikipedia.org/wiki/Gas_constant */
	constexpr static double Ru = 8.314462175; /* [J*mol^-1*K^-1] */

	double MolarMass() const {
		return 16.04 * 1.0e-3; 	/* [kg/mol] */
	}

	double AccentricityFactor()const {
		return 0.011 	 ;
	}

	double CriticalTemperature( ) const {
		return -82.7 + 273.15 ; /* [K] */
	}

	double CriticalPressure( ) const {
		return 45.96 * 1.0e5 ; /* [Pa] */
	}

	double Density(double T/*K*/, double Pg/*Pa*/, double z_CH4) const {

		double rho;
		/* rho: unit -> kg/m^3 */

		double R_CH4 = Ru/MolarMass();
		rho = Pg/( z_CH4 * R_CH4 * T);

		return rho/characteristicValue.density_c; /*ndim*/
	}

	double DynamicViscosity(double T/*K*/, double Pg/*Pa*/) const {

		double mu;
		/* mu: unit -> Pa.s */

		// Sutherland Correlation:
		// ref: http://portal.tpu.ru/SHARED/n/NATASHA/Material/Tab3/Glava_1.pdf
		double C = 162; // empirical constant
		double mu_0 = 1.0707e-5;
		// ref for mu_0 :http://www.pipeflowcalculations.com/tables/gas.php
		mu_0 *= (  1. - (1./(1.0707e-5)) * (  4.8134e-14 * Pg
											- 4.1719e-20 * Pg * Pg
											+ 7.3232e-28 * Pg * Pg * Pg )
				); // Pa.s -> ref: Frauenhofer Comsol Model
		mu = mu_0 * (273.15 + C ) * ( pow( (T/273.15), 1.5) / ( T + C ) ) ;

		return mu/characteristicValue.viscosity_c;/*ndim*/

	}

	double SolubilityCoefficient( double T/*K*/, double S ) const {

		double kHenry;
		/* [Pa^-1] */

		// REF: SUGAR TOOLBOX
		Water water;
		double Tc /*[K]*/   = water.CriticalTemperature();	//critical temperature of water
		double Pc /*[MPa]*/ = water.CriticalPressure()/1e6;	//critical pressure of water
		double tau = 1 - T/Tc;    							//dimensionless temperature
		double Ts = T/Tc;         							//reduced temperature

		// vapor pressure of water
		double a1 = -7.85951783;
		double a2 =  1.84408259;
		double a3 = -11.7866497;
		double a4 =  22.6807411;
		double a5 = -15.9618719;
		double a6 =  1.80122502;
		double pvap /*[MPa]*/ = Pc * exp( Tc/T * ( a1*tau + a2*pow(tau,1.5) + a3*pow(tau,3) + a4*pow(tau,3.5) + a5*pow(tau,4) + a6*pow(tau,7.5) ) );

		// Henry constant
		double A = -11.0094;
		double B = 4.8362;
		double C = 12.5220;
		kHenry /*MPa*/ = exp( log(pvap) + A/Ts + B*pow((1 - Ts),0.355)/Ts + C*exp(1 - Ts)*pow(Ts,(-0.41)) );
		kHenry *= 1.e6; /* [Pa] */

		return kHenry/characteristicValue.P_c; /*ndim*/
	}

};
