class BaseEoS{
private:
	constexpr static double eps = 1.e-6;
public:

	void printName(){
		std::cout<<"EOS: Base (z=1) " << std::endl;
	}

	double EvaluateCompressibilityFactor( double T/*K*/, double P/*Pa*/ )const{

		double compressibilityFactor = 1.0;
		return compressibilityFactor ;

	}
};


class FrauenhoferFunction{
private:
	constexpr static double eps = 1.e-6;
	constexpr static int iterMax = 10;
	Methane methane;
public:

	void printName(){
		std::cout<<"EOS: FRAUENHOFER COMSOL MODEL, empirical " << std::endl;
	}

	double EvaluateCompressibilityFactor( double T/*K*/, double P/*Pa*/ )const{

		double C0 = 1.0003;
		double C1 = 2.1173e-8;
		double C2 = 5.0207e-16;
		double C3 = 7.5669e-23;
		double C4 = 1.3259e-30;

		double compressibilityFactor = C0 - C1*P - C2*P*P + C3*P*P*P - C4*P*P*P*P;
		return compressibilityFactor ;

	}
};

class PengRobinson{
private:
	constexpr static double eps = 1.e-6;
	constexpr static int iterMax = 10;
	Methane methane;
public:

	void printName(){
		std::cout<<"EOS: PENG-ROBINSON " << std::endl;
	}

	std::vector<double> EvaluateEoSParams( double T/*K*/, double P/*Pa*/ )const{

		std::vector<double> PREoSParams(2);
		for(int i =0; i<PREoSParams.size();i++){
			PREoSParams[i] = 0.;
		}

		double omega 	= methane.AccentricityFactor();
		double Tc/*K*/	= methane.CriticalTemperature();
		double Pc/*Pa*/	= methane.CriticalPressure();
		double R_u/*J.mol^-1.K^-1*/ = methane.Ru;

		double kappa = 0.;
		if( omega <= 0.49 ){
			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
		}
		else{
			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
		}

		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );

		double a_T = 0.45724 * ( ( R_u * R_u * Tc * Tc ) / Pc ) * ac ;
		double b = 0.07780 * ( R_u * Tc / Pc );

		PREoSParams[0] = a_T;
		PREoSParams[1] = b  ;

		return PREoSParams;

	}

	std::vector<double> PolynomialCoefficients( double T/*K*/, double P/*Pa*/ )const{

		std::vector<double> EoSParameters(2);
		for(int i =0; i<EoSParameters.size();i++){
			EoSParameters[i] = 0.;
		}
		EoSParameters = EvaluateEoSParams( T,P );

		double a_T = EoSParameters[0];
		double b   = EoSParameters[1];
		double R_u = methane.Ru;

		double A = ( a_T * P ) / ( pow( ( R_u * T ) , 2 ) ) ;
		double B = ( b * P ) / ( R_u * T ) ;

		std::vector<double> Coeffs(4);
		for(int i =0; i<Coeffs.size();i++){
			Coeffs[i] = 0.;
		}
		// cubic equation: a0*z^3 + a1*z^2 + a2*z + a3

		Coeffs[0] = 1;
		Coeffs[1] = B - 1. ;
		Coeffs[2] = A - 2.*B - 3.*B*B ;
		Coeffs[3] = B*B*B + B*B - A*B ;

		return Coeffs;

	}

	double Method( std::vector<double> Coeffs )const{

//		std::cout <<"**NEWTON RAPHSON METHOD CALLED FOR Z_CH4 CALCULATION **"<< std::endl;

		int polynomialOrder = Coeffs.size() - 1;
//		std::cout<< "PolynomialOrder = " << polynomialOrder << std::endl;

		double z = 1.;
		double z_up = 1.;
		double defect = 0.;

		int iter = 0;
		do{
			iter += 1;
			double f_z = 0.;
			double df_z = 0.;
			for( int i = 0 ; i<=polynomialOrder ; i++ ){
				f_z += Coeffs[i] * pow( z , polynomialOrder - i );
				df_z += ( polynomialOrder - i ) * Coeffs[i] * pow( z , (polynomialOrder - i) -1 ) ;
			}
			z_up = z - (f_z / df_z );
			defect = z_up - z;

//			std::cout<< "iteration number = " << iter << "  , defect = " << defect << "  , z = " << z_up << std::endl;

			z = z_up;
		}
		while( std::abs(defect)>eps and iter<iterMax );

//		std::cout <<"z_CH4 CALCULATED = " << z_up << std::endl;

		return z_up;
	}

	double EvaluateCompressibilityFactor( double T/*K*/, double P/*Pa*/ )const{

		std::vector<double> Coeffs(4);
		for(int i =0; i<Coeffs.size();i++){
			Coeffs[i] = 0.;
		}
		Coeffs = PolynomialCoefficients(T,P);
		double compressibilityFactor = Method( Coeffs );
		return std::max(compressibilityFactor,0.1) ;
	}

};