template<typename PTree>
class Mixture
{
private:
	const PTree& ptree;
	Methane methane;
	Water water;
	Parameters<PTree> parameter;
	CharacteristicValues X_c;

public:

	  //! constructor
	  Mixture (const PTree& ptree_)
	  :ptree(ptree_),
	   parameter(ptree_)
	  {}

	/* MOLE FRACTIONS ( X -> liquid ; Y -> gas ) */

	std::vector<double> EquilibriumMoleFractions( double T/*K*/, double Pg/*Pa*/, double S, double z )const{

		double Y_H2O, Y_CH4, X_H2O, X_CH4;

		// NOTE: it is not necessary to check state dependence for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double Xc = parameter.ReferenceSaltConcentration();
		double f_CH4 /*ndim*/ = z*(Pg/X_c.P_c)/methane.SolubilityCoefficient(T,S);
		double f_H2O /*ndim*/ = (Pg/X_c.P_c)/water.SaturatedVaporPressure( T,S );

		Y_H2O = ((1.-Xc)-f_CH4)/(f_H2O-f_CH4);
		Y_CH4 = 1.-Y_H2O;
		X_H2O = Y_H2O * f_H2O;
		X_CH4 = 1. - Xc - X_H2O;

		std::vector<double> X(Indices::numOfConcs,0.);
		X[Indices::concId_XCH4] = X_CH4;
		X[Indices::concId_XH2O] = X_H2O;
		X[Indices::concId_YCH4] = Y_CH4;
		X[Indices::concId_YH2O] = Y_H2O;

		return X;
	}

	double YCH4( double X_CH4, double T, double Pg, double S, double z )const{

		double Y_CH4;

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		Y_CH4 = X_CH4 * methane.SolubilityCoefficient(T,S) / ( z * Pg/X_c.P_c ) ;

		return Y_CH4;
	}

	double XH2O( double Y_H2O, double T, double Pg, double S )const{

		double X_H2O;

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		X_H2O = Y_H2O * (Pg/X_c.P_c) / water.SaturatedVaporPressure( T,S );

		return X_H2O;
	}

	/* MASS TRANSFER COEFFICIENTS */

	double DiffCoeffH2OInGas( double T/*K*/, double Pg/*Pa*/ ) const {

		double D;
		/* m^2/s */

		double a0 = 0. ;
		double a1 = 2.26e-9;
		double a2 = 0.002554;

		D = ( a0 + a1*T + a2/Pg ) ;

		return D/X_c.dispersivity_c; /*ndim*/

	}

	double DiffCoeffCH4InLiquid( double T/*K*/, double Pw/*Pa*/ ) const {

		double D; /* m^2/s */

		//reference: https://arxiv.org/pdf/1710.06375.pdf
		double A = 15.13*1000./8.31 ; 	/* K */
		double D0 = 0.7e-6;		/* cm^2/s */

		D = D0 * exp(-A/T);

		return D/X_c.dispersivity_c; /*ndim*/
	}

};
