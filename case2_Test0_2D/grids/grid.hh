template<typename PTree>
class MeshParameters{
private:
	const PTree& ptree;
	constexpr static double eps = 1.0e-6;
	CharacteristicValues Xc;

	double Zmax;
	double Xmax;
#ifdef USE_YASP
	int nZ;
	int nX;
#endif

	double z_FGP;
	double dz_FGP;
	double x_FGP;
	double dx_FGP;

public:
	
	//! constructor
	MeshParameters (const PTree& ptree_)
	:ptree(ptree_)
	{
		std::cout<<"Construction MeshParameters...";

#ifdef USE_YASP
		Zmax = ptree.get("grid.yasp.LZ",(double)1.)/Xc.x_c; //m
		Xmax = ptree.get("grid.yasp.LX",(double)1.)/Xc.x_c; //m
		nZ = ptree.get("grid.yasp.NZ",(int)10)	;
		nX = ptree.get("grid.yasp.NX",(int)10)	;
#elif defined(USE_UG)
		Zmax = ptree.get("grid.ug.LZ",(double)1.)/Xc.x_c; //m
		Xmax = ptree.get("grid.ug.LX",(double)1.)/Xc.x_c; //m
#endif
		z_FGP  = ptree.get("free_gas_pocket.z", (double)0.);
		dz_FGP = ptree.get("free_gas_pocket.dz",(double)0.);
		x_FGP  = ptree.get("free_gas_pocket.x", (double)0.);
		dx_FGP = ptree.get("free_gas_pocket.dx",(double)0.);

		std::cout<< "Done." << std::endl;
	}

	/*
	 * 2D -> X and Z
	 */
	const static int dimension = 2;

	double get_Zmax() const {
		return Zmax;
	}

	double get_Xmax() const {
		return Xmax;
	}

#ifdef USE_YASP
	int get_Zcells() const {
		return nZ;
	}

	int get_Xcells() const {
		return nX;
	}
#endif

	bool isTop( Dune::FieldVector<double,dimension> globalPos ) const {
		if( globalPos[1] > Zmax - eps ){
			return true;
		}
		else return false;
	}

	bool isBottom( Dune::FieldVector<double,dimension> globalPos ) const {
		if( globalPos[1] < 0. + eps ){
			return true;
		}
		else return false;
	}
	
	bool isInitialFGP( Dune::FieldVector<double,dimension> globalPos ) const {
		if( 	globalPos[1]<z_FGP and globalPos[1]>(z_FGP-dz_FGP)
			and globalPos[0]<x_FGP and globalPos[0]>(x_FGP-dx_FGP)){
			return true;
		}
		else return false;
	}

};

