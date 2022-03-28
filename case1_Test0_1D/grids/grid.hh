template<typename PTree>
class MeshParameters{
private:
	const PTree& ptree;
	constexpr static double eps = 1.0e-6;

	CharacteristicValues Xc;
	double x;

	double z_domain;
	int nz;
	double z_FGP;
	double th_FGP;
	
public:
	
	 /****************************************************************
	 * z_domain : height of the computational domain [m]
	 * z_cells	: no. of cells along Z-axis
	 * z_FGP 	: location of the top of free gas pocket (FGP) [m]
	 * th_FGP	: thickness of the FGP [m]
	 *****************************************************************/

	//! constructor
	MeshParameters (const PTree& ptree_)
	:ptree(ptree_)
	{
		x = Xc.x_c;
		z_domain = ptree.get("grid.yasp.LZ",(double)1.)/x;
		nz = ptree.get("grid.yasp.NZ",(int)10)	;
		z_FGP = ptree.get("free_gas_pocket.z", (double)0.)/x;
		th_FGP = ptree.get("free_gas_pocket.dz",(double)0.)/x;
	}

	const static int dimension = 2;

	double get_Zmax() const {
		return z_domain; /*ndim*/
	}

	int get_nZ() const {
		return nz;
	}

	double get_Xmax() const {
		return z_domain/nz; /*ndim*/
	}

	int get_nX() const {
		return 1; /*ndim*/
	}

	bool isTop( Dune::FieldVector<double,dimension> globalPos /*ndim*/ ) const {
		if( globalPos[1] > z_domain - eps ){
			return true;
		}
		else return false;
	}

	bool isBottom( Dune::FieldVector<double,dimension> globalPos /*ndim*/ ) const {
		if( globalPos[1] < 0. + eps ){
			return true;
		}
		else return false;
	}
	
	bool isInitialFGP( Dune::FieldVector<double,dimension> globalPos /*ndim*/ ) const {
		if( globalPos[1]<z_FGP and globalPos[1]>(z_FGP-th_FGP) ){
			return true;
		}
		else return false;
	}

	bool isAboveInitialFGP( Dune::FieldVector<double,dimension> globalPos /*ndim*/ ) const {
		if( globalPos[1]>z_FGP ){
			return true;
		}
		else return false;
	}

	bool isBelowInitialFGP( Dune::FieldVector<double,dimension> globalPos /*ndim*/ ) const {
		if( globalPos[1]<(z_FGP-th_FGP) ){
			return true;
		}
		else return false;
	}

};
