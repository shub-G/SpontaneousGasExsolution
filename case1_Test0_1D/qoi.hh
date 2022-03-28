/*****************************************************
 * Quantities Of Interest
 *****************************************************
 *
 *****************************************************/

template< class GV, class PARAMS, class GFSOP, class UOP >
class QoI{

public:
	typedef Dune::PDELab::LocalFunctionSpace<GFSOP> LFS;

    typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
	typedef typename UOP::template LocalView<LFSCache> VectorView;

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;
    static const int dim = GV::dimension;

    QoI(std::string file_name_,
    	const GV& gv_,
    	const PARAMS& param_,
		GFSOP 	gfs_,
		UOP		*u_,
		double	*time_	,
		double	*dt_,
		int		*opnum_ )
    :file_name(file_name_),
	 gv(gv_),
	 param( param_ ),
	 gfs(gfs_),
	 u(u_),
	 time( time_ ),
	 dt( dt_ ),
	 opnum(opnum_)
    {
    	x = Xc.x_c;
    }

	virtual ~QoI(){}


	void generate_csv_output(){

		std::string name = file_name;
		name += ".";
		name += std::to_string((*opnum));
		name += ".csv";
		std::fstream result;
		result.open(name, std::fstream::out | std::fstream::trunc );
		if (!result.is_open()){
			std::cout << "Could not open file:" << name << std::endl;
			exit(0);
		}

		result << "z[m]" << " ";
		result << "Pg[Pa]" << " ";
		result << "Pw[Pa]" << " ";
		result << "Pc[Pa]" << " ";
		result << "Sg[-]" << " ";
		result << "Sw[-]" << " ";
		result << "XCH4[mol/m^3]" << " ";
		result << "XH2O[mol/m^3]" << " ";
		result << "YCH4[mol/m^3]" << " ";
		result << "YH2O[mol/m^3]" << " ";
		result << "rhog[kg/m^3]" << " ";
		result << "rhow[kg/m^3]" << " ";
		result << "K[m^2]" << " ";
		result << "krg[-]" << " ";
		result << "krw[-]" << " ";
		result << "mug[Pa.s]" << " ";
		result << "muw[Pa.s]" << " ";
		result << "zCH4[-]" << " ";
		result << "porosity[-]" << " ";
		result << "DCH4[m^2/s]" << " ";
		result << "DH2O[m^2/s]" << " ";
		result << "HCH4[1/Pa]" << " ";
		result << "pwsat[Pa]" << " ";
		result << "tortuosity[-]" << " ";
		result << std::endl;
		result.flush();

		LFS lfs(gfs);
		LFSCache lfs_cache(lfs);
		VectorView u_view( *u );

		// Iterate over each element
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			lfs.bind(*self);
			lfs_cache.update();
			u_view.bind(lfs_cache);
		    std::vector<double> ul(lfs.size());
		    u_view.read( ul );

			// Get Element geometry
			typedef typename LeafIterator::Entity::Geometry ElementGeometry;
			const auto& geo = (*self).geometry();
	        auto cell_center_local 	= geo.local(geo.center());
	        auto cell_center_global = geo.center();


	        // write values to file:
			result << (cell_center_global[dim-1]-param.mesh.get_Zmax()) * x << " ";
			result << ul[lfs.child(param.index.SVId_Pg	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_Pw	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_Pc	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_Sg	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_Sw	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_XCH4 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_XH2O ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_YCH4 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_YH2O ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_rhog ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_rhow ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_K	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_krg	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_krw	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_mug	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_muw	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_zCH4 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_por	 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_DCH4 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_DH2O ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_HCH4 ).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_Pwsat).localIndex(0)] << " ";
			result << ul[lfs.child(param.index.SVId_tau	 ).localIndex(0)] << " ";
			result << std::endl;
			result.flush();

		}//END: loop over volumes

		result.close();

	}//END: generate_csv_output()

private:
	std::string file_name;
    const GV& gv;
	const PARAMS& param;
	GFSOP gfs;
	UOP *u;
	double *time;
	double *dt;
	int *opnum;
	CharacteristicValues Xc;
	double x;
};
