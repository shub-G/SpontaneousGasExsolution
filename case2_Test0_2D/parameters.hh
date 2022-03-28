template<typename PTree>
class Parameters
{
private:
	const PTree& ptree;
	MeshParameters<PTree> mesh;
	const static int dim = MeshParameters<PTree>::dimension;

	double Sg_max_fgp;
	double nXCH4_t0;

	int numLayers;
	int numProps;
	std::vector<double> z;
	std::vector<std::vector<double> > prop;
	std::vector<bool> isRandom;

	double wave_H;
	double wave_storm_A;
	double wave_storm_T;
	double wave_storm_N;
	double wave_tide_A;
	double wave_tide_T;

	double ref_salinity;
	double ref_saltconcentration;
	double ref_temperature;

	bool gravity_flag;
	double g_magnitude;

public:

  //! constructor
  Parameters (const PTree& ptree_)
  :ptree(ptree_),
   mesh(ptree_)
  {
		numLayers = ptree.get("sediment.number_of_layers",(int)1);
		z = std::vector<double> (numLayers,0.);
		for(int n_layer=0; n_layer<numLayers; n_layer++ ){
			std::string name = "sediment.layer"+std::to_string(n_layer);
			z[n_layer] = ptree.get(name+".z",(double)100.);
		}

		numProps = ptree.get("sediment.number_of_properties",(int)1);
		prop = std::vector<std::vector<double> > (numLayers,std::vector<double>(numProps, 0.));
		isRandom = std::vector<bool> (numLayers,false);
		for(int n_layer=0; n_layer<numLayers; n_layer++ ){
			std::string name = "sediment.layer"+std::to_string(n_layer);
			isRandom[n_layer] = ptree.get(name+".random",(bool)false);
			prop[n_layer][0] = ptree.get(name+".por",	(double)0.5);
			prop[n_layer][1] = ptree.get(name+".K",		(double)1.e-12);
			prop[n_layer][2] = ptree.get(name+".pentry",(double)1.e4);
			prop[n_layer][3] = ptree.get(name+".lambda",(double)1.2);
			prop[n_layer][4] = ptree.get(name+".swr",	(double)0.);
			prop[n_layer][5] = ptree.get(name+".sgr",	(double)0.);
			prop[n_layer][6] = ptree.get(name+".beta",	(double)1.);
			prop[n_layer][7] = ptree.get(name+".compressibility_positivedP",(double)0.);
			prop[n_layer][8] = ptree.get(name+".compressibility_negativedP",(double)0.);
		}

		//wave
		wave_H = ptree.get("water_column.average_height",(double)50.);   // average water column height
		wave_storm_A = ptree.get("water_column.storm.wave_amplitude",(double)10.);	 // wave amplitude
		wave_storm_T = ptree.get("water_column.storm.wave_period",(double)10.);  	 // wave period -> time period of 1 wave
		wave_storm_N = ptree.get("water_column.storm.wave_number",(double)10.);  	 // wave number -> no. of waves per unit length
		wave_tide_A  = ptree.get("water_column.tide.wave_amplitude",(double)5.);	 // wave amplitude
		wave_tide_T  = ptree.get("water_column.tide.wave_period",(double)12.*3600.); // wave period -> time period of 1 wave

		//reference state
		ref_salinity = ptree.get("reference_state.salinity",(double)0.);
		ref_saltconcentration = ref_salinity * (18.0/58.4); /*MolarMass_H2O/MolarMass_salt*/
		ref_temperature = 273.15+ptree.get("reference_state.temperature",(double)0.);

		//initial free gas saturation in fgp
		Sg_max_fgp = ptree.get("free_gas_pocket.Sg_max",(double)0.0);

		//initial dissolved gas fatcor
		nXCH4_t0 = ptree.get("initial.nXCH4",(double)1.0);

		//gravity
		gravity_flag = ptree.get("gravity.flag",(bool)true);
		g_magnitude = ptree.get("gravity.magnitude",(double)9.81);
  }

	/**********************************************************************
	 * INPUTS
	 **********
	 * z_domain : height of the computational domain [m]
	 * z_cells	: no. of cells along Z-axis
	 * z_FGP 	: location of the top of free gas pocket (FGP) [m]
	 * th_FGP	: thickness of the FGP [m]
	 * Sg_FGP	: gas saturation in the FGP
	 * z_layers	: vector of top-locations of each sediment layer in the domain [m]
	 * layers 	: matrix of the properties of each sediment layer
	 * 			  Properties include: (0)porosity, (1)permeability, (2)pentry, (3)lambda, (4)Swr, (5)Sgr, (6)beta
	 * 			  matrix rows identify the layer number
	 * 			  matrix columns identify the property
	 **********************************************************************/

	//1. Initial percentage of dissolved gas (0-1)

	double nXCH4(Dune::FieldVector< double,dim > xglobal /*m*/) const {
		double nxch4 = nXCH4_t0;
		return nxch4;
	}

	//2. Free gas pocket (FGP)

	double Sg_FGP(Dune::FieldVector< double,dim > xglobal /*m*/) const {
		double Sg = Sg_max_fgp;
		return Sg;
	}

	//3. Sediment layer geometry and properties

	std::vector<double> z_layers() const {
		return z;
	}

	std::vector< std::vector<double> > layer_properties() const {
		return prop;
	}

	std::vector<bool> layer_property_distribution_flag() const {
		return isRandom;
	}


	// 4. Pw at sea floor
	double Pw_SF(Dune::FieldVector< double,dim > xglobal /*m*/, double time) const {
		double nu_storm = 1./wave_storm_T;  // storm frequency -> no. of waves per sec -> inverse of time period of 1 wave
		double nu_tide  = 1./wave_tide_T ;  // tide frequency  -> no. of waves per sec -> inverse of time period of 1 wave
		double h = wave_H;
		if(time>0){
			h += wave_storm_A*sin(2.*M_PI*(wave_storm_N*xglobal[0]-nu_storm*time)); // water depth at given time
			h += wave_tide_A*sin(2.*M_PI*nu_tide*time); // water depth at given time
		}
		double rhow = 1000.; // density of water in the water column
		double pwsf = rhow * g()[mesh.dimension-1] * h;
		return pwsf;
	}


	/**********************************************************************/
	/* REFERENCE STATE */
	double ReferenceSalinity() const {
		return ref_salinity; /*kg/kg*/
	}
	double ReferenceSaltConcentration() const {
		return ref_saltconcentration;
	}
	double ReferenceTemperature() const {
		return ref_temperature; /*K*/
	}

    /**********************************************************************/
    /* SEDIMENT BURIAL RATE */
	Dune::FieldVector<double,dim>
	SedimentationVelocity ( double time, double dt ) const {

		Dune::FieldVector<double,dim> vs( 0. );
		vs[dim-1] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}

	Dune::FieldVector<double,dim>
	SedimentVelocity ( double time, double dt ) const {

		Dune::FieldVector<double,dim> vs( 0. );
		vs[dim-1] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}

	/* GRAVITY VECTOR */
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = g_magnitude;
		gravity[dim-1] = g;
		gravity[0] = 0.;
		return gravity; /*N/kg*/
	}

	// COMPACTION
	double CompactionFunction( Dune::FieldVector<double,dim> globalpos ) const {

		double z = globalpos[dim-1]; /*m*/
		double beta = 1./3000.;
		double compaction_factor = std::exp( -beta*(mesh.Zmax-z) );

		return compaction_factor;
	}

	/**********************************************************************/

};
