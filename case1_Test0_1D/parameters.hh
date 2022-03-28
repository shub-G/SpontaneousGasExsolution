template<typename PTree>
class Parameters
{
private:
	const PTree& ptree;
	MeshParameters<PTree> mesh;
	const static int dim = MeshParameters<PTree>::dimension;

	double Sg_max_fgp = ptree.get("free_gas_pocket.Sg_max",(double)0.0);

	int numLayers;
	int numProps;
	std::vector<double> z;
	std::vector<std::vector<double> > prop;

	double nxch4;

	double wave_t0;
	double wave_H;
	double wave_A;
	double wave_T;

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
			z[n_layer] = ptree.get(name+".z",(double)0.);
		}

		numProps = ptree.get("sediment.number_of_properties",(int)9);
		prop = std::vector<std::vector<double> > (numLayers,std::vector<double>(numProps, 0.));
		for(int n_layer=0; n_layer<numLayers; n_layer++ ){
			std::string name = "sediment.layer"+std::to_string(n_layer);
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

		//initial dissolved methane fraction
		nxch4 = ptree.get("initial.nxch4",(double)1.);

		//wave
		wave_t0 = ptree.get("water_column.wave_start_time",(double)0.); // time at which storm starts
		wave_H = ptree.get("water_column.average_height",(double)50.);   // average water column height
		wave_A = ptree.get("water_column.wave_amplitude",(double)10.);	 // wave amplitude
		wave_T = ptree.get("water_column.wave_period",(double)10.);  	 // wave period -> time period of 1 wave

		//reference state
		ref_salinity = ptree.get("reference_state.salinity",(double)0.);
		ref_saltconcentration = ref_salinity * (18.0/58.4); /*MolarMass_H2O/MolarMass_salt*/
		ref_temperature = 273.15+ptree.get("reference_state.temperature",(double)0.);

		//gravity
		gravity_flag = ptree.get("gravity.flag",(bool)true);
		g_magnitude = ptree.get("gravity.magnitude",(double)9.81);
  }

	/**********************************************************************
	 * INPUTS
	 **********
	 * Sg_FGP	: gas saturation in the FGP
	 * z_layers	: vector of top-locations of each sediment layer in the domain [m] (counted from bottom)
	 * layers 	: matrix of the properties of each sediment layer
	 * 			  Properties include: (0)porosity, (1)permeability, (2)pentry, (3)lambda, (4)Swr, (5)Sgr, (6)beta, (7)sediment_compressibility (+dP), (8)sediment_compressibility (-dP)
	 * 			  matrix rows identify the layer number
	 * 			  matrix columns identify the property
	 **********************************************************************/

	//2. Free gas pocket (FGP)

	double Sg_FGP(Dune::FieldVector< double,dim > xglobal /*m*/) const {
		double Sg = Sg_max_fgp;
		return Sg;
	}

	//3. Sediment layer geometry and properties

	std::vector<double> z_layers() const {
		return z; /*m*/
	}

	std::vector< std::vector<double> > layer_properties() const {
		return prop; /*returned WITH dimensions*/
	}

	// 4. Pw at sea floor
	double Pw_SF(double time /*s*/) const {
		double nu = 1./wave_T;  // wave frequency -> no. of waves per sec -> inverse of time period of 1 wave
		double h = wave_H;
		if(time>wave_t0){
			h += wave_A*sin(2.*M_PI*nu*(time-wave_t0)); // water depth at given time
		}
		double rhow = 1000.; // density of water in the water column
		double pwsf = rhow * g()[mesh.dimension-1] * h;
		return pwsf; /*Pa*/
	}

	//5. initial dissolved gas fraction
	double nXCH4()const{
		return nxch4;
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

	/**********************************************************************/

};
