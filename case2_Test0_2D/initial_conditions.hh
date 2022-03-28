template<typename GV, typename Properties>
class ProblemInitialConditions
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;

public:

	  //! construct from grid view
	  ProblemInitialConditions (const GV& gv_, const Properties& property_)
	  : gv( gv_ ),
		property(property_)
	  {}

	  /* Initial Conditions */
	  std::vector< double >
	  evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
			  	const Dune::FieldVector<double,dim>& xlocal) const {

		  auto xglobal = element.geometry().global(xlocal);
		  
		  std::vector< double > icvalue(Indices::numOfPVs,0.);

		  /******************************************************************************/
		  // SATURATIONS
		  double Sg = 0.0;
		  if(property.mesh.isInitialFGP(xglobal)){
			  Sg = property.parameter.Sg_FGP(xglobal);
		  }
		  double Sw = 1.-Sg;

		  /******************************************************************************/
		  // PRESSURES
		  double porosity = property.soil.SedimentPorosity(element,xlocal);
		  double Pc = property.hydraulicProperty.CapillaryPressure(element,xlocal,Sw,porosity)
				    * property.characteristicValue.P_c; /*Pa*/
		  double Pw_top = property.parameter.Pw_SF(xglobal,0.);
		  double Pw = Pw_top + 1000.*property.parameter.g()[dim-1]*(property.mesh.get_Zmax()-xglobal[dim-1])*property.characteristicValue.x_c; 
		  double Pg = Pw + Pc; /*Pa*/
		  
		  /******************************************************************************/
		  // MOLE FRACTIONS
		  auto S = property.parameter.ReferenceSalinity();
		  auto T = property.parameter.ReferenceTemperature();
		  auto zCH4 = property.eos.EvaluateCompressibilityFactor( T,Pg );
		  auto Xf = property.mixture.EquilibriumMoleFractions( T, Pg, S, zCH4);

		  auto xch4 = property.parameter.nXCH4(xglobal)*Xf[Indices::concId_XCH4];
		  auto yh2o = 0.;
		  if(property.mesh.isInitialFGP(xglobal)){
			xch4 = Xf[Indices::concId_XCH4];
		  	yh2o = Xf[Indices::concId_YH2O];	
		  }

		  /******************************************************************************/
		  icvalue[Indices::PVId_Sg] = Sg ;
		  icvalue[Indices::PVId_Pw] = Pw/property.characteristicValue.P_c ;
		  icvalue[Indices::PVId_XCH4 ] = xch4;
		  icvalue[Indices::PVId_YH2O ] = yh2o;
		  /******************************************************************************/

		  return icvalue; /*ndim*/
	  }

	  //! get a reference to the grid view
	  inline const GV& getGridView () {return gv;}
};
