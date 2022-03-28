template<typename GV, typename Properties>
class ProblemInitialConditions
{
private:
	  const GV& gv;
	  const Properties& property;
	  double x_c;
	  double P_c;
	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;

public:

	  //! construct from grid view
	  ProblemInitialConditions (const GV& gv_, const Properties& property_)
	  : gv( gv_ ),
		property(property_)
	  {
		  x_c =property.characteristicValue.x_c;
		  P_c =property.characteristicValue.P_c;
	  }

	  /* Initial Conditions */
	  std::vector< double >
	  evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
			  	const Dune::FieldVector<double,dim>& xlocal) const {

		  auto xglobal/*ndim*/ = element.geometry().global(xlocal);
		  
		  std::vector< double > icvalue(Indices::numOfPVs,0.);

		  /******************************************************************************/
		  // SATURATIONS
		  double Sg = 0.0;
		  if(property.mesh.isInitialFGP(xglobal)){
			  Sg = property.parameter.Sg_FGP(xglobal*x_c);
		  }
		  double Sw = 1.-Sg;

		  /******************************************************************************/
		  // PRESSURES
		  double porosity = property.soil.SedimentPorosity(element,xlocal);
		  double Pc = property.hydraulicProperty.CapillaryPressure(element,xlocal,Sw,porosity) * P_c; /*Pa*/
		  double Pw_top = property.parameter.Pw_SF(0.); /*Pa*/
		  double Pw = Pw_top + 1000.*property.parameter.g()[dim-1]*(property.mesh.get_Zmax()-xglobal[dim-1])*x_c;
		  double Pg = Pw + Pc; /*Pa*/
		  
		  /******************************************************************************/
		  // MOLE FRACTIONS
		  auto S = property.parameter.ReferenceSalinity();/*kg/kg*/
		  auto T = property.parameter.ReferenceTemperature();/*K*/
		  auto zCH4 = property.eos.EvaluateCompressibilityFactor( T/*K*/,Pg/*Pa*/ );
		  auto Xf = property.mixture.EquilibriumMoleFractions( T/*K*/, Pg/*Pa*/, S/*kg/kg*/, zCH4);

		  auto xch4 = property.parameter.nXCH4() * Xf[Indices::concId_XCH4];
		  if( property.mesh.isInitialFGP(xglobal/*ndim*/)){
			  xch4 = Xf[Indices::concId_XCH4];
		  }
		  auto yh2o = Xf[Indices::concId_YH2O];

		  /******************************************************************************/
		  icvalue[Indices::PVId_Sg] = Sg ;
		  icvalue[Indices::PVId_Pw] = Pw/P_c ;
		  icvalue[Indices::PVId_XCH4 ] = xch4;
		  icvalue[Indices::PVId_YH2O ] = yh2o;
		  /******************************************************************************/

		  return icvalue; /*ndim*/
	  }

	  //! get a reference to the grid view
	  inline const GV& getGridView () {return gv;}
};
