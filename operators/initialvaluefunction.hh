/** \brief A function for initial values of Porosity
 */
template<typename GV, typename Properties>
class Por_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Por_Initial<GV,Properties> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;

public:

  //! construct from grid view
  Por_Initial (const GV& gv_, const Properties& property_)
  : gv( gv_ ),
	property(property_)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
                        const Dune::FieldVector<double,dim>& xlocal,
                        double& y) const
  {
    y /*ndim*/ = property.soil.SedimentPorosity( element,xlocal ) ;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Pw
 */
template<typename GV, typename Properties>
class Pw_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Pw_Initial<GV,Properties> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties> icvalue;

public:

  //! construct from grid view
  Pw_Initial (const GV& gv_, const Properties& property_)
  : gv( gv_ ),
	property(property_),
	icvalue(gv_,property_)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
                        const Dune::FieldVector<double,dim>& xlocal,
                        double& y) const
  {
    y /*ndim*/ = icvalue.evaluate(element,xlocal)[Indices::PVId_Pw] ;
//    std::cout<< "Pw_boundary = " << y << std::endl;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Sg
 */
template<typename GV, typename Properties>
class Sg_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Sg_Initial<GV,Properties> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties> icvalue;

public:

  //! construct from grid view
  Sg_Initial (const GV& gv_, const Properties& property_)
  : gv( gv_ ),
	property(property_),
	icvalue(gv_,property_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
	y = icvalue.evaluate(element,xlocal)[Indices::PVId_Sg] ;
//	std::cout<< "Sg_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/** \brief A function for initial values of XCH4
 */
template<typename GV, typename Properties>
class XCH4_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, XCH4_Initial<GV,Properties> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties> icvalue;

public:

  //! construct from grid view
  XCH4_Initial (const GV& gv_, const Properties& property_)
  : gv( gv_ ),
	property(property_),
    icvalue(gv_,property_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
	y = icvalue.evaluate(element,xlocal)[Indices::PVId_XCH4] ;
//	std::cout<< "XCH4_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of YH2O
 */
template<typename GV, typename Properties>
class YH2O_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, YH2O_Initial<GV,Properties> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties> icvalue;

public:

  //! construct from grid view
  YH2O_Initial (const GV& gv_, const Properties& property_)
  : gv( gv_ ),
	property(property_),
	icvalue(gv_,property_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
	y = icvalue.evaluate(element,xlocal)[Indices::PVId_YH2O] ;
//	std::cout<< "YH2O_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};
