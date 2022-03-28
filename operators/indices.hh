class Indices {

public:

	const static int numOfPVs 	= 4;
	const static int PVId_Pw	= 0;
	const static int PVId_Sg	= 1;
	const static int PVId_XCH4	= 2;
	const static int PVId_YH2O	= 3;

	const static int numOfConcs  = 4;
	const static int concId_XCH4 = 0;
	const static int concId_XH2O = 1;
	const static int concId_YCH4 = 2;
	const static int concId_YH2O = 3;

	const static int numOfSVs 	= 23;
	const static int SVId_Pg	= 0;  // gas phase pressure
	const static int SVId_Pw	= 1;  // water phase pressure
	const static int SVId_Pc	= 2;  // capillary pressure
	const static int SVId_Sg	= 3;  // gas saturation
	const static int SVId_Sw	= 4;  // water saturation
	const static int SVId_XCH4	= 5;  // CH4 mole fraction in water
	const static int SVId_XH2O	= 6;  // H2O mole fraction in water
	const static int SVId_YCH4	= 7;  // CH4 mole fraction in gas
	const static int SVId_YH2O	= 8;  // H2O mole fraction in gas
 	const static int SVId_rhow 	= 9;  // water density
 	const static int SVId_rhog 	= 10; // methane gas density
 	const static int SVId_K 	= 11; // absolute permeability
 	const static int SVId_krw 	= 12; // rel. water permeability
 	const static int SVId_krg 	= 13; // rel. gas permeability
 	const static int SVId_muw 	= 14; // dyn. water viscosity
 	const static int SVId_mug 	= 15; // dyn. gas viscosity
	const static int SVId_zCH4	= 16; // methane gas compressibility factor
	const static int SVId_por	= 17; // total porosity
 	const static int SVId_DH2O 	= 18; // H2O binary diffusion coeff. in methane
 	const static int SVId_DCH4 	= 19; // CH4 binary diffusion coeff. in water
 	const static int SVId_Pwsat = 20; // Saturation pressure for water vapor
 	const static int SVId_HCH4 	= 21; // Henry's constant
 	const static int SVId_tau	= 22; // tortuosity

	const static int numOfBCs = 2;
	const static int BC_w = 0;
	const static int BC_g = 1;
	const static int dirichlet 	= 0;
	const static int neumann	= 1;

};
