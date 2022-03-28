template<typename PTree>
class Sources
{
private:
	const PTree& ptree;
	Methane methane;
	Water water;
	CharacteristicValues X_c;
	Parameters<PTree> parameter;

public:
	
	//constructor
	Sources (const PTree& ptree_)
	: ptree(ptree_), parameter(ptree_)
	{}

	double GasGenerationRate( double T/*K*/, double Pg /*Pa*/ ) const {
		
		return 0.; //[kg/m^3/s]
	}
	
	double WaterGenerationRate( double T/*K*/, double Pw /*Pa*/ ) const {
		
		return 0.; //[kg/m^3/s]
	}
};
