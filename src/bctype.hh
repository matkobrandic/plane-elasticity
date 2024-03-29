#ifndef BCTYPE_HH
#define BCTYPE_HH

#include <iostream>

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>


/*	Klasa koja određuje vrijednost Dirichletovog rubnog uvjeta i 
	njegovo proširenje na čitavu domenu. 
	Template parametri:
		GV = GridView
		RF = Range Field Type (tip kojim su predstavljeni elementi slike funkcije)

		Ovo je rubni uvjet ua skalarnu komponentu.
*/
template<typename GV, typename RF>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, BCExtension<GV,RF> > {
	const GV& gv;
public:
	using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;
	BCExtension (const GV& gv_) : gv(gv_) {}

	// Izračunaj Dirichletovu vrijednost na elementu. Ako točka nije na 
	// Dirichletovoj granici, onda funkcija daje proširenje Dirichletovog rubnog
	// uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno. 
	// e     = element 
	// xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
	// y      = izračunata Dirichletova vrijednost
	inline void evaluate (const typename Traits::ElementType& e,
						  const typename Traits::DomainType& xlocal,
						  typename Traits::RangeType& y) const{
	return;
	}

	inline const GV& getGridView () {return gv;}
};

// Ovo je klasa za skalarnu komponentu
// što znači da pretpostavljamo isti tip uvjeta za sve komponente.
template <typename GV>
class BCTypeParam : public Dune::PDELab::DirichletConstraintsParameters{
	const GV& gv;
public:
	BCTypeParam(const GV& gv_) : gv(gv_) {}

public:
	// provjeravamo na kojem smo rubu domene
	template<typename I>
	bool isLeft(const I & intersection,
				const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
			   ) const{
	auto xg = intersection.geometry().global( coord );
	if(xg[0] < 1E-6)
		return true;
	return false;  
	}
	template<typename I>
	bool isDown(const I & intersection,
				const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
			   ) const{
	auto xg = intersection.geometry().global( coord );
	if(xg[1] < 1E-6)
		return true;
	return false;
	}
	template<typename I>
	bool isUp(const I & intersection,
			  const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
			 ) const{
	auto xg = intersection.geometry().global( coord );
	if(xg[1] > 1 - 1E-6)
		return true;
	return false;
	}
	template<typename I>
	bool isRight(const I & intersection,
				 const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
				) const{
	auto xg = intersection.geometry().global( coord );
	if(xg[0] > 1 - 1E-6)
		return true;
	return false;
	}
	inline const GV& getGridView(){ return gv; }
};
#endif
