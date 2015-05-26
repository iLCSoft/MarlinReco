#include "MCOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using IMPL::MCParticleImpl;
namespace TTbarAnalysis
{
	MCOperator:: MCOperator (LCCollection * col, LCCollection * rel)
	{
		myCollection = col;
		myRelCollection = rel;
		myPrimaryVertex[0] = 0.0;
		myPrimaryVertex[1] = 0.0;
		myPrimaryVertex[2] = 0.0;
		myAngleCut = 0.5; // \pi/4
	}
	DecayChain * MCOperator::RefineDecayChain(DecayChain * initial, vector<PDGTYPE> typeOfProducts)
	{
		if (!initial || initial->GetSize() < 2) 
		{
			return NULL;
		}
		int size = initial->GetSize();
		DecayChain * result = new DecayChain(initial->GetName(), initial->GetParentPDG());
		for (int i = 1; i < size; i++) 
		{
			MCParticle * particle = initial->Get(i);
			vector< MCParticle * > parents = particle->getParents();
			MCParticle * parent = NULL;
			if (parents.size() == 1 && CheckParticle(parents[0],typeOfProducts[i-1])) 
			{
				parent = parents[0];
			}
			else 
			{
				parent = initial->Get(i-1);
			}
			if (parent) 
			{
				result->Add(parent);
			}
		}

		result->Add(initial->Get(size-1));
		return result;
	}
	DecayChain * MCOperator::Construct(string name, int pdg, vector<PDGTYPE> typeOfProducts)
	{
		DecayChain * result = new DecayChain(name, pdg);
		int number = myCollection->getNumberOfElements();
		bool finished = false;
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) ) ;
			if (i > 100 && result->GetSize() == 0) 
			{
				streamlog_out(DEBUG) << "WARNING: EVENT RUN OUT OF TRIALS!\n";
				break;
			}
			if (particle->getPDG() == pdg) 
			{
				if (CheckForUnification(particle, pdg)) 
				{
		//			streamlog_out(DEBUG) << "WARNING: found a quark-antiquark merging!\n";
					continue;
				}
				for (unsigned int j = 0; j < typeOfProducts.size(); j++) 
				{
					PDGTYPE currentPDG = typeOfProducts[j];
					MCParticle * found = FindYoungestChild(particle, currentPDG);
					if (found) 
					{
						if (CheckForColorString(found, pdg)) 
						{
							MCParticle * service = found->getParents()[0];
							MCParticle * consistent = GetConsistentDaughter(particle, service, currentPDG);
							if (consistent) 
							{
								result->Add(consistent);
								found = consistent;
								finished = true;
							}
						}
						else 
						{
							result->Add(found);
							finished = true;
						}
					}
					else 
					{
						break;
					}
					particle =(finished)? found : NULL;
				}
			}
			if (finished)
			{
				break;
			}
		}
		return result;
	}
	bool MCOperator::CheckForUnification(MCParticle * particle, int pdg)
	{
		vector< MCParticle * > daughters = particle->getDaughters();
		if (daughters.size() == 0) 
		{
			return false;
		}
		if (daughters.size() == 1 && daughters[0]->getPDG() == pdg) 
		{
			return CheckForUnification(daughters[0], pdg);
		}
		if (daughters.size() == 1 && (daughters[0]->getPDG() == 92 || daughters[0]->getPDG() == 94)) 
		{
			vector< MCParticle * > grandDaughters = daughters[0]->getDaughters();
			bool foundParticle = false;
			bool foundAntiparticle = false;
			for (unsigned int i = 0; i < grandDaughters.size(); i++) 
			{
				if (grandDaughters[i]->getPDG() == pdg) 
				{
					foundParticle = true;
				}
				if (grandDaughters[i]->getPDG() == -pdg) 
				{
					foundAntiparticle = true;
				}
			}
			return foundParticle && foundAntiparticle;
		}
		return false;

	}
	bool MCOperator::CheckParticle(MCParticle * particle, PDGTYPE type)
	{
		bool result = false;
		vector<int> pdg =  ConstantStorage::GET_PDG(type);
		for (unsigned int i = 0; i < pdg.size(); i++) 
		{
			if (abs(particle->getPDG()) == pdg[i]) 
			{
				result = true;
				break;
			}
		}
		return result;
	}
	bool MCOperator::CheckForColorString(EVENT::MCParticle * daughter, int pdgOfParent)
	{
		vector< MCParticle * > parents = daughter->getParents();
		if (parents.size() == 1) 
		{
			MCParticle * parent = parents[0];
			if (parent->getPDG() == 92) 
			{
				vector< MCParticle * > grandParents = parent->getParents();
				int count = 0;
				for (unsigned int i = 0; i < grandParents.size(); i++) 
				{
					if (abs(grandParents[i]->getPDG()) == abs(pdgOfParent)) 
					{
						count++;
					}
					if (count > 1) 
					{
						streamlog_out(DEBUG) << "United color string found for PDG " << pdgOfParent <<"!\n";
						return true;
					}
				}
			}
		}
		return false;
	}

	MCParticle * MCOperator::GetConsistentDaughter(MCParticle * parent, MCParticle * service, PDGTYPE type)
	{
		vector< MCParticle * > daughters = SelectDaughtersOfType(service, type); //service->getDaughters();
		streamlog_out(DEBUG) << "Parent PDG: " << parent->getPDG() << ";\n";
		int size = daughters.size();
		vector<float> angle;
		if (size == 1) 
		{
			 //streamlog_out(DEBUG) << "We have only 1 meson of type " << type << '\n';
		}
		float dE = 10.0;
		for (int i = 0; i < size; i++) 
		{
			MCParticle * daughter = daughters[i];
			//streamlog_out(DEBUG) << "Daughter PDG: " << daughter->getPDG() << ";\n";
			if (MathOperator::getModule(daughter->getMomentum()) > MathOperator::getModule(parent->getMomentum())+dE)// || 
			{
				streamlog_out(DEBUG) << "Discarded by energy!\n";
				//continue;
			}
			if (daughter->getCharge()*parent->getCharge() > 0.0 
				&& (type == BOTTOM_MESONS || type ==  CHARMED_MESONS || type == STRANGE_MESONS))
				//&& type != BOTTOM_HADRONS 
				//&& type != BOTTOM_BARYONS 
				//&& type != CHARMED_HADRONS 
				//&& type != CHARMED_BARYONS)
			{
				//streamlog_out(DEBUG) << "Same sign of charge!\n";
				return daughter;
			}
			angle.push_back(MathOperator::getAngle(daughter->getMomentum(), parent->getMomentum()));
		}
		float minAngle = myAngleCut;
		int winner = -1;
		//streamlog_out(DEBUG) << "Checking angles...\n";
		for (unsigned int i = 0; i < angle.size(); i++) 
		{
			//streamlog_out(DEBUG) << "Angle " << i << ": " << angle[i] << '\n';
			if (angle[i] < minAngle)// && 
			   //daughters[i]->getEnergy() < parent->getEnergy()+dE &&
			   //(daughters[i]->getCharge()*parent->getCharge() > -0.0001 && ( type == BOTTOM_MESONS || type == CHARMED_MESONS || type == STRANGE_MESONS ))) 
			{
				minAngle = angle[i];
				winner = i;
			}
		}
		if (winner > -1) 
		{
			//streamlog_out(DEBUG) << "Choosing angle " << winner <<  "\n";
			if (daughters[winner]->getCharge()*parent->getCharge() < -0.0001  && type != BOTTOM_HADRONS && type != BOTTOM_BARYONS) 
			{
				streamlog_out(DEBUG) << "INFO: Charge is wrong!\n";
				
			}
			return daughters[winner];
		}
		//streamlog_out(DEBUG) << "Angle is wrong!\n";
		return NULL;
	
	}
	bool MCOperator::IsReconstructed(MCParticle * particle)
	{
		LCRelationNavigator navigator(myRelCollection);
		vector< LCObject * > obj = navigator.getRelatedFromObjects(particle);
		int nvtx = obj.size();
		if (nvtx < 1) 
		{
			streamlog_out(DEBUG) << "INFO: Particle not reconstructed\n";
			return false;
		}
		const vector< float > weights = navigator.getRelatedFromWeights(particle);
		if (nvtx > 0) 
		{
			int winner = -1;
			float maxweight = 0.0;
			for (unsigned int j = 0; j < obj.size(); j++) 
			{
				ReconstructedParticle * reco = dynamic_cast< ReconstructedParticle * >(obj[j]);
				if (weights[j] > maxweight && std::abs(reco->getCharge()) > 0.09)
				{
					maxweight = weights[j];
					winner = j;
				}
			}
			return winner > -1;
		}
		return nvtx > 0 && weights[0] > 0.5;

	}
	PDGTYPE MCOperator::GetParticleType(MCParticle * particle)
	{
		PDGTYPE result = EXCEPTIONAL_PDGTYPE;
		int count = _max_MESONS;
		for (int i = 1; i < count; i++) 
		{
			PDGTYPE candidate = static_cast<PDGTYPE>(i);
			vector<int> pdg =  ConstantStorage::GET_PDG(candidate);
			for (unsigned int j = 0; j < pdg.size(); j++) 
			{
				if (abs(particle->getPDG()) == pdg[j]) 
				{
					result = candidate;
					break;
				}
			}
		}
		return result;

	}
	MCParticle * MCOperator::FindExceptionalChild(MCParticle * parent, PDGTYPE parentType)
	{
		if (!parent) 
		{
			//streamlog_out(DEBUG) << "ERROR: Input particle is NULL!\n";
			return NULL;
		}
		vector< MCParticle * > daughters = parent->getDaughters();
		if (daughters.size() < 1) 
		{
			//streamlog_out(DEBUG)<<"ERROR: Found 0 daughters!\n";
			return NULL;
		}
		if (parentType == EXCEPTIONAL_PDGTYPE) 
		{
			parentType = GetParticleType(parent);
			if (parentType == EXCEPTIONAL_PDGTYPE) 
			{
				streamlog_out(DEBUG) << "ERROR: Probably, parent is not meson, to be continue....\n";
				return NULL;
			}
		}
		bool found = false;
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (CheckParticle(daughter,parentType)) 
			{
				//streamlog_out(DEBUG) << "Found PDG to exclude: " << daughter->getPDG() << "; Going recursion\n";
				found = true;
			}
		}
		if (!found) 
		{
			//streamlog_out(DEBUG) << "Found state without initial mesons!\n";
			return daughters[0];
		}
		//streamlog_out(DEBUG)<<"Not found in I gen. Checking next gen.\n";
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (daughter->getEnergy() < 1.0 && daughter->getMass() < 300.0) 
			{
				//streamlog_out(DEBUG)<<"Daughter does not pass the selection cuts\n";
				continue;
			}
			return FindExceptionalChild(daughter, parentType);
		}
		return NULL;
	}
	MCParticle * MCOperator::FindYoungestChild(MCParticle * parent, PDGTYPE type)
	{
		if (!parent) 
		{
			//streamlog_out(DEBUG) << "ERROR: Input particle is NULL!\n";
			return NULL;
		}
		if (type == EXCEPTIONAL_PDGTYPE) 
		{
			return FindExceptionalChild(parent, EXCEPTIONAL_PDGTYPE);
		}
		vector< MCParticle * > daughters = parent->getDaughters();
		if (daughters.size() < 1) 
		{
			//streamlog_out(DEBUG)<<"ERROR: Found 0 daughters!\n";
			return NULL;
		}
		//streamlog_out(DEBUG)<<"Checking " << daughters.size() <<" I gen of daughters. \n";
		MCParticle * result = NULL;
		int count = 0;
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (CheckParticle(daughter,type)) 
			{
				//streamlog_out(DEBUG) << "Found PDG: " << daughter->getPDG() << ";\n";
				count++;
				result = daughter;
			}
		}
		if (result) 
		{
			if (count > 1 && abs(parent->getPDG()) !=92  && (type == CHARMED_MESONS || type == CHARMED_HADRONS)) 
			{
				streamlog_out(DEBUG)<<"FATAL: Daughter multiplicity found for meson type " << type << " and parent " << parent->getPDG() <<"!\n";
				return NULL;
			}
			return result;
		}
		//streamlog_out(DEBUG)<<"Not found in I gen. Checking next gen.\n";
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (daughter->getEnergy() < 1.0 && daughter->getMass() < 300.0) 
			{
				//streamlog_out(DEBUG)<<"Daughter does not pass the selection cuts\n";
				//continue;
			}
			return FindYoungestChild(daughter, type);
		}
		return NULL;
	}
	EVENT::MCParticle * MCOperator::cureDoubleCharmDecay(std::vector< EVENT::MCParticle * > & selected)
	{
		MCParticleImpl * newparticle = new MCParticleImpl();
		newparticle->setMomentum(selected[0]->getMomentum());
		newparticle->setMass(selected[0]->getMass());
		newparticle->setCharge(selected[0]->getCharge());
		newparticle->setPDG(selected[0]->getPDG() + selected[1]->getPDG());
		return newparticle;
	}
	vector< MCParticle * > MCOperator::SelectDaughtersOfType(MCParticle * parent, PDGTYPE type)
	{
		vector< MCParticle * > daughters = parent->getDaughters();
		vector< MCParticle * > result;
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			if (CheckParticle(daughters[i], type)) 
			{
				result.push_back(daughters[i]);
			}
		}
		return result;
	}
	
	vector< MCParticle * > MCOperator::GetPairParticles(int pdg)
	{
		pdg = abs(pdg);
		vector< MCParticle * > pair;
		if (pdg < 1) 
		{
			return pair;
		}
		int number = myCollection->getNumberOfElements();
		MCParticle * b = NULL;
		MCParticle * bbar = NULL;
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (particle->getPDG() == pdg)// && countParticle == 0) 
			{
				b = particle;
			}
			if (particle->getPDG() == -pdg)// && countAntiparticle == 0) 
			{
				bbar =  particle;
			}
		}
		if (b) 
		{
			pair.push_back(new MCParticleImpl((const IMPL::MCParticleImpl&)(*b)));
			streamlog_out(DEBUG)<<"INFO: Found b!\n";
		}
		if (bbar) 
		{
			pair.push_back(new MCParticleImpl((const IMPL::MCParticleImpl&)(*bbar)));
			streamlog_out(DEBUG)<<"INFO: Found bbar!\n";
		}
		return pair;
	}
	vector< MCParticle * > MCOperator::GetPairParticles(PDGTYPE type)
	{
		vector< MCParticle * > pair;
		int number = myCollection->getNumberOfElements();
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (CheckParticle(particle, type) && particle->getParents()[0]->getPDG() == 92) 
			{
				pair.push_back(new MCParticleImpl((const IMPL::MCParticleImpl&)(*particle)));
				streamlog_out(DEBUG)<<"Found a suitable particle of type " << particle->getPDG() << '\n';
			}
			if (pair.size() > 1) 
			{
				break;
			}
		}
		return pair;
	}
	/*bool MCOperator::CheckProcessForPair(int pdg)
	{
		pdg = abs(pdg);
		if (pdg < 1) 
		{
			return false;
		}
		int countParticle = 0;
		int countAntiparticle = 0;
		int number = myCollection->getNumberOfElements();
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (particle->getPDG() == pdg) 
			{
				countParticle++;
			}
			if (particle->getPDG() == 0-pdg) 
			{
				countAntiparticle++;
			}
		}
		streamlog_out(DEBUG)<<"INFO: " << countParticle << " quarks and " << countAntiparticle << " antiquarks\n";
		return countParticle && countAntiparticle;
	}//*/
	vector< MCParticle * > * MCOperator::CheckProcessForPair(int pdg)
	{
		pdg = abs(pdg);
		if (pdg < 1) 
		{
		        return false;
		}
		MCParticle * particle = NULL;
		MCParticle * antiparticle = NULL;
		vector< MCParticle * > * result = NULL;
		for (int i = 0; i < 100; i++) 
		{
			MCParticle * candidate = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (candidate->getPDG() == pdg && !particle) 
			{
				particle = candidate;
			}
			if (candidate->getPDG() == -pdg && !antiparticle) 
			{
				antiparticle = candidate;
			}
		}
		if (particle && antiparticle) 
		{
			result = new vector< MCParticle * >();
			result->push_back(particle);
			result->push_back(antiparticle);
		}
		return result;
	}
	vector< MCParticle * > MCOperator::SelectStableCloseDaughters(MCParticle * parent,int excludePDG, bool selectReco, vector<MCParticle *> * misReconstructed)
	{
		vector< MCParticle * > result;
		if (!parent || 
		    CheckParticle(parent, NONTRACKABLE_PARTICLES) ||
		    parent->getPDG() == excludePDG)//(CheckParticle(parent, CHARMED_MESONS) && discardCharmedMesons)) 
		{
			return result;
		}
		vector< MCParticle * > daughters = parent->getDaughters();
		if (daughters.size()<1) 
		{
			return result;
		}
		if ( parent && MathOperator::getDistance(daughters[0]->getVertex(), parent->getVertex()) > 100.0) // Long distance!!!
		{
			streamlog_out(DEBUG)<<"WARNING: Long-lived noninteracting particle " << parent->getPDG() <<"!\n";
			return result;
		}
		bool finished = false;
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters[i];
			if (CheckParticle(daughter, TRACKABLE_PARTICLES))
			{
				//streamlog_out(DEBUG)<<"\tDaughter reached calorimeter!\n";
				if (MathOperator::getModule(daughter->getMomentum()) < 0.1) 
				{
					//streamlog_out(DEBUG)<<"CAUTION: Low momentum particle of " << MathOperator::getModule(daughter->getMomentum()) << " GeV!\n";
				}
				finished = true;
				if (selectReco && IsReconstructed(daughter)) 
				{
					result.push_back(daughter);
				}
				if (!selectReco) 
				{
					result.push_back(daughter);
				}
				if (misReconstructed && !IsReconstructed(daughter)) 
				{
					misReconstructed->push_back(daughter);
				}
			}
			else 
			{
				vector< MCParticle * > grannies = SelectStableCloseDaughters(daughter, excludePDG);
				for (unsigned int j = 0; j < grannies.size(); j++) 
				{
					result.push_back(grannies[j]);
				}
			}
		}
		return result;
	}
	float MCOperator::GetAccuracy(MCParticle * particle, float a, float b)
	{
		float p = MathOperator::getModule(particle->getMomentum());
		float m = particle->getMass();
		float gamma = sqrt( 1 + (p * p) / (m * m));
		vector<float> direction = MathOperator::getDirection(particle->getMomentum());
		vector<float> angles = MathOperator::getAngles(direction);
		float accuracy = sqrt(a*a + b*b /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
		return gamma*accuracy;
	}
	bool MCOperator::CheckCompatibility(const vector< MCParticle * > & daughters, MCParticle * parent, int plusCharge)
	{
		int charge = 0;
		float mass = 0.0;
		for (unsigned int i = 0; i < daughters.size(); i++)
		{
			MCParticle * daughter = daughters[i];
			mass += daughter->getMass();
			charge += daughter->getCharge();
		}
		if (charge + plusCharge != (int)(parent->getCharge())) 
		{
			streamlog_out(DEBUG)<< "CAUTION: Charge is not compatible: charges " << charge << " and " << parent->getCharge() << ", " << plusCharge <<"!\n";
			return false;
		}
		if (mass > parent->getMass()) 
		{
			streamlog_out(DEBUG)<< "CAUTION: Mass is not compatible!\n";
			return false;
		}
		return true;
	}
	vector< MCParticle * > MCOperator::CheckDaughterVisibility(const vector< MCParticle * > & daughters)
	{
		int size = daughters.size();
		vector< MCParticle * > filtered;
		for (int i = 0; i < size; i++) 
		{
			MCParticle * daughter = daughters[i];
			vector< float > direction = MathOperator::getDirection(daughter->getMomentum());
			float offset = MathOperator::getDistanceTo(myPrimaryVertex, direction, daughter->getVertex());
			float pt = MathOperator::getPt(daughter->getMomentum());
			float p = MathOperator::getModule(daughter->getMomentum());
			float eta = MathOperator::getAngles(direction)[1];
			streamlog_out(DEBUG) << "\tPDG: " << daughter->getPDG()
				  << " p: " << p
				  << " pt: " << pt
				  << " teta: " << eta
				  << " offset: " << offset
				  << '\n';
			if (offset > 0.01 && abs(eta) < 4 && p > 0.6) 
			{
				filtered.push_back(daughter);
			}
		}
		return filtered;
	}
} /*  */
