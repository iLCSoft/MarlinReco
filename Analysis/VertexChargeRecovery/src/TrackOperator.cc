#include "TrackOperator.hh"

using std::vector;
using std::string;
using std::abs;
using EVENT::Track;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
using IMPL::ReconstructedParticleImpl;
namespace TTbarAnalysis
{
	TrackOperator:: TrackOperator()
	{
	}
	ReconstructedParticle * TrackOperator::ReconstructParticle(EVENT::Track *  track)
	{
		float Bz = 3.5;
		float a = 3.0e-4;
		float omega = track->getOmega();
		float tanl = track->getTanLambda();
		float phi = track->getPhi();
		float pt = a * std::abs(Bz / omega);
		//float p = pt * std::sqrt(1 + tanl * tanl);
		
		float * momentum = new float[3];
		momentum[0] = pt * std::cos(phi);
		momentum[1] = pt * std::sin(phi);
		momentum[2] = pt * tanl;
		float mass = 0.140;
		float p = pt * sqrt( 1 + tanl * tanl );

		float charge = Bz / omega / std::abs(Bz / omega);
		ReconstructedParticleImpl * result = new ReconstructedParticleImpl();
		result->setCharge(charge);
		result->setMomentum(momentum);
		result->addTrack(track);
		result->setMass(mass);
		result->setEnergy(sqrt(p*p + mass*mass));
		return result;
	}
	float TrackOperator::GetOffset(ReconstructedParticle * particle)
	{
		if (!particle || particle->getTracks().size() < 1) 
		{
			std::cout << "The particle is null or 0 tracks!\n";
			return 0.0;
		}
		Track * track = particle->getTracks()[0];
		for (unsigned int i = 0; i < particle->getTracks().size(); i++) 
		{
			//PrintTrack(particle->getTracks()[i]);
		}
		float d0 = track->getD0();
		float z0 = track->getZ0();
		float offset = std::sqrt(d0*d0 + z0*z0);
		return offset;
	}
	float TrackOperator::GetDistanceBtw(const Vertex * sec, const Vertex * ip, const ReconstructedParticle * particle)
	{
		double * trackPos = GetStartPoint(particle);
		double * secPos = MathOperator::toDoubleArray(sec->getPosition(), 3);
		double * ipPos = MathOperator::toDoubleArray(ip->getPosition(), 3);
		double secDir[3];
		for (int i = 0; i < 3; i++) 
		{
			secDir[i] = secPos[i] - ipPos[i];
		}
		vector<float> secDirection = MathOperator::getDirection(secDir);
		vector<float> trackDirection = MathOperator::getDirection(particle->getMomentum());
		float distance = MathOperator::getDistanceBtw(trackPos, trackDirection, secPos, secDirection);
		return distance;
	}
	float TrackOperator::GetError(ReconstructedParticle * particle)
	{
		const vector<float> matrix = particle->getTracks()[0]->getCovMatrix();
		float sigma = std::sqrt( matrix[0]*matrix[0] + matrix[9]*matrix[9] );
		return sigma;
	}
	double * TrackOperator::GetStartPoint(const ReconstructedParticle * particle)
	{
		double * start = new double[3];
		Track * track = particle->getTracks()[0];
		float d0 = track->getD0();
		float z0 = track->getZ0();
		float phi0 = track->getPhi();
		start[0] =  - d0 * std::sin(phi0);
		start[1] =  d0 * std::cos(phi0);
		start[2] = z0;
		return start;

	}
	float TrackOperator::GetDprime(const EVENT::ReconstructedParticle * particle1, const EVENT::ReconstructedParticle * particle2, double * primaryPosition)
	{
		double * secPos = GetStartPoint(particle1);
		vector< float > secDir = MathOperator::getDirection(particle1->getMomentum());
		//float epsilon = MathOperator::getDistanceTo(primaryPosition, secDir, secPos);
		
		double * secPos1 = GetStartPoint(particle2);
		vector< float > secDir1 = MathOperator::getDirection(particle2->getMomentum());
		//float epsilon1 = MathOperator::getDistanceTo(primaryPosition, secDir1, secPos1);
		double side1[3];
		double side2[3];
		double vr[3];
		for (int i = 0; i < 3; i++) 
		{
			side1[i] = primaryPosition[i] - secPos[i];
			side2[i] = primaryPosition[i] - secPos1[i];
			vr[i] = secPos[i] - secPos1[i];
		}

		float r = MathOperator::getDistance(secPos1, secPos);
		//float tan = std::tan(MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum()) / 2.0);
		float sinGamma = std::sin(MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum()));
		float sinAlpha = std::sin(MathOperator::getAngle(particle1->getMomentum(), vr));
		float sinBeta = std::sin(MathOperator::getAngle(particle2->getMomentum(), vr));

		float dstar = r * sinAlpha * sinBeta / sinGamma;
		//float dstar = r / 2.0 / tan;
		//float median = std::sqrt(epsilon*epsilon + epsilon1*epsilon1 - r ) / 2.0;
		float median = std::sqrt(MathOperator::getModule(side1)*MathOperator::getModule(side1)*2+MathOperator::getModule(side2)*MathOperator::getModule(side2)*2 - r*r) /2.0;
		float dprime = std::sqrt(dstar * dstar + median * median);

		return dprime;
	}

	float TrackOperator::GetOffsetError(EVENT::ReconstructedParticle * particle, double * trackPosition, const EVENT::Vertex * ipVertex, double offset)
	{
		double p = MathOperator::getModule(particle->getMomentum());
		double m = offset / p;
		double * vec = new double[3];
		const float * ip = ipVertex->getPosition();
		for (int i = 0; i < 3; i++) 
		{
			vec[i] = trackPosition[i] - ip[i]; 
		}
		GConfig conf = {particle->getMomentum(), trackPosition, ip, m, p, vec};
		const vector< float > pcovMatrix = particle->getCovMatrix();
		const vector< float > ccovMatrix = getErrorPoint(particle);
		const vector< float > ipcovMatrix = ipVertex->getCovMatrix();
		/*std::cout <<"\n!!!cCovMatrix:\n";
		for (unsigned int i = 0; i < ccovMatrix.size(); i++) 
		{
			std::cout <<  i << ": " << ccovMatrix[i] << ' ';
		}
		std::cout <<"\n";*/
		return getError(conf, pcovMatrix, ccovMatrix, ipcovMatrix);
	}
	float TrackOperator::getError(GConfig & conf, const vector< float > pcovMatrix, const vector< float > ccovMatrix, const vector< float > ipcovMatrix)
	{
		float result = 0.0;
		printConf(conf);
		int index = -1;
		for (int i = 0; i < 3; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				if (j <= i) 
				{
					index++;
					//float ipError =(std::isnan(ipcovMatrix[index])) ? 0.0: ipcovMatrix[index];
					float add = doffsetdp(conf, i) *  doffsetdp(conf, j) * pcovMatrix[index];// + doffsetdC(conf, i) *  doffsetdC(conf, j) * ipError;
					result += (i != j)? 2*add : add;
					//result += doffsetdp(conf, i) *  doffsetdp(conf, j) * pcovMatrix[index];
					//result += doffsetdC(conf, i) *  doffsetdC(conf, j) * ipError;
				}
			}
		}
		result += doffsetdC(conf, 0) *  doffsetdC(conf, 0) * ccovMatrix[0];
		result += (doffsetdC(conf, 0) *  doffsetdC(conf, 1) + doffsetdC(conf, 1) *  doffsetdC(conf, 0)) * ccovMatrix[1];
		result += doffsetdC(conf, 1) *  doffsetdC(conf, 1) * ccovMatrix[2];
		result += doffsetdC(conf, 2) *  doffsetdC(conf, 2) * ccovMatrix[3];
		return result;
	}
	void TrackOperator::PrintTrack(Track * track)
	{
		if (!track) 
		{
			return;
		}
		std::cout << "|\t" << track->getD0() 
		       << "|\t" << track->getZ0() 
		       << "|\t" << track->getPhi() 
		       << "|\t" << track->getOmega() 
		       << "|\t" << track->getTanLambda()
		       << "|\t" << MathOperator::getModule(track->getReferencePoint ())

		       << '\n';
	}

	

	void TrackOperator::printConf(GConfig & conf)
	{
		/*std::cout << "Configuration print: \n";
		std::cout << "|\t" << conf.Momentum[0]
			<< "|\t" << conf.Momentum[1]
			<< "|\t" << conf.Momentum[2] << '\n'
			<< "|\t" << conf.C[0]
			<< "|\t" << conf.C[1]
			<< "|\t" << conf.C[2] << '\n'
			<< "|\t" << conf.IP[0]
			<< "|\t" << conf.IP[1]
			<< "|\t" << conf.IP[2]<< '\n'
			<< "|\t" << conf.m
			<< "|\t" << conf.p
			<< "|\t" << conf.v[0]
			<< "|\t" << conf.v[1]
			<< "|\t" << conf.v[2] << '\n';*/

	}
	void TrackOperator::test()
	{
		double * p = new double[3];
		p[0] = 45.0; p[1] = 78.2; p[2] = 0.5;
		double * p0 = new double[3];
		p0[0] = 5.0; p0[1] = 7.52; p0[2] = 10.5;
		float * ip = new float[3];
		ip[0] = 545.0; ip[1] = 547.52; ip[2] = 510.5;
		GConfig conf = {p, p0, ip,0.0,0.01, p};
		printConf(conf);
	}
	void TrackOperator::test(EVENT::ReconstructedParticle * particle)
	{
		double * C = GetStartPoint(particle);
		std::cout <<"\nPoint:\n";
		for (int i = 0; i < 3; i++) 
		{
			std::cout <<  i << ": " << C[i] << ' ';
		}
		std::cout <<"\nCovMatrix:\n";
		const vector<float> covMatrix = getErrorPoint(particle);	
		for (unsigned int i = 0; i < covMatrix.size(); i++) 
		{
			std::cout <<  i << ": " << std::sqrt(covMatrix[i]) << ' ';
		}
		std::cout <<"\n";
	}
	const vector<float> TrackOperator::getErrorPoint(const EVENT::ReconstructedParticle * particle) const
	{
		//double * C = GetStartPoint(particle);
		Track * track = particle->getTracks()[0];
		const vector<float> covMatrix = track->getCovMatrix();
		float d0 = track->getD0();
		//float z0 = track->getZ0(); //CHECK
		float phi0 = track->getPhi();
		
		float dxdd0 = - std::sin(phi0);
		float dydd0 = - std::cos(phi0);
		float dxdphi0 = - d0 * std::cos(phi0);
		float dydphi0 = - d0 * std::cos(phi0);

		vector<float> result;
		result.push_back(dxdd0 * dxdd0 * covMatrix[0] + 2 * dxdphi0 * dxdd0 * covMatrix[1] + dxdphi0 * dxdphi0 * covMatrix[2]); // xx
		result.push_back(dxdd0 * dydd0 * covMatrix[0] + dydphi0 * dxdd0 * covMatrix[1] + dxdphi0 * dydd0 * covMatrix[1] + dxdphi0 * dydphi0 * covMatrix[2]); // xy
		result.push_back(dydd0 * dydd0 * covMatrix[0] + 2 * dydphi0 * dydd0 * covMatrix[1] + dydphi0 * dydphi0 * covMatrix[2]); // yy
		result.push_back(covMatrix[9]); // zz

		return result;
	}
	
	double TrackOperator::doffsetdC(GConfig & conf, int i)
	{
		double result = 0.0;
		switch(i)
		{
			case 0:
				result = ( (conf.v[0]*conf.Momentum[1] - conf.v[1]*conf.Momentum[0])*conf.Momentum[1] - (conf.v[0]*conf.Momentum[2] - conf.v[2]*conf.Momentum[0])*conf.Momentum[2] ) / conf.m / conf.p;
				break;
			case 1:
				result = ( (conf.v[1]*conf.Momentum[2] - conf.v[2]*conf.Momentum[1])*conf.Momentum[2] - (conf.v[0]*conf.Momentum[1] - conf.v[1]*conf.Momentum[0])*conf.Momentum[0] ) / conf.m / conf.p;
				break;
			case 2:
				result = ( (conf.v[2]*conf.Momentum[0] - conf.v[0]*conf.Momentum[2])*conf.Momentum[0] - (conf.v[1]*conf.Momentum[2] - conf.v[2]*conf.Momentum[1])*conf.Momentum[1] ) / conf.m / conf.p;
				break;
		}
		result /= conf.p;
		return result;
	}
	double TrackOperator::doffsetdp(GConfig & conf, int i)
	{
		double result = 0.0;
		switch(i)
		{
			case 0:
				result = ( (conf.v[0]*conf.Momentum[1] - conf.v[1]*conf.Momentum[0])*conf.v[1] - (conf.v[0]*conf.Momentum[2] - conf.v[2]*conf.Momentum[0])*conf.v[2] ) / conf.m / conf.p;
				break;
			case 1:
				result = ( (conf.v[1]*conf.Momentum[2] - conf.v[2]*conf.Momentum[1])*conf.v[2] - (conf.v[0]*conf.Momentum[1] - conf.v[1]*conf.Momentum[0])*conf.v[0] ) / conf.m / conf.p;
				break;
			case 2:
				result = ( (conf.v[2]*conf.Momentum[0] - conf.v[0]*conf.Momentum[2])*conf.v[0] - (conf.v[1]*conf.Momentum[2] - conf.v[2]*conf.Momentum[1])*conf.v[1] ) / conf.m / conf.p;
				break;
		}
		result -= conf.m * conf.Momentum[i] / conf.p;
		result /= conf.p * conf.p;
		return result;
	}
}
