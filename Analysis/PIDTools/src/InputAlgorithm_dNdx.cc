#include <InputAlgorithm_dNdx.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DD4hep/DD4hepUnits.h>
#include <IMPL/TrackImpl.h>
#include <SimpleHelix.h>
#include <TRandom.h>
#include <TSpline.h>


using namespace lcio;

namespace cpid {

  InputAlgorithm_dNdx anInputAlgorithm_dNdx;

  InputAlgorithm_dNdx::InputAlgorithm_dNdx(): InputAlgorithm("dNdx")
  {
    _description = "Returns dE/dx value of PFO track";
  }

  std::vector<std::string> InputAlgorithm_dNdx::init(const std::vector<float>& inparF, const std::vector<std::string>& inparS)
  {
    print_init(inparF, inparS);

    if (inparF.size()>0) _gas = int(inparF[0]);
    else _gas = 0;
    if (inparF.size()>1) _scaling = inparF[1];
    else _scaling = 1;

    bool TPCfound = false;
    try
    {
      dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
      double bfieldV[3];
      lcdd.field().magneticField({0., 0., 0.}, bfieldV);
      _bField = bfieldV[2]/dd4hep::tesla;

      try
      {
        dd4hep::DetElement tpcDE = lcdd.detector("TPC");
        const dd4hep::rec::FixedPadSizeTPCData* tpc = tpcDE.extension<dd4hep::rec::FixedPadSizeTPCData>() ;
        _innRad = tpc->rMinReadout/dd4hep::mm;
        _outRad = tpc->rMaxReadout/dd4hep::mm;
        _zLen = tpc->driftLength/dd4hep::mm;
        TPCfound = true;
      }
      catch (...)
      {}
    }

    catch (...)
    {
      if (inparF.size()>2) _bField = inparF[2];
      else _bField = -3.5;
    }

    if (!TPCfound)
    {
      if (inparF.size()>5)
      {
        _innRad = inparF[3];
        _outRad = inparF[4];
        _zLen   = inparF[5];
      }
      else
      {
        _innRad =  372;
        _outRad = 1692;
        _zLen   = 2225;
      }
    }

    sloM << "init parameters used: gas = " << _gas << ", scaling = " << _scaling << ", B field = " << _bField << " T, track inner radius = " << _innRad << " mm, track outer radius = " << _outRad << " mm, tracks max z = " << _zLen << " mm" << std::endl;

    std::vector<std::string> obsNames{"elDis","muDis","piDis","kaDis","prDis"};
    return obsNames;
  }

  std::vector<std::pair<float,float> > InputAlgorithm_dNdx::extractObservables(ReconstructedParticleImpl* pfo, LCCollection* , int mcpdg)
  {
    std::vector<std::pair<float,float> > obsValues{};
    std::vector<std::pair<float,float> > failDefault{};
    for (int i=0; i<5; ++i) failDefault.push_back(std::pair<float,float>{-1e+10,-1});

    float masses[] = {0.000511, 0.10566, 0.13957, 0.49368, 0.93827};
    float mass = 0;
    switch (mcpdg)
    {
    case   11: mass = 0.000511; break;
    case   13: mass = 0.10566;  break;
    case  211: mass = 0.13957;  break;
    case  321: mass = 0.49368;  break;
    case 2212: mass = 0.93827;  break;
    }

    TrackVec tracks = pfo->getTracks();
    if (tracks.size() > 0 && mass)
    {
      Track* trk = tracks[0];

      const float* refPoint = trk->getReferencePoint();
      LCVector3D refPointV3D = LCVector3D(refPoint[0],refPoint[1],refPoint[2]);
      SimpleHelix trkHelix = SimpleHelix(trk->getD0(),trk->getPhi(),trk->getOmega(),trk->getZ0(),trk->getTanLambda(),refPointV3D);
      HelixClass trkHelixC;
      trkHelixC.Initialize_Canonical(trk->getPhi(),trk->getD0(),trk->getZ0(),trk->getOmega(),trk->getTanLambda(),_bField);

      const float* refP = trkHelixC.getReferencePoint();
      float ref[] = {refP[0], refP[1], refP[2]};
      float phi0 = computePointInTracker(trkHelixC, ref, _innRad, _zLen);
      if (phi0<-1e+10)
      {
        //obsValues.push_back(std::pair<float,float>{-1,-1});
        return failDefault;
      }
      float phi1 = computePointInTracker(trkHelixC, ref, _outRad, _zLen);
      if (phi1<-1e10)
      {
        const TrackState* tslh = trk->getTrackState(TrackState::AtLastHit);
        phi1 = tslh->getPhi();
      }

      const TrackState* tsIP = trk->getTrackState(TrackState::AtIP);
      float Omega = tsIP->getOmega();
      float tanL  = tsIP->getTanLambda();
      float len   = abs((phi0-phi1)*(1/Omega) * sqrt(1+tanL*tanL));

      double mom = sqrt(pow(trkHelixC.getMomentum()[0],2) + pow(trkHelixC.getMomentum()[1],2) + pow(trkHelixC.getMomentum()[2],2));
      double bg  = mom/mass;
      float  mCl = nClusters(bg, _gas)/10. * _scaling * len;
      float  nCl = gRandom->PoissonD(mCl);
      float dNdx = nCl/len;

      sloD << mass << " " << mom << " " << bg << " | " << phi0 << " " << phi1 << " " << len << " | " << mCl << " " << nCl << " " << dNdx << std::endl;

      for (int i=0; i<5; ++i)
      {
        double bgRef = mom/masses[i];
        float dNdxRef = nClusters(bgRef, _gas)/10. * _scaling;
        std::pair<float,float> dNdxValue {(dNdx-dNdxRef), sqrt(nCl)/len};
        if (nCl==0 || dNdxRef==0) dNdxValue.first = -1e+10;
        sloD << " " << dNdx-dNdxRef;
        obsValues.push_back(dNdxValue);
      }
      sloD << std::endl;
    }
    else return failDefault;

    return obsValues;
  }

  // returns the number of clusters per cm for a given gas at a given beta*gamma (for a pressure of 1 atm)
  float InputAlgorithm_dNdx::nClusters(float bg, int gas)
  {
    // gas = 0: He 90 - Isobutane 10
    //     = 1: pure He
    //     = 2: Argon 50 - Ethane 50
    //     = 3: pure Argon

    const int n = 18;
    double bgs[n] = {0.5, 0.8, 1., 2., 3., 4., 5., 8., 10., 12., 15., 20., 50., 100., 200., 500., 1000., 10000.};

    double nCl_He_Iso[n] = {42.94, 23.6, 18.97, 12.98, 12.2, 12.13, 12.24, 12.73, 13.03, 13.29, 13.63, 14.08, 15.56, 16.43, 16.8, 16.95, 16.98, 16.98};
    double nCl_He[n]     = {11.79, 6.5, 5.23, 3.59, 3.38, 3.37, 3.4, 3.54, 3.63, 3.7, 3.8, 3.92, 4.33, 4.61, 4.78, 4.87, 4.89, 4.89};
    double nCl_Ar_Eth[n] = {130.04, 71.55, 57.56, 39.44, 37.08, 36.9, 37.25, 38.76, 39.68, 40.49, 41.53, 42.91, 46.8, 48.09, 48.59, 48.85, 48.93, 48.93};
    double nCl_Ar[n]     = {88.69, 48.93, 39.41, 27.09, 25.51, 25.43, 25.69, 26.78, 27.44, 28.02, 28.77, 29.78, 32.67, 33.75, 34.24, 34.57, 34.68, 34.68};

    double nCl[n];
    switch (gas)
    {
    case 0: std::copy(nCl_He_Iso, nCl_He_Iso + n, nCl);
      break;
    case 1: std::copy(nCl_He, nCl_He + n, nCl);
      break;
    case 2: std::copy(nCl_Ar_Eth, nCl_Ar_Eth + n, nCl);
      break;
    case 3: std::copy(nCl_Ar, nCl_Ar + n, nCl);
      break;
    }

    float interp = 0.0;
    if (bg>1000) interp = nCl[17];
    else
    {
      TSpline3* sp3 = new TSpline3("sp3", bgs, nCl, n);
      if (bg > bgs[0] && bg < bgs[n - 1]) interp = sp3->Eval(bg);
    }
    return interp;
  }

  float InputAlgorithm_dNdx::computePointInTracker(HelixClass& helix, float* ref, float radius, float zmax)
  {
    float phi=-1e20;
    float points[6];
    float success = helix.getPointOnCircle(radius, ref, points);
    if (success>-1e20)
    {
      if (abs(points[2]>zmax && abs(points[5])>zmax))
      {
        helix.getPointInZ(zmax, ref, points);
        if (points[1]==0) phi = points[0]<0 ? M_PI*0.5 : M_PI*1.5;
        else phi = atan(-points[0]/points[1]);
      }
      else
        if (abs(points[2])<abs(points[5]))
          if (points[1]==0) phi = points[0]<0 ? M_PI*0.5 : M_PI*1.5;
          else phi = atan(-points[0]/points[1]);
        else
          if (points[4]==0) phi = points[3]<0 ? M_PI*0.5 : M_PI*1.5;
          else phi = atan(-points[3]/points[4]);
    }
    return phi;
  }

}
