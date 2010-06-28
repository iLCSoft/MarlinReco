// BcEnergyDensity class by A.Sapronov (sapronov@cern.ch)
// 18/07/2008

#include "iostream"

#include "TFile.h"
#include "TTree.h"

#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

#include "BcEnergyDensity.h"

using namespace ROOT::Math;

void BcEnergyDensity::Init(const char* inputfilename){
   //const char filename[100] = "bg_aver_LDC_4T_14mrad_AntiDID.root.root";
   TFile f(inputfilename, "READ");
   if ( f.IsZombie() ) {
      cerr << "Could not read data file. Exit." << endl;
      exit(1);
   }
   TTree *t = (TTree*) f.Get("tBcDensAverage");
   if ( t->IsZombie() ) {
      cerr << "Could not find data tree \"tBcDensAverage\". Exit." << endl;
      exit(1);
   }

   Int_t    cell_id[msSize];
   Double_t r_max;
   Double_t r_min;
   Double_t s_phi;
   Double_t d_phi;
   Double_t en_dens;
   Double_t en_dens_err;

   t->SetBranchAddress("sPos", cell_id);
   t->SetBranchAddress("sRin", &r_min);   // cm
   t->SetBranchAddress("sRout", &r_max);  // cm
   t->SetBranchAddress("sSphi", &s_phi);  // rad
   t->SetBranchAddress("sDphi", &d_phi);  // rad
   t->SetBranchAddress("sEnDens", &en_dens); // GeV/mm2
   t->SetBranchAddress("sEnDensErr", &en_dens_err);

   Int_t ne = t->GetEntries();

   for (int ie = 0; ie<ne; ie++){
      t->GetEntry(ie);
      r_min *= 10; // switch to mm
      r_max *= 10;
      CellType *cell = new CellType;
      for (int i = 0; i<msSize; i++){ cell->Id[i] = cell_id[i]; }
      cell->R = 0.5*( r_max + r_min );
      cell->Phi = s_phi + 0.5*d_phi;
      cell->EnDens = en_dens;
      cell->EnDensErr = en_dens_err;
      
      mCellStorage.push_back(cell);

      if ( (cell_id[0] == 1) && (cell_id[2] == 0) ) {
         mSegmentDeltaPhi.push_back(d_phi);
         if ( cell_id[1] == 0 ) {
	   mRMin = r_min;
	   mPhiMin = s_phi;
	   mRingDeltaR = r_max - r_min;
	 }
      }
   }

   mNumberOfRings = mSegmentDeltaPhi.size();
   mNumbersOfSegments.assign(mNumberOfRings, 0);

   for (int ie = 0; ie<ne; ie++){
      t->GetEntry(ie);
      r_min *= 10; // switch to mm
      r_max *= 10;

      if ( cell_id[1] == ( mNumberOfRings - 1 ) ){
         mRMax = r_max;
      }

      if ( cell_id[0] == 1 ) {
         mNumbersOfSegments.at(cell_id[1]) += 1;
      }
   }

   f.Close();
}

void BcEnergyDensity::Destroy(){
   vector<CellType*>::iterator it;
   for (it = mCellStorage.begin(); it!=mCellStorage.end(); it++){
      delete *it;
   }
}


Bool_t BcEnergyDensity::GetEnergyDensity(const Int_t& rLayer, 
                         const Double_t& rRadius, const Double_t& rPhi,
			 Double_t* pEnDens, Double_t* pEnDensError) const
{
  // initialize to zero such that 0 is returned in keyhole region!
  *pEnDens = 0.;
  *pEnDensError = 0.;
//   const Double_t TWO_PI = 2*TMath::Pi();
   Double_t the_phi = rPhi;
//   if ( the_phi < 0. ) the_phi += TWO_PI;
//   if ( the_phi > TWO_PI) the_phi -= TWO_PI;
   Int_t the_cell = 0;
   the_cell = GetCellNumber(rLayer, rRadius, the_phi); // central cell
   if ( the_cell < 0 ) return false;
   Int_t the_ring = mCellStorage.at(the_cell)->Id[1];

   const Int_t nsegs = 3;
   const Int_t nrings = 3;
   Int_t i_cell[nrings][nsegs];
   // nrings: 0 - down ring, 1 - original ring, 2 - up ring
   // nsegs:  0 - prev segment, 1 - original, 2 - next segment
   for ( int ir = 0; ir<nrings; ir++ ){
      for ( int is = 0; is<nsegs; is++ ){
         Double_t radius = rRadius + (ir - 1)*mRingDeltaR;
	 Int_t i_ring = the_ring + ir - 1;
	 if ( i_ring < 0 || i_ring >= mNumberOfRings ){
	    i_cell[ir][is] = -1;
	 } else {
	    Double_t phi = the_phi + (is - 1)*mSegmentDeltaPhi.at(i_ring);
            i_cell[ir][is] = GetCellNumber(rLayer, radius, phi);
	 }
      }
   }

   vector<double> phis[nrings];
   vector<ValErrType> en_dens[nrings];

   for (int is = 0; is<nsegs; is++){
      for (int ir = 0; ir<nrings; ir++){
         if ( !( i_cell[ir][is]<0 ) ){
            CellType *cell = mCellStorage.at(i_cell[ir][is]);
	    Double_t phi = cell->Phi;
	    if ( the_phi < mSegmentDeltaPhi.at(the_ring + ir - 1) && is == 0 ) {
	       phi -= 2*TMath::Pi();
	    }
	    if ( 2*TMath::Pi() - the_phi < mSegmentDeltaPhi.at(the_ring + ir - 1) 
	      && is == 2 ) {
	       phi += 2*TMath::Pi();
	    }
            phis[ir].push_back(phi);
	    ValErrType ed;
	    ed.Val = cell->EnDens;
	    Double_t ed_err = cell->EnDensErr;
	    ed.Lb = ed.Val - ed_err;
	    ed.Ub = ed.Val + ed_err;
            en_dens[ir].push_back(ed);
	 }
      }
   }

   // interpolate over phi in rings original, down and up 
   vector<double> r;
   vector<ValErrType> ed_int_by_phi; // values interpolated by phi in diff. rings
   for (int ir=0; ir<nrings; ir++){
      if ( en_dens[ir].size() != 0 && i_cell[ir][1] > 0 ){
         r.push_back(mCellStorage.at(i_cell[ir][1])->R);
         ValErrType ed_int;
	 Bool_t int_result = MakeInterpolation(the_phi, phis[ir], en_dens[ir], &ed_int);
	 if( int_result ){
	    ed_int_by_phi.push_back(ed_int);
	 } else {
	    cerr << "Interpolation failed." << endl;
	    return false;
	 }
      }
   }

   if (ed_int_by_phi.size()!=0){
      ValErrType ed_int;
      Bool_t int_result = MakeInterpolation(rRadius, r, ed_int_by_phi, &ed_int);
      if ( int_result ) {
         *pEnDens = ed_int.Val;
	 *pEnDensError = 0.5*(ed_int.Ub - ed_int.Lb);
      } else {
         cerr << "Interpolation failed." << endl;
         return false;
      }
   } else {
      return false;
   }

   return true;
}


// returns absolute number of cell which contains given coords
Int_t BcEnergyDensity::GetCellNumber(const Int_t& rLayer, const Double_t& rRadius, 
                                     const Double_t& rPhi) const
{
   unsigned int cell_cnt = 0; // cell counter

   if ( rRadius < mRMin || rRadius > mRMax ) return -1;

   Int_t layer_number =  TMath::Abs(rLayer) - 1;  // should start with 0
   Int_t layer_sign = ( rLayer>0 ) ? 0 : 1; // 0 - forward, 1 - backward
   Int_t ring_number = (Int_t) TMath::Floor( (rRadius - mRMin)/mRingDeltaR );
   double phi = ( rPhi < mPhiMin) ? rPhi + 2*TMath::Pi() : rPhi;
   Int_t segm_number = 
      (Int_t) TMath::Floor( (phi - mPhiMin)/mSegmentDeltaPhi.at(ring_number) );
   if ( segm_number >= mNumbersOfSegments.at(ring_number) ) return -1;

   Int_t n_in_layer = 0;
   for (int ir=0; ir<mNumberOfRings; ir++){
      n_in_layer += mNumbersOfSegments.at(ir);
   }

   for (int ir=0; ir<ring_number; ir++){
      cell_cnt += mNumbersOfSegments.at(ir);
   }
   cell_cnt += n_in_layer*layer_number + segm_number;
   cell_cnt *= 2;
   cell_cnt += layer_sign;

   if ( cell_cnt > mCellStorage.size() ) return -1;

   return cell_cnt;
}

Bool_t BcEnergyDensity::MakeInterpolation ( const Double_t& rX,
                                            const vector<Double_t>& rXData, 
                                            const vector<ValErrType>& rYData,
					    ValErrType* pResult) const
{					    
/* cout<<rXData.size()<<" "<<rYData.size()<<endl;
   cout << rX <<"\t";
   for (size_t i=0; i<rXData.size(); i++){
      cout << rXData.at(i) << "\t";
   }
   for (size_t i=0; i<rXData.size(); i++){
      cout << rYData.at(i).Val << "\t";
   }
   cout<<endl;
*/


#if ROOT_VERSION_CODE < ROOT_VERSION(5,20,0)
   Interpolation::Type interp_type = Interpolation::CSPLINE;
   if (rXData.size()<3) interp_type = Interpolation::LINEAR;
#else
   Interpolation::Type interp_type = Interpolation::kCSPLINE;
   if (rXData.size()<3) interp_type = Interpolation::kLINEAR;
#endif


//    Interpolator interp(rXData.size(), interp_type);
   Interpolator *interp=0;

   //cout<<interp->TypeGet()<<endl;
   vector<double> y_data;
   for (size_t i=0; i<rYData.size(); i++) y_data.push_back(rYData.at(i).Val);
   interp = new Interpolator (rXData, y_data, interp_type);
   //interp.SetData(rXData, y_data);
   if (interp > 0) {
     pResult->Val = interp->Eval(rX);
     delete interp;
   }
   y_data.clear();

   for (size_t i=0; i<rYData.size(); i++) y_data.push_back(rYData.at(i).Lb);
   interp = new Interpolator (rXData, y_data, interp_type);
   //interp.SetData(rXData, y_data);
   if (interp >0) {
     pResult->Lb = interp->Eval(rX);
     delete interp;
   }
   y_data.clear();

   for (size_t i=0; i<rYData.size(); i++) y_data.push_back(rYData.at(i).Ub);
   interp = new Interpolator (rXData, y_data, interp_type);
   //interp.SetData(rXData, y_data);
   if (interp > 0) {
   pResult->Ub = interp->Eval(rX);
   delete interp;
   }

   y_data.clear();

   return true;
}
