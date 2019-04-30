#include "../interface/analysisUtils.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1F.h"

double getDeltaPhi(double phi1, double phi2 ){
   const double PI = 3.14159265;
   double result = phi1 - phi2;

   if(result > PI) result = result - 2 * PI;
   if(result <= (-1 * PI)) result = result + 2 * PI;

   result = TMath::Abs(result);
   return result;
}

double deltaPhi(const double& phi1, const double& phi2)
{
   double deltaphi = fabs(phi1 - phi2);
   if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
   if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
   return deltaphi;
}
double deltaEta(const double& eta1, const double& eta2)
{
   double deltaeta = fabs(eta1 - eta2);
   return deltaeta;
}

double deltaR(const double& eta1, const double& phi1,
      const double& eta2, const double& phi2)
{
   double deltaphi = deltaPhi(phi1, phi2);
   double deltaeta = deltaEta(eta1, eta2);
   double deltar = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);
   return deltar;
}


//////////////////////////////////
//Ref: https://github.com/ram1123/LHEAnalyzer/blob/LHEanalyzer/LHEanalyzer.cpp
//////////////////////////////////
void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){

   //std::cout<<"After: "<<thep4H.Pt() << " " << thep4Z1.Pt() << " " << thep4M11.Pt() << " " <<thep4M12.Pt() << " "  << thep4Z2.Pt() << " " << thep4M21.Pt() << " " << thep4M22.Pt() << std::endl;
   ///////////////////////////////////////////////
   // check for z1/z2 convention, redefine all 4 vectors with convention
   ///////////////////////////////////////////////	
   TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
   p4H = thep4H;

   p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
   p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
   //// costhetastar
   TVector3 boostX = -(thep4H.BoostVector());
   TLorentzVector thep4Z1inXFrame( p4Z1 );
   TLorentzVector thep4Z2inXFrame( p4Z2 );	
   thep4Z1inXFrame.Boost( boostX );
   thep4Z2inXFrame.Boost( boostX );
   TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
   TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
   costhetastar = theZ1X_p3.CosTheta();

   //// --------------------------- costheta1
   TVector3 boostV1 = -(thep4Z1.BoostVector());
   TLorentzVector p4M11_BV1( p4M11 );
   TLorentzVector p4M12_BV1( p4M12 );	
   TLorentzVector p4M21_BV1( p4M21 );
   TLorentzVector p4M22_BV1( p4M22 );
   p4M11_BV1.Boost( boostV1 );
   p4M12_BV1.Boost( boostV1 );
   p4M21_BV1.Boost( boostV1 );
   p4M22_BV1.Boost( boostV1 );

   TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
   //// costheta1
   costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();

   if (boostV1.Mag()  > 1.0)
      std::cout<<"p4V2 = "<< p4V2_BV1.Vect().Mag() << "\t" << p4M11_BV1.Vect().Mag() << "\t Boost = " << boostV1.Mag() << std::endl;
   if (isnan( (float) costheta1) == 1)
      std::cout<<"p4V2 = "<< p4V2_BV1.Vect().Mag() << "\t" << p4M11_BV1.Vect().Mag() << "\t Boost = " << boostV1.Mag() << std::endl;

   //// --------------------------- costheta2
   TVector3 boostV2 = -(thep4Z2.BoostVector());
   TLorentzVector p4M11_BV2( p4M11 );
   TLorentzVector p4M12_BV2( p4M12 );	
   TLorentzVector p4M21_BV2( p4M21 );
   TLorentzVector p4M22_BV2( p4M22 );
   p4M11_BV2.Boost( boostV2 );
   p4M12_BV2.Boost( boostV2 );
   p4M21_BV2.Boost( boostV2 );
   p4M22_BV2.Boost( boostV2 );

   TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
   //// costheta2
   costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
   if (isnan( (float) costheta2) == 1)
      std::cout << "boostV2 = " << boostV2.Mag() << "\tp4M11_BV2 = " << p4M11_BV2.Mag() << "\tp4M11_BV2 = "<<p4M11_BV2.Mag() << "\tp4M12_BV2 = " << p4M12_BV2.Mag() << "\tp4M21_BV2 = " << p4M21_BV2.Mag() << "\tp4M22_BV2 = " << p4M22_BV2.Mag() << std::endl;

   //// --------------------------- Phi and Phi1
   //    TVector3 boostX = -(thep4H.BoostVector());
   TLorentzVector p4M11_BX( p4M11 );
   TLorentzVector p4M12_BX( p4M12 );	
   TLorentzVector p4M21_BX( p4M21 );
   TLorentzVector p4M22_BX( p4M22 );	

   p4M11_BX.Boost( boostX );
   p4M12_BX.Boost( boostX );
   p4M21_BX.Boost( boostX );
   p4M22_BX.Boost( boostX );

   TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
   TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    

   TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
   TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 

   //// Phi
   TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
   double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
   double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
   Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );


   //////////////

   TVector3 beamAxis(0,0,1);
   TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();

   TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
   TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
   TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );

   //// Phi1
   double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
   double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
   Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    

   //    std::cout << "extractAngles: " << std::endl;
   //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
   //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;    

}

double GetSFs_Lepton(double xAxis, double yAxis, TH1F *h1, TString posX, TString posY){
   double scaleFactor = 0.0;
   if (posX=="pt" && posY=="eta"){
      // Check if xAxis is not ouside upper limit; if so then to assign same SF as the last bin to high xAxis resacle tempPt to less then maximum range of xAxis defined in histo.
      if (xAxis > h1->GetXaxis()->GetXmax())  	
	 xAxis = h1->GetXaxis()->GetXmax() - 1.0;	
      if (xAxis < h1->GetXaxis()->GetXmin()) 
	 xAxis = h1->GetXaxis()->GetXmin() + 1.0;

      //std::cout << "====> pt = " <<xAxis<< "\t"<< h1->GetXaxis()->FindFixBin(xAxis) << "\t" << yAxis<< "\t" << h1->GetYaxis()->FindFixBin(yAxis) << "\n" << std::endl;

      scaleFactor = h1->GetBinContent(h1->GetXaxis()->FindFixBin(xAxis), h1->GetYaxis()->FindFixBin(yAxis));
   }
   else if (posX=="eta" && posY=="pt"){
      if (yAxis > h1->GetYaxis()->GetXmax())  	
	 yAxis = h1->GetYaxis()->GetXmax() - 1.0;	
      if (yAxis < h1->GetYaxis()->GetXmin()) 
	 yAxis = h1->GetYaxis()->GetXmin() + 1.0;
      //std::cout << "====> eta = " <<xAxis<< "\t"<< h1->GetXaxis()->FindFixBin(xAxis) << "\t" << yAxis<< "\t" << h1->GetYaxis()->FindFixBin(yAxis) << "\n" << std::endl;
      scaleFactor = h1->GetBinContent(h1->GetXaxis()->FindFixBin(xAxis), h1->GetYaxis()->FindFixBin(yAxis));
   }
   return scaleFactor;
}

double GetMin(double x, double y){
   if (x<y) return x;
   else	 return y;
}
double GetMax(double x, double y){
   if (x>y) return x;
   else	 return y;
}

float getPUPPIweight(TF1* puppisd_corrGEN, TF1* puppisd_corrRECO_cen, TF1* puppisd_corrRECO_for, float puppipt, float puppieta ){
   float genCorr  = 1.;
   float recoCorr = 1.;
   float totalWeight = 1.;

   genCorr =  puppisd_corrGEN->Eval( puppipt );
   if( fabs(puppieta)  <= 1.3 ){
      recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
   }
   else{
      recoCorr = puppisd_corrRECO_for->Eval( puppipt );
   }

   totalWeight = genCorr * recoCorr;


   return totalWeight;
}
double func( double pt, double eta, JetCorrectionUncertainty *fJetUnc) { 
   fJetUnc->setJetPt ( pt  );
   fJetUnc->setJetEta( eta );
   return fJetUnc->getUncertainty(true);
}
double GetPt_MET(double pfMET, double phi, double pz){

   double px = pfMET*TMath::Cos(phi);
   double py = pfMET*TMath::Sin(phi);
   return TMath::Sqrt(px*px + py*py); 
}
double GetEta_MET(double pfMET, double phi, double pz){
   double px = pfMET*TMath::Cos(phi);
   double py = pfMET*TMath::Sin(phi);
   double p = TMath::Sqrt(px*px + py*py + pz*pz);
   return 0.5*log((p+pz)/(p-pz));
}
