#define Calculate_CrossSection_cxx
#include "Calculate_CrossSection.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

void Calculate_CrossSection::Loop(std::string outFileName, float Total_Cross_Section)
{
//   In a ROOT session, you can do:
//      root> .L Calculate_CrossSection.C
//      root> Calculate_CrossSection t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   string Name[718] = {"FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS0", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FS1", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM0", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM1", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM6", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FM7", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT0", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT1", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2", "FT2"};
   double ParVal[718] = {-900.0, -880.0, -860.0, -840.0, -820.0, -800.0, -780.0, -760.0, -740.0, -720.0, -700.0, -680.0, -660.0, -640.0, -620.0, -600.0, -580.0, -560.0, -540.0, -520.0, -500.0, -480.0, -460.0, -440.0, -420.0, -400.0, -380.0, -360.0, -340.0, -320.0, -300.0, -280.0, -260.0, -240.0, -220.0, -200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0, 820.0, 840.0, 860.0, 880.0, 900.0, -330.0, -320.0, -310.0, -300.0, -290.0, -280.0, -270.0, -260.0, -250.0, -240.0, -230.0, -220.0, -210.0, -200.0, -190.0, -180.0, -170.0, -160.0, -150.0, -140.0, -130.0, -120.0, -110.0, -100.0, -90.0, -80.0, -70.0, -60.0, -50.0, -40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, -42.0, -41.0, -40.0, -39.0, -38.0, -37.0, -36.0, -35.0, -34.0, -33.0, -32.0, -31.0, -30.0, -29.0, -28.0, -27.0, -26.0, -25.0, -24.0, -23.0, -22.0, -21.0, -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, -165.0, -160.0, -155.0, -150.0, -145.0, -140.0, -135.0, -130.0, -125.0, -120.0, -115.0, -110.0, -105.0, -100.0, -95.0, -90.0, -85.0, -80.0, -75.0, -70.0, -65.0, -60.0, -55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, -84.0, -82.0, -80.0, -78.0, -76.0, -74.0, -72.0, -70.0, -68.0, -66.0, -64.0, -62.0, -60.0, -58.0, -56.0, -54.0, -52.0, -50.0, -48.0, -46.0, -44.0, -42.0, -40.0, -38.0, -36.0, -34.0, -32.0, -30.0, -28.0, -26.0, -24.0, -22.0, -20.0, -18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, 66.0, 68.0, 70.0, 72.0, 74.0, 76.0, 78.0, 80.0, 82.0, -300.0, -295.0, -290.0, -285.0, -280.0, -275.0, -270.0, -265.0, -260.0, -255.0, -250.0, -245.0, -240.0, -235.0, -230.0, -225.0, -220.0, -215.0, -210.0, -205.0, -200.0, -195.0, -190.0, -185.0, -180.0, -175.0, -170.0, -165.0, -160.0, -155.0, -150.0, -145.0, -140.0, -135.0, -130.0, -125.0, -120.0, -115.0, -110.0, -105.0, -100.0, -95.0, -90.0, -85.0, -80.0, -75.0, -70.0, -65.0, -60.0, -55.0, -50.0, -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, -6.8, -6.6, -6.4, -6.2, -6.0, -5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8, -3.6, -3.4, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, -12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, -20.5, -20.0, -19.5, -19.0, -18.5, -18.0, -17.5, -17.0, -16.5, -16.0, -15.5, -15.0, -14.5, -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5};


   ofstream file;
   if (outFileName.c_str() == "")
   	outFileName = "temp_data.txt";

   file.open (outFileName.c_str()); 
   file<<"#aQGC_ParName\tParaValue\tCross-section"<<endl;
   //file<<"SR.NO.\taQGC_ParName\tParaValue\tSumofWeight*Cross-Sec\tAvg.OFWeight*Cross-section"<<endl;

   fChain->SetBranchStatus("*",0); 
   fChain->SetBranchStatus("LHEWeight.weight",1);
   fChain->SetBranchStatus("LHEWeight.id",1);

   Long64_t nentries = fChain->GetEntries();
   cout<<"total number of entries = "<<nentries<<endl;

   Long64_t nbytes = 0, nb = 0;
   double SumWeight[718] = {0};
   double AvgWeight[718] = {0};
   double DevWeight[718] = {0};
   int countWgt ;

   cout<<"Total cross-section = "<<Total_Cross_Section<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //cout<<"Size = "<<sizeof(LHEWeight_id)/sizeof(LHEWeight_id[0])<<endl;
      //cout<<LHEWeight_id[0]<<"\t"<<LHEWeight_weight[0]<<endl;
      //cout<<LHEWeight_id[1]<<"\t"<<LHEWeight_weight[1]<<endl;
      //cout<<LHEWeight_id[2]<<"\t"<<LHEWeight_weight[2]<<endl;
      //cout<<"================="<<endl;
      // if (Cut(ientry) < 0) continue;
      countWgt = 0;
      //if (jentry<50)
      //	cout<<"LHEWeight_id[0] = "<<LHEWeight_weight[0]<<endl;
      // Ref: 
      for (int i =446; i<sizeof(LHEWeight_id)/sizeof(LHEWeight_id[0]);i++)
      {
      	SumWeight[countWgt] += (LHEWeight_weight[i])/(LHEWeight_weight[0]);
	//if (i==446+670)
	//cout<<i<<"\t"<<ParVal[i-446]<<"\tLHEWeight_weight["<<i<<"] = "<<LHEWeight_weight[i]<<"\t"<<LHEWeight_weight[0]<<"\t"<<(LHEWeight_weight[i])/(LHEWeight_weight[0])<<"\t"<<SumWeight[countWgt]<<endl;
	countWgt += 1;
      }
   }
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      countWgt = 0;
      for (int i =446; i<sizeof(LHEWeight_id)/sizeof(LHEWeight_id[0]);i++)
      {
      	AvgWeight[countWgt] = SumWeight[countWgt]/nentries;
	countWgt += 1;
      }
   }
   for(int i=0;i<718;i++){
   	//if (i==670)
   	cout<<i<<"\t"<<Name[i]<<"\t"<<ParVal[i]<<"\t"<<SumWeight[i]*Total_Cross_Section<<"\t"<<AvgWeight[i]*Total_Cross_Section<<endl;
   	file<<Name[i]<<"\t"<<ParVal[i]<<"\t"<<AvgWeight[i]*Total_Cross_Section<<endl;
   	//file<<i<<"\t"<<Name[i]<<"\t"<<ParVal[i]<<"\t"<<SumWeight[i]*Total_Cross_Section<<"\t"<<AvgWeight[i]*Total_Cross_Section<<endl;
   	//cout<<i<<"\t"<<Name[i]<<"\t"<<ParVal[i]<<"\t"<<SumWeight[i]<<"\t"<<AvgWeight[i]<<endl;
   	//file<<i<<"\t"<<Name[i]<<"\t"<<ParVal[i]<<"\t"<<SumWeight[i]<<"\t"<<AvgWeight[i]<<endl;
   }
file.close();
}
