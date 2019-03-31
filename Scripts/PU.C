/*
 * =====================================================================================
 *
 *       Filename:  PU.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Wednesday 13 July 2016 06:01.19  CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ramkrishna Sharma (Ram), ramkrishna.sharma71@gmail.com
 *   Organization:  University Of Delhi, Delhi, India
 *
 * =====================================================================================
 */

void PU() {
TGaxis::SetMaxDigits(3);
TFile f1("/uscms_data/d3/rasharma/aQGC_analysis/SecondStep/CMSSW_8_0_26_patch1/src/WWAnalysis/WWAnalysisRun2/output/WWTree_SingleElectron_Data.root");
TFile f2("/uscms_data/d3/rasharma/aQGC_analysis/SecondStep/CMSSW_8_0_26_patch1/src/WWAnalysis/WWAnalysisRun2/WWTree_TTToSemilepton_powheg_EleMu_EleMu.root");

TTree *T1 = (TTree*)f1.Get("otree");
TTree *T2 = (TTree*)f2.Get("otree");


TH1D *h1 = new TH1D("h1","Data",70,0,70);
h1->SetTitle("");
h1->SetStats(0);
h1->SetMarkerStyle(3);
T1->Draw("nPV>>h1");
//T1->Draw("nPV>>h1","pu_Weight");
h1->GetYaxis()->SetRange(0,90000);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitle("Normalized to Number of Events in Data");
h1->GetXaxis()->SetTitle("nPV");
//h1->GetXaxis()->SetTitle("nPV","pfMET>50");
//h1->Scale(1/h1->Integral());


TH1D *h2 = new TH1D("h2","WJets",70,0,70);
h2->SetTitle("");
h2->SetStats(0);
h2->SetMarkerStyle(2);
h2->SetLineColor(2);
h2->SetMarkerColor(2);
h2->GetYaxis()->SetTitleOffset(1.1);
h2->GetYaxis()->SetTitle("Normalized to Number of Events in Data");
h2->GetXaxis()->SetTitle("nPV");
//h2->SetFillColor(2);
h2->GetYaxis()->SetRange(0,90000);
T2->Draw("nPV>>h2");
//T2->Draw("nPV>>h2","pu_Weight");
//T2->Draw("nPV>>h2","totalEventWeight");
h2->Scale((h1->Integral())/(h2->Integral()));

TH1D *h3 = new TH1D("h3","WJets",70,0,70);
h3->SetTitle("");
h3->SetStats(0);
h3->SetMarkerStyle(2);
h3->SetLineColor(3);
h3->SetMarkerColor(3);
h3->GetYaxis()->SetTitleOffset(1.1);
h3->GetYaxis()->SetTitle("Normalized to Number of Events in Data");
h3->GetXaxis()->SetTitle("nPV");
//h3->SetFillColor(2);
h3->GetYaxis()->SetRange(0,90000);
//T2->Draw("nPV>>h3");
T2->Draw("nPV>>h3","pu_Weight");
//T2->Draw("nPV>>h3","totalEventWeight");
h3->Scale((h1->Integral())/(h3->Integral()));


TH1D *h4 = new TH1D("h4","WJets",70,0,70);
h4->SetTitle("");
h4->SetStats(0);
h4->SetMarkerStyle(2);
h4->SetLineColor(4);
h4->SetMarkerColor(4);
h4->GetYaxis()->SetTitleOffset(1.1);
h4->GetYaxis()->SetTitle("Normalized to Number of Events in Data");
h4->GetXaxis()->SetTitle("nPV");
//h4->SetFillColor(2);
h4->GetYaxis()->SetRange(0,90000);
//T2->Draw("nPV>>h4");
T2->Draw("nPV>>h4","pu_Weight*id_eff_Weight");
//T2->Draw("nPV>>h4","totalEventWeight");
h4->Scale((h1->Integral())/(h4->Integral()));

TH1D *h5 = new TH1D("h5","WJets",70,0,70);
h5->SetTitle("");
h5->SetStats(0);
h5->SetMarkerStyle(2);
h5->SetLineColor(5);
h5->SetMarkerColor(5);
h5->GetYaxis()->SetTitleOffset(1.1);
h5->GetYaxis()->SetTitle("Normalized to Number of Events in Data");
h5->GetXaxis()->SetTitle("nPV");
//h5->SetFillColor(2);
h5->GetYaxis()->SetRange(0,90000);
//T2->Draw("nPV>>h5");
T2->Draw("nPV>>h5","pu_Weight_2");
//T2->Draw("nPV>>h5","totalEventWeight");
h5->Scale((h1->Integral())/(h5->Integral()));


TCanvas *c1 = new TCanvas("c1");
c1->cd();

h1->Draw("hist");
h2->Draw("sames hist");
h3->Draw("sames hist");
h4->Draw("sames hist");
h5->Draw("sames hist");
TLegend * leg = new TLegend(0.5,0.7,0.90,0.9);
leg->AddEntry(h1,"Data","l");
leg->AddEntry(h2,"TTbar","l");
leg->AddEntry(h3,"TTbar with PU corr","l");
leg->AddEntry(h5,"TTbar with PU2 corr","l");
leg->AddEntry(h4,"TTbar with PU+IDIso","l");
leg->Draw();

c1->SaveAs("pileup.png");

//
//TFile f("/tmp/PileUpHisto.root","recreate");
//
// h1->Scale(1./h1->Integral());
// h2->Scale(1./h2->Integral());
// h1->Divide(h2);
//
//h1->Write();
////h2->Write();
//f.Close();
}

