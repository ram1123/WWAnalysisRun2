import ROOT as r

r.gROOT.SetBatch(True)

file1 = r.TFile("WWTree_WplusToLNuWminusTo2JJJ_EWK_LO_SM_MJJ100PTJ10_TuneCUETP8M1_13TeV-madgraph-pythia8_Save.root","READ")

treeIn = file1.Get("otree")

c1 = r.TCanvas("c1","c1",800,600);
c1.SetLogy();

h1 = r.TH1F("h1","",100,-70,70);
h1.SetLineColor(1)
h2 = r.TH1F("h2","",100,-70,70);
h2.SetLineColor(2)
h3 = r.TH1F("h3","",100,-70,70);
h3.SetLineColor(3)
h4 = r.TH1F("h4","",100,-70,70);
h4.SetLineColor(4)

#treeIn.Draw("nu_pz_gen>>h1")
treeIn.Draw("(nu_pz_type0-nu_pz_gen)/nu_pz_gen>>h2")
treeIn.Draw("(nu_pz_type2-nu_pz_gen)/nu_pz_gen>>h3")
treeIn.Draw("(nu_pz_run2-nu_pz_gen)/nu_pz_gen>>h4")


leg = r.TLegend(0.7,0.7,0.95,0.95)

h2.SetMaximum(13000)
h2.GetXaxis().SetTitle("(Sol - Truth)/Truth")
#h1.Draw()
#leg.AddEntry(h1,"Generator Value");
h2.Draw()
leg.AddEntry(h2,"type0")
h3.Draw("same")
leg.AddEntry(h3,"type2")
h4.Draw("same")
leg.AddEntry(h4,"run2")

leg.Draw();

c1.SaveAs("neutrino_pz.png")
c1.SaveAs("neutrino_pz.pdf")
