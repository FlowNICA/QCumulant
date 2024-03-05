#define PLOTV2EventPlane
#include "DrawTGraph.C"
#include "../constants.C"
void SaveGraphEventPlane(TString inputFileName = "GenericFrameworkPlotter/ToyModel_QC.root",
                      TString methodName = "MC",
                      Bool_t saveAsPNG = true,
                      TString outDirName = "pics")
// methodName can be "EtaSub", "SP", "FHCalEP", "LYZEP", "Eta3Sub", "MC"
{

  TFile *fi = TFile::Open(inputFileName.Data());
  TProfile2D *prV2CentEta = (TProfile2D *)fi->FindObjectAny(Form("prV2%svsEta",methodName.Data()));
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  TGraphErrors *gr[ncent][npid], *grRF[npid], *gr1040[npid];
  TProfile *tmp, *tmp1, *tmp2;
  const std::vector<TString> pidFancyNames = {"h^{+}", "#pi^{+}", "K^{+}", "p", "h^{-}", "#pi^{-}", "K^{-}", "#bar{p}", "h^{#pm}","#pi^{#pm}","K^{#pm}","p(#bar{p})"};
  for (Int_t i=0; i < npid; i++)
  {
    TProfile2D *prV2CentPt = (TProfile2D *)fi->FindObjectAny(Form("prV2%svsPt_pid%i",methodName.Data(),i));
    
    TProfile *tmp1 = PlotPtIntegratedV2(prV2CentPt,0.2,3.0);// v2 versus centrality, 0.2<pt<3.0 GeV/c
    grRF[i] = Converter(tmp1);
    grRF[i]->GetYaxis()-> SetTitle("v_{2}");
    grRF[i]->GetXaxis()-> SetTitle("Centrality, %");
    // grRF->GetXaxis()->SetRangeUser(0., 80);
    grRF[i]->SetMarkerStyle(20);
    grRF[i]->SetMarkerColor(kRed);

    tmp2 = PlotV2vsPt(prV2CentPt,10,40);// v2 versus pt, 10-40%
    gr1040[i] = Converter(tmp2);
    gr1040[i]->GetYaxis()-> SetTitle("v_{2}");
    gr1040[i]->SetMarkerStyle(20);
    gr1040[i]->SetMarkerColor(kRed);
    gr1040[i]->GetXaxis()-> SetTitle("p_{T}, GeV/c");

    for (Int_t j=0; j<ncent; j++)
    {
      tmp = PlotV2vsPt(prV2CentPt,bin_cent[j],bin_cent[j+1]);// v2 versus pt, 10-40%
      gr[j][i] = Converter(tmp);
      gr[j][i]->GetYaxis()-> SetTitle("v_{2}");
      // gr[j][i]->GetYaxis()->SetRangeUser(0., 0.2);
      gr[j][i]->SetMarkerStyle(20);
      gr[j][i]->SetMarkerColor(kRed);
      gr[j][i]->GetXaxis()-> SetTitle("p_{T}, GeV/c");
      // gr[j][i]->GetXaxis()->SetRangeUser(0., 3.5); 
    }
  }
  TFile *fo = new TFile(Form("graphs_v2_%s.root",methodName.Data()),"recreate");
  fo->cd();
  for (Int_t i=0; i < npid; i++)
  {
    grRF[i]->Write(Form("grRF_%i",i));
    gr1040[i]->Write(Form("gr_cent10-40_%i",i));
    for (Int_t j=0; j<ncent; j++)
    {
      gr[j][i]->Write(Form("gr_cent%i_%i",j,i));
    }
  }
  fo->Close();
}