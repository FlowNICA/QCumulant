#include "MultiparticleCorrelation.cxx"
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include "../DrawTGraph.C"
#include "../../constants.C"

#include <vector>
#include <iostream>
enum CovTerm {k24, k26, k46, kNCov};
enum DiffCovTerm {k22prime, k24prime, k26prime, k2prime4, k2prime4prime, k2prime6, k2prime6prime, k44prime, k46prime, k4prime6, k4prime6prime, k66prime, kNDiffCov};
enum Corr {k2, k4, k6, kNCorr};
enum DiffCorr {k2prime, k4prime, k6prime, kNDiffCorr};
enum Method {kStd, kGap, kNMethod};
// Flags
Bool_t saveAsPNG = true;
Int_t ratioToMethod = 0; // 2QC, 4QC, 6QC, 2QC-gapped, 4QC-gapped, 6QC-gapped
Int_t drawDifferentialFlowTill = 0; // Draw v2 vs pT (10% centrality cut) till: 0: no drawing; 1: till 10%; 2: till 20%; etc.
const std::vector<TString> pidFancyNames = {"h^{+}", "#pi^{+}", "K^{+}", "p", "h^{-}", "#pi^{-}", "K^{-}", "#bar{p}", "h^{#pm}","#pi^{#pm}","K^{#pm}","p(#bar{p})"};
const float eta_gap = 0.0;
TString model = "ToyModel";
// TString modelFancy = "AMPT SM, #sigma_{p}=1.5";
TString modelFancy = "ToyModel";
// TString modelFancy = model;
TString energy = "";
TString inputFileName = Form("./FirstRun.root");

const Int_t nmethod = 6; // 2QC, 4QC, 6QC, 2QC-gapped, 4QC-gapped, 6QC-gapped
const Int_t binMinPtRFP = 1;  // 0.2 GeV
const Int_t binMaxPtRFP = npt-1; // 2.8 GeV

const Double_t minptRFP = 0.2;
const Double_t maxptRFP = 3.0;

const Double_t maxpt = 3.6; // for v2 vs pt plotting
const Double_t minpt = 0.0;  // for v2 vs pt plotting

// const Int_t ncent = 9; // 0-80 %
const Double_t centrality_bin_center[ncent] = {2.5,7.5,15,25,35,45,55,65,75};
const Double_t bin_centE[ncent] = {0};
const std::vector<std::pair<Int_t,Int_t>> centrality = {{0,5},{5,10},{10,20},{20,30},{30,40},{40,50},{50,60},{60,70},{70,80}};


const Double_t mincent = 0.;  // for v2 vs centrality
const Double_t maxcent = 80.; // for v2 vs centrality

const Double_t minV2int = -0.005; // for v2 vs centrality plotting
const Double_t maxV2int = 0.1; // for v2 vs centrality plotting
const Double_t minV2dif = -0.01; // for v2 vs pt plotting
const Double_t maxV2dif = 0.25; // for v2 vs pt plotting


std::vector <Double_t> coordinateLeg = {0.18,0.63,0.45,0.889};
// std::vector<std::pair<Double_t,Double_t>> rangeRatio = {{0.84,1.16},{0.93,1.07}};
std::vector<std::pair<Double_t,Double_t>> rangeRatioRF ={{0.54,1.13},{0.54,1.13}};
std::vector<std::pair<Double_t,Double_t>> rangeRatio = {{0.67,1.13},{0.67,1.13}};
// std::vector<std::pair<Double_t,Double_t>> rangeRatioRF = {{0.78,1.22},{0.78,1.22}};
Int_t marker[nmethod]={kFullSquare,kFullCircle,kFullTriangleDown,kOpenSquare,kOpenCircle,kOpenTriangleDown}; // 2QC, 4QC, 6QC, 2QC-gapped, 4QC-gapped, 6QC-gapped

void CalStatErrCent1040(Double_t v2eDif1040[nmethod][npid][npt]){

  TFile *inFile = new TFile(inputFileName.Data(),"read");

  TProfile* pCov[ncent][kNMethod][kNCov];
  TProfile* pCovDiff[ncent][npt][npid][kNMethod][kNDiffCov];
  TProfile* pCorr[ncent][kNMethod][kNCorr];
  TProfile* pCorrDiff[ncent][npt][npid][kNMethod][kNDiffCorr];
  TProfile* pCorrWSquare[ncent][kNMethod][kNCorr];                     // for stat. err.
  TProfile* pCorrDiffWSquare[ncent][npt][npid][kNMethod][kNDiffCorr];  // for stat. err.

  for (Int_t icent = 0; icent < ncent; icent++)
  { // loop over centrality classes
    for(Int_t imeth(0);imeth<kNMethod;imeth++)
    {
      for (Int_t i(0);i<kNCorr;i++) pCorr[icent][imeth][i]        = (TProfile*) inFile->Get(Form("QCGF/pCorr_%i_%i_%i",icent,imeth,i));
      for (Int_t i(0);i<kNCorr;i++) pCorrWSquare[icent][imeth][i] = (TProfile*) inFile->Get(Form("QCGF/pCorrWSquare_%i_%i_%i",icent,imeth,i));
      for (Int_t i(0);i<kNCov;i++)  pCov[icent][imeth][i]         = (TProfile*) inFile->Get(Form("QCGF/pCov_%i_%i_%i",icent,imeth,i));
      for (Int_t id = 0; id < npid; id++)
      {
        for (Int_t ipt = 0; ipt < npt; ipt++)
        { // loop over pt bin
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]         = (TProfile*) inFile->Get(Form("QCGF/pCorrDiff_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i));
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]  = (TProfile*) inFile->Get(Form("QCGF/pCorrDiffWSquare_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i));
          for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][ipt][id][imeth][i]          = (TProfile*) inFile->Get(Form("QCGF/pCovDiff_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i));
        }
      }
    }
  } // end of loop over centrality classes


  // for (Int_t icent=0;icent<ncent;icent++) {
  //   for(Int_t imeth(0);imeth<kNMethod;imeth++) {
  //     for (Int_t ipt=0;ipt<npt;ipt++) {
  //       for (Int_t id=8;id<npid;id++) {
  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]         = (TProfile*) pCorrDiff[icent][ipt][id-8][imeth][i]         ->Clone();
  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]  = (TProfile*) pCorrDiffWSquare[icent][ipt][id-8][imeth][i]  ->Clone();
  //         for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][ipt][id][imeth][i]          = (TProfile*) pCovDiff[icent][ipt][id-8][imeth][i]          ->Clone();

  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]         ->Add(pCorrDiff[icent][ipt][id-4][imeth][i]);
  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]  ->Add(pCorrDiffWSquare[icent][ipt][id-4][imeth][i]);
  //         for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][ipt][id][imeth][i]          ->Add(pCovDiff[icent][ipt][id-4][imeth][i]);
  //       }
  //     }
  //   }
  // }


  // Add
  for (Int_t icent=3; icent<5; icent++){ // add 20-30% & 30-40% to 10-20%
    for(Int_t imeth(0);imeth<kNMethod;imeth++) {

      for (Int_t i(0);i<kNCorr;i++) pCorr[2][imeth][i]        ->Add(pCorr[icent][imeth][i]       );
      for (Int_t i(0);i<kNCorr;i++) pCorrWSquare[2][imeth][i] ->Add(pCorrWSquare[icent][imeth][i]);
      for (Int_t i(0);i<kNCov;i++)  pCov[2][imeth][i]         ->Add(pCov[icent][imeth][i]        );
      for(Int_t ipt=0; ipt<npt; ipt++){ // loop over pt bin
        for (Int_t id=0;id<npid;id++){ // loop over pid
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[2][ipt][id][imeth][i]         ->Add(pCorrDiff[icent][ipt][id][imeth][i]       );
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[2][ipt][id][imeth][i]  ->Add(pCorrDiffWSquare[icent][ipt][id][imeth][i]);
          for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[2][ipt][id][imeth][i]          ->Add(pCovDiff[icent][ipt][id][imeth][i]        );
        }
      } // end of loop over pt bin
    }
  }
  Int_t zero = 0;
  for (Int_t icent=2; icent<3; icent++){ // 10-40
    for (Int_t id=0;id<npid;id++){
      for(Int_t ipt=0; ipt<npt; ipt++){
        MultiparticleCorrelation *multCorr[kNMethod];
        for(Int_t imeth(0);imeth<kNMethod;imeth++) {
          multCorr[imeth] = new MultiparticleCorrelation(pCorr[icent][imeth][k2], pCorr[icent][imeth][k4], pCorr[icent][imeth][k6],
          pCorrDiff[icent][ipt][id][imeth][k2], pCorrDiff[icent][ipt][id][imeth][k4], pCorrDiff[icent][ipt][id][imeth][k6],
          pCorrWSquare[icent][imeth][k2], pCorrWSquare[icent][imeth][k4], pCorrWSquare[icent][imeth][k6],
          pCorrDiffWSquare[icent][ipt][id][imeth][k2], pCorrDiffWSquare[icent][ipt][id][imeth][k4], pCorrDiffWSquare[icent][ipt][id][imeth][k6],
          pCovDiff[icent][ipt][id][imeth][k22prime], pCov[icent][imeth][k24], pCovDiff[icent][ipt][id][imeth][k24prime],
          pCov[icent][imeth][k26], pCovDiff[icent][ipt][id][imeth][k26prime],
          pCovDiff[icent][ipt][id][imeth][k2prime4], pCovDiff[icent][ipt][id][imeth][k2prime4prime], pCovDiff[icent][ipt][id][imeth][k2prime6],
          pCovDiff[icent][ipt][id][imeth][k2prime6prime], pCovDiff[icent][ipt][id][imeth][k44prime], pCov[icent][imeth][k46], pCovDiff[icent][ipt][id][imeth][k46prime],
          pCovDiff[icent][ipt][id][imeth][k4prime6], pCovDiff[icent][ipt][id][imeth][k4prime6prime], pCovDiff[icent][ipt][id][imeth][k66prime],zero,zero,"");
        }
        v2eDif1040[0][id][ipt] = multCorr[0]->GetV22DifErr();
        v2eDif1040[1][id][ipt] = multCorr[0]->GetV24DifErr();
        v2eDif1040[2][id][ipt] = multCorr[0]->GetV26DifErr();
        v2eDif1040[3][id][ipt] = multCorr[1]->GetV22DifErr();
        v2eDif1040[4][id][ipt] = multCorr[1]->GetV24DifErr();
        v2eDif1040[5][id][ipt] = multCorr[1]->GetV26DifErr();
      } // end of loop for all pT bin
    } // end of loop for PID
  } // end of loop for centrality
  inFile->Close();
}

void PlotV2QCumulant(){
  TFile *outFile = new TFile(Form("./graphs_QCGF_%s.root",model.Data()),"recreate");
  TString outDirName=(TString)Form("%s_%s_eta_gap_%1.1f",model.Data(),energy.Data(),eta_gap*2);
  TString level= (TString) Form("%s, Au+Au at #sqrt{s_{NN}}=%s",modelFancy.Data(),energy.Data());
  Double_t v2eDif1040[nmethod][npid][npt];
  CalStatErrCent1040(v2eDif1040);

  TFile *inFile = new TFile(inputFileName.Data(),"read");

  // Temporary variables
  char hname[800]; // histogram hname

  // Input hist
  TProfile* pCov[ncent][kNMethod][kNCov];
  TProfile* pCovDiff[ncent][npt][npid][kNMethod][kNDiffCov];
  TProfile* pCorr[ncent][kNMethod][kNCorr];
  TProfile* pCorrDiff[ncent][npt][npid][kNMethod][kNDiffCorr];
  TProfile* pCorrWSquare[ncent][kNMethod][kNCorr];                     // for stat. err.
  TProfile* pCorrDiffWSquare[ncent][npt][npid][kNMethod][kNDiffCorr];  // for stat. err.

  TProfile *hcounter[ncent][npt][npid];

  // OUTPUT
  TGraphErrors *grDifFl[nmethod][ncent][npid];    // v2(pt); 3 = {2QC, 4QC, EP, gapped 2QC}
  TGraphErrors *grDifFl1040[nmethod][npid];


  for (Int_t icent = 0; icent < ncent; icent++)
  { // loop over centrality classes
    for (Int_t ipt=0;ipt<npt;ipt++) for (Int_t id=0;id<npid;id++) hcounter[icent][ipt][id] = (TProfile*) inFile->Get(Form("QCGF/hcounter_%i_%i_%i",icent,ipt,id));
    for (Int_t imeth(0);imeth<kNMethod;imeth++)
    {
      for (Int_t i(0);i<kNCorr;i++) pCorr[icent][imeth][i]        = (TProfile*) inFile->Get(Form("QCGF/pCorr_%i_%i_%i",icent,imeth,i));
      for (Int_t i(0);i<kNCorr;i++) pCorrWSquare[icent][imeth][i] = (TProfile*) inFile->Get(Form("QCGF/pCorrWSquare_%i_%i_%i",icent,imeth,i));
      for (Int_t i(0);i<kNCov;i++)  pCov[icent][imeth][i]         = (TProfile*) inFile->Get(Form("QCGF/pCov_%i_%i_%i",icent,imeth,i));
      for (Int_t id = 0; id < npid; id++)
      {
        for (Int_t ipt = 0; ipt < npt; ipt++)
        { // loop over pt bin
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]         = (TProfile*) inFile->Get(Form("QCGF/pCorrDiff_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i));
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]  = (TProfile*) inFile->Get(Form("QCGF/pCorrDiffWSquare_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i));
          for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][ipt][id][imeth][i]          = (TProfile*) inFile->Get(Form("QCGF/pCovDiff_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i));
        }
      }
    }
  } // end of loop over centrality classes

  for (Int_t icent = 0; icent < ncent; icent++)
  { // loop over centrality classes
    for (Int_t ipt=0;ipt<npt;ipt++) for (Int_t id=0;id<npid;id++) if (!hcounter[icent][ipt][id] ) { std::cout << "Err!!" << std::endl; }
    for (Int_t imeth(0);imeth<kNMethod;imeth++)
    {
      for (Int_t i(0);i<kNCorr;i++) if (!pCorr[icent][imeth][i]        ) { std::cout << "Err!!" << std::endl; }
      for (Int_t i(0);i<kNCorr;i++) if (!pCorrWSquare[icent][imeth][i] ) { std::cout << "Err!!" << std::endl; }
      for (Int_t i(0);i<kNCov;i++)  if (!pCov[icent][imeth][i]         ) { std::cout << "Err!!" << std::endl; }
      for (Int_t id = 0; id < npid; id++)
      {
        for (Int_t ipt = 0; ipt < npt; ipt++)
        { // loop over pt bin
          for (Int_t i(0);i<kNDiffCorr;i++) if (!pCorrDiff[icent][ipt][id][imeth][i]         ) { std::cout << "Err!!" << std::endl; }
          for (Int_t i(0);i<kNDiffCorr;i++) if (!pCorrDiffWSquare[icent][ipt][id][imeth][i]  ) { std::cout << "Err!!" << std::endl; }
          for (Int_t i(0);i<kNDiffCov;i++)  if (!pCovDiff[icent][ipt][id][imeth][i]          ) { std::cout << "Err!!" << std::endl; }
        }
      }
    }
  } // end of loop over centrality classes


  // for (Int_t icent=0;icent<ncent;icent++) {
  //   for (Int_t ipt=0;ipt<npt;ipt++) for (Int_t id=8;id<npid;id++)
  //   {
  //     hcounter[icent][ipt][id] = (TProfile*) hcounter[icent][ipt][id-8]->Clone();
  //     hcounter[icent][ipt][id] -> Add( hcounter[icent][ipt][id-4] );
  //   }
  //   for (Int_t imeth(0);imeth<kNMethod;imeth++) {
  //     for (Int_t ipt=0;ipt<npt;ipt++) {
  //       for (Int_t id=8;id<npid;id++) {
  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]         = (TProfile*) pCorrDiff[icent][ipt][id-8][imeth][i]         ->Clone();
  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]  = (TProfile*) pCorrDiffWSquare[icent][ipt][id-8][imeth][i]  ->Clone();
  //         for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][ipt][id][imeth][i]          = (TProfile*) pCovDiff[icent][ipt][id-8][imeth][i]          ->Clone();

  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]         ->Add(pCorrDiff[icent][ipt][id-4][imeth][i]);
  //         for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]  ->Add(pCorrDiffWSquare[icent][ipt][id-4][imeth][i]);
  //         for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][ipt][id][imeth][i]          ->Add(pCovDiff[icent][ipt][id-4][imeth][i]);
  //       }
  //     }
  //   }
  // }





  // Filling pT bin
  Double_t pt[npt];
  Double_t ept[npt]={0}; // error bin pT = 0.0
  for (Int_t ipt=0; ipt<npt; ipt++){
    // pt[icent][ipt][id] = hmult[icent][ipt][id] -> GetBinContent(1);
    pt[ipt] = ( pTBin[ipt] + pTBin[ipt+1] ) / 2.;
  }

  TProfile *prV2Dif1040[nmethod][npid];
  for (Int_t imeth=0; imeth<nmethod; imeth++){
    for (Int_t id=0; id<npid; id++){
      prV2Dif1040[imeth][id] = new TProfile(Form("prV2Dif1040_%i_%i",imeth,id),"",npt,0.,npt);
    }
  }
  // cross-check reference flow
  TProfile *prV2Ref[nmethod];
  for (Int_t imeth=0; imeth<nmethod; imeth++){
    prV2Ref[imeth] = new TProfile(Form("prV2Ref_%i",imeth),"",ncent,0,ncent);
  }
  Double_t v2Dif[nmethod][ncent][npid][npt], v2eDif[nmethod][ncent][npid][npt];
  Int_t zero = 0;
  for (Int_t icent=0; icent<ncent; icent++){ // 10-40
    for (Int_t id=0;id<npid;id++){
      for(Int_t ipt=0; ipt<npt; ipt++){
        MultiparticleCorrelation *multCorr[kNMethod];
        for(Int_t imeth(0);imeth<kNMethod;imeth++) {
          multCorr[imeth] = new MultiparticleCorrelation(pCorr[icent][imeth][k2], pCorr[icent][imeth][k4], pCorr[icent][imeth][k6],
          pCorrDiff[icent][ipt][id][imeth][k2], pCorrDiff[icent][ipt][id][imeth][k4], pCorrDiff[icent][ipt][id][imeth][k6],
          pCorrWSquare[icent][imeth][k2], pCorrWSquare[icent][imeth][k4], pCorrWSquare[icent][imeth][k6],
          pCorrDiffWSquare[icent][ipt][id][imeth][k2], pCorrDiffWSquare[icent][ipt][id][imeth][k4], pCorrDiffWSquare[icent][ipt][id][imeth][k6],
          pCovDiff[icent][ipt][id][imeth][k22prime], pCov[icent][imeth][k24], pCovDiff[icent][ipt][id][imeth][k24prime],
          pCov[icent][imeth][k26], pCovDiff[icent][ipt][id][imeth][k26prime],
          pCovDiff[icent][ipt][id][imeth][k2prime4], pCovDiff[icent][ipt][id][imeth][k2prime4prime], pCovDiff[icent][ipt][id][imeth][k2prime6],
          pCovDiff[icent][ipt][id][imeth][k2prime6prime], pCovDiff[icent][ipt][id][imeth][k44prime], pCov[icent][imeth][k46], pCovDiff[icent][ipt][id][imeth][k46prime],
          pCovDiff[icent][ipt][id][imeth][k4prime6], pCovDiff[icent][ipt][id][imeth][k4prime6prime], pCovDiff[icent][ipt][id][imeth][k66prime],zero,zero,"");
        }

        v2Dif[0][icent][id][ipt]  = multCorr[0]->GetV22Dif();
        v2eDif[0][icent][id][ipt] = multCorr[0]->GetV22DifErr();
        v2Dif[1][icent][id][ipt]  = multCorr[0]->GetV24Dif();
        v2eDif[1][icent][id][ipt] = multCorr[0]->GetV24DifErr();
        v2Dif[2][icent][id][ipt]  = multCorr[0]->GetV26Dif();
        v2eDif[2][icent][id][ipt] = multCorr[0]->GetV26DifErr();
        v2Dif[3][icent][id][ipt]  = multCorr[1]->GetV22Dif();
        v2eDif[3][icent][id][ipt] = multCorr[1]->GetV22DifErr();
        v2Dif[4][icent][id][ipt]  = multCorr[1]->GetV24Dif();
        v2eDif[4][icent][id][ipt] = multCorr[1]->GetV24DifErr();
        v2Dif[5][icent][id][ipt]  = multCorr[1]->GetV26Dif();
        v2eDif[5][icent][id][ipt] = multCorr[1]->GetV26DifErr();
        // if (id==8 && icent==3) std::cout << "v2{4}= " << v2Dif[1][icent][id][ipt] << "\t" << v2eDif[1][icent][id][ipt] << "v2{2,subevent}= " << v2Dif[2][icent][id][ipt] << "\t" << v2eDif[2][icent][id][ipt] << std::endl;

        if (icent>=2 && icent <=4) { // 10-40%

          prV2Dif1040[0][id] -> Fill(0.5+ipt,multCorr[0]->GetV22Dif(),hcounter[icent][ipt][id] -> GetBinEntries(1));
          prV2Dif1040[1][id] -> Fill(0.5+ipt,multCorr[0]->GetV24Dif(),hcounter[icent][ipt][id] -> GetBinEntries(1));
          prV2Dif1040[2][id] -> Fill(0.5+ipt,multCorr[0]->GetV26Dif(),hcounter[icent][ipt][id] -> GetBinEntries(1));
          prV2Dif1040[3][id] -> Fill(0.5+ipt,multCorr[1]->GetV22Dif(),hcounter[icent][ipt][id] -> GetBinEntries(2));
          prV2Dif1040[4][id] -> Fill(0.5+ipt,multCorr[1]->GetV24Dif(),hcounter[icent][ipt][id] -> GetBinEntries(2));
          prV2Dif1040[5][id] -> Fill(0.5+ipt,multCorr[1]->GetV26Dif(),hcounter[icent][ipt][id] -> GetBinEntries(2));
        }
        if (id == 8)
        {
          prV2Ref[0]->Fill(icent,multCorr[0]->GetV22Dif(),hcounter[icent][ipt][id] -> GetBinEntries(1));
          prV2Ref[1]->Fill(icent,multCorr[0]->GetV24Dif(),hcounter[icent][ipt][id] -> GetBinEntries(1));
          prV2Ref[2]->Fill(icent,multCorr[0]->GetV26Dif(),hcounter[icent][ipt][id] -> GetBinEntries(1));
          prV2Ref[3]->Fill(icent,multCorr[1]->GetV22Dif(),hcounter[icent][ipt][id] -> GetBinEntries(2));
          prV2Ref[4]->Fill(icent,multCorr[1]->GetV24Dif(),hcounter[icent][ipt][id] -> GetBinEntries(2));
          prV2Ref[5]->Fill(icent,multCorr[1]->GetV26Dif(),hcounter[icent][ipt][id] -> GetBinEntries(2));
        }

      } // end of loop for all pT bin

      for (Int_t i=0; i<nmethod; i++){
        grDifFl[i][icent][id] = new TGraphErrors(npt,pt,v2Dif[i][icent][id],ept,v2eDif[i][icent][id]);
        grDifFl[i][icent][id] -> SetMarkerStyle(marker[i]);
        grDifFl[i][icent][id] -> SetMarkerSize(1.5);
        grDifFl[i][icent][id] -> SetDrawOption("P");
      }

    } // end of loop for PID
  } // end of loop for centrality

  Double_t v2Dif1040[nmethod][npid][npt];
  for (Int_t imeth=0; imeth<nmethod; imeth++){
    for (Int_t id=0;id<npid;id++){
      for(Int_t ipt=0; ipt<npt; ipt++){
        v2Dif1040[imeth][id][ipt] = prV2Dif1040[imeth][id] -> GetBinContent(ipt+1);
      }
      grDifFl1040[imeth][id] = new TGraphErrors(npt,pt,v2Dif1040[imeth][id],ept,v2eDif1040[imeth][id]);
      grDifFl1040[imeth][id] -> SetMarkerStyle(marker[imeth]);
      grDifFl1040[imeth][id] -> SetMarkerSize(1.5);
      grDifFl1040[imeth][id] -> SetDrawOption("P");
    }
  }
  // for (Int_t ipt=0;ipt<npt;ipt++)
  // {
  //   std::cout << "v2{4} = " << v2eDif1040[1][8][ipt] << "\t" << "v2{2,sub-event} = " << v2eDif1040[2][8][ipt] << std::endl;
  // }

  const char *grTitle[nmethod]={"v_{2}{2}","v_{2}{4}","v_{2}{6}","v_{2}{2,|#Delta#eta|>0}","v_{2}{4,|#Delta#eta|>0}","v_{2}{6,|#Delta#eta|>0}"};
  outFile -> cd();
  for (Int_t imeth=0; imeth<nmethod; imeth++){
    for (Int_t id=0;id<npid;id++){
      grDifFl1040[imeth][id] -> SetTitle(grTitle[imeth]);
      grDifFl1040[imeth][id] -> GetYaxis()-> SetTitle("v_{2}");
      grDifFl1040[imeth][id] -> GetXaxis()-> SetTitle("p_{T}, GeV/c");
      grDifFl1040[imeth][id] -> Write(Form("gr_cent10-40_%i_%i",imeth,id));
      for (Int_t icent=0;icent<ncent;icent++){
        grDifFl[imeth][icent][id] -> SetTitle(grTitle[imeth]);
        grDifFl[imeth][icent][id] -> GetYaxis()-> SetTitle("v_{2}");
        grDifFl[imeth][icent][id] -> GetXaxis()-> SetTitle("p_{T}, GeV/c");
        grDifFl[imeth][icent][id] -> Write(Form("gr_cent%i_%i_%i",icent,imeth,id));
      }
    }
  }

  if (saveAsPNG) gSystem->Exec(Form("mkdir -p ./%s/",outDirName.Data()));
  TCanvas *cV2PT[ncent][npid];
  for (Int_t icent=0; icent<drawDifferentialFlowTill; icent++){
    for (Int_t id=0;id<npid;id++){
      std::vector<TGraphErrors*> vgrv2pt;
      vgrv2pt.push_back(grDifFl[ratioToMethod][icent][id]); // v2{gapped 2QC}
      for (Int_t i=0; i<nmethod; i++){
        if (i==ratioToMethod) continue;
        vgrv2pt.push_back(grDifFl[i][icent][id]);
      }
      sprintf(hname,"%s, %i-%i%%",pidFancyNames.at(id).Data(),centrality.at(icent).first,centrality.at(icent).second);
      cV2PT[icent][id] = (TCanvas*) DrawTGraph(vgrv2pt,"",rangeRatio.at(0).first, rangeRatio.at(0).second, minpt, maxpt, minV2dif, maxV2dif,
                                               coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                               level.Data(), hname, true, grTitle[ratioToMethod]);
      cV2PT[icent][id] -> SetName(hname);
      if (saveAsPNG) cV2PT[icent][id] -> SaveAs(Form("./%s/DifferentialFlow_Centrality%i-%i_%s.png",outDirName.Data(),centrality.at(icent).first,centrality.at(icent).second,pidNames.at(id).Data()));
    }
  }

  TCanvas *cV2PTMultPad[npid];
  // TString strCent[5] = {"0-5%","5-10%","10-20%","20-30%","30-40%"};
  TString strCent[5] = {"5-10%","10-20%","20-30%","30-40%","40-50%"};
  for (Int_t id=0;id<npid;id++)
  {
    std::vector<TGraphErrors*> vgrv2pt[5];
    for (Int_t icent=1; icent<6; icent++)
    {
      vgrv2pt[icent-1].push_back(grDifFl[ratioToMethod][icent][id]);
      for (Int_t imeth=0; imeth<nmethod; imeth++){
        if (imeth==ratioToMethod) continue;
        // if (imeth==0 || imeth==3 || imeth==4 || imeth==5) continue;
        // if (imeth==2 || imeth==4 || imeth==5) continue;
        vgrv2pt[icent-1].push_back(grDifFl[imeth][icent][id]);
      }
    }
    cV2PTMultPad[id] = (TCanvas*) DrawTGraph(vgrv2pt, 5,"",rangeRatio.at(0).first, rangeRatio.at(0).second, minpt, maxpt, minV2dif, maxV2dif,
                                              coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                              Form("%s, %s", level.Data(), pidFancyNames.at(id).Data()),
                                              strCent, true, grTitle[ratioToMethod]);
    cV2PTMultPad[id] -> SetName("");
    if (saveAsPNG) cV2PTMultPad[id] -> SaveAs(Form("./%s/DifferentialFlow_%s_Cent_0_40.png",outDirName.Data(),pidNames.at(id).Data()));
  }

  TCanvas *cV2PTMultPad2[npid];
  TString strCent2[4] = {"40-50%","50-60%","60-70%","70-80%"};
  for (Int_t id=0;id<npid;id++)
  {
    std::vector<TGraphErrors*> vgrv2pt[4];
    for (Int_t icent=5; icent<ncent; icent++)
    {
      vgrv2pt[icent-5].push_back(grDifFl[ratioToMethod][icent][id]);
      for (Int_t imeth=0; imeth<nmethod; imeth++){
        if (imeth==ratioToMethod) continue;
        // if (imeth==0 || imeth==3 || imeth==4 || imeth==5) continue;
        // if (imeth==2 || imeth==4 || imeth==5) continue;
        vgrv2pt[icent-5].push_back(grDifFl[imeth][icent][id]);
      }
    }
    cV2PTMultPad2[id] = (TCanvas*) DrawTGraph(vgrv2pt, 4,"",rangeRatio.at(0).first, rangeRatio.at(0).second, minpt, maxpt, minV2dif, maxV2dif,
                                              coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                              Form("%s, %s", level.Data(), pidFancyNames.at(id).Data()),
                                              strCent2, true, grTitle[ratioToMethod]);
    cV2PTMultPad2[id] -> SetName("");
    if (saveAsPNG) cV2PTMultPad2[id] -> SaveAs(Form("./%s/DifferentialFlow_%s_Cent_40_80.png",outDirName.Data(),pidNames.at(id).Data()));
  }

  TCanvas *cV2PT1040[npid];
  for (Int_t id=0;id<npid;id++){
    std::vector<TGraphErrors*> vgrv2pt1040;
    vgrv2pt1040.push_back(grDifFl1040[ratioToMethod][id]);
    for (Int_t imeth=0;imeth<nmethod;imeth++){
      if (imeth==ratioToMethod) continue;
      // if (imeth==0 || imeth==3 || imeth==4 || imeth==5) continue;
      // if (imeth==2 || imeth==4 || imeth==5) continue;
      vgrv2pt1040.push_back(grDifFl1040[imeth][id]);
    }
    sprintf(hname,"10-40%%, %s",pidFancyNames.at(id).Data());
    cV2PT1040[id] = (TCanvas*) DrawTGraph(vgrv2pt1040,"",rangeRatio.at(1).first, rangeRatio.at(1).second, minpt, 3.0, minV2dif, maxV2dif,
                                          coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                          level.Data(), hname, true, grTitle[ratioToMethod]);
    cV2PT1040[id] -> SetName(hname);
    if (saveAsPNG) cV2PT1040[id] -> SaveAs(Form("./%s/DifferentialFlow_Centrality10-40%%_%s.png",outDirName.Data(),pidNames.at(id).Data()));
    // if (saveAsPNG) cV2PT1040[id] -> SaveAs(Form("./%s/DifferentialFlow_Centrality10-40%%_%s.pdf",outDirName.Data(),pidNames.at(id).Data()));
  }
  //==========================================================================================================================

  TGraphErrors *grIntFlPID[nmethod][npid];
  TGraphErrors *grRefFl[nmethod];
  // v2 vs centrality for PID

  for (Int_t icent=0;icent<ncent;icent++) {
    // for (Int_t ipt=binMinPtRFP+1;ipt<binMaxPtRFP;ipt++) for (Int_t id=0;id<npid;id++) hcounter[icent][ipt][id] -> Add( hcounter[icent][ipt][id] );
    for (Int_t imeth(0);imeth<kNMethod;imeth++) {
      for (Int_t ipt=binMinPtRFP+1;ipt<binMaxPtRFP;ipt++){
        for (Int_t id=0;id<npid;id++) {
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][binMinPtRFP][id][imeth][i]         ->Add(pCorrDiff[icent][ipt][id][imeth][i]);
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][binMinPtRFP][id][imeth][i]  ->Add(pCorrDiffWSquare[icent][ipt][id][imeth][i]);
          for (Int_t i(0);i<kNDiffCov;i++)  pCovDiff[icent][binMinPtRFP][id][imeth][i]          ->Add(pCovDiff[icent][ipt][id][imeth][i]);
        }
      }
    }
  }

  Double_t v2[nmethod][npid][ncent], v2e[nmethod][npid][ncent];
  Double_t v2_RF[nmethod][ncent],    v2e_RF[nmethod][ncent];
  for (Int_t icent=0; icent<ncent; icent++){ // loop over centrality classes
    MultiparticleCorrelation *multCorrRef[kNMethod] = {nullptr};
    for(Int_t imeth(0);imeth<kNMethod;imeth++) {
      multCorrRef[imeth] = new MultiparticleCorrelation(pCorr[icent][imeth][k2],pCorr[icent][imeth][k4],pCorr[icent][imeth][k6],
                                                        pCorrWSquare[icent][imeth][k2],pCorrWSquare[icent][imeth][k4],pCorrWSquare[icent][imeth][k6],
                                                        pCov[icent][imeth][k24],pCov[icent][imeth][k26],pCov[icent][imeth][k46],zero,"");
    }
    v2_RF[0][icent]  = multCorrRef[0]->GetV22Ref();
    v2e_RF[0][icent] = multCorrRef[0]->GetV22RefErr();
    v2_RF[1][icent]  = multCorrRef[0]->GetV24Ref();
    v2e_RF[1][icent] = multCorrRef[0]->GetV24RefErr();
    v2_RF[2][icent]  = multCorrRef[0]->GetV26Ref();
    v2e_RF[2][icent] = multCorrRef[0]->GetV26RefErr();
    v2_RF[3][icent]  = multCorrRef[1]->GetV22Ref();
    v2e_RF[3][icent] = multCorrRef[1]->GetV22RefErr();
    v2_RF[4][icent]  = multCorrRef[1]->GetV24Ref();
    v2e_RF[4][icent] = multCorrRef[1]->GetV24RefErr();
    v2_RF[5][icent]  = multCorrRef[1]->GetV26Ref();
    v2e_RF[5][icent] = multCorrRef[1]->GetV26RefErr();
    // std::cout << "v22etasub = " << multCorrRef[1]->GetV26Ref() << "\t, v22 = " << multCorrRef[0]->GetV26Ref() << std::endl;
    // std::cout << "v22 = " << subeventQCRef->GetV22Ref() << "\t" << prV2Ref[2]->GetBinContent(icent+1) << "\t v24 = " << standardQCRef->GetV24Ref() << "\t" << prV2Ref[1]->GetBinContent(icent+1) << std::endl;

  } // end of loop over centrality classes
  // Differential flow calculation

  for (Int_t id = 0; id<npid; id++){
    for (Int_t icent=0; icent<ncent; icent++){ // loop over centrality classes
      for(Int_t ipt=binMinPtRFP; ipt<binMinPtRFP+1; ipt++){ // loop for all pT bin
        MultiparticleCorrelation *multCorr[kNMethod];
        for(Int_t imeth(0);imeth<kNMethod;imeth++) {
          multCorr[imeth] = new MultiparticleCorrelation(pCorr[icent][imeth][k2], pCorr[icent][imeth][k4], pCorr[icent][imeth][k6],
          pCorrDiff[icent][ipt][id][imeth][k2], pCorrDiff[icent][ipt][id][imeth][k4], pCorrDiff[icent][ipt][id][imeth][k6],
          pCorrWSquare[icent][imeth][k2], pCorrWSquare[icent][imeth][k4], pCorrWSquare[icent][imeth][k6],
          pCorrDiffWSquare[icent][ipt][id][imeth][k2], pCorrDiffWSquare[icent][ipt][id][imeth][k4], pCorrDiffWSquare[icent][ipt][id][imeth][k6],
          pCovDiff[icent][ipt][id][imeth][k22prime], pCov[icent][imeth][k24], pCovDiff[icent][ipt][id][imeth][k24prime],
          pCov[icent][imeth][k26], pCovDiff[icent][ipt][id][imeth][k26prime],
          pCovDiff[icent][ipt][id][imeth][k2prime4], pCovDiff[icent][ipt][id][imeth][k2prime4prime], pCovDiff[icent][ipt][id][imeth][k2prime6],
          pCovDiff[icent][ipt][id][imeth][k2prime6prime], pCovDiff[icent][ipt][id][imeth][k44prime], pCov[icent][imeth][k46], pCovDiff[icent][ipt][id][imeth][k46prime],
          pCovDiff[icent][ipt][id][imeth][k4prime6], pCovDiff[icent][ipt][id][imeth][k4prime6prime], pCovDiff[icent][ipt][id][imeth][k66prime],zero,zero,"");
        }

        v2[0][id][icent]  = multCorr[0]->GetV22Dif();
        v2e[0][id][icent] = multCorr[0]->GetV22DifErr();
        v2[1][id][icent]  = multCorr[0]->GetV24Dif();
        v2e[1][id][icent] = multCorr[0]->GetV24DifErr();
        v2[2][id][icent]  = multCorr[0]->GetV26Dif();
        v2e[2][id][icent] = multCorr[0]->GetV26DifErr();
        v2[3][id][icent]  = multCorr[1]->GetV22Dif();
        v2e[3][id][icent] = multCorr[1]->GetV22DifErr();
        v2[4][id][icent]  = multCorr[1]->GetV24Dif();
        v2e[4][id][icent] = multCorr[1]->GetV24DifErr();
        v2[5][id][icent]  = multCorr[1]->GetV26Dif();
        v2e[5][id][icent] = multCorr[1]->GetV26DifErr();

      } // end of loop for all pT bin
    } // end of loop over centrality classes
  } // end of loop over PID

  for (Int_t imeth=0; imeth<nmethod; imeth++){
    grRefFl[imeth] = new TGraphErrors(ncent,centrality_bin_center,v2_RF[imeth],bin_centE,v2e_RF[imeth]);
    grRefFl[imeth] -> SetMarkerStyle(marker[imeth]);
    grRefFl[imeth] -> SetMarkerSize(1.5);
    // grRefFl[imeth] -> RemovePoint(0);
    for (Int_t id=0; id<npid; id++){
      grIntFlPID[imeth][id] = new TGraphErrors(ncent,centrality_bin_center,v2[imeth][id],bin_centE,v2e[imeth][id]);
      grIntFlPID[imeth][id] -> SetMarkerStyle(marker[imeth]);
      grIntFlPID[imeth][id] -> SetMarkerSize(1.5);
      // grIntFlPID[imeth][id] -> RemovePoint(0);
    }
  }

  outFile -> cd();
  for (Int_t imeth=0; imeth<nmethod; imeth++){

    for (Int_t id=0;id<npid;id++){
      if (id==8) continue;
      grIntFlPID[imeth][id] -> SetTitle(grTitle[imeth]);
      grIntFlPID[imeth][id] -> GetYaxis()-> SetTitle("v_{2}");
      grIntFlPID[imeth][id] -> GetXaxis()-> SetTitle("Centrality, %");
      grIntFlPID[imeth][id] -> Write(Form("grRF_%i_%i",imeth,id));
    }
  
  }

  std::vector<TGraphErrors*> vgrv2cent[npid];
    for (Int_t id=0;id<npid;id++){
      if (id==8) continue;
      vgrv2cent[id].push_back(grIntFlPID[ratioToMethod][id]); // v2{gapped 2QC}
      for (Int_t imeth=0; imeth<nmethod; imeth++){
        if (imeth==ratioToMethod) continue;
        // if (imeth==0 || imeth==3 || imeth==4 || imeth==5) continue;
        // if (imeth==2 || imeth==4 || imeth==5) continue;
        vgrv2cent[id].push_back(grIntFlPID[imeth][id]);
      }
    }
  
  TCanvas *cV2Cent[npid];
  for (Int_t id=0;id<npid;id++){
    if (id==8) continue;
    cV2Cent[id] = (TCanvas*) DrawTGraph(vgrv2cent[id],"",rangeRatioRF.at(0).first, rangeRatioRF.at(0).second, mincent, maxcent, minV2int, maxV2int,
                                        coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                        level.Data(), Form("%s, %1.1f<p_{T}<%1.1f",pidFancyNames.at(id).Data(),pTBin[binMinPtRFP],pTBin[binMaxPtRFP]),true,grTitle[ratioToMethod]);

    cV2Cent[id] -> SetName(pidFancyNames.at(id).Data());
    if (saveAsPNG) cV2Cent[id] -> SaveAs(Form("./%s/IntegratedFlow_%s.png",outDirName.Data(),pidNames.at(id).Data()));
  }


  TCanvas *cV2CentRF;

  std::vector<TGraphErrors*> vgrv2cent_chargedHardons;
  for (Int_t imeth=0; imeth<nmethod; imeth++){
    grRefFl[imeth] -> SetTitle(grTitle[imeth]);
    grRefFl[imeth] -> GetYaxis()-> SetTitle("v_{2}");
    grRefFl[imeth] -> GetXaxis()-> SetTitle("Centrality, %");
    grRefFl[imeth] -> Write(Form("grRF_%i_8",imeth));
  }

  vgrv2cent_chargedHardons.push_back(grRefFl[ratioToMethod]);
  for (Int_t imeth=0;imeth<nmethod;imeth++){
    if (imeth==ratioToMethod) continue;
    // if (imeth==0 || imeth==3 || imeth==4 || imeth==5) continue;
    // if (imeth==2 || imeth==4 || imeth==5) continue;
    vgrv2cent_chargedHardons.push_back(grRefFl[imeth]);

  }

  cV2CentRF = (TCanvas*) DrawTGraph(vgrv2cent_chargedHardons,"",rangeRatioRF.at(1).first, rangeRatioRF.at(1).second,
                                    mincent, maxcent, minV2int, maxV2int,
                                    coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                    level.Data(), Form("h^{#pm}, %1.1f<p_{T}<%1.1f GeV/c",minptRFP,maxptRFP),true,grTitle[ratioToMethod]);
  cV2CentRF -> SetName("");
  if (saveAsPNG) cV2CentRF -> SaveAs(Form("./%s/IntegratedFlow_hadron.png",outDirName.Data()));

  inFile->Close();
  outFile->Close();
}

Int_t main()
{
  PlotV2QCumulant();
  return 0;
}
