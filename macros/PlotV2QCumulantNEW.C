// #include "MultiparticleCorrelation.cxx"
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include "DrawTGraph.C"
#include "../constants.C"

#include <vector>
#include <iostream>

struct term{ // structure for "Mean squared error of MEAN" calculation, using unbiased estimator for the root of variance
  term(){
    mVal = 0;
    mMSE = 0;
  }
  term(TProfile *const &pr, Int_t bin=0){
    // Eq. (15) http://arxiv.org/abs/2104.00588
    double Neff = pr -> GetBinEffectiveEntries(bin+1);
    mVal = pr -> GetBinContent(bin+1);
    pr -> SetErrorOption("s");
    double stdevW = pr -> GetBinError(bin+1);
    double S2 = stdevW*stdevW/(1-1./Neff);
    mMSE = S2/Neff;

    // Eq. (17) http://arxiv.org/abs/2104.00588
  };
  term(TProfile *const &pr1, TProfile *const &pr2, Int_t bin=0)
  {
    double Neff = pr1 -> GetBinEffectiveEntries(bin+1);
    // double N = pr1 -> GetEntries();
    mVal = pr1 -> GetBinContent(bin+1);
    pr2 -> SetErrorOption("s");
    double stdevW = pr2 -> GetBinError(bin+1);
    // mMSE = N/(N-1.)*stdevW*stdevW/Neff;
    mMSE = stdevW*stdevW/Neff;
  };
 public: 
  double mVal; // weithted mean value
  double mMSE; // Mean squared error of mean, https://en.wikipedia.org/wiki/Mean_squared_error

};

double Covariance(TProfile *const &hcovXY, TProfile *const &hX, TProfile *const &hY, Int_t binXY=0, Int_t binX=0, Int_t binY=0)
{
  double mSumWXY = hcovXY->GetBinEntries(binXY+1);
  double sumWX = hX->GetBinEntries(binX+1);
  double sumWY = hY->GetBinEntries(binY+1);

  double meanXY = hcovXY -> GetBinContent(binXY+1);
  double meanX = hX -> GetBinContent(binX+1);
  double meanY = hY -> GetBinContent(binY+1);
  double mVal = (meanXY-meanX*meanY)/(sumWX*sumWY/mSumWXY-1.); // Cov(x,y)/(sumWX*sumWY/sumWXY)
  return mVal;  
}

enum CovTerm {k24, k26, k46, kNCov};
enum DiffCovTerm {k22prime, k24prime, k26prime, k2prime4, k2prime4prime, k2prime6, k2prime6prime, k44prime, k46prime, k4prime6, k4prime6prime, k66prime, kNDiffCov};
enum Corr {k2, k4, k6, kNCorr};
enum DiffCorr {k2prime, k4prime, k6prime, kNDiffCorr};
enum Method {kStd, kGap, kNMethod};
// Flags
Bool_t saveAsPNG = true;
Int_t ratioToMethod = 3; // 2QC, 4QC, 6QC, 2QC-gapped, 4QC-gapped, 6QC-gapped
Int_t drawDifferentialFlowTill = 0; // Draw v2 vs pT (10% centrality cut) till: 0: no drawing; 1: till 10%; 2: till 20%; etc.
const std::vector<TString> pidFancyNames = {"h^{+}", "#pi^{+}", "K^{+}", "p", "h^{-}", "#pi^{-}", "K^{-}", "#bar{p}", "h^{#pm}","#pi^{#pm}","K^{#pm}","p(#bar{p})"};
const float eta_gap = 0.05;
TString model = "UrQMD";
TString modelFancy = "UrQMD";
TString energy = "4.5GeV";
TString inputFileName = Form("SecondRun_UrQMD_4GeV.root");
TFile *outFile = new TFile(Form("graphs_v2_QC_UrQMD_4GeV.root"),"recreate");

const Int_t nmethod = 6; // 2QC, 4QC, 6QC, 2QC-gapped, 4QC-gapped, 6QC-gapped
const Int_t binMinPtRFP = 0;      // 0.2 GeV
const Int_t binMaxPtRFP = npt-1;  // 3.0 GeV

const double minptRFP = 0.2;
const double maxptRFP = 3.0;

const double maxpt = 3.0; // for v2 vs pt plotting
const double minpt = 0.0;  // for v2 vs pt plotting

// const Int_t ncent = 19;
// const Double_t bin_cent[ncent + 1] = {0, 0.5, 1, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
                                      // 0,   1, 2, 3, 4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19

// const double wideCentBin[4] = {0,5,11,19};
const int bin10 = 5;
const int bin40 = 11;
const double wideCentBin[2] = {bin10,bin40};

const double mincent = 0.;  // for v2 vs centrality
const double maxcent = 80.; // for v2 vs centrality

const double minV2int = -0.005; // for v2 vs centrality plotting
const double maxV2int = 0.1; // for v2 vs centrality plotting
const double minV2dif = -0.01; // for v2 vs pt plotting
const double maxV2dif = 0.25; // for v2 vs pt plotting


std::vector <double> coordinateLeg = {0.18,0.63,0.45,0.889};
std::vector<std::pair<double,double>> rangeRatioRF ={{0.54,1.13},{0.54,1.13}};
std::vector<std::pair<double,double>> rangeRatio = {{0.67,1.13},{0.67,1.13}};
Int_t marker[nmethod]={kFullSquare,kFullCircle,kFullTriangleDown,kOpenSquare,kOpenCircle,kOpenTriangleDown}; // 2QC, 4QC, 6QC, 2QC-gapped, 4QC-gapped, 6QC-gapped

void CalStatErrCent1040(double v2eDif1040[nmethod][npid][npt]){

  TFile *inFile = new TFile(inputFileName.Data(),"read");

  TProfile* pCorr[kNCorr];
  TProfile* pCorrWSquare[kNCorr];
  TProfile* pCov[kNCov];
  TProfile* pCovDiff[kNDiffCov];
  TProfile* pCorrDiff[kNDiffCorr];
  TProfile* pCorrDiffWSquare[kNDiffCorr];
  // TProfile* tmp;
  for(Int_t imeth(0);imeth<kNMethod;imeth++) {
    for (Int_t i(0);i<kNCorr;i++) {
      pCorr[i]        = (TProfile*) inFile->Get(Form("QC/pCorr_%i_%i",imeth,i));
      pCorrWSquare[i] = (TProfile*) inFile->Get(Form("QC/pCorrWSquare_%i_%i",imeth,i));
      pCorr[i]        = (TProfile*) pCorr[i]       ->Rebin(1,"",&wideCentBin[0]);
      pCorrWSquare[i] = (TProfile*) pCorrWSquare[i]->Rebin(1,"",&wideCentBin[0]);
    }
    for (Int_t i(0);i<kNCov;i++) {
      pCov[i] = (TProfile*) inFile->Get(Form("QC/pCov_%i_%i",imeth,i));
      pCov[i] = (TProfile*) pCov[i]->Rebin(1,"",&wideCentBin[0]);
    }

    term cor2 = term(pCorr[k2],pCorrWSquare[k2]);
    double v22 = sqrt(cor2.mVal);
    double ev22 = sqrt(1./(4.*cor2.mVal)*cor2.mMSE);
    term cor4 = term(pCorr[k4],pCorrWSquare[k4]);
    double cov24 = Covariance(pCov[k24],pCorr[k2],pCorr[k4]);
    double v24 = pow(2*pow(cor2.mVal,2)-cor4.mVal,0.25);
    double ev24 = sqrt( 1./pow(v24,6)*(cor2.mVal*cor2.mVal*cor2.mMSE+1./16*cor4.mMSE-0.5*cor2.mVal*cov24) );
    term cor6 = term(pCorr[k6],pCorrWSquare[k6]);
    double cov26 = Covariance(pCov[k26],pCorr[k2],pCorr[k6]);
    double cov46 = Covariance(pCov[k46],pCorr[k4],pCorr[k6]);
    double v26 = pow(2,-1./3.)*pow(cor6.mVal-9*cor2.mVal*cor4.mVal+12*pow(cor2.mVal,3),1./6.);
    double ev26 = sqrt(1./(2*pow(2,4)*pow(v26,10))*(9./2.*pow(4*pow(cor2.mVal,2)-cor4.mVal,2)*cor2.mMSE
                      + 9./2.*pow(cor2.mVal,2)*cor4.mMSE+1./18.*cor6.mMSE
                      - 9*cor2.mVal*(4*pow(cor2.mVal,2)-cor4.mVal)*cov24
                      +(4*pow(cor2.mVal,2)-cor4.mVal)*cov26
                      - cor2.mVal*cov46));
    for (Int_t id=0;id<npid;id++){
      for(Int_t ipt=0; ipt<npt; ipt++){
        for (Int_t i(0);i<kNDiffCorr;i++) {
          pCorrDiff[i]         = (TProfile*) inFile->Get(Form("QC/pCorrDiff_%i_%i_%i_%i",ipt,id,imeth,i));
          pCorrDiffWSquare[i]  = (TProfile*) inFile->Get(Form("QC/pCorrDiffWSquare_%i_%i_%i_%i",ipt,id,imeth,i));
          pCorrDiff[i]         = (TProfile*) pCorrDiff[i]       ->Rebin(1,"",&wideCentBin[0]);
          pCorrDiffWSquare[i]  = (TProfile*) pCorrDiffWSquare[i]->Rebin(1,"",&wideCentBin[0]);
        }
        for (Int_t i(0);i<kNDiffCov;i++) {
          pCovDiff[i] = (TProfile*) inFile->Get(Form("QC/pCovDiff_%i_%i_%i_%i",ipt,id,imeth,i));
          pCovDiff[i] = (TProfile*) pCovDiff[i]->Rebin(1,"",&wideCentBin[0]);
        }        
        term cor2red = term(pCorrDiff[k2],pCorrDiffWSquare[k2]);
        double v22Dif = cor2red.mVal/v22;
        double cov22prime = Covariance(pCovDiff[k22prime],pCorr[k2],pCorrDiff[k2]);
        double ev22Dif = sqrt(0.25*pow(cor2.mVal,-3)*(pow(cor2red.mVal,2)*cor2.mMSE
                            + 4*pow(cor2.mVal,2)*cor2red.mMSE - 4*cor2.mVal*cor2red.mVal*cov22prime));
        term cor4red = term(pCorrDiff[k4],pCorrDiffWSquare[k4]);
        double cov24prime = Covariance(pCovDiff[k24prime],pCorr[k2],pCorrDiff[k4]);
        double cov42prime = Covariance(pCovDiff[k24prime],pCorr[k4],pCorrDiff[k2]);
        double cov44prime = Covariance(pCovDiff[k44prime],pCorr[k4],pCorrDiff[k4]);
        double cov2prime4prime = Covariance(pCovDiff[k2prime4prime],pCorrDiff[k2],pCorrDiff[k4]);
        double v24Dif = (2.*cor2.mVal*cor2red.mVal-cor4red.mVal)*pow(v24,-3);
        double ev24Dif = sqrt( pow(v24,-14)
            * (pow(2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal,2.)
            * cor2.mMSE
            + 9./16*pow(2.*cor2.mVal*cor2red.mVal-cor4red.mVal,2.)*cor4.mMSE
            + 4*pow(cor2.mVal,2)*pow(v24,8)*cor2red.mMSE
            + pow(v24,8)*cor4red.mMSE
            - 1.5*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
            * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
            * cov24
            - 4*cor2.mVal*pow(v24,4)
            * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
            * cov22prime
            + 2*pow(v24,4)
            * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
            * cov24prime
            + 3*cor2.mVal*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
            * cov42prime
            - 1.5*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
            * cov44prime
            - 4*cor2.mVal*pow(v24,8)*cov2prime4prime));
        term cor6red = term(pCorrDiff[k6],pCorrDiffWSquare[k6]);
        double cov26prime = Covariance(pCovDiff[k26prime],pCorr[k2],pCorrDiff[k6]);
        double cov2prime6 = Covariance(pCovDiff[k2prime6],pCorrDiff[k2],pCorr[k6]);
        double cov2prime6prime = Covariance(pCovDiff[k2prime6prime],pCorrDiff[k2],pCorrDiff[k6]);
        double cov46prime = Covariance(pCovDiff[k46prime],pCorr[k4],pCorrDiff[k6]);
        double cov4prime6 = Covariance(pCovDiff[k4prime6],pCorrDiff[k4],pCorr[k6]);
        double cov4prime6prime = Covariance(pCovDiff[k4prime6prime],pCorrDiff[k4],pCorrDiff[k6]);
        double cov66prime = Covariance(pCovDiff[k66prime],pCorr[k6],pCorrDiff[k6]);
        double dNumeratorDiffV26 = cor6red.mVal-9*cor2.mVal*cor4red.mVal+12*cor2red.mVal*pow(cor2.mVal,2);
        double v26Dif = dNumeratorDiffV26 / (4.*pow(v26,5));
        double dFourTimesV26RefToTheSixthPower = 4.*pow(v26,6.);
        double dPartialDerivativeOfDiffV26WRT2Corr    = (5.*pow(4.,1./3.)*(9.*cor4.mVal-36.*pow(cor2.mVal,2.))*dNumeratorDiffV26/(12.*pow(dFourTimesV26RefToTheSixthPower,11./6.)))
                                                        - (pow(4.,1./3.)*(9.*cor4red.mVal-24.*cor2.mVal*cor2red.mVal)/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.)));
        double dPartialDerivativeOfDiffV26WRT2CorrRed = 6.*pow(4.,1./3.)*pow(cor2.mVal,2.)/pow(dFourTimesV26RefToTheSixthPower,5./6.);
        double dPartialDerivativeOfDiffV26WRT4Corr    = 15.*pow(4.,1./3.)*cor2.mVal*dNumeratorDiffV26/(4.*pow(dFourTimesV26RefToTheSixthPower,11./6.));
        double dPartialDerivativeOfDiffV26WRT4CorrRed = (-1.)*9.*pow(4.,1./3.)*cor2.mVal/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.));
        double dPartialDerivativeOfDiffV26WRT6Corr    = (-1.)*5.*pow(4.,1./3.)*dNumeratorDiffV26/(12.*pow(dFourTimesV26RefToTheSixthPower,11./6.));
        double dPartialDerivativeOfDiffV26WRT6CorrRed = pow(4.,1./3.)/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.));
        double ev26Dif = sqrt(
          pow(dPartialDerivativeOfDiffV26WRT2Corr,    2.)*cor2.mMSE
        + pow(dPartialDerivativeOfDiffV26WRT2CorrRed, 2.)*cor2red.mMSE
        + pow(dPartialDerivativeOfDiffV26WRT4Corr,    2.)*cor4.mMSE
        + pow(dPartialDerivativeOfDiffV26WRT4CorrRed, 2.)*cor4red.mMSE
        + pow(dPartialDerivativeOfDiffV26WRT6Corr,    2.)*cor6.mMSE
        + pow(dPartialDerivativeOfDiffV26WRT6CorrRed, 2.)*cor6red.mMSE
        + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT2CorrRed) * cov22prime
        + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT4Corr)    * cov24
        + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov24prime
        + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT6Corr)    * cov26
        + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov26prime
        + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT4Corr   ) * cov42prime
        + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov2prime4prime
        + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT6Corr   ) * cov2prime6
        + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov2prime6prime
        + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov44prime
        + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT6Corr   ) * cov46
        + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov46prime
        + 2.*(dPartialDerivativeOfDiffV26WRT4CorrRed * dPartialDerivativeOfDiffV26WRT6Corr)     * cov4prime6
        + 2.*(dPartialDerivativeOfDiffV26WRT4CorrRed * dPartialDerivativeOfDiffV26WRT6CorrRed)  * cov4prime6prime
        + 2.*(dPartialDerivativeOfDiffV26WRT6Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov66prime
        );
        if (imeth == kStd) {
          v2eDif1040[0][id][ipt] = ev22Dif;
          v2eDif1040[1][id][ipt] = ev24Dif;
          v2eDif1040[2][id][ipt] = ev26Dif;
        }
        if (imeth == kGap) {
          v2eDif1040[3][id][ipt] = ev22Dif;
          v2eDif1040[4][id][ipt] = ev24Dif;
          v2eDif1040[5][id][ipt] = ev26Dif;
        }
      } // end of loop for all pT bin
    } // end of loop for PID
  } // for(Int_t imeth(0);imeth<kNMethod;imeth++)
  inFile->Close();
}

void PlotV2QCumulantNEW(){
  TString outDirName=(TString)Form("%s_%s_eta_gap_%1.1f",model.Data(),energy.Data(),eta_gap*2);
  TString level= (TString) Form("%s, Au+Au at #sqrt{s_{NN}}=%s",modelFancy.Data(),energy.Data());
  double v2eDif1040[nmethod][npid][npt];
  CalStatErrCent1040(v2eDif1040);
  TFile *inFile = new TFile(inputFileName.Data(),"read");
  // Temporary variables
  char hname[800]; // histogram hname
  TProfile* pCorr[kNCorr];
  TProfile* pCorrWSquare[kNCorr];
  TProfile* pCov[kNCov];
  TProfile* pCovDiff[npt][npid][kNDiffCov];
  TProfile* pCorrDiff[npt][npid][kNDiffCorr];
  TProfile* pCorrDiffWSquare[npt][npid][kNDiffCorr];
  TProfile *pMult[npt][npid];

  // OUTPUT
  TGraphErrors *grDifFl[nmethod][ncent][npid];    // v2(pt); 3 = {2QC, 4QC, EP, gapped 2QC}
  TGraphErrors *grDifFl1040[nmethod][npid];

  // Filling pT bin
  double pt[npt];
  const double ept[npt]={0}; // error bin pT = 0.0
  for (Int_t ipt=0; ipt<npt; ipt++){
    // pt[icent][ipt][id] = hmult[icent][ipt][id] -> GetBinContent(1);
    pt[ipt] = ( pTBin[ipt] + pTBin[ipt+1] ) / 2.;
  }
  double cent[ncent];
  const double centErr[ncent] = {0};
  for (Int_t icent=0; icent<ncent; icent++){
    cent[icent] = ( bin_cent[icent] + bin_cent[icent+1] ) / 2.;
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
  double v2Dif[nmethod][ncent][npid][npt], v2eDif[nmethod][ncent][npid][npt];
  for(Int_t imeth(0);imeth<kNMethod;imeth++) {
    for (Int_t i(0);i<kNCorr;i++) {
      pCorr[i]        = (TProfile*) inFile->Get(Form("QC/pCorr_%i_%i",imeth,i));
      pCorrWSquare[i] = (TProfile*) inFile->Get(Form("QC/pCorrWSquare_%i_%i",imeth,i));
    }
    for (Int_t i(0);i<kNCov;i++) {
      pCov[i] = (TProfile*) inFile->Get(Form("QC/pCov_%i_%i",imeth,i));
    }
    for (Int_t id=0;id<npid;id++) {
      for(Int_t ipt=0; ipt<npt; ipt++) {
        pMult[ipt][id] = (TProfile*) inFile->Get(Form("QC/pMult_%i_%i_%i",ipt,id,imeth));
        for (Int_t i(0);i<kNDiffCorr;i++) {
          pCorrDiff[ipt][id][i]         = (TProfile*) inFile->Get(Form("QC/pCorrDiff_%i_%i_%i_%i",ipt,id,imeth,i));
          pCorrDiffWSquare[ipt][id][i]  = (TProfile*) inFile->Get(Form("QC/pCorrDiffWSquare_%i_%i_%i_%i",ipt,id,imeth,i));
        }
        for (Int_t i(0);i<kNDiffCov;i++) {
          pCovDiff[ipt][id][i] = (TProfile*) inFile->Get(Form("QC/pCovDiff_%i_%i_%i_%i",ipt,id,imeth,i));
        }
      } // for(Int_t ipt=0; ipt<npt; ipt++)
    } // for (Int_t id=0;id<npid;id++)

    for (Int_t icent=0; icent<ncent; icent++) {
      term cor2 = term(pCorr[k2],pCorrWSquare[k2],icent);
      double v22 = sqrt(cor2.mVal);
      double ev22 = sqrt(1./(4.*cor2.mVal)*cor2.mMSE);
      term cor4 = term(pCorr[k4],pCorrWSquare[k4],icent);
      double cov24 = Covariance(pCov[k24],pCorr[k2],pCorr[k4],icent,icent,icent);
      double v24 = pow(2*pow(cor2.mVal,2)-cor4.mVal,0.25);
      double ev24 = sqrt( 1./pow(v24,6)*(cor2.mVal*cor2.mVal*cor2.mMSE+1./16*cor4.mMSE-0.5*cor2.mVal*cov24) );
      term cor6 = term(pCorr[k6],pCorrWSquare[k6],icent);
      double cov26 = Covariance(pCov[k26],pCorr[k2],pCorr[k6],icent,icent,icent);
      double cov46 = Covariance(pCov[k46],pCorr[k4],pCorr[k6],icent,icent,icent);
      double v26 = pow(2,-1./3.)*pow(cor6.mVal-9*cor2.mVal*cor4.mVal+12*pow(cor2.mVal,3),1./6.);
      double ev26 = sqrt(1./(2*pow(2,4)*pow(v26,10))*(9./2.*pow(4*pow(cor2.mVal,2)-cor4.mVal,2)*cor2.mMSE
                        + 9./2.*pow(cor2.mVal,2)*cor4.mMSE+1./18.*cor6.mMSE
                        - 9*cor2.mVal*(4*pow(cor2.mVal,2)-cor4.mVal)*cov24
                        +(4*pow(cor2.mVal,2)-cor4.mVal)*cov26
                        - cor2.mVal*cov46));
      for (Int_t id=0;id<npid;id++) {
        for(Int_t ipt=0; ipt<npt; ipt++) {
          term cor2red = term(pCorrDiff[ipt][id][k2],pCorrDiffWSquare[ipt][id][k2],icent);
          double v22Dif = cor2red.mVal/v22;
          double cov22prime = Covariance(pCovDiff[ipt][id][k22prime],pCorr[k2],pCorrDiff[ipt][id][k2],icent,icent,icent);
          double ev22Dif = sqrt(0.25*pow(cor2.mVal,-3)*(pow(cor2red.mVal,2)*cor2.mMSE
                              + 4*pow(cor2.mVal,2)*cor2red.mMSE - 4*cor2.mVal*cor2red.mVal*cov22prime));

          term cor4red = term(pCorrDiff[ipt][id][k4],pCorrDiffWSquare[ipt][id][k4],icent);
          double cov24prime = Covariance(pCovDiff[ipt][id][k24prime],pCorr[k2],pCorrDiff[ipt][id][k4],icent,icent,icent);
          double cov42prime = Covariance(pCovDiff[ipt][id][k24prime],pCorr[k4],pCorrDiff[ipt][id][k2],icent,icent,icent);
          double cov44prime = Covariance(pCovDiff[ipt][id][k44prime],pCorr[k4],pCorrDiff[ipt][id][k4],icent,icent,icent);
          double cov2prime4prime = Covariance(pCovDiff[ipt][id][k2prime4prime],pCorrDiff[ipt][id][k2],pCorrDiff[ipt][id][k4],icent,icent,icent);
          double v24Dif = (2.*cor2.mVal*cor2red.mVal-cor4red.mVal)*pow(v24,-3);
          double ev24Dif = sqrt( pow(v24,-14)
              * (pow(2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal,2.)
              * cor2.mMSE
              + 9./16*pow(2.*cor2.mVal*cor2red.mVal-cor4red.mVal,2.)*cor4.mMSE
              + 4*pow(cor2.mVal,2)*pow(v24,8)*cor2red.mMSE
              + pow(v24,8)*cor4red.mMSE
              - 1.5*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
              * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
              * cov24
              - 4*cor2.mVal*pow(v24,4)
              * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
              * cov22prime
              + 2*pow(v24,4)
              * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
              * cov24prime
              + 3*cor2.mVal*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
              * cov42prime
              - 1.5*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
              * cov44prime
              - 4*cor2.mVal*pow(v24,8)*cov2prime4prime));

          term cor6red = term(pCorrDiff[ipt][id][k6],pCorrDiffWSquare[ipt][id][k6],icent);
          double cov26prime = Covariance(pCovDiff[ipt][id][k26prime],pCorr[k2],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov2prime6 = Covariance(pCovDiff[ipt][id][k2prime6],pCorrDiff[ipt][id][k2],pCorr[k6],icent,icent,icent);
          double cov2prime6prime = Covariance(pCovDiff[ipt][id][k2prime6prime],pCorrDiff[ipt][id][k2],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov46prime = Covariance(pCovDiff[ipt][id][k46prime],pCorr[k4],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov4prime6 = Covariance(pCovDiff[ipt][id][k4prime6],pCorrDiff[ipt][id][k4],pCorr[k6],icent,icent,icent);
          double cov4prime6prime = Covariance(pCovDiff[ipt][id][k4prime6prime],pCorrDiff[ipt][id][k4],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov66prime = Covariance(pCovDiff[ipt][id][k66prime],pCorr[k6],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double dNumeratorDiffV26 = cor6red.mVal-9*cor2.mVal*cor4red.mVal+12*cor2red.mVal*pow(cor2.mVal,2);
          double v26Dif = dNumeratorDiffV26 / (4.*pow(v26,5));
          double dFourTimesV26RefToTheSixthPower = 4.*pow(v26,6.);
          double dPartialDerivativeOfDiffV26WRT2Corr    = (5.*pow(4.,1./3.)*(9.*cor4.mVal-36.*pow(cor2.mVal,2.))*dNumeratorDiffV26/(12.*pow(dFourTimesV26RefToTheSixthPower,11./6.)))
                                                          - (pow(4.,1./3.)*(9.*cor4red.mVal-24.*cor2.mVal*cor2red.mVal)/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.)));
          double dPartialDerivativeOfDiffV26WRT2CorrRed = 6.*pow(4.,1./3.)*pow(cor2.mVal,2.)/pow(dFourTimesV26RefToTheSixthPower,5./6.);
          double dPartialDerivativeOfDiffV26WRT4Corr    = 15.*pow(4.,1./3.)*cor2.mVal*dNumeratorDiffV26/(4.*pow(dFourTimesV26RefToTheSixthPower,11./6.));
          double dPartialDerivativeOfDiffV26WRT4CorrRed = (-1.)*9.*pow(4.,1./3.)*cor2.mVal/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.));
          double dPartialDerivativeOfDiffV26WRT6Corr    = (-1.)*5.*pow(4.,1./3.)*dNumeratorDiffV26/(12.*pow(dFourTimesV26RefToTheSixthPower,11./6.));
          double dPartialDerivativeOfDiffV26WRT6CorrRed = pow(4.,1./3.)/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.));
          double ev26Dif = sqrt(
            pow(dPartialDerivativeOfDiffV26WRT2Corr,    2.)*cor2.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT2CorrRed, 2.)*cor2red.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT4Corr,    2.)*cor4.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT4CorrRed, 2.)*cor4red.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT6Corr,    2.)*cor6.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT6CorrRed, 2.)*cor6red.mMSE
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT2CorrRed) * cov22prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT4Corr)    * cov24
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov24prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT6Corr)    * cov26
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov26prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT4Corr   ) * cov42prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov2prime4prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT6Corr   ) * cov2prime6
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov2prime6prime
          + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov44prime
          + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT6Corr   ) * cov46
          + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov46prime
          + 2.*(dPartialDerivativeOfDiffV26WRT4CorrRed * dPartialDerivativeOfDiffV26WRT6Corr)     * cov4prime6
          + 2.*(dPartialDerivativeOfDiffV26WRT4CorrRed * dPartialDerivativeOfDiffV26WRT6CorrRed)  * cov4prime6prime
          + 2.*(dPartialDerivativeOfDiffV26WRT6Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov66prime
          );
          if (imeth == kStd) {
            v2Dif[0][icent][id][ipt] = v22Dif;
            v2eDif[0][icent][id][ipt] = ev22Dif;
            v2Dif[1][icent][id][ipt] = v24Dif;
            v2eDif[1][icent][id][ipt] = ev24Dif;      
            v2Dif[2][icent][id][ipt] = v26Dif;
            v2eDif[2][icent][id][ipt] = ev26Dif;
          }
          if (imeth == kGap) {
            v2Dif[3][icent][id][ipt] = v22Dif;
            v2eDif[3][icent][id][ipt] = ev22Dif;
            v2Dif[4][icent][id][ipt] = v24Dif;
            v2eDif[4][icent][id][ipt] = ev24Dif;      
            v2Dif[5][icent][id][ipt] = v26Dif;
            v2eDif[5][icent][id][ipt] = ev26Dif;
          }


          if (icent>=bin10 && icent <bin40) { // 10-40%
            if (imeth == kStd) {
            prV2Dif1040[0][id] -> Fill(0.5+ipt,v22Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            prV2Dif1040[1][id] -> Fill(0.5+ipt,v24Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            prV2Dif1040[2][id] -> Fill(0.5+ipt,v26Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            } else {
            prV2Dif1040[3][id] -> Fill(0.5+ipt,v22Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            prV2Dif1040[4][id] -> Fill(0.5+ipt,v24Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            prV2Dif1040[5][id] -> Fill(0.5+ipt,v26Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            }
          }
          if (id == 8){
            if (imeth == kStd) {
              prV2Ref[0]->Fill(icent,v22Dif,pMult[ipt][id]->GetBinEntries(icent+1));
              prV2Ref[1]->Fill(icent,v24Dif,pMult[ipt][id]->GetBinEntries(icent+1));
              prV2Ref[2]->Fill(icent,v26Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            } else {
              prV2Ref[3]->Fill(icent,v22Dif,pMult[ipt][id]->GetBinEntries(icent+1));
              prV2Ref[4]->Fill(icent,v24Dif,pMult[ipt][id]->GetBinEntries(icent+1));
              prV2Ref[5]->Fill(icent,v26Dif,pMult[ipt][id]->GetBinEntries(icent+1));
            }
          }

        } // for(Int_t ipt=0; ipt<npt; ipt++)

        for (Int_t i=0; i<nmethod; i++) {
          grDifFl[i][icent][id] = new TGraphErrors(npt,pt,v2Dif[i][icent][id],ept,v2eDif[i][icent][id]);
          grDifFl[i][icent][id] -> SetMarkerStyle(marker[i]);
          grDifFl[i][icent][id] -> SetMarkerSize(1.5);
          grDifFl[i][icent][id] -> SetDrawOption("P");
        }
      } // end of loop for PID
    } // for (Int_t icent=0; icent<ncent; icent++)
  } // for(Int_t imeth(0);imeth<kNMethod;imeth++)

  double v2Dif1040[nmethod][npid][npt];
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
      cV2PT[icent][id] = (TCanvas*) DrawTGraph(vgrv2pt,"",rangeRatio.at(0).first, rangeRatio.at(0).second, minpt, maxpt, minV2dif, maxV2dif,
                                               coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                               level.Data(), Form("%s, %.0f-%.0f%%",pidFancyNames.at(id).Data(),bin_cent[icent],bin_cent[icent+1]), true, grTitle[ratioToMethod]);
      cV2PT[icent][id] -> SetName("");
      if (saveAsPNG) cV2PT[icent][id] -> SaveAs(Form("./%s/DifferentialFlow_Centrality%.0f-%.0f_%s.png",outDirName.Data(),bin_cent[icent],bin_cent[icent+1],pidNames.at(id).Data()));
    }
  }

  // TCanvas *cV2PTMultPad[npid];
  // TString strCent[5];
  // for (Int_t id=0;id<npid;id++)
  // {
  //   std::vector<TGraphErrors*> vgrv2pt[5];
  //   for (Int_t icent=1; icent<6; icent++)
  //   {
  //     vgrv2pt[icent-1].push_back(grDifFl[ratioToMethod][icent][id]);
  //     for (Int_t imeth=0; imeth<nmethod; imeth++){
  //       if (imeth==ratioToMethod) continue;
  //       vgrv2pt[icent-1].push_back(grDifFl[imeth][icent][id]);
  //     }
  //   }
  //   cV2PTMultPad[id] = (TCanvas*) DrawTGraph(vgrv2pt, 5,"",rangeRatio.at(0).first, rangeRatio.at(0).second, minpt, maxpt, minV2dif, maxV2dif,
  //                                             coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
  //                                             Form("%s, %s", level.Data(), pidFancyNames.at(id).Data()),
  //                                             strCent, true, grTitle[ratioToMethod]);
  //   cV2PTMultPad[id] -> SetName("");
  //   if (saveAsPNG) cV2PTMultPad[id] -> SaveAs(Form("./%s/DifferentialFlow_%s_Cent_0_40.png",outDirName.Data(),pidNames.at(id).Data()));
  // }

  TCanvas *cV2PT1040[npid];
  for (Int_t id=0;id<npid;id++){
    std::vector<TGraphErrors*> vgrv2pt1040;
    vgrv2pt1040.push_back(grDifFl1040[ratioToMethod][id]);
    for (Int_t imeth=0;imeth<nmethod;imeth++){
      if (imeth==ratioToMethod) continue;
      vgrv2pt1040.push_back(grDifFl1040[imeth][id]);
    }
    cV2PT1040[id] = (TCanvas*) DrawTGraph(vgrv2pt1040,"",rangeRatio.at(1).first, rangeRatio.at(1).second, minpt, 3.0, minV2dif, maxV2dif,
                                          coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                          level.Data(), Form("10-40%%, %s",pidFancyNames.at(id).Data()), true, grTitle[ratioToMethod]);
    cV2PT1040[id]->SetName("");
    if (saveAsPNG) cV2PT1040[id] -> SaveAs(Form("./%s/DifferentialFlow_Centrality10-40%%_%s.png",outDirName.Data(),pidNames.at(id).Data()));
  }
  //==========================================================================================================================

  TGraphErrors *grIntFlPID[nmethod][npid];
  TGraphErrors *grRefFl[nmethod];
  // v2 vs centrality for PID

  double v2[nmethod][npid][ncent];
  double v2e[nmethod][npid][ncent];
  double v2_RF[nmethod][ncent];
  double v2e_RF[nmethod][ncent];

  for(Int_t imeth(0);imeth<kNMethod;imeth++) {
    for (Int_t i(0);i<kNCorr;i++) {
      pCorr[i]        = (TProfile*) inFile->Get(Form("QC/pCorr_%i_%i",imeth,i));
      pCorrWSquare[i] = (TProfile*) inFile->Get(Form("QC/pCorrWSquare_%i_%i",imeth,i));
    }
    for (Int_t i(0);i<kNCov;i++) {
      pCov[i] = (TProfile*) inFile->Get(Form("QC/pCov_%i_%i",imeth,i));
    }
    for (Int_t id=0;id<npid;id++) {
      for(Int_t ipt=0; ipt<npt; ipt++) {
        for (Int_t i(0);i<kNDiffCorr;i++) {
          pCorrDiff[ipt][id][i]         = (TProfile*) inFile->Get(Form("QC/pCorrDiff_%i_%i_%i_%i",ipt,id,imeth,i));
          pCorrDiffWSquare[ipt][id][i]  = (TProfile*) inFile->Get(Form("QC/pCorrDiffWSquare_%i_%i_%i_%i",ipt,id,imeth,i));
        }
        for (Int_t i(0);i<kNDiffCov;i++) {
          pCovDiff[ipt][id][i] = (TProfile*) inFile->Get(Form("QC/pCovDiff_%i_%i_%i_%i",ipt,id,imeth,i));
        }
      } // for(Int_t ipt=0; ipt<npt; ipt++)
    } // for (Int_t id=0;id<npid;id++)
    for (Int_t icent=0; icent<ncent; icent++) {
      term cor2 = term(pCorr[k2],pCorrWSquare[k2],icent);
      double v22 = sqrt(cor2.mVal);
      double ev22 = sqrt(1./(4.*cor2.mVal)*cor2.mMSE);
      term cor4 = term(pCorr[k4],pCorrWSquare[k4],icent);
      double cov24 = Covariance(pCov[k24],pCorr[k2],pCorr[k4],icent,icent,icent);
      double v24 = pow(2*pow(cor2.mVal,2)-cor4.mVal,0.25);
      double ev24 = sqrt( 1./pow(v24,6)*(cor2.mVal*cor2.mVal*cor2.mMSE+1./16*cor4.mMSE-0.5*cor2.mVal*cov24) );
      term cor6 = term(pCorr[k6],pCorrWSquare[k6],icent);
      double cov26 = Covariance(pCov[k26],pCorr[k2],pCorr[k6],icent,icent,icent);
      double cov46 = Covariance(pCov[k46],pCorr[k4],pCorr[k6],icent,icent,icent);
      double v26 = pow(2,-1./3.)*pow(cor6.mVal-9*cor2.mVal*cor4.mVal+12*pow(cor2.mVal,3),1./6.);
      double ev26 = sqrt(1./(2*pow(2,4)*pow(v26,10))*(9./2.*pow(4*pow(cor2.mVal,2)-cor4.mVal,2)*cor2.mMSE
                        + 9./2.*pow(cor2.mVal,2)*cor4.mMSE+1./18.*cor6.mMSE
                        - 9*cor2.mVal*(4*pow(cor2.mVal,2)-cor4.mVal)*cov24
                        +(4*pow(cor2.mVal,2)-cor4.mVal)*cov26
                        - cor2.mVal*cov46));
      if (imeth == kStd) {
        v2_RF[0][icent]  = v22;
        v2e_RF[0][icent] = ev22;
        v2_RF[1][icent]  = v24;
        v2e_RF[1][icent] = ev24;
        v2_RF[2][icent]  = v26;
        v2e_RF[2][icent] = ev26;
      } else {
        v2_RF[3][icent]  = v22;
        v2e_RF[3][icent] = ev22;
        v2_RF[4][icent]  = v24;
        v2e_RF[4][icent] = ev24;
        v2_RF[5][icent]  = v26;
        v2e_RF[5][icent] = ev26;
      }
      // Differential flow calculation
      for (Int_t id = 0; id<npid; id++) {
        for(Int_t ipt=binMinPtRFP+1; ipt<binMaxPtRFP; ipt++) {
          for (Int_t i(0);i<kNDiffCorr;i++) {
            pCorrDiff[binMinPtRFP][id][i]->Add(pCorrDiff[ipt][id][i]);
            pCorrDiffWSquare[binMinPtRFP][id][i]->Add(pCorrDiffWSquare[ipt][id][i]);
          }
          for (Int_t i(0);i<kNDiffCov;i++) {
            pCovDiff[binMinPtRFP][id][i]->Add(pCovDiff[ipt][id][i]);
          }
        } // for(Int_t ipt=binMinPtRFP+1; ipt<binMaxPtRFP; ipt++)
        for(Int_t ipt=binMinPtRFP; ipt<binMinPtRFP+1; ipt++) {
          term cor2red = term(pCorrDiff[ipt][id][k2],pCorrDiffWSquare[ipt][id][k2],icent);
          double v22Dif = cor2red.mVal/v22;
          double cov22prime = Covariance(pCovDiff[ipt][id][k22prime],pCorr[k2],pCorrDiff[ipt][id][k2],icent,icent,icent);
          double ev22Dif = sqrt(0.25*pow(cor2.mVal,-3)*(pow(cor2red.mVal,2)*cor2.mMSE
                              + 4*pow(cor2.mVal,2)*cor2red.mMSE - 4*cor2.mVal*cor2red.mVal*cov22prime));

          term cor4red = term(pCorrDiff[ipt][id][k4],pCorrDiffWSquare[ipt][id][k4],icent);
          double cov24prime = Covariance(pCovDiff[ipt][id][k24prime],pCorr[k2],pCorrDiff[ipt][id][k4],icent,icent,icent);
          double cov42prime = Covariance(pCovDiff[ipt][id][k24prime],pCorr[k4],pCorrDiff[ipt][id][k2],icent,icent,icent);
          double cov44prime = Covariance(pCovDiff[ipt][id][k44prime],pCorr[k4],pCorrDiff[ipt][id][k4],icent,icent,icent);
          double cov2prime4prime = Covariance(pCovDiff[ipt][id][k2prime4prime],pCorrDiff[ipt][id][k2],pCorrDiff[ipt][id][k4],icent,icent,icent);
          double v24Dif = (2.*cor2.mVal*cor2red.mVal-cor4red.mVal)*pow(v24,-3);
          double ev24Dif = sqrt( pow(v24,-14)
              * (pow(2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal,2.)
              * cor2.mMSE
              + 9./16*pow(2.*cor2.mVal*cor2red.mVal-cor4red.mVal,2.)*cor4.mMSE
              + 4*pow(cor2.mVal,2)*pow(v24,8)*cor2red.mMSE
              + pow(v24,8)*cor4red.mMSE
              - 1.5*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
              * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
              * cov24
              - 4*cor2.mVal*pow(v24,4)
              * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
              * cov22prime
              + 2*pow(v24,4)
              * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
              * cov24prime
              + 3*cor2.mVal*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
              * cov42prime
              - 1.5*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
              * cov44prime
              - 4*cor2.mVal*pow(v24,8)*cov2prime4prime));

          term cor6red = term(pCorrDiff[ipt][id][k6],pCorrDiffWSquare[ipt][id][k6],icent);
          double cov26prime = Covariance(pCovDiff[ipt][id][k26prime],pCorr[k2],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov2prime6 = Covariance(pCovDiff[ipt][id][k2prime6],pCorrDiff[ipt][id][k2],pCorr[k6],icent,icent,icent);
          double cov2prime6prime = Covariance(pCovDiff[ipt][id][k2prime6prime],pCorrDiff[ipt][id][k2],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov46prime = Covariance(pCovDiff[ipt][id][k46prime],pCorr[k4],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov4prime6 = Covariance(pCovDiff[ipt][id][k4prime6],pCorrDiff[ipt][id][k4],pCorr[k6],icent,icent,icent);
          double cov4prime6prime = Covariance(pCovDiff[ipt][id][k4prime6prime],pCorrDiff[ipt][id][k4],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double cov66prime = Covariance(pCovDiff[ipt][id][k66prime],pCorr[k6],pCorrDiff[ipt][id][k6],icent,icent,icent);
          double dNumeratorDiffV26 = cor6red.mVal-9*cor2.mVal*cor4red.mVal+12*cor2red.mVal*pow(cor2.mVal,2);
          double v26Dif = dNumeratorDiffV26 / (4.*pow(v26,5));
          double dFourTimesV26RefToTheSixthPower = 4.*pow(v26,6.);
          double dPartialDerivativeOfDiffV26WRT2Corr    = (5.*pow(4.,1./3.)*(9.*cor4.mVal-36.*pow(cor2.mVal,2.))*dNumeratorDiffV26/(12.*pow(dFourTimesV26RefToTheSixthPower,11./6.)))
                                                          - (pow(4.,1./3.)*(9.*cor4red.mVal-24.*cor2.mVal*cor2red.mVal)/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.)));
          double dPartialDerivativeOfDiffV26WRT2CorrRed = 6.*pow(4.,1./3.)*pow(cor2.mVal,2.)/pow(dFourTimesV26RefToTheSixthPower,5./6.);
          double dPartialDerivativeOfDiffV26WRT4Corr    = 15.*pow(4.,1./3.)*cor2.mVal*dNumeratorDiffV26/(4.*pow(dFourTimesV26RefToTheSixthPower,11./6.));
          double dPartialDerivativeOfDiffV26WRT4CorrRed = (-1.)*9.*pow(4.,1./3.)*cor2.mVal/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.));
          double dPartialDerivativeOfDiffV26WRT6Corr    = (-1.)*5.*pow(4.,1./3.)*dNumeratorDiffV26/(12.*pow(dFourTimesV26RefToTheSixthPower,11./6.));
          double dPartialDerivativeOfDiffV26WRT6CorrRed = pow(4.,1./3.)/(2.*pow(dFourTimesV26RefToTheSixthPower,5./6.));
          double ev26Dif = sqrt(
            pow(dPartialDerivativeOfDiffV26WRT2Corr,    2.)*cor2.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT2CorrRed, 2.)*cor2red.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT4Corr,    2.)*cor4.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT4CorrRed, 2.)*cor4red.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT6Corr,    2.)*cor6.mMSE
          + pow(dPartialDerivativeOfDiffV26WRT6CorrRed, 2.)*cor6red.mMSE
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT2CorrRed) * cov22prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT4Corr)    * cov24
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov24prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT6Corr)    * cov26
          + 2.*(dPartialDerivativeOfDiffV26WRT2Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov26prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT4Corr   ) * cov42prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov2prime4prime
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT6Corr   ) * cov2prime6
          + 2.*(dPartialDerivativeOfDiffV26WRT2CorrRed * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov2prime6prime
          + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT4CorrRed) * cov44prime
          + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT6Corr   ) * cov46
          + 2.*(dPartialDerivativeOfDiffV26WRT4Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov46prime
          + 2.*(dPartialDerivativeOfDiffV26WRT4CorrRed * dPartialDerivativeOfDiffV26WRT6Corr)     * cov4prime6
          + 2.*(dPartialDerivativeOfDiffV26WRT4CorrRed * dPartialDerivativeOfDiffV26WRT6CorrRed)  * cov4prime6prime
          + 2.*(dPartialDerivativeOfDiffV26WRT6Corr * dPartialDerivativeOfDiffV26WRT6CorrRed) * cov66prime
          );

          if (imeth == kStd) {
            v2[0][id][icent]  = v22Dif;
            v2e[0][id][icent] = ev22Dif;
            v2[1][id][icent]  = v24Dif;
            v2e[1][id][icent] = ev24Dif;
            v2[2][id][icent]  = v26Dif;
            v2e[2][id][icent] = ev26Dif;
          } else {
            v2[3][id][icent]  = v22Dif;
            v2e[3][id][icent] = ev22Dif;
            v2[4][id][icent]  = v24Dif;
            v2e[4][id][icent] = ev24Dif;
            v2[5][id][icent]  = v26Dif;
            v2e[5][id][icent] = ev26Dif;
          }
        } // for(Int_t ipt=binMinPtRFP; ipt<binMinPtRFP+1; ipt++)
      } // for (Int_t id = 0; id<npid; id++)
    } // for (Int_t icent=0; icent<ncent; icent++)
  } // for(Int_t imeth(0);imeth<kNMethod;imeth++)

  for (Int_t imeth=0; imeth<nmethod; imeth++){
    grRefFl[imeth] = new TGraphErrors(ncent,cent,v2_RF[imeth],centErr,v2e_RF[imeth]);
    grRefFl[imeth] -> SetMarkerStyle(marker[imeth]);
    grRefFl[imeth] -> SetMarkerSize(1.5);
    // grRefFl[imeth] -> RemovePoint(0);
    for (Int_t id=0; id<npid; id++){
      grIntFlPID[imeth][id] = new TGraphErrors(ncent,cent,v2[imeth][id],centErr,v2e[imeth][id]);
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
  PlotV2QCumulantNEW();
  return 0;
}
