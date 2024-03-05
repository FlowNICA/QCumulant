#include <FlowAnalysisWithQCumulantGenericFramework.h>
ClassImp(FlowAnalysisWithQCumulantGenericFramework);
FlowAnalysisWithQCumulantGenericFramework::FlowAnalysisWithQCumulantGenericFramework() :
M(0),
fEtaGap(0.),
fUseWeight(kFALSE)
{
  Zero();
}

FlowAnalysisWithQCumulantGenericFramework::~FlowAnalysisWithQCumulantGenericFramework()
{
}

void FlowAnalysisWithQCumulantGenericFramework::Init()
{
  for (Int_t icent = 0; icent < ncent; icent++)
  { // loop over centrality classes
    for(Int_t ipt(0);ipt<npt;ipt++) for (Int_t id(0);id<npid;id++) hcounter[icent][ipt][id] = new TProfile(Form("hcounter_%i_%i_%i", icent, ipt, id), Form("hcounter_%i_%i_%i", icent, ipt, id), 2, 0., 2.);
    TString strCovTerm[kNCov] = {"Cov24","Cov26","Cov46"};
    TString strDiffCovTerm[kNDiffCov] = {"Cov22prime", "Cov24prime", "Cov26prime", "Cov2prime4", "Cov2prime4prime", "Cov2prime6", "Cov2prime6prime", "Cov44prime", "Cov46prime", "Cov4prime6", "Cov4prime6prime", "Cov66prime"};
    TString strCorr[kNCorr] = {"<2>","<4>","<6>"};
    TString strDiffCorr[kNDiffCorr] = {"<2'>","<4'>","<6'>"};
    TString strMethod[kNMethod] = {"Standard", "Sub-event"};
    for(Int_t imeth(0);imeth<kNMethod;imeth++)
    {
      for (Int_t i(0);i<kNCorr;i++) pCorr[icent][imeth][i] = new TProfile(Form("pCorr_%i_%i_%i",icent,imeth,i),Form("%s %s centrality %i",strMethod[imeth].Data(), strCorr[i].Data(), icent),1,0.,1.);
      for (Int_t i(0);i<kNCorr;i++) pCorrWSquare[icent][imeth][i] = new TProfile(Form("pCorrWSquare_%i_%i_%i",icent,imeth,i),Form("%s (squared weight) %s centrality %i",strMethod[imeth].Data(), strCorr[i].Data(), icent),1,0.,1.);
      for (Int_t i(0);i<kNCov;i++) pCov[icent][imeth][i] = new TProfile(Form("pCov_%i_%i_%i",icent,imeth,i),Form("%s %s centrality %i",strMethod[imeth].Data(), strCovTerm[i].Data(), icent),1,0.,1.);
      for (Int_t id = 0; id < npid; id++)
      {
        for (Int_t ipt = 0; ipt < npt; ipt++)
        { // loop over pt bin
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i] = new TProfile(Form("pCorrDiff_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i),Form("%s %s centrality %i; pt %i; %s",strMethod[imeth].Data(), strDiffCorr[i].Data(), icent,ipt,pidNames.at(id).Data()),1,0.,1.);
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i] = new TProfile(Form("pCorrDiffWSquare_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i),Form("%s (squared weight) %s centrality %i pt %i %s",strMethod[imeth].Data(), strDiffCorr[i].Data(), icent,ipt,pidNames.at(id).Data()),1,0.,1.);
          for (Int_t i(0);i<kNDiffCov;i++) pCovDiff[icent][ipt][id][imeth][i] = new TProfile(Form("pCovDiff_%i_%i_%i_%i_%i",icent,ipt,id,imeth,i),Form("%s %s centrality %i pt %i %s",strMethod[imeth].Data(), strDiffCovTerm[i].Data(), icent,ipt,pidNames.at(id).Data()),1,0.,1.);
        }
      }
    }
  } // end of loop over centrality classes

}

void FlowAnalysisWithQCumulantGenericFramework::Zero()
{
  M = 0;
  for(Int_t ipt(0);ipt<npt;ipt++){
    for(Int_t id(0);id<npid;id++){
      mp[ipt][id] = 0;
      for (Int_t ieta(0);ieta<2;ieta++){
        mpGap[ieta][ipt][id] = 0;
      }
    }
  }
        
  for(Int_t h=0;h<maxHarmonic;h++)
  {
    for(Int_t p=0;p<maxPower;p++)
    {
      fQvector[h][p](0.,0.);
      fQvectorGapPos[h][p](0.,0.);
      fQvectorGapNeg[h][p](0.,0.);
      for(Int_t ipt(0);ipt<npt;ipt++)
      {
        for(Int_t id(0);id<npid;id++)
        {
          fPvector[h][p][ipt][id](0.,0.);
          fSvector[h][p][ipt][id](0.,0.);
          fPvectorGapPos[h][p][ipt][id](0.,0.);
          fSvectorGapPos[h][p][ipt][id](0.,0.);
          fPvectorGapNeg[h][p][ipt][id](0.,0.);
          fSvectorGapNeg[h][p][ipt][id](0.,0.);        
        }
      }
    }
  }
}

void FlowAnalysisWithQCumulantGenericFramework::ProcessFirstTrackLoopRP(const Double_t &eta, const Double_t &phi)
{
  Double_t dWeight = 1.0;
  for(Int_t h=0;h<maxHarmonic;h++)
  {
    for(Int_t p=0;p<maxPower;p++)
    {
      fQvector[h][p] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
      if (eta <-fEtaGap)
        fQvectorGapNeg[h][p] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
      if (eta > fEtaGap)
        fQvectorGapPos[h][p] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
    }
  }
  M++;
}

void FlowAnalysisWithQCumulantGenericFramework::ProcessFirstTrackLoopPOI(const Double_t &eta, const Double_t &phi, const Double_t &pt, const Int_t &pid, const Double_t &charge, const Double_t &pTMinRef, const Double_t &pTMaxRef)
{
  Int_t ipt = -1;
  for (Int_t j(0);j<npt;j++) { if (pt>=pTBin[j] && pt<pTBin[j+1]) ipt=j; }
  if (ipt<0) return;
  Double_t dWeight = 1.0;
  if (charge!=0) mp[ipt][8]++;
  if (charge>0) mp[ipt][0]++;
  if (charge<0) mp[ipt][4]++;
  if (pid>0) mp[ipt][pid]++;
  if (pid==1 || pid ==5) mp[ipt][9]++;
  if (pid==2 || pid ==6) mp[ipt][10]++;
  if (pid==3 || pid ==7) mp[ipt][11]++;
  if (eta < -fEtaGap)
  {
    if (charge!=0) mpGap[0][ipt][8]++;
    if (charge>0) mpGap[0][ipt][0]++;
    if (charge<0) mpGap[0][ipt][4]++;
    if (pid>0) mpGap[0][ipt][pid]++;
    if (pid==1 || pid ==5) mpGap[0][ipt][9]++;
    if (pid==2 || pid ==6) mpGap[0][ipt][10]++;
    if (pid==3 || pid ==7) mpGap[0][ipt][11]++;
  }
  if (eta > fEtaGap)
  {
    if (charge!=0) mpGap[1][ipt][8]++;
    if (charge>0) mpGap[1][ipt][0]++;
    if (charge<0) mpGap[1][ipt][4]++;
    if (pid>0) mpGap[1][ipt][pid]++;
    if (pid==1 || pid ==5) mpGap[1][ipt][9]++;
    if (pid==2 || pid ==6) mpGap[1][ipt][10]++;
    if (pid==3 || pid ==7) mpGap[1][ipt][11]++;
  }  
  for(Int_t h=0;h<maxHarmonic;h++){
    for(Int_t p=0;p<maxPower;p++){
      if (charge != 0)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][8] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][8] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][8] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][8] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][8] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][8] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }

      if (charge > 0)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][0] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][0] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][0] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][0] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][0] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][0] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }
      if (charge < 0)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][4] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][4] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][4] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][4] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][4] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][4] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }
      if (pid > 0)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][pid] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][pid] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][pid] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][pid] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][pid] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][pid] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }
      if (pid == 1 || pid ==5)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][9] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][9] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][9] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][9] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][9] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][9] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }
      if (pid == 2 || pid ==6)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][10] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][10] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][10] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][10] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][10] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][10] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }
      if (pid == 3 || pid ==7)
      {
        if (pt > pTMinRef && pt < pTMaxRef) fSvector[h][p][ipt][11] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        fPvector[h][p][ipt][11] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
        if (eta < -fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapNeg[h][p][ipt][11] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapNeg[h][p][ipt][11] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
        if (eta > fEtaGap)
        {
          if (pt > pTMinRef && pt < pTMaxRef) fSvectorGapPos[h][p][ipt][11] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));
          fPvectorGapPos[h][p][ipt][11] += TMath::Power(dWeight,p)*TComplex(TMath::Cos(h*phi),TMath::Sin(h*phi));  
        }
      }      
    }
  }
} 

void FlowAnalysisWithQCumulantGenericFramework::ProcessEventAfterFirstTrackLoop(const Int_t &icent)
{
  Double_t w2 = Two(0,0).Re();
  Double_t w4 = Four(0,0,0,0).Re();
  Double_t w6 = Six(0,0,0,0,0,0).Re();
  if (w2>0. && w4>0. && w6>0.)
  {
    Double_t cor2 = Two(h1,h2).Re() / w2;
    Double_t cor4 = Four(h1,h2,h3,h4).Re() / w4;
    Double_t cor6 = Six(h1,h3,h5,h2,h4,h6).Re() / w6;

    pCorr[icent][kStd][k2]->Fill(0.5, cor2, w2); // <<2>>
    pCorr[icent][kStd][k4]->Fill(0.5, cor4, w4); // <<4>>
    pCorr[icent][kStd][k6]->Fill(0.5, cor6, w6); // <<6>>
    pCorrWSquare[icent][kStd][k2]->Fill(0.5, cor2, w2*w2);
    pCorrWSquare[icent][kStd][k4]->Fill(0.5, cor4, w4*w4);
    pCorrWSquare[icent][kStd][k6]->Fill(0.5, cor6, w6*w6);

    pCov[icent][kStd][k24]->Fill(0.5, cor2 * cor4, w2 * w4); // <2>*<4>
    pCov[icent][kStd][k26]->Fill(0.5, cor2 * cor6, w2 * w6); // <2>*<6>
    pCov[icent][kStd][k46]->Fill(0.5, cor4 * cor6, w4 * w6); // <4>*<6>
    for (Int_t ipt = 0; ipt < npt; ipt++)
    {
      for (Int_t id = 0; id < npid; id++)
      {
        Double_t wred2 = TwoDiff(0,0,ipt,id).Re();
        Double_t wred4 = FourDiff(0,0,0,0,ipt,id).Re();
        Double_t wred6 = SixDiff(0,0,0,0,0,0,ipt,id).Re();

        if (wred2>0. && wred4>0. && wred6>0.)
        {
          Double_t redCor2 = TwoDiff(h1,h2,ipt,id).Re() / wred2;
          Double_t redCor4 = FourDiff(h1,h2,h3,h4,ipt,id).Re() / wred4;
          Double_t redCor6 = SixDiff(h1,h3,h5,h2,h4,h6,ipt,id).Re() / wred6;

          pCorrDiff[icent][ipt][id][kStd][k2]->Fill(0.5, redCor2, wred2);
          pCorrDiff[icent][ipt][id][kStd][k4]->Fill(0.5, redCor4, wred4);
          pCorrDiff[icent][ipt][id][kStd][k6]->Fill(0.5, redCor6, wred6);
          pCorrDiffWSquare[icent][ipt][id][kStd][k2]->Fill(0.5, redCor2, wred2*wred2);
          pCorrDiffWSquare[icent][ipt][id][kStd][k4]->Fill(0.5, redCor4, wred4*wred4);
          pCorrDiffWSquare[icent][ipt][id][kStd][k6]->Fill(0.5, redCor6, wred6*wred6);


          pCovDiff[icent][ipt][id][kStd][k22prime]      ->Fill(0.5, cor2*redCor2,     w2*wred2);
          pCovDiff[icent][ipt][id][kStd][k24prime]      ->Fill(0.5, cor2*redCor4,     w2*wred4);
          pCovDiff[icent][ipt][id][kStd][k26prime]      ->Fill(0.5, cor2*redCor6,     w2*wred6);
          pCovDiff[icent][ipt][id][kStd][k2prime4]      ->Fill(0.5, redCor2*cor4,     wred2*w4);
          pCovDiff[icent][ipt][id][kStd][k2prime4prime] ->Fill(0.5, redCor2*redCor4,  wred2*wred4);
          pCovDiff[icent][ipt][id][kStd][k2prime6]      ->Fill(0.5, redCor2*cor6,     wred2*w6);
          pCovDiff[icent][ipt][id][kStd][k2prime6prime] ->Fill(0.5, redCor2*redCor6,  wred2*wred6);
          pCovDiff[icent][ipt][id][kStd][k44prime]      ->Fill(0.5, cor4*redCor4,     w4*wred4);
          pCovDiff[icent][ipt][id][kStd][k46prime]      ->Fill(0.5, cor4*redCor6,     w4*wred6);
          pCovDiff[icent][ipt][id][kStd][k4prime6]      ->Fill(0.5, redCor4*cor6,     wred4*w6);
          pCovDiff[icent][ipt][id][kStd][k4prime6prime] ->Fill(0.5, redCor4*redCor6,  wred4*wred6);
          pCovDiff[icent][ipt][id][kStd][k66prime]      ->Fill(0.5, cor6*redCor6,     w6*wred6);

          hcounter[icent][ipt][id]->Fill(0.5, 1, mp[ipt][id]);
        }
      }
    }
  }

  Double_t w2Gap = TwoGap(0,0).Re();
  Double_t w4Gap = FourGap(0,0,0,0).Re();
  Double_t w6Gap = SixGap(0,0,0,0,0,0).Re();

  if (w2Gap>0. && w4Gap>0. && w6Gap>0.)
  {
    Double_t cor2Gap = TwoGap(h1,h2).Re() / w2Gap;
    Double_t cor4Gap = FourGap(h1,h3,h2,h4).Re() / w4Gap;
    Double_t cor6Gap = SixGap(h1,h3,h5,h2,h4,h6).Re() / w6Gap;

    pCorr[icent][kGap][k2]->Fill(0.5, cor2Gap, w2Gap); // <<2>>
    pCorr[icent][kGap][k4]->Fill(0.5, cor4Gap, w4Gap); // <<4>>
    pCorr[icent][kGap][k6]->Fill(0.5, cor6Gap, w6Gap); // <<6>>
    pCorrWSquare[icent][kGap][k2]->Fill(0.5, cor2Gap, w2Gap*w2Gap);
    pCorrWSquare[icent][kGap][k4]->Fill(0.5, cor4Gap, w4Gap*w4Gap);
    pCorrWSquare[icent][kGap][k6]->Fill(0.5, cor6Gap, w6Gap*w6Gap);    
    pCov[icent][kGap][k24]->Fill(0.5, cor2Gap*cor4Gap, w2Gap*w4Gap); // <2>*<4>
    pCov[icent][kGap][k26]->Fill(0.5, cor2Gap*cor6Gap, w2Gap*w6Gap); // <2>*<6>
    pCov[icent][kGap][k46]->Fill(0.5, cor4Gap*cor6Gap, w4Gap*w6Gap); // <4>*<6>
    for (Int_t ipt = 0; ipt < npt; ipt++)
    {
      for (Int_t id = 0; id < npid; id++)
      {
        Double_t wred2GapPos = TwoDiffGapPos(0,0,ipt,id).Re();
        Double_t wred4GapPos = FourDiffGapPos(0,0,0,0,ipt,id).Re();
        Double_t wred6GapPos = SixDiffGapPos(0,0,0,0,0,0,ipt,id).Re();
        if (wred2GapPos>0. && wred4GapPos>0. && wred6GapPos>0.)
        {
          Double_t redCor2GapPos = TwoDiffGapPos(h1,h2,ipt,id).Re() / wred2GapPos;
          Double_t redCor4GapPos = FourDiffGapPos(h1,h3,h2,h4,ipt,id).Re() / wred4GapPos;
          Double_t redCor6GapPos = SixDiffGapPos(h1,h3,h5,h2,h4,h6,ipt,id).Re() / wred6GapPos;

          pCorrDiff[icent][ipt][id][kGap][k2]->Fill(0.5, redCor2GapPos, wred2GapPos);
          pCorrDiff[icent][ipt][id][kGap][k4]->Fill(0.5, redCor4GapPos, wred4GapPos);
          pCorrDiff[icent][ipt][id][kGap][k6]->Fill(0.5, redCor6GapPos, wred6GapPos);
          pCorrDiffWSquare[icent][ipt][id][kGap][k2]->Fill(0.5, redCor2GapPos, wred2GapPos*wred2GapPos);
          pCorrDiffWSquare[icent][ipt][id][kGap][k4]->Fill(0.5, redCor4GapPos, wred4GapPos*wred4GapPos);
          pCorrDiffWSquare[icent][ipt][id][kGap][k6]->Fill(0.5, redCor6GapPos, wred6GapPos*wred6GapPos);

          pCovDiff[icent][ipt][id][kGap][k22prime]      ->Fill(0.5, cor2Gap*redCor2GapPos,        w2Gap*wred2GapPos);
          pCovDiff[icent][ipt][id][kGap][k24prime]      ->Fill(0.5, cor2Gap*redCor4GapPos,        w2Gap*wred4GapPos);
          pCovDiff[icent][ipt][id][kGap][k26prime]      ->Fill(0.5, cor2Gap*redCor6GapPos,        w2Gap*wred6GapPos);
          pCovDiff[icent][ipt][id][kGap][k2prime4]      ->Fill(0.5, redCor2GapPos*cor4Gap,        wred2GapPos*w4Gap);
          pCovDiff[icent][ipt][id][kGap][k2prime4prime] ->Fill(0.5, redCor2GapPos*redCor4GapPos,  wred2GapPos*wred4GapPos);
          pCovDiff[icent][ipt][id][kGap][k2prime6]      ->Fill(0.5, redCor2GapPos*cor6Gap,        wred2GapPos*w6Gap);
          pCovDiff[icent][ipt][id][kGap][k2prime6prime] ->Fill(0.5, redCor2GapPos*redCor6GapPos,  wred2GapPos*wred6GapPos);
          pCovDiff[icent][ipt][id][kGap][k44prime]      ->Fill(0.5, cor4Gap*redCor4GapPos,        w4Gap*wred4GapPos);
          pCovDiff[icent][ipt][id][kGap][k46prime]      ->Fill(0.5, cor4Gap*redCor6GapPos,        w4Gap*wred6GapPos);
          pCovDiff[icent][ipt][id][kGap][k4prime6]      ->Fill(0.5, redCor4GapPos*cor6Gap,        wred4GapPos*w6Gap);
          pCovDiff[icent][ipt][id][kGap][k4prime6prime] ->Fill(0.5, redCor4GapPos*redCor6GapPos,  wred4GapPos*wred6GapPos);
          pCovDiff[icent][ipt][id][kGap][k66prime]      ->Fill(0.5, cor6Gap*redCor6GapPos,        w6Gap*wred6GapPos);

          hcounter[icent][ipt][id]->Fill(1.5, 1, mpGap[1][ipt][id]);
        }
        Double_t wred2GapNeg = TwoDiffGapNeg(0,0,ipt,id).Re();
        Double_t wred4GapNeg = FourDiffGapNeg(0,0,0,0,ipt,id).Re();
        Double_t wred6GapNeg = SixDiffGapNeg(0,0,0,0,0,0,ipt,id).Re();

        if (wred2GapNeg>0. && wred4GapNeg>0. && wred6GapNeg>0.)
        {
          Double_t redCor2GapNeg = TwoDiffGapNeg(h1,h2,ipt,id).Re() / wred2GapNeg;
          Double_t redCor4GapNeg = FourDiffGapNeg(h1,h3,h2,h4,ipt,id).Re() / wred4GapNeg;
          Double_t redCor6GapNeg = SixDiffGapNeg(h1,h3,h5,h2,h4,h6,ipt,id).Re() / wred6GapNeg;

          pCorrDiff[icent][ipt][id][kGap][k2]->Fill(0.5, redCor2GapNeg, wred2GapNeg);
          pCorrDiff[icent][ipt][id][kGap][k4]->Fill(0.5, redCor4GapNeg, wred4GapNeg);
          pCorrDiff[icent][ipt][id][kGap][k6]->Fill(0.5, redCor6GapNeg, wred6GapNeg);
          pCorrDiffWSquare[icent][ipt][id][kGap][k2]->Fill(0.5, redCor2GapNeg, wred2GapNeg*wred2GapNeg);
          pCorrDiffWSquare[icent][ipt][id][kGap][k4]->Fill(0.5, redCor4GapNeg, wred4GapNeg*wred4GapNeg);
          pCorrDiffWSquare[icent][ipt][id][kGap][k6]->Fill(0.5, redCor6GapNeg, wred6GapNeg*wred6GapNeg);

          pCovDiff[icent][ipt][id][kGap][k22prime]      ->Fill(0.5, cor2Gap*redCor2GapNeg,        w2Gap*wred2GapNeg);
          pCovDiff[icent][ipt][id][kGap][k24prime]      ->Fill(0.5, cor2Gap*redCor4GapNeg,        w2Gap*wred4GapNeg);
          pCovDiff[icent][ipt][id][kGap][k26prime]      ->Fill(0.5, cor2Gap*redCor6GapNeg,        w2Gap*wred6GapNeg);
          pCovDiff[icent][ipt][id][kGap][k2prime4]      ->Fill(0.5, redCor2GapNeg*cor4Gap,        wred2GapNeg*w4Gap);
          pCovDiff[icent][ipt][id][kGap][k2prime4prime] ->Fill(0.5, redCor2GapNeg*redCor4GapNeg,  wred2GapNeg*wred4GapNeg);
          pCovDiff[icent][ipt][id][kGap][k2prime6]      ->Fill(0.5, redCor2GapNeg*cor6Gap,        wred2GapNeg*w6Gap);
          pCovDiff[icent][ipt][id][kGap][k2prime6prime] ->Fill(0.5, redCor2GapNeg*redCor6GapNeg,  wred2GapNeg*wred6GapNeg);
          pCovDiff[icent][ipt][id][kGap][k44prime]      ->Fill(0.5, cor4Gap*redCor4GapNeg,        w4Gap*wred4GapNeg);
          pCovDiff[icent][ipt][id][kGap][k46prime]      ->Fill(0.5, cor4Gap*redCor6GapNeg,        w4Gap*wred6GapNeg);
          pCovDiff[icent][ipt][id][kGap][k4prime6]      ->Fill(0.5, redCor4GapNeg*cor6Gap,        wred4GapNeg*w6Gap);
          pCovDiff[icent][ipt][id][kGap][k4prime6prime] ->Fill(0.5, redCor4GapNeg*redCor6GapNeg,  wred4GapNeg*wred6GapNeg);
          pCovDiff[icent][ipt][id][kGap][k66prime]      ->Fill(0.5, cor6Gap*redCor6GapNeg,        w6Gap*wred6GapNeg);

          hcounter[icent][ipt][id]->Fill(1.5, 1, mpGap[0][ipt][id]);
        }
      }
    }
  }
  // Double_t w8Gap = EightGap(0,0,0,0,0,0,0,0).Re();
  // if (w8Gap!=0.)
  // {
  //   Double_t cor8Gap = EightGap(h1,h3,h5,h7,h2,h4,h6,h8).Re() / w8Gap;
  //   hv28Gap[icent]->Fill(0.5, cor8Gap, w8Gap);
  //   for (Int_t ipt = 0; ipt < npt; ipt++)
  //   {
  //     for (Int_t id = 0; id < npid; id++)
  //     {
  //       Double_t wred8GapPos = EightDiffGapPos(0,0,0,0,0,0,0,0,ipt,id).Re();
  //       if (wred8GapPos>0.)
  //       {
  //         Double_t redCor28GapPos = EightDiffGapPos(h1,h3,h5,h7,h2,h4,h6,h8,ipt,id).Re() / wred8GapPos;
  //         hv28ptGap[icent][ipt][id]->Fill(0.5, redCor28GapPos, wred8GapPos);
  //       }
  //       Double_t wred8GapNeg = EightDiffGapNeg(0,0,0,0,0,0,0,0,ipt,id).Re();
  //       if (wred8GapNeg>0.)
  //       {
  //         Double_t redCor28GapNeg = EightDiffGapNeg(h1,h3,h5,h7,h2,h4,h6,h8,ipt,id).Re() / wred8GapNeg;
  //         hv28ptGap[icent][ipt][id]->Fill(0.5, redCor28GapNeg, wred8GapNeg);
  //       }
  //     }
  //   }
  // }
}
void FlowAnalysisWithQCumulantGenericFramework::SaveHist()
{
  for (Int_t icent = 0; icent < ncent; icent++)
  {
    for(Int_t ipt(0);ipt<npt;ipt++) for (Int_t id(0);id<npid;id++) hcounter[icent][ipt][id]->Write();
    for(Int_t imeth(0);imeth<kNMethod;imeth++)
    {
      for (Int_t i(0);i<kNCorr;i++) pCorr[icent][imeth][i]->Write();
      for (Int_t i(0);i<kNCorr;i++) pCorrWSquare[icent][imeth][i]->Write();
      for (Int_t i(0);i<kNCov;i++) pCov[icent][imeth][i]->Write();
      for (Int_t id = 0; id < npid; id++)
      {
        for (Int_t ipt = 0; ipt < npt; ipt++)
        {
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiff[icent][ipt][id][imeth][i]->Write();
          for (Int_t i(0);i<kNDiffCorr;i++) pCorrDiffWSquare[icent][ipt][id][imeth][i]->Write();
          for (Int_t i(0);i<kNDiffCov;i++) pCovDiff[icent][ipt][id][imeth][i]->Write();
        }
      }
    }
  }
}

void FlowAnalysisWithQCumulantGenericFramework::SaveHist(TDirectoryFile *const &outputDir)
{
  for (Int_t icent = 0; icent < ncent; icent++)
  {
    for(Int_t ipt(0);ipt<npt;ipt++) for (Int_t id(0);id<npid;id++) outputDir->Add(hcounter[icent][ipt][id]);
    for(Int_t imeth(0);imeth<kNMethod;imeth++)
    {
      for (Int_t i(0);i<kNCorr;i++) outputDir->Add(pCorr[icent][imeth][i]);
      for (Int_t i(0);i<kNCorr;i++) outputDir->Add(pCorrWSquare[icent][imeth][i]);
      for (Int_t i(0);i<kNCov;i++) outputDir->Add(pCov[icent][imeth][i]);
      for (Int_t id = 0; id < npid; id++)
      {
        for (Int_t ipt = 0; ipt < npt; ipt++)
        {
          for (Int_t i(0);i<kNDiffCorr;i++) outputDir->Add(pCorrDiff[icent][ipt][id][imeth][i]);
          for (Int_t i(0);i<kNDiffCorr;i++) outputDir->Add(pCorrDiffWSquare[icent][ipt][id][imeth][i]);
          for (Int_t i(0);i<kNDiffCov;i++) outputDir->Add(pCovDiff[icent][ipt][id][imeth][i]);
        }
      }
    }
  }
  outputDir->Write();
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Q(const Int_t n, const Int_t p) const
{
  if (n < 0) return TComplex::Conjugate(fQvector[-n][p]);
  else return fQvector[n][p];
}

TComplex FlowAnalysisWithQCumulantGenericFramework::QGapPos(const Int_t n, const Int_t p) const
{
  if (n < 0) return TComplex::Conjugate(fQvectorGapPos[-n][p]);
  else return fQvectorGapPos[n][p];
}

TComplex FlowAnalysisWithQCumulantGenericFramework::QGapNeg(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fQvectorGapNeg[-n][p]);
  else return fQvectorGapNeg[n][p];
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Two(const Int_t n1, const Int_t n2) const
{
  TComplex two = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return two;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::TwoGap(const Int_t n1, const Int_t n2) const
{
  TComplex twoGap = QGapPos(n1,1)*QGapNeg(n2,1);
  return twoGap;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Three(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex three = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)-Q(n1,1)*Q(n2+n3,2)+2.0*Q(n1+n2+n3,3);
  return three; 
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Four(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex four = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
  return four;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::FourGap(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapPos(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)
                    -QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3+n4,2)+QGapPos(n1+n2,2)*QGapNeg(n3+n4,2);
  //same as
  //TComplex formula = TwoPos(n1,n2)*TwoNeg(n3,n4);
	return formula;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Five(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5) const
{
  TComplex five = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
  - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
  + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
  + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
  + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
  - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
  + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
  - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
  + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
  + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
  - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
  + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
  - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
  - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
  + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
  + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
  + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
  + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
  - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
  + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
  + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
  + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
  + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
  - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
  - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
  - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
  return five;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Six(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex six = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
  - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
  + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
  + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
  + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
  - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
  + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
  - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
  + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
  + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
  - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
  + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
  - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
  - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
  + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
  + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
  + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
  + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
  - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
  + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
  + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
  + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
  + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
  - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
  - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
  - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
  - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
  + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
  - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
  - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
  - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
  + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
  - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
  + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
  + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
  + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
  + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
  + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
  - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
  - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
  - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
  + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
  - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
  + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
  + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
  + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
  + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
  + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
  - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
  - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
  - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
  + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
  - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
  + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
  - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
  - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
  - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
  - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
  + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
  - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
  - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
  - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
  + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
  - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
  + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
  - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
  - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
  + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
  - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
  + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
  + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
  - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
  - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
  - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
  - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
  + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
  + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
  + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
  + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
  + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
  - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
  - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
  - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
  + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
  - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
  + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
  - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
  - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
  - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
  + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
  - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
  + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
  + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
  + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
  + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
  - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
  + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
  + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
  + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
  - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
  + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
  - 120.*Q(n1+n2+n3+n4+n5+n6,6);
  return six;
}
TComplex FlowAnalysisWithQCumulantGenericFramework::ThreePos(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)-QGapPos(n1+n2,2)*QGapPos(n3,1)-QGapPos(n2,1)*QGapPos(n1+n3,2)
 		                 - QGapPos(n1,1)*QGapPos(n2+n3,2)+2.*QGapPos(n1+n2+n3,3);
  return formula;
}
TComplex FlowAnalysisWithQCumulantGenericFramework::ThreeNeg(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)-QGapNeg(n1+n2,2)*QGapNeg(n3,1)-QGapNeg(n2,1)*QGapNeg(n1+n3,2)
 		                 - QGapNeg(n1,1)*QGapNeg(n2+n3,2)+2.*QGapNeg(n1+n2+n3,3);
  return formula;
}
TComplex FlowAnalysisWithQCumulantGenericFramework::SixGap(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = ThreePos(n1,n2,n3)*ThreeNeg(n4,n5,n6);
	return formula;
}
TComplex FlowAnalysisWithQCumulantGenericFramework::FourPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)-QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)-QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4,1)
                    - QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1)+2.0*QGapPos(n1+n2+n3,3)*QGapPos(n4,1)-QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4,2)
                    + QGapPos(n2+n3,2)*QGapPos(n1+n4,2)-QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)+QGapPos(n1+n3,2)*QGapPos(n2+n4,2)
                    + 2.0*QGapPos(n3,1)*QGapPos(n1+n2+n4,3)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2)+QGapPos(n1+n2,2)*QGapPos(n3+n4,2)
                    + 2.0*QGapPos(n2,1)*QGapPos(n1+n3+n4,3)+2.0*QGapPos(n1,1)*QGapPos(n2+n3+n4,3)-6.0*QGapPos(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::FourNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4,1)
                    - QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)+2.0*QGapNeg(n1+n2+n3,3)*QGapNeg(n4,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4,2)
                    + QGapNeg(n2+n3,2)*QGapNeg(n1+n4,2)-QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)+QGapNeg(n1+n3,2)*QGapNeg(n2+n4,2)
                    + 2.0*QGapNeg(n3,1)*QGapNeg(n1+n2+n4,3)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)+QGapNeg(n1+n2,2)*QGapNeg(n3+n4,2)
                    + 2.0*QGapNeg(n2,1)*QGapNeg(n1+n3+n4,3)+2.0*QGapNeg(n1,1)*QGapNeg(n2+n3+n4,3)-6.0*QGapNeg(n1+n2+n3+n4,4);
  return formula;
}
TComplex FlowAnalysisWithQCumulantGenericFramework::EightGap(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex formula = FourPos(n1,n2,n3,n4)*FourNeg(n5,n6,n7,n8);
	return formula;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::Recursion(Int_t n, Int_t* harmonic, Int_t mult /*= 1*/, Int_t skip /*= 0*/) 
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0)

TComplex FlowAnalysisWithQCumulantGenericFramework::P(const Int_t n, const Int_t p, const Int_t ipt, const Int_t pid) const
{
  if(n>=0) { return fPvector[n][p][ipt][pid]; }
  return TComplex::Conjugate(fPvector[-n][p][ipt][pid]);
}

TComplex FlowAnalysisWithQCumulantGenericFramework::PGapPos(const Int_t n, const Int_t p, const Int_t ipt, const Int_t pid) const
{
  if(n < 0) return TComplex::Conjugate(fPvectorGapPos[-n][p][ipt][pid]);
  else return fPvectorGapPos[n][p][ipt][pid];
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::PGapNeg(const Int_t n, const Int_t p, const Int_t ipt, const Int_t pid) const
{
  if(n < 0) return TComplex::Conjugate(fPvectorGapNeg[-n][p][ipt][pid]);
  else return fPvectorGapNeg[n][p][ipt][pid];
}

TComplex FlowAnalysisWithQCumulantGenericFramework::S(const Int_t n, const Int_t p, const Int_t ipt, const Int_t pid) const
{
  if(n>=0) { return fSvector[n][p][ipt][pid]; }
  return TComplex::Conjugate(fSvector[-n][p][ipt][pid]);
}

TComplex FlowAnalysisWithQCumulantGenericFramework::SGapPos(const Int_t n, const Int_t p, const Int_t ipt, const Int_t pid) const
{
  if(n < 0) return TComplex::Conjugate(fSvectorGapPos[-n][p][ipt][pid]);
  else return fSvectorGapPos[n][p][ipt][pid];
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::SGapNeg(const Int_t n, const Int_t p, const Int_t ipt, const Int_t pid) const
{
  if(n < 0) return TComplex::Conjugate(fSvectorGapNeg[-n][p][ipt][pid]);
  else return fSvectorGapNeg[n][p][ipt][pid];
}

TComplex FlowAnalysisWithQCumulantGenericFramework::TwoDiff(const Int_t n1, const Int_t n2, const Int_t ipt, const Int_t pid) const
{
  TComplex twoDiff = P(n1,1,ipt,pid)*Q(n2,1) - S(n1+n2,2,ipt,pid);
  return twoDiff;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::TwoDiffGapPos(const Int_t n1, const Int_t n2, const Int_t ipt, const Int_t pid) const
{
  TComplex twoDiffGapPos = PGapPos(n1,1,ipt,pid)*QGapNeg(n2,1);
  return twoDiffGapPos;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::TwoDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t ipt, const Int_t pid) const
{
  TComplex twoDiffGapNeg = PGapNeg(n1,1,ipt,pid)*QGapPos(n2,1);
  return twoDiffGapNeg;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::ThreeDiff(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t ipt, const Int_t pid) const
{
  TComplex threeDiff = P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)-S(n1+n2,2,ipt,pid)*Q(n3,1)-S(n1+n3,2,ipt,pid)*Q(n2,1)
 		                 - P(n1,1,ipt,pid)*Q(n2+n3,2)+2.0*S(n1+n2+n3,3,ipt,pid);
  return threeDiff;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::FourDiff(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t ipt, const Int_t pid) const
{
  TComplex fourDiff = P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2,ipt,pid)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2,ipt,pid)*Q(n4,1)
                    - P(n1,1,ipt,pid)*Q(n2+n3,2)*Q(n4,1)+2.0*S(n1+n2+n3,3,ipt,pid)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2,ipt,pid)
                    + Q(n2+n3,2)*S(n1+n4,2,ipt,pid)-P(n1,1,ipt,pid)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2,ipt,pid)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*S(n1+n2+n4,3,ipt,pid)-P(n1,1,ipt,pid)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2,ipt,pid)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*S(n1+n3+n4,3,ipt,pid)+2.0*P(n1,1,ipt,pid)*Q(n2+n3+n4,3)-6.0*S(n1+n2+n3+n4,4,ipt,pid);
  return fourDiff;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::FourDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = PGapPos(n1,1,ipt,pid)*QGapPos(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)
                      - SGapPos(n1+n2,2,ipt,pid)*QGapNeg(n3,1)*QGapNeg(n4,1)
                      - PGapPos(n1,1,ipt,pid)*QGapPos(n2,1)*QGapNeg(n3+n4,2)
                      + SGapPos(n1+n2,2,ipt,pid)*QGapNeg(n3+n4,2);
  // same as
  // TComplex formula = TwoDiffPos(n1,n2)*TwoNeg(n3,n4);
  return formula;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::FourDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = PGapNeg(n1,1,ipt,pid)*QGapNeg(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)
                      - SGapNeg(n1+n2,2,ipt,pid)*QGapPos(n3,1)*QGapPos(n4,1)
                      - PGapNeg(n1,1,ipt,pid)*QGapNeg(n2,1)*QGapPos(n3+n4,2)
                      + SGapNeg(n1+n2,2,ipt,pid)*QGapPos(n3+n4,2);
  // same as
  // TComplex formula = TwoDiffNeg(n1,n2)*TwoPos(n3,n4);
  return formula;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::ThreeDiffPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = PGapPos(n1,1,ipt,pid)*QGapPos(n2,1)*QGapPos(n3,1)-SGapPos(n1+n2,2,ipt,pid)*QGapPos(n3,1)-SGapPos(n1+n3,2,ipt,pid)*QGapPos(n2,1)
 		                 - PGapPos(n1,1,ipt,pid)*QGapPos(n2+n3,2)+2.0*SGapPos(n1+n2+n3,3,ipt,pid);
  return formula;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::ThreeDiffNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = PGapNeg(n1,1,ipt,pid)*QGapNeg(n2,1)*QGapNeg(n3,1)-SGapNeg(n1+n2,2,ipt,pid)*QGapNeg(n3,1)-SGapNeg(n1+n3,2,ipt,pid)*QGapNeg(n2,1)
 		                 - PGapNeg(n1,1,ipt,pid)*QGapNeg(n2+n3,2)+2.0*SGapNeg(n1+n2+n3,3,ipt,pid);
  return formula;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::SixDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = ThreeDiffPos(n1,n2,n3,ipt,pid)*ThreeNeg(n4,n5,n6);
  return formula;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::SixDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = ThreeDiffNeg(n1,n2,n3,ipt,pid)*ThreePos(n4,n5,n6);
  return formula;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::SixDiff(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t ipt, const Int_t pid) const
{
    TComplex formula = P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-S(n1+n2,2,ipt,pid)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    - Q(n2,1)*S(n1+n3,2,ipt,pid)*Q(n4,1)*Q(n5,1)*Q(n6,1)-P(n1,1,ipt,pid)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    + 2.*S(n1+n2+n3,3,ipt,pid)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2,ipt,pid)*Q(n5,1)*Q(n6,1)
    + Q(n2+n3,2)*S(n1+n4,2,ipt,pid)*Q(n5,1)*Q(n6,1)-P(n1,1,ipt,pid)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
    + S(n1+n3,2,ipt,pid)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*S(n1+n2+n4,3,ipt,pid)*Q(n5,1)*Q(n6,1)
    - P(n1,1,ipt,pid)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+S(n1+n2,2,ipt,pid)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n2,1)*S(n1+n3+n4,3,ipt,pid)*Q(n5,1)*Q(n6,1)+2.*P(n1,1,ipt,pid)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
    - 6.*S(n1+n2+n3+n4,4,ipt,pid)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*S(n1+n5,2,ipt,pid)*Q(n6,1)
    + Q(n2+n3,2)*Q(n4,1)*S(n1+n5,2,ipt,pid)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*S(n1+n5,2,ipt,pid)*Q(n6,1)
    + Q(n2,1)*Q(n3+n4,2)*S(n1+n5,2,ipt,pid)*Q(n6,1)-2.*Q(n2+n3+n4,3)*S(n1+n5,2,ipt,pid)*Q(n6,1)
    - P(n1,1,ipt,pid)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+S(n1+n3,2,ipt,pid)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
    + Q(n3,1)*S(n1+n4,2,ipt,pid)*Q(n2+n5,2)*Q(n6,1)+P(n1,1,ipt,pid)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
    - 2.*S(n1+n3+n4,3,ipt,pid)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*S(n1+n2+n5,3,ipt,pid)*Q(n6,1)
    - 2.*Q(n3+n4,2)*S(n1+n2+n5,3,ipt,pid)*Q(n6,1)-P(n1,1,ipt,pid)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
    + S(n1+n2,2,ipt,pid)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*S(n1+n4,2,ipt,pid)*Q(n3+n5,2)*Q(n6,1)
    + P(n1,1,ipt,pid)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*S(n1+n2+n4,3,ipt,pid)*Q(n3+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n4,1)*S(n1+n3+n5,3,ipt,pid)*Q(n6,1)-2.*Q(n2+n4,2)*S(n1+n3+n5,3,ipt,pid)*Q(n6,1)
    + 2.*P(n1,1,ipt,pid)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*S(n1+n4,2,ipt,pid)*Q(n2+n3+n5,3)*Q(n6,1)
    - 6.*Q(n4,1)*S(n1+n2+n3+n5,4,ipt,pid)*Q(n6,1)-P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
    + S(n1+n2,2,ipt,pid)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*S(n1+n3,2,ipt,pid)*Q(n4+n5,2)*Q(n6,1)
    + P(n1,1,ipt,pid)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*S(n1+n2+n3,3,ipt,pid)*Q(n4+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n3,1)*S(n1+n4+n5,3,ipt,pid)*Q(n6,1)-2.*Q(n2+n3,2)*S(n1+n4+n5,3,ipt,pid)*Q(n6,1)
    + 2.*P(n1,1,ipt,pid)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*S(n1+n3,2,ipt,pid)*Q(n2+n4+n5,3)*Q(n6,1)
    - 6.*Q(n3,1)*S(n1+n2+n4+n5,4,ipt,pid)*Q(n6,1)+2.*P(n1,1,ipt,pid)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
    - 2.*S(n1+n2,2,ipt,pid)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*S(n1+n3+n4+n5,4,ipt,pid)*Q(n6,1)
    - 6.*P(n1,1,ipt,pid)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*S(n1+n2+n3+n4+n5,5,ipt,pid)*Q(n6,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*S(n1+n6,2,ipt,pid)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*S(n1+n6,2,ipt,pid)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*S(n1+n6,2,ipt,pid)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*S(n1+n6,2,ipt,pid)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*S(n1+n6,2,ipt,pid)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*S(n1+n6,2,ipt,pid)
    - Q(n3+n4,2)*Q(n2+n5,2)*S(n1+n6,2,ipt,pid)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*S(n1+n6,2,ipt,pid)
    - Q(n2+n4,2)*Q(n3+n5,2)*S(n1+n6,2,ipt,pid)-2.*Q(n4,1)*Q(n2+n3+n5,3)*S(n1+n6,2,ipt,pid)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*S(n1+n6,2,ipt,pid)-Q(n2+n3,2)*Q(n4+n5,2)*S(n1+n6,2,ipt,pid)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*S(n1+n6,2,ipt,pid)-2.*Q(n2,1)*Q(n3+n4+n5,3)*S(n1+n6,2,ipt,pid)
    + 6.*Q(n2+n3+n4+n5,4)*S(n1+n6,2,ipt,pid)-P(n1,1,ipt,pid)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
    + S(n1+n3,2,ipt,pid)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*S(n1+n4,2,ipt,pid)*Q(n5,1)*Q(n2+n6,2)
    + P(n1,1,ipt,pid)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*S(n1+n3+n4,3,ipt,pid)*Q(n5,1)*Q(n2+n6,2)
    + Q(n3,1)*Q(n4,1)*S(n1+n5,2,ipt,pid)*Q(n2+n6,2)-Q(n3+n4,2)*S(n1+n5,2,ipt,pid)*Q(n2+n6,2)
    + P(n1,1,ipt,pid)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-S(n1+n4,2,ipt,pid)*Q(n3+n5,2)*Q(n2+n6,2)
    - 2.*Q(n4,1)*S(n1+n3+n5,3,ipt,pid)*Q(n2+n6,2)+P(n1,1,ipt,pid)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
    - S(n1+n3,2,ipt,pid)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*S(n1+n4+n5,3,ipt,pid)*Q(n2+n6,2)
    - 2.*P(n1,1,ipt,pid)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*S(n1+n3+n4+n5,4,ipt,pid)*Q(n2+n6,2)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*S(n1+n2+n6,3,ipt,pid)-2.*Q(n3+n4,2)*Q(n5,1)*S(n1+n2+n6,3,ipt,pid)
    - 2.*Q(n4,1)*Q(n3+n5,2)*S(n1+n2+n6,3,ipt,pid)-2.*Q(n3,1)*Q(n4+n5,2)*S(n1+n2+n6,3,ipt,pid)
    + 4.*Q(n3+n4+n5,3)*S(n1+n2+n6,3,ipt,pid)-P(n1,1,ipt,pid)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
    + S(n1+n2,2,ipt,pid)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*S(n1+n4,2,ipt,pid)*Q(n5,1)*Q(n3+n6,2)
    + P(n1,1,ipt,pid)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*S(n1+n2+n4,3,ipt,pid)*Q(n5,1)*Q(n3+n6,2)
    + Q(n2,1)*Q(n4,1)*S(n1+n5,2,ipt,pid)*Q(n3+n6,2)-Q(n2+n4,2)*S(n1+n5,2,ipt,pid)*Q(n3+n6,2)
    + P(n1,1,ipt,pid)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-S(n1+n4,2,ipt,pid)*Q(n2+n5,2)*Q(n3+n6,2)
    - 2.*Q(n4,1)*S(n1+n2+n5,3,ipt,pid)*Q(n3+n6,2)+P(n1,1,ipt,pid)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
    - S(n1+n2,2,ipt,pid)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*S(n1+n4+n5,3,ipt,pid)*Q(n3+n6,2)
    - 2.*P(n1,1,ipt,pid)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*S(n1+n2+n4+n5,4,ipt,pid)*Q(n3+n6,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*S(n1+n3+n6,3,ipt,pid)-2.*Q(n2+n4,2)*Q(n5,1)*S(n1+n3+n6,3,ipt,pid)
    - 2.*Q(n4,1)*Q(n2+n5,2)*S(n1+n3+n6,3,ipt,pid)-2.*Q(n2,1)*Q(n4+n5,2)*S(n1+n3+n6,3,ipt,pid)
    + 4.*Q(n2+n4+n5,3)*S(n1+n3+n6,3,ipt,pid)+2.*P(n1,1,ipt,pid)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
    - 2.*S(n1+n4,2,ipt,pid)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*S(n1+n5,2,ipt,pid)*Q(n2+n3+n6,3)
    - 2.*P(n1,1,ipt,pid)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*S(n1+n4+n5,3,ipt,pid)*Q(n2+n3+n6,3)
    - 6.*Q(n4,1)*Q(n5,1)*S(n1+n2+n3+n6,4,ipt,pid)+6.*Q(n4+n5,2)*S(n1+n2+n3+n6,4,ipt,pid)
    - P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+S(n1+n2,2,ipt,pid)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
    + Q(n2,1)*S(n1+n3,2,ipt,pid)*Q(n5,1)*Q(n4+n6,2)+P(n1,1,ipt,pid)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
    - 2.*S(n1+n2+n3,3,ipt,pid)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*S(n1+n5,2,ipt,pid)*Q(n4+n6,2)
    - Q(n2+n3,2)*S(n1+n5,2,ipt,pid)*Q(n4+n6,2)+P(n1,1,ipt,pid)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
    - S(n1+n3,2,ipt,pid)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*S(n1+n2+n5,3,ipt,pid)*Q(n4+n6,2)
    + P(n1,1,ipt,pid)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-S(n1+n2,2,ipt,pid)*Q(n3+n5,2)*Q(n4+n6,2)
    - 2.*Q(n2,1)*S(n1+n3+n5,3,ipt,pid)*Q(n4+n6,2)-2.*P(n1,1,ipt,pid)*Q(n2+n3+n5,3)*Q(n4+n6,2)
    + 6.*S(n1+n2+n3+n5,4,ipt,pid)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*S(n1+n4+n6,3,ipt,pid)
    - 2.*Q(n2+n3,2)*Q(n5,1)*S(n1+n4+n6,3,ipt,pid)-2.*Q(n3,1)*Q(n2+n5,2)*S(n1+n4+n6,3,ipt,pid)
    - 2.*Q(n2,1)*Q(n3+n5,2)*S(n1+n4+n6,3,ipt,pid)+4.*Q(n2+n3+n5,3)*S(n1+n4+n6,3,ipt,pid)
    + 2.*P(n1,1,ipt,pid)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*S(n1+n3,2,ipt,pid)*Q(n5,1)*Q(n2+n4+n6,3)
    - 2.*Q(n3,1)*S(n1+n5,2,ipt,pid)*Q(n2+n4+n6,3)-2.*P(n1,1,ipt,pid)*Q(n3+n5,2)*Q(n2+n4+n6,3)
    + 4.*S(n1+n3+n5,3,ipt,pid)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*S(n1+n2+n4+n6,4,ipt,pid)
    + 6.*Q(n3+n5,2)*S(n1+n2+n4+n6,4,ipt,pid)+2.*P(n1,1,ipt,pid)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
    - 2.*S(n1+n2,2,ipt,pid)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*S(n1+n5,2,ipt,pid)*Q(n3+n4+n6,3)
    - 2.*P(n1,1,ipt,pid)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*S(n1+n2+n5,3,ipt,pid)*Q(n3+n4+n6,3)
    - 6.*Q(n2,1)*Q(n5,1)*S(n1+n3+n4+n6,4,ipt,pid)+6.*Q(n2+n5,2)*S(n1+n3+n4+n6,4,ipt,pid)
    - 6.*P(n1,1,ipt,pid)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*S(n1+n5,2,ipt,pid)*Q(n2+n3+n4+n6,4)
    + 24.*Q(n5,1)*S(n1+n2+n3+n4+n6,5,ipt,pid)-P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
    + S(n1+n2,2,ipt,pid)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*S(n1+n3,2,ipt,pid)*Q(n4,1)*Q(n5+n6,2)
    + P(n1,1,ipt,pid)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*S(n1+n2+n3,3,ipt,pid)*Q(n4,1)*Q(n5+n6,2)
    + Q(n2,1)*Q(n3,1)*S(n1+n4,2,ipt,pid)*Q(n5+n6,2)-Q(n2+n3,2)*S(n1+n4,2,ipt,pid)*Q(n5+n6,2)
    + P(n1,1,ipt,pid)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-S(n1+n3,2,ipt,pid)*Q(n2+n4,2)*Q(n5+n6,2)
    - 2.*Q(n3,1)*S(n1+n2+n4,3,ipt,pid)*Q(n5+n6,2)+P(n1,1,ipt,pid)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
    - S(n1+n2,2,ipt,pid)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*S(n1+n3+n4,3,ipt,pid)*Q(n5+n6,2)
    - 2.*P(n1,1,ipt,pid)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*S(n1+n2+n3+n4,4,ipt,pid)*Q(n5+n6,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*S(n1+n5+n6,3,ipt,pid)-2.*Q(n2+n3,2)*Q(n4,1)*S(n1+n5+n6,3,ipt,pid)
    - 2.*Q(n3,1)*Q(n2+n4,2)*S(n1+n5+n6,3,ipt,pid)-2.*Q(n2,1)*Q(n3+n4,2)*S(n1+n5+n6,3,ipt,pid)
    + 4.*Q(n2+n3+n4,3)*S(n1+n5+n6,3,ipt,pid)+2.*P(n1,1,ipt,pid)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
    - 2.*S(n1+n3,2,ipt,pid)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*S(n1+n4,2,ipt,pid)*Q(n2+n5+n6,3)
    - 2.*P(n1,1,ipt,pid)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*S(n1+n3+n4,3,ipt,pid)*Q(n2+n5+n6,3)
    - 6.*Q(n3,1)*Q(n4,1)*S(n1+n2+n5+n6,4,ipt,pid)+6.*Q(n3+n4,2)*S(n1+n2+n5+n6,4,ipt,pid)
    + 2.*P(n1,1,ipt,pid)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*S(n1+n2,2,ipt,pid)*Q(n4,1)*Q(n3+n5+n6,3)
    - 2.*Q(n2,1)*S(n1+n4,2,ipt,pid)*Q(n3+n5+n6,3)-2.*P(n1,1,ipt,pid)*Q(n2+n4,2)*Q(n3+n5+n6,3)
    + 4.*S(n1+n2+n4,3,ipt,pid)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*S(n1+n3+n5+n6,4,ipt,pid)
    + 6.*Q(n2+n4,2)*S(n1+n3+n5+n6,4,ipt,pid)-6.*P(n1,1,ipt,pid)*Q(n4,1)*Q(n2+n3+n5+n6,4)
    + 6.*S(n1+n4,2,ipt,pid)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*S(n1+n2+n3+n5+n6,5,ipt,pid)
    + 2.*P(n1,1,ipt,pid)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*S(n1+n2,2,ipt,pid)*Q(n3,1)*Q(n4+n5+n6,3)
    - 2.*Q(n2,1)*S(n1+n3,2,ipt,pid)*Q(n4+n5+n6,3)-2.*P(n1,1,ipt,pid)*Q(n2+n3,2)*Q(n4+n5+n6,3)
    + 4.*S(n1+n2+n3,3,ipt,pid)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*S(n1+n4+n5+n6,4,ipt,pid)
    + 6.*Q(n2+n3,2)*S(n1+n4+n5+n6,4,ipt,pid)-6.*P(n1,1,ipt,pid)*Q(n3,1)*Q(n2+n4+n5+n6,4)
    + 6.*S(n1+n3,2,ipt,pid)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*S(n1+n2+n4+n5+n6,5,ipt,pid)
    - 6.*P(n1,1,ipt,pid)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*S(n1+n2,2,ipt,pid)*Q(n3+n4+n5+n6,4)
    + 24.*Q(n2,1)*S(n1+n3+n4+n5+n6,5,ipt,pid)+24.*P(n1,1,ipt,pid)*Q(n2+n3+n4+n5+n6,5)
    - 120.*S(n1+n2+n3+n4+n5+n6,6,ipt,pid);
    return formula;
}

TComplex FlowAnalysisWithQCumulantGenericFramework::FourDiffPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = PGapPos(n1,1,ipt,pid)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)-SGapPos(n1+n2,2,ipt,pid)*QGapPos(n3,1)*QGapPos(n4,1)-QGapPos(n2,1)*SGapPos(n1+n3,2,ipt,pid)*QGapPos(n4,1)
                    - PGapPos(n1,1,ipt,pid)*QGapPos(n2+n3,2)*QGapPos(n4,1)+2.0*SGapPos(n1+n2+n3,3,ipt,pid)*QGapPos(n4,1)-QGapPos(n2,1)*QGapPos(n3,1)*SGapPos(n1+n4,2,ipt,pid)
                    + QGapPos(n2+n3,2)*SGapPos(n1+n4,2,ipt,pid)-PGapPos(n1,1,ipt,pid)*QGapPos(n3,1)*QGapPos(n2+n4,2)+SGapPos(n1+n3,2,ipt,pid)*QGapPos(n2+n4,2)
                    + 2.0*QGapPos(n3,1)*SGapPos(n1+n2+n4,3,ipt,pid)-PGapPos(n1,1,ipt,pid)*QGapPos(n2,1)*QGapPos(n3+n4,2)+SGapPos(n1+n2,2,ipt,pid)*QGapPos(n3+n4,2)
                    + 2.0*QGapPos(n2,1)*SGapPos(n1+n3+n4,3,ipt,pid)+2.0*PGapPos(n1,1,ipt,pid)*QGapPos(n2+n3+n4,3)-6.0*SGapPos(n1+n2+n3+n4,4,ipt,pid);
  return formula;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::FourDiffNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = PGapNeg(n1,1,ipt,pid)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-SGapNeg(n1+n2,2,ipt,pid)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapNeg(n2,1)*SGapNeg(n1+n3,2,ipt,pid)*QGapNeg(n4,1)
                    - PGapNeg(n1,1,ipt,pid)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)+2.0*SGapNeg(n1+n2+n3,3,ipt,pid)*QGapNeg(n4,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*SGapNeg(n1+n4,2,ipt,pid)
                    + QGapNeg(n2+n3,2)*SGapNeg(n1+n4,2,ipt,pid)-PGapNeg(n1,1,ipt,pid)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)+SGapNeg(n1+n3,2,ipt,pid)*QGapNeg(n2+n4,2)
                    + 2.0*QGapNeg(n3,1)*SGapNeg(n1+n2+n4,3,ipt,pid)-PGapNeg(n1,1,ipt,pid)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)+SGapNeg(n1+n2,2,ipt,pid)*QGapNeg(n3+n4,2)
                    + 2.0*QGapNeg(n2,1)*SGapNeg(n1+n3+n4,3,ipt,pid)+2.0*PGapNeg(n1,1,ipt,pid)*QGapNeg(n2+n3+n4,3)-6.0*SGapNeg(n1+n2+n3+n4,4,ipt,pid);
  return formula;
}

// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::EightDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = FourDiffPos(n1,n2,n3,n4,ipt,pid)*FourNeg(n5,n6,n7,n8);
  return formula;
}
// ============================================================================
TComplex FlowAnalysisWithQCumulantGenericFramework::EightDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8, const Int_t ipt, const Int_t pid) const
{
  TComplex formula = FourDiffNeg(n1,n2,n3,n4,ipt,pid)*FourPos(n5,n6,n7,n8);
  return formula;
}