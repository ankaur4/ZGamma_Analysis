#!/bin/tcsh 
setenv pwd $PWD
#setenv leptonjetsDir /eos/uscms/store/user/leptonjets
#setenv leptonjetsDir /eos/uscms/store/group/leptonjets/ankaur/Ntuples_UL17/SingleMuon
foreach k ( Era_B Era_C Era_D Era_E Era_F )

#setenv sourceDir1 ${DatanTuplesDir1}/${k}/
#setenv Path1 ${Source1}/${k}/

#ankaur changes
#setenv DatanTuplesDir /uscms_data/d3/ankaur/ZGAMMA_2021/Analyser_Data_ZGamma_2017/CMSSW_10_2_10/src/Data
#setenv Source /uscms_data/d3/ankaur/ZGAMMA_2021/Analyser_Data_ZGamma_2017/CMSSW_10_2_10/src/Data


setenv DatanTuplesDir1 root://cmsxrootd.fnal.gov//store/group/leptonjets/ankaur/ZGamma/Ntuples_UL17/Data_2017UL 
setenv Source1 /store/group/leptonjets/ankaur/ZGamma/Ntuples_UL17/Data_2017UL
#setenv SinglePho SinglePhoton

foreach case ( EG )
set sampleIndex = 1

# *****************For DoublEGamma dataset*********
if ( ${case} == EG ) then
foreach i ( DoubleEGamma_v2_L1ECAL )
setenv sourceDir1 ${DatanTuplesDir1}/${i}/${k}/
setenv Path1 ${Source1}/${i}/${k}/
endif

#************** For SinglePhoton dataset***********
if ( ${case} == SP ) then
foreach i ( SinglePhoton_v2_L1ECAL ) 
setenv sourceDir1 ${DatanTuplesDir1}/${i}/${k}/
setenv Path1 ${Source1}/${i}/${k}/
endif

##*********************Making Text files for input ntuples**********
if ( -f Dataset1_${i}_${k}.txt ) then
echo "++++++++++++++ Deleting Dataset1_${i}_${k}.txt ++++++++++++++"
rm Dataset1_${i}_${k}.txt
endif


#@sampleIndex = ${sampleIndex} + 1
eosls ${Path1} > Dataset1_${i}_${k}.txt                  
setenv dataset1 Dataset1_${i}_${k}.txt


#************* Splitting the no. of files per job*******************
#Total no. of files per jobset r = 100   #change r only to change files per job
set r = 10
echo "Max Files per Job = ${r}"

set sf = 1   ## starting file from where to start

set maxf = ${r}

set Total_files =  `wc -l ${dataset1} | cut -c1-4`
#else
#set Total_files =  `wc -l ${dataset1} | cut -c1-2`
#endif

echo "Total number of files in Dataset1 ${k} are ${Total_files}" 
while ( ${sf} <= ${Total_files} )

if ( ${maxf} <= ${Total_files} ) then
setenv filenameTag1 ${i}_${k}_${sf}_${maxf}
else
setenv filenameTag1 ${i}_${k}_${sf}_${Total_files}
endif

##**************************************************
#while ( ${sf} <= ${Total_files} )


setenv destinationDir ${sourceDir1}


if ( -f PostAnalyzer_Data.C ) then
echo "++++++++++++++ Deleting PostAnalyzer_Data.C ++++++++++++++"
#rm PostAnalyzer_Data.C
endif

if ( -f PostAnalyzer_Data.h ) then
echo "++++++++++++++ Deleting PostAnalyzer_Data.h ++++++++++++++"
#rm PostAnalyzer_Data.h
endif

echo "Filename1= ${filenameTag1}"
echo "Source Dir1 = ${sourceDir1}"

#*************************************************
########## Making PostAnalyzer_Data.C ############
#*************************************************

cat > PostAnalyzer_Data.C <<EOF
#define PostAnalyzer_Data_cxx
#include "PostAnalyzer_Data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
void PostAnalyzer_Data::Loop()
{
   if (fChain == 0) return;

//***********************Defining various cuts here************
   
   Cut_Photon_pt = 15; 
   Cut_Photon_eta = 2.5; 
   Cut_Electron_eta = 2.5; 
   Cut_Electron_pt_Lead = 25;
   Cut_Electron_pt_SubLead = 25; 

//**************For choosing ID's*******************************
   //Cut_PhId = "tight";
   Cut_EleId = "medium"; //ZGamma
    //if (Pass_HLT_Ele ) {cout<<"Pass_HLT_Ele : "<<Pass_HLT_Ele<<endl;}


//****************Triggers using**************************

   SinglePhotonTriggers.push_back("HLT_Photon200_v");  //unprescaled, //ZGamma
   SinglePhotonTriggers.push_back("HLT_DoublePhoton70_v") ;   //ZGamma 2016-AN
   SinglePhotonTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");//2018 ZGamma Kyungwook

//triggers used in the study (efficiency of these three triggers is same as using all triggers)
   vector<ULong64_t> HLT_Pho;
   vector<ULong64_t> HLT_Ele;
   HLT_Ele.clear();
   HLT_Pho.clear();
   HLT_Ele.push_back(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v); 
   HLT_Pho.push_back(HLT_Photon200_v);
   //triggers for the denominator of trigger turn on
   vector<ULong64_t> HLT_deno;
   HLT_deno.clear();

//**************Define Output file here*********************
   TString OutputFile = "${filenameTag1}";
  //file = new TFile(OutputFile+".root", "RECREATE");
  
//******************Reduced************************
   Long64_t maxsize = 500000000;                 
   maxsize *= 100;                                
   TTree::SetMaxTreeSize(maxsize);

   string outFile ="${filenameTag1}";
   TFile *rfile = new TFile((outFile+"${i}.root").c_str(),"RECREATE");

   fChain->LoadTree(0);
   TTree *newtree = fChain->GetTree()->CloneTree(0);
   rfile->cd(); 

//************Define Histograms here*************************
   BookHistograms();
//************************Duplicate Event check*********************
  // map<int, set<int> > doublechecker;

//****************Event For loop starts from here *********************
   Long64_t nentries = fChain->GetEntries();
   //cout << "<Total entries" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;


   Int_t cachesize = 200000000;
   fChain->SetCacheSize(cachesize); 
   fChain->AddBranchToCache("*");
 
   // for (Long64_t jentry=0; jentry< 1000;jentry++) {
   for (Long64_t jentry=0; jentry< nentries;jentry++) {
   //cout << "*********Analyzing entry:*********" << jentry << endl;
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);
   nbytes += nb;



//***************Duplicate event removal from SinglePhoton Dataset***************************************
   
// Bool_t isitduplicate = read_from_file(run, event, lumis);
// if (isitduplicate) continue;

//*****************Initializing various parameters here************
   
   PC = -1;
   isZ = 0;
   EC1 = -1;
   EC2 = -1;
   GoodVertex = 0;  
   Pass_AllSelections = false;
//**********Initially making all bools to false********************
      
   double Zmass = 0;
   double ZGmass = 0;
   Pass_HLT = false;
   Pass_HLT_Ele = false; 
   Pass_HLT_Pho = false;
   HasPrimaryVtx = false;
   Pass_PhoPt = false;
   Pass_miniIso  = false;
   Pass_dR = false;
   Pass_dR_ll = false;
   Pass_dR1 =false;
   Pass_dR2 =false;
   Pass_PtByZGMiso = false;
   Pass_L1ECALPrefire = false;

//****************Running different functions *********************
  //Pass_HLT = PassHLT(HLT);
   Pass_HLT_Ele = PassHLT_Ele(HLT_Ele); 
   Pass_HLT_Pho = PassHLT_Pho(HLT_Pho);

   GoodVertex = nGoodVtx;
   if(GoodVertex > 0) HasPrimaryVtx = true;
  
   if(L1ECALPrefire < 1) Pass_L1ECALPrefire = true;
 
//Clear vector for next event
   GoodIsoPhotons.clear();
   GoodEtaPhotons.clear();
   GoodPtPhotons.clear();
   GoodIsoElectrons.clear();
   GoodEtaElectrons.clear();  //ZGamma
   GoodEtaPassElectrons.clear();
   GoodEtaElectrons = GoodElectrons(Cut_EleId);//ZGamma


//********************Electrons after passing Ele Id***********************

if (GoodEtaElectrons.size() != 0) {
    for (int ip = 0; ip < GoodEtaElectrons.size(); ip++) {
        if ((fabs((*eleEta)[GoodEtaElectrons[ip]]) < Cut_Electron_eta) && 
            (fabs((*eleEta)[GoodEtaElectrons[ip]]) < 1.444 || fabs((*eleEta)[GoodEtaElectrons[ip]]) > 1.566)) {
            GoodEtaPassElectrons.push_back(GoodEtaElectrons[ip]);
        }
    }
}
//(GoodEtaElectrons.size() check  ends here


if (GoodEtaPassElectrons.size() > 1) {
   for (int i = 0; i < GoodEtaPassElectrons.size(); ++i) {
   int electronIndex = GoodEtaPassElectrons[i];
   double electronPt = (*elePt)[electronIndex];
    
   if (electronPt > Cut_Electron_pt_Lead && EC1 == -1) {
      EC1 = electronIndex; // Set as leading electron if not already set
    } 
   else if (EC1 != -1 && EC2 == -1 && electronIndex != EC1) {
   double dR = GetdR((*eleEta)[EC1], (*eleEta)[electronIndex], (*elePhi)[EC1], (*elePhi)[electronIndex]);
   if (electronPt > Cut_Electron_pt_SubLead && dR > 0.1) {
      EC2 = electronIndex; // Set as subleading electron if dR condition is met
      break;
   }
 }
}
   if (EC1 > -1 && EC2 > -1 && (*eleCharge)[EC1] * (*eleCharge)[EC2] < 0 &&
        (GetZMass(EC1, EC2)) > 60.0 && (GetZMass(EC1, EC2)) < 120.0) {
        isZ = 1;
    }
}


//******** Photon selection********************************* 
//Loop through all photons
  
   int numPhotonsPerEvent = 0;


   for (Int_t i = 0; i < nPho; i++) {
          numPhotonsPerEvent++;      
           if ((((*phoIDbit)[i] >> 1) & 1) &&((*phoEleVeto)[i] == 1))
        {
          GoodIsoPhotons.push_back(i);
    }
  }
// h_PhotonCounts[1]->Fill(numPhotonsPerEvent);
 

      

if (GoodIsoPhotons.size() != 0) {
    for (int ip = 0; ip < GoodIsoPhotons.size(); ip++) {
       int photonIndex = GoodIsoPhotons[ip];
       double photonEta = fabs((*phoSCEta)[photonIndex]);
       
       if (photonEta < Cut_Photon_eta && (photonEta < 1.444 || photonEta > 1.566)) {
            GoodEtaPhotons.push_back(photonIndex);
        }
    
     }
    
 }
//h_PhotonCounts[2]->Fill(numPhotonsPerEvent);

// Check if the leading photon passes the pt, eta, and dR selections with both electrons

if (GoodEtaPhotons.size() != 0) {
      for (int j = 0; j < GoodEtaPhotons.size(); j++) {
            if ((*phoEt)[GoodEtaPhotons[j]] > Cut_Photon_pt ){
               GoodPtPhotons.push_back(GoodEtaPhotons[j]);
            }
        }

}


if (GoodPtPhotons.size() !=0 && EC1 > -1 && EC2 > -1){
       for (int l = 0; l < GoodPtPhotons.size(); l++){
        if(((GetdR((*phoSCEta)[GoodPtPhotons[l]], (*eleEta)[EC1], (*phoSCPhi)[GoodPtPhotons[l]], (*elePhi)[EC1])) > 0.4) &&
                ((GetdR((*phoSCEta)[GoodPtPhotons[l]], (*eleEta)[EC2], (*phoSCPhi)[GoodPtPhotons[l]], (*elePhi)[EC2])) > 0.4)) {
                PC = GoodPtPhotons[l];
                break;
// Exit the loop after finding the first photon that passes the selection   
           }
      }
 }

//h_PhotonCounts[3]->Fill(numPhotonsPerEvent);

    //Check on GoodEtaPhotons end here

//*******ZGAMMA Quantities and Kinematics******************************************************************************
             
        if (EC1 > -1 && EC2 > -1  && isZ ) Zmass = GetZMass(EC1,EC2);
        if (EC1 > -1 && EC2 > -1  && PC > -1 && isZ ) ZGmass = GetZGMass(EC1,EC2,PC);
        if ( EC1 > -1 && EC2 > -1 && isZ  && PC > -1 ) Pass_PtByZGMiso =(((*phoEt)[PC]/GetZGMass(EC1,EC2,PC)) > 0.115384615 ); 
            // if (EC1 > -1 ) cout <<"elePhi_EC1:"<< (*elePhi)[EC1]<<endl;
            // if (EC2 > -1 ) cout <<"elePhi_EC2:"<< (*elePhi)[EC2]<<endl;           
      

      h_CutFlow_ZGamma_Ele->Fill(0.5); //ZGamma
      h_event->Fill(event);
     //Primary vertex info for noCut only
//      h_trueNumofInt       ->Fill((*puTrue)[0]);
        h_goodPV_noWt        ->Fill(GoodVertex);
//      h_goodPV_LumiWt      ->Fill(GoodVertex, Lumi_EvtWt);
//      h_goodPV_PUWt        ->Fill(GoodVertex, PU_EvtWt);                           
   //------------------------------------------------------------
      //Photon distributions noCut
        if(nPho>0  &  nEle > 1 )
        {
	h_PhotonPt[0]               ->Fill((*phoEt)[0]);
        h_PhotonPtSublead[0]        ->Fill((*phoEt)[1]);
	h_PhotonCalibPt[0]          ->Fill((*phoCalibEt)[0]);
	h_PhotonEta[0]              ->Fill((*phoSCEta)[0]);
	h_PhotonPhi[0]              ->Fill((*phoSCPhi)[0]);
	h_Photon_SigmaIEtaIEta[0]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[0]);
	h_Photon_R9[0]              ->Fill((*phoR9)[0]);
	h_Photon_HoverE[0]          ->Fill((*phoHoverE)[0]);
	h_Photon_EleVeto[0]         ->Fill((*phoEleVeto)[0]);
	h_Photon_CorrPFChIso[0]     ->Fill(TMath::Max(((*phoPFChIso)[0] - rho*EAChargedHadrons((*phoSCEta)[0])), 0.0));
	h_Photon_CorrPFNeuIso[0]    ->Fill(TMath::Max(((*phoPFNeuIso)[0] - rho*EANeutralHadrons((*phoSCEta)[0])), 0.0));
	h_Photon_CorrPFPhoIso[0]    ->Fill(TMath::Max(((*phoPFPhoIso)[0] - rho*EAPhotons((*phoSCEta)[0])), 0.0));
        h_Photon_MVAID[0]           ->Fill((*phoIDMVA)[0]);		      
      //leptons distribution noCut
        h_leptonPt_Lead[0]               ->Fill((*elePt)[0]);
        h_leptonEta_Lead[0]              ->Fill((*eleEta)[0]);
        h_leptonSCEta_Lead[1]           ->Fill((*eleSCEta)[0]);
        h_leptondEtaAtVtx_Lead[0]        ->Fill((*eledEtaAtVtx)[0]);
        h_leptondPhiAtVtx_Lead[0]        ->Fill((*eledPhiAtVtx)[0]);
        h_leptonD0_Lead[0]               ->Fill((*eleD0)[0]);
        h_leptonDz_Lead[0]               ->Fill((*eleDz)[0]); 
        h_leptonPhi_Lead[0]              ->Fill((*elePhi)[0]);
        h_leptonCalibPt_Lead[0]          ->Fill((*eleCalibPt)[0]);
        h_lepton_R9_Lead[0]              ->Fill((*eleR9)[0]);
        h_lepton_HoverE_Lead[0]          ->Fill((*eleHoverE)[0]);
        h_lepton_SigmaIEtaIEta_Lead[0]   ->Fill((*eleSigmaIEtaIEtaFull5x5)[0]);
        h_lepton_charge_Lead[0]          ->Fill((*eleCharge)[0]);
        h_lepton_PFMiniIso_Lead[0]       ->Fill((*elePFMiniIso)[0]);
        h_lepton_ConvVeto_Lead[0]        ->Fill((*eleConvVeto)[0]);
        h_lepton_MiniIsoByPt_Lead[0]     ->Fill((*elePFMiniIso)[0]/(*elePt)[0]);
        h_lepton_MissHits_Lead[0]         ->Fill((*eleMissHits)[0]);
        h_lepton_EoverPInv_Lead[0]            ->Fill((*eleEoverPInv)[0]);
        
        h_leptonPt_SubLead[0]               ->Fill((*elePt)[1]); 
        h_leptonEta_SubLead[0]              ->Fill((*eleEta)[1]);
        h_leptonSCEta_SubLead[0]            ->Fill((*eleSCEta)[1]);
        h_leptondEtaAtVtx_SubLead[0]        ->Fill((*eledEtaAtVtx)[1]);
        h_leptondPhiAtVtx_SubLead[0]        ->Fill((*eledPhiAtVtx)[1]);
        h_leptonD0_SubLead[0]               ->Fill((*eleD0)[1]);
        h_leptonDz_SubLead[0]               ->Fill((*eleDz)[1]);
        h_leptonPhi_SubLead[0]              ->Fill((*elePhi)[1]);
        h_leptonCalibPt_SubLead[0]          ->Fill((*eleCalibPt)[1]);
        h_lepton_R9_SubLead[0]              ->Fill((*eleR9)[1]);
        h_lepton_HoverE_SubLead[0]          ->Fill((*eleHoverE)[1]);
        h_lepton_SigmaIEtaIEta_SubLead[0]   ->Fill((*eleSigmaIEtaIEtaFull5x5)[1]);
        h_lepton_charge_SubLead[0]          ->Fill((*eleCharge)[1]);
        h_lepton_PFMiniIso_SubLead[0]       ->Fill((*elePFMiniIso)[1]);
        h_lepton_ConvVeto_SubLead[0]        ->Fill((*eleConvVeto)[1]);
        h_lepton_MiniIsoByPt_SubLead[0]     ->Fill((*elePFMiniIso)[1]/(*elePt)[1]);
        h_lepton_MissHits_SubLead[0]         ->Fill((*eleMissHits)[1]);
        h_lepton_EoverPInv_SubLead[0]        ->Fill((*eleEoverPInv)[1]);       

        h_lepton_dPhi[0]                     ->Fill(GetdPhi((*elePhi)[0],(*elePhi)[1])); 
        h_lepton_dR[0]                      ->Fill(GetdR((*eleEta)[EC1],(*eleEta)[EC2],(*elePhi)[EC1],(*elePhi)[EC2]));
        h_lepton_dR1[0]              ->Fill(GetdR((*phoSCEta)[0],(*eleEta)[0],(*phoSCPhi)[0],(*elePhi)[0]));
        h_lepton_dR2[0]              ->Fill(GetdR((*phoSCEta)[0],(*eleEta)[1],(*phoSCPhi)[0],(*elePhi)[1]));
        h_Zmass[0]->Fill(GetZMass(0,1));
        h_ZGmass[0]->Fill(GetZGMass(0,1,0));


      //Primary vertex and number of photon and jets for noCut
	h_goodPV[0]                      ->Fill(GoodVertex);
	h_nIsoPhotons[0]                 ->Fill(GoodIsoPhotons.size());  // Tot # of isolated photons
        h_nGoodPhotons[0]                ->Fill(GoodIsoBarrelPhotons.size()); // Tot # of isolated photons with pt > cut and eta < cut 
        for(int ip = 0; ip < GoodIsoPhotons.size(); ip++){
	  h_IsoPhotonIdxVsPt[0]          ->Fill((*phoEt)[GoodIsoPhotons[ip]], ip+1);
         }
	for(int ii = 0; ii < GoodIsoBarrelPhotons.size(); ii++){
	  h_GoodPhotonIdxVsPt[0]         ->Fill((*phoEt)[GoodIsoBarrelPhotons[ii]], ii+1);
	}				     
      h_dRVsZMass[0]->Fill(GetdR((*eleEta)[EC1],(*eleEta)[EC2],(*elePhi)[EC1],(*elePhi)[EC2]),GetZMass(EC1,EC2)); 
      h_PtByZGMassVsZGMass[0]->Fill(GetZGMass(EC1,EC2,PC),((*phoEt)[PC]/GetZGMass(EC1,EC2,PC)));
      }//nPho> nele>0 check ends here

   // ------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//                       ZGAMMA_Ele
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if(Pass_HLT_Ele){
        h_CutFlow_ZGamma_Ele->Fill(1.5);
       
        // cout<< "Pass_HLT_Ele :" <<Pass_HLT_Ele<<endl;
        if(HasPrimaryVtx ) { 
          h_CutFlow_ZGamma_Ele->Fill(2.5);   
        
         if(GoodEtaElectrons.size() > 0 ) {
              h_CutFlow_ZGamma_Ele->Fill(3.5);
              h_GoodEtaElectrons[0]->Fill(GoodEtaElectrons.size());
            
              h_lepton_dR[1]      ->Fill(GetdR((*eleEta)[GoodEtaElectrons[0]],(*eleEta)[GoodEtaElectrons[1]],(*elePhi)[GoodEtaElectrons[0]],(*elePhi)[GoodEtaElectrons[1]]));     
     
             if(GoodEtaPassElectrons.size() > 1 ){
                h_CutFlow_ZGamma_Ele->Fill(4.5);
                 h_GoodEtaPassElectrons[0]->Fill(GoodEtaPassElectrons.size());

                if(EC1 > -1 && EC2 > -1) {
                  h_CutFlow_ZGamma_Ele->Fill(5.5);
          //          h_lepton_dR[2]       ->Fill(GetdR((*eleEta)[EC1],(*eleEta)[EC2],(*elePhi)[EC1],(*elePhi)[EC2])); 
              
                      if( isZ  ){ 
                       h_CutFlow_ZGamma_Ele->Fill(6.5);
                       h_Zmass[1]->Fill(GetZMass(EC1,EC2));
                       h_goodPV[1]  ->Fill(GoodVertex);        
                       h_dRVsZMass[1]->Fill(GetdR((*eleEta)[EC1],(*eleEta)[EC2],(*elePhi)[EC1],(*elePhi)[EC2]),GetZMass(EC1,EC2)); 
    
                          if(GoodIsoPhotons.size() > 0 ){
                           h_CutFlow_ZGamma_Ele->Fill(7.5);
                             h_PhotonCounts[0]->Fill(numPhotonsPerEvent); 
                             h_GoodIsoPhotons[0]->Fill(GoodIsoPhotons.size());
                             h_lepton_dR1[1]              ->Fill(GetdR((*phoSCEta)[GoodIsoPhotons[0]],(*eleEta)[EC1],(*phoSCPhi)[GoodIsoPhotons[0]],(*elePhi)[EC1]));
                             h_lepton_dR2[1]              ->Fill(GetdR((*phoSCEta)[GoodIsoPhotons[0]],(*eleEta)[EC2],(*phoSCPhi)[GoodIsoPhotons[0]],(*elePhi)[EC2]));


                            if(GoodEtaPhotons.size() > 0){
                             h_CutFlow_ZGamma_Ele->Fill(8.5);
                               h_PhotonCounts[1]->Fill(numPhotonsPerEvent);                                 
                                h_GoodEtaPhotons[0]->Fill(GoodEtaPhotons.size());

                            if(GoodPtPhotons.size() > 0){
                                h_CutFlow_ZGamma_Ele->Fill(9.5);
                                   h_PhotonCounts[2]->Fill(numPhotonsPerEvent);
                                   h_GoodPtPhotons[0]->Fill(GoodPtPhotons.size());              

                             if( PC > -1){
                               h_CutFlow_ZGamma_Ele->Fill(10.5);        
                         //      h_lepton_dR1[2]         ->Fill(GetdR((*phoSCEta)[PC],(*eleEta)[EC1],(*phoSCPhi)[PC],(*elePhi)[EC1]));
                         //      h_lepton_dR2[2]         ->Fill(GetdR((*phoSCEta)[PC],(*eleEta)[EC2],(*phoSCPhi)[PC],(*elePhi)[EC2]));
                               h_Zmass[2]   ->Fill(GetZMass(EC1,EC2));
                               h_ZGmass[1]   ->Fill(GetZGMass(EC1,EC2,PC));
                               h_PhotonCounts[3]->Fill(numPhotonsPerEvent);
                               /*
                               if((GetZGMass(EC1,EC2,PC)) > 130)
                                        {
                                         h_CutFlow_ZGamma_Ele->Fill(11.5);
                                   */
                                    

                                       h_PtByZGMiso[0]->Fill((*phoEt)[PC] / GetZGMass(EC1, EC2, PC));
                                       if(Pass_PtByZGMiso){
                                         Pass_AllSelections = true;
                                         h_CutFlow_ZGamma_Ele->Fill(12.5);
                                   
//    cout<<"run "<<run<<"event "<<event<<" (*phoEt)[PC] = "<<(*phoEt)[PC]<<"(*phoSCEta)[PC] = "<<(*phoSCEta)[PC]<<endl;
//    cout<<"(*elePt)[EC1] = "<<(*elePt)[EC1]<<"(*eleEta)[EC1] = "<<(*eleEta)[EC1]<<endl;
//    cout<<"(*elePt)[EC2] = "<<(*elePt)[EC2]<<"(*eleEta)[EC2] = "<<(*eleEta)[EC2]<<endl;                           
   
                                 
                     //cout<<"Pass_PtByZGMiso"<<Pass_PtByZGMiso<<"((*phoEt)[PC]/GetZGMass(EC1,EC2,PC)" <<((*phoEt)[PC]/GetZGMass(EC1,EC2,PC))<<endl;      
                                 //cout<<run<<":"<<event<<":"<<lumis<<endl; 

                h_PtByZGMiso[1]->Fill((*phoEt)[PC] / GetZGMass(EC1, EC2, PC));
              //----------------------------------------------------------
              //  [1]
                        h_lepton_dR[2]       ->Fill(GetdR((*eleEta)[EC1],(*eleEta)[EC2],(*elePhi)[EC1],(*elePhi)[EC2]));
                        h_lepton_dR1[2]         ->Fill(GetdR((*phoSCEta)[PC],(*eleEta)[EC1],(*phoSCPhi)[PC],(*elePhi)[EC1]));
                        h_lepton_dR2[2]         ->Fill(GetdR((*phoSCEta)[PC],(*eleEta)[EC2],(*phoSCPhi)[PC],(*elePhi)[EC2]));
                        h_GoodEtaElectrons[1]   ->Fill(GoodEtaElectrons.size());
                        h_GoodEtaPassElectrons[1]->Fill(GoodEtaPassElectrons.size());   
                        h_GoodIsoPhotons[1]     ->Fill(GoodIsoPhotons.size());
                        h_GoodEtaPhotons[1]     ->Fill(GoodEtaPhotons.size());
                        h_GoodPtPhotons[1]      ->Fill(GoodPtPhotons.size());
                 
                        h_goodPV[2]          ->Fill(GoodVertex);
                        h_PhotonPt[1]               ->Fill((*phoEt)[PC]);
                        h_PhotonCalibPt[1]          ->Fill((*phoCalibEt)[PC]);
                        h_PhotonEta[1]              ->Fill((*phoSCEta)[PC]);
                        h_PhotonPhi[1]              ->Fill((*phoSCPhi)[PC]);
                        h_Photon_SigmaIEtaIEta[1]   ->Fill((*phoSigmaIEtaIEtaFull5x5)[PC]);
                        h_Photon_R9[1]              ->Fill((*phoR9)[PC]);
                        h_Photon_HoverE[1]          ->Fill((*phoHoverE)[PC]);
                        h_Photon_EleVeto[1]         ->Fill((*phoEleVeto)[PC]);
                        h_Photon_CorrPFChIso[1]     ->Fill(TMath::Max(((*phoPFChIso)[PC] - rho*EAChargedHadrons((*phoSCEta)[PC])), 0.0));
                        h_Photon_CorrPFNeuIso[1]    ->Fill(TMath::Max(((*phoPFNeuIso)[PC] - rho*EANeutralHadrons((*phoSCEta)[PC])), 0.0));
                        h_Photon_CorrPFPhoIso[1]    ->Fill(TMath::Max(((*phoPFPhoIso)[PC] - rho*EAPhotons((*phoSCEta)[PC])), 0.0));
                        h_Photon_MVAID[1]           ->Fill((*phoIDMVA)[PC]);
      
                        h_leptonPt_Lead[1]               ->Fill((*elePt)[EC1]);
                        h_leptonEta_Lead[1]              ->Fill((*eleEta)[EC1]);
                        h_leptonSCEta_Lead[1]            ->Fill((*eleSCEta)[EC1]);
                        h_leptonPhi_Lead[1]              ->Fill((*elePhi)[EC1]);
                        h_leptonCalibPt_Lead[1]          ->Fill((*eleCalibPt)[EC1]);
                        h_lepton_R9_Lead[1]              ->Fill((*eleR9)[EC1]);
                        h_lepton_HoverE_Lead[1]          ->Fill((*eleHoverE)[EC1]);
                        h_lepton_SigmaIEtaIEta_Lead[1]   ->Fill((*eleSigmaIEtaIEtaFull5x5)[EC1]);
                        h_lepton_charge_Lead[1]          ->Fill((*eleCharge)[EC1]);
                        h_lepton_PFMiniIso_Lead[1]       ->Fill((*elePFMiniIso)[EC1]);
                        h_lepton_ConvVeto_Lead[1]        ->Fill((*eleConvVeto)[EC1]);
                        h_lepton_MiniIsoByPt_Lead[1]     ->Fill((*elePFMiniIso)[EC1]/(*elePt)[EC1]);
                        h_lepton_MissHits_Lead[1]        ->Fill((*eleMissHits)[EC1]);
                        h_leptondEtaAtVtx_Lead[1]         ->Fill((*eledEtaAtVtx)[EC1]);
                        h_leptondPhiAtVtx_Lead[1]        ->Fill((*eledPhiAtVtx)[EC1]);
                        h_leptonD0_Lead[1]               ->Fill((*eleD0)[EC1]);
                        h_leptonDz_Lead[1]               ->Fill((*eleDz)[EC1]); 
                        h_lepton_EoverPInv_Lead[1]       ->Fill((*eleEoverPInv)[EC1]);
                        h_lepton_relIsowithEA_Lead[1]    ->Fill((*elePFChIso)[EC1] + TMath::Max(0.0, (*elePFNeuIso)[EC1] + (*elePFPhoIso)[EC1] - rho*EANeutralHadrons((*eleSCEta)[EC1] ))); 
                        h_leptonPt_SubLead[1]               ->Fill((*elePt)[EC2]);
                        h_leptonEta_SubLead[1]              ->Fill((*eleEta)[EC2]);
                        h_leptonSCEta_SubLead[1]            ->Fill((*eleSCEta)[EC2]);
                        h_leptonPhi_SubLead[1]              ->Fill((*elePhi)[EC2]);
                        h_leptonCalibPt_SubLead[1]          ->Fill((*eleCalibPt)[EC2]);
                        h_lepton_R9_SubLead[1]              ->Fill((*eleR9)[EC2]);
                        h_lepton_HoverE_SubLead[1]          ->Fill((*eleHoverE)[EC2]);
                        h_lepton_SigmaIEtaIEta_SubLead[1]   ->Fill((*eleSigmaIEtaIEtaFull5x5)[EC2]);
                        h_lepton_charge_SubLead[1]          ->Fill((*eleCharge)[EC2]);
                        h_lepton_PFMiniIso_SubLead[1]       ->Fill((*elePFMiniIso)[EC2]);
                        h_lepton_ConvVeto_SubLead[1]        ->Fill((*eleConvVeto)[EC2]);
                        h_lepton_MiniIsoByPt_SubLead[1]     ->Fill((*elePFMiniIso)[EC2]/(*elePt)[EC2]);
                        h_lepton_MissHits_SubLead[1]        ->Fill((*eleMissHits)[EC2]);
                        h_leptondEtaAtVtx_SubLead[1]        ->Fill((*eledEtaAtVtx)[EC2]);
                        h_leptondPhiAtVtx_SubLead[1]        ->Fill((*eledPhiAtVtx)[EC2]);
                        h_leptonD0_SubLead[1]               ->Fill((*eleD0)[EC2]);
                        h_leptonDz_SubLead[1]               ->Fill((*eleDz)[EC2]); 
                        h_lepton_EoverPInv_SubLead[1]       ->Fill((*eleEoverPInv)[EC2]);                       
                        h_lepton_relIsowithEA_SubLead[1]    ->Fill((*elePFChIso)[EC2] + TMath::Max(0.0, (*elePFNeuIso)[EC2] + (*elePFPhoIso)[EC2] - rho*EANeutralHadrons((*eleSCEta)[EC2] )));
                        h_lepton_dPhi[1]                    ->Fill(GetdPhi((*elePhi)[EC1],(*elePhi)[EC2]) > 2.1 );
                        h_dRVsZMass[2]->Fill(GetdR((*eleEta)[EC1],(*eleEta)[EC2],(*elePhi)[EC1],(*elePhi)[EC2]),GetZMass(EC1,EC2)); 
                        h_PtByZGMassVsZGMass[2]->Fill(GetZGMass(EC1,EC2,PC),((*phoEt)[PC]/GetZGMass(EC1,EC2,PC)));
                        h_Zmass[3]   ->Fill(GetZMass(EC1,EC2));
                        h_ZGmass[2]   ->Fill(GetZGMass(EC1,EC2,PC));
                        }
                       }
                     }
                   }
                 }   
               }
             } 
           }
          }
        }
      }
  //  }
      //---------------------------------------------------------------------------------   
if(Pass_AllSelections)newtree->Fill();
   Pass_AllSelections = false;

 }//for jentry
   rfile->cd();
 
   newtree->Write();
 
   rfile->Write();
   rfile->Close();
}//Loop()

EOF



########## Making PostAnalyzer_Data.h ############
cat > PostAnalyzer_Data.h <<EOF
/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 11 00:35:43 2015 by ROOT version 6.02/05
// from TChain ggNtuplizer/EventTree/
//////////////////////////////////////////////////////////

#ifndef PostAnalyzer_Data_h
#define PostAnalyzer_Data_h

//ROOT include files
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h>
#include <TGraphAsymmErrors.h>
//#include "TRFIOFile.h"
//#include "TKDE.h" 

// Header file for the classes stored in the TTree if any.
#include "vector"

//c++ include files
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <map>
#include <set>


using namespace std;
using namespace ROOT;

class PostAnalyzer_Data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   //-----------------------
   //Variables defined by me
   //-----------------------
   //Outut root file
   TFile *file;
   
   //Bools
   Bool_t isZ; //ZGamma
   Bool_t Pass_HLT;
   Bool_t Pass_HLT_Ele; //ZGamma
   Bool_t Pass_HLT_Pho; //ZGamma
   Bool_t HasPrimaryVtx;
   Bool_t Pass_PhoPt;
   Bool_t Pass_PhoEtaEB;
   Bool_t Pass_JetPt;
   Bool_t Pass_JetEta;
   Bool_t Pass_GJdPhi;
   Bool_t Pass_GJdEta;
   Bool_t Pass_GJInvtMass;
   Bool_t Pass_AllSelections;
   Bool_t Pass_CSVv2bTag;
   Bool_t Fail_GJdEta;
   Bool_t Pass_dR_ll;
   Bool_t Pass_dR;
   Bool_t Pass_dR1;
   Bool_t Pass_dR2;  
   Bool_t Pass_miniIso;
   Bool_t Pass_dPhi;
   Bool_t Pass_PtByZGMiso;
   Bool_t Pass_L1ECALPrefire;
   //Ints, Doubles etc.
   Int_t GoodVertex;
   Int_t PC, JC, EC1, EC2;

   Double_t Cut_Vtx_z; //(cm)
   Double_t Cut_Vtx_ndof;
   Double_t Cut_Vtx_rho; //(cm)

   Double_t Cut_Photon_pt;
   Double_t Cut_Photon_eta;
   Double_t Cut_Electron_eta; 
   Double_t Cut_Electron_pt;
   Double_t Cut_Electron_pt_Lead;//ZGamma
   Double_t Cut_Electron_pt_SubLead;//ZGamma
   Double_t Cut_Photon_pt_EB; //(GeV) ZGamma
   Double_t Cut_Photon_eta_EB;
   Double_t Cut_Photon_pt_EE; //(GeV)
   Double_t Cut_Photon_eta_EE;
   Double_t Cut_Electron_MiniIsoByPt;// ZGamma
   Double_t Cut_Jet_pt; //(GeV)
   Double_t Cut_Jet_eta;

   Double_t Cut_GJdPhi;
   Double_t Cut_GJdEta;
   Double_t Cut_GJInvtMass;
   TString Cut_PhId;
   TString Cut_JetId;
   TString Cut_EleId; //ZGamma

   //Photon,Electron  and Jet Containers
   vector<Int_t> GoodIsoPhotons;
   vector<Int_t> GoodIsoSortedPhotons;
   vector<Int_t> GoodEtaPhotons;
   vector<Int_t> GoodPtPhotons;
   vector<Int_t> GoodIsoBarrelPhotons;
   vector<Int_t> GoodIsoJets;
   vector<Int_t> GoodIsoElectrons;
   vector<Int_t> GoodEtaElectrons;
   vector<Int_t> GoodEtaPassElectrons;
   vector<Int_t> GoodIsoMuons; 

//  Trigger Info (Required for data only)
   std::vector<std::string> SinglePhotonTriggers;
   //Unprescaled Triggers (Bits assigned according to ggAnalysis/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc)
   const ULong64_t HLT_Photon175_v = 7;
   const ULong64_t HLT_Photon250_NoHE_v = 8;
   const ULong64_t HLT_Photon300_NoHE_v = 9;
   const ULong64_t HLT_Photon200_v = 10; //ZGamma
   const ULong64_t HLT_DoublePhoton70_v  = 22 ;
//   const ULong64_t HLT_Photon600_v = 11;
   const ULong64_t HLT_Photon165_HE10_v = 12;
   const ULong64_t HLT_Photon165_R9Id90_HE10_IsoM_v = 34;
   //Prescaled Triggers (Bits assigned according to ggAnalysis/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc)
   const ULong64_t HLT_Photon75_v = 4;
   const ULong64_t HLT_Photon90_v = 5;
   const ULong64_t HLT_Photon120_v = 6;
   const ULong64_t  HLT_Photon33_v  = 2;
   const ULong64_t  HLT_Photon50_v  = 3;
   const ULong64_t  HLT_Photon30_v  = 1;
   const ULong64_t  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = 40; //ZGamma   
   const ULong64_t HLT_PFJet60_v   = 11;
  // const ULong64_t HLT_PFJet80_v   = 12; 
   const ULong64_t HLT_PFJet140_v   = 13; 
   //*****************************************************************
   //Histograms
   TH1F *h_PhotonPt[7];
   TH1F *h_PhotonPtSublead[7];
   TH1F *h_PhotonCalibPt[7];
   TH1F *h_PhotonEta[7];
   TH1F *h_PhotonPhi[7];
   TH1F *h_Photon_SigmaIEtaIEta[7];
   TH1F *h_Photon_R9[7];
   TH1F *h_Photon_HoverE[7];
   TH1F *h_Photon_EleVeto[7];
   TH1F *h_Photon_CorrPFChIso[7];
   TH1F *h_Photon_CorrPFNeuIso[7];
   TH1F *h_Photon_CorrPFPhoIso[7];
   TH1F *h_Photon_MVAID[7];
   TH1F *h_PhotonsPerEvent;
   TH1F *h_PhotonCounts[7];

   TH1F *h_leptonPt_Lead[7];
   TH1F *h_leptonCalibPt_Lead[7];
   TH1F *h_leptonEta_Lead[7];
   TH1F *h_leptonSCEta_Lead[7];
   TH1F *h_leptonPhi_Lead[7];
   TH1F *h_lepton_SigmaIEtaIEta_Lead[7];
   TH1F *h_lepton_R9_Lead[7];
   TH1F *h_lepton_HoverE_Lead[7];
   TH1F *h_lepton_charge_Lead[7];   
   TH1F *h_lepton_PFMiniIso_Lead[7];
   TH1F *h_lepton_ConvVeto_Lead[7];
   TH1F *h_lepton_MiniIsoByPt_Lead[7];
   TH1F *h_lepton_MissHits_Lead[7];
   TH1F *h_leptondEtaAtVtx_Lead[7];
   TH1F *h_leptondPhiAtVtx_Lead[7];
   TH1F *h_leptonD0_Lead[7];
   TH1F *h_leptonDz_Lead[7];   
   TH1F *h_lepton_EoverPInv_Lead[7];           
   TH1F *h_lepton_relIsowithEA_Lead[7];

   TH1F *h_leptonPt_SubLead[7];
   TH1F *h_leptonCalibPt_SubLead[7];
   TH1F *h_leptonEta_SubLead[7];
   TH1F *h_leptonSCEta_SubLead[7];
   TH1F *h_leptonPhi_SubLead[7];
   TH1F *h_lepton_SigmaIEtaIEta_SubLead[7];
   TH1F *h_lepton_R9_SubLead[7];
   TH1F *h_lepton_HoverE_SubLead[7];
   TH1F *h_lepton_charge_SubLead[7];   
   TH1F *h_lepton_PFMiniIso_SubLead[7];
   TH1F *h_lepton_ConvVeto_SubLead[7];
   TH1F *h_lepton_MiniIsoByPt_SubLead[7];
   TH1F *h_lepton_MissHits_SubLead[7];
   TH1F *h_leptondEtaAtVtx_SubLead[7];
   TH1F *h_leptondPhiAtVtx_SubLead[7];
   TH1F *h_leptonD0_SubLead[7];
   TH1F *h_leptonDz_SubLead[7];
   TH1F *h_lepton_EoverPInv_SubLead[7];  
   TH1F *h_lepton_relIsowithEA_SubLead[7];
   TH1F *h_lepton_dPhi[7];
   TH1F *h_lepton_dR[7];
   TH1F *h_lepton_dR1[7];
   TH1F *h_lepton_dR2[7]; 
   TH1F *h_GoodEtaElectrons[7];
   TH1F *h_GoodEtaPassElectrons[7];
   TH1F *h_GoodIsoPhotons[7];
   TH1F *h_GoodEtaPhotons[7];
   TH1F *h_GoodPtPhotons[7];
   TH1F *h_PtByZGMiso[7];
   TH1F  *h_Zmass[7];
   TH1F  *h_ZGmass[7];
   TH1F  *h_goodPV_noWt;
   TH1F *h_goodPV_WithCut;
   TH1F *h_TotalEvents;
   TH1F *h_bJetPt[7];
   TH1F *h_bJetEta[7];
   TH1F *h_bJetPhi[7];
   TH1F *h_bJet_Mt[7];
   TH1F *h_bJet_area[7];
   TH1F *h_bJet_Mass[7];
   TH1F *h_bJet_NHEF[7]; //Neutral Hadron energy fraction
   TH1F *h_bJet_NEEF[7]; //Neutral EM energy fraction
   TH1F *h_bJet_NConst[7]; //Number of constituents
   TH1F *h_bJet_CHEF[7];  //Charged Hadron energy fraction
   TH1F *h_bJet_ChMult[7]; //Charged Multiplicity
   TH1F *h_bJet_CEEF[7]; //Charged EM energy fraction
   TH1F *h_bJet_MUF[7]; //Muon energy fraction
   TH1F *h_bJet_NNP[7]; //Number of neutral particles

   TH1F *h_GbJetInvtMass_VarBin[7];
   TH1F *h_GbJetInvtMass_UnitBin[7];
   TH1F *h_GbJet_dEta[7];
   TH1F *h_GbJet_dPhi[7];
   TH1F *h_GbJet_dR[7];
   TH1F *h_cosThetaStar[7];

   TH1F *h_PFMet[7];
   TH1F *h_SumPFMet[7];
   TH1F *h_MetBySumMET[7];
   TH2F *h_PFMetVsGJmass[7];
   TH2F *h_PFMetOverSumEtVsGJmass[7];
   TH2F *h_PtByZGMassVsZGMass[7];
   TH2F *h_dRVsZMass[7];
   TH1F *h_MetByPhPt[7];

   TH2F *h_PhPt_vs_bJetPt[7];
   TH2F *h_PhEta_vs_bJetEta[7];

   TH1F *h_CSVv2Dist[7];
   TH2F *h_CSVv2_vs_bJetPt[7];
   TH2F *h_CSVv2_vs_bJetEta[7];

   //Pileup
   TH1F *h_goodPV[7];
   TH1F *h_event;

   //Number of photons and jets
   TH1F *h_nIsoPhotons[7];
   TH1F *h_nGoodPhotons[7];
   TH2F *h_IsoPhotonIdxVsPt[7];
   TH2F *h_GoodPhotonIdxVsPt[7];
   TH1F *h_nJets[7];
   TH2F *h_JetIdxVsPt[7];

   //Trigger Turn-on
   TH1F *h_TrigPhotonPt[2];
   TH1F *h_TrigEff;
   TH1F *h_TrigGJmass_Jetpt210;
   TH1F *h_TrigGJmass_Jetpt220;
   TH1F *h_TrigGJmass_Jetpt200;
   TH1F *h_TrigGJmass_Jetpt170;
   TH1F *h_TrigGJmass_Jetpt190;
   TH1F *h_LeadPt;

   TH1F *h_PC;
   TH1F *h_JC;

   //Cut Flow
   TH1F *h_CutFlow_qstar;
   TH1F *h_CutFlow_bstar;
   TH1F *h_CutFlow_ZGamma_Ele;   //ZGamma

/*   //Trigg Eff
   TH1F *h_PassHLT175_EBPassProbes;
   TH1F *h_PassHLT90OR120_EBAllProbes;
   TH1F *h_PassHLT175_EEPassProbes;
   TH1F *h_PassHLT90OR120_EEAllProbes;         */
   //***********************************************************

   //Variables of TTree::EventTree
   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Int_t           nTrksPV;
   Bool_t          isPVGood;
   Float_t         L1ECALPrefire;            
   Float_t         vtx;
   Float_t         vty;
   Float_t         vtz;
   Float_t         rho;
   Float_t         rhoCentral;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTPho;
   ULong64_t       HLTJet;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTJetIsPrescaled;
   vector<int>     *phoPrescale;
   Int_t           metFilters;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Float_t         pfMETPhi_T1JESUp;
   Float_t         pfMETPhi_T1JESDo;
   Float_t         pfMETPhi_T1UESUp;
   Float_t         pfMETPhi_T1UESDo;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoPx;
   vector<float>   *phoPy;
   vector<float>   *phoPz;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibEt;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEn;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoE1x3;
   vector<float>   *phoE1x5;
   vector<float>   *phoE2x2;
   vector<float>   *phoE2x5Max;
   vector<float>   *phoE5x5;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phomaxXtalenergyFull5x5;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoE1x3Full5x5;
   vector<float>   *phoE1x5Full5x5;
   vector<float>   *phoE2x2Full5x5;
   vector<float>   *phoE2x5MaxFull5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoSeedBCE;
   vector<float>   *phoSeedBCEta;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoCITKChIso;
   vector<float>   *phoCITKPhoIso;
   vector<float>   *phoCITKNeuIso;
   vector<float>   *phoIDMVA;
   vector<unsigned int> *phoFiredSingleTrgs;
   vector<unsigned int> *phoFiredDoubleTrgs;
   vector<unsigned int> *phoFiredL1Trgs;
//   vector<float>   *phoSeedTime;
 //  vector<float>   *phoSeedEnergy;
   vector<unsigned short> *phoxtalBits;
   vector<unsigned short> *phoIDbit;
   Int_t           npfPho;
   vector<float>   *pfphoEt;
   vector<float>   *pfphoEta;
   vector<float>   *pfphoPhi;
   Int_t           nEle;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleESEn;
   vector<float>   *eleESEnP1;
   vector<float>   *eleESEnP2;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleSIP;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>    *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eledEtaAtCalo;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *eleESEffSigmaRR;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<float>   *elePFMiniIso;
   vector<float>   *eleIDMVA;
   vector<float>   *eleIDMVAHZZ;
   vector<float>   *eledEtaseedAtVtx;
   vector<float>   *eleE1x5;
   vector<float>   *eleE2x5;
   vector<float>   *eleE5x5;
   vector<float>   *eleE1x5Full5x5;
   vector<float>   *eleE2x5Full5x5;
   vector<float>   *eleE5x5Full5x5;
   vector<float>   *eleR9Full5x5;
   vector<int>     *eleEcalDrivenSeed;
   vector<float>   *eleDr03EcalRecHitSumEt;
   vector<float>   *eleDr03HcalDepth1TowerSumEt;
   vector<float>   *eleDr03HcalDepth2TowerSumEt;
   vector<float>   *eleDr03HcalTowerSumEt;
   vector<float>   *eleDr03TkSumPt;
   vector<float>   *elecaloEnergy;
   vector<float>   *eleTrkdxy;
   vector<float>   *eleKFHits;
   vector<float>   *eleKFChi2;
   vector<float>   *eleGSFChi2;
   vector<vector<float> > *eleGSFPt;
   vector<vector<float> > *eleGSFEta;
   vector<vector<float> > *eleGSFPhi;
   vector<vector<float> > *eleGSFCharge;
   vector<vector<int> > *eleGSFHits;
   vector<vector<int> > *eleGSFMissHits;
   vector<vector<int> > *eleGSFNHitsMax;
   vector<vector<float> > *eleGSFVtxProb;
   vector<vector<float> > *eleGSFlxyPV;
   vector<vector<float> > *eleGSFlxyBS;
   vector<vector<float> > *eleBCEn;
   vector<vector<float> > *eleBCEta;
   vector<vector<float> > *eleBCPhi;
   vector<vector<float> > *eleBCS25;
   vector<vector<float> > *eleBCS15;
   vector<vector<float> > *eleBCSieie;
   vector<vector<float> > *eleBCSieip;
   vector<vector<float> > *eleBCSipip;
   vector<unsigned int> *eleFiredSingleTrgs;
   vector<unsigned int> *eleFiredDoubleTrgs;
   vector<unsigned int> *eleFiredL1Trgs;
   vector<unsigned short> *eleIDbit;
   Int_t           npfHF;
   vector<float>   *pfHFEn;
   vector<float>   *pfHFECALEn;
   vector<float>   *pfHFHCALEn;
   vector<float>   *pfHFPt;
   vector<float>   *pfHFEta;
   vector<float>   *pfHFPhi;
   vector<float>   *pfHFIso;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<unsigned short> *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muMatches;
   vector<int>     *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
  // vector<float>   *muPFMiniIso;
   vector<unsigned int> *muFiredTrgs;
   vector<unsigned int> *muFiredL1Trgs;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   Int_t           nJet;
   vector<float>   *jetPt;
   vector<float>   *jetEn;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetMt;
   vector<float>   *jetArea;
   vector<float>   *jetMass;
   vector<float>   *jetPx;
   vector<float>   *jetPy;
   vector<float>   *jetPz;
   vector<float>   *jetLeadTrackPt;
   vector<float>   *jetLeadTrackEta;
   vector<float>   *jetLeadTrackPhi;
   vector<int>     *jetLepTrackPID;
   vector<float>   *jetLepTrackPt;
   vector<float>   *jetLepTrackEta;
   vector<float>   *jetLepTrackPhi;
   vector<float>   *jetCSV2BJetTags;
   vector<float>   *jetJetProbabilityBJetTags;
   vector<float>   *jetpfCombinedMVAV2BJetTags;
   vector<float>   *jetDeepCSVTags_b;
   vector<float>   *jetDeepCSVTags_bb;
   vector<float>   *jetDeepCSVTags_c;
   vector<float>   *jetDeepCSVTags_udsg;
   vector<bool>    *jetPFLooseId;
   vector<int>     *jetID;
   vector<float>   *jetPUID;
   vector<float>   *jetJECUnc;
//#   vector<float>   *jetJERSmearing;
//#   vector<float>   *jetJERSmearingUp;
//#   vector<float>   *jetJERSmearingDown;
   vector<unsigned int> *jetFiredTrgs;
   vector<float>   *jetCHF; //chargedHadronEnergyFraction
   vector<float>   *jetNHF; //neutralHadronEnergyFraction
   vector<float>   *jetCEF; //chargedEmEnergyFraction
   vector<float>   *jetNEF; //neutralEmEnergyFraction
   vector<int>     *jetNCH; //chargedMultiplicity
   vector<int>     *jetNNP; //NumberofNeutralParticles or neutal Multiplicity
   vector<float>   *jetMUF; //Muon energy fraction
   vector<float>   *jetVtxPt;
   vector<float>   *jetVtxMass;
   vector<float>   *jetVtxNtrks;
   vector<float>   *jetVtx3DVal;
   vector<float>   *jetVtx3DSig;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConstituents; //Number of constituents
 
/*#   Int_t           nAK8Jet;
#   vector<float>   *AK8JetPt;
#   vector<float>   *AK8JetEn;
#   vector<float>   *AK8JetRawPt;
 #   vector<float>   *AK8JetRawEn;
   # vector<float>   *AK8JetEta;
#    vector<float>   *AK8JetPhi;
#    vector<float>   *AK8JetMass;
#    vector<float>   *AK8Jet_tau1;
#    vector<float>   *AK8Jet_tau2;
#    vector<float>   *AK8Jet_tau3;
#    vector<float>   *AK8JetCHF;
#    vector<float>   *AK8JetNHF;
#    vector<float>   *AK8JetCEF;
#    vector<float>   *AK8JetNEF;
#    vector<int>     *AK8JetNCH;
 #   vector<int>     *AK8JetNNP;
 #   vector<float>   *AK8JetMUF;
 #   vector<int>     *AK8Jetnconstituents# ;
#    vector<bool>    *AK8JetPFLooseId;
#    vector<bool>    *AK8JetPFTightLepVetoId;
#    vector<float>   *AK8JetSoftDropMass;
#    vector<float>   *AK8JetSoftDropMassCorr# ;
 #   vector<float>   *AK8JetPrunedMass;
 #   vector<float>   *AK8JetPrunedMassCorr;
 #   vector<float>   *AK8JetpfBoostedDSVBTag;
 #   vector<float>   *AK8JetDSVnewV4;
 #   vector<float>   *AK8JetCSV;
 #   vector<float>   *AK8JetJECUnc;
 #   vector<float>   *AK8JetL2L3corr;
 #   vector<float>   *AK8puppiPt;
 #   vector<float>   *AK8puppiMass;
 #   vector<float>   *AK8puppiEta;
 #   vector<float>   *AK8puppiPhi;
#    vector<float>   *AK8puppiTau1;
#    vector<float>   *AK8puppiTau2;
#    vector<float>   *AK8puppiTau3;
#    vector<float>   *AK8puppiSDL2L3corr;
   # vector<float>   *AK8puppiSDMass;
 #   vector<float>   *AK8puppiSDMassL2L3Corr;
#    vector<int>     *nAK8SDSJ;
#    vector<vector<float> > *AK8SDSJPt;
#    vector<vector<float> > *AK8SDSJEta;
#    vector<vector<float> > *AK8SDSJPhi;
#    vector<vector<float> > *AK8SDSJMass;
#    vector<vector<float> > *AK8SDSJE;
#    vector<vector<int> > *AK8SDSJCharge;
#    vector<vector<int> > *AK8SDSJFlavour;
#    vector<vector<float> > *AK8SDSJCSV;
#    vector<int>     *nAK8puppiSDSJ;
#    vector<vector<float> > *AK8puppiSDSJPt;
  #  vector<vector<float> > *AK8puppiSDSJEta;
  #  vector<vector<float> > *AK8puppiSDSJPhi;
  #  vector<vector<float> > *AK8puppiSDSJMass;
 #   vector<vector<float> > *AK8puppiSDSJE;
 #   vector<vector<int> > *AK8puppiSDSJCharge;
 #   vector<vector<int> > *AK8puppiSDSJFlavour;
#    vector<vector<float> > *AK8puppiSDSJCSV;
*/
   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_L1ECALPrefire;
   TBranch        *b_nTrksPV;   //!
   TBranch        *b_isPVGood;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_phoPrescale;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_pfMETPhi_T1JESUp;   //!
   TBranch        *b_pfMETPhi_T1JESDo;   //!
   TBranch        *b_pfMETPhi_T1UESUp;   //!
   TBranch        *b_pfMETPhi_T1UESDo;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phomaxXtalenergyFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoE1x3Full5x5;   //!
   TBranch        *b_phoE1x5Full5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   TBranch        *b_phoE2x5MaxFull5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoSeedBCE;   //!
   TBranch        *b_phoSeedBCEta;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoCITKChIso;   //!
   TBranch        *b_phoCITKPhoIso;   //!
   TBranch        *b_phoCITKNeuIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
//   TBranch        *b_phoSeedTime;   //!
//   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoxtalBits;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_npfPho;   //!
   TBranch        *b_pfphoEt;   //!
   TBranch        *b_pfphoEta;   //!
   TBranch        *b_pfphoPhi;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_eleESEnP1;   //!
   TBranch        *b_eleESEnP2;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleSIP;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eledEtaAtCalo;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_elePFMiniIso;   //!
   TBranch        *b_eleIDMVA;   //!
   TBranch        *b_eleIDMVAHZZ;   //!
   TBranch        *b_eledEtaseedAtVtx;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE2x5;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE1x5Full5x5;   //!
   TBranch        *b_eleE2x5Full5x5;   //!
   TBranch        *b_eleE5x5Full5x5;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_eleDr03EcalRecHitSumEt;   //!
   TBranch        *b_eleDr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_eleDr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_eleDr03HcalTowerSumEt;   //!
   TBranch        *b_eleDr03TkSumPt;   //!
   TBranch        *b_elecaloEnergy;   //!
   TBranch        *b_eleTrkdxy;   //!
   TBranch        *b_eleKFHits;   //!
   TBranch        *b_eleKFChi2;   //!
   TBranch        *b_eleGSFChi2;   //!
   TBranch        *b_eleGSFPt;   //!
   TBranch        *b_eleGSFEta;   //!
   TBranch        *b_eleGSFPhi;   //!
   TBranch        *b_eleGSFCharge;   //!
   TBranch        *b_eleGSFHits;   //!
   TBranch        *b_eleGSFMissHits;   //!
   TBranch        *b_eleGSFNHitsMax;   //!
   TBranch        *b_eleGSFVtxProb;   //!
   TBranch        *b_eleGSFlxyPV;   //!
   TBranch        *b_eleGSFlxyBS;   //!
   TBranch        *b_eleBCEn;   //!
   TBranch        *b_eleBCEta;   //!
   TBranch        *b_eleBCPhi;   //!
   TBranch        *b_eleBCS25;   //!
   TBranch        *b_eleBCS15;   //!
   TBranch        *b_eleBCSieie;   //!
   TBranch        *b_eleBCSieip;   //!
   TBranch        *b_eleBCSipip;   //!
   TBranch        *b_eleFiredSingleTrgs;   //!
   TBranch        *b_eleFiredDoubleTrgs;   //!
   TBranch        *b_eleFiredL1Trgs;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_npfHF;   //!
   TBranch        *b_pfHFEn;   //!
   TBranch        *b_pfHFECALEn;   //!
   TBranch        *b_pfHFHCALEn;   //!
   TBranch        *b_pfHFPt;   //!
   TBranch        *b_pfHFEta;   //!
   TBranch        *b_pfHFPhi;   //!
   TBranch        *b_pfHFIso;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
  // TBranch        *b_muPFMiniIso;   //!
   TBranch        *b_muFiredTrgs;   //!
   TBranch        *b_muFiredL1Trgs;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetPx;   //!
   TBranch        *b_jetPy;   //!
   TBranch        *b_jetPz;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetLeadTrackEta;   //!
   TBranch        *b_jetLeadTrackPhi;   //!
   TBranch        *b_jetLepTrackPID;   //!
   TBranch        *b_jetLepTrackPt;   //!
   TBranch        *b_jetLepTrackEta;   //!
   TBranch        *b_jetLepTrackPhi;   //!
   TBranch        *b_jetCSV2BJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetpfCombinedMVAV2BJetTags;   //!
   TBranch        *b_jetDeepCSVTags_b;
   TBranch        *b_jetDeepCSVTags_bb;
   TBranch        *b_jetDeepCSVTags_c;
   TBranch        *b_jetDeepCSVTags_udsg;

   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_jetPUID;   //!
   TBranch        *b_jetJECUnc;   //!
//#   TBranch        *b_jetJERSmearing;   //!
//#   TBranch        *b_jetJERSmearingUp;   //!
//#   TBranch        *b_jetJERSmearingDown;   //!
   TBranch        *b_jetFiredTrgs;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetNNP;   //!
   TBranch        *b_jetMUF;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtxNtrks;   //!
   TBranch        *b_jetVtx3DVal;   //!
   TBranch        *b_jetVtx3DSig;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!

/*#    TBranch        *b_nAK8Jet;   //!
#    TBranch        *b_AK8JetPt;   //!
#    TBranch        *b_AK8JetEn;   //!
#    TBranch        *b_AK8JetRawPt;   //!
#    TBranch        *b_AK8JetRawEn;   //!
#    TBranch        *b_AK8JetEta;   //!
#    TBranch        *b_AK8JetPhi;   //!
#    TBranch        *b_AK8JetMass;   //!
#    TBranch        *b_AK8Jet_tau1;   //!
#    TBranch        *b_AK8Jet_tau2;   //!
#    TBranch        *b_AK8Jet_tau3;   //!
 #   TBranch        *b_AK8JetCHF;   //!
 #   TBranch        *b_AK8JetNHF;   //!
 #   TBranch        *b_AK8JetCEF;   //!
#    TBranch        *b_AK8JetNEF;   //!
#    TBranch        *b_AK8JetNCH;   //!
#    TBranch        *b_AK8JetNNP;   //!
#    TBranch        *b_AK8JetMUF;   //!
#    TBranch        *b_AK8Jetnconstituents;   //!
#    TBranch        *b_AK8JetPFLooseId;   //!
#    TBranch        *b_AK8JetPFTightLepVetoId;   //!
#    TBranch        *b_AK8JetSoftDropMass;   //!
#    TBranch        *b_AK8JetSoftDropMassCorr;   //!
#    TBranch        *b_AK8JetPrunedMass;   //!
#    TBranch        *b_AK8JetPrunedMassCorr;   //!
#    TBranch        *b_AK8JetpfBoostedDSVBTag;   //!
#    TBranch        *b_AK8JetDSVnewV4;   //!
#    TBranch        *b_AK8JetCSV;   //!
#    TBranch        *b_AK8JetJECUnc;   //!
#   TBranch        *b_AK8JetL2L3corr;   //!
#   TBranch        *b_AK8puppiPt;   //!
#   TBranch        *b_AK8puppiMass;   //!
#   TBranch        *b_AK8puppiEta;   //!
#   TBranch        *b_AK8puppiPhi;   //!
#   TBranch        *b_AK8puppiTau1;   //!
#   TBranch        *b_AK8puppiTau2;   //!
#   TBranch        *b_AK8puppiTau3;   //!
#   TBranch        *b_AK8puppiSDL2L3corr;   //!
 #  TBranch        *b_AK8puppiSDMass;   //!
 #  TBranch        *b_AK8puppiSDMassL2L3Corr;   //!
 #  TBranch        *b_nAK8SDSJ;   //!
 #  TBranch        *b_AK8SDSJPt;   //!
 #  TBranch        *b_AK8SDSJEta;   //!
 #  TBranch        *b_AK8SDSJPhi;   //!
 #  TBranch        *b_AK8SDSJMass;   //!
 #  TBranch        *b_AK8SDSJE;   //!
 #  TBranch        *b_AK8SDSJCharge;   //!
 #  TBranch        *b_AK8SDSJFlavour;   //!
#   TBranch        *b_AK8SDSJCSV;   //!
 #  TBranch        *b_nAK8puppiSDSJ;   //!
 #  TBranch        *b_AK8puppiSDSJPt;   //!
 #  TBranch        *b_AK8puppiSDSJEta;   //!
 #  TBranch        *b_AK8puppiSDSJPhi;   //!
 #  TBranch        *b_AK8puppiSDSJMass;   //!
 #  TBranch        *b_AK8puppiSDSJE;   //!
 #  TBranch        *b_AK8puppiSDSJCharge;   //!
 #  TBranch        *b_AK8puppiSDSJFlavour;   //!
 #  TBranch        *b_AK8puppiSDSJCSV;   //!
*/
   PostAnalyzer_Data(TTree *tree=0);
   virtual ~PostAnalyzer_Data();
//   virtual void     Efficiency(TGraphAsymmErrors * eff);
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   //User-Defined functions
//   virtual Bool_t   PassHLT(vector<ULong64_t> trigBits);
   virtual Bool_t   PassHLT_Ele(vector<ULong64_t> trigBits);
   virtual Bool_t   PassHLT_Pho(vector<ULong64_t> trigBits);
   virtual Bool_t   PassHLT_Prescale(ULong64_t trigBit);
//   virtual Bool_t   PassHLT();
//   virtual Bool_t   PassHLT_Photon90();
//   virtual Bool_t   PassHLT_Photon120();
   //virtual Bool_t   GoodPrimaryVtx(Int_t &GoodVertex);
//   virtual Bool_t   read_from_file(int R, Long64_t E, int L); 
   virtual Bool_t   CutBasedPhotonID(Int_t ipho, TString phoWP);
   virtual Bool_t  CutBasedElectronID(Int_t iele, TString eleWP);
   virtual Bool_t   MVABasedPhotonID(Int_t ipho, TString phoWP);
   virtual Double_t EAChargedHadrons(Double_t eta);
   virtual Double_t EANeutralHadrons(Double_t eta);
   virtual Double_t EAPhotons(Double_t eta);
   virtual Int_t    FirstGoodPhoton(TString phoWP); 
   virtual vector<Int_t> GoodPhotons(TString phoWP); 
   virtual vector<Int_t> GoodMuons();
   virtual vector<Int_t>  GoodElectrons(TString eleWP);
   virtual Bool_t   ResSpikes(Int_t);
   virtual Bool_t   JetId(Int_t iJet, TString jetWP);
   virtual Int_t    FirstGoodJet(TString jetWP);
   virtual vector<Int_t> GoodJets(TString jetWP);

   //For 80X (WP cuts need to change for 76X)
   virtual Bool_t   CSVv2bTag(Int_t ijet, TString WP);

   virtual Double_t GetdEta(Double_t eta1, Double_t eta2);
   virtual Double_t GetdPhi(Double_t phi1, Double_t phi2);
   virtual Double_t GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
   virtual Double_t GetCosThetaStar(Double_t eta1, Double_t eta2);
  // virtual Double_t GetInvtMass(Int_t ph, Int_t jet);
   virtual Double_t GetInvtMass(Int_t ph, TLorentzVector wide);
   virtual Double_t GetZMass(Int_t first, Int_t second);//ZGamma
   virtual Double_t GetZGMass(Int_t first, Int_t second, Int_t third); //ZGamma
  std::vector<float> leadPt;

   virtual void     BookHistograms();

};

#endif

#ifdef PostAnalyzer_Data_cxx
PostAnalyzer_Data::PostAnalyzer_Data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("ggNtuplizer/EventTree",tree);
       // f->GetObject("EventTree",tree);
#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
         TChain * chain = new TChain("ggNtuplizer/EventTree","");
         // TChain * chain = new TChain("EventTree","");
       //     chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/ReReco-BCDEFG_PromptReco-H/Run2016B-23Sep2016_ReReco/Data_Run2016B-23Sep2016_ReReco_4151.root/ggNtuplizer/EventTree");
     // chain->Add("/eos/uscms/store/user/lpcqstar/13TeV/Ntuples/80X/Data/reMiniAOD03Feb2017/singlePhotonRun2016Bv2/Data_newSinglePhotonRun2016Bv2_1.root/ggNtuplizer/EventTree");
   // chain->Add("Data_ntuple.root");
     
      //Uncomment this part in script
/*
      ///--------------------Use this part while submitting job in one go for a dataset--------------///
      TString main_path = "${sourceDir1}";

      TSystemDirectory sourceDir1("sysDir",main_path);
      TList* fileList = sourceDir1.GetListOfFiles();
      TIter next(fileList);
      TSystemFile* fileName;

      int fileNumber = 1;
      int maxFiles = -1;

      while ((fileName = (TSystemFile*)next()) && fileNumber > maxFiles){
        if(TString(fileName->GetName()) == "." || TString(fileName->GetName()) == ".."){continue;}

        TString FullPathInputFile = (main_path+fileName->GetName());

        //cout << FullPathInputFile << endl;

        chain->Add(FullPathInputFile+"/ggNtuplizer/EventTree");

        fileNumber++;

      }

      //cout << "Total files in this set = " << fileNumber - 1 << endl;
*/
      ///-----------Only change Tstring part with this part for job su10mission in parts for a dataset-----///
      ifstream Dataset;
     Dataset.open("${dataset1}", ifstream::in);
     //   Dataset.open("/uscms_data/d3/jbabbar/trial/CMSSW_9_4_6_patch1/src/PA_Main//${dataset1}", ifstream::in);
      char datafilename[500];

     // if(Dataset) cout << "${dataset1} Opened" << endl;

      int a = 0;
      for(Int_t i = 1; i <= ${maxf} && i <= ${Total_files}; i++){      
	Dataset >> datafilename;   ////dataset >> datafilename always starts from the 1st line if the file is just opened and will start from the 
                                   //// next to last line if already opened.
	string fname(datafilename);
	string main_path = "${sourceDir1}";
     	string myevt = "/ggNtuplizer/EventTree";
     //   string myevt = "/EventTree";        
        if(i >= ${sf} && i <= ${Total_files}){
          chain->Add((main_path+fname+myevt).c_str());
          //cout<<(main_path+fname+myevt).c_str()<<endl;
          a++;
        }
      }
     
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

PostAnalyzer_Data::~PostAnalyzer_Data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
//   file->cd();

/*   TGraphAsymmErrors *eff = new TGraphAsymmErrors();
   eff->Divide(h_TrigPhotonPt[1], h_TrigPhotonPt[0]);
  eff->SetName("eff");
   eff->Write();*/
//   file->Write();
//   file->Close();
 
}

/*void PostAnalyzer_Data::Efficiency(TGraphAsymmErrors * eff)
{
//eff->SetName(Efficiency);
eff->Divide(h_TrigPhotonPt[1], h_TrigPhotonPt[0]);
eff->SetName("eff");
eff->Write();
h_TrigEff->Divide(h_TrigPhotonPt[1],h_TrigPhotonPt[0]);
}*/

Int_t PostAnalyzer_Data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PostAnalyzer_Data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PostAnalyzer_Data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   phoPrescale = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   phoCalibE = 0;
   phoCalibEt = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEn = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoE1x3 = 0;
   phoE1x5 = 0;
   phoE2x2 = 0;
   phoE2x5Max = 0;
   phoE5x5 = 0;
   phoESEffSigmaRR = 0;
   phomaxXtalenergyFull5x5 = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE1x3Full5x5 = 0;
   phoE1x5Full5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE2x5MaxFull5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoSeedBCE = 0;
   phoSeedBCEta = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   phoCITKChIso = 0;
   phoCITKPhoIso = 0;
   phoCITKNeuIso = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredL1Trgs = 0;
//   phoSeedTime = 0;
 //  phoSeedEnergy = 0;
   phoxtalBits = 0;
   phoIDbit = 0;
   pfphoEt = 0;
   pfphoEta = 0;
   pfphoPhi = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleESEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaAtCalo = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   elePFMiniIso = 0;
   eleIDMVA = 0;
   eleIDMVAHZZ = 0;
   eledEtaseedAtVtx = 0;
   eleE1x5 = 0;
   eleE2x5 = 0;
   eleE5x5 = 0;
   eleE1x5Full5x5 = 0;
   eleE2x5Full5x5 = 0;
   eleE5x5Full5x5 = 0;
   eleR9Full5x5 = 0;
   eleEcalDrivenSeed = 0;
   eleDr03EcalRecHitSumEt = 0;
   eleDr03HcalDepth1TowerSumEt = 0;
   eleDr03HcalDepth2TowerSumEt = 0;
   eleDr03HcalTowerSumEt = 0;
   eleDr03TkSumPt = 0;
   elecaloEnergy = 0;
   eleTrkdxy = 0;
   eleKFHits = 0;
   eleKFChi2 = 0;
   eleGSFChi2 = 0;
   eleGSFPt = 0;
   eleGSFEta = 0;
   eleGSFPhi = 0;
   eleGSFCharge = 0;
   eleGSFHits = 0;
   eleGSFMissHits = 0;
   eleGSFNHitsMax = 0;
   eleGSFVtxProb = 0;
   eleGSFlxyPV = 0;
   eleGSFlxyBS = 0;
   eleBCEn = 0;
   eleBCEta = 0;
   eleBCPhi = 0;
   eleBCS25 = 0;
   eleBCS15 = 0;
   eleBCSieie = 0;
   eleBCSieip = 0;
   eleBCSipip = 0;
   eleFiredSingleTrgs = 0;
   eleFiredDoubleTrgs = 0;
   eleFiredL1Trgs = 0;
   eleIDbit = 0;
   pfHFEn = 0;
   pfHFECALEn = 0;
   pfHFHCALEn = 0;
   pfHFPt = 0;
   pfHFEta = 0;
   pfHFPhi = 0;
   pfHFIso = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   //muPFMiniIso = 0;
   muFiredTrgs = 0;
   muFiredL1Trgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   jetPt = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetMt = 0;
   jetArea = 0;
   jetMass = 0;
   jetPx = 0;
   jetPy = 0;
   jetPz = 0;
   jetLeadTrackPt = 0;
   jetLeadTrackEta = 0;
   jetLeadTrackPhi = 0;
   jetLepTrackPID = 0;
   jetLepTrackPt = 0;
   jetLepTrackEta = 0;
   jetLepTrackPhi = 0;
   jetCSV2BJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetpfCombinedMVAV2BJetTags = 0;
   jetDeepCSVTags_b = 0;
   jetDeepCSVTags_bb = 0;
   jetDeepCSVTags_c = 0;
   jetDeepCSVTags_udsg = 0;

   jetPFLooseId = 0;
   jetID = 0;
   jetPUID = 0;
   jetJECUnc = 0;
 // # jetJERSmearing = 0;
//  # jetJERSmearingUp = 0;
// #  jetJERSmearingDown = 0;
   jetFiredTrgs = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetNNP = 0;
   jetMUF = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtxNtrks = 0;
   jetVtx3DVal = 0;
   jetVtx3DSig = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConstituents = 0;
   
 
/*  # AK8JetPt = 0;
  # AK8JetEn = 0;
  # AK8JetRawPt = 0;
  # AK8JetRawEn = 0;
  # AK8JetEta = 0;
  # AK8JetPhi = 0;
  # AK8JetMass = 0;
  # AK8Jet_tau1 = 0;
  # AK8Jet_tau2 = 0;
  # AK8Jet_tau3 = 0;
  # AK8JetCHF = 0;
  # AK8JetNHF = 0;
  # AK8JetCEF = 0;
 #  AK8JetNEF = 0;
 #  AK8JetNCH = 0;
  # AK8JetNNP = 0;
  # AK8JetMUF = 0;
 #  AK8Jetnconstituents = 0;
 #  AK8JetPFLooseId = 0;
 #  AK8JetPFTightLepVetoId = 0;
 #  AK8JetSoftDropMass = 0;
 #  AK8JetSoftDropMassCorr = 0;
 #  AK8JetPrunedMass = 0;
 #  AK8JetPrunedMassCorr = 0;
#   AK8JetpfBoostedDSVBTag = 0;
 #  AK8JetDSVnewV4 = 0;
#   AK8JetCSV = 0;
#   AK8JetJECUnc = 0;
#   AK8JetL2L3corr = 0;
#   AK8puppiPt = 0;
#   AK8puppiMass = 0;
#   AK8puppiEta = 0;
#   AK8puppiPhi = 0;
#   AK8puppiTau1 = 0;
#   AK8puppiTau2 = 0;
#   AK8puppiTau3 = 0;
#   AK8puppiSDL2L3corr = 0;
#   AK8puppiSDMass = 0;
#   AK8puppiSDMassL2L3Corr = 0;
#   nAK8SDSJ = 0;
#   AK8SDSJPt = 0;
#   AK8SDSJEta = 0;
#   AK8SDSJPhi = 0;
#   AK8SDSJMass = 0;
#   AK8SDSJE = 0;
#   AK8SDSJCharge = 0;
#   AK8SDSJFlavour = 0;
#   AK8SDSJCSV = 0;
 #  nAK8puppiSDSJ = 0;
  # AK8puppiSDSJPt = 0;
#   AK8puppiSDSJEta = 0;
 #  AK8puppiSDSJPhi = 0;
  # AK8puppiSDSJMass = 0;
 #  AK8puppiSDSJE = 0;
#   AK8puppiSDSJCharge = 0;
#   AK8puppiSDSJFlavour = 0;
#   AK8puppiSDSJCSV = 0;  */
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("L1ECALPrefire", &L1ECALPrefire, &b_L1ECALPrefire);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("phoPrescale", &phoPrescale, &b_phoPrescale);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
   fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
   fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
   fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phomaxXtalenergyFull5x5", &phomaxXtalenergyFull5x5, &b_phomaxXtalenergyFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE1x3Full5x5", &phoE1x3Full5x5, &b_phoE1x3Full5x5);
   fChain->SetBranchAddress("phoE1x5Full5x5", &phoE1x5Full5x5, &b_phoE1x5Full5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE2x5MaxFull5x5", &phoE2x5MaxFull5x5, &b_phoE2x5MaxFull5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoSeedBCE", &phoSeedBCE, &b_phoSeedBCE);
   fChain->SetBranchAddress("phoSeedBCEta", &phoSeedBCEta, &b_phoSeedBCEta);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoCITKChIso", &phoCITKChIso, &b_phoCITKChIso);
   fChain->SetBranchAddress("phoCITKPhoIso", &phoCITKPhoIso, &b_phoCITKPhoIso);
   fChain->SetBranchAddress("phoCITKNeuIso", &phoCITKNeuIso, &b_phoCITKNeuIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
//   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
//   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoxtalBits", &phoxtalBits, &b_phoxtalBits);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("npfPho", &npfPho, &b_npfPho);
   fChain->SetBranchAddress("pfphoEt", &pfphoEt, &b_pfphoEt);
   fChain->SetBranchAddress("pfphoEta", &pfphoEta, &b_pfphoEta);
   fChain->SetBranchAddress("pfphoPhi", &pfphoPhi, &b_pfphoPhi);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   fChain->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("elePFMiniIso", &elePFMiniIso, &b_elePFMiniIso);
   fChain->SetBranchAddress("eleIDMVA", &eleIDMVA, &b_eleIDMVA);
   fChain->SetBranchAddress("eleIDMVAHZZ", &eleIDMVAHZZ, &b_eleIDMVAHZZ);
   fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   fChain->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleDr03EcalRecHitSumEt", &eleDr03EcalRecHitSumEt, &b_eleDr03EcalRecHitSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt, &b_eleDr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt, &b_eleDr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("eleDr03HcalTowerSumEt", &eleDr03HcalTowerSumEt, &b_eleDr03HcalTowerSumEt);
   fChain->SetBranchAddress("eleDr03TkSumPt", &eleDr03TkSumPt, &b_eleDr03TkSumPt);
   fChain->SetBranchAddress("elecaloEnergy", &elecaloEnergy, &b_elecaloEnergy);
   fChain->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   fChain->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   fChain->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   fChain->SetBranchAddress("eleGSFChi2", &eleGSFChi2, &b_eleGSFChi2);
   fChain->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   fChain->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   fChain->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   fChain->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   fChain->SetBranchAddress("eleGSFHits", &eleGSFHits, &b_eleGSFHits);
   fChain->SetBranchAddress("eleGSFMissHits", &eleGSFMissHits, &b_eleGSFMissHits);
   fChain->SetBranchAddress("eleGSFNHitsMax", &eleGSFNHitsMax, &b_eleGSFNHitsMax);
   fChain->SetBranchAddress("eleGSFVtxProb", &eleGSFVtxProb, &b_eleGSFVtxProb);
   fChain->SetBranchAddress("eleGSFlxyPV", &eleGSFlxyPV, &b_eleGSFlxyPV);
   fChain->SetBranchAddress("eleGSFlxyBS", &eleGSFlxyBS, &b_eleGSFlxyBS);
   fChain->SetBranchAddress("eleBCEn", &eleBCEn, &b_eleBCEn);
   fChain->SetBranchAddress("eleBCEta", &eleBCEta, &b_eleBCEta);
   fChain->SetBranchAddress("eleBCPhi", &eleBCPhi, &b_eleBCPhi);
   fChain->SetBranchAddress("eleBCS25", &eleBCS25, &b_eleBCS25);
   fChain->SetBranchAddress("eleBCS15", &eleBCS15, &b_eleBCS15);
   fChain->SetBranchAddress("eleBCSieie", &eleBCSieie, &b_eleBCSieie);
   fChain->SetBranchAddress("eleBCSieip", &eleBCSieip, &b_eleBCSieip);
   fChain->SetBranchAddress("eleBCSipip", &eleBCSipip, &b_eleBCSipip);
   fChain->SetBranchAddress("eleFiredSingleTrgs", &eleFiredSingleTrgs, &b_eleFiredSingleTrgs);
   fChain->SetBranchAddress("eleFiredDoubleTrgs", &eleFiredDoubleTrgs, &b_eleFiredDoubleTrgs);
   fChain->SetBranchAddress("eleFiredL1Trgs", &eleFiredL1Trgs, &b_eleFiredL1Trgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("npfHF", &npfHF, &b_npfHF);
   fChain->SetBranchAddress("pfHFEn", &pfHFEn, &b_pfHFEn);
   fChain->SetBranchAddress("pfHFECALEn", &pfHFECALEn, &b_pfHFECALEn);
   fChain->SetBranchAddress("pfHFHCALEn", &pfHFHCALEn, &b_pfHFHCALEn);
   fChain->SetBranchAddress("pfHFPt", &pfHFPt, &b_pfHFPt);
   fChain->SetBranchAddress("pfHFEta", &pfHFEta, &b_pfHFEta);
   fChain->SetBranchAddress("pfHFPhi", &pfHFPhi, &b_pfHFPhi);
   fChain->SetBranchAddress("pfHFIso", &pfHFIso, &b_pfHFIso);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
//   fChain->SetBranchAddress("muPFMiniIso", &muPFMiniIso, &b_muPFMiniIso);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetPx", &jetPx, &b_jetPx);
   fChain->SetBranchAddress("jetPy", &jetPy, &b_jetPy);
   fChain->SetBranchAddress("jetPz", &jetPz, &b_jetPz);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
   fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
   fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
   fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
   fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
   fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
   fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetpfCombinedMVAV2BJetTags", &jetpfCombinedMVAV2BJetTags, &b_jetpfCombinedMVAV2BJetTags);
   fChain->SetBranchAddress("jetDeepCSVTags_b", &jetDeepCSVTags_b, &b_jetDeepCSVTags_b);
   fChain->SetBranchAddress("jetDeepCSVTags_bb", &jetDeepCSVTags_bb, &b_jetDeepCSVTags_bb);
   fChain->SetBranchAddress("jetDeepCSVTags_c", &jetDeepCSVTags_c, &b_jetDeepCSVTags_c);
   fChain->SetBranchAddress("jetDeepCSVTags_udsg", &jetDeepCSVTags_udsg, &b_jetDeepCSVTags_udsg);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
//#   fChain->SetBranchAddress("jetJERSmearing", &jetJERSmearing, &b_jetJERSmearing);
//#   fChain->SetBranchAddress("jetJERSmearingUp", &jetJERSmearingUp, &b_jetJERSmearingUp);
//#   fChain->SetBranchAddress("jetJERSmearingDown", &jetJERSmearingDown, &b_jetJERSmearingDown);
   fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetNNP", &jetNNP, &b_jetNNP);
   fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtxNtrks", &jetVtxNtrks, &b_jetVtxNtrks);
   fChain->SetBranchAddress("jetVtx3DVal", &jetVtx3DVal, &b_jetVtx3DVal);
   fChain->SetBranchAddress("jetVtx3DSig", &jetVtx3DSig, &b_jetVtx3DSig);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);

/*  # fChain->SetBranchAddress("nAK8Jet", &nAK8Jet, &b_nAK8Jet);
  # fChain->SetBranchAddress("AK8JetPt", &AK8JetPt, &b_AK8JetPt);
  # fChain->SetBranchAddress("AK8JetEn", &AK8JetEn, &b_AK8JetEn);
  # fChain->SetBranchAddress("AK8JetRawPt", &AK8JetRawPt, &b_AK8JetRawPt);
  # fChain->SetBranchAddress("AK8JetRawEn", &AK8JetRawEn, &b_AK8JetRawEn);
  # fChain->SetBranchAddress("AK8JetEta", &AK8JetEta, &b_AK8JetEta);
  # fChain->SetBranchAddress("AK8JetPhi", &AK8JetPhi, &b_AK8JetPhi);
  # fChain->SetBranchAddress("AK8JetMass", &AK8JetMass, &b_AK8JetMass);
  # fChain->SetBranchAddress("AK8Jet_tau1", &AK8Jet_tau1, &b_AK8Jet_tau1);
  # fChain->SetBranchAddress("AK8Jet_tau2", &AK8Jet_tau2, &b_AK8Jet_tau2);
  # fChain->SetBranchAddress("AK8Jet_tau3", &AK8Jet_tau3, &b_AK8Jet_tau3);
  # fChain->SetBranchAddress("AK8JetCHF", &AK8JetCHF, &b_AK8JetCHF);
 #  fChain->SetBranchAddress("AK8JetNHF", &AK8JetNHF, &b_AK8JetNHF);
 #  fChain->SetBranchAddress("AK8JetCEF", &AK8JetCEF, &b_AK8JetCEF);
 #  fChain->SetBranchAddress("AK8JetNEF", &AK8JetNEF, &b_AK8JetNEF);
 #  fChain->SetBranchAddress("AK8JetNCH", &AK8JetNCH, &b_AK8JetNCH);
 #  fChain->SetBranchAddress("AK8JetNNP", &AK8JetNNP, &b_AK8JetNNP);
 #  fChain->SetBranchAddress("AK8JetMUF", &AK8JetMUF, &b_AK8JetMUF);
 #  fChain->SetBranchAddress("AK8Jetnconstituents", &AK8Jetnconstituents, &b_AK8Jetnconstituents);
 #  fChain->SetBranchAddress("AK8JetPFLooseId", &AK8JetPFLooseId, &b_AK8JetPFLooseId);
 #  fChain->SetBranchAddress("AK8JetPFTightLepVetoId", &AK8JetPFTightLepVetoId, &b_AK8JetPFTightLepVetoId);
 ##  fChain->SetBranchAddress("AK8JetSoftDropMass", &AK8JetSoftDropMass, &b_AK8JetSoftDropMass);
 #  fChain->SetBranchAddress("AK8JetSoftDropMassCorr", &AK8JetSoftDropMassCorr, &b_AK8JetSoftDropMassCorr);
 #  fChain->SetBranchAddress("AK8JetPrunedMass", &AK8JetPrunedMass, &b_AK8JetPrunedMass);
 #  fChain->SetBranchAddress("AK8JetPrunedMassCorr", &AK8JetPrunedMassCorr, &b_AK8JetPrunedMassCorr);
 #  fChain->SetBranchAddress("AK8JetpfBoostedDSVBTag", &AK8JetpfBoostedDSVBTag, &b_AK8JetpfBoostedDSVBTag);
 #  fChain->SetBranchAddress("AK8JetDSVnewV4", &AK8JetDSVnewV4, &b_AK8JetDSVnewV4);
 #  fChain->SetBranchAddress("AK8JetCSV", &AK8JetCSV, &b_AK8JetCSV);
 #  fChain->SetBranchAddress("AK8JetJECUnc", &AK8JetJECUnc, &b_AK8JetJECUnc);
 #  fChain->SetBranchAddress("AK8JetL2L3corr", &AK8JetL2L3corr, &b_AK8JetL2L3corr);
 #  fChain->SetBranchAddress("AK8puppiPt", &AK8puppiPt, &b_AK8puppiPt);
 #  fChain->SetBranchAddress("AK8puppiMass", &AK8puppiMass, &b_AK8puppiMass);
 #  fChain->SetBranchAddress("AK8puppiEta", &AK8puppiEta, &b_AK8puppiEta);
 #  fChain->SetBranchAddress("AK8puppiPhi", &AK8puppiPhi, &b_AK8puppiPhi);
 #  fChain->SetBranchAddress("AK8puppiTau1", &AK8puppiTau1, &b_AK8puppiTau1);
 #  fChain->SetBranchAddress("AK8puppiTau2", &AK8puppiTau2, &b_AK8puppiTau2);
 #  fChain->SetBranchAddress("AK8puppiTau3", &AK8puppiTau3, &b_AK8puppiTau3);
 #  fChain->SetBranchAddress("AK8puppiSDL2L3corr", &AK8puppiSDL2L3corr, &b_AK8puppiSDL2L3corr);
 #  fChain->SetBranchAddress("AK8puppiSDMass", &AK8puppiSDMass, &b_AK8puppiSDMass);
 #  fChain->SetBranchAddress("AK8puppiSDMassL2L3Corr", &AK8puppiSDMassL2L3Corr, &b_AK8puppiSDMassL2L3Corr);
 #  fChain->SetBranchAddress("nAK8SDSJ", &nAK8SDSJ, &b_nAK8SDSJ);
 #  fChain->SetBranchAddress("AK8SDSJPt", &AK8SDSJPt, &b_AK8SDSJPt#);
 
 #  fChain->SetBranchAddress("AK8SDSJEta", &AK8SDSJEta, &b_AK8SDSJEta);
#   fChain->SetBranchAddress("AK8SDSJPhi", &AK8SDSJPhi, &b_AK8SDSJPhi);
 #  fChain->SetBranchAddress("AK8SDSJMass", &AK8SDSJMass, &b_AK8SDSJMass);
  # fChain->SetBranchAddress("AK8SDSJE", &AK8SDSJE, &b_AK8SDSJE);
 #  fChain->SetBranchAddress("AK8SDSJCharge", &AK8SDSJCharge, &b_AK8SDSJCharge);
 #  fChain->SetBranchAddress("AK8SDSJFlavour", &AK8SDSJFlavour, &b_AK8SDSJFlavour);
 #  fChain->SetBranchAddress("AK8SDSJCSV", &AK8SDSJCSV, &b_AK8SDSJCSV);
 #  fChain->SetBranchAddress("nAK8puppiSDSJ", &nAK8puppiSDSJ, &b_nAK8puppiSDSJ);
 #  fChain->SetBranchAddress("AK8puppiSDSJPt", &AK8puppiSDSJPt, &b_AK8puppiSDSJPt);
 #  fChain->SetBranchAddress("AK8puppiSDSJEta", &AK8puppiSDSJEta, &b_AK8puppiSDSJEta);
 #  fChain->SetBranchAddress("AK8puppiSDSJPhi", &AK8puppiSDSJPhi, &b_AK8puppiSDSJPhi);
 #  fChain->SetBranchAddress("AK8puppiSDSJMass", &AK8puppiSDSJMass, &b_AK8puppiSDSJMass);
 #  fChain->SetBranchAddress("AK8puppiSDSJE", &AK8puppiSDSJE, &b_AK8puppiSDSJE);
 #  fChain->SetBranchAddress("AK8puppiSDSJCharge", &AK8puppiSDSJCharge, &b_AK8puppiSDSJCharge);
 #  fChain->SetBranchAddress("AK8puppiSDSJFlavour", &AK8puppiSDSJFlavour, &b_AK8puppiSDSJFlavour);
 #  fChain->SetBranchAddress("AK8puppiSDSJCSV", &AK8puppiSDSJCSV, &b_AK8puppiSDSJCSV);*/
    Notify();
}

Bool_t PostAnalyzer_Data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
 // user if needed. The return value is currently not used.

   return kTRUE;
}

void PostAnalyzer_Data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PostAnalyzer_Data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Bool_t PostAnalyzer_Data::PassHLT_Ele(vector<ULong64_t> trigBits){

  int count = 0;
  bool trigFired = false;

  for(Int_t i = 0; i < trigBits.size(); i++){
    if(HLTEleMuX>>(trigBits[i]) & 1){            //HLTPho is a 64 bit integer having only those bits as one whose corresponding triggers have been fired
      count++;     
           }
}                                              //by this particular event. see the correspondence between triggers and bits in ggAnalysis/ggNtuplizer/
                                             //plugins/ggNtuplizer_globalEvent.cc. trigBit in input represent the bit for trigger we want to check                                                                                 //So if a trigger has been fired by evt, then its bit in HLTPho is 1 and right shift operator (>>) will
  if(count > 0) trigFired = true;                 //shift this bit to the first place and its bitwise AND (&) with 1 will be 1, otherwise not.
  return trigFired;
}


Bool_t PostAnalyzer_Data::PassHLT_Pho(vector<ULong64_t> trigBits){

  int count = 0;
  bool trigFired = false;

  for(Int_t i = 0; i < trigBits.size(); i++){
    if(HLTPho>>(trigBits[i]) & 1){
      count++;
           }
}

  if(count > 0) trigFired = true;
  return trigFired;
}



Bool_t PostAnalyzer_Data::PassHLT_Prescale(ULong64_t trigBit){ //HLTPhoIsPrescaled is 1 for trig bit whose prescale value > 1 i.e trig is prescaled
                                                               //So PassHLT_Prescale() will retrun true for the trigger with prescale > 1.
  bool trigPre = false;                                        //As for HLT_Photon165_HE10, prescale = 1, so this fun will return 0 for that.
   if((HLTPhoIsPrescaled >> trigBit) & 1){
    trigPre = true;
  }

                          
                                         
  return trigPre;
}

/*
Bool_t PostAnalyzer_Data::GoodPrimaryVtx(Int_t &GoodVertex){

  Bool_t passVtx = false;
  GoodVertex = 0;

  for(Int_t i=0; i < nVtx; ++i){
    if( (fabs((*vtz)[i])) <= Cut_Vtx_z &&
        (*vndof)[i] >= Cut_Vtx_ndof    &&
        !((*isFake)[i])                &&
        (fabs((*vrho)[i])) <= Cut_Vtx_rho )
      GoodVertex++;
  }

  if(GoodVertex > 0) passVtx = true;

  return passVtx;

}
*/

//Cut Based Ph ID for 13 TeV 2016data(https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Recommended_Working_points_for_2)
//Latest ID for full 2016 Rereco data (Spring16 selection)
Bool_t PostAnalyzer_Data::CutBasedPhotonID(Int_t ipho, TString phoWP){

  Bool_t PhID = false;

if(phoWP == "loose"){ //loose
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
        ((*phoHoverE)[ipho] < 0.04596)                   &&
        ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0106)    &&
        ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 1.694)                    &&
        ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 24.032 + 0.01512*(*phoEt)[ipho] + 2.259e-05*(*phoEt)[ipho]*(*phoEt)[ipho])                                      &&
        ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 2.876+0.004017*(*phoEt)[ipho]);

    }
if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.0590 )                    &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0272)      &&
      ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 2.089)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 19.722 + 0.0117*(*phoEt)[ipho] + 2.3e-05*(*phoEt)[ipho]*(*phoEt)[ipho])                                         &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 4.162 + 0.0037*(*phoEt)[ipho]);

    }
  }

if(phoWP == "medium"){ //medium
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
        ((*phoHoverE)[ipho] < 0.02197)                   &&
        ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.01015)    &&
        ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 1.141)                    &&
        ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 1.189 + 0.01512*(*phoEt)[ipho] + 2.259e-05*(*phoEt)[ipho]*(*phoEt)[ipho])                                       &&
        ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 2.08 + 0.004017*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.0326)                     &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0272)      &&
      ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 1.051)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 2.718 + 0.0117*(*phoEt)[ipho] + 2.3e-05*(*phoEt)[ipho]*(*phoEt)[ipho])                                         &&
      ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 3.867 + 0.0037*(*phoEt)[ipho]);

    }
  }
 
if(phoWP == "tight"){ //tight
    if((fabs((*phoSCEta)[ipho])) <= 1.4442){ //Barrel

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
        ((*phoHoverE)[ipho] < 0.02148)                   &&
        ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.00996)    &&
        ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 0.65)                    &&
        ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 0.317 + 0.01512*(*phoEt)[ipho] + 2.259e-05*(*phoEt)[ipho]*(*phoEt)[ipho])                                          &&
        ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 2.044 + 0.004017*(*phoEt)[ipho]);

    }
    if((fabs((*phoSCEta)[ipho])) > 1.4442){ //Endcap

      PhID = ((*phoEleVeto)[ipho] == 1)                 &&
      ((*phoHoverE)[ipho] < 0.0321)                     &&
      ((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.0271)      &&
      ((TMath::Max(((*phoPFChIso)[ipho] - rho*EAChargedHadrons((*phoSCEta)[ipho])), 0.0)) < 0.517)                      &&
      ((TMath::Max(((*phoPFNeuIso)[ipho] - rho*EANeutralHadrons((*phoSCEta)[ipho])), 0.0)) < 2.716 + 0.0117*(*phoEt)[ipho] + 2.3e-05*(*phoEt)[ipho]*(*phoEt)[ipho])                                           &&
 
     ((TMath::Max(((*phoPFPhoIso)[ipho] - rho*EAPhotons((*phoSCEta)[ipho])), 0.0)) < 3.032 + 0.0037*(*phoEt)[ipho]);

    }
  }
  return PhID;
}

//Electron ID 
//ZGamma
Bool_t PostAnalyzer_Data::CutBasedElectronID(Int_t iele, TString eleWP){

  Bool_t EleID = false;

if(eleWP == "loose"){ //loose
    if((fabs((*eleSCEta)[iele])) <= 1.479){ //Barrel

        EleID = ((*eleConvVeto)[iele] == 1)                 &&
        (((*elePFChIso)[iele] + TMath::Max(0.0, (*elePFNeuIso)[iele] + (*elePFPhoIso)[iele] - rho*EANeutralHadrons((*eleSCEta)[iele] )))/(*elePt)[iele] < ( 0.112 + 0.506 /(*elePt)[iele])) &&
        ((*eleHoverE)[iele] < 0.05+ 1.16/(*eleSCEn)[iele] + 0.0324*rho/(*eleSCEn)[iele])                   &&
        (fabs((*eledEtaAtVtx)[iele]) < 0.00377  ) && 
        (fabs((*eledPhiAtVtx)[iele])  < 0.0884)  && 
        ((*eleSigmaIEtaIEtaFull5x5)[iele] < 0.0112)    &&
        ((*eleMissHits)[iele] <=1 ) &&
        ((*eleEoverPInv)[iele] < 0.193 ) &&
        ((*eleD0)[iele] < 0.05) &&
        ((*eleDz)[iele] < 0.10) ;
    }
if((fabs((*eleSCEta)[iele])) > 1.479){ //Endcap

      EleID = ((*eleConvVeto)[iele] == 1)                 &&
      (((*elePFChIso)[iele] + TMath::Max(0.0, (*elePFNeuIso)[iele] + (*elePFPhoIso)[iele] - rho*EANeutralHadrons((*eleSCEta)[iele] )))/(*elePt)[iele] < ( 0.108 + 0.963 /(*elePt)[iele])) && 
      ((*eleHoverE)[iele] < 0.0441 + 2.54/(*eleSCEn)[iele] + 0.0183*rho/(*eleSCEn)[iele])                    &&
      ((*eleSigmaIEtaIEtaFull5x5)[iele] < 0.0425)      &&
      (fabs((*eledEtaAtVtx)[iele]) < 0.00674  ) &&
      (fabs((*eledPhiAtVtx)[iele])  < 0.169 )  &&
      ((*eleMissHits)[iele] <=1 ) &&
      ((*eleEoverPInv)[iele] < 0.111) &&
      ((*eleD0)[iele] < 0.10) &&
      ((*eleDz)[iele] < 0.20) ;
     }
    }

if(eleWP == "medium"){
if((fabs((*eleSCEta)[iele])) <= 1.479){
EleID = ((*eleConvVeto)[iele] == 1)                 &&
        ( ((*elePFChIso)[iele] + TMath::Max(0.0, (*elePFNeuIso)[iele] + (*elePFPhoIso)[iele] - rho*EANeutralHadrons((*eleSCEta)[iele] )))/(*elePt)[iele] < ( 0.0478 + 0.506 /(*elePt)[iele]) ) &&
        ((*eleHoverE)[iele] < 0.046 + 1.16/(*eleSCEn)[iele] + 0.0324*rho/(*eleSCEn)[iele])                   &&
        (fabs((*eledEtaAtVtx)[iele]) < 0.0032  ) &&
        (fabs((*eledPhiAtVtx)[iele])  < 0.0547)  &&
        ((*eleSigmaIEtaIEtaFull5x5)[iele] < 0.0106)    &&
        ((*eleMissHits)[iele] <=1 ) &&
        ((*eleEoverPInv)[iele] < 0.184 ) &&
        ((*eleD0)[iele] < 0.05) &&
        ((*eleDz)[iele] < 0.10) ;
 }
if((fabs((*eleSCEta)[iele])) > 1.479){
 EleID = ((*eleConvVeto)[iele] == 1)                 &&
        (((*elePFChIso)[iele] + TMath::Max(0.0, (*elePFNeuIso)[iele] + (*elePFPhoIso)[iele] - rho*EANeutralHadrons((*eleSCEta)[iele] )))/(*elePt)[iele] < ( 0.0658 + 0.963 /(*elePt)[iele])) &&
      ((*eleHoverE)[iele] < 0.0275 + 2.52/(*eleSCEn)[iele] + 0.0183*rho/(*eleSCEn)[iele])                    &&
      ((*eleSigmaIEtaIEtaFull5x5)[iele] < 0.0387)      &&
      (fabs((*eledEtaAtVtx)[iele]) < 0.00632  ) &&
      (fabs((*eledPhiAtVtx)[iele])  < 0.0394 )  &&
      ((*eleMissHits)[iele] <=1 ) &&
      ((*eleEoverPInv)[iele] < 0.0721) &&
      ((*eleD0)[iele] < 0.10) &&
      ((*eleDz)[iele] < 0.20) ;
     }
    }


if(eleWP == "tight"){
if((fabs((*eleSCEta)[iele])) <= 1.479){
EleID = ((*eleConvVeto)[iele] == 1)                 &&
        ( ((*elePFChIso)[iele] + TMath::Max(0.0, (*elePFNeuIso)[iele] + (*elePFPhoIso)[iele] - rho*EANeutralHadrons((*eleSCEta)[iele] )))/(*elePt)[iele] < ( 0.0287 + 0.506 /(*elePt)[iele]) ) &&
        ((*eleHoverE)[iele] < 0.026 + 1.15/(*eleSCEn)[iele] + 0.0324*rho/(*eleSCEn)[iele])                   &&
        (fabs((*eledEtaAtVtx)[iele]) < 0.00255  ) &&
        (fabs((*eledPhiAtVtx)[iele])  < 0.022)  &&
        ((*eleSigmaIEtaIEtaFull5x5)[iele] < 0.0104)    &&
        ((*eleMissHits)[iele] <=1 ) &&
        ((*eleEoverPInv)[iele] < 0.159 ) &&
        ((*eleD0)[iele] < 0.05) &&
        ((*eleDz)[iele] < 0.10) ;
 }
if((fabs((*eleSCEta)[iele])) > 1.479){
 EleID = ((*eleConvVeto)[iele] == 1)                 &&
        (((*elePFChIso)[iele] + TMath::Max(0.0, (*elePFNeuIso)[iele] + (*elePFPhoIso)[iele] - rho*EANeutralHadrons((*eleSCEta)[iele] )))/(*elePt)[iele] < ( 0.0445 + 0.963 /(*elePt)[iele])) &&
      ((*eleHoverE)[iele] < 0.0188 + 2.06/(*eleSCEn)[iele] + 0.0183*rho/(*eleSCEn)[iele])                    &&
      ((*eleSigmaIEtaIEtaFull5x5)[iele] < 0.0353)      &&
      (fabs((*eledEtaAtVtx)[iele]) < 0.00501  ) &&
      (fabs((*eledPhiAtVtx)[iele])  < 0.0236 )  &&
      ((*eleMissHits)[iele] <=1 ) &&
      ((*eleEoverPInv)[iele] < 0.0197) &&
      ((*eleD0)[iele] < 0.10) &&
      ((*eleDz)[iele] < 0.20) ;
     }
    }
  return EleID;
  }
Bool_t PostAnalyzer_Data::MVABasedPhotonID(Int_t ipho, TString phoWP){

Bool_t PhID = false;
if(phoWP == "tight"){
PhID = (((fabs((*phoSCEta)[ipho]) < 1.444 && (*phoIDMVA)[ipho] > 0.42 ) ||  (fabs((*phoSCEta)[ipho]) > 1.5 && fabs((*phoSCEta)[ipho])) < 2.5  &&  ((*phoIDMVA)[ipho] > 0.14)) && (*phoEt)[ipho] > 15  && (*phoEleVeto)[ipho] == 1) ;
}
return PhID;
}


//(https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Selection implementation details for SPRING16)
Double_t PostAnalyzer_Data::EAChargedHadrons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.0112;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.0108;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0106;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.01002;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.0098;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.0089;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.0087;

  return EffArea;

}

Double_t PostAnalyzer_Data::EANeutralHadrons(Double_t eta){

  Double_t EffArea = 0;
/*
  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.0668;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1054;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0786;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0233;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.0078;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.0028;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.0137;
*/
//ZGamma included EA For CutBasedElectron ID
  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.1440;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1562;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.1032;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0859;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.1116;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.1321;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.1654;
 
 return EffArea;

}

Double_t PostAnalyzer_Data::EAPhotons(Double_t eta){

  Double_t EffArea = 0;

  if( fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffArea = 0.1113;
  if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.0953;
  if( fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffArea = 0.0619;
  if( fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffArea = 0.0837;
  if( fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffArea = 0.1070;
  if( fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffArea = 0.1212;
  if( fabs(eta) >= 2.4                        ) EffArea = 0.1466;

  return EffArea;

}

Int_t PostAnalyzer_Data::FirstGoodPhoton(TString phoWP){

  Int_t pc = -1;
  Bool_t ID = false;
  for(Int_t i = 0; i < nPho; i++){
    ID = MVABasedPhotonID(i, phoWP);
    if(ID){
      pc = i;
      break;
    }
  }
  return pc;
}

/*
  Bool_t PostAnalyzer_Data::read_from_file(int R, Long64_t E, int L){
        string RunN;
        string LumiN;
        string EvtN;

       int run,lum;
       Long64_t evt;
       bool duplicateEvt= false;

       fstream file("DoubleEGamma_v2_L1ECAL_21_2_24.txt");
       for(int i=1;i<4467;++i)
       {
            getline(file, RunN,':');
            //cout<<"RunN:"<<RunN<<endl;
            run = std::stoi(RunN);
            getline(file, EvtN,':');
            evt = std::stol(EvtN);
            //cout<<"EvtN:"<<EvtN<<endl;
            getline(file, LumiN);
            lum = std::stoi(LumiN);
            //cout<<"LumiN:"<<LumiN<<endl;
            if(R == run && L == lum && E == evt) duplicateEvt = true;
       }
       return  duplicateEvt;
   }
*/
/*
Bool_t PostAnalyzer_Data::read_from_file(int R, Long64_t E, int L){
        string RunN;
        string LumiN;
        string EvtN;

       int run,lum;
       Long64_t evt;
       bool duplicateEvt= false;

       fstream file("DoubleEGamma_v2_L1ECAL.txt");
       for(int i=1;i<117;++i)
       {
            getline(file, RunN,':');
            cout<<"RunN:"<<RunN<<endl;
            run = std::stoi(RunN);
            getline(file, EvtN,':');
            cout<<"EvtN:"<<EvtN<<endl;
            evt = std::stol(EvtN);
            getline(file, LumiN);
            cout<<"LumiN:"<<LumiN<<endl;
            lum = std::stoi(LumiN);
            if(R == run && L == lum && E == evt) duplicateEvt = true;
        cout<<"R:"<< R <<"L:"<< L << "E:"<< E <<endl;
         }
        cout<<"R:"<< R <<"L:"<< L << "E:"<< E <<endl;

       return  duplicateEvt;
}
*/
/*
vector<Int_t> PostAnalyzer_Data::GoodPhotons(TString phoWP){

  vector<Int_t> goodphs;
//  goodphs.clear();

  for(Int_t i = 0; i < nPho; i++){
    if(CutBasedPhotonID(i, phoWP) && ResSpikes(i) && (*phoEt)[i] > 30.0){
      goodphs.push_back(i);
    }
  }
  return goodphs;
}
*/

//ZGamma:MVA ID ~80% + phoEleVeto , Eta selected here so not included in PC selection
vector<Int_t> PostAnalyzer_Data::GoodPhotons(TString phoWP){

  vector<Int_t> goodphs;
  goodphs.clear();

  for(Int_t i = 0; i < nPho; i++){
  if(MVABasedPhotonID(i, phoWP)){
    goodphs.push_back(i);
   } 
 }
  return goodphs;
 }
 
                 
vector<Int_t> PostAnalyzer_Data::GoodMuons(){
vector<Int_t> goodmu;
goodmu.clear();
 for(Int_t i = 0; i < nMu; i++){
       if((*muIDbit)[i] >> 1 & 1){
            if((fabs((*muEta)[i])) < 2.5 && ((*muPt)[i]) > 26){
       goodmu.push_back(i);
       }      
      }
     }
   return goodmu;
} 

vector<Int_t> PostAnalyzer_Data::GoodElectrons(TString eleWP){

  vector<Int_t> goodele;
  goodele.clear();
 for(Int_t i = 0; i < nEle; i++){
    if(CutBasedElectronID(i, eleWP)){
      goodele.push_back(i);
    }
  }
  return goodele;
}


Bool_t PostAnalyzer_Data::ResSpikes(Int_t i){
  Bool_t spikes = false;
  if(// fabs((*phoSeedTime)[i]) < 3.0    &&  //time of arrival of ith photon at seed crystal
      (*phoSigmaIEtaIEtaFull5x5)[i] > 0.001   &&
      (*phoSigmaIPhiIPhiFull5x5)[i] > 0.001   &&
      //fabs(GetLICTD(i)) < 5.0               &&   //LICTD is the largest time difference between the seed crystal and the any other crystal 
      (*phoR9Full5x5)[i] < 1.0){
    spikes = true;
  }
  return spikes;
}

//Recommended JetID for 13 TeV 2016 data(https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016)
Bool_t PostAnalyzer_Data::JetId(Int_t iJet, TString jetWP){

  Bool_t JetID = false;

  if(fabs((*jetEta)[iJet]) <= 2.7){
    if(jetWP == "loose"){

      JetID = ((*jetNHF)[iJet] < 0.99 && (*jetNEF)[iJet] < 0.99 && (*jetNConstituents)[iJet] > 1 ) &&
	((fabs((*jetEta)[iJet]) <= 2.4 && (*jetCHF)[iJet] > 0 && (*jetNCH)[iJet] > 0 && (*jetCEF)[iJet] < 0.99) || fabs((*jetEta)[iJet]) > 2.4);
	 
    }
    if(jetWP == "tight"){

      JetID = ((*jetNHF)[iJet] < 0.90 && (*jetNEF)[iJet] < 0.90 && (*jetNConstituents)[iJet] > 1 ) &&
	((fabs((*jetEta)[iJet]) <= 2.4 && (*jetCHF)[iJet] > 0 && (*jetNCH)[iJet] > 0  && (*jetCEF)[iJet] < 0.99) || fabs((*jetEta)[iJet]) > 2.4);

    }  
    if(jetWP == "tightLepVeto"){

      JetID = ((*jetNHF)[iJet] < 0.90 && (*jetNEF)[iJet] < 0.90 && (*jetNConstituents)[iJet] > 1 && (*jetMUF)[iJet] < 0.8) &&
	((fabs((*jetEta)[iJet]) <= 2.4 && (*jetCHF)[iJet] > 0 && (*jetNCH)[iJet] > 0 && (*jetCEF)[iJet] < 0.80) || fabs((*jetEta)[iJet]) > 2.4);

    }
  }

  if(fabs((*jetEta)[iJet]) > 2.7 && fabs((*jetEta)[iJet]) <= 3.0){
    if(jetWP == "loose"){
      JetID = (*jetNEF)[iJet] > 0.01 && (*jetNHF)[iJet] < 0.98 && (*jetNNP)[iJet] > 2;
    }
    if(jetWP == "tight"){
      JetID = (*jetNEF)[iJet] > 0.01 && (*jetNHF)[iJet] < 0.98 && (*jetNNP)[iJet] > 2;
    }
  }
  
  if(fabs((*jetEta)[iJet]) > 3.0){
    if(jetWP == "loose"){
      JetID = (*jetNEF)[iJet] < 0.90 && (*jetNNP)[iJet] > 10;
    }
    if(jetWP == "tight"){
      JetID = (*jetNEF)[iJet] < 0.90 && (*jetNNP)[iJet]> 10;
    }
  }

  return JetID;
} 


Int_t PostAnalyzer_Data::FirstGoodJet(TString jetWP){

  Int_t jc = -1;
  for(Int_t i = 0; i < nJet; i++){
    Bool_t ID = false;
    ID = JetId(i, jetWP);
    if(ID && (*jetPt)[i] > 30.0){
      double minDR = 99.0;
      for(Int_t ph = 0; ph < nPho; ph++){
        if(CutBasedPhotonID(ph, "loose") && ResSpikes(ph) && (*phoEt)[ph] > 30.0){
          double temp_dR = GetdR((*phoSCEta)[ph], (*jetEta)[i], (*phoSCPhi)[ph], (*jetPhi)[i]);
          if(temp_dR < minDR) minDR = temp_dR;
        }
      }
      if(minDR > 0.5 && minDR < 99.0){
	jc = i;
	break;
      }
    }
  }
  return jc;
}

vector<Int_t> PostAnalyzer_Data::GoodJets(TString jetWP){

  vector<Int_t> goodjets;
  goodjets.clear();

  for(Int_t i = 0; i < nJet; i++){
    Bool_t ID = false;
    ID = JetId(i, jetWP);
    if(ID && (*jetPt)[i] > 30.0){
      double minDR = 99.0;
      for(Int_t ph = 0; ph < nPho; ph++){
        if(CutBasedPhotonID(ph, "loose") && ResSpikes(ph) && (*phoEt)[ph] > 30.0){
          double temp_dR = GetdR((*phoSCEta)[ph], (*jetEta)[i], (*phoSCPhi)[ph], (*jetPhi)[i]);
          if(temp_dR < minDR) minDR = temp_dR;
        }
      }
      if(minDR > 0.5 && minDR < 99.0) goodjets.push_back(i);
    }
  }
  return goodjets;
}


//80XReReco recommendations (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco)
Bool_t PostAnalyzer_Data::CSVv2bTag(Int_t ijet, TString WP){

  Bool_t passTag = false;
/*  0216 CSV 
  if(WP == "L"){ //loose
    if((*jetCSV2BJetTags)[ijet] > 0.5426) passTag = true;
  }
  if(WP == "M"){ //medium
    if((*jetCSV2BJetTags)[ijet] > 0.8484) passTag = true;
  }
  if(WP == "T"){ //tight
    if((*jetCSV2BJetTags)[ijet] > 0.9535) passTag = true;
  }
*/
if(WP == "L"){ //loose
    if(((*jetDeepCSVTags_b)[ijet] + (*jetDeepCSVTags_bb)[ijet]) > 0.2217) passTag = true;
  }
  if(WP == "M"){ //medium
    if(((*jetDeepCSVTags_b)[ijet] + (*jetDeepCSVTags_bb)[ijet]) > 0.6321) passTag = true;
  }
  if(WP == "T"){ //tight
    if(((*jetDeepCSVTags_b)[ijet] + (*jetDeepCSVTags_bb)[ijet]) > 0.8953) passTag = true;
  }

  return passTag;
}


Double_t PostAnalyzer_Data::GetdEta(Double_t eta1, Double_t eta2){

  Double_t dEta = fabs(eta1 - eta2);
  return dEta;
}

Double_t PostAnalyzer_Data::GetdPhi(Double_t phi1, Double_t phi2){

  Double_t dphi = (phi1 - phi2);
  Double_t twoPi = 2.0*(TMath::Pi());

  if(dphi < 0) dphi = - dphi;
  if(dphi >= (twoPi - dphi)) dphi = twoPi - dphi;

  return dphi;
}

  Double_t PostAnalyzer_Data::GetdR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){

  Double_t dEta = GetdEta(eta1, eta2);
  Double_t dPhi = GetdPhi(phi1, phi2);

  Double_t dR = 0.0;
  dR = sqrt(dEta*dEta + dPhi*dPhi);

  return dR;
}

//----------------------                                                          
// Compute cosThetaStar                                                                                                                            
//---------------------
Double_t PostAnalyzer_Data::GetCosThetaStar(Double_t eta1, Double_t eta2){ 
  Double_t theta = tanh( GetdEta(eta1,eta2)/2.0 );                                                                                                  
  return theta;                                                                                                                                     
}      

/*
Double_t PostAnalyzer_Data::GetInvtMass(Int_t pho, Int_t jet){

  Double_t mass = 0.0;

  TLorentzVector Pho;
  Pho.SetPtEtaPhiE((*phoEt)[pho], (*phoSCEta)[pho], (*phoSCPhi)[pho], (*phoE)[pho]);

  TLorentzVector Jet;
  Jet.SetPtEtaPhiE((*jetPt)[jet], (*jetEta)[jet], (*jetPhi)[jet], (*jetEn)[jet] );

  mass = (Pho+Jet).M();

  return mass;
}
*/
//ZGamma
Double_t PostAnalyzer_Data::GetZMass(Int_t  first, Int_t second){

  Double_t mass = 0.0;

  TLorentzVector FirstE;
  FirstE.SetPtEtaPhiE((*elePt)[first], (*eleEta)[first], (*elePhi)[first], (*eleEn)[first]);

  TLorentzVector SecondE;
  SecondE.SetPtEtaPhiE((*elePt)[second], (*eleEta)[second], (*elePhi)[second], (*eleEn)[second] );

  mass = fabs((FirstE+SecondE).M());
  //cout<<" mass value in functions "<< mass <<endl;
  return mass;
}
//ZGamma
Double_t PostAnalyzer_Data::GetZGMass(Int_t  first, Int_t second, Int_t third){

  Double_t mass = 0.0;

  TLorentzVector FirstE;
  FirstE.SetPtEtaPhiE((*elePt)[first], (*eleEta)[first], (*elePhi)[first], (*eleEn)[first]);

  TLorentzVector SecondE;
  SecondE.SetPtEtaPhiE((*elePt)[second], (*eleEta)[second], (*elePhi)[second], (*eleEn)[second] );

  TLorentzVector ThirdE;
  ThirdE.SetPtEtaPhiE((*phoEt)[third], (*phoSCEta)[third], (*phoSCPhi)[third], (*phoE)[third] );
  mass = fabs((FirstE+SecondE+ThirdE).M());
 // cout<<" mass value in functions "<< mass <<endl;
  return mass;
}


Double_t PostAnalyzer_Data::GetInvtMass(Int_t pho, TLorentzVector wide){

  Double_t mass = 0.0;
  TLorentzVector Pho;
  Pho.SetPtEtaPhiE((*phoEt)[pho], (*phoSCEta)[pho], (*phoSCPhi)[pho], (*phoE)[pho]);

  TLorentzVector Jet;
  Jet.SetPtEtaPhiE(wide.Pt(), wide.Eta(), wide.Phi(), wide.E());
    mass = (Pho+Jet).M();
      return mass;
      }

void PostAnalyzer_Data::BookHistograms(){
  //file->cd();
  
  char name[100];
  /*
  const Int_t nMassBins_qstar = 119;
  const Double_t MassBin_qstar[nMassBins_qstar+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};
  */
 /* const Int_t nMassBins = 109;
  const Double_t MassBin[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};

const Double_t MassBin[nMassBins+1] = {1, 15, 29, 44, 59, 74, 90, 107, 124, 142, 160, 179, 199, 219, 240, 262, 284, 307, 331, 356, 381, 407, 434, 462, 491, 521, 552, 584, 617, 651, 686, 723, 761, 800, 839, 880, 922, 965, 1010, 1056, 1104, 1154, 1205, 1258, 1313, 1370, 1428, 1488, 1550, 1615, 1682, 1751, 1822, 1896, 1972, 2051, 2132, 2216, 2303, 2393, 2486, 2582, 2681, 2784, 2890, 3000, 3113, 3230, 3351, 3476, 3605, 3739, 3877, 4020, 4168, 4321, 4479, 4642, 4811, 4985, 5165, 5351, 5543, 5742, 5948, 6161, 6381, 6608, 6843, 7086, 7337, 7596, 7864, 8141, 8427, 8723, 9029, 9345, 9672, 10010, 10359, 10720, 11093, 11479, 11878, 12290, 12716, 13156, 13611, 14082};
*/
const Int_t nMassBins = 108;
const Double_t MassBin[nMassBins+1] = {1, 15, 29, 44, 59, 74, 90, 107, 124, 142, 160, 179, 199, 219, 240, 262, 284, 307, 331, 356, 381, 407, 434, 462, 491, 521, 552, 584, 617, 651, 686, 723, 761,800, 840, 881, 924, 968, 1014, 1061, 1110, 1161, 1213, 1267, 1323, 1381, 1441, 1503, 1567, 1633, 1701, 1772, 1845, 1921, 1999, 2080, 2164, 2251, 2341, 2434, 2530, 2629, 2732, 2838, 2948, 3061, 3178, 3299, 3425, 3555, 3689, 3828, 3972, 4121, 4275, 4434, 4598, 4768, 4944, 5126, 5314, 5509, 5710, 5918, 6133, 6356, 6586, 6824, 7070, 7325, 7588, 7860, 8142, 8433, 8734, 9046, 9368, 9701, 10046, 10403, 10772, 11153, 11547, 11955, 12377, 12813, 13264, 13731, 14214 };

  std::string cut0[7] = {"noCut", "noMassCut", "MassCut", "1BTag_noMassCut", "1BTag_MassCut", "0BTag_noMassCut", "0BTag_MassCut"};
  for(Int_t hist0 = 0; hist0 < 7; ++hist0){

    sprintf(name, "h_PhotonPt_%s",cut0[hist0].c_str());
    h_PhotonPt[hist0] = new TH1F(name,"Pt distribution of photons",60,0.0,1200.0);
    h_PhotonPt[hist0]->GetYaxis()->SetTitle("Events/20 GeV");  h_PhotonPt[hist0]->GetYaxis()->CenterTitle();
    h_PhotonPt[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_PhotonPt[hist0]->GetXaxis()->CenterTitle();
    h_PhotonPt[hist0]->Sumw2();
  
    sprintf(name, "h_PhotonPtSublead_%s",cut0[hist0].c_str());
    h_PhotonPtSublead[hist0] = new TH1F(name,"Pt distribution of photons",60,0.0,1200.0);
    h_PhotonPtSublead[hist0]->GetYaxis()->SetTitle("Events/20 GeV");  h_PhotonPtSublead[hist0]->GetYaxis()->CenterTitle();
    h_PhotonPtSublead[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_PhotonPtSublead[hist0]->GetXaxis()->CenterTitle();
    h_PhotonPtSublead[hist0]->Sumw2();


     sprintf(name, "h_Photon_MVAID_%s",cut0[hist0].c_str());
      h_Photon_MVAID[hist0] = new TH1F(name,",MVA photons",6,0.0,3.0);
      h_Photon_MVAID[hist0]->GetYaxis()->SetTitle("Events");  h_Photon_MVAID[hist0]->GetYaxis()->CenterTitle();
      h_Photon_MVAID[hist0]->GetXaxis()->SetTitle(" (GeV)");  h_Photon_MVAID[hist0]->GetXaxis()->CenterTitle();
      h_Photon_MVAID[hist0]->Sumw2(); 
 
    sprintf(name, "h_PhotonCalibPt_%s",cut0[hist0].c_str());
    h_PhotonCalibPt[hist0] = new TH1F(name,"Pt distribution of Calibrated photons",150,0.0,6000.0);
    h_PhotonCalibPt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_PhotonCalibPt[hist0]->GetYaxis()->CenterTitle();
    h_PhotonCalibPt[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_PhotonCalibPt[hist0]->GetXaxis()->CenterTitle();
    h_PhotonCalibPt[hist0]->Sumw2();

    sprintf(name, "h_PhotonEta_%s",cut0[hist0].c_str());
    h_PhotonEta[hist0] = new TH1F(name,"Eta distribution of photons",25,-2.5,2.5);
    h_PhotonEta[hist0]->GetYaxis()->SetTitle("Events/0.2");  h_PhotonEta[hist0]->GetYaxis()->CenterTitle();
    h_PhotonEta[hist0]->GetXaxis()->SetTitle("#eta^{#gamma}");  h_PhotonEta[hist0]->GetXaxis()->CenterTitle();
    h_PhotonEta[hist0]->Sumw2();

    sprintf(name, "h_PhotonPhi_%s",cut0[hist0].c_str());
    h_PhotonPhi[hist0] = new TH1F(name,"Phi distribution of photons",100,-4.0,4.0);
    h_PhotonPhi[hist0]->GetYaxis()->SetTitle("Events");  h_PhotonPhi[hist0]->GetYaxis()->CenterTitle();
    h_PhotonPhi[hist0]->GetXaxis()->SetTitle("#phi^{#gamma}");   h_PhotonPhi[hist0]->GetXaxis()->CenterTitle();
    h_PhotonPhi[hist0]->Sumw2();

    sprintf(name, "h_Photon_SigmaIEtaIEta_%s",cut0[hist0].c_str());
    h_Photon_SigmaIEtaIEta[hist0] = new TH1F(name,"Photon SigmaIetaIeta Distribution",100,0.0,0.05);
    h_Photon_SigmaIEtaIEta[hist0]->GetYaxis()->SetTitle("Events");  h_Photon_SigmaIEtaIEta[hist0]->GetYaxis()->CenterTitle();
    h_Photon_SigmaIEtaIEta[hist0]->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");  h_Photon_SigmaIEtaIEta[hist0]->GetXaxis()->CenterTitle();
    h_Photon_SigmaIEtaIEta[hist0]->Sumw2();

    sprintf(name, "h_Photon_R9_%s",cut0[hist0].c_str());
    h_Photon_R9[hist0] = new TH1F(name,"Photon R9 Distribution",100,0.0,10.0);
    h_Photon_R9[hist0]->GetYaxis()->SetTitle("Events");       h_Photon_R9[hist0]->GetYaxis()->CenterTitle();
    h_Photon_R9[hist0]->GetXaxis()->SetTitle("Photon_r9");    h_Photon_R9[hist0]->GetXaxis()->CenterTitle();
    h_Photon_R9[hist0]->Sumw2();

    sprintf(name, "h_Photon_HoverE_%s",cut0[hist0].c_str());
    h_Photon_HoverE[hist0] = new TH1F(name,"Photon HoverE Distribution",50,0.0,0.1);
    h_Photon_HoverE[hist0]->GetYaxis()->SetTitle("Events");   h_Photon_HoverE[hist0]->GetYaxis()->CenterTitle();
    h_Photon_HoverE[hist0]->GetXaxis()->SetTitle("H/E");      h_Photon_HoverE[hist0]->GetXaxis()->CenterTitle();
    h_Photon_HoverE[hist0]->Sumw2();

    sprintf(name, "h_Photon_EleVeto_%s",cut0[hist0].c_str());
    h_Photon_EleVeto[hist0] = new TH1F(name,"Photon ElectronVeto",3,0,3);
    h_Photon_EleVeto[hist0]->GetYaxis()->SetTitle("Events");     h_Photon_EleVeto[hist0]->GetYaxis()->CenterTitle();
    h_Photon_EleVeto[hist0]->GetXaxis()->SetTitle("EleVeto");    h_Photon_EleVeto[hist0]->GetXaxis()->CenterTitle();
    h_Photon_EleVeto[hist0]->Sumw2();

    sprintf(name, "h_Photon_CorrPFChIso_%s",cut0[hist0].c_str());
    h_Photon_CorrPFChIso[hist0] = new TH1F(name,"Rho Corrected PF Charged Isolation",100,0,5);
    h_Photon_CorrPFChIso[hist0]->GetYaxis()->SetTitle("Events");        h_Photon_CorrPFChIso[hist0]->GetYaxis()->CenterTitle();
    h_Photon_CorrPFChIso[hist0]->GetXaxis()->SetTitle("CorrPFChIso");   h_Photon_CorrPFChIso[hist0]->GetXaxis()->CenterTitle();
    h_Photon_CorrPFChIso[hist0]->Sumw2();

    sprintf(name, "h_Photon_CorrPFNeuIso_%s",cut0[hist0].c_str());
    h_Photon_CorrPFNeuIso[hist0] = new TH1F(name,"Rho Corrected PF Neutral Isolation",60,0,30);
    h_Photon_CorrPFNeuIso[hist0]->GetYaxis()->SetTitle("Events");    h_Photon_CorrPFNeuIso[hist0]->GetYaxis()->CenterTitle();
    h_Photon_CorrPFNeuIso[hist0]->GetXaxis()->SetTitle("CorrPFNeuIso");    h_Photon_CorrPFNeuIso[hist0]->GetXaxis()->CenterTitle();
    h_Photon_CorrPFNeuIso[hist0]->Sumw2();

    sprintf(name, "h_Photon_CorrPFPhoIso_%s",cut0[hist0].c_str());
    h_Photon_CorrPFPhoIso[hist0] = new TH1F(name,"Rho Corrected PF Photon Isolation",20,0,10);
    h_Photon_CorrPFPhoIso[hist0]->GetYaxis()->SetTitle("Events");    h_Photon_CorrPFPhoIso[hist0]->GetYaxis()->CenterTitle();
    h_Photon_CorrPFPhoIso[hist0]->GetXaxis()->SetTitle("CorrPFPhoIso"); h_Photon_CorrPFPhoIso[hist0]->GetXaxis()->CenterTitle();
    h_Photon_CorrPFPhoIso[hist0]->Sumw2();
   

//    h_PhotonCounts = new TH2F("h_PhotonCounts", "Number of Photons per Event vs. Total Number of Photons;Total Number of Photons;Number of Photons per Event",
  //                                 5, 0, 5, 12, 0, 12);
      //h_PhotonCounts = new TH1F("h_PhotonCounts", "Number of Photons per Event;Number of Photons per Event;Events", 10, 0, 10);
 
     
      sprintf(name, "h_PhotonCounts_%s",cut0[hist0].c_str());
      h_PhotonCounts[hist0] = new TH1F(name,"Number of Photons per Event",5,0.5,5.5);
      h_PhotonCounts[hist0]->GetYaxis()->SetTitle("Events");  h_PhotonCounts[hist0]->GetYaxis()->CenterTitle();
      h_PhotonCounts[hist0]->GetXaxis()->SetTitle("Number of Photons per Event");  h_PhotonCounts[hist0]->GetXaxis()->CenterTitle();
      h_PhotonCounts[hist0]->Sumw2();

      sprintf(name, "h_GoodEtaElectrons_%s",cut0[hist0].c_str());
      h_GoodEtaElectrons[hist0] = new TH1F(name,"Number of GoodEta Electrons per Event",5,0.5,5.5);
      h_GoodEtaElectrons[hist0]->GetYaxis()->SetTitle("Events");  h_GoodEtaElectrons[hist0]->GetYaxis()->CenterTitle();
      h_GoodEtaElectrons[hist0]->GetXaxis()->SetTitle(" GoodEtaElectrons");  h_GoodEtaElectrons[hist0]->GetXaxis()->CenterTitle();
      h_GoodEtaElectrons[hist0]->Sumw2();


      sprintf(name, "h_GoodEtaPassElectrons_%s",cut0[hist0].c_str());
      h_GoodEtaPassElectrons[hist0] = new TH1F(name,"Number of GoodEtaPass Electrons per Event",5,0.5,5.5);
      h_GoodEtaPassElectrons[hist0]->GetYaxis()->SetTitle("Events");  h_GoodEtaPassElectrons[hist0]->GetYaxis()->CenterTitle();
      h_GoodEtaPassElectrons[hist0]->GetXaxis()->SetTitle("GoodEtaPassElectrons");  h_GoodEtaPassElectrons[hist0]->GetXaxis()->CenterTitle();
      h_GoodEtaPassElectrons[hist0]->Sumw2();

      sprintf(name, "h_GoodIsoPhotons_%s",cut0[hist0].c_str());
      h_GoodIsoPhotons[hist0] = new TH1F(name,"Number of GoodEtaPass Photons per Event",5,0.5,5.5);
      h_GoodIsoPhotons[hist0]->GetYaxis()->SetTitle("Events");  h_GoodIsoPhotons[hist0]->GetYaxis()->CenterTitle();
      h_GoodIsoPhotons[hist0]->GetXaxis()->SetTitle("GoodIsoPhotons");  h_GoodIsoPhotons[hist0]->GetXaxis()->CenterTitle();
      h_GoodIsoPhotons[hist0]->Sumw2();
     
      sprintf(name, "h_GoodEtaPhotons_%s",cut0[hist0].c_str());
      h_GoodEtaPhotons[hist0] = new TH1F(name,"Number of GoodEtaPass Photons per Event",5,0.5,5.5);
      h_GoodEtaPhotons[hist0]->GetYaxis()->SetTitle("Events");  h_GoodEtaPhotons[hist0]->GetYaxis()->CenterTitle();
      h_GoodEtaPhotons[hist0]->GetXaxis()->SetTitle("GoodEtaPhotons");  h_GoodEtaPhotons[hist0]->GetXaxis()->CenterTitle();
      h_GoodEtaPhotons[hist0]->Sumw2();

      sprintf(name, "h_GoodPtPhotons_%s",cut0[hist0].c_str());
      h_GoodPtPhotons[hist0] = new TH1F(name,"Number of GoodPt Photons per Event",5,0.5,5.5);
      h_GoodPtPhotons[hist0]->GetYaxis()->SetTitle("Events");  h_GoodPtPhotons[hist0]->GetYaxis()->CenterTitle();
      h_GoodPtPhotons[hist0]->GetXaxis()->SetTitle("GoodPtPhotons");  h_GoodPtPhotons[hist0]->GetXaxis()->CenterTitle();
      h_GoodPtPhotons[hist0]->Sumw2();  


 
      sprintf(name, "h_leptonPt_Lead_%s",cut0[hist0].c_str());
      h_leptonPt_Lead[hist0] = new TH1F(name,"Pt distribution of leading leptons",100,0.0,1000.0);
      h_leptonPt_Lead[hist0]->GetYaxis()->SetTitle("Events/10 GeV");  h_leptonPt_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonPt_Lead[hist0]->GetXaxis()->SetTitle("P_{T} (GeV)");  h_leptonPt_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonPt_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptonCalibPt_Lead_%s",cut0[hist0].c_str());
      h_leptonCalibPt_Lead[hist0] = new TH1F(name,"Pt distribution of Calibrated  leading leptons",100,0.0,1200.0);
      h_leptonCalibPt_Lead[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_leptonCalibPt_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonCalibPt_Lead[hist0]->GetXaxis()->SetTitle("P_{T} (GeV)");  h_leptonCalibPt_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonCalibPt_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptonEta_Lead_%s",cut0[hist0].c_str());
      h_leptonEta_Lead[hist0] = new TH1F(name,"Eta distribution of  leading leptons",25,-2.5,2.5);
      h_leptonEta_Lead[hist0]->GetYaxis()->SetTitle("Events/0.2");  h_leptonEta_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonEta_Lead[hist0]->GetXaxis()->SetTitle("#eta");  h_leptonEta_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonEta_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptonSCEta_Lead_%s",cut0[hist0].c_str());
      h_leptonSCEta_Lead[hist0] = new TH1F(name," SC Eta distribution of  leading leptons",25,-2.5,2.5);
      h_leptonSCEta_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonSCEta_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonSCEta_Lead[hist0]->GetXaxis()->SetTitle("SC#eta");  h_leptonSCEta_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonSCEta_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptonPhi_Lead_%s",cut0[hist0].c_str());
      h_leptonPhi_Lead[hist0] = new TH1F(name,"Phi distribution of  leading leptons",100,-4.0,4.0);
      h_leptonPhi_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonPhi_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonPhi_Lead[hist0]->GetXaxis()->SetTitle("#phi");   h_leptonPhi_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonPhi_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_SigmaIEtaIEta_Lead_%s",cut0[hist0].c_str());
      h_lepton_SigmaIEtaIEta_Lead[hist0] = new TH1F(name,"leading leptons SigmaIetaIeta Distribution",10,0.0,0.05);
      h_lepton_SigmaIEtaIEta_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_SigmaIEtaIEta_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_SigmaIEtaIEta_Lead[hist0]->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");  h_lepton_SigmaIEtaIEta_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_SigmaIEtaIEta_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_R9_Lead_%s",cut0[hist0].c_str());
      h_lepton_R9_Lead[hist0] = new TH1F(name,"leading leptons R9 Distribution",100,0.0,5.0);
      h_lepton_R9_Lead[hist0]->GetYaxis()->SetTitle("Events");       h_lepton_R9_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_R9_Lead[hist0]->GetXaxis()->SetTitle("lepton_r9");    h_lepton_R9_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_R9_Lead[hist0]->Sumw2();
   
      sprintf(name, "h_lepton_HoverE_Lead_%s",cut0[hist0].c_str());
      h_lepton_HoverE_Lead[hist0] = new TH1F(name,"leading leptons HoverE Distribution",50,0.0,0.4);
      h_lepton_HoverE_Lead[hist0]->GetYaxis()->SetTitle("Events/0.01");   h_lepton_HoverE_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_HoverE_Lead[hist0]->GetXaxis()->SetTitle("H/E");      h_lepton_HoverE_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_HoverE_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_charge_Lead_%s",cut0[hist0].c_str());
      h_lepton_charge_Lead[hist0] = new TH1F(name,"leading leptons charge",4,-2.0,2.0);
      h_lepton_charge_Lead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_charge_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_charge_Lead[hist0]->GetXaxis()->SetTitle("lepton charge");    h_lepton_charge_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_charge_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_PFMiniIso_Lead_%s",cut0[hist0].c_str());
      h_lepton_PFMiniIso_Lead[hist0] = new TH1F(name,"leading leptons MiniIso",100,0.0,200.0);
      h_lepton_PFMiniIso_Lead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_PFMiniIso_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_PFMiniIso_Lead[hist0]->GetXaxis()->SetTitle("lepton Mini Isolation");    h_lepton_PFMiniIso_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_PFMiniIso_Lead[hist0]->Sumw2();
      
      sprintf(name, "h_lepton_MissHits_Lead_%s",cut0[hist0].c_str());
      h_lepton_MissHits_Lead[hist0] = new TH1F(name,"leading leptons expected Missing Hits ",5,0.0,5.0);
      h_lepton_MissHits_Lead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_MissHits_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_MissHits_Lead[hist0]->GetXaxis()->SetTitle("leptons Missing Hits");    h_lepton_MissHits_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_MissHits_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_ConvVeto_Lead_%s",cut0[hist0].c_str());
      h_lepton_ConvVeto_Lead[hist0] = new TH1F(name,"leading leptons conversion Veto",2,0.0,2.0);
      h_lepton_ConvVeto_Lead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_ConvVeto_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_ConvVeto_Lead[hist0]->GetXaxis()->SetTitle("lepton Conversion Veto");    h_lepton_ConvVeto_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_ConvVeto_Lead[hist0]->Sumw2(); 

      sprintf(name, "h_lepton_MiniIsoByPt_Lead_%s",cut0[hist0].c_str());
      h_lepton_MiniIsoByPt_Lead[hist0] = new TH1F(name,"leading leptons MiniIso to Pt ratio",1000,0.0,100);
      h_lepton_MiniIsoByPt_Lead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_MiniIsoByPt_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_MiniIsoByPt_Lead[hist0]->GetXaxis()->SetTitle("lepton MiniIso/Pt ");    h_lepton_MiniIsoByPt_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_MiniIsoByPt_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptondEtaAtVtx_Lead_%s",cut0[hist0].c_str());
      h_leptondEtaAtVtx_Lead[hist0] = new TH1F(name,"dEta distribution of  leading leptons",20,-0.2,0.2);
      h_leptondEtaAtVtx_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_leptondEtaAtVtx_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptondEtaAtVtx_Lead[hist0]->GetXaxis()->SetTitle("#eta");  h_leptondEtaAtVtx_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptondEtaAtVtx_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptondPhiAtVtx_Lead_%s",cut0[hist0].c_str());
      h_leptondPhiAtVtx_Lead[hist0] = new TH1F(name,"dPhi distribution of subleading  leptons",20,-1.0,1.0);
      h_leptondPhiAtVtx_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_leptondPhiAtVtx_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptondPhiAtVtx_Lead[hist0]->GetXaxis()->SetTitle("#phi");   h_leptondPhiAtVtx_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptondPhiAtVtx_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptonD0_Lead_%s",cut0[hist0].c_str());
      h_leptonD0_Lead[hist0] = new TH1F(name,"D0 distribution of leading  leptons",20,-1.0,1.0);
      h_leptonD0_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonD0_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonD0_Lead[hist0]->GetXaxis()->SetTitle("D0");   h_leptonD0_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonD0_Lead[hist0]->Sumw2();

      sprintf(name, "h_leptonDz_Lead_%s",cut0[hist0].c_str());
      h_leptonDz_Lead[hist0] = new TH1F(name,"Dz distribution of leading  leptons",100,-20.0,20.0);
      h_leptonDz_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonDz_Lead[hist0]->GetYaxis()->CenterTitle();
      h_leptonDz_Lead[hist0]->GetXaxis()->SetTitle("Dz");   h_leptonDz_Lead[hist0]->GetXaxis()->CenterTitle();
      h_leptonDz_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_EoverPInv_Lead_%s",cut0[hist0].c_str());
      h_lepton_EoverPInv_Lead[hist0] = new TH1F(name,"EoverPInv distribution of leading  leptons",25,-1.0,1.0);
      h_lepton_EoverPInv_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_EoverPInv_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_EoverPInv_Lead[hist0]->GetXaxis()->SetTitle("EoverP");   h_lepton_EoverPInv_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_EoverPInv_Lead[hist0]->Sumw2();

      sprintf(name, "h_lepton_relIsowithEA_Lead_%s",cut0[hist0].c_str());
      h_lepton_relIsowithEA_Lead[hist0] = new TH1F(name,"relIso distribution of leading  leptons",40,0.0,8.0);
      h_lepton_relIsowithEA_Lead[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_relIsowithEA_Lead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_relIsowithEA_Lead[hist0]->GetXaxis()->SetTitle("relIso");   h_lepton_relIsowithEA_Lead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_relIsowithEA_Lead[hist0]->Sumw2();
     
      sprintf(name, "h_lepton_relIsowithEA_SubLead_%s",cut0[hist0].c_str());
      h_lepton_relIsowithEA_SubLead[hist0] = new TH1F(name,"relIso distribution of subleading  leptons",40,0.0,8.0);
      h_lepton_relIsowithEA_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_relIsowithEA_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_relIsowithEA_SubLead[hist0]->GetXaxis()->SetTitle("relIso");   h_lepton_relIsowithEA_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_relIsowithEA_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_EoverPInv_SubLead_%s",cut0[hist0].c_str());
      h_lepton_EoverPInv_SubLead[hist0] = new TH1F(name,"EoverPInv distribution of leading  leptons",25,-1.0,1.0);
      h_lepton_EoverPInv_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_EoverPInv_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_EoverPInv_SubLead[hist0]->GetXaxis()->SetTitle("EoverP");   h_lepton_EoverPInv_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_EoverPInv_SubLead[hist0]->Sumw2();

      sprintf(name, "h_leptonD0_SubLead_%s",cut0[hist0].c_str());
      h_leptonD0_SubLead[hist0] = new TH1F(name,"D0 distribution of subleading  leptons",40,-4.0,4.0);
      h_leptonD0_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonD0_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonD0_SubLead[hist0]->GetXaxis()->SetTitle("D0");   h_leptonD0_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonD0_SubLead[hist0]->Sumw2();

      sprintf(name, "h_leptonDz_SubLead_%s",cut0[hist0].c_str());
      h_leptonDz_SubLead[hist0] = new TH1F(name,"Dz distribution of subleading  leptons",100,-20.0,20.0);
      h_leptonDz_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonDz_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonDz_SubLead[hist0]->GetXaxis()->SetTitle("Dz");   h_leptonDz_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonDz_SubLead[hist0]->Sumw2();


      sprintf(name, "h_leptondPhiAtVtx_SubLead_%s",cut0[hist0].c_str());
      h_leptondPhiAtVtx_SubLead[hist0] = new TH1F(name,"dPhi distribution of subleading  leptons",20,-1.0,1.0);
      h_leptondPhiAtVtx_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptondPhiAtVtx_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptondPhiAtVtx_SubLead[hist0]->GetXaxis()->SetTitle("#phi");   h_leptondPhiAtVtx_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptondPhiAtVtx_SubLead[hist0]->Sumw2();


      sprintf(name, "h_leptondEtaAtVtx_SubLead_%s",cut0[hist0].c_str());
      h_leptondEtaAtVtx_SubLead[hist0] = new TH1F(name,"dEta distribution of  Subleading leptons",20,-0.2,0.2);
      h_leptondEtaAtVtx_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptondEtaAtVtx_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptondEtaAtVtx_SubLead[hist0]->GetXaxis()->SetTitle("#eta");  h_leptondEtaAtVtx_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptondEtaAtVtx_SubLead[hist0]->Sumw2();


      sprintf(name, "h_leptonPt_SubLead_%s",cut0[hist0].c_str());
      h_leptonPt_SubLead[hist0] = new TH1F(name,"Pt distribution of subleading leptons",30,0.0,300.0);
      h_leptonPt_SubLead[hist0]->GetYaxis()->SetTitle("Events GeV");  h_leptonPt_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonPt_SubLead[hist0]->GetXaxis()->SetTitle("P_{T} (GeV)");  h_leptonPt_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonPt_SubLead[hist0]->Sumw2();

      sprintf(name, "h_leptonCalibPt_SubLead_%s",cut0[hist0].c_str());
      h_leptonCalibPt_SubLead[hist0] = new TH1F(name,"Pt distribution of Calibrated subleading leptons",100,0.0,1200.0);
      h_leptonCalibPt_SubLead[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_leptonCalibPt_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonCalibPt_SubLead[hist0]->GetXaxis()->SetTitle("P_{T} (GeV)");  h_leptonCalibPt_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonCalibPt_SubLead[hist0]->Sumw2();

      sprintf(name, "h_leptonEta_SubLead_%s",cut0[hist0].c_str());
      h_leptonEta_SubLead[hist0] = new TH1F(name,"Eta distribution of subleading  leptons",25,-2.5,2.5);
      h_leptonEta_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonEta_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonEta_SubLead[hist0]->GetXaxis()->SetTitle("#eta");  h_leptonEta_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonEta_SubLead[hist0]->Sumw2();


      sprintf(name, "h_leptonSCEta_SubLead_%s",cut0[hist0].c_str());
      h_leptonSCEta_SubLead[hist0] = new TH1F(name," SC Eta distribution of subleading  leptons",25,-2.5,2.5);
      h_leptonSCEta_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonSCEta_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonSCEta_SubLead[hist0]->GetXaxis()->SetTitle("#eta");  h_leptonSCEta_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonSCEta_SubLead[hist0]->Sumw2();

      sprintf(name, "h_leptonPhi_SubLead_%s",cut0[hist0].c_str());
      h_leptonPhi_SubLead[hist0] = new TH1F(name,"Phi distribution of subleading  leptons",100,-4.0,4.0);
      h_leptonPhi_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_leptonPhi_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_leptonPhi_SubLead[hist0]->GetXaxis()->SetTitle("#phi");   h_leptonPhi_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_leptonPhi_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_SigmaIEtaIEta_SubLead_%s",cut0[hist0].c_str());
      h_lepton_SigmaIEtaIEta_SubLead[hist0] = new TH1F(name,"subleading leptons SigmaIetaIeta Distribution",10,0.0,0.05);
      h_lepton_SigmaIEtaIEta_SubLead[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_SigmaIEtaIEta_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_SigmaIEtaIEta_SubLead[hist0]->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");  h_lepton_SigmaIEtaIEta_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_SigmaIEtaIEta_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_R9_SubLead_%s",cut0[hist0].c_str());
      h_lepton_R9_SubLead[hist0] = new TH1F(name,"subleading leptons R9 Distribution",100,0.0,5.0);
      h_lepton_R9_SubLead[hist0]->GetYaxis()->SetTitle("Events");       h_lepton_R9_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_R9_SubLead[hist0]->GetXaxis()->SetTitle("lepton_r9");    h_lepton_R9_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_R9_SubLead[hist0]->Sumw2();

     sprintf(name, "h_lepton_HoverE_SubLead_%s",cut0[hist0].c_str());
      h_lepton_HoverE_SubLead[hist0] = new TH1F(name,"subleading leptons HoverE Distribution",50,0.0,0.4);
      h_lepton_HoverE_SubLead[hist0]->GetYaxis()->SetTitle("Events");   h_lepton_HoverE_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_HoverE_SubLead[hist0]->GetXaxis()->SetTitle("H/E");      h_lepton_HoverE_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_HoverE_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_charge_SubLead_%s",cut0[hist0].c_str());
      h_lepton_charge_SubLead[hist0] = new TH1F(name,"subleading leptons charge",4,-2.0,2.0);
      h_lepton_charge_SubLead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_charge_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_charge_SubLead[hist0]->GetXaxis()->SetTitle("lepton charge");    h_lepton_charge_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_charge_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_PFMiniIso_SubLead_%s",cut0[hist0].c_str());
      h_lepton_PFMiniIso_SubLead[hist0] = new TH1F(name,"subleading leptons MiniIso",100,0.0,200.0);
      h_lepton_PFMiniIso_SubLead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_PFMiniIso_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_PFMiniIso_SubLead[hist0]->GetXaxis()->SetTitle("lepton Mini Isolation");    h_lepton_PFMiniIso_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_PFMiniIso_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_MissHits_SubLead_%s",cut0[hist0].c_str());
      h_lepton_MissHits_SubLead[hist0] = new TH1F(name,"leptons expected Missing Hits ",5,0.0,5.0);
      h_lepton_MissHits_SubLead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_MissHits_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_MissHits_SubLead[hist0]->GetXaxis()->SetTitle("leptons Missing Hits");    h_lepton_MissHits_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_MissHits_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_ConvVeto_SubLead_%s",cut0[hist0].c_str());
      h_lepton_ConvVeto_SubLead[hist0] = new TH1F(name,"subleading leptons conversion Veto",2,0.0,2.0);
      h_lepton_ConvVeto_SubLead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_ConvVeto_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_ConvVeto_SubLead[hist0]->GetXaxis()->SetTitle("lepton Conversion Veto");    h_lepton_ConvVeto_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_ConvVeto_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_MiniIsoByPt_SubLead_%s",cut0[hist0].c_str());
      h_lepton_MiniIsoByPt_SubLead[hist0] = new TH1F(name,"subleading leptons MiniIso to Pt ratio",1000,0.0,100);
      h_lepton_MiniIsoByPt_SubLead[hist0]->GetYaxis()->SetTitle("Events");     h_lepton_MiniIsoByPt_SubLead[hist0]->GetYaxis()->CenterTitle();
      h_lepton_MiniIsoByPt_SubLead[hist0]->GetXaxis()->SetTitle("lepton MiniIso/Pt ");    h_lepton_MiniIsoByPt_SubLead[hist0]->GetXaxis()->CenterTitle();
      h_lepton_MiniIsoByPt_SubLead[hist0]->Sumw2();

      sprintf(name, "h_lepton_dPhi_%s",cut0[hist0].c_str());
      h_lepton_dPhi[hist0] = new TH1F(name,"Phi difference of  leptons",100,-4.0,4.0);
      h_lepton_dPhi[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_dPhi[hist0]->GetYaxis()->CenterTitle();
      h_lepton_dPhi[hist0]->GetXaxis()->SetTitle("#phi");   h_lepton_dPhi[hist0]->GetXaxis()->CenterTitle();
      h_lepton_dPhi[hist0]->Sumw2();

      sprintf(name, "h_lepton_dR_%s",cut0[hist0].c_str());
      h_lepton_dR[hist0] = new TH1F(name," Delta R of  leptons",100,0.0,10);
      h_lepton_dR[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_dR[hist0]->GetYaxis()->CenterTitle();
      h_lepton_dR[hist0]->GetXaxis()->SetTitle("Delta R");   h_lepton_dR[hist0]->GetXaxis()->CenterTitle();
      h_lepton_dR[hist0]->Sumw2();     

      sprintf(name, "h_lepton_dR1_%s",cut0[hist0].c_str());
      h_lepton_dR1[hist0] = new TH1F(name," Delta R between leading lepton and Photon",100,0.0,10);
      h_lepton_dR1[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_dR1[hist0]->GetYaxis()->CenterTitle();
      h_lepton_dR1[hist0]->GetXaxis()->SetTitle("Delta R");   h_lepton_dR1[hist0]->GetXaxis()->CenterTitle();
      h_lepton_dR1[hist0]->Sumw2();


      sprintf(name, "h_lepton_dR2_%s",cut0[hist0].c_str());
      h_lepton_dR2[hist0] = new TH1F(name," Delta R subleading leptons",100,0.0,10);
      h_lepton_dR2[hist0]->GetYaxis()->SetTitle("Events");  h_lepton_dR2[hist0]->GetYaxis()->CenterTitle();
      h_lepton_dR2[hist0]->GetXaxis()->SetTitle("Delta R");   h_lepton_dR2[hist0]->GetXaxis()->CenterTitle();
      h_lepton_dR2[hist0]->Sumw2();




      sprintf(name, "h_Zmass_%s",cut0[hist0].c_str());
      h_Zmass[hist0] = new TH1F(name," Z boson ",60,60,120);
      h_Zmass[hist0]->GetYaxis()->SetTitle("Events");     h_Zmass[hist0]->GetYaxis()->CenterTitle();
      h_Zmass[hist0]->GetXaxis()->SetTitle("Z (GeV) ");    h_Zmass[hist0]->GetXaxis()->CenterTitle();
      h_Zmass[hist0]->Sumw2();
 
      sprintf(name, "h_ZGmass_%s",cut0[hist0].c_str());
      h_ZGmass[hist0] = new TH1F(name," ZG ",280,0.0,7000);
      h_ZGmass[hist0]->GetYaxis()->SetTitle("Events");     h_ZGmass[hist0]->GetYaxis()->CenterTitle();
      h_ZGmass[hist0]->GetXaxis()->SetTitle("ZG (GeV) ");    h_ZGmass[hist0]->GetXaxis()->CenterTitle();
      h_ZGmass[hist0]->Sumw2();

      h_goodPV_noWt = new TH1F("h_goodPV_noWt","No. of Good Primary Vertices", 70, 0, 70);
      h_goodPV_noWt->GetYaxis()->SetTitle("Events");   h_goodPV_noWt->GetYaxis()->CenterTitle();
      h_goodPV_noWt->GetXaxis()->SetTitle("nPV");      h_goodPV_noWt->GetXaxis()->CenterTitle();
      h_goodPV_noWt->Sumw2();
/*  
      h_goodPV_WithCut = new TH1F("h_goodPV_WithCut","No. of Good Primary Vertices after cut ", 70, 0, 70);
      h_goodPV_WithCut->GetYaxis()->SetTitle("Events");   h_goodPV_WithCut->GetYaxis()->CenterTitle();
      h_goodPV_WithCut->GetXaxis()->SetTitle("nPV");      h_goodPV_WithCut->GetXaxis()->CenterTitle();
      h_goodPV_WithCut->Sumw2();
*/
    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_JetPt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJetPt_%s",cut0[hist0].c_str());
    h_bJetPt[hist0] = new TH1F(name,"Pt distribution of jets",150,0.0,6000.0);
    h_bJetPt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");  h_bJetPt[hist0]->GetYaxis()->CenterTitle();
    h_bJetPt[hist0]->GetXaxis()->SetTitle("P_{T}^{Jet} (GeV)");  h_bJetPt[hist0]->GetXaxis()->CenterTitle();
    h_bJetPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_JetEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJetEta_%s",cut0[hist0].c_str());
    h_bJetEta[hist0] = new TH1F(name,"Eta distribution of jets",120,-3.0,3.0);
    h_bJetEta[hist0]->GetYaxis()->SetTitle("Events");      h_bJetEta[hist0]->GetYaxis()->CenterTitle();
    h_bJetEta[hist0]->GetXaxis()->SetTitle("#eta^{Jet}");  h_bJetEta[hist0]->GetXaxis()->CenterTitle();
    h_bJetEta[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_JetPhi_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJetPhi_%s",cut0[hist0].c_str());
    h_bJetPhi[hist0] = new TH1F(name,"Phi distribution of jets",100,-4.0,4.0);
    h_bJetPhi[hist0]->GetYaxis()->SetTitle("Events");      h_bJetPhi[hist0]->GetYaxis()->CenterTitle();
    h_bJetPhi[hist0]->GetXaxis()->SetTitle("#phi^{Jet}");  h_bJetPhi[hist0]->GetXaxis()->CenterTitle();
    h_bJetPhi[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_Mt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_Mt_%s",cut0[hist0].c_str());
    h_bJet_Mt[hist0] = new TH1F(name,"Transverse mass distribution of jets",150,0.0,6000.0);
    h_bJet_Mt[hist0]->GetYaxis()->SetTitle("Events/40 GeV");      h_bJet_Mt[hist0]->GetYaxis()->CenterTitle();
    h_bJet_Mt[hist0]->GetXaxis()->SetTitle("M_{T}^{Jet}");  h_bJet_Mt[hist0]->GetXaxis()->CenterTitle();
    h_bJet_Mt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_area_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_area_%s",cut0[hist0].c_str());
    h_bJet_area[hist0] = new TH1F(name,"Jets area",25,0.0,1000.0);
    h_bJet_area[hist0]->GetYaxis()->SetTitle("Events/40 GeV");      h_bJet_area[hist0]->GetYaxis()->CenterTitle();
    h_bJet_area[hist0]->GetXaxis()->SetTitle("JetArea");  h_bJet_area[hist0]->GetXaxis()->CenterTitle();
    h_bJet_area[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_Mass_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_Mass_%s",cut0[hist0].c_str());
    h_bJet_Mass[hist0] = new TH1F(name,"Mass distribution of jets",150,0.0,6000.0);
    h_bJet_Mass[hist0]->GetYaxis()->SetTitle("Events/40 GeV");      h_bJet_Mass[hist0]->GetYaxis()->CenterTitle();
    h_bJet_Mass[hist0]->GetXaxis()->SetTitle("M^{Jet}");  h_bJet_Mass[hist0]->GetXaxis()->CenterTitle();
    h_bJet_Mass[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NHEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NHEF_%s",cut0[hist0].c_str());
    h_bJet_NHEF[hist0] = new TH1F(name,"Neutral Hadron Energy Fraction",25,0,1);
    h_bJet_NHEF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_NHEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NHEF[hist0]->GetXaxis()->SetTitle("Jet_NHEF");    h_bJet_NHEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NHEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NEEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NEEF_%s",cut0[hist0].c_str());
    h_bJet_NEEF[hist0] = new TH1F(name,"Neutral Em Energy Fraction",25,0,1);
    h_bJet_NEEF[hist0]->GetYaxis()->SetTitle("Events");       h_bJet_NEEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NEEF[hist0]->GetXaxis()->SetTitle("Jet_NEEF");     h_bJet_NEEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NEEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NConst_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NConst_%s",cut0[hist0].c_str());
    h_bJet_NConst[hist0] = new TH1F(name,"No. of Constituents",100,0,100);
    h_bJet_NConst[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_NConst[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NConst[hist0]->GetXaxis()->SetTitle("Jet_NConst");  h_bJet_NConst[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NConst[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_CHEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_CHEF_%s",cut0[hist0].c_str());
    h_bJet_CHEF[hist0] = new TH1F(name,"Charged Hadron Energy Fraction",25,0,1);
    h_bJet_CHEF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_CHEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_CHEF[hist0]->GetXaxis()->SetTitle("Jet_CHEF");    h_bJet_CHEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_CHEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_ChMult_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_ChMult_%s",cut0[hist0].c_str());
    h_bJet_ChMult[hist0] = new TH1F(name,"Charged Multiplicity",100,0,100);
    h_bJet_ChMult[hist0]->GetYaxis()->SetTitle("Events");       h_bJet_ChMult[hist0]->GetYaxis()->CenterTitle();
    h_bJet_ChMult[hist0]->GetXaxis()->SetTitle("Jet_ChMult");   h_bJet_ChMult[hist0]->GetXaxis()->CenterTitle();
    h_bJet_ChMult[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_CEEF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_CEEF_%s",cut0[hist0].c_str());
    h_bJet_CEEF[hist0] = new TH1F(name,"Charged Em Energy Fraction",25,0,1);
    h_bJet_CEEF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_CEEF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_CEEF[hist0]->GetXaxis()->SetTitle("Jet_CEEF");    h_bJet_CEEF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_CEEF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_MUF_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_MUF_%s",cut0[hist0].c_str());
    h_bJet_MUF[hist0] = new TH1F(name,"Muon Fraction",25,0,1);
    h_bJet_MUF[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_MUF[hist0]->GetYaxis()->CenterTitle();
    h_bJet_MUF[hist0]->GetXaxis()->SetTitle("Jet_MUF");    h_bJet_MUF[hist0]->GetXaxis()->CenterTitle();
    h_bJet_MUF[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_Jet_NNP_%s",cut0[hist0].c_str());
    else sprintf(name, "h_bJet_NNP_%s",cut0[hist0].c_str());
    h_bJet_NNP[hist0] = new TH1F(name,"Number of neutral particles",100,0,100);
    h_bJet_NNP[hist0]->GetYaxis()->SetTitle("Events");      h_bJet_NNP[hist0]->GetYaxis()->CenterTitle();
    h_bJet_NNP[hist0]->GetXaxis()->SetTitle("Jet_NNP");    h_bJet_NNP[hist0]->GetXaxis()->CenterTitle();
    h_bJet_NNP[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJetInvtMass_VarBin_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJetInvtMass_VarBin_%s",cut0[hist0].c_str());
    h_GbJetInvtMass_VarBin[hist0] = new TH1F(name, "Invt mass distribution", nMassBins, MassBin);
    h_GbJetInvtMass_VarBin[hist0]->GetYaxis()->SetTitle("Events/VarBin"); h_GbJetInvtMass_VarBin[hist0]->GetYaxis()->CenterTitle();
    h_GbJetInvtMass_VarBin[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)");    h_GbJetInvtMass_VarBin[hist0]->GetXaxis()->CenterTitle();
    h_GbJetInvtMass_VarBin[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJetInvtMass_UnitBin_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJetInvtMass_UnitBin_%s",cut0[hist0].c_str());
    h_GbJetInvtMass_UnitBin[hist0] = new TH1F(name,"Invt mass distribution", 14000, 0.0, 14000.0);
    h_GbJetInvtMass_UnitBin[hist0]->GetYaxis()->SetTitle("Events/UnitBin");   h_GbJetInvtMass_UnitBin[hist0]->GetYaxis()->CenterTitle();
    h_GbJetInvtMass_UnitBin[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)"); h_GbJetInvtMass_UnitBin[hist0]->GetXaxis()->CenterTitle();
    h_GbJetInvtMass_UnitBin[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJet_dEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJet_dEta_%s",cut0[hist0].c_str());
    h_GbJet_dEta[hist0] = new TH1F(name, "delta Eta dist", 120, 0, 6);
    h_GbJet_dEta[hist0]->GetYaxis()->SetTitle("Events");         h_GbJet_dEta[hist0]->GetYaxis()->CenterTitle();
    h_GbJet_dEta[hist0]->GetXaxis()->SetTitle("#Delta #eta");    h_GbJet_dEta[hist0]->GetXaxis()->CenterTitle();
    h_GbJet_dEta[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJet_dPhi_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJet_dPhi_%s",cut0[hist0].c_str());
    h_GbJet_dPhi[hist0] = new TH1F(name, "delta Phi dist", 64, 0, 3.2);
    h_GbJet_dPhi[hist0]->GetYaxis()->SetTitle("Events");       h_GbJet_dPhi[hist0]->GetYaxis()->CenterTitle();
    h_GbJet_dPhi[hist0]->GetXaxis()->SetTitle("#Delta #phi");  h_GbJet_dPhi[hist0]->GetXaxis()->CenterTitle();
    h_GbJet_dPhi[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_GJet_dR_%s",cut0[hist0].c_str());
    else sprintf(name, "h_GbJet_dR_%s",cut0[hist0].c_str());
    h_GbJet_dR[hist0] = new TH1F(name, "delta R dist", 100, 0.0, 10.0);
    h_GbJet_dR[hist0]->GetYaxis()->SetTitle("Events");          h_GbJet_dR[hist0]->GetYaxis()->CenterTitle();
    h_GbJet_dR[hist0]->GetXaxis()->SetTitle("#Delta R");        h_GbJet_dR[hist0]->GetXaxis()->CenterTitle();
    h_GbJet_dR[hist0]->Sumw2();
 
    sprintf(name, "h_cosThetaStar_%s",cut0[hist0].c_str());
    h_cosThetaStar[hist0] = new TH1F(name,"Cos theta star of photon+jet",50,0,1);
    h_cosThetaStar[hist0]->GetYaxis()->SetTitle("Events");        h_cosThetaStar[hist0]->GetYaxis()->CenterTitle();
    h_cosThetaStar[hist0]->GetXaxis()->SetTitle("cos#theta*");    h_cosThetaStar[hist0]->GetXaxis()->CenterTitle();
    h_cosThetaStar[hist0]->Sumw2();

    sprintf(name, "h_PFMet_%s",cut0[hist0].c_str());
    h_PFMet[hist0] = new TH1F(name, "PFMet distribution", 200,0.0,1000.0);
    h_PFMet[hist0]->GetYaxis()->SetTitle("Events/5 GeV");           h_PFMet[hist0]->GetYaxis()->CenterTitle();
    h_PFMet[hist0]->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");    h_PFMet[hist0]->GetXaxis()->CenterTitle();
    h_PFMet[hist0]->Sumw2();

    sprintf(name, "h_SumPFMet_%s",cut0[hist0].c_str());
    h_SumPFMet[hist0]  = new TH1F(name,"SumET PF Met distribution",80,0.0,4000.0);
    h_SumPFMet[hist0]->GetYaxis()->SetTitle("Events/50 GeV");            h_SumPFMet[hist0]->GetYaxis()->CenterTitle();
    h_SumPFMet[hist0]->GetXaxis()->SetTitle("#sum#slash{E}_{T} (GeV)");  h_SumPFMet[hist0]->GetXaxis()->CenterTitle();
    h_SumPFMet[hist0]->Sumw2();

    sprintf(name, "h_MetBySumMET_%s",cut0[hist0].c_str());
    h_MetBySumMET[hist0]  = new TH1F(name,"MET / SumET PF Met",50,0.0,1.0);
    h_MetBySumMET[hist0]->GetYaxis()->SetTitle("Events");        h_MetBySumMET[hist0]->GetYaxis()->CenterTitle();
    h_MetBySumMET[hist0]->GetXaxis()->SetTitle("#slash{E}_{T}/#sum#slash{E}_{T}"); h_MetBySumMET[hist0]->GetXaxis()->CenterTitle();
    h_MetBySumMET[hist0]->Sumw2();

    sprintf(name, "h_PFMetVsGJmass_%s",cut0[hist0].c_str());
    h_PFMetVsGJmass[hist0] = new TH2F(name, "PFMet Vs InvtMass", 1400, 0.0, 14000.0, 150, 0.0, 6000.0);
    h_PFMetVsGJmass[hist0]->GetYaxis()->SetTitle("PFMET (GeV)");    h_PFMetVsGJmass[hist0]->GetYaxis()->CenterTitle();
    h_PFMetVsGJmass[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)");  h_PFMetVsGJmass[hist0]->GetXaxis()->CenterTitle();
    h_PFMetVsGJmass[hist0]->Sumw2();

    sprintf(name, "h_PFMetOverSumEtVsGJmass_%s",cut0[hist0].c_str());
    h_PFMetOverSumEtVsGJmass[hist0] = new TH2F(name, "PFMet/SumEt Vs InvtMass", 1400, 0.0, 14000.0, 1000, 0.0, 10.0);
    h_PFMetOverSumEtVsGJmass[hist0]->GetYaxis()->SetTitle("PFMET/SumEt (GeV)");  h_PFMetOverSumEtVsGJmass[hist0]->GetYaxis()->CenterTitle();
    h_PFMetOverSumEtVsGJmass[hist0]->GetXaxis()->SetTitle("Invt Mass (GeV)");  h_PFMetOverSumEtVsGJmass[hist0]->GetXaxis()->CenterTitle();
    h_PFMetOverSumEtVsGJmass[hist0]->Sumw2();

    sprintf(name, "h_PtByZGMassVsZGMass_%s",cut0[hist0].c_str());
    h_PtByZGMassVsZGMass[hist0] = new TH2F(name, "PtPho/InvtMass Vs InvtMass", 280, 0.0, 7000.0, 100, 0.0, 10.0);
    h_PtByZGMassVsZGMass[hist0]->GetYaxis()->SetTitle("PtPho/InvtMass (GeV) ");    h_PtByZGMassVsZGMass[hist0]->GetYaxis()->CenterTitle();
    h_PtByZGMassVsZGMass[hist0]->GetXaxis()->SetTitle(" ZG Invt Mass (GeV)");  h_PtByZGMassVsZGMass[hist0]->GetXaxis()->CenterTitle();
    h_PtByZGMassVsZGMass[hist0]->Sumw2();

    sprintf(name, "h_dRVsZMass_%s",cut0[hist0].c_str());
    h_dRVsZMass[hist0] = new TH2F(name, "Delta R Vs Invt Z_Mass", 50, 50.0, 150.0, 100, 0.0, 10.0);
    h_dRVsZMass[hist0]->GetYaxis()->SetTitle(" Delta R");    h_dRVsZMass[hist0]->GetYaxis()->CenterTitle();
    h_dRVsZMass[hist0]->GetXaxis()->SetTitle(" Z Invt Mass (GeV)");  h_dRVsZMass[hist0]->GetXaxis()->CenterTitle();
    h_dRVsZMass[hist0]->Sumw2();

    sprintf(name, "h_MetByPhPt_%s",cut0[hist0].c_str());
    h_MetByPhPt[hist0]  = new TH1F(name,"MET / PTgamma",100,0.0,2.0);
    h_MetByPhPt[hist0]->GetYaxis()->SetTitle("Events");        h_MetByPhPt[hist0]->GetYaxis()->CenterTitle();
    h_MetByPhPt[hist0]->GetXaxis()->SetTitle("MET/P_{T}^{#gamma}"); h_MetByPhPt[hist0]->GetXaxis()->CenterTitle();
    h_MetByPhPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_PhPt_vs_JetPt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_PhPt_vs_bJetPt_%s",cut0[hist0].c_str());
    h_PhPt_vs_bJetPt[hist0] = new TH2F(name,"Pt of Photon vs Jet ",120,0.0,4800.0,120,0.0,4800.0);
    h_PhPt_vs_bJetPt[hist0]->GetYaxis()->SetTitle("P_{T}^{Jet}");   h_PhPt_vs_bJetPt[hist0]->GetYaxis()->CenterTitle();
    h_PhPt_vs_bJetPt[hist0]->GetXaxis()->SetTitle("P_{T}^{#gamma}");  h_PhPt_vs_bJetPt[hist0]->GetXaxis()->CenterTitle();
    h_PhPt_vs_bJetPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_PhEta_vs_JetEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_PhEta_vs_bJetEta_%s",cut0[hist0].c_str());
    h_PhEta_vs_bJetEta[hist0] = new TH2F(name,"eta of Photon vs Jet ",100,-2.5,2.5,100,-2.5,2.5);
    h_PhEta_vs_bJetEta[hist0]->GetYaxis()->SetTitle("#eta^{Jet}");     h_PhEta_vs_bJetEta[hist0]->GetYaxis()->CenterTitle();
    h_PhEta_vs_bJetEta[hist0]->GetXaxis()->SetTitle("#eta^{#gamma}");  h_PhEta_vs_bJetEta[hist0]->GetXaxis()->CenterTitle();
    h_PhEta_vs_bJetEta[hist0]->Sumw2();

    sprintf(name, "h_CSVv2Dist_%s",cut0[hist0].c_str());
    h_CSVv2Dist[hist0] = new TH1F(name, "CSVv2 bTagger Distribution", 500, 0.0, 5.0);
    h_CSVv2Dist[hist0]->GetYaxis()->SetTitle("Events");      h_CSVv2Dist[hist0]->GetYaxis()->CenterTitle();
    h_CSVv2Dist[hist0]->GetXaxis()->SetTitle("CSVv2");       h_CSVv2Dist[hist0]->GetXaxis()->CenterTitle();
    h_CSVv2Dist[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_CSVv2_vs_JetPt_%s",cut0[hist0].c_str());
    else sprintf(name, "h_CSVv2_vs_bJetPt_%s",cut0[hist0].c_str());
    h_CSVv2_vs_bJetPt[hist0] = new TH2F(name,"CSVv2 vs Jet Pt", 120, 0.0, 4800.0, 500, 0.0, 5.0);
    h_CSVv2_vs_bJetPt[hist0]->GetYaxis()->SetTitle("CSVv2");         h_CSVv2_vs_bJetPt[hist0]->GetYaxis()->CenterTitle();
    h_CSVv2_vs_bJetPt[hist0]->GetXaxis()->SetTitle("P_{T}^{Jet}");   h_CSVv2_vs_bJetPt[hist0]->GetXaxis()->CenterTitle();
    h_CSVv2_vs_bJetPt[hist0]->Sumw2();

    if(hist0 == 0 || hist0 == 1 || hist0 == 2) sprintf(name, "h_CSVv2_vs_JetEta_%s",cut0[hist0].c_str());
    else sprintf(name, "h_CSVv2_vs_bJetEta_%s",cut0[hist0].c_str());
    h_CSVv2_vs_bJetEta[hist0] = new TH2F(name,"CSVv2 vs Jet Eta", 100, -2.5, 2.5, 500, 0.0, 5.0);
    h_CSVv2_vs_bJetEta[hist0]->GetYaxis()->SetTitle("CSVv2");         h_CSVv2_vs_bJetEta[hist0]->GetYaxis()->CenterTitle();
    h_CSVv2_vs_bJetEta[hist0]->GetXaxis()->SetTitle("#eta^{Jet}");    h_CSVv2_vs_bJetEta[hist0]->GetXaxis()->CenterTitle();
    h_CSVv2_vs_bJetEta[hist0]->Sumw2();

    sprintf(name, "h_goodPV_%s",cut0[hist0].c_str());
    h_goodPV[hist0] = new TH1F(name,"No. of Good Primary Vertices", 70, 0, 70);
    h_goodPV[hist0]->GetYaxis()->SetTitle("Events");   h_goodPV[hist0]->GetYaxis()->CenterTitle();
    h_goodPV[hist0]->GetXaxis()->SetTitle("nPV");      h_goodPV[hist0]->GetXaxis()->CenterTitle();
    h_goodPV[hist0]->Sumw2();
       

    sprintf(name, "h_PtByZGMiso_%s",cut0[hist0].c_str());
    h_PtByZGMiso[hist0] = new TH1F(name,"photon Pt/zgamma mass", 20, 0, 2);
    h_PtByZGMiso[hist0]->GetYaxis()->SetTitle("Events");   h_PtByZGMiso[hist0]->GetYaxis()->CenterTitle();
    h_PtByZGMiso[hist0]->GetXaxis()->SetTitle("ratio");      h_PtByZGMiso[hist0]->GetXaxis()->CenterTitle();
    h_PtByZGMiso[hist0]->Sumw2();

    sprintf(name, "h_nIsoPhotons_%s",cut0[hist0].c_str());
    h_nIsoPhotons[hist0] = new TH1F(name,"No. of Isolated Photons", 400, 0, 200);
    h_nIsoPhotons[hist0]->GetYaxis()->SetTitle("Events");         h_nIsoPhotons[hist0]->GetYaxis()->CenterTitle();
    h_nIsoPhotons[hist0]->GetXaxis()->SetTitle("nIsoPhotons");    h_nIsoPhotons[hist0]->GetXaxis()->CenterTitle();
    h_nIsoPhotons[hist0]->Sumw2();

    sprintf(name, "h_nGoodPhotons_%s",cut0[hist0].c_str());
    h_nGoodPhotons[hist0] = new TH1F(name,"No. of Good Isolated Photons", 400, 0, 200);
    h_nGoodPhotons[hist0]->GetYaxis()->SetTitle("Events");          h_nGoodPhotons[hist0]->GetYaxis()->CenterTitle();
    h_nGoodPhotons[hist0]->GetXaxis()->SetTitle("nGoodPhotons");    h_nGoodPhotons[hist0]->GetXaxis()->CenterTitle();
    h_nGoodPhotons[hist0]->Sumw2();

    sprintf(name, "h_IsoPhotonIdxVsPt_%s",cut0[hist0].c_str());
    h_IsoPhotonIdxVsPt[hist0] = new TH2F(name,"Isolated Photons Index Vs Pt", 120, 0.0, 4800.0, 10, 0.0, 10.0);
    h_IsoPhotonIdxVsPt[hist0]->GetYaxis()->SetTitle("IsoPhotonIdx");      h_IsoPhotonIdxVsPt[hist0]->GetYaxis()->CenterTitle();
    h_IsoPhotonIdxVsPt[hist0]->GetXaxis()->SetTitle("P^{T}_{#gamma}");    h_IsoPhotonIdxVsPt[hist0]->GetXaxis()->CenterTitle();
    h_IsoPhotonIdxVsPt[hist0]->Sumw2();

    sprintf(name, "h_GoodPhotonIdxVsPt_%s",cut0[hist0].c_str());
    h_GoodPhotonIdxVsPt[hist0] = new TH2F(name,"Good Isolated Photons Idx Vs Pt", 120, 0.0, 4800.0, 10, 0.0, 10.0);
    h_GoodPhotonIdxVsPt[hist0]->GetYaxis()->SetTitle("GoodPhotonIdx");         h_GoodPhotonIdxVsPt[hist0]->GetYaxis()->CenterTitle();
    h_GoodPhotonIdxVsPt[hist0]->GetXaxis()->SetTitle("P^{T}_{#gamma}");       h_GoodPhotonIdxVsPt[hist0]->GetXaxis()->CenterTitle();
    h_GoodPhotonIdxVsPt[hist0]->Sumw2();

    sprintf(name, "h_nJets_%s",cut0[hist0].c_str());
    h_nJets[hist0] = new TH1F(name,"No. of Jets", 400, 0, 200);
    h_nJets[hist0]->GetYaxis()->SetTitle("Events");      h_nJets[hist0]->GetYaxis()->CenterTitle();
    h_nJets[hist0]->GetXaxis()->SetTitle("nJets");       h_nJets[hist0]->GetXaxis()->CenterTitle();
    h_nJets[hist0]->Sumw2();

    sprintf(name, "h_JetIdxVsPt_%s",cut0[hist0].c_str());
    h_JetIdxVsPt[hist0] = new TH2F(name,"Jet Idx vs Pt", 120, 0.0, 4800.0, 20, 0.0, 20.0);
    h_JetIdxVsPt[hist0]->GetYaxis()->SetTitle("Jet Idx");           h_JetIdxVsPt[hist0]->GetYaxis()->CenterTitle();
    h_JetIdxVsPt[hist0]->GetXaxis()->SetTitle("P^{T}_{Jet}");       h_JetIdxVsPt[hist0]->GetXaxis()->CenterTitle();
    h_JetIdxVsPt[hist0]->Sumw2();

  }

  h_event = new TH1F("h_event", "h_event", 1000000, 1, 1000000000);
  h_event->GetYaxis()->SetTitle("Events");    h_event->GetYaxis()->CenterTitle();
  h_event->GetXaxis()->SetTitle("Evts");      h_event->GetXaxis()->CenterTitle();
  h_event->Sumw2();


  //Position of photon and jet
  h_PC = new TH1F("h_PC", "Photon Candidate", 10, 0, 10);
  h_PC->GetYaxis()->SetTitle("Events");                       h_PC->GetYaxis()->CenterTitle();
  h_PC->GetXaxis()->SetTitle("Position of Photon");           h_PC->GetXaxis()->CenterTitle();
  h_PC->Sumw2();

  h_JC = new TH1F("h_JC", "Jet Candidate", 20, 0, 20);
  h_JC->GetYaxis()->SetTitle("Events");                       h_JC->GetYaxis()->CenterTitle();
  h_JC->GetXaxis()->SetTitle("Position of Jet");              h_JC->GetXaxis()->CenterTitle();
  h_JC->Sumw2();
 
/* 
  h_PassHLT175_EBPassProbes = new TH1F("h_PassHLT175_EBPassProbes", "Pt Distribution of Photons for events passing HLT_Photon175 in EB", 300, 0.0, 1500.0);
  h_PassHLT175_EBPassProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT175_EBPassProbes->GetYaxis()->CenterTitle();
  h_PassHLT175_EBPassProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT175_EBPassProbes->GetXaxis()->CenterTitle();
  h_PassHLT175_EBPassProbes->Sumw2();

  h_PassHLT90OR120_EBAllProbes = new TH1F("h_PassHLT90OR120_EBAllProbes", "Pt Distribution of Photons for events passing HLT_Photon90 OR HLT_Photon120 in EB", 300, 0.0, 1500.0);
  h_PassHLT90OR120_EBAllProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT90OR120_EBAllProbes->GetYaxis()->CenterTitle();
  h_PassHLT90OR120_EBAllProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT90OR120_EBAllProbes->GetXaxis()->CenterTitle();
  h_PassHLT90OR120_EBAllProbes->Sumw2();

  h_PassHLT175_EEPassProbes = new TH1F("h_PassHLT175_EEPassProbes", "Pt Distribution of Photons for events passing HLT_Photon175 in EE", 300, 0.0, 1500.0);
  h_PassHLT175_EEPassProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT175_EEPassProbes->GetYaxis()->CenterTitle();
  h_PassHLT175_EEPassProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT175_EEPassProbes->GetXaxis()->CenterTitle();
  h_PassHLT175_EEPassProbes->Sumw2();

  h_PassHLT90OR120_EEAllProbes = new TH1F("h_PassHLT90OR120_EEAllProbes", "Pt Distribution of Photons for events passing HLT_Photon90 OR HLT_Photon120 in EE", 300, 0.0, 1500.0);
  h_PassHLT90OR120_EEAllProbes->GetYaxis()->SetTitle("Events/5 GeV");          h_PassHLT90OR120_EEAllProbes->GetYaxis()->CenterTitle();
  h_PassHLT90OR120_EEAllProbes->GetXaxis()->SetTitle("P_{T}^{#gamma}");        h_PassHLT90OR120_EEAllProbes->GetXaxis()->CenterTitle();
  h_PassHLT90OR120_EEAllProbes->Sumw2();
*/
 
  //Trigger Turn-on
  const Int_t nTrigPtBins = 28; 
  const Double_t TrigPtBins[nTrigPtBins+1] = { 0, 50, 100, 120, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 300, 400, 500, 700, 1000};// 1500, 2000, 3000, 4000, 5000};       

  std::string cut1[2] = {"deno", "num"};
  for(Int_t hist1 =0; hist1 < 2; ++hist1){
  
    sprintf(name, "h_TrigPhotonPt_%s",cut1[hist1].c_str());
    h_TrigPhotonPt[hist1] = new TH1F(name, "pt of trigger photon", nTrigPtBins, TrigPtBins);
    h_TrigPhotonPt[hist1]->GetYaxis()->SetTitle("");      h_TrigPhotonPt[hist1]->GetYaxis()->CenterTitle();                      
    h_TrigPhotonPt[hist1]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");   h_TrigPhotonPt[hist1]->GetXaxis()->CenterTitle();
    h_TrigPhotonPt[hist1]->Sumw2();
  }
  h_TrigEff = new TH1F("h_TrigEff", "Trigger Efficiency" , nTrigPtBins , TrigPtBins);
  h_TrigEff->GetYaxis()->SetTitle("");   h_TrigEff->GetYaxis()->CenterTitle();
  h_TrigEff->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_TrigEff->GetXaxis()->CenterTitle();
  h_TrigEff->Sumw2(); 

/*  h_PassPhotonPt = new TH1F("h_PassPhotonPt", "Mass Plot" ,150,0.0,6000.0);
  h_PassPhotonPt->GetYaxis()->SetTitle("");  h_PassPhotonPt->GetYaxis()->CenterTitle();
  h_PassPhotonPt->GetXaxis()->SetTitle("Pt (Gev)"); h_PassPhotonPt->GetXaxis()->CenterTitle();
  h_PassPhotonPt->Sumw2();
  
  h_PassJetPt = new TH1F("h_PassJetPt", "Mass Plot" ,150,0.0,6000.0);
  h_PassJetPt->GetYaxis()->SetTitle("");  h_PassJetPt->GetYaxis()->CenterTitle();
  h_PassJetPt->GetXaxis()->SetTitle("Pt (Gev)"); h_PassJetPt->GetXaxis()->CenterTitle();
  h_PassJetPt->Sumw2();
*/
  h_TrigGJmass_Jetpt210 = new TH1F("h_TrigGJmass_Jetpt60", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt210->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt210->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt210->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt210->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt210->Sumw2();

  h_TrigGJmass_Jetpt220 = new TH1F("h_TrigGJmass_Jetpt100", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt220->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt220->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt220->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt220->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt220->Sumw2();

  h_TrigGJmass_Jetpt200 = new TH1F("h_TrigGJmass_Jetpt120", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt200->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt200->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt200->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt200->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt200->Sumw2();

  h_TrigGJmass_Jetpt170 = new TH1F("h_TrigGJmass_Jetpt150", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt170->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt170->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt170->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt170->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt170->Sumw2();

  h_TrigGJmass_Jetpt190 = new TH1F("h_TrigGJmass_Jetpt190", "Invt Mass", nMassBins, MassBin);
  h_TrigGJmass_Jetpt190->GetYaxis()->SetTitle("");   h_TrigGJmass_Jetpt190->GetYaxis()->CenterTitle();
  h_TrigGJmass_Jetpt190->GetXaxis()->SetTitle("Mass (GeV)");  h_TrigGJmass_Jetpt190->GetXaxis()->CenterTitle();
  h_TrigGJmass_Jetpt190->Sumw2();
  
  h_LeadPt = new TH1F("h_LeadPt", "Leading Pt" , nTrigPtBins , TrigPtBins);
  h_LeadPt->GetYaxis()->SetTitle(""); h_LeadPt->GetYaxis()->CenterTitle();
  h_LeadPt->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");  h_TrigEff->GetXaxis()->CenterTitle();
  h_LeadPt->Sumw2();

  //Defining the histogram from filling number of events after various cuts
  const int nbins_qstar = 11;
  const int nbins_bstar = 16;
  const int nbins_zgamma = 13;//ZGamma

  TString CutFlowLabels_qstar[nbins_qstar] = {"Total", "HLT", "PrimaryVtx_&_METFilters", "PhotonID", "PhotonPtEta", "JetID", "JetPt", "JetEta", "DPhi", "DEta", "MassCut"};
  TString CutFlowLabels_bstar[nbins_bstar] = {"Total", "HLT", "PrimaryVtx_&_METFilters", "PhotonID", "PhotonPtEta", "JetID", "JetPt", "JetEta", "1BTag", "1BTag_DPhi", "1BTag_DEta", "1BTag_MassCut", "0BTag", "0BTag_DPhi", "0BTag_DEta", "0BTag_MassCut"};
  TString CutFlowLabels_zgamma[nbins_zgamma] = {"Total", "HLT", "PrimaryVtx_", "EleId" ,"ElePt&Eta" , "EC1+EC2" ,"isZ","PhotonId" , "PhotonCandidate" ,"PixelSeedVeto", "DeltaR1&DeltaR2","ZGmassCut", "ptByMZGamma" };//ZGamma

  h_CutFlow_qstar = new TH1F("h_CutFlow_qstar", "Events Passing Various Cuts for qstar", nbins_qstar, 0, nbins_qstar);
  h_CutFlow_qstar->GetYaxis()->SetTitle("Events");         h_CutFlow_qstar->GetYaxis()->CenterTitle();
  for(int i = 0; i < nbins_qstar; i++){
    h_CutFlow_qstar->GetXaxis()->SetBinLabel(i+1, CutFlowLabels_qstar[i]);
  }

  h_CutFlow_ZGamma_Ele = new TH1F("h_CutFlow_ZGamma_Ele", "Events Passing Various Cuts for ZGamma", nbins_zgamma, 0, nbins_zgamma);
  h_CutFlow_ZGamma_Ele->GetYaxis()->SetTitle("Events");         h_CutFlow_ZGamma_Ele->GetYaxis()->CenterTitle();
  for(int i = 0; i < nbins_zgamma; i++){
    h_CutFlow_ZGamma_Ele->GetXaxis()->SetBinLabel(i+1, CutFlowLabels_zgamma[i]);
  }

  h_CutFlow_bstar = new TH1F("h_CutFlow_bstar", "Events Passing Various Cuts for bstar", nbins_bstar, 0, nbins_bstar);
  h_CutFlow_bstar->GetYaxis()->SetTitle("Events");         h_CutFlow_bstar->GetYaxis()->CenterTitle();
  for(int i = 0; i < nbins_bstar; i++){
    h_CutFlow_bstar->GetXaxis()->SetBinLabel(i+1, CutFlowLabels_bstar[i]);
  }

}
#endif // #ifdef PostAnalyzer_Data_cxx
 
EOF

cat > analysis_${filenameTag1}.C <<EOF
#include "PostAnalyzer_Data.C"
#include "TROOT.h"
int main(){
    PostAnalyzer_Data a;     
    a.Loop();
return 0;
}
EOF


###Compilation

g++ -Wno-deprecated analysis_${filenameTag1}.C -o ${filenameTag1}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs` 
####Execution
#./${filenameTag}.exe

###Submit jobs


chmod 775 dup_MakeCondorFiles_Selections.csh

source dup_MakeCondorFiles_Selections.csh ${filenameTag1}
##**************************************************
@ sf = ${sf} + ${r}
@ maxf = ${maxf} + ${r}

end ##end of while loop
##**************************************************
@ sampleIndex = ${sampleIndex} + 1
end ##end of for loop##
end ##end of k for loop
end ##end of case for loop
