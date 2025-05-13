#include "PhysicsTools/XGBoost/interface/XGBooster.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include <math.h>
#include "TMath.h"
#include <limits>
#include "TSystem.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "Math/GenVector/LorentzVector.h"

using namespace std;

//  *Notes on model implementation from Anagha*
//  Following the --> are the branch names in the ntuple, *new* are added by this module
//
//  mu-tau 1 bjet
//  'mT_thMET': massT(taus4V[:,0], events.MET), --> mtMET_1
//  'mT_mMET':  massT(muons4V[:,0], events.MET), --> mtMET_2
//  'd_zeta': d_zeta_ll(muons2V[:,0], taus2V[:,0], MET2V) --> d_zeta
//  While training for mu-tau channel I used my own ntuples, where the below method was used for mutau_pt. For etau and emu, pt_vis_nominal branch was used directly. But definition is the same.
//  'mutau_pt':(muons4V[:,0] + taus4V[:,0]).pt,  --> pt_vis_nominal
//  'm_b1mt':((taus4V[:,0] + bjet4VLV[:,0] + muons4V[:,0]).mass), --> m_btautau_vis_nominal
//  'bmt_dR':(bjets4V[:,0].delta_r(muons4V[:,0] + taus4V[:,0])),  --> btautau_dR   *new*
//  'b1th_dR':(bjets4V[:,0].delta_r(taus4V[:,0])), --> btau2_dR   *new*
//  'njets'
//  
//  mu-tau >1 bjet
//  'm_b2mt':(taus4V[:,0] + bjet4VLV[:,1] + muons4V[:,0]).mass, --> m_b2tautau_vis_nominal   *new*
//  'm_bbmt': (bjet4VLV[:,0] + bjet4VLV[:,1] + muons4V[:,0] + taus4V[:,0]).mass, --> m_bbtautau_vis_nominal   *new*
//  'b2th_dR': bjet4VLV[:,1].delta_r(taus4V[:,0]),  --> b2tau2_dR   *new*
//  'd_ma':((events.m_bb - events.m_mutau)/events.m_mutau) --> d_ma   *new*
//  
//  For e-tau and e-mu, here, I have only listed variables I computed, not the pre-existing branches in the ntuples
//  etau 1 bjet
//  'b1e_dR': delta_r(events.beta_deepflavour_1.to_numpy(), events.bphi_deepflavour_1.to_numpy(),  --> btau1_dR   *new*
//                         events.eta_1.to_numpy(), events.phi_1.to_numpy())
//  'b1th_dR': delta_r(events.beta_deepflavour_1.to_numpy(), events.bphi_deepflavour_1.to_numpy(), --> btau2_dR   *new*
//                          events.eta_2.to_numpy(), events.phi_2.to_numpy())
//  'mT_b1MET': MT between the 1st b jet and MET --> mtMET_b   *new*
//  
//  etau >1 bjet
//  ## BUG! The following variable should be: 'b2e_dR' since its the dR between e and sub-leading b
//  'b2th_dR':delta_r(events.beta_deepflavour_2.to_numpy(), events.bphi_deepflavour_2.to_numpy(),  --> b2tau1_dR   *new*
//                         events.eta_1.to_numpy(), events.phi_1.to_numpy())
//  'd_ma':(events['m_b1b2']- events['m_vis_nominal'])/events['m_vis_nominal']  --> d_ma   *new*
//  'mbb': invariant mass of bb --> m_bb   *new*
//  
//  emu 1 bjet
//  'b1emu_dR':delta_r(em_eta.to_numpy(), em_phi.to_numpy(),
//                          events.beta_deepflavour_1.to_numpy(), events.bphi_deepflavour_1.to_numpy()) --> btautau_dR   *new*
//  
//  emu >1 bjet
//  'b2emu_dR':delta_r(em_eta.to_numpy(), em_phi.to_numpy(),
//                          events.beta_deepflavour_2.to_numpy(), events.bphi_deepflavour_2.to_numpy()) --> b2tautau_dR   *new*
//  'd_ma':((b1b2_m- events.m_vis_nominal)/(events.m_vis_nominal)) --> d_ma   *new*

struct CalculatedFeatures {
  float btautau_dR = -10;
  float btau1_dR = -10;
  float btau2_dR = -10;
  float mtMET_b = -10;
  float m_bb = -10;
  float m_bbtautau_vis = -10;
  float m_b2tautau_vis = -10;
  float d_ma = -10;
  float b2tau1_dR = -10;
  float b2tau2_dR = -10;
  float b2tautau_dR = -10;
};

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], pat::XGBooster boosters[], bool isMC, bool isEmbedded);
void runBDTeval(float& BDTval, CalculatedFeatures& features_to_save,
  int booster_idx, pat::XGBooster booster, float jet_1_pt, float jet_1_eta, 
  float jet_1_phi, float jet_1_mass, float jet_2_pt, float jet_2_eta, 
  float jet_2_phi, float jet_2_mass, float pt_1, float eta_1, 
  float phi_1, float mass_1, float pt_2, float eta_2, 
  float phi_2, float mass_2, float met, float met_phi, 
  float d_zeta, int njets, float m_btautau_vis, float mtMET_1, 
  float mtMET_2, float pt_vis
);
void CopyDir(TDirectory *source, optutl::CommandLineParser parser);
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void copyFiles(optutl::CommandLineParser parser, TFile* fOld, TFile* fNew);

float compute_deltaPhi(float v1, float v2) {
  float r = std::fmod(v2 - v1, 2*M_PI);
  if (r < -M_PI) {
    r += 2*M_PI;
  }
  else if (r >= M_PI) {
    r -= 2*M_PI;
  }
  return r;
}

float compute_deltaR(float e1, float e2, float p1, float p2) {
  float deta = e1 - e2;
  float dphi = compute_deltaPhi(p1, p2);
  return std::sqrt(deta*deta + dphi*dphi);
}

float compute_mt(float pt_1, float phi_1, float pt_met, float phi_met) {
  const auto dphi = compute_deltaPhi(phi_1, phi_met);
  return std::sqrt(2.0 * pt_1 * pt_met * (1.0 - std::cos(dphi)));
}

int main (int argc, char** argv) {
  optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
  parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
  parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
  parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);
  parser.addOption("isMC",optutl::CommandLineParser::kBool,"isMC",false);
  parser.addOption("isEmbedded",optutl::CommandLineParser::kBool,"isEmbedded",false);

  parser.parseArguments (argc, argv);

  pat::XGBooster mutau_1bjet_booster("/afs/cern.ch/work/p/pdas/public/Ha1a2_bdtfiles/mutau-1bjet.model");
  std::cout<<"mutau-1bjet.model loaded"<<std::endl;
  mutau_1bjet_booster.addFeature("mT_thMET");
  mutau_1bjet_booster.addFeature("mT_mMET");
  mutau_1bjet_booster.addFeature("d_zeta");
  mutau_1bjet_booster.addFeature("mutau_pt");
  mutau_1bjet_booster.addFeature("m_b1mt");
  mutau_1bjet_booster.addFeature("bmt_dR");
  mutau_1bjet_booster.addFeature("b1th_dR");
  mutau_1bjet_booster.addFeature("njets");

  pat::XGBooster etau_1bjet_booster("/afs/cern.ch/work/p/pdas/public/Ha1a2_bdtfiles/etau-1bjet.model");
  std::cout<<"etau-1bjet.model loaded"<<std::endl;
  etau_1bjet_booster.addFeature("pt_vis_nominal");
  etau_1bjet_booster.addFeature("D_zeta_nominal");
  etau_1bjet_booster.addFeature("b1e_dR");
  etau_1bjet_booster.addFeature("m_btautau_vis_nominal");
  etau_1bjet_booster.addFeature("mtMET_1_nominal");
  etau_1bjet_booster.addFeature("b1th_dR");
  etau_1bjet_booster.addFeature("mtMET_2_nominal");
  etau_1bjet_booster.addFeature("njets");

  pat::XGBooster emu_1bjet_booster("/afs/cern.ch/work/p/pdas/public/Ha1a2_bdtfiles/emu-1bjet.model");
  std::cout<<"emu-1bjet.model loaded"<<std::endl;
  emu_1bjet_booster.addFeature("pt_vis_nominal");
  emu_1bjet_booster.addFeature("D_zeta_nominal");
  emu_1bjet_booster.addFeature("mT_b1MET");
  emu_1bjet_booster.addFeature("pt_1_nominal");
  emu_1bjet_booster.addFeature("m_btautau_vis_nominal");
  emu_1bjet_booster.addFeature("mtMET_1_nominal");
  emu_1bjet_booster.addFeature("b1emu_dR");
  emu_1bjet_booster.addFeature("njets");

  pat::XGBooster mutau_2bjet_booster("/afs/cern.ch/work/p/pdas/public/Ha1a2_bdtfiles/mutau_morethan1b.model");
  std::cout<<"mutau_morethan1b.model loaded"<<std::endl;
  mutau_2bjet_booster.addFeature("m_b2mt");
  mutau_2bjet_booster.addFeature("bmt_dR");
  mutau_2bjet_booster.addFeature("b2th_dR");
  mutau_2bjet_booster.addFeature("mT_mMET");
  mutau_2bjet_booster.addFeature("m_bbmt");
  mutau_2bjet_booster.addFeature("d_ma");

  pat::XGBooster etau_2bjet_booster("/afs/cern.ch/work/p/pdas/public/Ha1a2_bdtfiles/etau-morethan1b.model");
  std::cout<<"etau-morethan1b.model loaded"<<std::endl;
  etau_2bjet_booster.addFeature("mbb");
  etau_2bjet_booster.addFeature("m_btautau_vis_nominal");
  etau_2bjet_booster.addFeature("mtMET_1_nominal");
  etau_2bjet_booster.addFeature("b1th_dR");
  etau_2bjet_booster.addFeature("b1e_dR");
  etau_2bjet_booster.addFeature("b2th_dR");

  pat::XGBooster emu_2bjet_booster("/afs/cern.ch/work/p/pdas/public/Ha1a2_bdtfiles/emu-morethan1b.model");
  std::cout<<"emu-morethan1b.model loaded"<<std::endl;
  emu_2bjet_booster.addFeature("pt_vis_nominal");
  emu_2bjet_booster.addFeature("mT_b1MET");
  emu_2bjet_booster.addFeature("mtMET_1_nominal");
  emu_2bjet_booster.addFeature("b2emu_dR");
  emu_2bjet_booster.addFeature("d_ma");

  pat::XGBooster boosters[6] = {mutau_1bjet_booster, etau_1bjet_booster, emu_1bjet_booster, mutau_2bjet_booster, etau_2bjet_booster, emu_2bjet_booster};

  char TreeToUse[80]="first" ;

  TFile *fProduce;

  TFile *f = TFile::Open(parser.stringValue("inputFile").c_str());
  std::cout<<"Creating new outputfile"<<std::endl;
  std::string newFileName = parser.stringValue("newFile");
  bool isMC = parser.boolValue("isMC");
  bool isEmbedded = parser.boolValue("isEmbedded");

  fProduce = new TFile(newFileName.c_str(),"RECREATE");
  copyFiles(parser, f, fProduce);
  fProduce = new TFile(newFileName.c_str(),"UPDATE");
  std::cout<<"listing the directories================="<<std::endl;
  fProduce->ls();
  readdir(fProduce,parser,TreeToUse,boosters,isMC,isEmbedded);
  fProduce->Close();
  f->Close();
}

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], pat::XGBooster boosters[], bool isMC, bool isEmbedded) {
  TDirectory *dirsav = gDirectory;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  char stringA[80]="first";
  dir->cd();
  int k=0;
  bool got_mutau_tree(false);
  bool got_etau_tree(false);
  bool got_emu_tree(false);
  while ((key = (TKey*)next())) {
    printf("Found key=%s \n",key->GetName());

    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      sprintf(TreeToUse,"%s",key->GetName());
      if (std::string(key->GetName()).find("_tree")) {
         readdir(subdir,parser,TreeToUse,boosters,isMC,isEmbedded);
      }
      dirsav->cd();
    }
    // else if(k<1 && obj->IsA()->InheritsFrom(TTree::Class())) {
    // Check all items in the tree: do not do k<1 criteria
    else if (obj->IsA()->InheritsFrom(TTree::Class())) {
      k++;
      sprintf(TreeToUse,"%s",obj->GetName());
      if (got_mutau_tree && std::string(key->GetName()).find("mutau") != std::string::npos) continue;
      if (got_etau_tree && std::string(key->GetName()).find("etau") != std::string::npos) continue;
      if (got_emu_tree && std::string(key->GetName()).find("emu") != std::string::npos) continue;
      TTree *t = (TTree*)obj;

      // add the new branches
      float BDTval_nominal = -10;
      float BDTval_es1Up = -10;
      float BDTval_es1Down = -10;
      float BDTval_es2Up = -10;
      float BDTval_es2Down = -10;
      float BDTval_UESUp = -10;
      float BDTval_UESDown = -10;
      float BDTval_ResponseUp = -10;
      float BDTval_ResponseDown = -10;
      float BDTval_ResolutionUp = -10;
      float BDTval_ResolutionDown = -10;
      float BDTval_JetAbsoluteUp = -10;
      float BDTval_JetAbsoluteDown = -10;
      float BDTval_JetAbsoluteyearUp = -10;
      float BDTval_JetAbsoluteyearDown = -10;
      float BDTval_JetBBEC1Up = -10;
      float BDTval_JetBBEC1Down = -10;
      float BDTval_JetBBEC1yearUp = -10;
      float BDTval_JetBBEC1yearDown = -10;
      float BDTval_JetEC2Up = -10;
      float BDTval_JetEC2Down = -10;
      float BDTval_JetEC2yearUp = -10;
      float BDTval_JetEC2yearDown = -10;
      float BDTval_JetFlavorQCDUp = -10;
      float BDTval_JetFlavorQCDDown = -10;
      float BDTval_JetHFUp = -10;
      float BDTval_JetHFDown = -10;
      float BDTval_JetHFyearUp = -10;
      float BDTval_JetHFyearDown = -10;
      float BDTval_JetRelativeBalUp = -10;
      float BDTval_JetRelativeBalDown = -10;
      float BDTval_JetRelativeSampleUp = -10;
      float BDTval_JetRelativeSampleDown = -10;
      float BDTval_JERUp = -10;
      float BDTval_JERDown = -10;

      CalculatedFeatures features_nominal;
      CalculatedFeatures features_es1Up;
      CalculatedFeatures features_es1Down;
      CalculatedFeatures features_es2Up;
      CalculatedFeatures features_es2Down;
      CalculatedFeatures features_UESUp;
      CalculatedFeatures features_UESDown;
      CalculatedFeatures features_responseUp;
      CalculatedFeatures features_responseDown;
      CalculatedFeatures features_resolutionUp;
      CalculatedFeatures features_resolutionDown;
      CalculatedFeatures features_JetAbsoluteUp;
      CalculatedFeatures features_JetAbsoluteDown;
      CalculatedFeatures features_JetAbsoluteyearUp;
      CalculatedFeatures features_JetAbsoluteyearDown;
      CalculatedFeatures features_JetBBEC1Up;
      CalculatedFeatures features_JetBBEC1Down;
      CalculatedFeatures features_JetBBEC1yearUp;
      CalculatedFeatures features_JetBBEC1yearDown;
      CalculatedFeatures features_JetEC2Up;
      CalculatedFeatures features_JetEC2Down;
      CalculatedFeatures features_JetEC2yearUp;
      CalculatedFeatures features_JetEC2yearDown;
      CalculatedFeatures features_JetFlavorQCDUp;
      CalculatedFeatures features_JetFlavorQCDDown;
      CalculatedFeatures features_JetHFUp;
      CalculatedFeatures features_JetHFDown;
      CalculatedFeatures features_JetHFyearUp;
      CalculatedFeatures features_JetHFyearDown;
      CalculatedFeatures features_JetRelativeBalUp;
      CalculatedFeatures features_JetRelativeBalDown;
      CalculatedFeatures features_JetRelativeSampleUp;
      CalculatedFeatures features_JetRelativeSampleDown;
      CalculatedFeatures features_JERUp;
      CalculatedFeatures features_JERDown;
      
      TBranch *newBranch1 = t->Branch("bdtscore", &BDTval_nominal, "bdtscore/F");
      TBranch *newBranch1_es1U = t->Branch("bdtscore_es1Up",   &BDTval_es1Up,   "bdtscore_es1Up/F");
      TBranch *newBranch1_es1D = t->Branch("bdtscore_es1Down",   &BDTval_es1Down,   "bdtscore_es1Down/F");
      TBranch *newBranch1_es2U = t->Branch("bdtscore_es2Up", &BDTval_es2Up, "bdtscore_es2Up/F");
      TBranch *newBranch1_es2D = t->Branch("bdtscore_es2Down",   &BDTval_es2Down,   "bdtscore_es2Down/F");
      TBranch *newBranch1_UESU = t->Branch("bdtscore_UESUp", &BDTval_UESUp, "bdtscore_UESUp/F");
      TBranch *newBranch1_UESD = t->Branch("bdtscore_UESDown", &BDTval_UESDown, "bdtscore_UESDown/F");
      TBranch *newBranch1_ResponseU = t->Branch("bdtscore_responseUp", &BDTval_ResponseUp, "bdtscore_ResponseUp/F");
      TBranch *newBranch1_ResponseD = t->Branch("bdtscore_responseDown", &BDTval_ResponseDown, "bdtscore_ResponseDown/F");
      TBranch *newBranch1_ResolutionU = t->Branch("bdtscore_resolutionUp", &BDTval_ResolutionUp, "bdtscore_ResolutionUp/F");
      TBranch *newBranch1_ResolutionD = t->Branch("bdtscore_resolutionDown", &BDTval_ResolutionDown, "bdtscore_ResolutionDown/F");
      TBranch *newBranch1_JetAbsoluteU = t->Branch("bdtscore_JetAbsoluteUp", &BDTval_JetAbsoluteUp, "bdtscore_JetAbsoluteUp/F");
      TBranch *newBranch1_JetAbsoluteD = t->Branch("bdtscore_JetAbsoluteDown", &BDTval_JetAbsoluteDown, "bdtscore_JetAbsoluteDown/F");
      TBranch *newBranch1_JetAbsoluteyearU = t->Branch("bdtscore_JetAbsoluteyearUp", &BDTval_JetAbsoluteyearUp, "bdtscore_JetAbsoluteyearUp/F");
      TBranch *newBranch1_JetAbsoluteyearD = t->Branch("bdtscore_JetAbsoluteyearDown", &BDTval_JetAbsoluteyearDown, "bdtscore_JetAbsoluteyearDown/F");
      TBranch *newBranch1_JetBBEC1U = t->Branch("bdtscore_JetBBEC1Up", &BDTval_JetBBEC1Up, "bdtscore_JetBBEC1Up/F");
      TBranch *newBranch1_JetBBEC1D = t->Branch("bdtscore_JetBBEC1Down", &BDTval_JetBBEC1Down, "bdtscore_JetBBEC1Down/F");
      TBranch *newBranch1_JetBBEC1yearU = t->Branch("bdtscore_JetBBEC1yearUp", &BDTval_JetBBEC1yearUp, "bdtscore_JetBBEC1yearUp/F");
      TBranch *newBranch1_JetBBEC1yearD = t->Branch("bdtscore_JetBBEC1yearDown", &BDTval_JetBBEC1yearDown, "bdtscore_JetBBEC1yearDown/F");
      TBranch *newBranch1_JetEC2U = t->Branch("bdtscore_JetEC2Up", &BDTval_JetEC2Up, "bdtscore_JetEC2Up/F");
      TBranch *newBranch1_JetEC2D = t->Branch("bdtscore_JetEC2Down", &BDTval_JetEC2Down, "bdtscore_JetEC2Down/F");
      TBranch *newBranch1_JetEC2yearU = t->Branch("bdtscore_JetEC2yearUp", &BDTval_JetEC2yearUp, "bdtscore_JetEC2yearUp/F");
      TBranch *newBranch1_JetEC2yearD = t->Branch("bdtscore_JetEC2yearDown", &BDTval_JetEC2yearDown, "bdtscore_JetEC2yearDown/F");
      TBranch *newBranch1_JetFlavorQCDU = t->Branch("bdtscore_JetFlavorQCDUp", &BDTval_JetFlavorQCDUp, "bdtscore_JetFlavorQCDUp/F");
      TBranch *newBranch1_JetFlavorQCDD = t->Branch("bdtscore_JetFlavorQCDDown", &BDTval_JetFlavorQCDDown, "bdtscore_JetFlavorQCDDown/F");
      TBranch *newBranch1_JetHFU = t->Branch("bdtscore_JetHFUp", &BDTval_JetHFUp, "bdtscore_JetHFUp/F");
      TBranch *newBranch1_JetHFD = t->Branch("bdtscore_JetHFDown", &BDTval_JetHFDown, "bdtscore_JetHFDown/F");
      TBranch *newBranch1_JetHFyearU = t->Branch("bdtscore_JetHFyearUp", &BDTval_JetHFyearUp, "bdtscore_JetHFyearUp/F");
      TBranch *newBranch1_JetHFyearD = t->Branch("bdtscore_JetHFyearDown", &BDTval_JetHFyearDown, "bdtscore_JetHFyearDown/F");
      TBranch *newBranch1_JetRelativeBalU = t->Branch("bdtscore_JetRelativeBalUp", &BDTval_JetRelativeBalUp, "bdtscore_JetRelativeBalUp/F");
      TBranch *newBranch1_JetRelativeBalD = t->Branch("bdtscore_JetRelativeBalDown", &BDTval_JetRelativeBalDown, "bdtscore_JetRelativeBalDown/F");
      TBranch *newBranch1_JetRelativeSampleU = t->Branch("bdtscore_JetRelativeSampleUp", &BDTval_JetRelativeSampleUp, "bdtscore_JetRelativeSampleUp/F");
      TBranch *newBranch1_JetRelativeSampleD = t->Branch("bdtscore_JetRelativeSampleDown", &BDTval_JetRelativeSampleDown, "bdtscore_JetRelativeSampleDown/F");
      TBranch *newBranch1_JERU = t->Branch("bdtscore_JERUp", &BDTval_JERUp, "bdtscore_JERUp/F");
      TBranch *newBranch1_JERD = t->Branch("bdtscore_JERDown", &BDTval_JERDown, "bdtscore_JERDown/F");


      // nominal features
      TBranch *newBranch1_btautau_dR = t->Branch("btautau_dR", &features_nominal.btautau_dR, "btautau_dR/F");
      TBranch *newBranch1_btau1_dR = t->Branch("btau1_dR", &features_nominal.btau1_dR, "btau1_dR/F");
      TBranch *newBranch1_btau2_dR = t->Branch("btau2_dR", &features_nominal.btau2_dR, "btau2_dR/F");
      TBranch *newBranch1_mtMET_b_nominal = t->Branch("mtMET_b_nominal", &features_nominal.mtMET_b, "mtMET_b_nominal/F");
      TBranch *newBranch1_m_bb_nominal = t->Branch("m_bb_nominal", &features_nominal.m_bb, "m_bb_nominal/F");
      TBranch *newBranch1_m_bbtautau_vis_nominal = t->Branch("m_bbtautau_vis_nominal", &features_nominal.m_bbtautau_vis, "m_bbtautau_vis_nominal/F");
      TBranch *newBranch1_m_b2tautau_vis_nominal = t->Branch("m_b2tautau_vis_nominal", &features_nominal.m_b2tautau_vis, "m_b2tautau_vis_nominal/F");
      TBranch *newBranch1_d_ma_nominal= t->Branch("d_ma_nominal", &features_nominal.d_ma, "d_ma_nominal/F");
      TBranch *newBranch1_b2tau1_dR = t->Branch("b2tau1_dR", &features_nominal.b2tau1_dR, "b2tau1_dR/F");
      TBranch *newBranch1_b2tau2_dR = t->Branch("b2tau2_dR", &features_nominal.b2tau2_dR, "b2tau2_dR/F");
      TBranch *newBranch1_b2tautau_dR = t->Branch("b2tautau_dR", &features_nominal.b2tautau_dR, "b2tautau_dR/F");
      // es1
      TBranch *newBranch1_m_bbtautau_vis_es1Up = t->Branch("m_bbtautau_vis_es1Up", &features_es1Up.m_bbtautau_vis, "m_bbtautau_vis_es1Up/F");
      TBranch *newBranch1_m_b2tautau_vis_es1Up = t->Branch("m_b2tautau_vis_es1Up", &features_es1Up.m_b2tautau_vis, "m_b2tautau_vis_es1Up/F");
      TBranch *newBranch1_d_ma_es1Up = t->Branch("d_ma_es1Up", &features_es1Up.d_ma, "d_ma_es1Up/F");
      TBranch *newBranch1_m_bbtautau_vis_es1Down = t->Branch("m_bbtautau_vis_es1Down", &features_es1Down.m_bbtautau_vis, "m_bbtautau_vis_es1Down/F");
      TBranch *newBranch1_m_b2tautau_vis_es1Down = t->Branch("m_b2tautau_vis_es1Down", &features_es1Down.m_b2tautau_vis, "m_b2tautau_vis_es1Down/F");
      TBranch *newBranch1_d_ma_es1Down = t->Branch("d_ma_es1Down", &features_es1Down.d_ma, "d_ma_es1Down/F");
      // es1
      TBranch *newBranch1_m_bbtautau_vis_es2Up = t->Branch("m_bbtautau_vis_es2Up", &features_es2Up.m_bbtautau_vis, "m_bbtautau_vis_es2Up/F");
      TBranch *newBranch1_m_b2tautau_vis_es2Up = t->Branch("m_b2tautau_vis_es2Up", &features_es2Up.m_b2tautau_vis, "m_b2tautau_vis_es2Up/F");
      TBranch *newBranch1_d_ma_es2Up = t->Branch("d_ma_es2Up", &features_es2Up.d_ma, "d_ma_es2Up/F");
      TBranch *newBranch1_m_bbtautau_vis_es2Down = t->Branch("m_bbtautau_vis_es2Down", &features_es2Down.m_bbtautau_vis, "m_bbtautau_vis_es2Down/F");
      TBranch *newBranch1_m_b2tautau_vis_es2Down = t->Branch("m_b2tautau_vis_es2Down", &features_es2Down.m_b2tautau_vis, "m_b2tautau_vis_es2Down/F");
      TBranch *newBranch1_d_ma_es2Down = t->Branch("d_ma_es2Down", &features_es2Down.d_ma, "d_ma_es2Down/F");
      // MET systs
      TBranch *newBranch1_mtMET_b_UESUp = t->Branch("mtMET_b_UESUp", &features_UESUp.mtMET_b, "mtMET_b_UESUp/F");
      TBranch *newBranch1_mtMET_b_UESDown = t->Branch("mtMET_b_UESDown", &features_UESDown.mtMET_b, "mtMET_b_UESDown/F");
      TBranch *newBranch1_mtMET_b_responseUp = t->Branch("mtMET_b_responseUp", &features_responseUp.mtMET_b, "mtMET_b_responseUp/F");
      TBranch *newBranch1_mtMET_b_responseDown = t->Branch("mtMET_b_responseDown", &features_responseDown.mtMET_b, "mtMET_b_responseDown/F");
      TBranch *newBranch1_mtMET_b_resolutionUp = t->Branch("mtMET_b_resolutionUp", &features_resolutionUp.mtMET_b, "mtMET_b_resolutionUp/F");
      TBranch *newBranch1_mtMET_b_resolutionDown = t->Branch("mtMET_b_resolutionDown", &features_resolutionDown.mtMET_b, "mtMET_b_resolutionDown/F");
      // JetAbsolute
      TBranch *newBranch1_mtMET_b_JetAbsoluteUp = t->Branch("mtMET_b_JetAbsoluteUp", &features_JetAbsoluteUp.mtMET_b, "mtMET_b_JetAbsoluteUp/F");
      TBranch *newBranch1_m_bb_JetAbsoluteUp = t->Branch("m_bb_JetAbsoluteUp", &features_JetAbsoluteUp.m_bb, "m_bb_JetAbsoluteUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetAbsoluteUp = t->Branch("m_bbtautau_vis_JetAbsoluteUp", &features_JetAbsoluteUp.m_bbtautau_vis, "m_bbtautau_vis_JetAbsoluteUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetAbsoluteUp = t->Branch("m_b2tautau_vis_JetAbsoluteUp", &features_JetAbsoluteUp.m_b2tautau_vis, "m_b2tautau_vis_JetAbsoluteUp/F");
      TBranch *newBranch1_d_ma_JetAbsoluteUp = t->Branch("d_ma_JetAbsoluteUp", &features_JetAbsoluteUp.d_ma, "d_ma_JetAbsoluteUp/F");
      TBranch *newBranch1_mtMET_b_JetAbsoluteDown = t->Branch("mtMET_b_JetAbsoluteDown", &features_JetAbsoluteDown.mtMET_b, "mtMET_b_JetAbsoluteDown/F");
      TBranch *newBranch1_m_bb_JetAbsoluteDown = t->Branch("m_bb_JetAbsoluteDown", &features_JetAbsoluteDown.m_bb, "m_bb_JetAbsoluteDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetAbsoluteDown = t->Branch("m_bbtautau_vis_JetAbsoluteDown", &features_JetAbsoluteDown.m_bbtautau_vis, "m_bbtautau_vis_JetAbsoluteDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetAbsoluteDown = t->Branch("m_b2tautau_vis_JetAbsoluteDown", &features_JetAbsoluteDown.m_b2tautau_vis, "m_b2tautau_vis_JetAbsoluteDown/F");
      TBranch *newBranch1_d_ma_JetAbsoluteDown = t->Branch("d_ma_JetAbsoluteDown", &features_JetAbsoluteDown.d_ma, "d_ma_JetAbsoluteDown/F");
      // JetAbsoluteyear
      TBranch *newBranch1_mtMET_b_JetAbsoluteyearUp = t->Branch("mtMET_b_JetAbsoluteyearUp", &features_JetAbsoluteyearUp.mtMET_b, "mtMET_b_JetAbsoluteyearUp/F");
      TBranch *newBranch1_m_bb_JetAbsoluteyearUp = t->Branch("m_bb_JetAbsoluteyearUp", &features_JetAbsoluteyearUp.m_bb, "m_bb_JetAbsoluteyearUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetAbsoluteyearUp = t->Branch("m_bbtautau_vis_JetAbsoluteyearUp", &features_JetAbsoluteyearUp.m_bbtautau_vis, "m_bbtautau_vis_JetAbsoluteyearUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetAbsoluteyearUp = t->Branch("m_b2tautau_vis_JetAbsoluteyearUp", &features_JetAbsoluteyearUp.m_b2tautau_vis, "m_b2tautau_vis_JetAbsoluteyearUp/F");
      TBranch *newBranch1_d_ma_JetAbsoluteyearUp = t->Branch("d_ma_JetAbsoluteyearUp", &features_JetAbsoluteyearUp.d_ma, "d_ma_JetAbsoluteyearUp/F");
      TBranch *newBranch1_mtMET_b_JetAbsoluteyearDown = t->Branch("mtMET_b_JetAbsoluteyearDown", &features_JetAbsoluteyearDown.mtMET_b, "mtMET_b_JetAbsoluteyearDown/F");
      TBranch *newBranch1_m_bb_JetAbsoluteyearDown = t->Branch("m_bb_JetAbsoluteyearDown", &features_JetAbsoluteyearDown.m_bb, "m_bb_JetAbsoluteyearDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetAbsoluteyearDown = t->Branch("m_bbtautau_vis_JetAbsoluteyearDown", &features_JetAbsoluteyearDown.m_bbtautau_vis, "m_bbtautau_vis_JetAbsoluteyearDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetAbsoluteyearDown = t->Branch("m_b2tautau_vis_JetAbsoluteyearDown", &features_JetAbsoluteyearDown.m_b2tautau_vis, "m_b2tautau_vis_JetAbsoluteyearDown/F");
      TBranch *newBranch1_d_ma_JetAbsoluteyearDown = t->Branch("d_ma_JetAbsoluteyearDown", &features_JetAbsoluteyearDown.d_ma, "d_ma_JetAbsoluteyearDown/F");
      // JetBBEC1
      TBranch *newBranch1_mtMET_b_JetBBEC1Up = t->Branch("mtMET_b_JetBBEC1Up", &features_JetBBEC1Up.mtMET_b, "mtMET_b_JetBBEC1Up/F");
      TBranch *newBranch1_m_bb_JetBBEC1Up = t->Branch("m_bb_JetBBEC1Up", &features_JetBBEC1Up.m_bb, "m_bb_JetBBEC1Up/F");
      TBranch *newBranch1_m_bbtautau_vis_JetBBEC1Up = t->Branch("m_bbtautau_vis_JetBBEC1Up", &features_JetBBEC1Up.m_bbtautau_vis, "m_bbtautau_vis_JetBBEC1Up/F");
      TBranch *newBranch1_m_b2tautau_vis_JetBBEC1Up = t->Branch("m_b2tautau_vis_JetBBEC1Up", &features_JetBBEC1Up.m_b2tautau_vis, "m_b2tautau_vis_JetBBEC1Up/F");
      TBranch *newBranch1_d_ma_JetBBEC1Up = t->Branch("d_ma_JetBBEC1Up", &features_JetBBEC1Up.d_ma, "d_ma_JetBBEC1Up/F");
      TBranch *newBranch1_mtMET_b_JetBBEC1Down = t->Branch("mtMET_b_JetBBEC1Down", &features_JetBBEC1Down.mtMET_b, "mtMET_b_JetBBEC1Down/F");
      TBranch *newBranch1_m_bb_JetBBEC1Down = t->Branch("m_bb_JetBBEC1Down", &features_JetBBEC1Down.m_bb, "m_bb_JetBBEC1Down/F");
      TBranch *newBranch1_m_bbtautau_vis_JetBBEC1Down = t->Branch("m_bbtautau_vis_JetBBEC1Down", &features_JetBBEC1Down.m_bbtautau_vis, "m_bbtautau_vis_JetBBEC1Down/F");
      TBranch *newBranch1_m_b2tautau_vis_JetBBEC1Down = t->Branch("m_b2tautau_vis_JetBBEC1Down", &features_JetBBEC1Down.m_b2tautau_vis, "m_b2tautau_vis_JetBBEC1Down/F");
      TBranch *newBranch1_d_ma_JetBBEC1Down = t->Branch("d_ma_JetBBEC1Down", &features_JetBBEC1Down.d_ma, "d_ma_JetBBEC1Down/F");
      // JetBBEC1year
      TBranch *newBranch1_mtMET_b_JetBBEC1yearUp = t->Branch("mtMET_b_JetBBEC1yearUp", &features_JetBBEC1yearUp.mtMET_b, "mtMET_b_JetBBEC1yearUp/F");
      TBranch *newBranch1_m_bb_JetBBEC1yearUp = t->Branch("m_bb_JetBBEC1yearUp", &features_JetBBEC1yearUp.m_bb, "m_bb_JetBBEC1yearUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetBBEC1yearUp = t->Branch("m_bbtautau_vis_JetBBEC1yearUp", &features_JetBBEC1yearUp.m_bbtautau_vis, "m_bbtautau_vis_JetBBEC1yearUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetBBEC1yearUp = t->Branch("m_b2tautau_vis_JetBBEC1yearUp", &features_JetBBEC1yearUp.m_b2tautau_vis, "m_b2tautau_vis_JetBBEC1yearUp/F");
      TBranch *newBranch1_d_ma_JetBBEC1yearUp = t->Branch("d_ma_JetBBEC1yearUp", &features_JetBBEC1yearUp.d_ma, "d_ma_JetBBEC1yearUp/F");
      TBranch *newBranch1_mtMET_b_JetBBEC1yearDown = t->Branch("mtMET_b_JetBBEC1yearDown", &features_JetBBEC1yearDown.mtMET_b, "mtMET_b_JetBBEC1yearDown/F");
      TBranch *newBranch1_m_bb_JetBBEC1yearDown = t->Branch("m_bb_JetBBEC1yearDown", &features_JetBBEC1yearDown.m_bb, "m_bb_JetBBEC1yearDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetBBEC1yearDown = t->Branch("m_bbtautau_vis_JetBBEC1yearDown", &features_JetBBEC1yearDown.m_bbtautau_vis, "m_bbtautau_vis_JetBBEC1yearDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetBBEC1yearDown = t->Branch("m_b2tautau_vis_JetBBEC1yearDown", &features_JetBBEC1yearDown.m_b2tautau_vis, "m_b2tautau_vis_JetBBEC1yearDown/F");
      TBranch *newBranch1_d_ma_JetBBEC1yearDown = t->Branch("d_ma_JetBBEC1yearDown", &features_JetBBEC1yearDown.d_ma, "d_ma_JetBBEC1yearDown/F");
      // JetEC2
      TBranch *newBranch1_mtMET_b_JetEC2Up = t->Branch("mtMET_b_JetEC2Up", &features_JetEC2Up.mtMET_b, "mtMET_b_JetEC2Up/F");
      TBranch *newBranch1_m_bb_JetEC2Up = t->Branch("m_bb_JetEC2Up", &features_JetEC2Up.m_bb, "m_bb_JetEC2Up/F");
      TBranch *newBranch1_m_bbtautau_vis_JetEC2Up = t->Branch("m_bbtautau_vis_JetEC2Up", &features_JetEC2Up.m_bbtautau_vis, "m_bbtautau_vis_JetEC2Up/F");
      TBranch *newBranch1_m_b2tautau_vis_JetEC2Up = t->Branch("m_b2tautau_vis_JetEC2Up", &features_JetEC2Up.m_b2tautau_vis, "m_b2tautau_vis_JetEC2Up/F");
      TBranch *newBranch1_d_ma_JetEC2Up = t->Branch("d_ma_JetEC2Up", &features_JetEC2Up.d_ma, "d_ma_JetEC2Up/F");
      TBranch *newBranch1_mtMET_b_JetEC2Down = t->Branch("mtMET_b_JetEC2Down", &features_JetEC2Down.mtMET_b, "mtMET_b_JetEC2Down/F");
      TBranch *newBranch1_m_bb_JetEC2Down = t->Branch("m_bb_JetEC2Down", &features_JetEC2Down.m_bb, "m_bb_JetEC2Down/F");
      TBranch *newBranch1_m_bbtautau_vis_JetEC2Down = t->Branch("m_bbtautau_vis_JetEC2Down", &features_JetEC2Down.m_bbtautau_vis, "m_bbtautau_vis_JetEC2Down/F");
      TBranch *newBranch1_m_b2tautau_vis_JetEC2Down = t->Branch("m_b2tautau_vis_JetEC2Down", &features_JetEC2Down.m_b2tautau_vis, "m_b2tautau_vis_JetEC2Down/F");
      TBranch *newBranch1_d_ma_JetEC2Down = t->Branch("d_ma_JetEC2Down", &features_JetEC2Down.d_ma, "d_ma_JetEC2Down/F");
      // JetEC2year
      TBranch *newBranch1_mtMET_b_JetEC2yearUp = t->Branch("mtMET_b_JetEC2yearUp", &features_JetEC2yearUp.mtMET_b, "mtMET_b_JetEC2yearUp/F");
      TBranch *newBranch1_m_bb_JetEC2yearUp = t->Branch("m_bb_JetEC2yearUp", &features_JetEC2yearUp.m_bb, "m_bb_JetEC2yearUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetEC2yearUp = t->Branch("m_bbtautau_vis_JetEC2yearUp", &features_JetEC2yearUp.m_bbtautau_vis, "m_bbtautau_vis_JetEC2yearUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetEC2yearUp = t->Branch("m_b2tautau_vis_JetEC2yearUp", &features_JetEC2yearUp.m_b2tautau_vis, "m_b2tautau_vis_JetEC2yearUp/F");
      TBranch *newBranch1_d_ma_JetEC2yearUp = t->Branch("d_ma_JetEC2yearUp", &features_JetEC2yearUp.d_ma, "d_ma_JetEC2yearUp/F");
      TBranch *newBranch1_mtMET_b_JetEC2yearDown = t->Branch("mtMET_b_JetEC2yearDown", &features_JetEC2yearDown.mtMET_b, "mtMET_b_JetEC2yearDown/F");
      TBranch *newBranch1_m_bb_JetEC2yearDown = t->Branch("m_bb_JetEC2yearDown", &features_JetEC2yearDown.m_bb, "m_bb_JetEC2yearDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetEC2yearDown = t->Branch("m_bbtautau_vis_JetEC2yearDown", &features_JetEC2yearDown.m_bbtautau_vis, "m_bbtautau_vis_JetEC2yearDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetEC2yearDown = t->Branch("m_b2tautau_vis_JetEC2yearDown", &features_JetEC2yearDown.m_b2tautau_vis, "m_b2tautau_vis_JetEC2yearDown/F");
      TBranch *newBranch1_d_ma_JetEC2yearDown = t->Branch("d_ma_JetEC2yearDown", &features_JetEC2yearDown.d_ma, "d_ma_JetEC2yearDown/F");
      // JetFlavorQCD
      TBranch *newBranch1_mtMET_b_JetFlavorQCDUp = t->Branch("mtMET_b_JetFlavorQCDUp", &features_JetFlavorQCDUp.mtMET_b, "mtMET_b_JetFlavorQCDUp/F");
      TBranch *newBranch1_m_bb_JetFlavorQCDUp = t->Branch("m_bb_JetFlavorQCDUp", &features_JetFlavorQCDUp.m_bb, "m_bb_JetFlavorQCDUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetFlavorQCDUp = t->Branch("m_bbtautau_vis_JetFlavorQCDUp", &features_JetFlavorQCDUp.m_bbtautau_vis, "m_bbtautau_vis_JetFlavorQCDUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetFlavorQCDUp = t->Branch("m_b2tautau_vis_JetFlavorQCDUp", &features_JetFlavorQCDUp.m_b2tautau_vis, "m_b2tautau_vis_JetFlavorQCDUp/F");
      TBranch *newBranch1_d_ma_JetFlavorQCDUp = t->Branch("d_ma_JetFlavorQCDUp", &features_JetFlavorQCDUp.d_ma, "d_ma_JetFlavorQCDUp/F");
      TBranch *newBranch1_mtMET_b_JetFlavorQCDDown = t->Branch("mtMET_b_JetFlavorQCDDown", &features_JetFlavorQCDDown.mtMET_b, "mtMET_b_JetFlavorQCDDown/F");
      TBranch *newBranch1_m_bb_JetFlavorQCDDown = t->Branch("m_bb_JetFlavorQCDDown", &features_JetFlavorQCDDown.m_bb, "m_bb_JetFlavorQCDDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetFlavorQCDDown = t->Branch("m_bbtautau_vis_JetFlavorQCDDown", &features_JetFlavorQCDDown.m_bbtautau_vis, "m_bbtautau_vis_JetFlavorQCDDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetFlavorQCDDown = t->Branch("m_b2tautau_vis_JetFlavorQCDDown", &features_JetFlavorQCDDown.m_b2tautau_vis, "m_b2tautau_vis_JetFlavorQCDDown/F");
      TBranch *newBranch1_d_ma_JetFlavorQCDDown = t->Branch("d_ma_JetFlavorQCDDown", &features_JetFlavorQCDDown.d_ma, "d_ma_JetFlavorQCDDown/F");
      // JetHF
      TBranch *newBranch1_mtMET_b_JetHFUp = t->Branch("mtMET_b_JetHFUp", &features_JetHFUp.mtMET_b, "mtMET_b_JetHFUp/F");
      TBranch *newBranch1_m_bb_JetHFUp = t->Branch("m_bb_JetHFUp", &features_JetHFUp.m_bb, "m_bb_JetHFUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetHFUp = t->Branch("m_bbtautau_vis_JetHFUp", &features_JetHFUp.m_bbtautau_vis, "m_bbtautau_vis_JetHFUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetHFUp = t->Branch("m_b2tautau_vis_JetHFUp", &features_JetHFUp.m_b2tautau_vis, "m_b2tautau_vis_JetHFUp/F");
      TBranch *newBranch1_d_ma_JetHFUp = t->Branch("d_ma_JetHFUp", &features_JetHFUp.d_ma, "d_ma_JetHFUp/F");
      TBranch *newBranch1_mtMET_b_JetHFDown = t->Branch("mtMET_b_JetHFDown", &features_JetHFDown.mtMET_b, "mtMET_b_JetHFDown/F");
      TBranch *newBranch1_m_bb_JetHFDown = t->Branch("m_bb_JetHFDown", &features_JetHFDown.m_bb, "m_bb_JetHFDown/F"); 
      TBranch *newBranch1_m_bbtautau_vis_JetHFDown = t->Branch("m_bbtautau_vis_JetHFDown", &features_JetHFDown.m_bbtautau_vis, "m_bbtautau_vis_JetHFDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetHFDown = t->Branch("m_b2tautau_vis_JetHFDown", &features_JetHFDown.m_b2tautau_vis, "m_b2tautau_vis_JetHFDown/F");
      TBranch *newBranch1_d_ma_JetHFDown = t->Branch("d_ma_JetHFDown", &features_JetHFDown.d_ma, "d_ma_JetHFDown/F"); 
      // JetHFyear
      TBranch *newBranch1_mtMET_b_JetHFyearUp = t->Branch("mtMET_b_JetHFyearUp", &features_JetHFyearUp.mtMET_b, "mtMET_b_JetHFyearUp/F");
      TBranch *newBranch1_m_bb_JetHFyearUp = t->Branch("m_bb_JetHFyearUp", &features_JetHFyearUp.m_bb, "m_bb_JetHFyearUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetHFyearUp = t->Branch("m_bbtautau_vis_JetHFyearUp", &features_JetHFyearUp.m_bbtautau_vis, "m_bbtautau_vis_JetHFyearUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetHFyearUp = t->Branch("m_b2tautau_vis_JetHFyearUp", &features_JetHFyearUp.m_b2tautau_vis, "m_b2tautau_vis_JetHFyearUp/F");
      TBranch *newBranch1_d_ma_JetHFyearUp = t->Branch("d_ma_JetHFyearUp", &features_JetHFyearUp.d_ma, "d_ma_JetHFyearUp/F");
      TBranch *newBranch1_mtMET_b_JetHFyearDown = t->Branch("mtMET_b_JetHFyearDown", &features_JetHFyearDown.mtMET_b, "mtMET_b_JetHFyearDown/F");
      TBranch *newBranch1_m_bb_JetHFyearDown = t->Branch("m_bb_JetHFyearDown", &features_JetHFyearDown.m_bb, "m_bb_JetHFyearDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetHFyearDown = t->Branch("m_bbtautau_vis_JetHFyearDown", &features_JetHFyearDown.m_bbtautau_vis, "m_bbtautau_vis_JetHFyearDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetHFyearDown = t->Branch("m_b2tautau_vis_JetHFyearDown", &features_JetHFyearDown.m_b2tautau_vis, "m_b2tautau_vis_JetHFyearDown/F");
      TBranch *newBranch1_d_ma_JetHFyearDown = t->Branch("d_ma_JetHFyearDown", &features_JetHFyearDown.d_ma, "d_ma_JetHFyearDown/F");
      // JetRelativeBal
      TBranch *newBranch1_mtMET_b_JetRelativeBalUp = t->Branch("mtMET_b_JetRelativeBalUp", &features_JetRelativeBalUp.mtMET_b, "mtMET_b_JetRelativeBalUp/F");
      TBranch *newBranch1_m_bb_JetRelativeBalUp = t->Branch("m_bb_JetRelativeBalUp", &features_JetRelativeBalUp.m_bb, "m_bb_JetRelativeBalUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetRelativeBalUp = t->Branch("m_bbtautau_vis_JetRelativeBalUp", &features_JetRelativeBalUp.m_bbtautau_vis, "m_bbtautau_vis_JetRelativeBalUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetRelativeBalUp = t->Branch("m_b2tautau_vis_JetRelativeBalUp", &features_JetRelativeBalUp.m_b2tautau_vis, "m_b2tautau_vis_JetRelativeBalUp/F");
      TBranch *newBranch1_d_ma_JetRelativeBalUp = t->Branch("d_ma_JetRelativeBalUp", &features_JetRelativeBalUp.d_ma, "d_ma_JetRelativeBalUp/F");
      TBranch *newBranch1_mtMET_b_JetRelativeBalDown = t->Branch("mtMET_b_JetRelativeBalDown", &features_JetRelativeBalDown.mtMET_b, "mtMET_b_JetRelativeBalDown/F");
      TBranch *newBranch1_m_bb_JetRelativeBalDown = t->Branch("m_bb_JetRelativeBalDown", &features_JetRelativeBalDown.m_bb, "m_bb_JetRelativeBalDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetRelativeBalDown = t->Branch("m_bbtautau_vis_JetRelativeBalDown", &features_JetRelativeBalDown.m_bbtautau_vis, "m_bbtautau_vis_JetRelativeBalDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetRelativeBalDown = t->Branch("m_b2tautau_vis_JetRelativeBalDown", &features_JetRelativeBalDown.m_b2tautau_vis, "m_b2tautau_vis_JetRelativeBalDown/F");
      TBranch *newBranch1_d_ma_JetRelativeBalDown = t->Branch("d_ma_JetRelativeBalDown", &features_JetRelativeBalDown.d_ma, "d_ma_JetRelativeBalDown/F");
      // JetRelativeSample
      TBranch *newBranch1_mtMET_b_JetRelativeSampleUp = t->Branch("mtMET_b_JetRelativeSampleUp", &features_JetRelativeSampleUp.mtMET_b, "mtMET_b_JetRelativeSampleUp/F");
      TBranch *newBranch1_m_bb_JetRelativeSampleUp = t->Branch("m_bb_JetRelativeSampleUp", &features_JetRelativeSampleUp.m_bb, "m_bb_JetRelativeSampleUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JetRelativeSampleUp = t->Branch("m_bbtautau_vis_JetRelativeSampleUp", &features_JetRelativeSampleUp.m_bbtautau_vis, "m_bbtautau_vis_JetRelativeSampleUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JetRelativeSampleUp = t->Branch("m_b2tautau_vis_JetRelativeSampleUp", &features_JetRelativeSampleUp.m_b2tautau_vis, "m_b2tautau_vis_JetRelativeSampleUp/F");
      TBranch *newBranch1_d_ma_JetRelativeSampleUp = t->Branch("d_ma_JetRelativeSampleUp", &features_JetRelativeSampleUp.d_ma, "d_ma_JetRelativeSampleUp/F");
      TBranch *newBranch1_mtMET_b_JetRelativeSampleDown = t->Branch("mtMET_b_JetRelativeSampleDown", &features_JetRelativeSampleDown.mtMET_b, "mtMET_b_JetRelativeSampleDown/F");
      TBranch *newBranch1_m_bb_JetRelativeSampleDown = t->Branch("m_bb_JetRelativeSampleDown", &features_JetRelativeSampleDown.m_bb, "m_bb_JetRelativeSampleDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JetRelativeSampleDown = t->Branch("m_bbtautau_vis_JetRelativeSampleDown", &features_JetRelativeSampleDown.m_bbtautau_vis, "m_bbtautau_vis_JetRelativeSampleDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JetRelativeSampleDown = t->Branch("m_b2tautau_vis_JetRelativeSampleDown", &features_JetRelativeSampleDown.m_b2tautau_vis, "m_b2tautau_vis_JetRelativeSampleDown/F");
      TBranch *newBranch1_d_ma_JetRelativeSampleDown = t->Branch("d_ma_JetRelativeSampleDown", &features_JetRelativeSampleDown.d_ma, "d_ma_JetRelativeSampleDown/F");
      // JER
      TBranch *newBranch1_mtMET_b_JERUp = t->Branch("mtMET_b_JERUp", &features_JERUp.mtMET_b, "mtMET_b_JERUp/F");
      TBranch *newBranch1_m_bb_JERUp = t->Branch("m_bb_JERUp", &features_JERUp.m_bb, "m_bb_JERUp/F");
      TBranch *newBranch1_m_bbtautau_vis_JERUp = t->Branch("m_bbtautau_vis_JERUp", &features_JERUp.m_bbtautau_vis, "m_bbtautau_vis_JERUp/F");
      TBranch *newBranch1_m_b2tautau_vis_JERUp = t->Branch("m_b2tautau_vis_JERUp", &features_JERUp.m_b2tautau_vis, "m_b2tautau_vis_JERUp/F");
      TBranch *newBranch1_d_ma_JERUp = t->Branch("d_ma_JERUp", &features_JERUp.d_ma, "d_ma_JERUp/F");
      TBranch *newBranch1_mtMET_b_JERDown = t->Branch("mtMET_b_JERDown", &features_JERDown.mtMET_b, "mtMET_b_JERDown/F");
      TBranch *newBranch1_m_bb_JERDown = t->Branch("m_bb_JERDown", &features_JERDown.m_bb, "m_bb_JERDown/F");
      TBranch *newBranch1_m_bbtautau_vis_JERDown = t->Branch("m_bbtautau_vis_JERDown", &features_JERDown.m_bbtautau_vis, "m_bbtautau_vis_JERDown/F");
      TBranch *newBranch1_m_b2tautau_vis_JERDown = t->Branch("m_b2tautau_vis_JERDown", &features_JERDown.m_b2tautau_vis, "m_b2tautau_vis_JERDown/F");
      TBranch *newBranch1_d_ma_JERDown = t->Branch("d_ma_JERDown", &features_JERDown.d_ma, "d_ma_JERDown/F");

      // read the branches needed for input to BDT
      Int_t channel;
      Float_t pt_1_nominal;
      Float_t pt_1_es1Up;
      Float_t pt_1_es1Down;
      Float_t eta_1;
      Float_t phi_1;
      Float_t m_1_nominal;
      Float_t m_1_es1Up;
      Float_t m_1_es1Down;
      Float_t pt_2_nominal;
      Float_t pt_2_es2Up;
      Float_t pt_2_es2Down;
      Float_t eta_2;
      Float_t phi_2;
      Float_t m_2_nominal;
      Float_t m_2_es2Up;
      Float_t m_2_es2Down;
      Float_t bpt_deepflavour_1;
      Float_t bpt_deepflavour_JERUp_1;
      Float_t bpt_deepflavour_JERDown_1;
      Float_t bpt_deepflavour_JetAbsoluteUp_1;
      Float_t bpt_deepflavour_JetAbsoluteDown_1;
      Float_t bpt_deepflavour_JetAbsoluteyearUp_1;
      Float_t bpt_deepflavour_JetAbsoluteyearDown_1;
      Float_t bpt_deepflavour_JetBBEC1Up_1;
      Float_t bpt_deepflavour_JetBBEC1Down_1;
      Float_t bpt_deepflavour_JetBBEC1yearUp_1;
      Float_t bpt_deepflavour_JetBBEC1yearDown_1;
      Float_t bpt_deepflavour_JetEC2Up_1;
      Float_t bpt_deepflavour_JetEC2Down_1;
      Float_t bpt_deepflavour_JetEC2yearUp_1;
      Float_t bpt_deepflavour_JetEC2yearDown_1;
      Float_t bpt_deepflavour_JetFlavorQCDUp_1;
      Float_t bpt_deepflavour_JetFlavorQCDDown_1;
      Float_t bpt_deepflavour_JetHFUp_1;
      Float_t bpt_deepflavour_JetHFDown_1;
      Float_t bpt_deepflavour_JetHFyearUp_1;
      Float_t bpt_deepflavour_JetHFyearDown_1;
      Float_t bpt_deepflavour_JetRelativeBalUp_1;
      Float_t bpt_deepflavour_JetRelativeBalDown_1;
      Float_t bpt_deepflavour_JetRelativeSampleUp_1;
      Float_t bpt_deepflavour_JetRelativeSampleDown_1;
      Float_t beta_deepflavour_1;
      Float_t bphi_deepflavour_1;
      Float_t bm_deepflavour_1;
      Float_t bm_deepflavour_JERUp_1;
      Float_t bm_deepflavour_JERDown_1;
      Float_t bm_deepflavour_JetAbsoluteUp_1;
      Float_t bm_deepflavour_JetAbsoluteDown_1;
      Float_t bm_deepflavour_JetAbsoluteyearUp_1;
      Float_t bm_deepflavour_JetAbsoluteyearDown_1;
      Float_t bm_deepflavour_JetBBEC1Up_1;
      Float_t bm_deepflavour_JetBBEC1Down_1;
      Float_t bm_deepflavour_JetBBEC1yearUp_1;
      Float_t bm_deepflavour_JetBBEC1yearDown_1;
      Float_t bm_deepflavour_JetEC2Up_1;
      Float_t bm_deepflavour_JetEC2Down_1;
      Float_t bm_deepflavour_JetEC2yearUp_1;
      Float_t bm_deepflavour_JetEC2yearDown_1;
      Float_t bm_deepflavour_JetFlavorQCDUp_1;
      Float_t bm_deepflavour_JetFlavorQCDDown_1;
      Float_t bm_deepflavour_JetHFUp_1;
      Float_t bm_deepflavour_JetHFDown_1;
      Float_t bm_deepflavour_JetHFyearUp_1;
      Float_t bm_deepflavour_JetHFyearDown_1;
      Float_t bm_deepflavour_JetRelativeBalUp_1;
      Float_t bm_deepflavour_JetRelativeBalDown_1;
      Float_t bm_deepflavour_JetRelativeSampleUp_1;
      Float_t bm_deepflavour_JetRelativeSampleDown_1;
      Float_t met_nominal;
      Float_t met_es1Up;
      Float_t met_es1Down;
      Float_t met_es2Up;
      Float_t met_es2Down;
      Float_t met_UESUp;
      Float_t met_UESDown;
      Float_t met_responseUp;
      Float_t met_responseDown;
      Float_t met_resolutionUp;
      Float_t met_resolutionDown;
      Float_t met_JERUp;
      Float_t met_JERDown;
      Float_t met_JetAbsoluteUp;
      Float_t met_JetAbsoluteDown;
      Float_t met_JetAbsoluteyearUp;
      Float_t met_JetAbsoluteyearDown;
      Float_t met_JetBBEC1Up;
      Float_t met_JetBBEC1Down;
      Float_t met_JetBBEC1yearUp;
      Float_t met_JetBBEC1yearDown;
      Float_t met_JetEC2Up;
      Float_t met_JetEC2Down;
      Float_t met_JetEC2yearUp;
      Float_t met_JetEC2yearDown;
      Float_t met_JetFlavorQCDUp;
      Float_t met_JetFlavorQCDDown;
      Float_t met_JetHFUp;
      Float_t met_JetHFDown;
      Float_t met_JetHFyearUp;
      Float_t met_JetHFyearDown;
      Float_t met_JetRelativeBalUp;
      Float_t met_JetRelativeBalDown;
      Float_t met_JetRelativeSampleUp;
      Float_t met_JetRelativeSampleDown;
      Float_t metphi_nominal;
      Float_t metphi_es1Up;
      Float_t metphi_es1Down;
      Float_t metphi_es2Up;
      Float_t metphi_es2Down;
      Float_t metphi_UESUp;
      Float_t metphi_UESDown;
      Float_t metphi_responseUp;
      Float_t metphi_responseDown;
      Float_t metphi_resolutionUp;
      Float_t metphi_resolutionDown;
      Float_t metphi_JERUp;
      Float_t metphi_JERDown;
      Float_t metphi_JetAbsoluteUp;
      Float_t metphi_JetAbsoluteDown;
      Float_t metphi_JetAbsoluteyearUp;
      Float_t metphi_JetAbsoluteyearDown;
      Float_t metphi_JetBBEC1Up;
      Float_t metphi_JetBBEC1Down;
      Float_t metphi_JetBBEC1yearUp;
      Float_t metphi_JetBBEC1yearDown;
      Float_t metphi_JetEC2Up;
      Float_t metphi_JetEC2Down;
      Float_t metphi_JetEC2yearUp;
      Float_t metphi_JetEC2yearDown;
      Float_t metphi_JetFlavorQCDUp;
      Float_t metphi_JetFlavorQCDDown;
      Float_t metphi_JetHFUp;
      Float_t metphi_JetHFDown;
      Float_t metphi_JetHFyearUp;
      Float_t metphi_JetHFyearDown;
      Float_t metphi_JetRelativeBalUp;
      Float_t metphi_JetRelativeBalDown;
      Float_t metphi_JetRelativeSampleUp;
      Float_t metphi_JetRelativeSampleDown;
      Float_t D_zeta_nominal;
      Float_t D_zeta_es1Up;
      Float_t D_zeta_es1Down;
      Float_t D_zeta_es2Up;
      Float_t D_zeta_es2Down;
      Float_t D_zeta_UESUp;
      Float_t D_zeta_UESDown;
      Float_t D_zeta_responseUp;
      Float_t D_zeta_responseDown;
      Float_t D_zeta_resolutionUp;
      Float_t D_zeta_resolutionDown;
      Float_t D_zeta_JERUp;
      Float_t D_zeta_JERDown;
      Float_t D_zeta_JetAbsoluteUp;
      Float_t D_zeta_JetAbsoluteDown;
      Float_t D_zeta_JetAbsoluteyearUp;
      Float_t D_zeta_JetAbsoluteyearDown;
      Float_t D_zeta_JetBBEC1Up;
      Float_t D_zeta_JetBBEC1Down;
      Float_t D_zeta_JetBBEC1yearUp;
      Float_t D_zeta_JetBBEC1yearDown;
      Float_t D_zeta_JetEC2Up;
      Float_t D_zeta_JetEC2Down;
      Float_t D_zeta_JetEC2yearUp;
      Float_t D_zeta_JetEC2yearDown;
      Float_t D_zeta_JetFlavorQCDUp;
      Float_t D_zeta_JetFlavorQCDDown;
      Float_t D_zeta_JetHFUp;
      Float_t D_zeta_JetHFDown;
      Float_t D_zeta_JetHFyearUp;
      Float_t D_zeta_JetHFyearDown;
      Float_t D_zeta_JetRelativeBalUp;
      Float_t D_zeta_JetRelativeBalDown;
      Float_t D_zeta_JetRelativeSampleUp;
      Float_t D_zeta_JetRelativeSampleDown;

      Float_t bpt_deepflavour_2;
      Float_t bpt_deepflavour_JERUp_2;
      Float_t bpt_deepflavour_JERDown_2;
      Float_t bpt_deepflavour_JetAbsoluteUp_2;
      Float_t bpt_deepflavour_JetAbsoluteDown_2;
      Float_t bpt_deepflavour_JetAbsoluteyearUp_2;
      Float_t bpt_deepflavour_JetAbsoluteyearDown_2;
      Float_t bpt_deepflavour_JetBBEC1Up_2;
      Float_t bpt_deepflavour_JetBBEC1Down_2;
      Float_t bpt_deepflavour_JetBBEC1yearUp_2;
      Float_t bpt_deepflavour_JetBBEC1yearDown_2;
      Float_t bpt_deepflavour_JetEC2Up_2;
      Float_t bpt_deepflavour_JetEC2Down_2;
      Float_t bpt_deepflavour_JetEC2yearUp_2;
      Float_t bpt_deepflavour_JetEC2yearDown_2;
      Float_t bpt_deepflavour_JetFlavorQCDUp_2;
      Float_t bpt_deepflavour_JetFlavorQCDDown_2;
      Float_t bpt_deepflavour_JetHFUp_2;
      Float_t bpt_deepflavour_JetHFDown_2;
      Float_t bpt_deepflavour_JetHFyearUp_2;
      Float_t bpt_deepflavour_JetHFyearDown_2;
      Float_t bpt_deepflavour_JetRelativeBalUp_2;
      Float_t bpt_deepflavour_JetRelativeBalDown_2;
      Float_t bpt_deepflavour_JetRelativeSampleUp_2;
      Float_t bpt_deepflavour_JetRelativeSampleDown_2;
      Float_t beta_deepflavour_2;
      Float_t bphi_deepflavour_2;
      Float_t bm_deepflavour_2;
      Float_t bm_deepflavour_JERUp_2;
      Float_t bm_deepflavour_JERDown_2;
      Float_t bm_deepflavour_JetAbsoluteUp_2;
      Float_t bm_deepflavour_JetAbsoluteDown_2;
      Float_t bm_deepflavour_JetAbsoluteyearUp_2;
      Float_t bm_deepflavour_JetAbsoluteyearDown_2;
      Float_t bm_deepflavour_JetBBEC1Up_2;
      Float_t bm_deepflavour_JetBBEC1Down_2;
      Float_t bm_deepflavour_JetBBEC1yearUp_2;
      Float_t bm_deepflavour_JetBBEC1yearDown_2;
      Float_t bm_deepflavour_JetEC2Up_2;
      Float_t bm_deepflavour_JetEC2Down_2;
      Float_t bm_deepflavour_JetEC2yearUp_2;
      Float_t bm_deepflavour_JetEC2yearDown_2;
      Float_t bm_deepflavour_JetFlavorQCDUp_2;
      Float_t bm_deepflavour_JetFlavorQCDDown_2;
      Float_t bm_deepflavour_JetHFUp_2;
      Float_t bm_deepflavour_JetHFDown_2;
      Float_t bm_deepflavour_JetHFyearUp_2;
      Float_t bm_deepflavour_JetHFyearDown_2;
      Float_t bm_deepflavour_JetRelativeBalUp_2;
      Float_t bm_deepflavour_JetRelativeBalDown_2;
      Float_t bm_deepflavour_JetRelativeSampleUp_2;
      Float_t bm_deepflavour_JetRelativeSampleDown_2;
      
      Int_t njets;
      Int_t njets_JERUp;
      Int_t njets_JERDown;
      Int_t njets_JetAbsoluteUp;
      Int_t njets_JetAbsoluteDown;
      Int_t njets_JetAbsoluteyearUp;
      Int_t njets_JetAbsoluteyearDown;
      Int_t njets_JetBBEC1Up;
      Int_t njets_JetBBEC1Down;
      Int_t njets_JetBBEC1yearUp;
      Int_t njets_JetBBEC1yearDown;
      Int_t njets_JetEC2Up;
      Int_t njets_JetEC2Down;
      Int_t njets_JetEC2yearUp;
      Int_t njets_JetEC2yearDown;
      Int_t njets_JetFlavorQCDUp;
      Int_t njets_JetFlavorQCDDown;
      Int_t njets_JetHFUp;
      Int_t njets_JetHFDown;
      Int_t njets_JetHFyearUp;
      Int_t njets_JetHFyearDown;
      Int_t njets_JetRelativeBalUp;
      Int_t njets_JetRelativeBalDown;
      Int_t njets_JetRelativeSampleUp;
      Int_t njets_JetRelativeSampleDown;
      Float_t pt_vis_nominal;
      Float_t pt_vis_es1Up;
      Float_t pt_vis_es1Down;
      Float_t pt_vis_es2Up;
      Float_t pt_vis_es2Down;
      Float_t m_btautau_vis_nominal;
      Float_t m_btautau_vis_es1Up;
      Float_t m_btautau_vis_es1Down;
      Float_t m_btautau_vis_es2Up;
      Float_t m_btautau_vis_es2Down;
      Float_t m_btautau_vis_JetAbsoluteUp_1;
      Float_t m_btautau_vis_JetAbsoluteDown_1;
      Float_t m_btautau_vis_JetAbsoluteyearUp_1;
      Float_t m_btautau_vis_JetAbsoluteyearDown_1;
      Float_t m_btautau_vis_JetBBEC1Up_1;
      Float_t m_btautau_vis_JetBBEC1Down_1;
      Float_t m_btautau_vis_JetBBEC1yearUp_1;
      Float_t m_btautau_vis_JetBBEC1yearDown_1;
      Float_t m_btautau_vis_JetEC2Up_1;
      Float_t m_btautau_vis_JetEC2Down_1;
      Float_t m_btautau_vis_JetEC2yearUp_1;
      Float_t m_btautau_vis_JetEC2yearDown_1;
      Float_t m_btautau_vis_JetFlavorQCDUp_1;
      Float_t m_btautau_vis_JetFlavorQCDDown_1;
      Float_t m_btautau_vis_JetHFUp_1;
      Float_t m_btautau_vis_JetHFDown_1;
      Float_t m_btautau_vis_JetHFyearUp_1;
      Float_t m_btautau_vis_JetHFyearDown_1;
      Float_t m_btautau_vis_JetRelativeBalUp_1;
      Float_t m_btautau_vis_JetRelativeBalDown_1;
      Float_t m_btautau_vis_JetRelativeSampleUp_1;
      Float_t m_btautau_vis_JetRelativeSampleDown_1;
      Float_t m_btautau_vis_JERUp_1;
      Float_t m_btautau_vis_JERDown_1;
      Float_t mtMET_1_nominal;
      Float_t mtMET_1_es1Up;
      Float_t mtMET_1_es1Down;
      Float_t mtMET_1_responseUp;
      Float_t mtMET_1_responseDown;
      Float_t mtMET_1_resolutionUp;
      Float_t mtMET_1_resolutionDown;
      Float_t mtMET_1_JetAbsoluteUp;
      Float_t mtMET_1_JetAbsoluteDown;
      Float_t mtMET_1_JetAbsoluteyearUp;
      Float_t mtMET_1_JetAbsoluteyearDown;
      Float_t mtMET_1_JetBBEC1Up;
      Float_t mtMET_1_JetBBEC1Down;
      Float_t mtMET_1_JetBBEC1yearUp;
      Float_t mtMET_1_JetBBEC1yearDown;
      Float_t mtMET_1_JetEC2Up;
      Float_t mtMET_1_JetEC2Down;
      Float_t mtMET_1_JetEC2yearUp;
      Float_t mtMET_1_JetEC2yearDown;
      Float_t mtMET_1_JetFlavorQCDUp;
      Float_t mtMET_1_JetFlavorQCDDown;
      Float_t mtMET_1_JetHFUp;
      Float_t mtMET_1_JetHFDown;
      Float_t mtMET_1_JetHFyearUp;
      Float_t mtMET_1_JetHFyearDown;
      Float_t mtMET_1_JetRelativeBalUp;
      Float_t mtMET_1_JetRelativeBalDown;
      Float_t mtMET_1_JetRelativeSampleUp;
      Float_t mtMET_1_JetRelativeSampleDown;
      Float_t mtMET_1_JERUp;
      Float_t mtMET_1_JERDown;
      Float_t mtMET_2_nominal;
      Float_t mtMET_2_es2Up;
      Float_t mtMET_2_es2Down;
      Float_t mtMET_2_responseUp;
      Float_t mtMET_2_responseDown;
      Float_t mtMET_2_resolutionUp;
      Float_t mtMET_2_resolutionDown;
      Float_t mtMET_2_JetAbsoluteUp;
      Float_t mtMET_2_JetAbsoluteDown;
      Float_t mtMET_2_JetAbsoluteyearUp;
      Float_t mtMET_2_JetAbsoluteyearDown;
      Float_t mtMET_2_JetBBEC1Up;
      Float_t mtMET_2_JetBBEC1Down;
      Float_t mtMET_2_JetBBEC1yearUp;
      Float_t mtMET_2_JetBBEC1yearDown;
      Float_t mtMET_2_JetEC2Up;
      Float_t mtMET_2_JetEC2Down;
      Float_t mtMET_2_JetEC2yearUp;
      Float_t mtMET_2_JetEC2yearDown;
      Float_t mtMET_2_JetFlavorQCDUp;
      Float_t mtMET_2_JetFlavorQCDDown;
      Float_t mtMET_2_JetHFUp;
      Float_t mtMET_2_JetHFDown;
      Float_t mtMET_2_JetHFyearUp;
      Float_t mtMET_2_JetHFyearDown;
      Float_t mtMET_2_JetRelativeBalUp;
      Float_t mtMET_2_JetRelativeBalDown;
      Float_t mtMET_2_JetRelativeSampleUp;
      Float_t mtMET_2_JetRelativeSampleDown;
      Float_t mtMET_2_JERUp;
      Float_t mtMET_2_JERDown;

      t->SetBranchAddress("channel", &channel);
      t->SetBranchAddress("pt_1_nominal", &pt_1_nominal);
      t->SetBranchAddress("pt_1_es1Up",&pt_1_es1Up);
      t->SetBranchAddress("pt_1_es1Down",&pt_1_es1Down);
      t->SetBranchAddress("eta_1", &eta_1);
      t->SetBranchAddress("phi_1", &phi_1);
      t->SetBranchAddress("m_1_nominal", &m_1_nominal);
      t->SetBranchAddress("m_1_es1Up", &m_1_es1Up);
      t->SetBranchAddress("m_1_es1Down", &m_1_es1Down);
      t->SetBranchAddress("pt_2_nominal", &pt_2_nominal);
      t->SetBranchAddress("pt_2_es2Up", &pt_2_es2Up);
      t->SetBranchAddress("pt_2_es2Down",&pt_2_es2Down);
      t->SetBranchAddress("eta_2", &eta_2);
      t->SetBranchAddress("phi_2", &phi_2);
      t->SetBranchAddress("m_2_nominal", &m_2_nominal);
      t->SetBranchAddress("m_2_es2Up",  &m_2_es2Up);
      t->SetBranchAddress("m_2_es2Down",&m_2_es2Down);
      t->SetBranchAddress("bpt_deepflavour_1", &bpt_deepflavour_1);
      t->SetBranchAddress("bpt_deepflavour_JERUp_1", &bpt_deepflavour_JERUp_1);
      t->SetBranchAddress("bpt_deepflavour_JERDown_1", &bpt_deepflavour_JERDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteUp_1", &bpt_deepflavour_JetAbsoluteUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteDown_1", &bpt_deepflavour_JetAbsoluteDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearUp_1", &bpt_deepflavour_JetAbsoluteyearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearDown_1", &bpt_deepflavour_JetAbsoluteyearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Up_1", &bpt_deepflavour_JetBBEC1Up_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Down_1", &bpt_deepflavour_JetBBEC1Down_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearUp_1", &bpt_deepflavour_JetBBEC1yearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearDown_1", &bpt_deepflavour_JetBBEC1yearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Up_1", &bpt_deepflavour_JetEC2Up_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Down_1", &bpt_deepflavour_JetEC2Down_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearUp_1", &bpt_deepflavour_JetEC2yearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearDown_1", &bpt_deepflavour_JetEC2yearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDUp_1", &bpt_deepflavour_JetFlavorQCDUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDDown_1", &bpt_deepflavour_JetFlavorQCDDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFUp_1", &bpt_deepflavour_JetHFUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFDown_1", &bpt_deepflavour_JetHFDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearUp_1", &bpt_deepflavour_JetHFyearUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearDown_1", &bpt_deepflavour_JetHFyearDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalUp_1", &bpt_deepflavour_JetRelativeBalUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalDown_1", &bpt_deepflavour_JetRelativeBalDown_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleUp_1", &bpt_deepflavour_JetRelativeSampleUp_1);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleDown_1", &bpt_deepflavour_JetRelativeSampleDown_1);
      t->SetBranchAddress("beta_deepflavour_1", &beta_deepflavour_1);
      t->SetBranchAddress("bphi_deepflavour_1", &bphi_deepflavour_1);
      t->SetBranchAddress("bm_deepflavour_1", &bm_deepflavour_1);
      t->SetBranchAddress("bm_deepflavour_JERUp_1", &bm_deepflavour_JERUp_1);
      t->SetBranchAddress("bm_deepflavour_JERDown_1", &bm_deepflavour_JERDown_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteUp_1", &bm_deepflavour_JetAbsoluteUp_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteDown_1", &bm_deepflavour_JetAbsoluteDown_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearUp_1", &bm_deepflavour_JetAbsoluteyearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearDown_1", &bm_deepflavour_JetAbsoluteyearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Up_1", &bm_deepflavour_JetBBEC1Up_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Down_1", &bm_deepflavour_JetBBEC1Down_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearUp_1", &bm_deepflavour_JetBBEC1yearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearDown_1", &bm_deepflavour_JetBBEC1yearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2Up_1", &bm_deepflavour_JetEC2Up_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2Down_1", &bm_deepflavour_JetEC2Down_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearUp_1", &bm_deepflavour_JetEC2yearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearDown_1", &bm_deepflavour_JetEC2yearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDUp_1", &bm_deepflavour_JetFlavorQCDUp_1);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDDown_1", &bm_deepflavour_JetFlavorQCDDown_1);
      t->SetBranchAddress("bm_deepflavour_JetHFUp_1", &bm_deepflavour_JetHFUp_1);
      t->SetBranchAddress("bm_deepflavour_JetHFDown_1", &bm_deepflavour_JetHFDown_1);
      t->SetBranchAddress("bm_deepflavour_JetHFyearUp_1", &bm_deepflavour_JetHFyearUp_1);
      t->SetBranchAddress("bm_deepflavour_JetHFyearDown_1", &bm_deepflavour_JetHFyearDown_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalUp_1", &bm_deepflavour_JetRelativeBalUp_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalDown_1", &bm_deepflavour_JetRelativeBalDown_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleUp_1", &bm_deepflavour_JetRelativeSampleUp_1);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleDown_1", &bm_deepflavour_JetRelativeSampleDown_1);
      t->SetBranchAddress("met_nominal", &met_nominal);
      t->SetBranchAddress("met_es1Up", &met_es1Up);
      t->SetBranchAddress("met_es1Down", &met_es1Down);
      t->SetBranchAddress("met_es2Up", &met_es2Up);
      t->SetBranchAddress("met_es2Down", &met_es2Down);
      t->SetBranchAddress("met_UESUp", &met_UESUp);
      t->SetBranchAddress("met_UESDown", &met_UESDown);
      t->SetBranchAddress("met_responseUp", &met_responseUp);
      t->SetBranchAddress("met_responseDown", &met_responseDown);
      t->SetBranchAddress("met_resolutionUp", &met_resolutionUp);
      t->SetBranchAddress("met_resolutionDown", &met_resolutionDown);
      t->SetBranchAddress("met_JetAbsoluteUp", &met_JetAbsoluteUp);
      t->SetBranchAddress("met_JetAbsoluteDown", &met_JetAbsoluteDown);
      t->SetBranchAddress("met_JetAbsoluteyearUp", &met_JetAbsoluteyearUp);
      t->SetBranchAddress("met_JetAbsoluteyearDown", &met_JetAbsoluteyearDown);
      t->SetBranchAddress("met_JetBBEC1Up", &met_JetBBEC1Up);
      t->SetBranchAddress("met_JetBBEC1Down", &met_JetBBEC1Down);
      t->SetBranchAddress("met_JetBBEC1yearUp", &met_JetBBEC1yearUp);
      t->SetBranchAddress("met_JetBBEC1yearDown", &met_JetBBEC1yearDown);
      t->SetBranchAddress("met_JetEC2Up", &met_JetEC2Up);
      t->SetBranchAddress("met_JetEC2Down", &met_JetEC2Down);
      t->SetBranchAddress("met_JetEC2yearUp", &met_JetEC2yearUp);
      t->SetBranchAddress("met_JetEC2yearDown", &met_JetEC2yearDown);
      t->SetBranchAddress("met_JetFlavorQCDUp", &met_JetFlavorQCDUp);
      t->SetBranchAddress("met_JetFlavorQCDDown", &met_JetFlavorQCDDown);
      t->SetBranchAddress("met_JetHFUp", &met_JetHFUp);
      t->SetBranchAddress("met_JetHFDown", &met_JetHFDown);
      t->SetBranchAddress("met_JetHFyearUp", &met_JetHFyearUp);
      t->SetBranchAddress("met_JetHFyearDown", &met_JetHFyearDown);
      t->SetBranchAddress("met_JetRelativeBalUp", &met_JetRelativeBalUp);
      t->SetBranchAddress("met_JetRelativeBalDown", &met_JetRelativeBalDown);
      t->SetBranchAddress("met_JetRelativeSampleUp", &met_JetRelativeSampleUp);
      t->SetBranchAddress("met_JetRelativeSampleDown", &met_JetRelativeSampleDown);
      t->SetBranchAddress("met_JERUp", &met_JERUp);
      t->SetBranchAddress("met_JERDown", &met_JERDown);
      t->SetBranchAddress("metphi_nominal", &metphi_nominal);
      t->SetBranchAddress("metphi_es1Up", &metphi_es1Up);
      t->SetBranchAddress("metphi_es1Down",&metphi_es1Down);
      t->SetBranchAddress("metphi_es2Up", &metphi_es2Up);
      t->SetBranchAddress("metphi_es2Down", &metphi_es2Down);
      t->SetBranchAddress("metphi_UESUp", &metphi_UESUp);
      t->SetBranchAddress("metphi_UESDown", &metphi_UESDown);
      t->SetBranchAddress("metphi_responseUp", &metphi_responseUp);
      t->SetBranchAddress("metphi_responseDown", &metphi_responseDown);
      t->SetBranchAddress("metphi_resolutionUp", &metphi_resolutionUp);
      t->SetBranchAddress("metphi_resolutionDown", &metphi_resolutionDown);
      t->SetBranchAddress("metphi_JetAbsoluteUp", &metphi_JetAbsoluteUp);
      t->SetBranchAddress("metphi_JetAbsoluteDown", &metphi_JetAbsoluteDown);
      t->SetBranchAddress("metphi_JetAbsoluteyearUp", &metphi_JetAbsoluteyearUp);
      t->SetBranchAddress("metphi_JetAbsoluteyearDown", &metphi_JetAbsoluteyearDown);
      t->SetBranchAddress("metphi_JetBBEC1Up", &metphi_JetBBEC1Up);
      t->SetBranchAddress("metphi_JetBBEC1Down", &metphi_JetBBEC1Down);
      t->SetBranchAddress("metphi_JetBBEC1yearUp", &metphi_JetBBEC1yearUp);
      t->SetBranchAddress("metphi_JetBBEC1yearDown", &metphi_JetBBEC1yearDown);
      t->SetBranchAddress("metphi_JetEC2Up", &metphi_JetEC2Up);
      t->SetBranchAddress("metphi_JetEC2Down", &metphi_JetEC2Down);
      t->SetBranchAddress("metphi_JetEC2yearUp", &metphi_JetEC2yearUp);
      t->SetBranchAddress("metphi_JetEC2yearDown", &metphi_JetEC2yearDown);
      t->SetBranchAddress("metphi_JetFlavorQCDUp", &metphi_JetFlavorQCDUp);
      t->SetBranchAddress("metphi_JetFlavorQCDDown", &metphi_JetFlavorQCDDown);
      t->SetBranchAddress("metphi_JetHFUp", &metphi_JetHFUp);
      t->SetBranchAddress("metphi_JetHFDown", &metphi_JetHFDown);
      t->SetBranchAddress("metphi_JetHFyearUp", &metphi_JetHFyearUp);
      t->SetBranchAddress("metphi_JetHFyearDown", &metphi_JetHFyearDown);
      t->SetBranchAddress("metphi_JetRelativeBalUp", &metphi_JetRelativeBalUp);
      t->SetBranchAddress("metphi_JetRelativeBalDown", &metphi_JetRelativeBalDown);
      t->SetBranchAddress("metphi_JetRelativeSampleUp", &metphi_JetRelativeSampleUp);
      t->SetBranchAddress("metphi_JetRelativeSampleDown", &metphi_JetRelativeSampleDown);
      t->SetBranchAddress("metphi_JERUp", &metphi_JERUp);
      t->SetBranchAddress("metphi_JERDown", &metphi_JERDown);
      t->SetBranchAddress("D_zeta_nominal", &D_zeta_nominal);
      t->SetBranchAddress("D_zeta_es1Up", &D_zeta_es1Up);
      t->SetBranchAddress("D_zeta_es1Down", &D_zeta_es1Down);
      t->SetBranchAddress("D_zeta_es2Up", &D_zeta_es2Up);
      t->SetBranchAddress("D_zeta_es2Down", &D_zeta_es2Down);
      t->SetBranchAddress("D_zeta_responseUp", &D_zeta_responseUp);
      t->SetBranchAddress("D_zeta_responseDown", &D_zeta_responseDown);
      t->SetBranchAddress("D_zeta_resolutionUp", &D_zeta_resolutionUp);
      t->SetBranchAddress("D_zeta_resolutionDown", &D_zeta_resolutionDown);
      t->SetBranchAddress("D_zeta_UESUp", &D_zeta_UESUp);
      t->SetBranchAddress("D_zeta_UESDown", &D_zeta_UESDown);
      t->SetBranchAddress("D_zeta_JERUp", &D_zeta_JERUp);
      t->SetBranchAddress("D_zeta_JERDown", &D_zeta_JERDown);
      t->SetBranchAddress("D_zeta_JetAbsoluteUp", &D_zeta_JetAbsoluteUp);
      t->SetBranchAddress("D_zeta_JetAbsoluteDown", &D_zeta_JetAbsoluteDown);
      t->SetBranchAddress("D_zeta_JetAbsoluteyearUp", &D_zeta_JetAbsoluteyearUp);
      t->SetBranchAddress("D_zeta_JetAbsoluteyearDown", &D_zeta_JetAbsoluteyearDown);
      t->SetBranchAddress("D_zeta_JetBBEC1Up", &D_zeta_JetBBEC1Up);
      t->SetBranchAddress("D_zeta_JetBBEC1Down", &D_zeta_JetBBEC1Down);
      t->SetBranchAddress("D_zeta_JetBBEC1yearUp", &D_zeta_JetBBEC1yearUp);
      t->SetBranchAddress("D_zeta_JetBBEC1yearDown", &D_zeta_JetBBEC1yearDown);
      t->SetBranchAddress("D_zeta_JetEC2Up", &D_zeta_JetEC2Up);
      t->SetBranchAddress("D_zeta_JetEC2Down", &D_zeta_JetEC2Down);
      t->SetBranchAddress("D_zeta_JetEC2yearUp", &D_zeta_JetEC2yearUp);
      t->SetBranchAddress("D_zeta_JetEC2yearDown", &D_zeta_JetEC2yearDown);
      t->SetBranchAddress("D_zeta_JetFlavorQCDUp", &D_zeta_JetFlavorQCDUp);
      t->SetBranchAddress("D_zeta_JetFlavorQCDDown", &D_zeta_JetFlavorQCDDown);
      t->SetBranchAddress("D_zeta_JetHFUp", &D_zeta_JetHFUp);
      t->SetBranchAddress("D_zeta_JetHFDown", &D_zeta_JetHFDown);
      t->SetBranchAddress("D_zeta_JetHFyearUp", &D_zeta_JetHFyearUp);
      t->SetBranchAddress("D_zeta_JetHFyearDown", &D_zeta_JetHFyearDown);
      t->SetBranchAddress("D_zeta_JetRelativeBalUp", &D_zeta_JetRelativeBalUp);
      t->SetBranchAddress("D_zeta_JetRelativeBalDown", &D_zeta_JetRelativeBalDown);
      t->SetBranchAddress("D_zeta_JetRelativeSampleUp", &D_zeta_JetRelativeSampleUp);
      t->SetBranchAddress("D_zeta_JetRelativeSampleDown", &D_zeta_JetRelativeSampleDown);

      // New Change
      t->SetBranchAddress("bpt_deepflavour_2", &bpt_deepflavour_2);
      t->SetBranchAddress("bpt_deepflavour_JERUp_2", &bpt_deepflavour_JERUp_2);
      t->SetBranchAddress("bpt_deepflavour_JERDown_2", &bpt_deepflavour_JERDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteUp_2", &bpt_deepflavour_JetAbsoluteUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteDown_2", &bpt_deepflavour_JetAbsoluteDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearUp_2", &bpt_deepflavour_JetAbsoluteyearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetAbsoluteyearDown_2", &bpt_deepflavour_JetAbsoluteyearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Up_2", &bpt_deepflavour_JetBBEC1Up_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1Down_2", &bpt_deepflavour_JetBBEC1Down_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearUp_2", &bpt_deepflavour_JetBBEC1yearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetBBEC1yearDown_2", &bpt_deepflavour_JetBBEC1yearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Up_2", &bpt_deepflavour_JetEC2Up_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2Down_2", &bpt_deepflavour_JetEC2Down_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearUp_2", &bpt_deepflavour_JetEC2yearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetEC2yearDown_2", &bpt_deepflavour_JetEC2yearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDUp_2", &bpt_deepflavour_JetFlavorQCDUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetFlavorQCDDown_2", &bpt_deepflavour_JetFlavorQCDDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFUp_2", &bpt_deepflavour_JetHFUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFDown_2", &bpt_deepflavour_JetHFDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearUp_2", &bpt_deepflavour_JetHFyearUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetHFyearDown_2", &bpt_deepflavour_JetHFyearDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalUp_2", &bpt_deepflavour_JetRelativeBalUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeBalDown_2", &bpt_deepflavour_JetRelativeBalDown_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleUp_2", &bpt_deepflavour_JetRelativeSampleUp_2);
      t->SetBranchAddress("bpt_deepflavour_JetRelativeSampleDown_2", &bpt_deepflavour_JetRelativeSampleDown_2);
      t->SetBranchAddress("beta_deepflavour_2", &beta_deepflavour_2);
      t->SetBranchAddress("bphi_deepflavour_2", &bphi_deepflavour_2);
      t->SetBranchAddress("bm_deepflavour_2", &bm_deepflavour_2);
      t->SetBranchAddress("bm_deepflavour_JERUp_2", &bm_deepflavour_JERUp_2);
      t->SetBranchAddress("bm_deepflavour_JERDown_2", &bm_deepflavour_JERDown_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteUp_2", &bm_deepflavour_JetAbsoluteUp_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteDown_2", &bm_deepflavour_JetAbsoluteDown_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearUp_2", &bm_deepflavour_JetAbsoluteyearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetAbsoluteyearDown_2", &bm_deepflavour_JetAbsoluteyearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Up_2", &bm_deepflavour_JetBBEC1Up_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1Down_2", &bm_deepflavour_JetBBEC1Down_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearUp_2", &bm_deepflavour_JetBBEC1yearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetBBEC1yearDown_2", &bm_deepflavour_JetBBEC1yearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2Up_2", &bm_deepflavour_JetEC2Up_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2Down_2", &bm_deepflavour_JetEC2Down_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearUp_2", &bm_deepflavour_JetEC2yearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetEC2yearDown_2", &bm_deepflavour_JetEC2yearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDUp_2", &bm_deepflavour_JetFlavorQCDUp_2);
      t->SetBranchAddress("bm_deepflavour_JetFlavorQCDDown_2", &bm_deepflavour_JetFlavorQCDDown_2);
      t->SetBranchAddress("bm_deepflavour_JetHFUp_2", &bm_deepflavour_JetHFUp_2);
      t->SetBranchAddress("bm_deepflavour_JetHFDown_2", &bm_deepflavour_JetHFDown_2);
      t->SetBranchAddress("bm_deepflavour_JetHFyearUp_2", &bm_deepflavour_JetHFyearUp_2);
      t->SetBranchAddress("bm_deepflavour_JetHFyearDown_2", &bm_deepflavour_JetHFyearDown_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalUp_2", &bm_deepflavour_JetRelativeBalUp_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeBalDown_2", &bm_deepflavour_JetRelativeBalDown_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleUp_2", &bm_deepflavour_JetRelativeSampleUp_2);
      t->SetBranchAddress("bm_deepflavour_JetRelativeSampleDown_2", &bm_deepflavour_JetRelativeSampleDown_2);
      t->SetBranchAddress("njets", &njets);
      t->SetBranchAddress("njets_JERUp", &njets_JERUp);
      t->SetBranchAddress("njets_JERDown", &njets_JERDown);
      t->SetBranchAddress("njets_JetAbsoluteUp", &njets_JetAbsoluteUp);
      t->SetBranchAddress("njets_JetAbsoluteDown", &njets_JetAbsoluteDown);
      t->SetBranchAddress("njets_JetAbsoluteyearUp", &njets_JetAbsoluteyearUp);
      t->SetBranchAddress("njets_JetAbsoluteyearDown", &njets_JetAbsoluteyearDown);
      t->SetBranchAddress("njets_JetBBEC1Up", &njets_JetBBEC1Up);
      t->SetBranchAddress("njets_JetBBEC1Down", &njets_JetBBEC1Down);
      t->SetBranchAddress("njets_JetBBEC1yearUp", &njets_JetBBEC1yearUp);
      t->SetBranchAddress("njets_JetBBEC1yearDown", &njets_JetBBEC1yearDown);
      t->SetBranchAddress("njets_JetEC2Up", &njets_JetEC2Up);
      t->SetBranchAddress("njets_JetEC2Down", &njets_JetEC2Down);
      t->SetBranchAddress("njets_JetEC2yearUp", &njets_JetEC2yearUp);
      t->SetBranchAddress("njets_JetEC2yearDown", &njets_JetEC2yearDown);
      t->SetBranchAddress("njets_JetFlavorQCDUp", &njets_JetFlavorQCDUp);
      t->SetBranchAddress("njets_JetFlavorQCDDown", &njets_JetFlavorQCDDown);
      t->SetBranchAddress("njets_JetHFUp", &njets_JetHFUp);
      t->SetBranchAddress("njets_JetHFDown", &njets_JetHFDown);
      t->SetBranchAddress("njets_JetHFyearUp", &njets_JetHFyearUp);
      t->SetBranchAddress("njets_JetHFyearDown", &njets_JetHFyearDown);
      t->SetBranchAddress("njets_JetRelativeBalUp", &njets_JetRelativeBalUp);
      t->SetBranchAddress("njets_JetRelativeBalDown", &njets_JetRelativeBalDown);
      t->SetBranchAddress("njets_JetRelativeSampleUp", &njets_JetRelativeSampleUp);
      t->SetBranchAddress("njets_JetRelativeSampleDown", &njets_JetRelativeSampleDown);
      t->SetBranchAddress("pt_vis_nominal", &pt_vis_nominal);
      t->SetBranchAddress("pt_vis_es1Up", &pt_vis_es1Up);
      t->SetBranchAddress("pt_vis_es1Down", &pt_vis_es1Down);
      t->SetBranchAddress("pt_vis_es2Up", &pt_vis_es2Up);
      t->SetBranchAddress("pt_vis_es2Down", &pt_vis_es2Down);
      t->SetBranchAddress("m_btautau_vis_nominal", &m_btautau_vis_nominal);
      t->SetBranchAddress("m_btautau_vis_es1Up", &m_btautau_vis_es1Up);
      t->SetBranchAddress("m_btautau_vis_es1Down", &m_btautau_vis_es1Down);
      t->SetBranchAddress("m_btautau_vis_es2Up", &m_btautau_vis_es2Up);
      t->SetBranchAddress("m_btautau_vis_es2Down", &m_btautau_vis_es2Down);
      t->SetBranchAddress("m_btautau_vis_JetAbsoluteUp_1", &m_btautau_vis_JetAbsoluteUp_1);
      t->SetBranchAddress("m_btautau_vis_JetAbsoluteDown_1", &m_btautau_vis_JetAbsoluteDown_1);
      t->SetBranchAddress("m_btautau_vis_JetAbsoluteyearUp_1", &m_btautau_vis_JetAbsoluteyearUp_1);
      t->SetBranchAddress("m_btautau_vis_JetAbsoluteyearDown_1", &m_btautau_vis_JetAbsoluteyearDown_1);
      t->SetBranchAddress("m_btautau_vis_JetBBEC1Up_1", &m_btautau_vis_JetBBEC1Up_1);
      t->SetBranchAddress("m_btautau_vis_JetBBEC1Down_1", &m_btautau_vis_JetBBEC1Down_1);
      t->SetBranchAddress("m_btautau_vis_JetBBEC1yearUp_1", &m_btautau_vis_JetBBEC1yearUp_1);
      t->SetBranchAddress("m_btautau_vis_JetBBEC1yearDown_1", &m_btautau_vis_JetBBEC1yearDown_1);
      t->SetBranchAddress("m_btautau_vis_JetEC2Up_1", &m_btautau_vis_JetEC2Up_1);
      t->SetBranchAddress("m_btautau_vis_JetEC2Down_1", &m_btautau_vis_JetEC2Down_1);
      t->SetBranchAddress("m_btautau_vis_JetEC2yearUp_1", &m_btautau_vis_JetEC2yearUp_1);
      t->SetBranchAddress("m_btautau_vis_JetEC2yearDown_1", &m_btautau_vis_JetEC2yearDown_1);
      t->SetBranchAddress("m_btautau_vis_JetFlavorQCDUp_1", &m_btautau_vis_JetFlavorQCDUp_1);
      t->SetBranchAddress("m_btautau_vis_JetFlavorQCDDown_1", &m_btautau_vis_JetFlavorQCDDown_1);
      t->SetBranchAddress("m_btautau_vis_JetHFUp_1", &m_btautau_vis_JetHFUp_1);
      t->SetBranchAddress("m_btautau_vis_JetHFDown_1", &m_btautau_vis_JetHFDown_1);
      t->SetBranchAddress("m_btautau_vis_JetHFyearUp_1", &m_btautau_vis_JetHFyearUp_1);
      t->SetBranchAddress("m_btautau_vis_JetHFyearDown_1", &m_btautau_vis_JetHFyearDown_1);
      t->SetBranchAddress("m_btautau_vis_JetRelativeBalUp_1", &m_btautau_vis_JetRelativeBalUp_1);
      t->SetBranchAddress("m_btautau_vis_JetRelativeBalDown_1", &m_btautau_vis_JetRelativeBalDown_1);
      t->SetBranchAddress("m_btautau_vis_JetRelativeSampleUp_1", &m_btautau_vis_JetRelativeSampleUp_1);
      t->SetBranchAddress("m_btautau_vis_JetRelativeSampleDown_1", &m_btautau_vis_JetRelativeSampleDown_1);
      t->SetBranchAddress("m_btautau_vis_JERUp_1", &m_btautau_vis_JERUp_1);
      t->SetBranchAddress("m_btautau_vis_JERDown_1", &m_btautau_vis_JERDown_1);
      t->SetBranchAddress("mtMET_1_nominal", &mtMET_1_nominal);
      t->SetBranchAddress("mtMET_1_es1Up", &mtMET_1_es1Up);
      t->SetBranchAddress("mtMET_1_es1Down", &mtMET_1_es1Down);
      t->SetBranchAddress("mtMET_1_responseUp", &mtMET_1_responseUp);
      t->SetBranchAddress("mtMET_1_responseDown", &mtMET_1_responseDown);
      t->SetBranchAddress("mtMET_1_resolutionUp", &mtMET_1_resolutionUp);
      t->SetBranchAddress("mtMET_1_resolutionDown", &mtMET_1_resolutionDown);
      t->SetBranchAddress("mtMET_1_JetAbsoluteUp", &mtMET_1_JetAbsoluteUp);
      t->SetBranchAddress("mtMET_1_JetAbsoluteDown", &mtMET_1_JetAbsoluteDown);
      t->SetBranchAddress("mtMET_1_JetAbsoluteyearUp", &mtMET_1_JetAbsoluteyearUp);
      t->SetBranchAddress("mtMET_1_JetAbsoluteyearDown", &mtMET_1_JetAbsoluteyearDown);
      t->SetBranchAddress("mtMET_1_JetBBEC1Up", &mtMET_1_JetBBEC1Up);
      t->SetBranchAddress("mtMET_1_JetBBEC1Down", &mtMET_1_JetBBEC1Down);
      t->SetBranchAddress("mtMET_1_JetBBEC1yearUp", &mtMET_1_JetBBEC1yearUp);
      t->SetBranchAddress("mtMET_1_JetBBEC1yearDown", &mtMET_1_JetBBEC1yearDown);
      t->SetBranchAddress("mtMET_1_JetEC2Up", &mtMET_1_JetEC2Up);
      t->SetBranchAddress("mtMET_1_JetEC2Down",&mtMET_1_JetEC2Down);
      t->SetBranchAddress("mtMET_1_JetEC2yearUp", &mtMET_1_JetEC2yearUp);
      t->SetBranchAddress("mtMET_1_JetEC2yearDown", &mtMET_1_JetEC2yearDown);
      t->SetBranchAddress("mtMET_1_JetFlavorQCDUp", &mtMET_1_JetFlavorQCDUp);
      t->SetBranchAddress("mtMET_1_JetFlavorQCDDown", &mtMET_1_JetFlavorQCDDown);
      t->SetBranchAddress("mtMET_1_JetHFUp", &mtMET_1_JetHFUp);
      t->SetBranchAddress("mtMET_1_JetHFDown", &mtMET_1_JetHFDown);
      t->SetBranchAddress("mtMET_1_JetHFyearUp", &mtMET_1_JetHFyearUp);
      t->SetBranchAddress("mtMET_1_JetHFyearDown", &mtMET_1_JetHFyearDown);
      t->SetBranchAddress("mtMET_1_JetRelativeBalUp", &mtMET_1_JetRelativeBalUp);
      t->SetBranchAddress("mtMET_1_JetRelativeBalDown", &mtMET_1_JetRelativeBalDown);
      t->SetBranchAddress("mtMET_1_JetRelativeSampleUp", &mtMET_1_JetRelativeSampleUp);
      t->SetBranchAddress("mtMET_1_JetRelativeSampleDown", &mtMET_1_JetRelativeSampleDown);
      t->SetBranchAddress("mtMET_1_JERUp", &mtMET_1_JERUp);
      t->SetBranchAddress("mtMET_1_JERDown", &mtMET_1_JERDown);
      t->SetBranchAddress("mtMET_2_nominal", &mtMET_2_nominal);
      t->SetBranchAddress("mtMET_2_es2Up", &mtMET_2_es2Up);
      t->SetBranchAddress("mtMET_2_es2Down", &mtMET_2_es2Down);
      t->SetBranchAddress("mtMET_2_responseUp", &mtMET_2_responseUp);
      t->SetBranchAddress("mtMET_2_responseDown", &mtMET_2_responseDown);
      t->SetBranchAddress("mtMET_2_resolutionUp", &mtMET_2_resolutionUp);
      t->SetBranchAddress("mtMET_2_resolutionDown", &mtMET_2_resolutionDown);
      t->SetBranchAddress("mtMET_2_JetAbsoluteUp", &mtMET_2_JetAbsoluteUp);
      t->SetBranchAddress("mtMET_2_JetAbsoluteDown", &mtMET_2_JetAbsoluteDown);
      t->SetBranchAddress("mtMET_2_JetAbsoluteyearUp", &mtMET_2_JetAbsoluteyearUp);
      t->SetBranchAddress("mtMET_2_JetAbsoluteyearDown", &mtMET_2_JetAbsoluteyearDown);
      t->SetBranchAddress("mtMET_2_JetBBEC1Up", &mtMET_2_JetBBEC1Up);
      t->SetBranchAddress("mtMET_2_JetBBEC1Down", &mtMET_2_JetBBEC1Down);
      t->SetBranchAddress("mtMET_2_JetBBEC1yearUp", &mtMET_2_JetBBEC1yearUp);
      t->SetBranchAddress("mtMET_2_JetBBEC1yearDown", &mtMET_2_JetBBEC1yearDown);
      t->SetBranchAddress("mtMET_2_JetEC2Up", &mtMET_2_JetEC2Up);
      t->SetBranchAddress("mtMET_2_JetEC2Down",&mtMET_2_JetEC2Down);
      t->SetBranchAddress("mtMET_2_JetEC2yearUp", &mtMET_2_JetEC2yearUp);
      t->SetBranchAddress("mtMET_2_JetEC2yearDown", &mtMET_2_JetEC2yearDown);
      t->SetBranchAddress("mtMET_2_JetFlavorQCDUp", &mtMET_2_JetFlavorQCDUp);
      t->SetBranchAddress("mtMET_2_JetFlavorQCDDown", &mtMET_2_JetFlavorQCDDown);
      t->SetBranchAddress("mtMET_2_JetHFUp", &mtMET_2_JetHFUp);
      t->SetBranchAddress("mtMET_2_JetHFDown", &mtMET_2_JetHFDown);
      t->SetBranchAddress("mtMET_2_JetHFyearUp", &mtMET_2_JetHFyearUp);
      t->SetBranchAddress("mtMET_2_JetHFyearDown", &mtMET_2_JetHFyearDown);
      t->SetBranchAddress("mtMET_2_JetRelativeBalUp", &mtMET_2_JetRelativeBalUp);
      t->SetBranchAddress("mtMET_2_JetRelativeBalDown", &mtMET_2_JetRelativeBalDown);
      t->SetBranchAddress("mtMET_2_JetRelativeSampleUp", &mtMET_2_JetRelativeSampleUp);
      t->SetBranchAddress("mtMET_2_JetRelativeSampleDown", &mtMET_2_JetRelativeSampleDown);
      t->SetBranchAddress("mtMET_2_JERUp", &mtMET_2_JERUp);
      t->SetBranchAddress("mtMET_2_JERDown", &mtMET_2_JERDown);


      printf("Found tree -> weighting\n");
      for (Int_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        BDTval_nominal = -10;
        BDTval_es1Up = -10;
        BDTval_es1Down = -10;
        BDTval_es2Up = -10;
        BDTval_es2Down = -10;
        BDTval_UESUp = -10;
        BDTval_UESDown = -10;
        BDTval_ResponseUp = -10;
        BDTval_ResponseDown = -10;
        BDTval_ResolutionUp = -10;
        BDTval_ResolutionDown = -10;
        BDTval_JetAbsoluteUp = -10;
        BDTval_JetAbsoluteDown = -10;
        BDTval_JetAbsoluteyearUp = -10;
        BDTval_JetAbsoluteyearDown = -10;
        BDTval_JetBBEC1Up = -10;
        BDTval_JetBBEC1Down = -10;
        BDTval_JetBBEC1yearUp = -10;
        BDTval_JetBBEC1yearDown = -10;
        BDTval_JetEC2Up = -10;
        BDTval_JetEC2Down = -10;
        BDTval_JetEC2yearUp = -10;
        BDTval_JetEC2yearDown = -10;
        BDTval_JetFlavorQCDUp = -10;
        BDTval_JetFlavorQCDDown = -10;
        BDTval_JetHFUp = -10;
        BDTval_JetHFDown = -10;
        BDTval_JetHFyearUp = -10;
        BDTval_JetHFyearDown = -10;
        BDTval_JetRelativeBalUp = -10;
        BDTval_JetRelativeBalDown = -10;
        BDTval_JetRelativeSampleUp = -10;
        BDTval_JetRelativeSampleDown = -10;
        BDTval_JERUp = -10;
        BDTval_JERDown = -10;

        int booster_idx = channel;
        if (bpt_deepflavour_1 > 0 && bpt_deepflavour_2 > 0) booster_idx += 3;

	if (bpt_deepflavour_1 > 0) {
	  runBDTeval(BDTval_nominal, features_nominal, booster_idx, boosters[booster_idx], bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1, bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2, pt_1_nominal, eta_1, phi_1, m_1_nominal, pt_2_nominal, eta_2, phi_2, m_2_nominal, met_nominal, metphi_nominal, D_zeta_nominal, njets, m_btautau_vis_nominal, mtMET_1_nominal, mtMET_2_nominal, pt_vis_nominal);
	  if (isMC) {
            runBDTeval(BDTval_es1Up, features_es1Up, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_es1Up, eta_1, phi_1, m_1_es1Up,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_es1Up, metphi_es1Up, D_zeta_es1Up, njets,
            		m_btautau_vis_es1Up, mtMET_1_es1Up, mtMET_2_nominal, pt_vis_es1Up);
            runBDTeval(BDTval_es1Down, features_es1Down, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_es1Down, eta_1, phi_1, m_1_es1Down,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_es1Down, metphi_es1Down, D_zeta_es1Down, njets,
            		m_btautau_vis_es1Down, mtMET_1_es1Down, mtMET_2_nominal, pt_vis_es1Down);
            runBDTeval(BDTval_es2Up, features_es2Up, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_es2Up, eta_2, phi_2, m_2_es2Up,
            		met_es2Up, metphi_es2Up, D_zeta_es2Up, njets,
            		m_btautau_vis_es2Up, mtMET_1_nominal, mtMET_2_es2Up, pt_vis_es2Up);
            runBDTeval(BDTval_es2Down, features_es2Down, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_es2Down, eta_2, phi_2, m_2_es2Down,
            		met_es2Down, metphi_es2Down, D_zeta_es2Down, njets,
            		m_btautau_vis_es2Down, mtMET_1_nominal, mtMET_2_es2Down, pt_vis_es2Down);
            runBDTeval(BDTval_UESUp, features_UESUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_UESUp, metphi_UESUp, D_zeta_UESUp, njets,
            		m_btautau_vis_nominal, mtMET_1_nominal, mtMET_2_nominal, pt_vis_nominal);
            runBDTeval(BDTval_UESDown, features_UESDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_UESDown, metphi_UESDown, D_zeta_UESDown, njets,
            		m_btautau_vis_nominal, mtMET_1_nominal, mtMET_2_nominal, pt_vis_nominal);
            runBDTeval(BDTval_ResponseUp, features_responseUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_responseUp, metphi_responseUp, D_zeta_responseUp, njets,
            		m_btautau_vis_nominal, mtMET_1_responseUp, mtMET_2_responseUp, pt_vis_nominal);
            runBDTeval(BDTval_ResponseDown, features_responseDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_responseDown, metphi_responseDown, D_zeta_responseDown, njets,
            		m_btautau_vis_nominal, mtMET_1_responseDown, mtMET_2_responseDown, pt_vis_nominal);
            runBDTeval(BDTval_ResolutionUp, features_resolutionUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_resolutionUp, metphi_resolutionUp, D_zeta_resolutionUp, njets,
            		m_btautau_vis_nominal, mtMET_1_resolutionUp, mtMET_2_resolutionUp, pt_vis_nominal);
            runBDTeval(BDTval_ResolutionDown, features_resolutionDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
            		bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_resolutionDown, metphi_resolutionDown, D_zeta_resolutionDown, njets,
            		m_btautau_vis_nominal, mtMET_1_resolutionDown, mtMET_2_resolutionDown, pt_vis_nominal);
            runBDTeval(BDTval_JetAbsoluteUp, features_JetAbsoluteUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetAbsoluteUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteUp_1,
            		bpt_deepflavour_JetAbsoluteUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetAbsoluteUp, metphi_JetAbsoluteUp, D_zeta_JetAbsoluteUp, njets_JetAbsoluteUp,
            		m_btautau_vis_JetAbsoluteUp_1, mtMET_1_JetAbsoluteUp, mtMET_2_JetAbsoluteUp, pt_vis_nominal);
            runBDTeval(BDTval_JetAbsoluteDown, features_JetAbsoluteDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetAbsoluteDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteDown_1,
            		bpt_deepflavour_JetAbsoluteDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetAbsoluteDown, metphi_JetAbsoluteDown, D_zeta_JetAbsoluteDown, njets_JetAbsoluteDown,
            		m_btautau_vis_JetAbsoluteDown_1, mtMET_1_JetAbsoluteDown, mtMET_2_JetAbsoluteDown, pt_vis_nominal);
            runBDTeval(BDTval_JetAbsoluteyearUp, features_JetAbsoluteyearUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetAbsoluteyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearUp_1,
            		bpt_deepflavour_JetAbsoluteyearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteyearUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetAbsoluteyearUp, metphi_JetAbsoluteyearUp, D_zeta_JetAbsoluteyearUp, njets_JetAbsoluteyearUp,
            		m_btautau_vis_JetAbsoluteyearUp_1, mtMET_1_JetAbsoluteyearUp, mtMET_2_JetAbsoluteyearUp, pt_vis_nominal);
            runBDTeval(BDTval_JetAbsoluteyearDown, features_JetAbsoluteyearDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetAbsoluteyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetAbsoluteyearDown_1,
            		bpt_deepflavour_JetAbsoluteyearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetAbsoluteyearDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetAbsoluteyearDown, metphi_JetAbsoluteyearDown, D_zeta_JetAbsoluteyearDown, njets_JetAbsoluteyearDown,
            		m_btautau_vis_JetAbsoluteyearDown_1, mtMET_1_JetAbsoluteyearDown, mtMET_2_JetAbsoluteyearDown, pt_vis_nominal);
            runBDTeval(BDTval_JetBBEC1Up, features_JetBBEC1Up, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetBBEC1Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Up_1,
            		bpt_deepflavour_JetBBEC1Up_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1Up_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetBBEC1Up, metphi_JetBBEC1Up, D_zeta_JetBBEC1Up, njets_JetBBEC1Up,
            		m_btautau_vis_JetBBEC1Up_1, mtMET_1_JetBBEC1Up, mtMET_2_JetBBEC1Up, pt_vis_nominal);
            runBDTeval(BDTval_JetBBEC1Down, features_JetBBEC1Down, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetBBEC1Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1Down_1,
            		bpt_deepflavour_JetBBEC1Down_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1Down_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetBBEC1Down, metphi_JetBBEC1Down, D_zeta_JetBBEC1Down, njets_JetBBEC1Down,
            		m_btautau_vis_JetBBEC1Down_1, mtMET_1_JetBBEC1Down, mtMET_2_JetBBEC1Down, pt_vis_nominal);
            runBDTeval(BDTval_JetBBEC1yearUp, features_JetBBEC1yearUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetBBEC1yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearUp_1,
            		bpt_deepflavour_JetBBEC1yearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1yearUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetBBEC1yearUp, metphi_JetBBEC1yearUp, D_zeta_JetBBEC1yearUp, njets_JetBBEC1yearUp,
            		m_btautau_vis_JetBBEC1yearUp_1, mtMET_1_JetBBEC1yearUp, mtMET_2_JetBBEC1yearUp, pt_vis_nominal);
            runBDTeval(BDTval_JetBBEC1yearDown, features_JetBBEC1yearDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetBBEC1yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetBBEC1yearDown_1,
            		bpt_deepflavour_JetBBEC1yearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetBBEC1yearDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetBBEC1yearDown, metphi_JetBBEC1yearDown, D_zeta_JetBBEC1yearDown, njets_JetBBEC1yearDown,
            		m_btautau_vis_JetBBEC1yearDown_1, mtMET_1_JetBBEC1yearDown, mtMET_2_JetBBEC1yearDown, pt_vis_nominal);
            runBDTeval(BDTval_JetEC2Up, features_JetEC2Up, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetEC2Up_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Up_1,
            		bpt_deepflavour_JetEC2Up_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2Up_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetEC2Up, metphi_JetEC2Up, D_zeta_JetEC2Up, njets_JetEC2Up,
            		m_btautau_vis_JetEC2Up_1, mtMET_1_JetEC2Up, mtMET_2_JetEC2Up, pt_vis_nominal);
            runBDTeval(BDTval_JetEC2Down, features_JetEC2Down, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetEC2Down_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2Down_1,
            		bpt_deepflavour_JetEC2Down_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2Down_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetEC2Down, metphi_JetEC2Down, D_zeta_JetEC2Down, njets_JetEC2Down,
            		m_btautau_vis_JetEC2Down_1, mtMET_1_JetEC2Down, mtMET_2_JetEC2Down, pt_vis_nominal);
            runBDTeval(BDTval_JetEC2yearUp, features_JetEC2yearUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetEC2yearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearUp_1,
            		bpt_deepflavour_JetEC2yearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2yearUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetEC2yearUp, metphi_JetEC2yearUp, D_zeta_JetEC2yearUp, njets_JetEC2yearUp,
            		m_btautau_vis_JetEC2yearUp_1, mtMET_1_JetEC2yearUp, mtMET_2_JetEC2yearUp, pt_vis_nominal);
            runBDTeval(BDTval_JetEC2yearDown, features_JetEC2yearDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetEC2yearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetEC2yearDown_1,
            		bpt_deepflavour_JetEC2yearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetEC2yearDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetEC2yearDown, metphi_JetEC2yearDown, D_zeta_JetEC2yearDown, njets_JetEC2yearDown,
            		m_btautau_vis_JetEC2yearDown_1, mtMET_1_JetEC2yearDown, mtMET_2_JetEC2yearDown, pt_vis_nominal);
            runBDTeval(BDTval_JetFlavorQCDUp, features_JetFlavorQCDUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetFlavorQCDUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDUp_1,
            		bpt_deepflavour_JetFlavorQCDUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetFlavorQCDUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetFlavorQCDUp, metphi_JetFlavorQCDUp, D_zeta_JetFlavorQCDUp, njets_JetFlavorQCDUp,
            		m_btautau_vis_JetFlavorQCDUp_1, mtMET_1_JetFlavorQCDUp, mtMET_2_JetFlavorQCDUp, pt_vis_nominal);
            runBDTeval(BDTval_JetFlavorQCDDown, features_JetFlavorQCDDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetFlavorQCDDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetFlavorQCDDown_1,
            		bpt_deepflavour_JetFlavorQCDDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetFlavorQCDDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetFlavorQCDDown, metphi_JetFlavorQCDDown, D_zeta_JetFlavorQCDDown, njets_JetFlavorQCDDown,
            		m_btautau_vis_JetFlavorQCDDown_1, mtMET_1_JetFlavorQCDDown, mtMET_2_JetFlavorQCDDown, pt_vis_nominal);
            runBDTeval(BDTval_JetHFUp, features_JetHFUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetHFUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFUp_1,
            		bpt_deepflavour_JetHFUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetHFUp, metphi_JetHFUp, D_zeta_JetHFUp, njets_JetHFUp,
            		m_btautau_vis_JetHFUp_1, mtMET_1_JetHFUp, mtMET_2_JetHFUp, pt_vis_nominal);
            runBDTeval(BDTval_JetHFDown, features_JetHFDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetHFDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFDown_1,
            		bpt_deepflavour_JetHFDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetHFDown, metphi_JetHFDown, D_zeta_JetHFDown, njets_JetHFDown,
            		m_btautau_vis_JetHFDown_1, mtMET_1_JetHFDown, mtMET_2_JetHFDown, pt_vis_nominal);
            runBDTeval(BDTval_JetHFyearUp, features_JetHFyearUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetHFyearUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearUp_1,
            		bpt_deepflavour_JetHFyearUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFyearUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetHFyearUp, metphi_JetHFyearUp, D_zeta_JetHFyearUp, njets_JetHFyearUp,
            		m_btautau_vis_JetHFyearUp_1, mtMET_1_JetHFyearUp, mtMET_2_JetHFyearUp, pt_vis_nominal);
            runBDTeval(BDTval_JetHFyearDown, features_JetHFyearDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetHFyearDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetHFyearDown_1,
            		bpt_deepflavour_JetHFyearDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetHFyearDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetHFyearDown, metphi_JetHFyearDown, D_zeta_JetHFyearDown, njets_JetHFyearDown,
            		m_btautau_vis_JetHFyearDown_1, mtMET_1_JetHFyearDown, mtMET_2_JetHFyearDown, pt_vis_nominal);
            runBDTeval(BDTval_JetRelativeBalUp, features_JetRelativeBalUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetRelativeBalUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalUp_1,
            		bpt_deepflavour_JetRelativeBalUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeBalUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetRelativeBalUp, metphi_JetRelativeBalUp, D_zeta_JetRelativeBalUp, njets_JetRelativeBalUp,
            		m_btautau_vis_JetRelativeBalUp_1, mtMET_1_JetRelativeBalUp, mtMET_2_JetRelativeBalUp, pt_vis_nominal);
            runBDTeval(BDTval_JetRelativeBalDown, features_JetRelativeBalDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetRelativeBalDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeBalDown_1,
            		bpt_deepflavour_JetRelativeBalDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeBalDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetRelativeBalDown, metphi_JetRelativeBalDown, D_zeta_JetRelativeBalDown, njets_JetRelativeBalDown,
            		m_btautau_vis_JetRelativeBalDown_1, mtMET_1_JetRelativeBalDown, mtMET_2_JetRelativeBalDown, pt_vis_nominal);
            runBDTeval(BDTval_JetRelativeSampleUp, features_JetRelativeSampleUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetRelativeSampleUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleUp_1,
            		bpt_deepflavour_JetRelativeSampleUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeSampleUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetRelativeSampleUp, metphi_JetRelativeSampleUp, D_zeta_JetRelativeSampleUp, njets_JetRelativeSampleUp,
            		m_btautau_vis_JetRelativeSampleUp_1, mtMET_1_JetRelativeSampleUp, mtMET_2_JetRelativeSampleUp, pt_vis_nominal);
            runBDTeval(BDTval_JetRelativeSampleDown, features_JetRelativeSampleDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JetRelativeSampleDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JetRelativeSampleDown_1,
            		bpt_deepflavour_JetRelativeSampleDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JetRelativeSampleDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JetRelativeSampleDown, metphi_JetRelativeSampleDown, D_zeta_JetRelativeSampleDown, njets_JetRelativeSampleDown,
            		m_btautau_vis_JetRelativeSampleDown_1, mtMET_1_JetRelativeSampleDown, mtMET_2_JetRelativeSampleDown, pt_vis_nominal);
            runBDTeval(BDTval_JERUp, features_JERUp, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JERUp_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERUp_1,
            		bpt_deepflavour_JERUp_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JERUp_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JERUp, metphi_JERUp, D_zeta_JERUp, njets_JERUp,
            		m_btautau_vis_JERUp_1, mtMET_1_JERUp, mtMET_2_JERUp, pt_vis_nominal);
            runBDTeval(BDTval_JERDown, features_JERDown, booster_idx, boosters[booster_idx],
            		bpt_deepflavour_JERDown_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_JERDown_1,
            		bpt_deepflavour_JERDown_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_JERDown_2,
            		pt_1_nominal, eta_1, phi_1, m_1_nominal,
            		pt_2_nominal, eta_2, phi_2, m_2_nominal,
            		met_JERDown, metphi_JERDown, D_zeta_JERDown, njets_JERDown,
            		m_btautau_vis_JERDown_1, mtMET_1_JERDown, mtMET_2_JERDown, pt_vis_nominal);
	  }
	  if (isEmbedded) {
            runBDTeval(BDTval_es1Up, features_es1Up, booster_idx, boosters[booster_idx],
                        bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
                        bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
                        pt_1_es1Up, eta_1, phi_1, m_1_es1Up,
                        pt_2_nominal, eta_2, phi_2, m_2_nominal,
                        met_es1Up, metphi_es1Up, D_zeta_es1Up, njets,
                        m_btautau_vis_es1Up, mtMET_1_es1Up, mtMET_2_nominal, pt_vis_es1Up);
            runBDTeval(BDTval_es1Down, features_es1Down, booster_idx, boosters[booster_idx],
                        bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
                        bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
                        pt_1_es1Down, eta_1, phi_1, m_1_es1Down,
                        pt_2_nominal, eta_2, phi_2, m_2_nominal,
                        met_es1Down, metphi_es1Down, D_zeta_es1Down, njets,
                        m_btautau_vis_es1Down, mtMET_1_es1Down, mtMET_2_nominal, pt_vis_es1Down);
            runBDTeval(BDTval_es2Up, features_es2Up, booster_idx, boosters[booster_idx],
                        bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
                        bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
                        pt_1_nominal, eta_1, phi_1, m_1_nominal,
                        pt_2_es2Up, eta_2, phi_2, m_2_es2Up,
                        met_es2Up, metphi_es2Up, D_zeta_es2Up, njets,
                        m_btautau_vis_es2Up, mtMET_1_nominal, mtMET_2_es2Up, pt_vis_es2Up);
            runBDTeval(BDTval_es2Down, features_es2Down, booster_idx, boosters[booster_idx],
                        bpt_deepflavour_1, beta_deepflavour_1, bphi_deepflavour_1, bm_deepflavour_1,
                        bpt_deepflavour_2, beta_deepflavour_2, bphi_deepflavour_2, bm_deepflavour_2,
                        pt_1_nominal, eta_1, phi_1, m_1_nominal,
                        pt_2_es2Down, eta_2, phi_2, m_2_es2Down,
                        met_es2Down, metphi_es2Down, D_zeta_es2Down, njets,
                        m_btautau_vis_es2Down, mtMET_1_nominal, mtMET_2_es2Down, pt_vis_es2Down);
	  }
        } 

        newBranch1->Fill();
        newBranch1_es1U->Fill();
        newBranch1_es1D->Fill();
        newBranch1_es2U->Fill();
        newBranch1_es2D->Fill();
        newBranch1_UESU->Fill();
        newBranch1_UESD->Fill();
        newBranch1_ResponseU->Fill();
        newBranch1_ResponseD->Fill();
        newBranch1_ResolutionU->Fill();
        newBranch1_ResolutionD->Fill();
        newBranch1_JetAbsoluteU->Fill();
        newBranch1_JetAbsoluteD->Fill();
        newBranch1_JetAbsoluteyearU->Fill();
        newBranch1_JetAbsoluteyearD->Fill();
        newBranch1_JetBBEC1U->Fill();
        newBranch1_JetBBEC1D->Fill();
        newBranch1_JetBBEC1yearU->Fill();
        newBranch1_JetBBEC1yearD->Fill();
        newBranch1_JetEC2U->Fill();
        newBranch1_JetEC2D->Fill();
        newBranch1_JetEC2yearU->Fill();
        newBranch1_JetEC2yearD->Fill();
        newBranch1_JetFlavorQCDU->Fill();
        newBranch1_JetFlavorQCDD->Fill();
        newBranch1_JetHFU->Fill();
        newBranch1_JetHFD->Fill();
        newBranch1_JetHFyearU->Fill();
        newBranch1_JetHFyearD->Fill();
        newBranch1_JetRelativeBalU->Fill();
        newBranch1_JetRelativeBalD->Fill();
        newBranch1_JetRelativeSampleU->Fill();
        newBranch1_JetRelativeSampleD->Fill();
        newBranch1_JERU->Fill();
        newBranch1_JERD->Fill();

	
	newBranch1_btautau_dR->Fill();
        newBranch1_btau1_dR->Fill();
        newBranch1_btau2_dR->Fill();
        newBranch1_mtMET_b_nominal->Fill();
        newBranch1_m_bb_nominal->Fill();
        newBranch1_m_bbtautau_vis_nominal->Fill();
        newBranch1_m_b2tautau_vis_nominal->Fill();
        newBranch1_d_ma_nominal->Fill();
        newBranch1_b2tau1_dR->Fill();
        newBranch1_b2tau2_dR->Fill();
        newBranch1_b2tautau_dR->Fill();
	newBranch1_m_bbtautau_vis_es1Up->Fill();
	newBranch1_m_b2tautau_vis_es1Up->Fill();
	newBranch1_d_ma_es1Up->Fill();
	newBranch1_m_bbtautau_vis_es1Down->Fill();
	newBranch1_m_b2tautau_vis_es1Down->Fill();
	newBranch1_d_ma_es1Down->Fill();
	newBranch1_m_bbtautau_vis_es2Up->Fill();
        newBranch1_m_b2tautau_vis_es2Up->Fill();
        newBranch1_d_ma_es2Up->Fill();
        newBranch1_m_bbtautau_vis_es2Down->Fill();
        newBranch1_m_b2tautau_vis_es2Down->Fill();
        newBranch1_d_ma_es2Down->Fill();
	newBranch1_mtMET_b_UESUp->Fill();
	newBranch1_mtMET_b_UESDown->Fill();
	newBranch1_mtMET_b_responseUp->Fill();
	newBranch1_mtMET_b_responseDown->Fill();
	newBranch1_mtMET_b_resolutionUp->Fill();
	newBranch1_mtMET_b_resolutionDown->Fill();
        newBranch1_mtMET_b_JetAbsoluteUp->Fill();
        newBranch1_m_bb_JetAbsoluteUp->Fill();
        newBranch1_m_bbtautau_vis_JetAbsoluteUp->Fill();
        newBranch1_m_b2tautau_vis_JetAbsoluteUp->Fill();
        newBranch1_d_ma_JetAbsoluteUp->Fill();
        newBranch1_mtMET_b_JetAbsoluteDown->Fill();
        newBranch1_m_bb_JetAbsoluteDown->Fill();
        newBranch1_m_bbtautau_vis_JetAbsoluteDown->Fill();
        newBranch1_m_b2tautau_vis_JetAbsoluteDown->Fill();
        newBranch1_d_ma_JetAbsoluteDown->Fill();
        newBranch1_mtMET_b_JetAbsoluteyearUp->Fill();
        newBranch1_m_bb_JetAbsoluteyearUp->Fill();
        newBranch1_m_bbtautau_vis_JetAbsoluteyearUp->Fill();
        newBranch1_m_b2tautau_vis_JetAbsoluteyearUp->Fill();
        newBranch1_d_ma_JetAbsoluteyearUp->Fill();
        newBranch1_mtMET_b_JetAbsoluteyearDown->Fill();
        newBranch1_m_bb_JetAbsoluteyearDown->Fill();
        newBranch1_m_bbtautau_vis_JetAbsoluteyearDown->Fill();
        newBranch1_m_b2tautau_vis_JetAbsoluteyearDown->Fill();
        newBranch1_d_ma_JetAbsoluteyearDown->Fill();
        newBranch1_mtMET_b_JetBBEC1Up->Fill();
        newBranch1_m_bb_JetBBEC1Up->Fill();
        newBranch1_m_bbtautau_vis_JetBBEC1Up->Fill();
        newBranch1_m_b2tautau_vis_JetBBEC1Up->Fill();
        newBranch1_d_ma_JetBBEC1Up->Fill();
        newBranch1_mtMET_b_JetBBEC1Down->Fill();
        newBranch1_m_bb_JetBBEC1Down->Fill();
        newBranch1_m_bbtautau_vis_JetBBEC1Down->Fill();
        newBranch1_m_b2tautau_vis_JetBBEC1Down->Fill();
        newBranch1_d_ma_JetBBEC1Down->Fill();
        newBranch1_mtMET_b_JetBBEC1yearUp->Fill();
        newBranch1_m_bb_JetBBEC1yearUp->Fill();
        newBranch1_m_bbtautau_vis_JetBBEC1yearUp->Fill();
        newBranch1_m_b2tautau_vis_JetBBEC1yearUp->Fill();
        newBranch1_d_ma_JetBBEC1yearUp->Fill();
        newBranch1_mtMET_b_JetBBEC1yearDown->Fill();
        newBranch1_m_bb_JetBBEC1yearDown->Fill();
        newBranch1_m_bbtautau_vis_JetBBEC1yearDown->Fill();
        newBranch1_m_b2tautau_vis_JetBBEC1yearDown->Fill();
        newBranch1_d_ma_JetBBEC1yearDown->Fill();
        newBranch1_mtMET_b_JetEC2Up->Fill();
        newBranch1_m_bb_JetEC2Up->Fill();
        newBranch1_m_bbtautau_vis_JetEC2Up->Fill();
        newBranch1_m_b2tautau_vis_JetEC2Up->Fill();
        newBranch1_d_ma_JetEC2Up->Fill();
        newBranch1_mtMET_b_JetEC2Down->Fill();
        newBranch1_m_bb_JetEC2Down->Fill();
        newBranch1_m_bbtautau_vis_JetEC2Down->Fill();
        newBranch1_m_b2tautau_vis_JetEC2Down->Fill();
        newBranch1_d_ma_JetEC2Down->Fill();
        newBranch1_mtMET_b_JetEC2yearUp->Fill();
        newBranch1_m_bb_JetEC2yearUp->Fill();
        newBranch1_m_bbtautau_vis_JetEC2yearUp->Fill();
        newBranch1_m_b2tautau_vis_JetEC2yearUp->Fill();
        newBranch1_d_ma_JetEC2yearUp->Fill();
        newBranch1_mtMET_b_JetEC2yearDown->Fill();
        newBranch1_m_bb_JetEC2yearDown->Fill();
        newBranch1_m_bbtautau_vis_JetEC2yearDown->Fill();
        newBranch1_m_b2tautau_vis_JetEC2yearDown->Fill();
        newBranch1_d_ma_JetEC2yearDown->Fill();
        newBranch1_mtMET_b_JetFlavorQCDUp->Fill();
        newBranch1_m_bb_JetFlavorQCDUp->Fill();
        newBranch1_m_bbtautau_vis_JetFlavorQCDUp->Fill();
        newBranch1_m_b2tautau_vis_JetFlavorQCDUp->Fill();
        newBranch1_d_ma_JetFlavorQCDUp->Fill();
        newBranch1_mtMET_b_JetFlavorQCDDown->Fill();
        newBranch1_m_bb_JetFlavorQCDDown->Fill();
        newBranch1_m_bbtautau_vis_JetFlavorQCDDown->Fill();
        newBranch1_m_b2tautau_vis_JetFlavorQCDDown->Fill();
        newBranch1_d_ma_JetFlavorQCDDown->Fill();
        newBranch1_mtMET_b_JetHFUp->Fill();
        newBranch1_m_bb_JetHFUp->Fill();
        newBranch1_m_bbtautau_vis_JetHFUp->Fill();
        newBranch1_m_b2tautau_vis_JetHFUp->Fill();
        newBranch1_d_ma_JetHFUp->Fill();
        newBranch1_mtMET_b_JetHFDown->Fill();
        newBranch1_m_bb_JetHFDown->Fill();
        newBranch1_m_bbtautau_vis_JetHFDown->Fill();
        newBranch1_m_b2tautau_vis_JetHFDown->Fill();
        newBranch1_d_ma_JetHFDown->Fill();
        newBranch1_mtMET_b_JetHFyearUp->Fill();
        newBranch1_m_bb_JetHFyearUp->Fill();
        newBranch1_m_bbtautau_vis_JetHFyearUp->Fill();
        newBranch1_m_b2tautau_vis_JetHFyearUp->Fill();
        newBranch1_d_ma_JetHFyearUp->Fill();
        newBranch1_mtMET_b_JetHFyearDown->Fill();
        newBranch1_m_bb_JetHFyearDown->Fill();
        newBranch1_m_bbtautau_vis_JetHFyearDown->Fill();
        newBranch1_m_b2tautau_vis_JetHFyearDown->Fill();
        newBranch1_d_ma_JetHFyearDown->Fill();
        newBranch1_mtMET_b_JetRelativeBalUp->Fill();
        newBranch1_m_bb_JetRelativeBalUp->Fill();
        newBranch1_m_bbtautau_vis_JetRelativeBalUp->Fill();
        newBranch1_m_b2tautau_vis_JetRelativeBalUp->Fill();
        newBranch1_d_ma_JetRelativeBalUp->Fill();
        newBranch1_mtMET_b_JetRelativeBalDown->Fill();
        newBranch1_m_bb_JetRelativeBalDown->Fill();
        newBranch1_m_bbtautau_vis_JetRelativeBalDown->Fill();
        newBranch1_m_b2tautau_vis_JetRelativeBalDown->Fill();
        newBranch1_d_ma_JetRelativeBalDown->Fill();
        newBranch1_mtMET_b_JetRelativeSampleUp->Fill();
        newBranch1_m_bb_JetRelativeSampleUp->Fill();
        newBranch1_m_bbtautau_vis_JetRelativeSampleUp->Fill();
        newBranch1_m_b2tautau_vis_JetRelativeSampleUp->Fill();
        newBranch1_d_ma_JetRelativeSampleUp->Fill();
        newBranch1_mtMET_b_JetRelativeSampleDown->Fill();
        newBranch1_m_bb_JetRelativeSampleDown->Fill();
        newBranch1_m_bbtautau_vis_JetRelativeSampleDown->Fill();
        newBranch1_m_b2tautau_vis_JetRelativeSampleDown->Fill();
        newBranch1_d_ma_JetRelativeSampleDown->Fill();
        newBranch1_mtMET_b_JERUp->Fill();
        newBranch1_m_bb_JERUp->Fill();
        newBranch1_m_bbtautau_vis_JERUp->Fill();
        newBranch1_m_b2tautau_vis_JERUp->Fill();
        newBranch1_d_ma_JERUp->Fill();
        newBranch1_mtMET_b_JERDown->Fill();
        newBranch1_m_bb_JERDown->Fill();
        newBranch1_m_bbtautau_vis_JERDown->Fill();
        newBranch1_m_b2tautau_vis_JERDown->Fill();
        newBranch1_d_ma_JERDown->Fill();


      }
      dir->cd();
      t->Write("",TObject::kOverwrite);
      strcpy(TreeToUse,stringA);
      if (std::string(t->GetName()).find("mutau") != std::string::npos) got_mutau_tree = true;
      if (std::string(t->GetName()).find("etau") != std::string::npos) got_etau_tree = true;
      if (std::string(t->GetName()).find("emu") != std::string::npos) got_emu_tree = true;
    }
  }
}

void runBDTeval(float& BDTval, CalculatedFeatures& features_to_save,
  int booster_idx, pat::XGBooster booster, float jet_1_pt, float jet_1_eta, 
  float jet_1_phi, float jet_1_mass, float jet_2_pt, float jet_2_eta, 
  float jet_2_phi, float jet_2_mass, float pt_1, float eta_1, 
  float phi_1, float mass_1, float pt_2, float eta_2, 
  float phi_2, float mass_2, float met, float met_phi, 
  float d_zeta, int njets, float m_btautau_vis, float mtMET_1, 
  float mtMET_2, float pt_vis
) {
  ROOT::Math::PtEtaPhiMVector bjet1(jet_1_pt, jet_1_eta, jet_1_phi, jet_1_mass);
  ROOT::Math::PtEtaPhiMVector bjet2(jet_2_pt, jet_2_eta, jet_2_phi, jet_2_mass);
  ROOT::Math::PtEtaPhiMVector bjet12 = bjet1 + bjet2;
  ROOT::Math::PtEtaPhiMVector tau1(pt_1, eta_1, phi_1, mass_1);
  ROOT::Math::PtEtaPhiMVector tau2(pt_2, eta_2, phi_2, mass_2);
  ROOT::Math::PtEtaPhiMVector tau12 = tau1 + tau2;

  features_to_save.btautau_dR = compute_deltaR(jet_1_eta, tau12.Eta(), jet_1_phi, tau12.Phi());
  features_to_save.btau1_dR = compute_deltaR(jet_1_eta, eta_1, jet_1_phi, phi_1);
  features_to_save.btau2_dR = compute_deltaR(jet_1_eta, eta_2, jet_1_phi, phi_2);
  features_to_save.mtMET_b = compute_mt(jet_1_pt, jet_1_phi, met, met_phi);

  if (booster_idx > 2) {
    features_to_save.m_bb = bjet12.M();
    features_to_save.m_bbtautau_vis = (bjet12 + tau12).M();
    features_to_save.m_b2tautau_vis = (bjet2 + tau12).M();
    features_to_save.d_ma = (bjet12.M() - tau12.M())/tau12.M();
    features_to_save.b2tau1_dR = compute_deltaR(jet_2_eta, eta_1, jet_2_phi, phi_1);
    features_to_save.b2tau2_dR = compute_deltaR(jet_2_eta, eta_2, jet_2_phi, phi_2);
    features_to_save.b2tautau_dR = compute_deltaR(jet_2_eta, tau12.Eta(), jet_2_phi, tau12.Phi());
  }

  if (booster_idx == 0) {
    float in_mT_thMET, in_mT_mMET, in_d_zeta, in_mutau_pt, in_m_b1mt, in_bmt_dR, in_b1th_dR, in_njets;
    in_mT_thMET = mtMET_2;
    in_mT_mMET = mtMET_1;
    in_d_zeta = d_zeta;
    in_mutau_pt = pt_vis;
    in_m_b1mt = m_btautau_vis;
    in_bmt_dR = features_to_save.btautau_dR;
    in_b1th_dR = features_to_save.btau2_dR;
    in_njets = njets;

    booster.set("mT_thMET", in_mT_thMET);
    booster.set("mT_mMET", in_mT_mMET);
    booster.set("d_zeta", in_d_zeta);
    booster.set("mutau_pt", in_mutau_pt);
    booster.set("m_b1mt", in_m_b1mt);
    booster.set("bmt_dR", in_bmt_dR);
    booster.set("b1th_dR", in_b1th_dR);
    booster.set("njets", in_njets);
    BDTval = booster.predict();
  } 
  else if (booster_idx == 1) {
    float in_pt_vis, in_d_zeta, in_b1e_dR, in_m_btautau_vis, in_mtMET_1, in_b1th_dR, in_mtMET_2, in_njets;
    in_pt_vis = pt_vis;
    in_d_zeta = d_zeta;
    in_b1e_dR = features_to_save.btau1_dR;
    in_m_btautau_vis = m_btautau_vis;
    in_mtMET_1 = mtMET_1;
    in_b1th_dR = features_to_save.btau2_dR;
    in_mtMET_2 = mtMET_2;
    in_njets = njets;

    booster.set("pt_vis_nominal", in_pt_vis);
    booster.set("D_zeta_nominal", in_d_zeta);
    booster.set("b1e_dR", in_b1e_dR);
    booster.set("m_btautau_vis_nominal", in_m_btautau_vis);
    booster.set("mtMET_1_nominal", in_mtMET_1);
    booster.set("b1th_dR", in_b1th_dR);
    booster.set("mtMET_2_nominal", in_mtMET_2);
    booster.set("njets", in_njets);
    BDTval = booster.predict();
  } 
  else if (booster_idx == 2) {
    float in_pt_vis, in_d_zeta, in_mT_b1MET, in_pt_1, in_m_btautau_vis, in_mtMET_1, in_b1emu_dR, in_njets;
    in_pt_vis = pt_vis;
    in_d_zeta = d_zeta;
    in_mT_b1MET = features_to_save.mtMET_b;
    in_pt_1 = pt_1;
    in_m_btautau_vis = m_btautau_vis;
    in_mtMET_1 = mtMET_1;
    in_b1emu_dR = features_to_save.btautau_dR;
    in_njets = njets;

    booster.set("pt_vis_nominal", in_pt_vis);
    booster.set("D_zeta_nominal", in_d_zeta);
    booster.set("mT_b1MET", in_mT_b1MET);
    booster.set("pt_1_nominal", in_pt_1);
    booster.set("m_btautau_vis_nominal", in_m_btautau_vis);
    booster.set("mtMET_1_nominal", in_mtMET_1);
    booster.set("b1emu_dR", in_b1emu_dR);
    booster.set("njets", in_njets);
    BDTval = booster.predict();
  } 
  else if (booster_idx == 3) {
    float in_m_b2mt, in_bmt_dR, in_b2th_dR, in_mT_mMET, in_m_bbmt, in_d_ma;
    in_m_b2mt = features_to_save.m_b2tautau_vis;
    in_bmt_dR = features_to_save.btautau_dR;
    in_b2th_dR = features_to_save.b2tau2_dR;
    in_mT_mMET = mtMET_1;
    in_m_bbmt = features_to_save.m_bbtautau_vis;
    in_d_ma = features_to_save.d_ma;

    booster.set("m_b2mt", in_m_b2mt);
    booster.set("bmt_dR", in_bmt_dR);
    booster.set("b2th_dR", in_b2th_dR);
    booster.set("mT_mMET", in_mT_mMET);
    booster.set("m_bbmt", in_m_bbmt);
    booster.set("d_ma", in_d_ma);
    BDTval = booster.predict();
  } 
  else if (booster_idx == 4) {
    float in_mbb, in_m_btautau_vis, in_mtMET_1, in_b1th_dR, in_b1e_dR, in_b2th_dR;
    in_mbb = features_to_save.m_bb;
    in_m_btautau_vis = m_btautau_vis;
    in_mtMET_1 = mtMET_1;
    in_b1th_dR = features_to_save.btau2_dR;
    in_b1e_dR = features_to_save.btau1_dR;
    in_b2th_dR = features_to_save.b2tau1_dR; //BUG in naming at model level

    booster.set("mbb", in_mbb);
    booster.set("m_btautau_vis_nominal", in_m_btautau_vis);
    booster.set("mtMET_1_nominal", in_mtMET_1);
    booster.set("b1th_dR", in_b1th_dR);
    booster.set("b1e_dR", in_b1e_dR);
    booster.set("b2th_dR", in_b2th_dR);
    BDTval = booster.predict();
  }
  else if (booster_idx == 5) {
    float in_pt_vis, in_mT_b1MET, in_mtMET_1, in_b2emu_dR, in_d_ma;
    in_pt_vis = pt_vis;
    in_mT_b1MET = features_to_save.mtMET_b;
    in_mtMET_1 = mtMET_1;
    in_b2emu_dR = features_to_save.b2tautau_dR;
    in_d_ma = features_to_save.d_ma;

    booster.set("pt_vis_nominal", in_pt_vis);
    booster.set("mT_b1MET", in_mT_b1MET);
    booster.set("mtMET_1_nominal", in_mtMET_1);
    booster.set("b2emu_dR", in_b2emu_dR);
    booster.set("d_ma", in_d_ma);
    BDTval = booster.predict();
  } 
  else BDTval = -10.;
}

void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir;
  if(source->GetName()!=parser.stringValue("inputFile")){
    adir = savdir->mkdir(source->GetName());
    std::cout<<"Source name is not outputfile name"<<std::endl;
    adir->cd();
  }
  else{
    adir->cd();
  }

  // loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    std::cout<<"My key is: "<<key->GetName()<<std::endl;
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir,parser);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyFile(const char *fname, optutl::CommandLineParser parser) {
  // copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
}

void copyFiles(optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) {
  // prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();
}

