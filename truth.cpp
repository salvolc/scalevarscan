#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "RConfig.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h" 
#include "TMultiGraph.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include <iostream>
#include "TLorentzVector.h"
#include "math.h"
#include "ExRootClasses.h"
#include "ExRootTreeReader.h"
#include "TLeaf.h"
#include "classes/DelphesClasses.h"


#include "/home/salv/eigen/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using std::ofstream;
using std::cout;
using std::endl;
using std::string;
using std::cin;

int get_nevents(string fileName);
MatrixXd get_eventdisplay(string fileName, int event);
void ladebalken(int i, int max);
void speichere(std::string name, MatrixXd data);
void speichere(std::string name, VectorXd data);
VectorXd get_event_MET(string fileName, int event);
TLorentzVector get_event_MET_neutrino(string fileName, int event);
MatrixXd get_eventdisplay_particle(string fileName, int event, int PID);
int get_numberOfPhotons(string fileName);
int get_numberOfJets(string fileName);
int get_numberOfbJets(string fileName);
TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[], double mass);
TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[]);
int get_nTruth(string fileName, int PID);
VectorXd vertex_match(TClonesArray* TCP, int iEvent);
double get_weight(TBranch* branchE, int iEvent);
int get_highestpt_particle(TClonesArray* TCP,int iPID);
int follow_particle(TClonesArray* TCP, int ID);
bool is_uc(int i){i = abs(i);if (i == 2 || i == 4){return true;}else return false;}
bool is_ucg(int i){i = abs(i);if (i == 2 || i == 4 || i == 21 || i == 9){return true;}else return false;}
bool is_lep(int i){i = abs(i);if (i > 10 && i < 17){return true;}else return false;}
bool is_top(int i){i = abs(i);if (i == 6){return true;}else return false;}
bool is_col(int i){i = abs(i);if (i == 21 || (i > 0 && i < 6)){return true;}else return false;}



int error=0;
bool debug=true;

int main(int argc, char const *argv[])
{


	gSystem->Load("/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/libExRootAnalysis.so");
	gSystem->Load("/home/salv/Dokumente/Masterarbeit/Delphes/libDelphes.so");

	string fileNames[6];
	fileNames[0] = "samples/tua_LH_production_LO_multileg_ta_taj_2_0.1_PY8_halftop_scale_DELCMS_50.root";
	fileNames[1] = "samples/tua_LH_production_LO_multileg_ta_taj_2_0.1_PY8_twotop_scale_DELCMS_50.root";
	fileNames[2] = "samples/tua_LH_production_LO_multileg_ta_taj_2_0.1_PY8_dynscale2_DELCMS_50.root";
	fileNames[3] = "samples/tua_LH_production_LO_multileg_ta_taj_2_0.1_PY8_top_scale_DELCMS_50.root";
	fileNames[4] = "samples/tua_LH_production_LO_multileg_ta_taj_2_0.1_PY8_tripeltop_scale_DELCMS_50.root";
	fileNames[5] = "samples/tua_LH_production_LO_multileg_ta_taj_2_0.1_PY8_owndyn_scale_DELCMS_50.root";
	//fileNames[0] = "/home/salv/fusessh/samplescan/tua_LH_decay_LO_multileg_tt_wbua_2_0.001_PY8_DELATL_50.root";
	//fileNames[1] = "/home/salv/fusessh/samplescan/tua_LH_interference_LO_multileg_tja_ta_2_0.001_PY8_DELATL_50.root";
	//fileNames[2] = "/home/salv/fusessh/samplescan/tua_LH_production_LO_multileg_ta_taj_2_0.001_PY8_DELATL_50.root";

	string fileTPPTNames[6];	fileTPPTNames[0] = "data/htop_Photon_PT_truth";		fileTPPTNames[1] = "data/ttop_Photon_PT_truth";		fileTPPTNames[2] = "data/dyn_Photon_PT_truth";		fileTPPTNames[3] = "data/top_Photon_PT_truth";		fileTPPTNames[4] = "data/tritop_Photon_PT_truth";		fileTPPTNames[5] = "data/owntop_Photon_PT_truth";
	string fileTPEtaNames[6];	fileTPEtaNames[0] = "data/htop_Photon_Eta_truth";	fileTPEtaNames[1] = "data/ttop_Photon_Eta_truth";	fileTPEtaNames[2] = "data/dyn_Photon_Eta_truth";	fileTPEtaNames[3] = "data/top_Photon_Eta_truth";	fileTPEtaNames[4] = "data/tritop_Photon_Eta_truth";		fileTPEtaNames[5] = "data/owntop_Photon_Eta_truth";
	string fileTPPhiNames[6];	fileTPPhiNames[0] = "data/htop_Photon_Phi_truth";	fileTPPhiNames[1] = "data/ttop_Photon_Phi_truth";	fileTPPhiNames[2] = "data/dyn_Photon_Phi_truth";	fileTPPhiNames[3] = "data/top_Photon_Phi_truth";	fileTPPhiNames[4] = "data/tritop_Photon_Phi_truth";		fileTPPhiNames[5] = "data/owntop_Photon_Phi_truth";
	
	string fileTWPTNames[6];	fileTWPTNames[0] = "data/htop_WBoson_PT_truth";		fileTWPTNames[1] = "data/ttop_WBoson_PT_truth";		fileTWPTNames[2] = "data/dyn_WBoson_PT_truth";		fileTWPTNames[3] = "data/top_WBoson_PT_truth";		fileTWPTNames[4] = "data/tritop_WBoson_PT_truth";		fileTWPTNames[5] = "data/owntop_WBoson_PT_truth";
	string fileTWEtaNames[6];	fileTWEtaNames[0] = "data/htop_WBoson_Eta_truth";	fileTWEtaNames[1] = "data/ttop_WBoson_Eta_truth";	fileTWEtaNames[2] = "data/dyn_WBoson_Eta_truth";	fileTWEtaNames[3] = "data/top_WBoson_Eta_truth";	fileTWEtaNames[4] = "data/tritop_WBoson_Eta_truth";		fileTWEtaNames[5] = "data/owntop_WBoson_Eta_truth";
	string fileTWPhiNames[6];	fileTWPhiNames[0] = "data/htop_WBoson_Phi_truth";	fileTWPhiNames[1] = "data/ttop_WBoson_Phi_truth";	fileTWPhiNames[2] = "data/dyn_WBoson_Phi_truth";	fileTWPhiNames[3] = "data/top_WBoson_Phi_truth";	fileTWPhiNames[4] = "data/tritop_WBoson_Phi_truth";		fileTWPhiNames[5] = "data/owntop_WBoson_Phi_truth";
	string fileTWMNames[6];		fileTWMNames[0] = "data/htop_WBoson_M_truth";		fileTWMNames[1] = "data/ttop_WBoson_M_truth";		fileTWMNames[2] = "data/dyn_WBoson_M_truth";		fileTWMNames[3] = "data/top_WBoson_M_truth";		fileTWMNames[4] = "data/tritop_WBoson_M_truth";			fileTWMNames[5] = "data/owntop_WBoson_M_truth";
	
	string fileTtPTNames[6];	fileTtPTNames[0] = "data/htop_TopQuark_PT_truth";	fileTtPTNames[1] = "data/ttop_TopQuark_PT_truth";	fileTtPTNames[2] = "data/dyn_TopQuark_PT_truth";	fileTtPTNames[3] = "data/top_TopQuark_PT_truth";	fileTtPTNames[4] = "data/tritop_TopQuark_PT_truth";		fileTtPTNames[5] = "data/owntop_TopQuark_PT_truth";
	string fileTtEtaNames[6];	fileTtEtaNames[0] = "data/htop_TopQuark_Eta_truth";	fileTtEtaNames[1] = "data/ttop_TopQuark_Eta_truth";	fileTtEtaNames[2] = "data/dyn_TopQuark_Eta_truth";	fileTtEtaNames[3] = "data/top_TopQuark_Eta_truth";	fileTtEtaNames[4] = "data/tritop_TopQuark_Eta_truth";	fileTtEtaNames[5] = "data/owntop_TopQuark_Eta_truth";
	string fileTtPhiNames[6];	fileTtPhiNames[0] = "data/htop_TopQuark_Phi_truth";	fileTtPhiNames[1] = "data/ttop_TopQuark_Phi_truth";	fileTtPhiNames[2] = "data/dyn_TopQuark_Phi_truth";	fileTtPhiNames[3] = "data/top_TopQuark_Phi_truth";	fileTtPhiNames[4] = "data/tritop_TopQuark_Phi_truth";	fileTtPhiNames[5] = "data/owntop_TopQuark_Phi_truth";
	string fileTtMNames[6];		fileTtMNames[0] = "data/htop_TopQuark_M_truth";		fileTtMNames[1] = "data/ttop_TopQuark_M_truth";		fileTtMNames[2] = "data/dyn_TopQuark_M_truth";		fileTtMNames[3] = "data/top_TopQuark_M_truth";		fileTtMNames[4] = "data/tritop_TopQuark_M_truth";		fileTtMNames[5] = "data/owntop_TopQuark_M_truth";
	
	string fileTbPTNames[6];	fileTbPTNames[0] = "data/htop_BQuark_PT_truth";		fileTbPTNames[1] = "data/ttop_BQuark_PT_truth";		fileTbPTNames[2] = "data/dyn_BQuark_PT_truth";		fileTbPTNames[3] = "data/top_BQuark_PT_truth";		fileTbPTNames[4] = "data/tritop_BQuark_PT_truth";		fileTbPTNames[5] = "data/owntop_BQuark_PT_truth";
	string fileTbEtaNames[6];	fileTbEtaNames[0] = "data/htop_BQuark_Eta_truth";	fileTbEtaNames[1] = "data/ttop_BQuark_Eta_truth";	fileTbEtaNames[2] = "data/dyn_BQuark_Eta_truth";	fileTbEtaNames[3] = "data/top_BQuark_Eta_truth";	fileTbEtaNames[4] = "data/tritop_BQuark_Eta_truth";		fileTbEtaNames[5] = "data/owntop_BQuark_Eta_truth";
	string fileTbPhiNames[6];	fileTbPhiNames[0] = "data/htop_BQuark_Phi_truth";	fileTbPhiNames[1] = "data/ttop_BQuark_Phi_truth";	fileTbPhiNames[2] = "data/dyn_BQuark_Phi_truth";	fileTbPhiNames[3] = "data/top_BQuark_Phi_truth";	fileTbPhiNames[4] = "data/tritop_BQuark_Phi_truth";		fileTbPhiNames[5] = "data/owntop_BQuark_Phi_truth";
	string fileTbMNames[6];		fileTbMNames[0] = "data/htop_BQuark_M_truth";		fileTbMNames[1] = "data/ttop_BQuark_M_truth";		fileTbMNames[2] = "data/dyn_BQuark_M_truth";		fileTbMNames[3] = "data/top_BQuark_M_truth";		fileTbMNames[4] = "data/tritop_BQuark_M_truth";			fileTbMNames[5] = "data/owntop_BQuark_M_truth";
	
	string fileTuPTNames[6];	fileTuPTNames[0] = "data/htop_UQuark_PT_truth";		fileTuPTNames[1] = "data/ttop_UQuark_PT_truth";		fileTuPTNames[2] = "data/dyn_UQuark_PT_truth";		fileTuPTNames[3] = "data/top_UQuark_PT_truth";		fileTuPTNames[4] = "data/tritop_UQuark_PT_truth";		fileTuPTNames[5] = "data/owntop_UQuark_PT_truth";
	string fileTuEtaNames[6];	fileTuEtaNames[0] = "data/htop_UQuark_Eta_truth";	fileTuEtaNames[1] = "data/ttop_UQuark_Eta_truth";	fileTuEtaNames[2] = "data/dyn_UQuark_Eta_truth";	fileTuEtaNames[3] = "data/top_UQuark_Eta_truth";	fileTuEtaNames[4] = "data/tritop_UQuark_Eta_truth";		fileTuEtaNames[5] = "data/owntop_UQuark_Eta_truth";
	string fileTuPhiNames[6];	fileTuPhiNames[0] = "data/htop_UQuark_Phi_truth";	fileTuPhiNames[1] = "data/ttop_UQuark_Phi_truth";	fileTuPhiNames[2] = "data/dyn_UQuark_Phi_truth";	fileTuPhiNames[3] = "data/top_UQuark_Phi_truth";	fileTuPhiNames[4] = "data/tritop_UQuark_Phi_truth";		fileTuPhiNames[5] = "data/owntop_UQuark_Phi_truth";
	string fileTuMNames[6];		fileTuMNames[0] = "data/htop_UQuark_M_truth";		fileTuMNames[1] = "data/ttop_UQuark_M_truth";		fileTuMNames[2] = "data/dyn_UQuark_M_truth";		fileTuMNames[3] = "data/top_UQuark_M_truth";		fileTuMNames[4] = "data/tritop_UQuark_M_truth";			fileTuMNames[5] = "data/owntop_UQuark_M_truth";
	
	string fileWeightNames[6];	fileWeightNames[0] = "data/htop_weight_truth";		fileWeightNames[1] = "data/ttop_weight_truth";		fileWeightNames[2] = "data/dyn_weight_truth";		fileWeightNames[3] = "data/top_weight_truth";		fileWeightNames[4] = "data/tritop_weight_truth";		fileWeightNames[5] = "data/owntop_weight_truth";

	string fileTtPRNames[6];	fileTtPRNames[0] = "data/htop_TopQuark_Photon_R_truth";	fileTtPRNames[1] = "data/ttop_TopQuark_Photon_R_truth";	fileTtPRNames[2] = "data/dyn_TopQuark_Photon_R_truth";	fileTtPRNames[3] = "data/top_TopQuark_Photon_R_truth";	fileTtPRNames[4] = "data/tritop_TopQuark_Photon_R_truth";	fileTtPRNames[5] = "data/owntop_TopQuark_Photon_R_truth";
	string fileTtbRNames[6];	fileTtbRNames[0] = "data/htop_TopQuark_BQuark_R_truth";	fileTtbRNames[1] = "data/ttop_TopQuark_BQuark_R_truth";	fileTtbRNames[2] = "data/dyn_TopQuark_BQuark_R_truth";	fileTtbRNames[3] = "data/top_TopQuark_BQuark_R_truth";	fileTtbRNames[4] = "data/tritop_TopQuark_BQuark_R_truth";	fileTtbRNames[5] = "data/owntop_TopQuark_BQuark_R_truth";
	string fileTtWRNames[6];	fileTtWRNames[0] = "data/htop_TopQuark_WBoson_R_truth";	fileTtWRNames[1] = "data/ttop_TopQuark_WBoson_R_truth";	fileTtWRNames[2] = "data/dyn_TopQuark_WBoson_R_truth";	fileTtWRNames[3] = "data/top_TopQuark_WBoson_R_truth";	fileTtWRNames[4] = "data/tritop_TopQuark_WBoson_R_truth";	fileTtWRNames[5] = "data/owntop_TopQuark_WBoson_R_truth";
	string fileTtuRNames[6];	fileTtuRNames[0] = "data/htop_TopQuark_UQuark_R_truth";	fileTtuRNames[1] = "data/ttop_TopQuark_UQuark_R_truth";	fileTtuRNames[2] = "data/dyn_TopQuark_UQuark_R_truth";	fileTtuRNames[3] = "data/top_TopQuark_UQuark_R_truth";	fileTtuRNames[4] = "data/tritop_TopQuark_UQuark_R_truth";	fileTtuRNames[5] = "data/owntop_TopQuark_UQuark_R_truth";

	string fileTtPMNames[6];	fileTtPMNames[0] = "data/htop_TopQuark_Photon_M_truth";	fileTtPMNames[1] = "data/ttop_TopQuark_Photon_M_truth";	fileTtPMNames[2] = "data/dyn_TopQuark_Photon_M_truth";	fileTtPMNames[3] = "data/top_TopQuark_Photon_M_truth";	fileTtPMNames[4] = "data/tritop_TopQuark_Photon_M_truth";	fileTtPMNames[5] = "data/owntop_TopQuark_Photon_M_truth";
	string fileTtbMNames[6];	fileTtbMNames[0] = "data/htop_TopQuark_BQuark_M_truth";	fileTtbMNames[1] = "data/ttop_TopQuark_BQuark_M_truth";	fileTtbMNames[2] = "data/dyn_TopQuark_BQuark_M_truth";	fileTtbMNames[3] = "data/top_TopQuark_BQuark_M_truth";	fileTtbMNames[4] = "data/tritop_TopQuark_BQuark_M_truth";	fileTtbMNames[5] = "data/owntop_TopQuark_BQuark_M_truth";
	string fileTtWMNames[6];	fileTtWMNames[0] = "data/htop_TopQuark_WBoson_M_truth";	fileTtWMNames[1] = "data/ttop_TopQuark_WBoson_M_truth";	fileTtWMNames[2] = "data/dyn_TopQuark_WBoson_M_truth";	fileTtWMNames[3] = "data/top_TopQuark_WBoson_M_truth";	fileTtWMNames[4] = "data/tritop_TopQuark_WBoson_M_truth";	fileTtWMNames[5] = "data/owntop_TopQuark_WBoson_M_truth";
	string fileTtuMNames[6];	fileTtuMNames[0] = "data/htop_TopQuark_UQuark_M_truth";	fileTtuMNames[1] = "data/ttop_TopQuark_UQuark_M_truth";	fileTtuMNames[2] = "data/dyn_TopQuark_UQuark_M_truth";	fileTtuMNames[3] = "data/top_TopQuark_UQuark_M_truth";	fileTtuMNames[4] = "data/tritop_TopQuark_UQuark_M_truth";	fileTtuMNames[5] = "data/owntop_TopQuark_UQuark_M_truth";

	string fileTPtRNames[6];	fileTPtRNames[0] = "data/htop_Photon_TopQuark_R_truth";	fileTPtRNames[1] = "data/ttop_Photon_TopQuark_R_truth";	fileTPtRNames[2] = "data/dyn_Photon_TopQuark_R_truth";	fileTPtRNames[3] = "data/top_Photon_TopQuark_R_truth";	fileTPtRNames[4] = "data/tritop_Photon_TopQuark_R_truth";	fileTPtRNames[5] = "data/owntop_Photon_TopQuark_R_truth";
	string fileTPbRNames[6];	fileTPbRNames[0] = "data/htop_Photon_BQuark_R_truth";	fileTPbRNames[1] = "data/ttop_Photon_BQuark_R_truth";	fileTPbRNames[2] = "data/dyn_Photon_BQuark_R_truth";	fileTPbRNames[3] = "data/top_Photon_BQuark_R_truth";	fileTPbRNames[4] = "data/tritop_Photon_BQuark_R_truth";		fileTPbRNames[5] = "data/owntop_Photon_BQuark_R_truth";
	string fileTPWRNames[6];	fileTPWRNames[0] = "data/htop_Photon_WBoson_R_truth";	fileTPWRNames[1] = "data/ttop_Photon_WBoson_R_truth";	fileTPWRNames[2] = "data/dyn_Photon_WBoson_R_truth";	fileTPWRNames[3] = "data/top_Photon_WBoson_R_truth";	fileTPWRNames[4] = "data/tritop_Photon_WBoson_R_truth";		fileTPWRNames[5] = "data/owntop_Photon_WBoson_R_truth";
	string fileTPuRNames[6];	fileTPuRNames[0] = "data/htop_Photon_UQuark_R_truth";	fileTPuRNames[1] = "data/ttop_Photon_UQuark_R_truth";	fileTPuRNames[2] = "data/dyn_Photon_UQuark_R_truth";	fileTPuRNames[3] = "data/top_Photon_UQuark_R_truth";	fileTPuRNames[4] = "data/tritop_Photon_UQuark_R_truth";		fileTPuRNames[5] = "data/owntop_Photon_UQuark_R_truth";
	
	string fileTPtMNames[6];	fileTPtMNames[0] = "data/htop_Photon_TopQuark_M_truth";	fileTPtMNames[1] = "data/ttop_Photon_TopQuark_M_truth";	fileTPtMNames[2] = "data/dyn_Photon_TopQuark_M_truth";	fileTPtMNames[3] = "data/top_Photon_TopQuark_M_truth";	fileTPtMNames[4] = "data/tritop_Photon_TopQuark_M_truth";	fileTPtMNames[5] = "data/owntop_Photon_TopQuark_M_truth";
	string fileTPbMNames[6];	fileTPbMNames[0] = "data/htop_Photon_BQuark_M_truth";	fileTPbMNames[1] = "data/ttop_Photon_BQuark_M_truth";	fileTPbMNames[2] = "data/dyn_Photon_BQuark_M_truth";	fileTPbMNames[3] = "data/top_Photon_BQuark_M_truth";	fileTPbMNames[4] = "data/tritop_Photon_BQuark_M_truth";		fileTPbMNames[5] = "data/owntop_Photon_BQuark_M_truth";
	string fileTPWMNames[6];	fileTPWMNames[0] = "data/htop_Photon_WBoson_M_truth";	fileTPWMNames[1] = "data/ttop_Photon_WBoson_M_truth";	fileTPWMNames[2] = "data/dyn_Photon_WBoson_M_truth";	fileTPWMNames[3] = "data/top_Photon_WBoson_M_truth";	fileTPWMNames[4] = "data/tritop_Photon_WBoson_M_truth";		fileTPWMNames[5] = "data/owntop_Photon_WBoson_M_truth";
	string fileTPuMNames[6];	fileTPuMNames[0] = "data/htop_Photon_UQuark_M_truth";	fileTPuMNames[1] = "data/ttop_Photon_UQuark_M_truth";	fileTPuMNames[2] = "data/dyn_Photon_UQuark_M_truth";	fileTPuMNames[3] = "data/top_Photon_UQuark_M_truth";	fileTPuMNames[4] = "data/tritop_Photon_UQuark_M_truth";		fileTPuMNames[5] = "data/owntop_Photon_UQuark_M_truth";


	TFile* files[6];

	for (int iFile = 5; iFile < 6; ++iFile)
	{
		files[iFile] = new TFile(fileNames[iFile].c_str(),"READ");
		TTree* tree = (TTree*)files[iFile]->Get("Delphes");

		TBranch *bP 		= tree->GetBranch("Particle");
		TBranch *bE 		= tree->GetBranch("Event");
		TClonesArray *TCP 	= 0;
		TClonesArray *TCE 	= 0;
		bP->SetAddress(&TCP);
		bE->SetAddress(&TCE);
		bP->GetEntry(0);

		int nEvents = get_nevents(fileNames[iFile].c_str());


		VectorXd VTPhotonPT = VectorXd::Zero(nEvents);VectorXd VTPhotonEta = VectorXd::Zero(nEvents);VectorXd VTPhotonPhi = VectorXd::Zero(nEvents);
		VectorXd VTTopQuarkPT = VectorXd::Zero(nEvents);VectorXd VTTopQuarkEta = VectorXd::Zero(nEvents);VectorXd VTTopQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTTopQuarkM = VectorXd::Zero(nEvents);
		VectorXd VTWBosonPT = VectorXd::Zero(nEvents);VectorXd VTWBosonEta = VectorXd::Zero(nEvents);VectorXd VTWBosonPhi = VectorXd::Zero(nEvents);VectorXd VTWBosonM = VectorXd::Zero(nEvents);
		VectorXd VTBQuarkPT = VectorXd::Zero(nEvents);VectorXd VTBQuarkEta = VectorXd::Zero(nEvents);VectorXd VTBQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTBQuarkM = VectorXd::Zero(nEvents);
		VectorXd VTUQuarkPT = VectorXd::Zero(nEvents);VectorXd VTUQuarkEta = VectorXd::Zero(nEvents);VectorXd VTUQuarkPhi = VectorXd::Zero(nEvents);VectorXd VTUQuarkM = VectorXd::Zero(nEvents);

		VectorXd VWeight = VectorXd::Zero(nEvents);

		VectorXd VTTopQuark_Photon_R = VectorXd::Zero(nEvents);VectorXd VTTopQuark_bQuark_R = VectorXd::Zero(nEvents);VectorXd VTTopQuark_WBoson_R = VectorXd::Zero(nEvents);VectorXd VTTopQuark_UQuark_R = VectorXd::Zero(nEvents);
		VectorXd VTTopQuark_Photon_M = VectorXd::Zero(nEvents);VectorXd VTTopQuark_bQuark_M = VectorXd::Zero(nEvents);VectorXd VTTopQuark_WBoson_M = VectorXd::Zero(nEvents);VectorXd VTTopQuark_UQuark_M = VectorXd::Zero(nEvents);

		VectorXd VTPhoton_TopQuark_R = VectorXd::Zero(nEvents);VectorXd VTPhoton_bQuark_R = VectorXd::Zero(nEvents);VectorXd VTPhoton_WBoson_R = VectorXd::Zero(nEvents);VectorXd VTPhoton_UQuark_R = VectorXd::Zero(nEvents);
		VectorXd VTPhoton_TopQuark_M = VectorXd::Zero(nEvents);VectorXd VTPhoton_bQuark_M = VectorXd::Zero(nEvents);VectorXd VTPhoton_WBoson_M = VectorXd::Zero(nEvents);VectorXd VTPhoton_UQuark_M = VectorXd::Zero(nEvents);

		int i_unmatched = 0;
		//nEvents = 10;
		for (int iEvent = 0; iEvent < nEvents; ++iEvent)
		{
			//cout << get_eventdisplay(fileNames[iFile],iEvent).block(0,0,10,10) << endl;
			//cout << get_eventdisplay_particle(fileNames[iFile],iEvent,22).block(0,0,10,10) << endl;
			//cout << get_eventdisplay_particle(fileNames[iFile],iEvent,6) << endl;
			//cout << get_eventdisplay_particle(fileNames[iFile],iEvent,-6) << endl;

			bP->GetEntry(iEvent);
			bE->GetEntry(iEvent);

			//cout << "\n" << vertex_match(TCP,iFile) << endl;
			//if(debug){cout<<"Input Number for next Event"<<endl;string a;cin >> a;}

			HepMCEvent *E_Event = (HepMCEvent*)TCE->At(0);
			VWeight(iEvent) = E_Event->Weight;
			if (E_Event->Weight == 0){continue;}

			int nParticles = TCP->GetEntries();

			VectorXd i_matched = VectorXd::Zero(10); i_matched = vertex_match(TCP,iFile);

			if(i_matched(9)==0){i_unmatched++;continue;}

			GenParticle* P_Particle = (GenParticle*)TCP->At(i_matched(1));
			TLorentzVector 	TVPhoton(P_Particle->Px,P_Particle->Py,P_Particle->Pz,P_Particle->E);
			VTPhotonPT(iEvent) = P_Particle->PT;
			VTPhotonEta(iEvent) = P_Particle->Eta;
			VTPhotonPhi(iEvent) = P_Particle->Phi;

			P_Particle = (GenParticle*)TCP->At(i_matched(3));
			TLorentzVector	TVTopQuark(P_Particle->Px,P_Particle->Py,P_Particle->Pz,P_Particle->E);
			VTTopQuarkPT(iEvent) = P_Particle->PT;
			VTTopQuarkEta(iEvent) = P_Particle->Eta;
			VTTopQuarkPhi(iEvent) = P_Particle->Phi;
			VTTopQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(get_highestpt_particle(TCP,2));
			TLorentzVector	TVuQuark(P_Particle->Px,P_Particle->Py,P_Particle->Pz,P_Particle->E);
			VTUQuarkPT(iEvent) = P_Particle->PT;
			VTUQuarkEta(iEvent) = P_Particle->Eta;
			VTUQuarkPhi(iEvent) = P_Particle->Phi;
			VTUQuarkM(iEvent) = P_Particle->Mass;
	
			P_Particle = (GenParticle*)TCP->At(i_matched(4));
			TLorentzVector	TVbQuark(P_Particle->Px,P_Particle->Py,P_Particle->Pz,P_Particle->E);
			VTBQuarkPT(iEvent) = P_Particle->PT;
			VTBQuarkEta(iEvent) = P_Particle->Eta;
			VTBQuarkPhi(iEvent) = P_Particle->Phi;
			VTBQuarkM(iEvent) = P_Particle->Mass;

			P_Particle = (GenParticle*)TCP->At(i_matched(5));
			TLorentzVector	TVWBoson(P_Particle->Px,P_Particle->Py,P_Particle->Pz,P_Particle->E);
			VTWBosonPT(iEvent) = P_Particle->PT;
			VTWBosonEta(iEvent) = P_Particle->Eta;
			VTWBosonPhi(iEvent) = P_Particle->Phi;
			VTWBosonM(iEvent) = P_Particle->Mass;

			VTTopQuark_Photon_R(iEvent) = TVTopQuark.DeltaR(TVPhoton);
			VTTopQuark_bQuark_R(iEvent) = TVTopQuark.DeltaR(TVbQuark);
			VTTopQuark_WBoson_R(iEvent) = TVTopQuark.DeltaR(TVWBoson);
			VTTopQuark_UQuark_R(iEvent) = TVTopQuark.DeltaR(TVuQuark);

			VTTopQuark_Photon_M(iEvent) = (TVTopQuark+TVPhoton).M();
			VTTopQuark_bQuark_M(iEvent) = (TVTopQuark+TVbQuark).M();
			VTTopQuark_WBoson_M(iEvent) = (TVTopQuark+TVWBoson).M();
			VTTopQuark_UQuark_M(iEvent) = (TVTopQuark+TVuQuark).M();

			VTPhoton_TopQuark_R(iEvent) = TVPhoton.DeltaR(TVTopQuark);
			VTPhoton_bQuark_R(iEvent) = TVPhoton.DeltaR(TVbQuark);
			VTPhoton_WBoson_R(iEvent) = TVPhoton.DeltaR(TVWBoson);
			VTPhoton_UQuark_R(iEvent) = TVPhoton.DeltaR(TVuQuark);

			VTPhoton_TopQuark_M(iEvent) = (TVPhoton+TVTopQuark).M();
			VTPhoton_bQuark_M(iEvent) = (TVPhoton+TVbQuark).M();
			VTPhoton_WBoson_M(iEvent) = (TVPhoton+TVWBoson).M();
			VTPhoton_UQuark_M(iEvent) = (TVPhoton+TVuQuark).M();


			//if(iEvent%10000==0){ladebalken(iEvent*(iFile+1),nEvents*3);}
			//if(debug)cout << "Event " << iEvent << " Photon PT " << VTPhotonPT(iEvent) << " TopQuarkPT " << VTTopQuarkPT(iEvent) << endl;
			//if(debug)cout << get_eventdisplay_particle(fileNames[iFile],iEvent,22).block(0,0,5,12)<< endl;
			//if(debug)cout << get_eventdisplay_particle(fileNames[iFile],iEvent,6)<< endl;
			//if(debug)cout << get_eventdisplay_particle(fileNames[iFile],iEvent,-6)<< endl;
		}

		speichere(fileTPPTNames[iFile],VTPhotonPT);
		speichere(fileTPEtaNames[iFile],VTPhotonEta);
		speichere(fileTPPhiNames[iFile],VTPhotonPhi);

		speichere(fileTtPTNames[iFile],VTTopQuarkPT);
		speichere(fileTtEtaNames[iFile],VTTopQuarkEta);
		speichere(fileTtPhiNames[iFile],VTTopQuarkPhi);
		speichere(fileTtMNames[iFile],VTTopQuarkM);

		speichere(fileTbPTNames[iFile],VTBQuarkPT);
		speichere(fileTbEtaNames[iFile],VTBQuarkEta);
		speichere(fileTbPhiNames[iFile],VTBQuarkPhi);
		speichere(fileTbMNames[iFile],VTBQuarkM);

		speichere(fileTuPTNames[iFile],VTBQuarkPT);
		speichere(fileTuEtaNames[iFile],VTBQuarkEta);
		speichere(fileTuPhiNames[iFile],VTBQuarkPhi);
		speichere(fileTuMNames[iFile],VTBQuarkM);

		speichere(fileTWPTNames[iFile],VTWBosonPT);
		speichere(fileTWEtaNames[iFile],VTWBosonEta);
		speichere(fileTWPhiNames[iFile],VTWBosonPhi);
		speichere(fileTWMNames[iFile],VTWBosonM);

		speichere(fileTtPRNames[iFile],VTTopQuark_Photon_R);
		speichere(fileTtbRNames[iFile],VTTopQuark_bQuark_R);
		speichere(fileTtWRNames[iFile],VTTopQuark_WBoson_R);
		speichere(fileTtuRNames[iFile],VTTopQuark_UQuark_R);

		speichere(fileTtPMNames[iFile],VTTopQuark_Photon_M);
		speichere(fileTtbMNames[iFile],VTTopQuark_bQuark_M);
		speichere(fileTtWMNames[iFile],VTTopQuark_WBoson_M);
		speichere(fileTtuMNames[iFile],VTTopQuark_UQuark_M);

		speichere(fileTPtRNames[iFile],VTPhoton_TopQuark_R);
		speichere(fileTPbRNames[iFile],VTPhoton_bQuark_R);
		speichere(fileTPWRNames[iFile],VTPhoton_WBoson_R);
		speichere(fileTPuRNames[iFile],VTPhoton_UQuark_R);

		speichere(fileTPtMNames[iFile],VTPhoton_TopQuark_M);
		speichere(fileTPbMNames[iFile],VTPhoton_bQuark_M);
		speichere(fileTPWMNames[iFile],VTPhoton_WBoson_M);
		speichere(fileTPuMNames[iFile],VTPhoton_UQuark_M);


		cout << i_unmatched << endl;

		files[iFile]->Close();
		//file->Close();
	}

	return 0;

}

VectorXd vertex_match(TClonesArray* TCP,int iFile){
	int nPart = TCP->GetEntries();
	VectorXd i_vertex_matched = VectorXd::Zero(10);
	int nh = 999999;

	for (int i_p = 0; i_p < nPart; ++i_p)
	{
		GenParticle* P = (GenParticle*)TCP->At(i_p);
		GenParticle* PH1 = (GenParticle*)TCP->At(i_p);
		GenParticle* PH2 = (GenParticle*)TCP->At(i_p);

		if (abs(P->PID)==6 && i_vertex_matched(8)==0)
		{
			nh = follow_particle(TCP,i_p);
			P = (GenParticle*)TCP->At(nh);
			GenParticle* PH1 = (GenParticle*)TCP->At(P->D1);
			GenParticle* PH2 = (GenParticle*)TCP->At(P->D2);
			if((abs(PH1->PID)==22 && is_uc(PH2->PID)) || (abs(PH2->PID)==22 && is_uc(PH1->PID)))
			{
				i_vertex_matched(0) = nh;
				i_vertex_matched(8) = 1;
			}
		}

		if (abs(P->PID)==6 && i_vertex_matched(7)==0)
		{
			nh = follow_particle(TCP,i_p);
			P = (GenParticle*)TCP->At(nh);
			GenParticle* PH1 = (GenParticle*)TCP->At(P->D1);
			GenParticle* PH2 = (GenParticle*)TCP->At(P->D2);
			if((abs(PH1->PID)==24 && abs(PH2->PID)==5))
			{
				i_vertex_matched(3) = nh;
				i_vertex_matched(4) = follow_particle(TCP,P->D2);
				i_vertex_matched(5) = follow_particle(TCP,P->D1);
				i_vertex_matched(7) = 1;
			}
			if((abs(PH2->PID)==24 && abs(PH1->PID)==5))
			{
				i_vertex_matched(3) = nh;
				i_vertex_matched(4) = follow_particle(TCP,P->D1);
				i_vertex_matched(5) = follow_particle(TCP,P->D2);
				i_vertex_matched(7) = 1;
			}
		}

		if(abs(P->PID)==22 && i_vertex_matched(6)==0)
		{
			P = (GenParticle*)TCP->At(i_p);
			nh = follow_particle(TCP,i_p);
			//No Mother no FCNC
			if (P->M1 == -1 && P->M2 == -1){continue;
			} else if ((P->M1 == -1 || P->M2 == -1) && iFile != 2)
			{	
				//One mother = Top then FCNC
				if (P->M1 != -1)
					{
						GenParticle* PH3 = (GenParticle*)TCP->At(P->M1);
						if (is_top(PH3->PID))
						{
							i_vertex_matched(1) = nh;
							i_vertex_matched(6) = 1;
						}
					}
				if (P->M2 != -1)
					{
						GenParticle* PH4 = (GenParticle*)TCP->At(P->M2);
						if (is_top(PH4->PID))
						{
							i_vertex_matched(1) = nh;
							i_vertex_matched(6) = 1;
						}
					}
			} else if(iFile != 100 && P->M1 != -1 && P->M2 != -1){
				//Or two mothers bsp. 21 2 means in primary vertex (FCNC) dont do for decay
				GenParticle* PH3 = (GenParticle*)TCP->At(P->M1);
				GenParticle* PH4 = (GenParticle*)TCP->At(P->M2);
				if (is_col(PH3->PID) && is_col(PH4->PID))
				{
					i_vertex_matched(1) = nh;
					i_vertex_matched(6) = 1;
				}
			}
		}
	}
	
	if (i_vertex_matched(6) == 1 && i_vertex_matched(7) == 1) {i_vertex_matched(9)=1;}
	
	return i_vertex_matched;
}

/*VectorXd vertex_match(TClonesArray* TCP,int iFile){
	return VectorXd::Zero(12);
}*/


MatrixXd get_eventdisplay(string fileName, int event){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");

	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();
	MatrixXd display = MatrixXd::Zero(numberOfParticles,10);


    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{

		GenParticle *P = (GenParticle*)TCP->At(ipart);
		GenParticle *Ph = (GenParticle*)TCP->At(ipart);
		display(ipart,0)= ipart;
		display(ipart,1)=P->PID;
		//display(ipart,1)=P->Px;
		//display(ipart,2)=P->Py;
		//display(ipart,3)=P->Pz;
		if (P->D1 == -1){
			display(ipart,2)= 0;
		}else{	
			Ph = (GenParticle*)TCP->At(P->D1);
			display(ipart,2)= Ph->PID;}

		if (P->D2 == -1){
			display(ipart,3)= 0;
		}else{	
			Ph = (GenParticle*)TCP->At(P->D2);
			display(ipart,3)= Ph->PID;}

		display(ipart,4)= P->D1;
		display(ipart,5)= P->D2;
		if (P->M1 == -1){
			display(ipart,6)= 0;
		}else{	
			Ph = (GenParticle*)TCP->At(P->M1);
			display(ipart,6)= Ph->PID;}

		if (P->M2 == -1){
			display(ipart,7)= 0;
		}else{	
			Ph = (GenParticle*)TCP->At(P->M2);
			display(ipart,7)= Ph->PID;}

		display(ipart,8)= P->M1;
		display(ipart,9)= P->M2;
	}

	file->Close();
	return display;
}

MatrixXd get_eventdisplay_particle(string fileName, int event, int PID){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	int nPIDs = 0;
	for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if (P->PID == PID)nPIDs++;
	}
	if(nPIDs == 0){cout << "Kein Teilchen gefunden!" << endl; return MatrixXd::Zero(1,1);}
	MatrixXd display = MatrixXd::Zero(nPIDs,10);
	nPIDs = 0;
	GenParticle *Ph;
    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if (P->PID == PID){
			display(nPIDs,0)=ipart;
			display(nPIDs,1)=P->PID;
			if (P->D1 == -1){
				display(nPIDs,2)= 0;
			}else{	
				Ph = (GenParticle*)TCP->At(P->D1);
				display(nPIDs,2)= Ph->PID;}

			if (P->D2 == -1){
				display(nPIDs,3)= 0;
			}else{	
				Ph = (GenParticle*)TCP->At(P->D2);
				display(nPIDs,3)= Ph->PID;}

			display(nPIDs,4)= P->D1;
			display(nPIDs,5)= P->D2;
			if (P->M1 == -1){
				display(nPIDs,6)= 0;
			}else{	
				Ph = (GenParticle*)TCP->At(P->M1);
				display(nPIDs,6)= Ph->PID;}

			if (P->M2 == -1){
				display(nPIDs,7)= 0;
			}else{	
				Ph = (GenParticle*)TCP->At(P->M2);
				display(nPIDs,7)= Ph->PID;}

			display(nPIDs,8)= P->M1;
			display(nPIDs,9)= P->M2;

		/*display(nPIDs,0)=P->PID;
		display(nPIDs,1)=P->Px;
		display(nPIDs,2)=P->Py;
		display(nPIDs,3)=P->Pz;
		display(nPIDs,4)=P->PT;
		display(nPIDs,5)=P->Mass;
		display(nPIDs,6)=P->E;
		display(nPIDs,7)=P->Eta;
		display(nPIDs,8)=P->Phi;
		display(nPIDs,9)=ipart;
		display(nPIDs,10)=P->D1;
		display(nPIDs,11)=P->D2;*/
		nPIDs++;
		}
	}

	file->Close();
	return display;
}

int follow_particle(TClonesArray* TCP, int ID){
	GenParticle* P = (GenParticle*)TCP->At(ID);
	int part_num = ID;
	int myPID = P->PID;
	if(P->D1 == -1){return part_num;}//cerr << "Teilchen im Endzustand!" << endl;
	if(P->D2 == -1){return part_num;}//cerr << "Teilchen im Endzustand!" << endl;
	GenParticle* PD1 = (GenParticle*)TCP->At(P->D1);
	GenParticle* PD2 = (GenParticle*)TCP->At(P->D2);
	do
	{
		if(P->D1 == -1){return part_num;}//cerr << "Teilchen im Endzustand!" << endl;
		if(P->D2 == -1){return part_num;}//cerr << "Teilchen im Endzustand!" << endl;
		PD1 = (GenParticle*)TCP->At(P->D1);
		PD2 = (GenParticle*)TCP->At(P->D2);
		if (PD1->PID == P->PID && PD2->PID == P->PID){
			if (P->D1 == P->D2){part_num = P->D1;}
			if (P->D1 != P->D2 && !is_ucg(P->PID)){cout << "Gleiche Teilchentyp mit unterschiedlichen indizies, breche ab.." << endl;return part_num;}
		}
		if(PD1->PID == P->PID){
			part_num = P->D1;
		}else if(PD2->PID == P->PID){
			part_num = P->D2;
		}
		P = (GenParticle*)TCP->At(part_num);
	}while(PD1->PID == P->PID || PD2->PID == P->PID);

	return part_num;
}

int get_highestpt_particle(TClonesArray* TCP,int iPID){
	int numberOfParticles = TCP->GetEntries();

	int particle_number = -99999;
	double pt = 0;
    for (int ipart = 0; ipart < numberOfParticles; ++ipart)
	{
	   	GenParticle *P = (GenParticle*)TCP->At(ipart);
		if ((abs(P->PID) == iPID) && (P->PT > pt))
		{pt = P->PT; particle_number = ipart;}
	}
	if (particle_number == -99999){cout << "Warning highest pt light quark not found!!" << endl;}
	return particle_number;
}

int get_nTruth(string fileName, int PID){
  	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(0);  
	int nEvents = branchP->GetEntries();

	int nPIDs = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		branchP->GetEntry(iEvent);  
		int nParticles = TCP->GetEntries();
		for (int ipart = 0; ipart < nParticles; ++ipart)
		{
		   	GenParticle *P = (GenParticle*)TCP->At(ipart);
			if(abs(P->PID) == PID)nPIDs++;
		}
	}

	file->Close();
	return nPIDs;
}


int get_nevents(string fileName){
	TFile* h_file = new TFile(fileName.c_str(),"READ");
	TTree* h_tree = (TTree*)h_file->Get("Delphes");
	int numberOfEntries = h_tree->GetEntries();
	h_file->Close();
	return numberOfEntries;

}

double get_weight(TBranch* branchE, int iEvent){
	TClonesArray *TCE = 0;
	branchE->SetAddress(&TCE);
	branchE->GetEntry(iEvent);  
	HepMCEvent *E = (HepMCEvent*)TCE->At(0);
	return E->Weight;
}


VectorXd get_event_MET(string fileName, int event){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	VectorXd MET = VectorXd::Zero(3);

	for (int iparticle = 0; iparticle < numberOfParticles; ++iparticle)
	{
		GenParticle *P = (GenParticle*)TCP->At(iparticle);
		MET(0)=MET(0)-P->Px;MET(1)=MET(1)-P->Py;
	}
	MET(2)=sqrt(pow(MET(0),2)+pow(MET(1),2));

	file->Close();
	return MET;
}


TLorentzVector get_event_MET_neutrino(string fileName, int event){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch *branchP = tree->GetBranch("Particle");
	TClonesArray *TCP = 0;
	branchP->SetAddress(&TCP);
	branchP->GetEntry(event);  
	Long64_t numberOfParticles = TCP->GetEntries();

	TLorentzVector MET;

	for (int iparticle = 0; iparticle < numberOfParticles; ++iparticle)
	{
		GenParticle *P = (GenParticle*)TCP->At(iparticle);
		if (P->PID == 12 || P->PID == -12 || P->PID == 14 || P->PID == -14 || P->PID == 18 || P->PID == -18)
		{
			TLorentzVector nu(P->Px,P->Py,P->Pz,P->E);
			MET = MET+nu;
		}
	}
	MET.SetE(MET.Vect().Mag());

	file->Close();
	return MET;
}



int get_numberOfPhotons(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bPhoton = tree->GetBranch("Photon");
	TClonesArray* TCPhoton = 0;
	bPhoton->SetAddress(&TCPhoton);
	int nEvents = bPhoton->GetEntries();

	int nPhotons = 0;

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bPhoton->GetEntry(iEvent);
		nPhotons += TCPhoton->GetEntries();
	}
	file->Close();
	return nPhotons;
}

int get_numberOfJets(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bJet = tree->GetBranch("Jet");
	TClonesArray* TCJet = 0;
	bJet->SetAddress(&TCJet);
	int nEvents = bJet->GetEntries();

	int nJet = 0;

	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bJet->GetEntry(iEvent);
		nJet += TCJet->GetEntries();
	}
	file->Close();
	return nJet;
}

int get_numberOfbJets(string fileName){
	TFile* file = new TFile(fileName.c_str(),"READ");
	TTree* tree = (TTree*)file->Get("Delphes");
	TBranch* bJet = tree->GetBranch("Jet");
	TClonesArray* TCJet = 0;
	bJet->SetAddress(&TCJet);
	int nEvents = bJet->GetEntries();
	int nbJet = 0;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent)
	{
		bJet->GetEntry(iEvent);
		int nJet = TCJet->GetEntries();
		for (int iJet = 0; iJet < nJet; ++iJet)
		{
			Jet* P_Jet = (Jet*)TCJet->At(iJet);
			if(P_Jet->BTag == 1)nbJet++;
		}
	}
	file->Close();
	return nbJet;
}


void ladebalken(int i, int max){
	double progress = (1.*i)/(1.*max)*100;
	#pragma omp critical
	std::cout << "\rSchritt " << i << " von " << max << " geschafft! " << "Also " << progress << " %";
	if(i == max-1 || i==max){std::cout << "\rSchritt " << max << " von " << max << " geschafft. Fertig!" << std::endl;}
	
	return;
}


void speichere(std::string name, MatrixXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}



void speichere(std::string name, VectorXd data){
	cout << ("\nSpeichere Datei " +name+ ".txt ab:").c_str() << endl;
	ofstream dat((name+".txt").c_str());
	dat.is_open();
	dat << data << "\n";
	dat.close();
	cout << "Datei " +name+ ".txt abgespeichert!\n" << endl;
	return;
}




TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[], double mass)
{
	if(nVec == 1)
	{
		return vecs[0];
	}
	else if(nVec == 2)
	{
		return vecs[0]+vecs[1];
	}
	else if(nVec == 3)
	{
		TLorentzVector perms[3];
		double perm_m[3];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();
		perms[2] = vecs[1]+vecs[2];perm_m[2]=perms[2].M();
		if (fabs(perm_m[0]-mass) > fabs(perm_m[1]))
		{
			if (fabs(perm_m[1]-mass) > fabs(perm_m[2]))
			{
				return perms[2];
			}else
			{
				return perms[1];
			}
		}
		if (fabs(perm_m[0]-mass) > fabs(perm_m[2]))
		{
			return perms[2];
		}
		return perms[0];
	}
	else if(nVec == 4)
	{
		TLorentzVector perms[6];
		double perm_m[6];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();perm_m[0]=fabs(perm_m[0]-mass);
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();perm_m[1]=fabs(perm_m[1]-mass);
		perms[2] = vecs[0]+vecs[3];perm_m[2]=perms[2].M();perm_m[2]=fabs(perm_m[2]-mass);
		perms[3] = vecs[1]+vecs[2];perm_m[3]=perms[3].M();perm_m[3]=fabs(perm_m[3]-mass);
		perms[4] = vecs[1]+vecs[3];perm_m[4]=perms[4].M();perm_m[4]=fabs(perm_m[4]-mass);
		perms[5] = vecs[2]+vecs[3];perm_m[5]=perms[5].M();perm_m[5]=fabs(perm_m[5]-mass);
		for (int i = 0; i < 6; ++i)
		{
			if((perm_m[i] < perm_m[(i+1)%6]) && (perm_m[i] < perm_m[(i+2)%6]) && (perm_m[i] < perm_m[(i+3)%6]) && (perm_m[i] < perm_m[(i+4)%6]) && (perm_m[i] < perm_m[(i+5)%6]))
			{
				return perms[i];
			}
		}
	}
	else if(nVec == 5)
	{
		TLorentzVector perms[10];
		double perm_m[10];
		perms[0] = vecs[0]+vecs[1];perm_m[0]=perms[0].M();perm_m[0]=fabs(perm_m[0]-mass);
		perms[1] = vecs[0]+vecs[2];perm_m[1]=perms[1].M();perm_m[1]=fabs(perm_m[1]-mass);
		perms[2] = vecs[0]+vecs[3];perm_m[2]=perms[2].M();perm_m[2]=fabs(perm_m[2]-mass);
		perms[3] = vecs[0]+vecs[4];perm_m[3]=perms[3].M();perm_m[3]=fabs(perm_m[3]-mass);
		perms[4] = vecs[1]+vecs[2];perm_m[4]=perms[4].M();perm_m[4]=fabs(perm_m[4]-mass);
		perms[5] = vecs[1]+vecs[3];perm_m[5]=perms[5].M();perm_m[5]=fabs(perm_m[5]-mass);
		perms[6] = vecs[1]+vecs[4];perm_m[6]=perms[6].M();perm_m[6]=fabs(perm_m[6]-mass);
		perms[7] = vecs[2]+vecs[3];perm_m[7]=perms[7].M();perm_m[7]=fabs(perm_m[7]-mass);
		perms[8] = vecs[2]+vecs[4];perm_m[8]=perms[8].M();perm_m[8]=fabs(perm_m[8]-mass);
		perms[9] = vecs[3]+vecs[4];perm_m[9]=perms[9].M();perm_m[9]=fabs(perm_m[9]-mass);
		for (int i = 0; i < 10; ++i)
		{
			if((perm_m[i] < perm_m[(i+1)%10]) && (perm_m[i] < perm_m[(i+2)%10]) && (perm_m[i] < perm_m[(i+3)%10]) && (perm_m[i] < perm_m[(i+4)%10]) && (perm_m[i] < perm_m[(i+5)%10]) && (perm_m[i] < perm_m[(i+6)%10]) && (perm_m[i] < perm_m[(i+7)%10]) && (perm_m[i] < perm_m[(i+8)%10]) && (perm_m[i] < perm_m[(i+9)%10]))
			{
				return perms[i];
			}
		}
	}
	else{
		cout << "ERROR: Not supported number of Jets" << endl;
		return vecs[0];
	}
}


TLorentzVector permute_to_mass_reco(int nVec, TLorentzVector vecs[])
{
	if(error < 2){error++;
	cout << "INFO: No Mass was given choosing highest PT" << endl;}
	if (nVec < 2)
	{	
		if(error < 2){error++;
		cout << "INFO: Only one Vector given -> Outout = Input" << endl;}
		return vecs[0];
	}

	TLorentzVector reco[2];
	TLorentzVector swap;
	reco[0] = vecs[0];
	reco[1] = vecs[1];
	for (int i = 2; i < nVec; ++i)
	{
		if (reco[0].Pt() < reco[1].Pt())
		{
			swap = reco[1]; reco[1]=reco[0]; reco[0]=swap;
		}
		if (vecs[i].Pt() > reco[1].Pt())
		{
			reco[1] = vecs[i];
		}
	}
	swap = reco[0]+reco[1];
	return swap;
}


