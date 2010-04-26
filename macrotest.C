#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include "TStopwatch.h"
#define RUNS 20

AliESDEvent* load_event() {
	// Load .root
	// http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/ttree.html
	TFile *esdfile = new TFile("AliESDs.root","READ");  //read root-file
	TTree *esdTree = (TTree *) esdfile->Get("esdTree"); //read tree

	AliESDEvent *event = new AliESDEvent();
	event->ReadFromTree(esdTree);
	esdTree->GetEvent();
	return event;
}

time_vertexing(AliESDEvent* event, Float_t* array, Int_t RUNS) {
	AliV0vertexer v0vertexer;
	TStopwatch *tstopw = new TStopwatch();
	Int_t i;

  for(i=0;i<RUNS;i++){
  tstopw->Start();
  Int_t N_v0 = v0vertexer.Tracks2V0vertices(event);
  tstopw->Stop();
  array[i] = tstopw->RealTime();
 }

}

TH1F* draw_histrogram(Float_t array[]) {
// Draw histogram from data
	TCanvas *c1;
	// Create canvas with grid
	c1 = new TCanvas("c1","Blopp",600,400);
	c1->SetGrid();
	// create histogram
	hist = new TH1F("times","V0Vertexer time distribution",100,1.4,1.6);
	hist->SetFillColor(16);
	hist->SetMarkerStyle(21);
	hist->SetMarkerSize(0.7);
	hist->Draw("hmm");

for(int i; i<RUNS; i++){ 
	hist->Fill(array[i]); }

  c1->Modified();
  c1->Update();
	return hist;
}


void PrintArray(Float_t array[], Int_t N)
{
// Print results
for(i=0;i<N;i++)
 printf("%f",array[i]);
}

