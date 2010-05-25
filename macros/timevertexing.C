#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TXMLEngine.h"
#include "TXMLSetup.h"
#include "TCanvas.h"
#include "TStorage.h"
#include <stdio.h>
#include <stdlib.h>

AliESDEvent* load_event() {
	// Load .root
	// http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/root_crib/ttree.html
	TFile *esdfile = new TFile("AliESDs.root","READ");
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
	AliESDEvent *eventcopy = new AliESDEvent();
	for(i=0;i<RUNS;i++){
	 // AliESDEvent *eventcopy = new AliESDEvent();
	  (*eventcopy) = (*event);
	  tstopw->Start();
	  Int_t N_v0 = v0vertexer.Tracks2V0vertices(eventcopy);
	  tstopw->Stop();
	  array[i] = tstopw->RealTime();
	  //delete eventcopy;
	  eventcopy->Reset();
	}

}

TH1F* draw_histogram(Float_t array[], Int_t N, Int_t nbinsx) {
	
	Float_t M = 0;
	Float_t xlow = 0;
	Float_t xup = 0;
	Float_t mean = 0;
	for(i=0;i<N;i++){
		M = M + array[i];
	}	
	mean = M / N;
	xlow = mean * 0.96;
	xup = mean * 1.04;

	// Draw histogram from data
	TCanvas *c1;
	// Create canvas with grid
	c1 = new TCanvas("c1","Blopp",600,400);
	c1->SetGrid();
	// create histogram
	hist = new TH1F("times","V0Vertexer time distribution", nbinsx, xlow, xup);
	hist->SetFillColor(16);
	hist->SetMarkerStyle(20);
	hist->SetMarkerSize(0.7);
	hist->Draw("hmm");

for(int i=0; i<N; i++){ 
	hist->Fill(array[i]); }

  c1->Modified();
  c1->Update();
  c1->Print("histogram.pdf","pdfLandscape");
	return hist;
}


void PrintArray(Float_t array[], Int_t N)
{
// Print results
for(i=0;i<N;i++)
 printf("%f",array[i]);
}
void save_data(Float_t array[], Int_t N, Char_t filename[]){
	//TFile *file = new TFile("Data.xml","RECREATE");  //create file
	//hist->Write();
	//delete file;
	FILE *file;
	file = fopen(filename,"a+");
	for(i=0;i<N;i++)	
		fprintf(file, "%f\n",array[i]);
	fclose(file);
	//return 0;
}

void load_data(Float_t array[], Char_t filename[], Int_t N){
        FILE *ifile;
        ifile = fopen(filename,"r");
        Int_t i = 0;
        //array = new Float_t[N];
        for(i = 0; i<N; i++){		
		Int_t ret = fscanf(ifile, "%f\n", &array[i]);
		if(ret == EOF) break;	
	}
        fclose(ifile);
        //return 0;
}

