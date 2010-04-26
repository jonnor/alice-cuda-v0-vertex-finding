
#include "common.h"
//#include "aliexternaltrackparam.cu"
#include "aliv0vertexer.h"


//NOTE: AliESDEvent event had to be replaced
// Information this function needs
// - All the tracks
// - The number of tracks
// - The magnetic field
// - The primary vertex
__global__ void Tracks2V0vertices_kernel(struct privertex *vtxT3D, 
                                            struct trackparam *tracks, 
                                            Int_t nentr, Double_t b) {
    ;
    //FIXME: when enabled, this call makes cuda very angry and causes a syntax error in the .ptx file
    //Int_t nvrtx = Tracks2V0vertices(vtxT3D, tracks, nentr, b);
}

#define GetXv() fPosition[0]
#define GetYv() fPosition[1]
#define GetZv() fPosition[2]
__device__ __host__ Int_t Tracks2V0vertices(struct privertex *vtxT3D, 
                                            struct trackparam *tracks, 
                                            Int_t nentr, Double_t b) {

    //TODO: put this outside function
    // http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliV0vertexer.cxx.html#27
    Double_t fgChi2max=33.; //max chi2
    Double_t fgDNmin=0.05;  //min imp parameter for the 1st daughter
    Double_t fgDPmin=0.05;  //min imp parameter for the 2nd daughter
    Double_t fgDCAmax=0.5;  //max DCA between the daughter tracks
    Double_t fgCPAmin=0.99; //min cosine of V0's pointing angle
    Double_t fgRmin=0.2;    //min radius of the fiducial volume
    Double_t fgRmax=100.;   //max radius of the fiducial volume

    // http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliV0vertexer.h.html#46
    Double_t fChi2max = fgChi2max; 
    Double_t fDNmin = fgDNmin;
    Double_t fDPmin = fgDPmin;
    Double_t fDCAmax = fgDCAmax;
    Double_t fCPAmin = fgCPAmin;
    Double_t fRmin = fgRmin;
    Double_t fRmax = fgRmax;

   Double_t xPrimaryVertex=vtxT3D->GetXv();
   Double_t yPrimaryVertex=vtxT3D->GetYv();
   Double_t zPrimaryVertex=vtxT3D->GetZv();

   if (nentr<2) return 0; 

   Int_t neg[nentr];
   Int_t pos[nentr];

   Int_t nneg=0, npos=0, nvtx=0;

   //Sorts all tracks that meet certain cuts into positive and negative
   Int_t i;
   for (i=0; i<nentr; i++) {
     struct trackparam *esdTrack=&tracks[i]; // AliESDtrack *esdTrack=event->GetTrack(i);
//     ULong_t status=esdTrack->GetStatus();

//      if ((status&AliESDtrack::kITSrefit)==0)
//         if ((status&AliESDtrack::kTPCrefit)==0) continue;

     Double_t d=GetD(esdTrack,xPrimaryVertex,yPrimaryVertex,b);
     if (abs(d)<fDPmin) continue;
     if (abs(d)>fRmax) continue;

     if (GetSign(esdTrack) < 0.) neg[nneg++]=i;
     else pos[npos++]=i;
   }

//TODO: not implemented
//    //Tries to match negative tracks with positive to find v0s
//    for (i=0; i<nneg; i++) {
//       Int_t nidx=neg[i];
//       AliESDtrack *ntrk=event->GetTrack(nidx);
//
//       for (Int_t k=0; k<npos; k++) {
//          Int_t pidx=pos[k];
// 	 AliESDtrack *ptrk=event->GetTrack(pidx);
//
//          if (TMath::Abs(ntrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin)
// 	   if (TMath::Abs(ptrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin) continue;
//
//          Double_t xn, xp, dca=ntrk->GetDCA(ptrk,b,xn,xp);
//          if (dca > fDCAmax) continue;
//          if ((xn+xp) > 2*fRmax) continue;
//          if ((xn+xp) < 2*fRmin) continue;
//    
//          AliExternalTrackParam nt(*ntrk), pt(*ptrk);
//          Bool_t corrected=kFALSE;
//          if ((nt.GetX() > 3.) && (xn < 3.)) {
// 	   //correct for the beam pipe material
//            corrected=kTRUE;
//          }
//          if ((pt.GetX() > 3.) && (xp < 3.)) {
// 	   //correct for the beam pipe material
//            corrected=kTRUE;
//          }
//          if (corrected) {
// 	   dca=nt.GetDCA(&pt,b,xn,xp);
//            if (dca > fDCAmax) continue;
//            if ((xn+xp) > 2*fRmax) continue;
//            if ((xn+xp) < 2*fRmin) continue;
// 	 }
//
//          nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);

// TODO: not implemented
//          AliESDv0 vertex(nt,nidx,pt,pidx);
//          if (vertex.GetChi2V0() > fChi2max) continue;
// 	 
// 	 Float_t cpa=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
// 	 if (cpa < fCPAmin) continue;
// 	 vertex.SetDcaV0Daughters(dca);
//          vertex.SetV0CosineOfPointingAngle(cpa);
//          vertex.ChangeMassHypothesis(kK0Short);
//
//          event->AddV0(&vertex);

//          nvtx++;
//       }
//    }

//    Info("Tracks2V0vertices","Number of reconstructed V0 vertices: %d",nvtx);

    return nvtx;
}
