
#include "AliCudaDefs.h"
//#include "aliexternaltrackparam.cu"
#include "aliv0vertexer.h"

/*
// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDv0.cxx.html#pUv_1E
__device__ __host__ v0_contructor(struct trackparam* t1, Int_t i1,
                                    struct trackparam* t2,  Int_t i2) {
    struct v0vertex v0;
    v0->fParamN = *t1;
    v0->fParamP = *t2;
    v0->fNidx = i1;
    v0->fPidx = i2;

    //TODO: find out how much of this is actually neccesary and implement. 
    //If we're gonna make this into a ordinary object on the CPU side, its 
    //probably just the dynamic/non-standard parts?
  fEffMass(TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass()),
  fDcaV0Daughters(0),
  fChi2V0(0.),
  fRr(0),
  fDistSigma(0),
  fChi2Before(0),
  fChi2After(0),
  fPointAngleFi(0),
  fPointAngleTh(0),
  fPointAngle(0),
  fPdgCode(kK0Short),
  fNidx(i1),
  fPidx(i2),
  fStatus(0),
  fNBefore(0),
  fNAfter(0),
  fOnFlyStatus(kFALSE)
{
  //--------------------------------------------------------------------
  // Main constructor  (K0s)
  //--------------------------------------------------------------------

  for (Int_t i=0; i<6; i++) {
    fPosCov[i]= 0.;
  }

  //Trivial estimation of the vertex parameters
  Double_t alpha=t1.GetAlpha(), cs=cos(alpha), sn=sin(alpha);
  Double_t tmp[3];
  t1.GetPxPyPz(tmp);
  Double_t px1=tmp[0], py1=tmp[1], pz1=tmp[2];
  t1.GetXYZ(tmp);
  Double_t  x1=tmp[0],  y1=tmp[1],  z1=tmp[2];
  const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
  Double_t sx1=sn*sn*t1.GetSigmaY2()+ss, sy1=cs*cs*t1.GetSigmaY2()+ss; 


  alpha=t2.GetAlpha(); cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  t2.GetPxPyPz(tmp);
  Double_t px2=tmp[0], py2=tmp[1], pz2=tmp[2];
  t2.GetXYZ(tmp);
  Double_t  x2=tmp[0],  y2=tmp[1],  z2=tmp[2];
  Double_t sx2=sn*sn*t2.GetSigmaY2()+ss, sy2=cs*cs*t2.GetSigmaY2()+ss; 
    
  Double_t sz1=t1.GetSigmaZ2(), sz2=t2.GetSigmaZ2();
  Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
  Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
  Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
  fPos[0]=wx1*x1 + wx2*x2; fPos[1]=wy1*y1 + wy2*y2; fPos[2]=wz1*z1 + wz2*z2;

  //fPos[0]=0.5*(x1+x2); fPos[1]=0.5*(y1+y2); fPos[2]=0.5*(z1+z2);
  fNmom[0]=px1; fNmom[1]=py1; fNmom[2]=pz1; 
  fPmom[0]=px2; fPmom[1]=py2; fPmom[2]=pz2;

  for (Int_t i=0;i<6;i++){fClusters[0][i]=0; fClusters[1][i]=0;}
  fNormDCAPrim[0]=fNormDCAPrim[1]=0;
  for (Int_t i=0;i<3;i++){fAngle[i]=0;}
  for (Int_t i=0;i<4;i++){fCausality[i]=0;}

}
*/

//NOTE: AliESDEvent event had to be replaced
// Information this function needs
// - All the tracks
// - The number of tracks
// - The magnetic field
// - The primary vertex

#define GetXv() fPosition[0]
#define GetYv() fPosition[1]
#define GetZv() fPosition[2]

#define MAXTRACKS 4000
__device__ Int_t neg[MAXTRACKS];
__device__ Int_t pos[MAXTRACKS];
__global__ void Tracks2V0vertices_kernel(struct privertex *vtxT3D, 
                                            struct trackparam *tracks, 
                                            Int_t nentr, Double_t b) {
   // ;
    //FIXME: when enabled, this call makes cuda very angry and causes a syntax error in the .ptx file
    //Int_t nvrtx = Tracks2V0vertices(vtxT3D, tracks, nentr, b);
// }


/*
// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliV0vertexer.cxx.html#gwutwE
__device__ __host__ Int_t Tracks2V0vertices(struct privertex *vtxT3D, 
                                            struct trackparam *tracks, 
                                            Int_t nentr, Double_t b) {
*/
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

//  if (nentr<2) return 0; 

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

   //Tries to match negative tracks with positive to find v0s
   for (i=0; i<nneg; i++) {
      Int_t nidx=neg[i];
      struct trackparam *ntrk=&tracks[nidx]; // AliESDtrack *ntrk=event->GetTrack(nidx);

      for (Int_t k=0; k<npos; k++) {
         Int_t pidx=pos[k];
	 struct trackparam *ptrk=&tracks[pidx]; // AliESDtrack *ptrk=event->GetTrack(pidx);

         if (abs(GetD(ntrk,xPrimaryVertex,yPrimaryVertex,b))<fDNmin)
	   if (abs(GetD(ptrk,xPrimaryVertex,yPrimaryVertex,b))<fDNmin) continue;

         Double_t xn, xp, dca=GetDCA(ntrk,ptrk,b,xn,xp);
         if (dca > fDCAmax) continue;
         if ((xn+xp) > 2*fRmax) continue;
         if ((xn+xp) < 2*fRmin) continue;

         // XXX: might need a cast instead, and should probably be copies
         struct trackparam nt(*ntrk), pt(*ptrk); 
         Bool_t corrected=kFALSE;
         if ((GetX(&nt) > 3.) && (xn < 3.)) {
	   //correct for the beam pipe material
           corrected=kTRUE;
         }
         if ((GetX(&pt) > 3.) && (xp < 3.)) {
	   //correct for the beam pipe material
           corrected=kTRUE;
         }
         if (corrected) {
	   dca=GetDCA(&nt,&pt,b,xn,xp);
           if (dca > fDCAmax) continue;
           if ((xn+xp) > 2*fRmax) continue;
           if ((xn+xp) < 2*fRmin) continue;
	 }

         PropagateTo(&nt,xn,b); PropagateTo(&pt,xp,b);

/*
         struct v0vertex* vertex = v0vertex_contructor(nt, nidx, pt, pidx);
         if (GetChi2V0(vertex) > fChi2max) continue;
	 
	 Float_t cpa=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
	 if (cpa < fCPAmin) continue;
	 SetDcaV0Daughters(vertex,dca);
         SetV0CosineOfPointingAngle(vertex,cpa);
         ChangeMassHypothesis(vertex,kK0Short);

          //TODO: find and implement an equivalent way to do this. Just use an array?
         //event->AddV0(&vertex); http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDEvent.cxx.html#qqk9g
*/
         nvtx++;
      }
    }

//    Info("Tracks2V0vertices","Number of reconstructed V0 vertices: %d",nvtx);

//   return nvtx;
}
