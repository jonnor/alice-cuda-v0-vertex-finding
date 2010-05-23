

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDv0.cxx.html#HTK9SE
void __device__ __host__ GetPxPyPz(struct v0vertex* v0, Double_t &px, Double_t &py, Double_t &pz) {
  //--------------------------------------------------------------------
  // This function returns V0's momentum (global)
  //--------------------------------------------------------------------
  px=v0->fNmom[0]+v0->fPmom[0]; 
  py=v0->fNmom[1]+v0->fPmom[1]; 
  pz=v0->fNmom[2]+v0->fPmom[2]; 
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDv0.cxx.html#qIMmkE
Float_t __device__ __host__ GetV0CosineOfPointingAngle(struct v0vertex* v0, 
                Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
  // calculates the pointing angle of the V0 wrt a reference point

  Double_t momV0[3]; //momentum of the V0
  GetPxPyPz(v0, momV0[0],momV0[1],momV0[2]);

  Double_t deltaPos[3]; //vector between the reference point and the V0 vertex
  deltaPos[0] = v0->fPos[0] - refPointX;
  deltaPos[1] = v0->fPos[1] - refPointY;
  deltaPos[2] = v0->fPos[2] - refPointZ;

  Double_t momV02    = momV0[0]*momV0[0] + momV0[1]*momV0[1] + momV0[2]*momV0[2];
  Double_t deltaPos2 = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];

  Double_t cosinePointingAngle = (deltaPos[0]*momV0[0] +
				  deltaPos[1]*momV0[1] +
				  deltaPos[2]*momV0[2] ) / sqrt(momV02 * deltaPos2);
  
  return cosinePointingAngle;
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliESDv0.html#AliESDv0:GetChi2V0
Double_t __device__ __host__ GetChi2V0(struct v0vertex* v0) {
    return v0->fChi2V0;
}

#define fPos v0->fPos
#define fNmom v0->fNmom
#define fPmom v0->fPmom
// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDv0.cxx.html#pUv_1E
struct v0vertex* __device__ __host__ v0vertex_contructor(struct trackparam* t1, Int_t i1,
                                    struct trackparam* t2,  Int_t i2) {
    struct v0vertex v0vtx;
    struct v0vertex* v0 = &v0vtx;
    v0->fParamN = *t1;
    v0->fParamP = *t2;
    v0->fNidx = i1;
    v0->fPidx = i2;

    //TODO: find out how much of this is actually neccesary and implement. 
    //If we're gonna make this into a ordinary object on the CPU side, its 
    //probably just the dynamic/non-standard parts?
/*
  fEffMass(TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass()),

  fPdgCode(kK0Short),
  fNidx(i1),
  fPidx(i2),

  fOnFlyStatus(kFALSE)
*/
  //Trivial estimation of the vertex parameters
  Double_t alpha=t1->fAlpha, cs=cos(alpha), sn=sin(alpha);
  Double_t tmp[3];
  GetPxPyPz(t1, tmp);
  Double_t px1=tmp[0], py1=tmp[1], pz1=tmp[2];
  GetXYZ(t1, tmp);
  Double_t  x1=tmp[0],  y1=tmp[1],  z1=tmp[2];
  const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
  Double_t sx1=sn*sn*t1->fC[0]+ss, sy1=cs*cs*t1->fC[0]+ss; 


  alpha=t2->fAlpha; cs=cos(alpha); sn=sin(alpha);
  GetPxPyPz(t2, tmp);
  Double_t px2=tmp[0], py2=tmp[1], pz2=tmp[2];
  GetXYZ(t2, tmp);
  Double_t  x2=tmp[0],  y2=tmp[1],  z2=tmp[2];
  Double_t sx2=sn*sn*t2->fC[0]+ss, sy2=cs*cs*t2->fC[0]+ss; 
    
  Double_t sz1=t1->fC[2], sz2=t2->fC[2];
  Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
  Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
  Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
  fPos[0]=wx1*x1 + wx2*x2; fPos[1]=wy1*y1 + wy2*y2; fPos[2]=wz1*z1 + wz2*z2;

  fNmom[0]=px1; fNmom[1]=py1; fNmom[2]=pz1; 
  fPmom[0]=px2; fPmom[1]=py2; fPmom[2]=pz2;

  return v0;
}
#undef fPos
#undef fPmom
#undef fNmom


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliV0vertexer.cxx.html#27
// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliV0vertexer.h.html#46
#define fChi2max 33. //max chi2
#define fDNmin 0.05  //min imp parameter for the 1st daughter
#define fDPmin 0.05  //min imp parameter for the 2nd daughter
#define fDCAmax 0.5  //max DCA between the daughter tracks
#define fCPAmin 0.99 //min cosine of V0's pointing angle
#define fRmin 0.2    //min radius of the fiducial volume
#define fRmax 100.   //max radius of the fiducial volume


__device__ __host__ Int_t CompareTracks(struct trackparam* ntrk, struct trackparam* ptrk, 
                                    Int_t nidx, Int_t pidx, Double_t b,
                                    Double_t xPrimaryVertex, Double_t yPrimaryVertex, Double_t zPrimaryVertex) {

     if ( (abs(GetD(ntrk,xPrimaryVertex,yPrimaryVertex,b))<fDNmin) && 
            (abs(GetD(ptrk,xPrimaryVertex,yPrimaryVertex,b))<fDNmin) ) return 0;

     Double_t xn, xp, dca;
     dca=GetDCA(ntrk,ptrk,b,xn,xp);
     if (dca > fDCAmax) return 0;
     if ((xn+xp) > 2*fRmax) return 0;
     if ((xn+xp) < 2*fRmin) return 0;

     struct trackparam nt, pt;
     nt=(*ntrk); pt=(*ptrk); 

     PropagateTo(&nt,xn,b); PropagateTo(&pt,xp,b);

     struct v0vertex* vertex = v0vertex_contructor(&nt, nidx, &pt, pidx);
 
     Float_t cpa=GetV0CosineOfPointingAngle(vertex, 
                        xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
     if (cpa < fCPAmin) return 0;

    return 1; // This is a match
}


#define GetXv() fPosition[0]
#define GetYv() fPosition[1]
#define GetZv() fPosition[2]

__global__ void SortTracks_kernel(struct privertex *vtxT3D, 
                                    struct trackparam *tracks, 
                                    Int_t *pos, Int_t *neg, 
                                    Int_t *npos_ptr, Int_t *nneg_ptr, 
                                    Int_t nentr, Int_t b) {

if ((blockIdx.x * blockDim.x + threadIdx.x) == 1) { 
   Double_t xPrimaryVertex=vtxT3D->GetXv();
   Double_t yPrimaryVertex=vtxT3D->GetYv();
   Double_t zPrimaryVertex=vtxT3D->GetZv();

   (*nneg_ptr)=0; (*npos_ptr)=0;
   //Sorts all tracks that meet certain cuts into positive and negative
   Int_t i;
   for (i=0; i<nentr; i++) {
     struct trackparam *esdTrack=&tracks[i]; // AliESDtrack *esdTrack=event->GetTrack(i);
    ULong_t status=esdTrack->fFlags;
     if ((status&kITSrefit)==0)
        if ((status&kTPCrefit)==0) continue;
     Double_t d=GetD(esdTrack,xPrimaryVertex,yPrimaryVertex,b);
     if (abs(d)<fDPmin) continue;
     if (abs(d)>fRmax) continue;

     if (GetSign(esdTrack) < 0.) neg[(*nneg_ptr)++]=i;
     else pos[(*npos_ptr)++]=i;
   }
}

}


__global__ void Tracks2V0vertices_kernel(struct privertex *vtxT3D, 
                                            struct trackparam *tracks,
                                            Int_t *pos, Int_t *neg,
                                            Int_t *npos_ptr, Int_t *nneg_ptr, 
                                            Int_t *nv0s_ptr, Double_t b) {

   Double_t xPrimaryVertex=vtxT3D->GetXv();
   Double_t yPrimaryVertex=vtxT3D->GetYv();
   Double_t zPrimaryVertex=vtxT3D->GetZv();

   Int_t npos = *npos_ptr;
   Int_t nneg = *nneg_ptr;
   (*nv0s_ptr)=0;

   Int_t id=(blockIdx.x * blockDim.x) + threadIdx.x;
   //Tries to match negative tracks with positive to find v0s
   if (id < nneg) {
      Int_t nidx=neg[id];
      struct trackparam *ntrk=&tracks[nidx];

      Int_t k;
      for (k=0; k<npos; k++) {
         Int_t pidx=pos[k];
	 struct trackparam *ptrk=&tracks[pidx];

        if (CompareTracks(ntrk, ptrk, nidx, pidx, b, 
            xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) atomicAdd(nv0s_ptr,1);
      }
    }


}

/*
	 SetDcaV0Daughters(vertex,dca);
         SetV0CosineOfPointingAngle(vertex,cpa);
         ChangeMassHypothesis(vertex,kK0Short);

          //TODO: find and implement an equivalent way to do this. Just use an array?
         //event->AddV0(&vertex); http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDEvent.cxx.html#qqk9g
*/       
