
#ifndef ALIEXTERNALTRACKPARAM_H
#define ALIEXTERNALTRACKPARAM_H

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#178
// We could also just use a simple array instead of the struct
// Double_t trackparam[7]
// When we start caring about multiple tracks, they can either be represented by
// - An array of structs/arrays
// - An array where the params are separated by a rowstride equal to number of tracks
// The latter might perform better, but be less readable
struct trackparam
{
	Double_t fP[5];
	Double_t fAlpha;
	Double_t fX;
        Double_t fC[15];
};

Double_t GetSign(struct trackparam *tp);
__host__ __device__ Double_t GetC(struct trackparam *tp, Double_t b);
Bool_t GetPxPyPz(struct trackparam *tp, Double_t p[3]);
Bool_t GetXYZ(struct trackparam *tp, Double_t *r);
__global__ void GetHelixParameters(struct trackparam *tp, Double_t hlx[6], Double_t b);
void Evaluate(const Double_t *h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]); //FIXME: should be static?
Bool_t PropagateTo(struct trackparam *tp, Double_t xk, Double_t b);
Double_t GetLinearD(struct trackparam *tp, Double_t xv,Double_t yv);
Double_t GetD(struct trackparam *tp, Double_t x,Double_t y,Double_t b);


#endif
