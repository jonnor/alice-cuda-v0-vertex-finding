
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
        ULong_t fFlags;
};

__device__ __host__ Double_t GetDCA(struct trackparam *tp1, struct trackparam *tp2,
                                    Double_t b, Double_t &xthis, Double_t &xp);
__host__ __device__ Double_t GetSign(struct trackparam *tp);
__host__ __device__ Double_t GetC(struct trackparam *tp, Double_t b);
__host__ __device__ Bool_t GetPxPyPz(struct trackparam *tp, Double_t p[3]);
__host__ __device__ Bool_t GetXYZ(struct trackparam *tp, Double_t *r);
__host__ __device__ void GetHelixParameters(struct trackparam *tp, Double_t hlx[6], Double_t b);
__host__ __device__ void Evaluate(const Double_t *h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]); //FIXME: should be static?
__host__ __device__ Bool_t PropagateTo(struct trackparam *tp, Double_t xk, Double_t b);
__host__ __device__ Double_t GetLinearD(struct trackparam *tp, Double_t xv,Double_t yv);
__host__ __device__ Double_t GetD(struct trackparam *tp, Double_t x,Double_t y,Double_t b);


#endif
