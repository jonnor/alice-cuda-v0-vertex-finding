// AliCudaV0vertexer.h

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

// Datastructure for v0 secondary vertices
struct v0vertex {
    struct trackparam *fParamN;
    struct trackparam *fParamP;
    Int_t *fNidx;
    Int_t *fPidx;
    Double_t fPos[3];

    Double_t fDcaV0Daughters;
    Double_t fChi2V0;
};

#define MAX_VERTICES 1000

// Datastructure for the primary vertex
struct privertex {
    Double_t fPosition[3];
    struct v0vertex fVertices[MAX_VERTICES];
};

Int_t cuda_v0_vertexer(struct privertex* vtx, struct trackparam* tracks, 
                        Int_t ntrks, Double_t b);

int test_cuda_v0_vertexer();
