#ifndef ALIV0VERTEXER_H
#define ALIV0VERTEXER_H

#include "aliexternaltrackparam.h"

// Datastructure for v0 secondary vertices
struct v0vertex {
    void *fParamN;
    void *fParamP;
    Double_t fPos[3];

    Double_t fDcaV0Daughters;
    Double_t fChi2V0;
};

const int MAX_VERTICES = 1000;

// Datastructure for the primary vertex
struct privertex {
    Double_t fPosition[3];
    struct v0vertex *fVertices[MAX_VERTICES];
};


__global__ void Tracks2V0vertices_kernel(struct privertex *vtxT3d,
                                            struct trackparam *tracks,
                                            Int_t nentr, Double_t b);
__device__ __host__ Int_t Tracks2V0vertices(struct privertex *vtxT3d,
                                            struct trackparam *tracks,
                                            Int_t nentr, Double_t b);

#endif
