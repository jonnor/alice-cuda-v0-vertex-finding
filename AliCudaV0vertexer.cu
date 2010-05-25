/* This file is part of a project that implements GPU based 
 * v0 vertex finding for use with AliROOT in ALICE
 * http://ri-pro.hive.no/prosjekter/EN2010-01/ 
 * Code at http://gitorious.org/cuda-alice-vertex-finding
 */

/**************************************************************************
 * Copyright(c) 2010  Vestfold University College, All rights reserved.   *
 *                                                                        *
 * Authors: Jon Nordby, Lars Bratrud                                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cuda.h>

#include "AliCudaDefs.h"
#include "AliCudaV0vertexer.h"
#include "aliexternaltrackparam.cu"
#include "aliv0vertexer.cu"

const int TRACK_SIZE = sizeof(struct trackparam);
const int HELIX_SIZE = sizeof(Double_t)*6;
const int VERTEX_SIZE = sizeof(struct privertex);

Int_t cuda_v0_vertexer(struct privertex* vtx, struct trackparam* tracks, 
                        Int_t ntrks, Double_t b) {

    if (ntrks < 2) return 0;

    // Host data
    Int_t nv0s=0;

    // Declare and allocate device data
    struct trackparam* tracks_d;
    struct privertex* vtx_d;
    cudaMalloc((void**)&vtx_d, VERTEX_SIZE);
    cudaMalloc((void**)&tracks_d, TRACK_SIZE*ntrks);

    Int_t *nv0s_d;
    cudaMalloc((void**)&nv0s_d, sizeof(Int_t));

    Int_t *ptracks_d;
    Int_t *ntracks_d;
    cudaMalloc((void**)&ntracks_d, sizeof(Int_t)*ntrks);
    cudaMalloc((void**)&ptracks_d, sizeof(Int_t)*ntrks);

    Int_t *npos_d; Int_t *nneg_d;
    cudaMalloc((void**)&npos_d, sizeof(Int_t));
    cudaMalloc((void**)&nneg_d, sizeof(Int_t));
    printf("Malloc DONE\n");

    // Copy data to device
    cudaMemcpy(tracks_d, tracks, TRACK_SIZE*ntrks, cudaMemcpyHostToDevice);
    cudaMemcpy(vtx_d, vtx, VERTEX_SIZE, cudaMemcpyHostToDevice);
    cudaMemcpy(nv0s_d, &nv0s, sizeof(Int_t), cudaMemcpyHostToDevice);
    printf("Memcpy HostToDevice DONE\n");

    // Execute
    const Int_t warpsize = 32;
    dim3  grid(ceil(ntrks/warpsize), 1, 1); // x,y,z
    dim3  block(warpsize, 1, 1); // x,y,z

    SortTracks_kernel<<<grid, block>>>(vtx_d, tracks_d, 
                                        ptracks_d, ntracks_d, 
                                        npos_d, nneg_d,
                                        ntrks, b);
    cudaThreadSynchronize();
    printf("SortTracks kernel execution DONE\n");


    // DEBUG
    Int_t npos, nneg;
    cudaMemcpy(&npos, npos_d, sizeof(Int_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(&nneg, nneg_d, sizeof(Int_t), cudaMemcpyDeviceToHost);
    printf("SortTracks: npos=%d, nneg=%d\n", npos, nneg);

    Tracks2V0vertices_kernel<<<grid, block>>>(vtx_d, tracks_d,
                                                ptracks_d, ntracks_d, 
                                                npos_d, nneg_d,
                                                nv0s_d, b);
    cudaThreadSynchronize();
    printf("Tracks2V0vertices kernel execution DONE\n");

    // Copy data back and clean up
    cudaMemcpy(vtx, vtx_d, VERTEX_SIZE, cudaMemcpyDeviceToHost);
    cudaFree(vtx_d); cudaFree(tracks_d);

    cudaMemcpy(&nv0s, nv0s_d, sizeof(Int_t), cudaMemcpyDeviceToHost);
    cudaFree(nv0s_d);

    return nv0s;
}

int test_cuda_v0_vertexer()
{

    Double_t b = -5.00668;
    struct trackparam *tp;
    Double_t *hp;

    // Allocate memory
    tp = (struct trackparam*)malloc(TRACK_SIZE);
    hp = (Double_t *)malloc(HELIX_SIZE);

    // Initialize data
    tp->fP[0] = -0.00429969;
    tp->fP[1] = -4.56162;
    tp->fP[2] = 2.38928e-09;
    tp->fP[3] = 0.855637;
    tp->fP[4] = -1.96397;
    tp->fAlpha = 1.90909;
    tp->fX = -0.00462971;
    for(int i=0; i<6;i++) hp[i] = 0;

    printf("GetHelixParameters\n");
    GetHelixParameters(tp, hp, b);
    for (int i=0; i<6; i++) printf("%f\n", hp[i]);

    // TODO: find real inputdata and test
    Double_t t = 3; 
    Double_t rv[3], d[3], dd[3];
    Evaluate(hp, t, rv, d, dd);

    // TODO: find real inputdata and test
    Double_t xk = 50.0;
    PropagateTo(tp, xk, b);
    printf("PropagateTo\n");
    for (int i=0; i<5; i++) printf("%f\n",tp->fP[i]);

    // TODO: find real inputdata and test
    Double_t xv = 1.0, yv = 1.0;
    Double_t d_lin = GetLinearD(tp, xv, yv);
    printf("GetLinearD = %f\n", d_lin);
    Double_t d_ = GetD(tp, xv, yv, b);
    printf("GetD = %f\n", d_);


    // TODO: find real inputdata and test. Must include covariance matrix
    struct trackparam *tp2;
    tp2 = (struct trackparam*)malloc(TRACK_SIZE);

    tp2->fP[0] = -0.0315877;
    tp2->fP[1] = -4.54952;
    tp2->fP[2] = 3.74249e-09;
    tp2->fP[3] = 1.15249;
    tp2->fP[4] = 1.67247;
    tp2->fAlpha = 0.107172;
    tp2->fX = 0.000891429;

    Double_t xp=1.0, xn=1.0;
    Double_t dca = 1.0;
    dca = GetDCA(tp, tp2, b, xn, xp);
    printf("GetDCA = %f\n", dca);

    // v0 vertexer
    struct privertex *vtxT3d;
    vtxT3d = (struct privertex *) malloc(sizeof(struct privertex));
    vtxT3d->fPosition[0] = 0.0;
    vtxT3d->fPosition[1] = 0.0;
    vtxT3d->fPosition[2] = 0.0;

    const int NTRACKS=2;
    struct trackparam *tracks;
    tracks = (struct trackparam*)malloc(sizeof(struct trackparam)*NTRACKS);
    tracks[0] = *tp;
    tracks[1] = *tp2;

    printf("Tracks2V0vertices\n");
//    Tracks2V0vertices(vtxT3d, tracks, NTRACKS, b);
//     Tracks2V0vertices_kernel<<<1,1>>>(vtxT3d, tracks, NTRACKS, b);
//     cudaThreadSynchronize();

    cuda_v0_vertexer(vtxT3d, tracks, NTRACKS, b);

    // Cleanup
    free(tp); free(hp); free(tp2); free(vtxT3d);

    return 1;
}

int main() {
    return test_cuda_v0_vertexer();
}
