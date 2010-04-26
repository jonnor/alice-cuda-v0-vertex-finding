
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cuda.h>

#include "common.h"
#include "aliexternaltrackparam.cu"
#include "aliv0vertexer.cu"

int main()
{

    Double_t b = -5.00668;
    const int TRACK_SIZE = sizeof(struct trackparam);
    const int HELIX_SIZE = sizeof(Double_t)*6;

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

    //FIXME: should be pointer to struct trackparam instead??
    const int NTRACKS=2;
    struct trackparam *tracks[NTRACKS];
    tracks[0] = tp;
    tracks[1] = tp2;

    printf("Tracks2V0vertices\n");
    Tracks2V0vertices(vtxT3d, *tracks, NTRACKS, b);

    // Cleanup
    free(tp); free(hp); free(tp2); free(vtxT3d);

    return 1;
}
