
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cuda.h>

#include "common.h"
#include "aliexternaltrackparam.h"

int main()
{

    Double_t b = -5.00668;
    const int TRACK_SIZE = sizeof(struct trackparam);
    const int HELIX_SIZE = sizeof(Double_t)*6;

    struct trackparam *tp;
    struct trackparam *tp_d;

    Double_t *hp;
    Double_t *hp_d;

    // Allocate memory
    tp = (struct trackparam*)malloc(TRACK_SIZE);
    hp = (Double_t *)malloc(HELIX_SIZE);
    cudaMalloc((void **) &tp_d, TRACK_SIZE);
    cudaMalloc((void **) &hp_d, HELIX_SIZE);

    // Initialize data
    tp->fP[0] = -0.00429969;
    tp->fP[1] = -4.56162;
    tp->fP[2] = 2.38928e-09;
    tp->fP[3] = 0.855637;
    tp->fP[4] = -1.96397;
    tp->fAlpha = 1.90909;
    tp->fX = -0.00462971;
    for(int i=0; i<6;i++) hp[i] = 0;

    // Transfer data and do computation
    cudaMemcpy(tp_d, tp, TRACK_SIZE, cudaMemcpyHostToDevice);
    printf("GetHelixParameters\n");
    GetHelixParameters <<<1,1>>> (tp_d, hp_d, b);

    // Retrieve data and check results
    cudaMemcpy(hp, hp_d, HELIX_SIZE, cudaMemcpyDeviceToHost);
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

    // Cleanup
    free(tp); free(hp);
    cudaFree(tp_d); cudaFree(hp_d);

    return 1;
}
