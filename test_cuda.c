
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <cuda.h>

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliVParticle.h.html#18
const float kB2C=-0.299792458e-3;
const float kAlmost0=FLT_MIN;

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#178
// We could also just use a simple array instead of the struct
// float trackparam[7]
// When we start caring about multiple tracks, they can either be represented by
// - An array of structs/arrays
// - An array where the params are separated by a rowstride equal to number of tracks
// The latter might perform better, but be less readable
struct trackparam
{
	float fP[5];
	float fAlpha;
	float fX;
};

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#XB.FNC
__device__ float GetC(struct trackparam *tp, float b) {
    return tp->fP[4]*b*kB2C;
}


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
__global__ void GetHelixParameters(struct trackparam *tp, float hlx[6], float b) {

    float cs=cos(tp->fAlpha), sn=sin(tp->fAlpha); 

    hlx[0]=tp->fP[0]; hlx[1]=tp->fP[1]; hlx[2]=tp->fP[2]; hlx[3]=tp->fP[3];

    hlx[5]=tp->fX*cs - hlx[0]*sn;               // x0
    hlx[0]=tp->fX*sn + hlx[0]*cs;               // y0
    //hlx[1]=                                 // z0
    hlx[2]=sin(hlx[2]) + tp->fAlpha;    // phi0
    //hlx[3]=                                 // tgl
    hlx[4]=GetC(tp, b);                         // C
} 


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
static void Evaluate(const float *h, float t,
                     float r[3],  //radius vector
                     float g[3],  //first defivatives
                     float gg[3]) //second derivatives
{

  float phase=h[4]*t+h[2];
  float sn=sin(phase), cs=cos(phase);

  r[0] = h[5];
  r[1] = h[0];
  if (fabs(h[4])>kAlmost0) {
     r[0] += (sn - h[6])/h[4];
     r[1] -= (cs - h[7])/h[4];  
  }
  r[2] = h[1] + h[3]*t;

  g[0] = cs; g[1]=sn; g[2]=h[3];
 
  gg[0]=-h[4]*sn; gg[1]=h[4]*cs; gg[2]=0.;
}



int main()
{

    float b = -5.00668;
    const int TRACK_SIZE = sizeof(struct trackparam);
    const int HELIX_SIZE = sizeof(float)*6;

    struct trackparam *tp;
    struct trackparam *tp_d;

    float *hp;
    float *hp_d;

    // Allocate memory
    tp = (struct trackparam*)malloc(TRACK_SIZE);
    hp = (float *)malloc(HELIX_SIZE);
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

    // Cleanup
    free(tp); free(hp);
    cudaFree(tp_d); cudaFree(hp_d);

    // TODO: find real inputdata and test
    float t = 3; 
    float rv[3], d[3], dd[3];
    Evaluate(hp, t, rv, d, dd);


    return 1;
}
