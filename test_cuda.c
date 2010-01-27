

#include <stdio.h>
#include <math.h>


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliVParticle.h.html#18
const double kB2C=-0.299792458e-3;


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#178
// We could also just use a simple array instead of the struct
// double trackparam[7]
// When we start caring about multiple tracks, they can either be represented by
// - An array of structs/arrays
// - An array where the params are separated by a rowstride equal to number of tracks
// The latter might perform better, but be less readable
struct trackparam
{
	double fP[5];
	double fAlpha;
	double fX;
};

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#XB.FNC
double GetC(struct trackparam *tp, double b) {
    return tp->fP[4]*b*kB2C;
}


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
void GetHelixParameters(struct trackparam *tp, double hlx[6], double b) {

    double cs=cos(tp->fAlpha), sn=sin(tp->fAlpha); 

    hlx[0]=tp->fP[0]; hlx[1]=tp->fP[1]; hlx[2]=tp->fP[2]; hlx[3]=tp->fP[3];

    hlx[5]=tp->fX*cs - hlx[0]*sn;               // x0
    hlx[0]=tp->fX*sn + hlx[0]*cs;               // y0
    //hlx[1]=                                 // z0
    hlx[2]=sin(hlx[2]) + tp->fAlpha;    // phi0
    //hlx[3]=                                 // tgl
    hlx[4]=GetC(tp, b);                         // C
} 


int main()
{

    double b = -5.00668;
    struct trackparam *tp;
    tp->fP[0] = -0.00429969;
    tp->fP[1] = -4.56162;
    tp->fP[2] = 2.38928e-09;
    tp->fP[3] = 0.855637;
    tp->fP[4] = -1.96397;
    tp->fAlpha = 1.90909;
    tp->fX = -0.00462971;

    double hp[6];
    GetHelixParameters(tp, hp, b);
    printf("GetHelixParameters\n");
    for (int i; i<6; i++) printf("%f\n", hp[i]);

    return 1;
}
