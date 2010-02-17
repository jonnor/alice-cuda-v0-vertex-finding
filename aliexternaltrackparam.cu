
#include "common.h"
#include "aliexternaltrackparam.h"

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#nauVnC
Double_t GetSign(struct trackparam *tp) {return (tp->fP[4]>0) ? 1 : -1;}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#XB.FNC
__host__ __device__ Double_t GetC(struct trackparam *tp, Double_t b) {
    return tp->fP[4]*b*kB2C;
}


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
__global__ void GetHelixParameters(struct trackparam *tp, Double_t hlx[6], Double_t b) {

    Double_t cs=cos(tp->fAlpha), sn=sin(tp->fAlpha); 

    hlx[0]=tp->fP[0]; hlx[1]=tp->fP[1]; hlx[2]=tp->fP[2]; hlx[3]=tp->fP[3];

    hlx[5]=tp->fX*cs - hlx[0]*sn;               // x0
    hlx[0]=tp->fX*sn + hlx[0]*cs;               // y0
    //hlx[1]=                                 // z0
    hlx[2]=sin(hlx[2]) + tp->fAlpha;    // phi0
    //hlx[3]=                                 // tgl
    hlx[4]=GetC(tp, b);                         // C
} 


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
void Evaluate(const Double_t *h, Double_t t,
                     Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
                     Double_t gg[3]) //second derivatives
{

  Double_t phase=h[4]*t+h[2];
  Double_t sn=sin(phase), cs=cos(phase);

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

//FIXME: move these to the top of file
#define fP tp->fP
#define fC tp->fC
#define fX tp->fX
#define fAlpha tp->fAlpha

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#u.xhAD
Bool_t PropagateTo(struct trackparam *tp, Double_t xk, Double_t b) {
  Double_t dx=xk-fX;
  if (abs(dx)<=kAlmost0)  return kTRUE;

  Double_t crv=GetC(tp, b);
  if (abs(b) < kAlmost0Field) crv=0.;

  Double_t f1=fP[2], f2=f1 + crv*dx;
  if (abs(f1) >= kAlmost1) return kFALSE;
  if (abs(f2) >= kAlmost1) return kFALSE;

  Double_t &fP0=fP[0], &fP1=fP[1], &fP2=fP[2], &fP3=fP[3], &fP4=fP[4];
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];

  Double_t r1=sqrt((1.-f1)*(1.+f1)), r2=sqrt((1.-f2)*(1.+f2));

  fX=xk;
  fP0 += dx*(f1+f2)/(r1+r2);
  fP1 += dx*(r2 + f2*(f1+f2)/(r1+r2))*fP3;  // Many thanks to P.Hristov !
  fP2 += dx*crv;

  //f = F - 1

  Double_t f02=    dx/(r1*r1*r1);            Double_t cc=crv/fP4;
  Double_t f04=0.5*dx*dx/(r1*r1*r1);         f04*=cc;
  Double_t f12=    dx*fP3*f1/(r1*r1*r1);
  Double_t f14=0.5*dx*dx*fP3*f1/(r1*r1*r1);  f14*=cc;
  Double_t f13=    dx/r1;
  Double_t f24=    dx;                       f24*=cc;

  //b = C*ft
  Double_t b00=f02*fC20 + f04*fC40, b01=f12*fC20 + f14*fC40 + f13*fC30;
  Double_t b02=f24*fC40;
  Double_t b10=f02*fC21 + f04*fC41, b11=f12*fC21 + f14*fC41 + f13*fC31;
  Double_t b12=f24*fC41;
  Double_t b20=f02*fC22 + f04*fC42, b21=f12*fC22 + f14*fC42 + f13*fC32;
  Double_t b22=f24*fC42;
  Double_t b40=f02*fC42 + f04*fC44, b41=f12*fC42 + f14*fC44 + f13*fC43;
  Double_t b42=f24*fC44;
  Double_t b30=f02*fC32 + f04*fC43, b31=f12*fC32 + f14*fC43 + f13*fC33;
  Double_t b32=f24*fC43;

  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f24*b42;

  //F*C*Ft = C + (b + bt + a)
  fC00 += b00 + b00 + a00;
  fC10 += b10 + b01 + a01; 
  fC20 += b20 + b02 + a02;
  fC30 += b30;
  fC40 += b40;
  fC11 += b11 + b11 + a11;
  fC21 += b21 + b12 + a12;
  fC31 += b31; 
  fC41 += b41;
  fC22 += b22 + b22 + a22;
  fC32 += b32;
  fC42 += b42;

  // CheckCovariance(); //TODO: implement

  return kTRUE;
}

Double_t GetLinearD(struct trackparam *tp, Double_t xv,Double_t yv) {
  Double_t sn=sin(fAlpha), cs=cos(fAlpha);
  Double_t x= xv*cs + yv*sn;
  Double_t y=-xv*sn + yv*cs;

  Double_t d = (fX-x)*fP[2] - (fP[0]-y)*sqrt((1.-fP[2])*(1.+fP[2]));

  return -d;
}

Double_t GetD(struct trackparam *tp, Double_t x,Double_t y,Double_t b) {
  if (abs(b) < kAlmost0Field) return GetLinearD(tp,x,y);
  Double_t rp4=GetC(tp, b);

  Double_t xt=fX, yt=fP[0];

  Double_t sn=sin(fAlpha), cs=cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  sn=rp4*xt - fP[2]; cs=rp4*yt + sqrt((1.- fP[2])*(1.+fP[2]));
  a=2*(xt*fP[2] - yt*sqrt((1.-fP[2])*(1.+fP[2])))-rp4*(xt*xt + yt*yt);
  return  -a/(1 + sqrt(sn*sn + cs*cs));
}

//FIXME: put functions above into separate file so that we dont have to do this
#undef fP
#undef fC
#undef fAlpha
#undef fX
