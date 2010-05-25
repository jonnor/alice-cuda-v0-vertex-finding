/* This file contains code derived from AliROOT, ported to the CUDA platform
 * Part of a project that implements GPU based v0 vertex finding  
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

// Forward declaration
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


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#1175
// NOTE: I am not sure about the function qualifier
__device__ __host__ Double_t GetDCA(struct trackparam *tp1, struct trackparam *tp2, Double_t b, Double_t &xthis, Double_t &xp) {

  Double_t dy2 = tp1->fC[0] + tp2->fC[0]; //GetSigmaY2
  Double_t dz2 = tp1->fC[2] + tp2->fC[2]; //GetSigmaZ2
  Double_t dx2 = dy2; 

  Double_t p1[8]; GetHelixParameters(tp1,p1,b);
  p1[6] = sin(p1[2]); p1[7] = cos(p1[2]);
  Double_t p2[8]; GetHelixParameters(tp2,p2,b);
  p2[6] = sin(p2[2]); p2[7] = cos(p2[2]);

  Double_t r1[3], g1[3], gg1[3]; Double_t t1 = 0.;
  Evaluate(p1, t1, r1, g1, gg1);//calculates some derivatives and such from p1
  Double_t r2[3], g2[3], gg2[3]; Double_t t2 = 0.;
  Evaluate(p2, t2, r2, g2, gg2);//calculates some derivatives and such from p2

  // some partial derivatives
  Double_t dx = r2[0]-r1[0], dy = r2[1]-r1[1], dz = r2[2]-r1[2];
  Double_t dm = dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;

  Int_t max=27;
  while (max--) {
     Double_t gt1 = -(dx*g1[0]/dx2 + dy*g1[1]/dy2 + dz*g1[2]/dz2);
     Double_t gt2 = +(dx*g2[0]/dx2 + dy*g2[1]/dy2 + dz*g2[2]/dz2);
	 // define hessian matrix
     Double_t h11 = (g1[0]*g1[0] - dx*gg1[0])/dx2 + 
                 (g1[1]*g1[1] - dy*gg1[1])/dy2 +
                 (g1[2]*g1[2] - dz*gg1[2])/dz2;
     Double_t h22 = (g2[0]*g2[0] + dx*gg2[0])/dx2 + 
                 (g2[1]*g2[1] + dy*gg2[1])/dy2 +
                 (g2[2]*g2[2] + dz*gg2[2])/dz2;
     Double_t h12 = -(g1[0]*g2[0]/dx2 + g1[1]*g2[1]/dy2 + g1[2]*g2[2]/dz2);
	//calculate determinant of hessian matrix
     Double_t det = h11*h22-h12*h12;

     Double_t dt1,dt2;
     if ( abs(det)<1.e-33 ) {
        //(quasi)singular Hessian
        dt1=-gt1; dt2=-gt2;
     } else {
        dt1=-(gt1*h22 - gt2*h12)/det; 
        dt2=-(h11*gt2 - h12*gt1)/det;
     }

     if ((dt1*gt1+dt2*gt2)>0) {dt1=-dt1; dt2=-dt2;}

     //check delta(phase1) ?
     //check delta(phase2) ?

     if ( abs(dt1)/(abs(t1)+1.e-3 ) < 1.e-4)
     if ( abs(dt2)/(abs(t2)+1.e-3 ) < 1.e-4) {
        if ( (gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2) 
          ; // 	  printf(" stopped at not a stationary point !"); //AliDebug(1," stopped at not a stationary point !");
        Double_t lmb=h11+h22; lmb=lmb-sqrt(lmb*lmb-4*det);
        if (lmb < 0.) 
          ; // 	  printf(" stopped at not a minimum !"); //AliDebug(1," stopped at not a minimum !");
        break;
     }

     Double_t dd = dm;
     for (Int_t div=1 ; ; div*=2) {
	Evaluate(p1,t1+dt1,r1,g1,gg1);
        Evaluate(p2,t2+dt2,r2,g2,gg2);

        dx = r2[0]-r1[0]; dy = r2[1]-r1[1]; dz = r2[2]-r1[2];
        dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
	if (dd<dm) break;
        dt1*=0.5; dt2*=0.5;
        if (div>512) {
            break; //printf(" overshoot !"); AliDebug(1," overshoot !"); break;
        }
     }
     dm=dd;

     t1+=dt1;
     t2+=dt2;

  } // end of while

  if (max<=0) ; // printf(" too many iterations !"); //AliDebug(1," too many iterations !");

  Double_t cs = cos(tp1->fAlpha);
  Double_t sn = sin(tp1->fAlpha);
  xthis=(r1[0] * cs) + (r1[1] * sn);

  cs=cos(tp2->fAlpha);
  sn=sin(tp2->fAlpha);
  xp = (r2[0] * cs) + (r2[1] * sn);

  return sqrt( dm*sqrt(dy2*dz2) );
}

#define fP tp->fP
#define fC tp->fC
#define fX tp->fX
#define fAlpha tp->fAlpha

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliVParticle.cxx.html#EkM_t
__host__ __device__ Bool_t Local2GlobalMomentum(Double_t p[3], Double_t alpha) {
  if (abs(p[0])<=kAlmost0) return kFALSE;
  if (abs(p[1])> kAlmost1) return kFALSE;

  Double_t pt=1./abs(p[0]);
  Double_t cs=cos(alpha), sn=sin(alpha);
  Double_t r=sqrt((1. - p[1])*(1. + p[1]));
  p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];

  return kTRUE;
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliVParticle.cxx.html#EkM_t
__host__ __device__ Bool_t Local2GlobalPosition(Double_t r[3], Double_t alpha) {
  Double_t cs=cos(alpha), sn=sin(alpha), x=r[0];
  r[0]=x*cs - r[1]*sn; r[1]=x*sn + r[1]*cs;

  return kTRUE;
}


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#WGE_nE
__host__ __device__ Bool_t GetPxPyPz(struct trackparam *tp, Double_t p[3]) {
  p[0]=fP[4]; p[1]=fP[2]; p[2]=fP[3];
  return Local2GlobalMomentum(p,fAlpha);
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#JBAezD
__host__ __device__ Double_t GetX(struct trackparam *tp) {
    return fX;
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#zlsQ3B
__host__ __device__ Bool_t GetXYZ(struct trackparam *tp, Double_t *r) {
  r[0]=fX; r[1]=fP[0]; r[2]=fP[1];
  return Local2GlobalPosition(r,fAlpha);
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#nauVnC
__host__ __device__ Double_t GetSign(struct trackparam *tp) {return (fP[4]>0) ? 1 : -1;}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#XB.FNC
__host__ __device__ Double_t GetC(struct trackparam *tp, Double_t b) {
    return fP[4]*b*kB2C;
}

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
__host__ __device__ void GetHelixParameters(struct trackparam *tp, Double_t hlx[6], Double_t b) {

    Double_t cs=cos(fAlpha), sn=sin(fAlpha); 

    hlx[0]=fP[0]; hlx[1]=fP[1]; hlx[2]=fP[2]; hlx[3]=fP[3];

    hlx[5]=fX*cs - hlx[0]*sn;               // x0
    hlx[0]=fX*sn + hlx[0]*cs;               // y0
    //hlx[1]=                                 // z0
    hlx[2]=sin(hlx[2]) + fAlpha;    // phi0
    //hlx[3]=                                 // tgl
    hlx[4]=GetC(tp, b);                         // C
} 


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#RJz9EE
__host__ __device__ void Evaluate(const Double_t *h, Double_t t,
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


// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.cxx.html#u.xhAD
__host__ __device__ Bool_t PropagateTo(struct trackparam *tp, Double_t xk, Double_t b) {
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

__host__ __device__ Double_t GetLinearD(struct trackparam *tp, Double_t xv,Double_t yv) {
  Double_t sn=sin(fAlpha), cs=cos(fAlpha);
  Double_t x= xv*cs + yv*sn;
  Double_t y=-xv*sn + yv*cs;

  Double_t d = (fX-x)*fP[2] - (fP[0]-y)*sqrt((1.-fP[2])*(1.+fP[2]));

  return -d;
}

__host__ __device__ Double_t GetD(struct trackparam *tp, Double_t x,Double_t y,Double_t b) {
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

#undef fP
#undef fC
#undef fX
#undef fAlpha
