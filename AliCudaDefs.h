
#ifndef COMMON_H
#define COMMON_H

#include <float.h>
#define sin sinf
#define cos cosf

#define Float_t float
#define Double_t double
#define Int_t int
#define Bool_t int

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliVParticle.h.html#18
const Double_t kB2C=-0.299792458e-3;
const Double_t kAlmost0=(Double_t)FLT_MIN;
const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
const Double_t kAlmost0Field=1.e-13;

const int kTRUE=1;
const int kFALSE=0;

#endif
