
#ifndef COMMON_H
#define COMMON_H

#include <float.h>
#define sin sin
#define cos cos

#define Float_t float
#define Double_t double
#define Int_t int
#define Bool_t int
#define ULong_t unsigned long

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliVParticle.h.html#18
const Double_t kB2C=-0.299792458e-3;
const Double_t kAlmost0=(Double_t)FLT_MIN;
const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
const Double_t kAlmost0Field=1.e-13;

const int kTRUE=1;
const int kFALSE=0;

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/AliESDtrack.html#AliESDtrack:kTPCrefit
#define kITSrefit 0x0004
#define kTPCrefit 0x0040


#endif
