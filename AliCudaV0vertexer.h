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

// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliExternalTrackParam.h.html#178
struct trackparam
{
	Double_t fP[5];
	Double_t fAlpha;
	Double_t fX;
        Double_t fC[15];
        ULong_t fFlags;
};

// Datastructure for v0 secondary vertices
struct v0vertex {
    struct trackparam fParamN;
    struct trackparam fParamP;
    Int_t fNidx;
    Int_t fPidx;
    Double_t fNmom[3];
    Double_t fPmom[3];
    Double_t fPos[3];

    Double_t fDcaV0Daughters;
    Double_t fChi2V0;
};

#define MAX_VERTICES 1000

// Datastructure for the primary vertex
struct privertex {
    Double_t fPosition[3];
    struct v0vertex* fVertices[MAX_VERTICES];
};

Int_t cuda_v0_vertexer(struct privertex* vtx, struct trackparam* tracks, 
                        Int_t ntrks, Double_t b);

int test_cuda_v0_vertexer();
