
#include "aliexternaltrackparam.h"

/*
// TODO. implement
// TODO: Needs to go into a file in our cuda project code, which we canz link against
Int_t cuda_v0vertexer(struct trackparam* tracks) {
    pass;
}

// TODO: implement, lower priority
void add_v0s_to_event(AliESDEvent *event) {
    ;
}
*/

/*
// TODO: Needs to go into AliV0vertexer.cxx
// http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliV0vertexer.cxx.html#gwutwE
Int_t Tracks2v0vertices_gpu(AliESDEvent *event) {

    // Extract relevant data and convert to a cuda friendly format
    Double_t b=event->GetMagneticField();

    Double_t pos[3];
    AliESDVertex *priVertex = event->GetPrimaryVertex();
    priVertex->GetXYZ(pos);

    Int_t ntrks;
    struct trackparam* tracks = tracks_from_event(event, ntrks);

    // Do the actual vertexfinding
    cuda_v0vertexer(tracks);
    free(tracks);

    // Put the data into the event
    add_v0s_to_event(event);
    Int_t nv0s;
    return nv0s
}
*/

struct trackparam* tracks_from_event(AliESDEvent *event, Int_t &ntrks) {
    ntrks = event->GetNumberOfTracks();
    struct trackparam* tracks;
    tracks = (struct trackparam*)malloc(sizeof(struct trackparam)*ntrks);

    // Fill in the data
    for(Int_t i=0; i<ntrks; i++) {
        struct trackparam* tp=&tracks[i];
        AliESDtrack *esdTrack=event->GetTrack(i);

        tp->fAlpha = esdTrack->GetAlpha();
        esdTrack->GetExternalCovariance(tp->fC);
        esdTrack->GetExternalParameters(tp->fX, tp->fP);
    }
    return tracks;
}

