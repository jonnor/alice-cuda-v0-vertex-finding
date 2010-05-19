
// Global constants
Double_t b = -5.00668; //Magnetic field

// Initialize a track
Double_t x = -0.00462971;
Double_t alpha = 1.90909;
Double_t param[5] = {-0.00429969, -4.56162, 2.38928e-09, 0.855637, -1.96397};
Double_t covar[15] =    {1.0, 2.0, 3.0, 4.0, 5.0, 
                        1.0, 2.0, 3.0, 4.0, 5.0, 
                        1.0, 2.0, 3.0, 4.0, 5.0};
for (int i=0; i<15; i++) covar[i] = 10;

track = new AliExternalTrackParam(x, alpha, param, covar);


// GetHelixParameters()
Double_t helixparam[6];
track->GetHelixParameters(helixparam, b)
for (int i; i<6; i++) printf("%f ", helixparam[i]);


// Attributes
track.GetParameter()
track.GetAlpha()
track.GetX()
track.GetCovariance()


// PropagateTo()
Double_t x = 5.0;
Bool_t success = track->PropagateTo(x, b); // NB: Modifies the track

// GetDCA()
track1 = new AliExternalTrackParam(x, alpha, param, covar);
track2 = new AliExternalTrackParam(x, alpha, param, covar);

Double_t xthis;
Double_t xp;
Double_t dca = track1.GetDCA(track2, b, xthis, xp);


// GetD()
track = new AliExternalTrackParam(x, alpha, param, covar);
Double_t traverse_imp_param = track.GetD(x, y, b);
