
esdTree
Int_t nevents = (Int_t)esdTree->GetEntries()
for (Int_t i=0; i<nevents; i++) esdTree->Show(i,1e9); >> events.txt

.L ../alice1/macrotest.C
vertexer = new AliV0vertexer()
event = load_event("/home/lbratrud/Jon/1/AliESDs.root",0)
vertexer->Tracks2V0vertices(event)
vertexer->Tracks2V0verticesCuda(event)
