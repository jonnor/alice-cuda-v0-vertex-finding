outfile = new TXMLFile("esd_tracks.xml", "recreate")
infile = new TFile("AliESDs.root")


TTree *esdtree 
esdtree = (TTree)infile.Get("esdTree")

TBranch *tracks = esdtree.GetBranch("Tracks")
TLeaf *fP = esdtree.GetLeaf("Tracks.fP")

//Prints out all the info about entry nr 1 in the tree
esdtree.Show(1) 

esdtree.Scan("Tracks.fP")
esdtree.Scan("Tracks.fP[0]")

// Gives a segfault :(
esdtree->SetScanField(0)
esdtree->Scan("*"); > esdtree.txt

// Does not work with XML files? :(
outfile.cd()
esdtree->Write()
