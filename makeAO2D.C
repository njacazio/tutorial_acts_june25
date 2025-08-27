#include "TreeReader.h"
#include "dm_O2track_iu.h"
#include "ReconstructionDataFormats/Track.h"

using O2Track = o2::track::TrackParCov;

void makeO2Track(O2Track &o2track, std::array<float, 3> xyz, std::array<float, 3> ptetaphi, int charge)
{
    std::array<float, 5> params;
    std::array<float, 15> covm = {
        0.,
        0., 0.,
        0., 0., 0.,
        0., 0., 0., 0.,
        0., 0., 0., 0., 0.};
    float s, c, x;
    o2::math_utils::sincos(ptetaphi[2], s, c);
    o2::math_utils::rotateZInv(xyz[0], xyz[1], x, params[0], s, c);
    params[1] = xyz[2];
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * atan(exp(-ptetaphi[1]));
    params[3] = 1. / tan(theta);
    params[4] = charge / ptetaphi[0];

    new (&o2track) O2Track(x, ptetaphi[2], params, covm);
}

struct Track{
    float X;
    float Y;
    float Z;
    float Xd;
    float Yd;
    float Zd;
    float PT;
    float Eta;
    float Phi;
    int Charge;
};

void convertTrackToO2Track(const Track &track, O2Track &o2track, bool atDCA)
{
    std::array<float, 3> xyz = {0.1f * static_cast<float>(atDCA ? track.Xd : track.X),
                                0.1f * static_cast<float>(atDCA ? track.Yd : track.Y),
                                0.1f * static_cast<float>(atDCA ? track.Zd : track.Z)};

    std::array<float, 3> ptetaphi = {static_cast<float>(track.PT), static_cast<float>(track.Eta), static_cast<float>(track.Phi)};
    int charge = track.Charge;
    makeO2Track(o2track, xyz, ptetaphi, charge);
}

void makeAO2D(std::string inputFile = "../reco_output_pythia/tracksummary_ambi.root")
{
    O2track_iu track_iu;
    TFile *inFile = TFile::Open(inputFile.c_str());
    TFile *outFile = TFile::Open("AO2D.root", "RECREATE");
    outFile->mkdir("DF_001");
    outFile->cd("DF_001");
    TTree *tree = (TTree *)inFile->Get("tracksummary");
    TrackSummaryReader tsReader(tree, false);

    Track ttt;
    O2Track o2ttt;
    for (unsigned iEntry = 0; iEntry < tree->GetEntries(); iEntry++)
    {
        tsReader.getEntry(iEntry);

        // Select only fitted tracks
        if (not tsReader.hasFittedParams)
            continue;

        // Loop over tracks in the event
        for (unsigned ip = 0; ip < tsReader.eQOP_fit->size(); ip++)
        {
            track_iu.fIndexCollisions = iEntry;
            track_iu.fTrackType = 1;
            ttt.X = tsReader.t_d0->at(ip);
            ttt.Y = tsReader.t_d0->at(ip);
            ttt.Z = tsReader.t_z0->at(ip);
            ttt.PT = tsReader.t_pT->at(ip);
            ttt.Eta = tsReader.t_eta->at(ip);
            ttt.Phi = tsReader.t_phi->at(ip);
            ttt.Charge = tsReader.t_charge->at(ip);
            convertTrackToO2Track(ttt, o2ttt, false);
            track_iu.fX = o2ttt.getX();
            track_iu.fAlpha = o2ttt.getAlpha();
            track_iu.fZ = o2ttt.getZ();
            track_iu.fSnp = o2ttt.getSnp();
            track_iu.fTgl = o2ttt.getTgl();
            track_iu.fSigned1Pt = o2ttt.getCharge2Pt();
            track_iu.fill();
        }
    }

    track_iu.write();
    outFile->Close();
}
