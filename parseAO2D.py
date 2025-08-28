#! /usr/bin/env python3

"""Script to parse the AO2D data model for the interesting trees and write a header with the correct data model to fill the TTrees"""

DATA_MODEL = """
DECLARE_SOA_INDEX_COLUMN(Collision, collision);    //! Collision to which this track belongs
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t); //! Type of track. See enum TrackTypeEnum. This cannot be used to decide which detector has contributed to this track. Use hasITS, hasTPC, etc.
DECLARE_SOA_COLUMN(X, x, float);                   //!
DECLARE_SOA_COLUMN(Alpha, alpha, float);           //!
DECLARE_SOA_COLUMN(Y, y, float);                   //!
DECLARE_SOA_COLUMN(Z, z, float);                   //!
DECLARE_SOA_COLUMN(Snp, snp, float);               //!
DECLARE_SOA_COLUMN(Tgl, tgl, float);               //!
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);   //! (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float,     //! Phi of the track, in radians within [0, 2pi)
                              ifnode(nasin(aod::track::snp) + aod::track::alpha < 0.0f, nasin(aod::track::snp) + aod::track::alpha + o2::constants::math::TwoPI,
                                     ifnode(nasin(aod::track::snp) + aod::track::alpha >= o2::constants::math::TwoPI, nasin(aod::track::snp) + aod::track::alpha - o2::constants::math::TwoPI,
                                            nasin(aod::track::snp) + aod::track::alpha)));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity
                              -1.f * nlog(ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::track::tgl))));
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum of the track in GeV/c
                              ifnode(nabs(aod::track::signed1Pt) <= o2::constants::math::Almost0, o2::constants::math::VeryBig, nabs(1.f / aod::track::signed1Pt)));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Momentum in Gev/c
                              ifnode(nabs(aod::track::signed1Pt) <= o2::constants::math::Almost0, o2::constants::math::VeryBig, 0.5f * (ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::track::tgl)) + 1.f / ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::track::tgl))) / nabs(aod::track::signed1Pt)));


DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);             //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);               //! Impact parameter in Z of the track to the primary vertex
DECLARE_SOA_COLUMN(SigmaDcaXY2, sigmaDcaXY2, float); //! Impact parameter sigma^2 in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(SigmaDcaZ2, sigmaDcaZ2, float);   //! Impact parameter sigma^2 in Z of the track to the primary vertex

DECLARE_SOA_COLUMN(SigmaY, sigmaY, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaZ, sigmaZ, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaSnp, sigmaSnp, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaTgl, sigmaTgl, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(Sigma1Pt, sigma1Pt, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(RhoZY, rhoZY, int8_t);         //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoSnpY, rhoSnpY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoSnpZ, rhoSnpZ, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglY, rhoTglY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglZ, rhoTglZ, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglSnp, rhoTglSnp, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtY, rho1PtY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtZ, rho1PtZ, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtSnp, rho1PtSnp, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtTgl, rho1PtTgl, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_EXPRESSION_COLUMN(CYY, cYY, float, //!
                              aod::track::sigmaY* aod::track::sigmaY);
DECLARE_SOA_EXPRESSION_COLUMN(CZY, cZY, float, //!
                              (aod::track::rhoZY / 128.f) * (aod::track::sigmaZ * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CZZ, cZZ, float, //!
                              aod::track::sigmaZ* aod::track::sigmaZ);
DECLARE_SOA_EXPRESSION_COLUMN(CSnpY, cSnpY, float, //!
                              (aod::track::rhoSnpY / 128.f) * (aod::track::sigmaSnp * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CSnpZ, cSnpZ, float, //!
                              (aod::track::rhoSnpZ / 128.f) * (aod::track::sigmaSnp * aod::track::sigmaZ));
DECLARE_SOA_EXPRESSION_COLUMN(CSnpSnp, cSnpSnp, float, //!
                              aod::track::sigmaSnp* aod::track::sigmaSnp);
DECLARE_SOA_EXPRESSION_COLUMN(CTglY, cTglY, float, //!
                              (aod::track::rhoTglY / 128.f) * (aod::track::sigmaTgl * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CTglZ, cTglZ, float, //!
                              (aod::track::rhoTglZ / 128.f) * (aod::track::sigmaTgl * aod::track::sigmaZ));
DECLARE_SOA_EXPRESSION_COLUMN(CTglSnp, cTglSnp, float, //!
                              (aod::track::rhoTglSnp / 128.f) * (aod::track::sigmaTgl * aod::track::sigmaSnp));
DECLARE_SOA_EXPRESSION_COLUMN(CTglTgl, cTglTgl, float, //!
                              aod::track::sigmaTgl* aod::track::sigmaTgl);
DECLARE_SOA_EXPRESSION_COLUMN(C1PtY, c1PtY, float, //!
                              (aod::track::rho1PtY / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtZ, c1PtZ, float, //!
                              (aod::track::rho1PtZ / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaZ));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtSnp, c1PtSnp, float, //!
                              (aod::track::rho1PtSnp / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaSnp));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtTgl, c1PtTgl, float, //!
                              (aod::track::rho1PtTgl / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaTgl));
DECLARE_SOA_EXPRESSION_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float, //!
                              aod::track::sigma1Pt* aod::track::sigma1Pt);


DECLARE_SOA_COLUMN(IsReconstructed, isReconstructed, bool); //! is reconstructed or not
DECLARE_SOA_COLUMN(NSiliconHits, nSiliconHits, int);        //! number of silicon hits
DECLARE_SOA_COLUMN(NTPCHits, nTPCHits, int);                //! number of tpc hits

DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);     //! Bit mask to indicate detector mismatches (bit ON means mismatch). Bit 0-6: mismatch at ITS layer. Bit 7-9: # of TPC mismatches in the ranges 0, 1, 2-3, 4-7, 8-15, 16-31, 32-63, >64. Bit 10: TRD, bit 11: TOF, bit 15: indicates negative label

DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);       //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label

DECLARE_SOA_INDEX_COLUMN(BC, bc);                              //! Most probably BC to where this collision has occured
DECLARE_SOA_COLUMN(PosX, posX, float);                         //! X Vertex position in cm
DECLARE_SOA_COLUMN(PosY, posY, float);                         //! Y Vertex position in cm
DECLARE_SOA_COLUMN(PosZ, posZ, float);                         //! Z Vertex position in cm
DECLARE_SOA_COLUMN(CovXX, covXX, float);                       //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovXY, covXY, float);                       //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovXZ, covXZ, float);                       //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovYY, covYY, float);                       //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovYZ, covYZ, float);                       //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovZZ, covZZ, float);                       //! Vertex covariance matrix
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);                    //! Run 2: see CollisionFlagsRun2 | Run 3: see Vertex::Flags
DECLARE_SOA_COLUMN(Chi2, chi2, float);                         //! Chi2 of vertex fit
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t);          //! Number of tracks used for the vertex
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);       //! Collision time in ns relative to BC stored in bc()
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float); //! Resolution of collision time

DECLARE_SOA_COLUMN(MultDensity, multDensity, float); //! midrapidity Nch density


DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);                                      //! Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(Flags, flags, uint32_t);                                                   //! Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t);                               //! Clusters sizes, four bits per a layer, starting from the innermost
DECLARE_SOA_COLUMN(ITSClusterMap, itsClusterMap, uint8_t);                                    //! Old cluster ITS cluster map, kept for version 0 compatibility
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, uint8_t);                                //! Findable TPC clusters for this track geometry
DECLARE_SOA_COLUMN(TPCNClsFindableMinusFound, tpcNClsFindableMinusFound, int8_t);             //! TPC Clusters: Findable - Found
DECLARE_SOA_COLUMN(TPCNClsFindableMinusPID, tpcNClsFindableMinusPID, int8_t);                 //! TPC Clusters: Findable - Found clusters used for PID
DECLARE_SOA_COLUMN(TPCNClsFindableMinusCrossedRows, tpcNClsFindableMinusCrossedRows, int8_t); //! TPC Clusters: Findable - crossed rows
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t);                                    //! Number of shared TPC clusters
DECLARE_SOA_COLUMN(ITSSharedClusterMap, itsSharedClusterMap, uint8_t);                        //! shared ITS cluster map (Run 2)
DECLARE_SOA_COLUMN(TRDPattern, trdPattern, uint8_t);                                          //! Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);                                            //! Chi2 / cluster for the ITS track segment
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);                                            //! Chi2 / cluster for the TPC track segment
DECLARE_SOA_COLUMN(TRDChi2, trdChi2, float);                                                  //! Chi2 for the TRD track segment
DECLARE_SOA_COLUMN(TOFChi2, tofChi2, float);                                                  //! Chi2 for the TOF track segment
DECLARE_SOA_COLUMN(ITSSignal, itsSignal, float);                                              //! dE/dx signal in the ITS (Run 2)
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                              //! dE/dx signal in the TPC
DECLARE_SOA_COLUMN(TRDSignal, trdSignal, float);                                              //! PID signal in the TRD
DECLARE_SOA_COLUMN(Length, length, float);                                                    //! Track length
DECLARE_SOA_COLUMN(TOFExpMom, tofExpMom, float);                                              //! TOF expected momentum obtained in tracking, used to compute the expected times
DECLARE_SOA_COLUMN(TrackEtaEMCAL, trackEtaEmcal, float);                                      //!
DECLARE_SOA_COLUMN(TrackPhiEMCAL, trackPhiEmcal, float);                                      //!
DECLARE_SOA_COLUMN(TrackTime, trackTime, float);                                              //! Estimated time of the track in ns wrt collision().bc() or ambiguoustrack.bcSlice()[0]
DECLARE_SOA_COLUMN(TrackTimeRes, trackTimeRes, float);                                        //! Resolution of the track time in ns (see TrackFlags::TrackTimeResIsRange)

// expression columns changing between versions have to be declared in different namespaces

DECLARE_SOA_EXPRESSION_COLUMN(DetectorMap, detectorMap, uint8_t, //! Detector map: see enum DetectorMapEnum
                              ifnode(aod::track::itsClusterMap > (uint8_t)0, static_cast<uint8_t>(o2::aod::track::ITS), (uint8_t)0x0) |
                                ifnode(aod::track::tpcNClsFindable > (uint8_t)0, static_cast<uint8_t>(o2::aod::track::TPC), (uint8_t)0x0) |
                                ifnode(aod::track::trdPattern > (uint8_t)0, static_cast<uint8_t>(o2::aod::track::TRD), (uint8_t)0x0) |
                                ifnode((aod::track::tofChi2 >= 0.f) && (aod::track::tofExpMom > 0.f), static_cast<uint8_t>(o2::aod::track::TOF), (uint8_t)0x0));

                                
DECLARE_SOA_INDEX_COLUMN(BC, bc);                            //! BC index
DECLARE_SOA_COLUMN(GeneratorsID, generatorsID, short);       //! disentangled generator IDs should be accessed using getGeneratorId, getSubGeneratorId and getSourceId
DECLARE_SOA_COLUMN(PosX, posX, float);                       //! X vertex position in cm
DECLARE_SOA_COLUMN(PosY, posY, float);                       //! Y vertex position in cm
DECLARE_SOA_COLUMN(PosZ, posZ, float);                       //! Z vertex position in cm
DECLARE_SOA_COLUMN(T, t, float);                             //! Collision time relative to given bc in ns
DECLARE_SOA_COLUMN(Weight, weight, float);                   //! MC weight
DECLARE_SOA_COLUMN(ImpactParameter, impactParameter, float); //! Impact parameter for A-A
DECLARE_SOA_COLUMN(EventPlaneAngle, eventPlaneAngle, float); //! Event plane angle for A-A

DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);                                     //! MC collision of this particle
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);                                              //! PDG code
DECLARE_SOA_COLUMN(StatusCode, statusCode, int);                                        //! Generators status code or physics process. Do not use directly. Use dynamic columns getGenStatusCode() or getProcess()
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);                                              //! ALICE specific flags, see MCParticleFlags. Do not use directly. Use the dynamic columns, e.g. producedByGenerator()
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother0, mother0, int, "McParticles_Mother0");       //! Track index of the first mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother1, mother1, int, "McParticles_Mother1");       //! Track index of the last mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter0, daughter0, int, "McParticles_Daughter0"); //! Track index of the first daugther
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter1, daughter1, int, "McParticles_Daughter1"); //! Track index of the last daugther
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);                                  //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);                              //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_COLUMN(Weight, weight, float);                                              //! MC weight
DECLARE_SOA_COLUMN(Px, px, float);                                                      //! Momentum in x in GeV/c
DECLARE_SOA_COLUMN(Py, py, float);                                                      //! Momentum in y in GeV/c
DECLARE_SOA_COLUMN(Pz, pz, float);                                                      //! Momentum in z in GeV/c
DECLARE_SOA_COLUMN(E, e, float);                                                        //! Energy
DECLARE_SOA_COLUMN(Vx, vx, float);                                                      //! X production vertex in cm
DECLARE_SOA_COLUMN(Vy, vy, float);                                                      //! Y production vertex in cm
DECLARE_SOA_COLUMN(Vz, vz, float);                                                      //! Z production vertex in cm
DECLARE_SOA_COLUMN(Vt, vt, float);                                                      //! Production time
DECLARE_SOA_DYNAMIC_COLUMN(ProducedByGenerator, producedByGenerator,                    //! True if particle produced by the generator (==TMCProcess::kPrimary); False if by the transport code
                           [](uint8_t flags) -> bool { return (flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0; });
DECLARE_SOA_DYNAMIC_COLUMN(FromBackgroundEvent, fromBackgroundEvent, //! Particle from background event
                           [](uint8_t flags) -> bool { return (flags & o2::aod::mcparticle::enums::FromBackgroundEvent) == o2::aod::mcparticle::enums::FromBackgroundEvent; });
DECLARE_SOA_DYNAMIC_COLUMN(GetProcess, getProcess, //! The VMC physics code (as int) that generated this particle (see header TMCProcess.h in ROOT)
                           [](uint8_t flags, int statusCode) -> int { if ((flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0) { return 0 /*TMCProcess::kPrimary*/; } else { return statusCode; } });
DECLARE_SOA_DYNAMIC_COLUMN(GetGenStatusCode, getGenStatusCode, //! The native status code put by the generator, or -1 if a particle produced during transport
                           [](uint8_t flags, int statusCode) -> int { if ((flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0) { return o2::mcgenstatus::getGenStatusCode(statusCode); } else { return -1; } });
DECLARE_SOA_DYNAMIC_COLUMN(GetHepMCStatusCode, getHepMCStatusCode, //! The HepMC status code put by the generator, or -1 if a particle produced during transport
                           [](uint8_t flags, int statusCode) -> int { if ((flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0) { return o2::mcgenstatus::getHepMCStatusCode(statusCode); } else { return -1; } });
DECLARE_SOA_DYNAMIC_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, //! True if particle is considered a physical primary according to the ALICE definition
                           [](uint8_t flags) -> bool { return (flags & o2::aod::mcparticle::enums::PhysicalPrimary) == o2::aod::mcparticle::enums::PhysicalPrimary; });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector, //! Momentum vector in x,y,z-directions in GeV/c
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array<float, 3>{px, py, pz}; });

DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float, //! Phi in the range [0, 2pi)
                              o2::constants::math::PI + natan2(-1.0f * aod::mcparticle::py, -1.0f * aod::mcparticle::px));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity, conditionally defined to avoid FPEs
                              ifnode((nsqrt(aod::mcparticle::px * aod::mcparticle::px +
                                            aod::mcparticle::py * aod::mcparticle::py +
                                            aod::mcparticle::pz * aod::mcparticle::pz) -
                                      aod::mcparticle::pz) < static_cast<float>(1e-7),
                                     ifnode(aod::mcparticle::pz < 0.f, -100.f, 100.f),
                                     0.5f * nlog((nsqrt(aod::mcparticle::px * aod::mcparticle::px +
                                                        aod::mcparticle::py * aod::mcparticle::py +
                                                        aod::mcparticle::pz * aod::mcparticle::pz) +
                                                  aod::mcparticle::pz) /
                                                 (nsqrt(aod::mcparticle::px * aod::mcparticle::px +
                                                        aod::mcparticle::py * aod::mcparticle::py +
                                                        aod::mcparticle::pz * aod::mcparticle::pz) -
                                                  aod::mcparticle::pz))));
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum in GeV/c
                              nsqrt(aod::mcparticle::px* aod::mcparticle::px +
                                    aod::mcparticle::py * aod::mcparticle::py));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Total momentum in GeV/c
                              nsqrt(aod::mcparticle::px* aod::mcparticle::px +
                                    aod::mcparticle::py * aod::mcparticle::py +
                                    aod::mcparticle::pz * aod::mcparticle::pz));
DECLARE_SOA_EXPRESSION_COLUMN(Y, y, float, //! Particle rapidity, conditionally defined to avoid FPEs
                              ifnode((aod::mcparticle::e - aod::mcparticle::pz) < static_cast<float>(1e-7),
                                     ifnode(aod::mcparticle::pz < 0.f, -100.f, 100.f),
                                     0.5f * nlog((aod::mcparticle::e + aod::mcparticle::pz) /
                                                 (aod::mcparticle::e - aod::mcparticle::pz))));


"""


def main(input_ao2d_name: str, tree_list: list):

    # Cleanup DATA_MODEL, keep only lines with DECLARE_SOA_INDEX_COLUMN or DECLARE_SOA_COLUMN
    cleaned_data_model = "\n".join(
        line for line in DATA_MODEL.splitlines()
        if "DECLARE_SOA_INDEX_COLUMN" in line or "DECLARE_SOA_COLUMN" in line or "DECLARE_SOA_EXPRESSION_COLUMN" in line or "DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN" in line or "DECLARE_SOA_SELF_SLICE_INDEX_COLUMN" in line
    )
    # Strip every line after ;
    cleaned_data_model = [line.split(";")[0] for line in cleaned_data_model.splitlines()]
    types = {}
    for i in cleaned_data_model:
        name = None
        t = None
        if "DECLARE_SOA_COLUMN(" in i:
            name = i.split("(")[1].split(",")[0].strip()
            t = i.split(",")[-1].strip(")").strip()
        elif "DECLARE_SOA_INDEX_COLUMN(" in i:
            name = i.split("(")[1].split(",")[0].strip()
            name = "Index" + name + "s"
            t = "uint64_t"
        elif "DECLARE_SOA_EXPRESSION_COLUMN(" in i:
            name = i.split("(")[1].split(",")[0].strip()
            t = None
        elif "DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(" in i:
            name = i.split("(")[1].split(",")[0].strip()
            name = f"IndexArray_{name}"
            t = "std::vector<uint64_t>"
            types[f"f{name}_size"] = "int"
        elif "DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(" in i:
            name = i.split("(")[1].split(",")[0].strip()
            name = f"IndexSlice_{name}"
            t = "std::pair<uint64_t, uint64_t>"
        types[f"f{name}"] = t

    print(types)

    from ROOT import TFile
    input_ao2d = TFile.Open(input_ao2d_name, "READ")
    dfname = input_ao2d.GetListOfKeys().At(0).GetName()
    print(dfname)
    if "DF_" not in dfname:
        raise ValueError("Input file is not a DataFrame")
    df = input_ao2d.Get(dfname)
    df.ls()
    for i in tree_list:
        tree = df.Get(i)
        if not tree:
            raise ValueError(f"Tree {i} not found in DataFrame")
        print(f"Found tree {i} with {tree.GetEntries()} entries")
        # We now write a class with the same data members
        with open(f"dm_{i}.h", "w") as f:
            f.write("#include <TTree.h>\n")
            f.write("struct "+i+"{\n")
            for branch in tree.GetListOfBranches():
                if branch.GetName() not in types:
                    print("->>>>>", branch.GetName(), branch.GetTitle())
                if types[branch.GetName()] is not None:
                    f.write(f"    {types[branch.GetName()]} {branch.GetName()};\n")
            f.write("    TTree *fTree;\n")
            # Constructor
            f.write("    "+i+"() {\n")
            f.write("        fTree = new TTree(\""+i+"\", \"\");\n")
            for branch in tree.GetListOfBranches():
                t = types[branch.GetName()]
                if t is None:
                    continue
                f.write(
                    f"        fTree->Branch(\"{branch.GetName()}\", &{branch.GetName()}, \"{branch.GetTitle()}\");\n")
            f.write("    }\n")
            # Destructor
            f.write("    ~"+i+"() {\n")
            f.write("        delete fTree;\n")
            f.write("    }\n")
            # Fill method
            f.write("    void fill() {\n")
            f.write("        fTree->Fill();\n")
            f.write("    }\n")
            # Write method
            f.write("    void write() {\n")
            f.write("        fTree->Write();\n")
            f.write("    }\n")
            f.write("};\n")


if 0:
    main(input_ao2d_name="/tmp/AnalysisResults_trees.root", tree_list=["O2trackdca",
                                                                       "O2extrackcov",
                                                                       "O2extrack",
                                                                       "O2tracksextraa3",
                                                                       "O2tracksalice3",
                                                                       "O2trackdcacov",
                                                                       "O2trackcov",
                                                                       "O2track",
                                                                       "O2mctracklabel",
                                                                       "O2mccollislabel",
                                                                       "O2collision_001",
                                                                       "O2collalice3",
                                                                       "O2extrackextra_002"])


if 1:
    main(input_ao2d_name="/tmp/AnalysisResults_trees_pp.root", tree_list=["O2mccollision_001",
                                                                          "O2mcparticle_001"])
