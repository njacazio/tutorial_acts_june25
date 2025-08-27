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
"""


def main(input_ao2d_name="/home/njacazio/cernbox/pidtest/AO2D.root", tree_list=["O2track_iu"]):

    # Cleanup DATA_MODEL, keep only lines with DECLARE_SOA_INDEX_COLUMN or DECLARE_SOA_COLUMN
    cleaned_data_model = "\n".join(
        line for line in DATA_MODEL.splitlines()
        if "DECLARE_SOA_INDEX_COLUMN" in line or "DECLARE_SOA_COLUMN" in line
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
            name = "Index"+name + "s"
            t = "uint64_t"
        types[f"f{name}"] = t

    print(types)

    from ROOT import TFile
    input_ao2d = TFile.Open(input_ao2d_name, "READ")
    dfname = input_ao2d.GetListOfKeys().At(0).GetName()
    print(dfname)
    if "DF_" not in dfname:
        raise ValueError("Input file is not a DataFrame")
    df = input_ao2d.Get(dfname)
    # df.ls()
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
                f.write(f"    {types[branch.GetName()]} {branch.GetName()};\n")
            f.write("    TTree *fTree;\n")
            # Constructor
            f.write("    "+i+"() {\n")
            f.write("        fTree = new TTree(\""+i+"\", \"\");\n")
            for branch in tree.GetListOfBranches():
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


main()
