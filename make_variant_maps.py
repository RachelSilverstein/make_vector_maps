import Bio.Seq
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


inputs = pd.read_csv("inputs.csv")

def format_AAs(AAs):
    AAs = "".join(AAs.split())
    AAs = AAs.split(",")
    return AAs

def format_positions(positions):
    positions = "".join(positions.split())
    positions = positions.split(",")
    positions = [int(x) for x in positions]
    return positions

def mutate_codon(sequence, annotation_start_pos, codon, codon_pos):
    nt_start = annotation_start_pos + (3 * codon_pos) - 3
    for i in range(0, 3):
        sequence[nt_start + i] = codon[i]

def get_codon(AA, codon_table):
    codon = codon_table["codon"][AA]
    return codon


for index, row in inputs.iterrows():
    name = row["name"]
    print("Working on %s" % name)
    AAs = format_AAs(row["AAs"])
    positions = format_positions(row["positions"])
    ref_annot = row["reference_annotation_name"]
    new_annot = row["new_annotation_name"]
    ref_filename = row["reference_map"]
    codon_table = pd.read_csv(row["codon_table"], sep="	", index_col=0)
    ref_record = SeqIO.read(ref_filename, "genbank")
    feature_of_interest = None
    new_features = []

    for feat in ref_record.features:
        try:
            if feat.qualifiers["standard_name"] == [ref_annot]:
                feature_of_interest = feat
        except:
            pass
        if not feature_of_interest == feat:
            new_features.append(feat)

    if feature_of_interest is None:
        print("Couldn't find feature '%s'. Quitting..." % ref_annot)
        exit()

    annot_location = feature_of_interest.location
    if not feature_of_interest.location.strand == 1:
        print("Feature of interest must be on the + strand. Quitting...")
        exit()
    sequence = ref_record.seq

    # mutate the codons and make new feature for each one
    new_sequence = sequence.tomutable()
    assert(len(AAs) == len(positions))
    for i in range(len(AAs)):
        codon = get_codon(AAs[i], codon_table)
        codon_pos = positions[i]
        codon_pos_in_annot = positions[i]
        mutate_codon(new_sequence, annot_location.start, codon, codon_pos_in_annot)
        feat_nt_start = annot_location.start + (3 * codon_pos) - 3
        feat = SeqFeature(FeatureLocation(feat_nt_start, feat_nt_start+3), type = "motif")
        feat.qualifiers["standard_name"] =  AAs[i] + str(codon_pos)
        new_features.append(feat)

    # make a new feature for the edited CDS
    new_feature_of_interest = feature_of_interest
    new_feature_of_interest.qualifiers['standard_name'] = new_annot
    new_features.append(new_feature_of_interest)

    new_record = SeqRecord(new_sequence,
                           id = name,
                           name = name,
                           features = new_features,
                           annotations = ref_record.annotations)
    with open(name + '.gb', 'w') as handle:
        SeqIO.write(new_record, handle, 'genbank')





