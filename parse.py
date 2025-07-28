from Bio import SeqIO
from Bio import motifs
import SBOL_exporter

consensus_list = ""
motif_columns = []
conserved_sequences = {}
sequences_amount = 0


def parse_nucleotides(fasta_file):
    global consensus_list
    global motif
    instances = []

    index = 1
    biggest_len = 0
    for seqRecord in SeqIO.parse(fasta_file, "fasta"):
        name = "Comp_" + str(index)
        instances.append(seqRecord.seq.upper())
        current_len = len(seqRecord)
        if current_len > biggest_len:
            biggest_len = current_len
        index = index + 1

    dashes = ""
    for seq in range(len(instances)):
        if len(instances[seq]) < biggest_len:
            num_of_dashes = biggest_len - len(instances[seq])
            for d in range(num_of_dashes):
                dashes = dashes + "-"
            instances[seq] = instances[seq] + dashes
            dashes = ""

    motif = motifs.create(instances)
    consensus_list = motif.consensus

    return instances


def get_updated_top_row(instances):
    global motif
    global consensus_list
    motif = motifs.create(instances)
    consensus_list = motif.consensus
    top_row = construct_top_row(instances)
    return top_row


def construct_top_row(instances):
    global motif_columns
    global conserved_sequences
    global sequences_amount
    motif_columns = []
    conserved_sequences = {}
    top_row = []
    sequences_amount = len(instances)

    # Putting the motifs into a list, so it is easily searchable
    for i in range(motif.length):
        array = [motif.counts[0][i], motif.counts[1][i], motif.counts[2][i], motif.counts[3][i]]
        motif_columns.append(array)

    # Construct a top row with the symbols
    # Give a dictionary for all conserved sequences and their positions
    count = 0
    for nucleotide in range(motif.length):
        if consensus_list[nucleotide] == "A" and 1.0 > (motif_columns[nucleotide][0] / sequences_amount) >= 0.9:
            top_row.append("+")
            count = count + 1
        elif consensus_list[nucleotide] == "A" and (motif_columns[nucleotide][0] / sequences_amount) == 1.0:
            top_row.append("*")
            count = count + 1
        elif consensus_list[nucleotide] == "C" and 1.0 > (motif_columns[nucleotide][1] / sequences_amount) >= 0.9:
            top_row.append("+")
            count = count + 1
        elif consensus_list[nucleotide] == "C" and (motif_columns[nucleotide][1] / sequences_amount) == 1.0:
            top_row.append("*")
            count = count + 1
        elif consensus_list[nucleotide] == "G" and 1.0 > (motif_columns[nucleotide][2] / sequences_amount) >= 0.9:
            top_row.append("+")
            count = count + 1
        elif consensus_list[nucleotide] == "G" and (motif_columns[nucleotide][2] / sequences_amount) == 1.0:
            top_row.append("*")
            count = count + 1
        elif consensus_list[nucleotide] == "T" and 1.0 > (motif_columns[nucleotide][3] / sequences_amount) >= 0.9:
            top_row.append("+")
            count = count + 1
        elif consensus_list[nucleotide] == "T" and (motif_columns[nucleotide][3] / sequences_amount) == 1.0:
            top_row.append("*")
            count = count + 1
        else:
            top_row.append(".")
            if count != 0:
                pos1 = nucleotide - count + 1
                pos2 = nucleotide
                sequence = consensus_list[pos1 - 1: pos2]
                position = pos1, pos2
                conserved_sequences.update({sequence: position})
            count = 0
        if nucleotide == motif.length - 1 and count != 0:
            pos1 = nucleotide + 1 - count + 1
            pos2 = nucleotide + 1
            sequence = consensus_list[pos1 - 1: pos2]
            position = pos1, pos2
            conserved_sequences.update({sequence: position})
    temp_dict = conserved_sequences.copy()
    for seq in temp_dict:
        if len(seq) < 5:
            del conserved_sequences[seq]
    return top_row

