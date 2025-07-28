from Bio.Seq import Seq, MutableSeq
import parse


def create_linker_strand(target_strand):
    return Seq(target_strand).complement()


def create_incumbent_strand(target_strand):
    incumbent_strand = MutableSeq(target_strand)
    count = 0

    for nuc in range(len(target_strand)):
        count = count + 1
        if count == 4:
            current_nuc = incumbent_strand[nuc]
            if current_nuc == "A":
                incumbent_strand[nuc] = "C"
            elif current_nuc == "C":
                incumbent_strand[nuc] = "A"
            elif current_nuc == "T":
                incumbent_strand[nuc] = "G"
            elif current_nuc == "G":
                incumbent_strand[nuc] = "T"
            count = 0

    return incumbent_strand


def get_longest_conserved():
    conserved_sequences = parse.conserved_sequences
    longest_conserved_sequences = []

    biggest_len = 0
    for seq in conserved_sequences:
        current_len = len(seq)
        if current_len > biggest_len:
            biggest_len = current_len

    for seq in conserved_sequences:
        if len(seq) >= biggest_len:
            longest_conserved_sequences.append([len(seq), seq, conserved_sequences[seq]])

    return longest_conserved_sequences


def calculate_most_conserved_seq(user_length):
    consensus_list = parse.consensus_list
    motif_columns = parse.motif_columns
    sequences_amount = parse.sequences_amount

    # Use the consensus list and motifs to give like ratings for each column
    ratings = []
    for nuc in range(len(consensus_list)):
        if consensus_list[nuc] == "A":
            ratings.append(motif_columns[nuc][0]/sequences_amount)
        elif consensus_list[nuc] == "C":
            ratings.append(motif_columns[nuc][1]/sequences_amount)
        elif consensus_list[nuc] == "G":
            ratings.append(motif_columns[nuc][2]/sequences_amount)
        elif consensus_list[nuc] == "T":
            ratings.append(motif_columns[nuc][3]/sequences_amount)

    # Check the length for each position and check the average
    #   - Get the number of different combinations with the users length and the length of the consensus or rating
    #   - Loop through this length then add each one and divide them by the users length
    biggest_avg = 0
    position = []
    for sec in range(len(ratings) - user_length + 1):
        current_sec = ratings[sec: sec + user_length]
        avg = 0
        for decimal in current_sec:
            avg = avg + decimal
        avg = avg/user_length
        if avg >= biggest_avg:
            biggest_avg = avg
            pos1 = sec + 1
            pos2 = sec + user_length + 1
            position = pos1, pos2

    most_conserved_seq = [biggest_avg*100, position, consensus_list[position[0] - 1: position[1] - 1]]
    return most_conserved_seq

