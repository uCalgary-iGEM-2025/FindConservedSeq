from Bio.Seq import Seq, MutableSeq
import parse


def shift_mutations(instances):
    conserved_sequences = parse.conserved_sequences
    # dict {(sequence, pos1, pos2, index)}
    mutations = []

    for seq in conserved_sequences:
        position = conserved_sequences[seq]
        for index in range(len(instances)):
            line = instances[index]
            if line[position[0] - 1: position[1]] != seq and line.find(seq) != -1:
                pos1 = line.find(seq) + 1
                pos2 = pos1 + len(seq) - 1
                mutation_info = seq, pos1, pos2, index
                mutations.append(mutation_info)

    for m in mutations:
        mutated_seq = m[0]
        mutated_pos1 = m[1]
        conserved_pos1, conserved_pos2 = conserved_sequences[mutated_seq]
        if conserved_pos1 > mutated_pos1:
            instances = deletion(instances, [m[0], m[1], m[2], m[3]], conserved_pos1)
        elif conserved_pos1 < mutated_pos1:
            instances, conserved_sequences = insertion(instances, conserved_sequences, [m[0], m[1], m[2], m[3]], conserved_pos1, conserved_pos2)

    top_row = parse.get_updated_top_row(instances)
    return instances, top_row


def deletion(instances, mutation_info, conserved_pos1):
    sequence = mutation_info[0]
    mutation_pos1 = mutation_info[1]
    index = mutation_info[3]
    mutation_seq = MutableSeq(instances[index])

    if (conserved_pos1 - mutation_pos1) <= mutation_seq.count("-"):
        end_seq = mutation_seq[mutation_pos1 - 1:]
        dashes = ""
        for i in range(conserved_pos1 - mutation_pos1):
            dashes = dashes + "-"
        new_seq = end_seq[:len(end_seq) - len(dashes)]
        new_seq = dashes + new_seq
        mutation_seq[mutation_pos1 - 1:] = new_seq

    instances[index] = Seq(mutation_seq)
    return instances


def insertion(instances, conserved_sequences, mutation_info, conserved_pos1, conserved_pos2):
    sequence = mutation_info[0]
    mutation_pos1 = mutation_info[1]
    mutation_pos2 = mutation_info[2]
    indexes = []
    for index in range(len(instances)):
        inst = instances[index]
        if inst[conserved_pos1 - 1: conserved_pos2] == sequence:
            indexes.append(index)

    for num in indexes:
        mutation_seq = MutableSeq(instances[num])
        if (mutation_pos1 - conserved_pos1) <= mutation_seq.count("-"):
            end_seq = mutation_seq[conserved_pos1 - 1:]
            dashes = ""
            for i in range(mutation_pos1 - conserved_pos1):
                dashes = dashes + "-"
            new_seq = end_seq[:len(end_seq) - len(dashes)]
            new_seq = dashes + new_seq
            mutation_seq[conserved_pos1 - 1:] = new_seq
        instances[num] = Seq(mutation_seq)

    conserved_sequences[sequence] = mutation_pos1, mutation_pos2
    return instances, conserved_sequences

