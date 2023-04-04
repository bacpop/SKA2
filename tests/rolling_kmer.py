def main():
    sequence = "NNAGACAANCCATNA"
    k = 3
    rolling_kmers(sequence, k)

def rolling_kmers(sequence, k):
    pos = 0
    current_kmer_dict = {}
    current_kmer_dict = make_the_mers(pos, sequence, k, current_kmer_dict)
    print(current_kmer_dict)


def init(pos, k, sequence):
    current_kmer = sequence[pos:pos + k]
    return current_kmer


def make_the_mers(pos, sequence, k, current_kmer_dict):
    # pos is the position of the beginning of the current k-mer
    # current k-mer
    current_kmer = init(pos, k, sequence)

    if 'N' in current_kmer:
        pos_of_N = current_kmer.index('N')
        pos = pos + pos_of_N + 1  # start with new k-mer directly after the N
        if pos < len(sequence) - k + 1:
            make_the_mers(pos, sequence, k, current_kmer_dict)
    else:
        split_kmer = current_kmer[0] + current_kmer[2]
        split_base = current_kmer[1]
        # add split k-mer and split base to dictionary
        current_kmer_dict[split_kmer] = split_base

        # add to position to get next k-mer
        pos += 1

        # getting next k-mer from new pos
        for i in range(pos, len(sequence) - k + 1):
            new_base = sequence[i + k - 1]
            if new_base == 'N':  # the new base is an N
                # update position
                pos = i + k + 1
                if pos < len(sequence) - k + 1:
                    return make_the_mers(pos, sequence, k, current_kmer_dict)
                else:
                    return current_kmer_dict
            else:
                # TODO: needs to be refined for longer k-mers

                split_base1 = split_kmer[1]  # last base for k=3
                split_kmer = split_base + new_base

                # add split k-mer and split base to dictionary
                current_kmer_dict[split_kmer] = split_base1
                split_base = split_base1

    return current_kmer_dict


if __name__ == "__main__":
    main()
