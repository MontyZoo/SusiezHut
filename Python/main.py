from Bio.Seq import Seq
from collections import OrderedDict


_fasta_file = '/users/karl/Downloads/dna1.fasta'
#_fasta_file = '/users/karl/Downloads/dna.example.fasta'
_seq_dict = OrderedDict()

def read_seqs(file):
    name = ''
    for line in iter(file):
        line = line.rstrip()
        if line.startswith('>'):
            words = line.split()
            name = words[0][1:]
            _seq_dict[name] = {}
            _seq_dict[name]["seq"] = ''
        else:
            _seq_dict[name]["seq"] = _seq_dict[name]["seq"] + line


def get_frame_orfs(seq):
    orfs = []
    orf_start = -1
    i = 0
    while i < len(seq):
        codon = seq[i:i+3]
        if orf_start < 0:
            if codon == "ATG":
                orf_start = i
        else:
            if codon == "TAA" or codon == "TAG" or codon == "TGA":
                orfs.append((seq[orf_start:i+3], orf_start+1))
                orf_start = -1
        i += 3
    return orfs

def find_orfs_forward(frames, frame_filter = None):
    if not frame_filter:
        frame_filter = (1,2,3)
    orfs = []
    for frame_id in frame_filter:
        orfs += get_frame_orfs(frames[frame_id])
    return orfs

def find_orfs(seq, frame_filter = None, check_reverse = True):
    frames = OrderedDict()
    frames[1] = seq
    frames[2] = seq[1:]
    frames[3] = seq[2:]
    orfs = find_orfs_forward(frames, frame_filter)

    if (check_reverse):
        reverse_seq = Seq(seq).reverse_complement()
        r_frames = OrderedDict()
        r_frames[1] = reverse_seq
        r_frames[2] = reverse_seq[1:]
        r_frames[3] = reverse_seq[2:]
        orfs += find_orfs_forward(reverse_seq)

    return orfs


def find_repeat_in_seq(seq, dict, len):
    seq_len = seq.__len__()
    for index, val in enumerate(seq):
        if index + len <= seq_len:
            repeat = seq[index:index+len]
            if dict.has_key(repeat):
                dict[repeat] += 1
            else:
                dict[repeat] = 1


def find_repeats(len):
    repeat_dict = {}
    for title, entry in _seq_dict.iteritems():
        find_repeat_in_seq(entry["seq"], repeat_dict, len)

    max_repeats = []
    max_repeat_count = 0
    for repeat, count in repeat_dict.iteritems():
        if max_repeat_count < count:
            max_repeat_count = count
            max_repeats = [repeat]
        elif max_repeat_count == count:
            max_repeats.append(repeat)

    print "Max Repeat Length of length %d: %d" % (len, max_repeat_count)
    print "Number of max repeat strings: %d" % max_repeats.__len__()


def run():
    file = None
    try:
        file = open(_fasta_file)
        read_seqs(file)
    except Exception, e:
        print e
    finally:
        if file:
            file.close()

    print "Counts: " + str(_seq_dict.__len__())

    longest = -1
    shortest = 999999999
    orf_count = 0
    longest_orf = -1
    orf_start = -1
    for title, entry in _seq_dict.iteritems():
        seq = entry["seq"]
        length = len(seq)
        if length > longest:
            longest = length
        if shortest > length:
            shortest = length

        if title != 'gi|142022655|gb|EQ086233.1|97':
            continue
        orfs = find_orfs(seq, [1,2,3])
        for orf in orfs:
            orf_count += 1
            orf_len = len(orf[0])
            if orf_len > longest_orf:
                longest_orf = orf_len
                orf_start = orf[1]

    print "Longest: " + str(longest)
    print "Shortest: " + str(shortest)
    print "ORF Count: " + str(orf_count)
    print "Longest ORF: " + str(longest_orf)
    print "Longest ORF Start Pos: " + str(orf_start)

    find_repeats(7)

run()
