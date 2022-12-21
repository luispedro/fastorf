import fasta
import genetic_code
from collections import namedtuple
ORFPos = namedtuple('ORFPos', ['start', 'end', 'rc'])

MIN_LEN = 90
def find_orfs(seq, accept_incomplete=False):
    if accept_incomplete:
        active = [0, 1, 2]
    else:
        active = [-1, -1, -1]
    orfs = []
    for i in range(len(seq)):
        ix = i % 3
        if genetic_code.is_start(seq[i:i+3]):
            if active[ix] == -1:
                active[ix] = i
        if genetic_code.is_stop(seq[i:i+3]):
            if active[ix] != -1:
                if i - active[ix] > MIN_LEN:
                    orfs.append(ORFPos(active[ix], i+3, False))
                active[ix] = -1
    if accept_incomplete:
        for ix in range(3):
            if active[ix] != -1:
                if len(seq) - active[ix] > MIN_LEN:
                    orfs.append(ORFPos(active[ix], len(seq) - (len(seq)-ix)%3, False))
    return orfs


def find_orfs_reverse_complement(seq, accept_incomplete=False):
    orfs = find_orfs(genetic_code.reverse_complement(seq), accept_incomplete)
    return [ORFPos(len(seq)-b, len(seq)-a, True) for a,b,_ in orfs][::-1]


def write_orfs(seq, header, ix, orfs, rc, out, coords_out):
    for i, (a, b, _) in enumerate(orfs):
        out.write(f'>{header}_{i+ix}\n'
            f'{genetic_code.translate(seq[a:b] if not rc else genetic_code.reverse_complement(seq[a:b]))}\n')
        coords_out.write(f'{header}_{i+ix}\t{a}\t{b}\t{b-a}\t{rc}\n')

