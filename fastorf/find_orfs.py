from . import genetic_code
from collections import namedtuple
ORFPos = namedtuple('ORFPos', ['start', 'end', 'rc'])

def findall(seq, pats):
    '''Finds all the matches to the set of patterns given
    '''
    matches = []
    for pat in pats:
        next_start = seq.find(pat)
        while next_start >= 0:
            matches.append(next_start)
            next_start = seq.find(pat, next_start+1)
    matches.sort()
    return matches


MIN_LEN = 90
def find_orfs_fwd(seq, accept_incomplete=False):
    '''Find ORFs in the forward strand

    accept_incomplete : bool, optional
    '''
    if accept_incomplete:
        active = [0, 1, 2]
    else:
        active = [-1, -1, -1]
    orfs = []
    positions = findall(seq, genetic_code.START_CODONS + genetic_code.STOP_CODONS)
    for i in positions:
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


def find_orfs_rev(seq, accept_incomplete=False):
    orfs = find_orfs_fwd(genetic_code.reverse_complement(seq), accept_incomplete)
    return [ORFPos(len(seq)-b, len(seq)-a, True) for a,b,_ in orfs][::-1]

def find_orfs(seq, accept_incomplete=False):
    '''Find ORFs in nucletide sequence'''
    return find_orfs_fwd(seq, accept_incomplete) + find_orfs_rev(seq, accept_incomplete)

def extract(seq, orf):
    '''Extract ORF sequence'''
    seq = seq[orf.start: orf.end]
    if orf.rc:
        seq = genetic_code.reverse_complement(seq)
    return seq

def write_orfs(seq, header, orfs, faa_out, coords_out):
    '''
    seq : DNA sequence
    header : Each orf will be named {header}_{ix}
    orfs : list of ORFs
    faa_out : output file object (optional)
    coords_out : output file object (optional)
    '''
    for i,orf in enumerate(orfs):
        if faa_out is not None:
            faa_out.write(f'>{header}_{i} {orf.start}-{orf.end} {"-1" if orf.rc else "+1"}\n')
            faa_out.write(f'{genetic_code.translate(extract(seq, orf))}\n')
        if coords_out is not None:
            coords_out.write(f'{header}_{i}\t{orf.start}\t{orf.end}\t{orf.rc}\n')

