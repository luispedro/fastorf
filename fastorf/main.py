import argparse

from . import find_orfs
from . import fastorf_version
from .fasta import fasta_iter

def main(args=None):
    if args is None:
        import sys
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(description='Find ORFs in a FASTA file using a fast and naïve algorithm.')
    parser.add_argument('fasta', help='FASTA file')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + fastorf_version.__version__)
    parser.add_argument('-a', '--faa', dest='faa_out', help='Output file for amino acid sequences')
    args = parser.parse_args(args)

    with open(args.faa_out, 'wt') as faa:
        for h, seq in fasta_iter(args.fasta):
            orfs = find_orfs.find_orfs(seq)
            find_orfs.write_orfs(seq, h, orfs, faa, None)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
