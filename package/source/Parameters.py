#!/usr/bin/python
# Copyright (c) 2013, Gao Wang <gaow@bcm.edu>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
from argparse import ArgumentParser, SUPPRESS
from SEQLinco.Utils import env, getColumn
from SEQLinco.Core import main

def checkParams(args):
    '''set default arguments or make warnings'''
    args.vcf = os.path.abspath(os.path.expanduser(args.vcf))
    if args.output:
        env.output = os.path.split(args.output)[-1]
    env.missing = args.missing
    #
    if len([x for x in set(getColumn(args.tfam, 6)) if x.lower() not in env.ped_missing]) > 2:
        env.trait = 'quantitative'
    env.log('{} trait detected in [{}]'.format(env.trait.capitalize(), args.tfam))
    if not args.blueprint:
        args.blueprint = os.path.join(env.resource_dir, 'genemap.txt')
    # pop plink/mega2 format to first
    args.format = [x.lower() for x in set(args.format)]
    for item in ['mega2', 'plink']:
        if item in args.format:
            args.format.insert(0, args.format.pop(args.format.index(item)))
    return True
    
class Args:
    def __init__(self):
        self.parser = ArgumentParser(
        description = '''{}, encode sequence data for linkage analysis'''.format(env.proj),
        prog = env.prog,
        fromfile_prefix_chars = '@',
        epilog = '''Copyright (c) 2013 Gao Wang <gaow@bcm.edu> and Di Zhang <di.zhang@bcm.edu> GNU General Public License''')
        self.parser.add_argument('--version',
                                 action='version',
                                 version='%(prog)s version {0}'.format(env.version))
        self.getEncoderArguments(self.parser)
        self.getIOArguments(self.parser)
        self.getRuntimeArguments(self.parser)
        self.parser.set_defaults(func=main)

    def get(self):
        return self.parser.parse_args()

    def getEncoderArguments(self, parser):
        vargs = parser.add_argument_group('Encoder arguments')
        vargs.add_argument('--size', default = 1, type = int,
                           help='''Number of variants to collapse into a bin. Set to -1 to collapse
        all variants into one bin. Default to 1 (1 variant per bin).''')
        vargs.add_argument('--blueprint', metavar = 'FILE',
                           help='''Blueprint file for super marker
        (format: "chr startpos endpos name distance").''')
        vargs.add_argument('--missing', metavar = 'STRING', default="NA",
                           help='''Coding for missing allele. Default to "NA".''')

    def getIOArguments(self, parser):
        vargs = parser.add_argument_group('Input / output options')
        vargs.add_argument('--vcf', metavar='FILE.gz', required=True, help='''Input VCF file, bgzipped.''')
        vargs.add_argument('--tfam', metavar='FILE', required=True,
                           help='''Input pedigree and phenotype information in TFAM format.''')
        vargs.add_argument('-o', '--output', metavar='FILE', help='''Output file name prefix.''')
        vargs.add_argument('--format', metavar = 'FORMAT', nargs='*', default=['plink'],
                           help='''Output format, choose from {} (multiple choices allowed).
        Default to PLINK format.'''.format(repr(env.formats.keys())))
        
    def getRuntimeArguments(self, parser):
        vargs = parser.add_argument_group('Runtime arguments')
        vargs.add_argument('-j', '--jobs', metavar='N', type = int, default = 2,
                           help='''Number of CPUs to use.''')
        vargs.add_argument('--vanilla', action='store_true',
                           help='''Encode data from scratch without loading cache.''')
