from libencoder import * #

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

def main(args):
    '''the main encoder function'''
    status = checkParams(args)
    status = downloadResources([('http://tigerwang.org/uploads/genemap.txt', env.resource_dir),
                               ('http://tigerwang.org/uploads/tabix', env.resource_bin),
                               ('http://tigerwang.org/uploads/plink', env.resource_bin)])
    cache = Cache(env.cache_dir, env.output)
    # STEP 1: write encoded data to TPED format
    if not args.vanilla and cache.check():
        env.log('Loading data from archive ...')
        cache.load()
    else:
        status = indexVCF(args.vcf)
        samples_vcf = extractSamplenames(args.vcf)
        status = checkSamples(samples_vcf, getColumn(args.tfam, 2))
        families = extractFamilies(args.tfam)
        if not families:
            env.error('No family found in [{}]!'.format(args.tfam), exit = True)
        else:
            # samples have to be in tfam file and be in the order of their appearance in vcf file
            samples = [x for x in samples_vcf if x in flatten([flatten(item) for item in families.values()])]
            env.log('{:,d} families with a total of {:,d} samples will be processed'.\
                    format(len(families), len(samples)))
        status = rewriteFamfile(args.tfam, samples)
        #
        with open(args.blueprint, 'r') as f:
            regions = [x.strip().split() for x in f.readlines()]
        env.log('Scanning [{}] for {:,d} pre-defined units'.format(args.vcf, len(regions)))
        env.jobs = max(min(args.jobs, len(regions)), 1)
        regions.extend([None] * env.jobs)
        queue = Queue()
        try:
            for i in regions:
                queue.put(i)
            jobs = [EncoderWorker(
                i,
                queue,
                RData(samples, families),
                RegionExtractor(args.vcf, [idx + 9 for idx, x in enumerate(samples_vcf) if x not in samples]),
                MendelianErrorChecker(),
                GenoEncoder(args.size),
                LinkageWritter()
                ) for i in range(env.jobs)]
            for j in jobs:
                j.start()
            for j in jobs:
                j.join()
        except KeyboardInterrupt:
            # FIXME: need to properly close all jobs
            raise ValueError("Use 'killall {}' to properly terminate all processes!".format(env.prog))
        else:
            env.log('\nArchiving to directory [{}]'.format(env.cache_dir))
            cache.write(pre = env.output, ext = '.tped', otherfiles = [args.vcf, args.tfam, args.blueprint])
    # STEP 2: write to PLINK or mega2 format
    tpeds = [os.path.join(env.cache_dir, item) for item in os.listdir(env.cache_dir) if item.startswith(env.output) and item.endswith('.tped')]
    if 'plink' in args.format:
        env.log('Saving data to directory [PLINK] ...')
        formatPlink(tpeds, [args.tfam] * len(tpeds), 'PLINK')
    if 'mega2' in args.format:
        env.log('Saving data to directory [MEGA2] ...')
        formatMega2('MEGA2/{}*'.format(env.output))
    if 'mlink' in args.format:
        env.log('Saving data to directory [MLINK] ...')
        formatMlink(tpeds, [args.tfam] * len(tpeds), 'MLINK')
    runMlink()
    plotMlink()
    #please implement this section
    # STEP final: clean up unwanted files
    #for item in args.format:
    #    removeFiles(item.upper(), exclude = env.formats[item])
    #removeFiles(env.cache_dir, exclude = ['.cache'])
    
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

if __name__ == '__main__':
    # check platform
    import platform
    system = platform.system()
    arch = platform.architecture()[0]
    if not (system == 'Linux' and arch == '64bit'):
        env.error("{} cannot be executed on {}-{} computers "\
                  "(only 64bit Linux environment is supported)".format(env.prog, system, arch), exit = True)
    # encode
    try:
        args = Args().get()
        args.func(args)
    except Exception as e:
        raise
        env.error("{}".format(e))
