#########The structure of the main script is modified from bin3C########
from Script.raw_contact import ContactMatrix, ContactMatrix_LC
from Script.normalized_contact import NormCCMap, NormCCMap_LC
from Script.predict_species_number import gen_bestk
from Script.cluster import ClusterBin
from Script.post_processing import Postprocess, merge_bin
from Script.exceptions import ApplicationException
from Script.utils import load_object, save_object, make_dir, gen_bins, gen_sub_bins, make_random_seed
import scipy.sparse as scisp
import argparse
import warnings
import logging
import shutil
import sys
import os


##Ignore the warning information of package deprecation##
warnings.filterwarnings("ignore")

__version__ = '1.1.0, released at 04/2022'

if __name__ == '__main__':
    
    def mk_version():
        return 'MetaCC v{}'.format(__version__)

    def out_name(base, suffix):
        return '{}{}'.format(base, suffix)

    def ifelse(arg, default):
        if arg is None:
            return default
        else:
            return arg

    runtime_defaults = {
        'min_len': 1000,
        'min_signal': 2,
        'min_mapq': 30,
        'min_match': 30,
        'hic_len':100,
        'thres': 0.05,
        'min_binsize':150000
    }

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-V', '--version', default=False, action='store_true', help='Show the application version')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('--cover', default=False, action='store_true', help='Cover existing files')
    global_parser.add_argument('--log', help='Log file path [OUTDIR/MetaCC.log]')


    parser = argparse.ArgumentParser(description='MetaCC: a scalable and integrative analysis framework for both short-read and long-read metagenomic Hi-C datasets')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')

    cmd_norm = subparsers.add_parser('norm', parents=[global_parser],
                                      description='Normalize contacts.')
                                      
    cmd_cl = subparsers.add_parser('cluster', parents=[global_parser],
                                      description='Do the binning.')
                                      
    cmd_pp = subparsers.add_parser('recluster', parents=[global_parser],
                                      description='post-processing step on partially containminated bins.')

    #cmd_test = subparsers.add_parser('test', parents=[global_parser],
                                        #description='pipeline testing.')

    '''
    Normalization subparser input
    '''
    cmd_norm.add_argument('--min-len', type=int,
                           help='Minimum acceptable reference length [1000]')
    cmd_norm.add_argument('--min-signal', type=int,
                           help='Minimum acceptable signal [2]')
    cmd_norm.add_argument('--min-mapq', type=int,
                           help='Minimum acceptable mapping quality [30]')
    cmd_norm.add_argument('--min-match', type=int,
                           help='Accepted alignments must being N matches [30]')
    cmd_norm.add_argument('-hic-len' , type=int,
                           help='Length of Hi-C reads, used as the scaling factor of CovCC abundances [100]')
    cmd_norm.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                           help='Case-sensitive enzyme name. Use multiple times for multiple enzymes')
    cmd_norm.add_argument('--thres', type=float,
                           help='acceptable fraction of incorrectly identified valid contacts [0.05]')
    cmd_norm.add_argument('FASTA', help='Reference fasta sequence')
    cmd_norm.add_argument('BAM', help='Input bam file in query order')
    cmd_norm.add_argument('OUTDIR', help='Output directory')
         
    
    '''
    Clutering subsparser input
    '''
    cmd_cl.add_argument('--min-binsize', type=int,
                               help='Minimum bin size used in output [150000]')
    cmd_cl.add_argument('--num-gene', type=int,
                               help='Number of maker genes detected, automatically detected if not input')
    cmd_cl.add_argument('-s', '--seed', default=None,
                               help='Random seed')
    cmd_cl.add_argument('FASTA', help='Reference fasta sequence')
    cmd_cl.add_argument('OUTDIR', help='Output directory of sub bins')
    
    
    '''
    Post-processing step
    '''
    cmd_pp.add_argument('--min-binsize', type=int,
                               help='Minimum bin size used in output [150000]')
    cmd_pp.add_argument('FASTA', help='Reference fasta sequence')
    cmd_pp.add_argument('CHECKM',  help='CheckM result')
    cmd_pp.add_argument('OUTDIR', help='Output directory of sub bins')
    
    
    '''
    Testing of NormCC software
    '''
    #cmd_test.add_argument('OUTDIR', help='Output directory of testing results')

    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    try:
        make_dir(args.OUTDIR, args.cover)
    except IOError:
        print('Error: cannot find out directory or the directory already exists')
        sys.exit(1)
        
    temp_folder = os.path.join(args.OUTDIR , 'tmp')
    os.mkdir(temp_folder)
        
    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    # File log listens to all levels from root
    if args.log is not None:
        log_path = args.log
    else:
        log_path = os.path.join(args.OUTDIR, 'MetaCC.log')
    fh = logging.FileHandler(log_path, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:
        if args.command == 'norm':
            if args.enzyme is not None:
                logger.info('Begin constructing raw contact matrix...')
                cm = ContactMatrix(args.BAM,
                                args.enzyme,
                                args.FASTA,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_signal=ifelse(args.min_signal, runtime_defaults['min_signal']),
                                hic_len=ifelse(args.hic_len, runtime_defaults['hic_len']))

                logger.info('Raw contact matrix construction finished')

                logger.info('Begin normalizing raw contacts by NormCC...')
                
                from rpy2 import robjects
                r = robjects.r
                r.source('NormCC/normcc.R')
                contig_file = os.path.join(temp_folder , 'contig_info.csv')
                norm_result = r.normcc(contig_file);logger.info('Coefficient is {}'.format(norm_result))
                
                ######Construct normalized matrix of Hi-C interaction maps#############
                hzmap = NormCCMap(args.OUTDIR,
                                cm.seq_info,
                                cm.seq_map,
                                norm_result,
                                thres = ifelse(args.thres, runtime_defaults['thres']))
                                
                logger.info('NormCC normalization finished')
                            
            else:
                logger.info('Begin constructing raw contact matrix...')
                cm = ContactMatrix_LC(args.BAM,
                                args.FASTA,
                                args.OUTDIR,
                                min_mapq=ifelse(args.min_mapq, runtime_defaults['min_mapq']),
                                min_len=ifelse(args.min_len, runtime_defaults['min_len']),
                                min_match=ifelse(args.min_match, runtime_defaults['min_match']),
                                min_signal=ifelse(args.min_signal, runtime_defaults['min_signal']),
                                hic_len=ifelse(args.hic_len, runtime_defaults['hic_len']))
                                
                logger.info('Raw contact matrix construction finished')

                logger.info('Begin normalizing raw contacts by site-free NormCC due to no enzyme input detected...')
                
                from rpy2 import robjects
                r = robjects.r
                r.source('NormCC/normcc_site_free.R')
                contig_file = os.path.join(temp_folder , 'contig_info.csv')
                norm_result = r.normcc(contig_file);logger.info('Coefficient is {}'.format(norm_result))
                
                ######Construct normalized matrix of Hi-C interaction maps#############
                hzmap = NormCCMap_LC(args.OUTDIR,
                                cm.seq_info,
                                cm.seq_map,
                                norm_result,
                                thres = ifelse(args.thres, runtime_defaults['thres']))
                                
                logger.info('Site-free NormCC normalization finished')
        
            shutil.rmtree(temp_folder) ######Remove all intermediate files#######
            scisp.save_npz(os.path.join(args.OUTDIR, 'Normalized_contact_matrix.npz'), hzmap.seq_map.tocsr())
            save_object(os.path.join(args.OUTDIR, 'NormCC_normalized_contact'), hzmap)
            logger.info('Normalization results have been saved')
            
        if args.command == 'cluster':
            if not os.path.exists(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz')):
                raise IOError('Please run the NormCC normalization step before MagCC binning')
            
            ###########Load the normalization instance to get access to the normalized Hi-C contact maps##########
            logger.info('Loading normalized contact maps by NormCC from: {}'.format(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz')))
            hzmap = load_object(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz'))
            
            #########Scan the marker gene to determine the hyperparameter in the Leiden clustering#########
            if args.num_gene is None:
                logger.info('Begin scanning marker genes...')
                args.num_gene = gen_bestk(args.OUTDIR , args.FASTA)
                if args.num_gene == 0:
                    logger.warning('No marker gene is detected from the assembled contigs!')
                
            logger.info('There are {} marker genes in the assembled contigs'.format(args.num_gene))
            
            if not args.seed:
                args.seed = make_random_seed()
                logger.info('Generated random seed for clustering: {}'.format(args.seed))
            else:
                logger.info('User set random seed for clustering: {}'.format(args.seed))
                
            cluster_process = ClusterBin(args.OUTDIR , hzmap.name , hzmap.len , hzmap.seq_map ,
                                            ifelse(args.min_binsize, runtime_defaults['min_binsize']), args.num_gene, args.seed)
            logger.info('Writing bins...')
            gen_bins(args.FASTA , os.path.join(temp_folder , 'cluster.txt') , os.path.join(args.OUTDIR ,'BIN'))
            shutil.rmtree(temp_folder) ######Remove all intermediate files#######
            logger.info('MetaCC binning fininshed.')


        if args.command == 'recluster':
            if not os.path.exists(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz')):
                logger.error('File NormCC_normalized_contact.gz is not found')
                raise IOError('Please run the NormCC normalization step before postprocessing the contaminated bins')
            if not os.path.exists(os.path.join(args.OUTDIR , 'BIN')):
                logger.error('BIN directory is not found')
                raise IOError('Please run the MagCC binning step before postprocessing the contaminated bins')
            logger.info('Loading normalized contact map instance from: {}'.format(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz')))
            hzmap = load_object(os.path.join(args.OUTDIR , 'NormCC_normalized_contact.gz'))
            post = Postprocess(args.OUTDIR , args.CHECKM , hzmap.name , hzmap.len , hzmap.seq_map, ifelse(args.min_binsize, runtime_defaults['min_binsize']))
            logger.info('Writing sub bins...')
            gen_sub_bins(args.FASTA , os.path.join(args.OUTDIR ,'tmp','cluster_sub.txt') , os.path.join(args.OUTDIR ,'tmp','SUB_BIN'))
            logger.info('Merging bins and sub bins...')
            merge_bin(args.OUTDIR , args.CHECKM)
            shutil.rmtree(temp_folder) ######Remove all intermediate files#######
            logger.info('Post-processing finished.')


    except ApplicationException:
        logger.error('ApplicationException Error')
        sys.exit(1)
