import os
import sys
import logging

# package logger
logger = logging.getLogger(__name__)

metacc_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#######This function is used to scan the number of marker genes from assembled contigs#######
def gen_bestk(dir_path , contig_file):
    ##dir_path path of contig file
    fragScanURL = 'run_FragGeneScan.pl'
    hmmExeURL = 'hmmsearch'
    markerExeURL = os.path.join(metacc_directory, 'Auxiliary' , 'test_getmarker.pl')
    markerURL = os.path.join(metacc_directory, 'Auxiliary' , 'marker.hmm')
    
    temp_folder = os.path.join(dir_path , 'tmp')
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)

    seedURL = os.path.join(temp_folder , 'contigs.seed')
    fragResultURL = os.path.join(temp_folder , 'contigs.frag.faa')
    hmmResultURL = os.path.join(temp_folder , 'contigs.hmmout')

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + os.path.join(temp_folder , 'contigs.frag') + " -complete=0 -train=complete -thread=10 1>" + os.path.join(temp_folder , 'contigs.frag.out') + " 2>" + os.path.join(temp_folder , 'contigs.frag.err')
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu 48 " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1000 " + seedURL
                os.system(markerCmd)

            if os.path.exists(seedURL):
                candK = file_len(seedURL)
                return(candK)
            else:
                return(0)
        else:
            logger.error("Hmmsearch failed! Not exist: " + hmmResultURL)
            sys.exit()
    else:
        logger.error("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()



def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
