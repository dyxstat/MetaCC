import os
import sys
import logging

# package logger
logger = logging.getLogger(__name__)

script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

#######This function is used to scan the number of marker genes from assembled contigs#######
def gen_bestk(dir_path , contig_file):
    ##dir_path path of contig file
    fragURL = os.path.join(f'{script_directory}','Auxiliary', 'FragGeneScan', 'FragGeneScan')
    os.system("chmod 777 " + fragURL)
    fragScanURL = os.path.join(f'{script_directory}','Auxiliary', 'FragGeneScan', 'run_FragGeneScan.pl')
    os.system("chmod 777 " + fragScanURL)
    hmmExeURL = os.path.join(f'{script_directory}','Auxiliary', 'hmmer-3.3.2', 'bin', 'hmmsearch')
    os.system("chmod 777 " + hmmExeURL)
    markerExeURL = os.path.join(f'{script_directory}','Auxiliary' , 'test_getmarker.pl')
    os.system("chmod 777 " + markerExeURL)
    markerURL = os.path.join(f'{script_directory}','Auxiliary' , 'marker.hmm')
    temp_folder = os.path.join(dir_path , 'tmp')
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
