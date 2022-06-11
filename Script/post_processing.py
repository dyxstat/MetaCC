import pandas as pd
import numpy as np
import scipy.sparse as scisp
import tqdm
import Bio.SeqIO as SeqIO
import os
import logging
import igraph as ig
import leidenalg
from utils import open_input, count_fasta_sequences

# package logger
logger = logging.getLogger(__name__)

class Postprocess:

    def __init__(self , path , checkm_result , contig_name , contig_len , seq_map , min_binsize):
        #######Do the recursive Leiden Algorithm##################
        cm = pd.read_csv(checkm_result , sep = ',' , header=None)
        ##ind.csv is the name of fasta file
        cm = cm.values[:,0]

        ##ref help you determine the location of contigs in the fasta file#####
        ref = {}
        numnode = 0
        rang = []
        dist_cluster={}
        self.path = path
        self.name = contig_name
        self.len = contig_len
        self.map = seq_map
        self.binsize = min_binsize

        for counter, value in enumerate(list(self.name)):
            ref[value] = counter

        for k in range(cm.shape[0]):
            fi = str(cm[k]) + '.fa'
            logger.info('Handling fasta file {}'.format(fi))
            name = []
            with open_input(os.path.join(self.path, 'BIN' ,fi)) as multi_fasta:
                fasta_count = count_fasta_sequences(os.path.join(self.path, 'BIN' ,fi))
                for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count):
                    name.append(seqrec.id) 

            index_sub = []
            for i in name:
                index_sub.append(ref[i])

            map_sub = self.map.tocsr()
            map_sub = map_sub[index_sub , :]
            map_sub = map_sub.tocsc()
            map_sub = map_sub[: , index_sub]
            map_sub = map_sub.tocoo()

            len_sub = self.len[index_sub]
            name_sub = self.name[index_sub]
            logger.debug('There are {} contigs with total length {} in {}'.format(len(len_sub),sum(len_sub),fi))

            vcount = map_sub.shape[0]
            sources = map_sub.row
            targets = map_sub.col
            wei = map_sub.data
            index = sources>targets
            sources = sources[index]
            targets = targets[index]
            wei = wei[index]
            edgelist = list(zip(sources, targets))
            g = ig.Graph(vcount, edgelist)

            part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition  , weights=wei , n_iterations=-1)
            part = list(part)

            # dict of communities
            for ci in range(len(part)):
                if np.sum(len_sub[part[ci]]) >= self.binsize:
                    rang.append(ci)
                    numnode = numnode+len(part[ci])
                    for id in part[ci]:
                        dist_cluster[name_sub[id]] = str(cm[k])+str(ci)
        logger.debug('There are {} contigs in {} sub bins'.format(numnode,len(rang)))

        ########create file for checkm################
        with open(os.path.join(self.path ,'tmp','cluster_sub.txt'),'w') as out:
            for key , value in dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n')
                
        
    
def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)

    folder_list_sorte = sorted(folder_list)
    return folder_list_sorte



def merge_bin(pwd_folder , contam_file):
    bin_path = os.path.join(pwd_folder , 'tmp' , 'UNPROCESSED_BIN')
    mv = 'mv' + ' ' + os.path.join(pwd_folder , 'BIN') + ' ' + bin_path
    os.system(mv)
    sub_path = os.path.join(pwd_folder , 'tmp' , 'SUB_BIN')
    output_dir = os.path.join(pwd_folder , 'BIN')
    mk = 'mkdir' + ' ' + output_dir
    os.system(mk)

    bin_list = get_no_hidden_folder_list(bin_path)
    sub_list = get_no_hidden_folder_list(sub_path)

    cm = pd.read_csv(contam_file , sep = ',' , header=None)
    cm = cm.values[:,0]
    cm = list(cm)
    skip = 0
    global_index = 1
    for index1 in range(len(bin_list)):
        if index1 < 10:
            bin_name = 'BIN'+ '000' + str(index1)
        elif index1 >= 10 and index1 < 100:
            bin_name = 'BIN'+ '00' + str(index1)
        elif index1 >= 100 and index1 < 1000:
            bin_name = 'BIN'+ '0' + str(index1)
        else:
            bin_name = 'BIN'+str(index1)

        if bin_name in cm:
            skip += 1
            continue
        else:
            bin = bin_path  + '/' +bin_name + '.fa'
            output = output_dir + '/BIN' + str(global_index) + '.fa'
            mv = 'cp' + ' ' + bin + ' ' + output
            global_index += 1
            os.system(mv)
    for index2 in range(len(sub_list)):
        if index2 < 10:
            bin_name = 'SUB'+ '000' + str(index2)
        elif index2 >= 10 and index2 < 100:
            bin_name = 'SUB'+ '00' + str(index2)
        elif index2 >= 100 and index2 < 1000:
            bin_name = 'SUB'+ '0' + str(index2)
        else:
            bin_name = 'SUB'+str(index2)

        bin = sub_path + '/' + bin_name + '.fa'
        output = output_dir + '/BIN' + str(global_index) + '.fa'
        mv = 'cp' + ' ' + bin + ' ' + output
        global_index += 1
        os.system(mv)
