#!/usr/bin/env python
# coding: utf-8

import numpy as np
import logging
import igraph as ig
import leidenalg
import os


# package logger
logger = logging.getLogger(__name__)


class ClusterBin:
    def __init__(self, path , contig_name , contig_len , seq_map  , min_binsize , num_gene , random_seed):
        '''
        min_binsize: minimum bin size of output bins
        '''
        self.path = path
        self.name = contig_name
        self.len = contig_len
        self.seq_map = seq_map
        self.binsize = min_binsize
        self.num_gene = num_gene
        self.dist_cluster={}
        self.random_seed = random_seed

        logger.info('Run Leiden Algorithm')
        self.leiden()
        self._write_cluster()
        


    def leiden(self):
        #########Use Leiden Algorithm to do clustering########
        _map_del = self.seq_map.tocoo()
        _vcount = _map_del.shape[0]
        _sources = _map_del.row
        _targets = _map_del.col
        _wei = _map_del.data
        _index = _sources<_targets
        _sources = _sources[_index]
        _targets = _targets[_index]
        _wei = _wei[_index]
        _edgelist = list(zip(_sources, _targets))
        g = ig.Graph(_vcount, _edgelist)
        
        #############determine the resolution parameter###########
        res_option = [1,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500]
        for res in res_option:
            part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=_wei , resolution_parameter = res , n_iterations = -1 , seed = self.random_seed)
            part = list(part)
            numcom = 0
            for ci in range(len(part)):
                if np.sum(self.len[part[ci]]) >= self.binsize:
                    numcom += 1
            if numcom > self.num_gene:
                break
                

        # dict of communities
        numnode = 0
        rang = []
        for ci in range(len(part)):
            if np.sum(self.len[part[ci]]) >= self.binsize:
                rang.append(ci)
                numnode = numnode+len(part[ci])
                for id in part[ci]:
                    self.dist_cluster[self.name[id]] = 'group'+str(ci)

        logger.debug('The optimal resolution is {}'.format(res))
        logger.debug('There are {} contigs in {} bins'.format(numnode , len(rang)))


    def _write_cluster(self):
        ########create file for checkm################
        with open(os.path.join(self.path , 'tmp' , 'cluster.txt'),'w') as out:
            for key , value in self.dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n') 


