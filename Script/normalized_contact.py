#!/usr/bin/env python
# coding: utf-8
import numpy as np
from math import log,exp,sqrt
import logging

# package logger
logger = logging.getLogger(__name__)


class NormCCMap:
    def __init__(self, path , contig_info , seq_map , norm_result , thres):
        '''
        perc: threshold of spurious contacts
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.thres = thres
        self.name = []
        self.site = []
        self.len = []
        self.covcc = []

        for i in range(len(contig_info)):
            temp = contig_info[i]
            self.name.append(temp.name)
            self.site.append(temp.sites)
            self.len.append(temp.length)
            self.covcc.append(temp.covcc)

        del contig_info
        
        ####transfer the list to numpy array to do slicing#####
        self.name = np.array(self.name)
        self.site = np.array(self.site)
        self.len = np.array(self.len)
        self.covcc = np.array(self.covcc)
        
        #####Normalize raw contacts######
        self.norm()

        

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _index = _map_row<_map_col
        _map_row = _map_row[_index]
        _map_col = _map_col[_index]
        _map_data = _map_data[_index]
        
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result
        
        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(np.float)
        
        mu_vector = []
        for contig_feature in zip(self.site, self.len, self.covcc):
            mu_vector.append(exp(coeff[0] + coeff[1]*log(contig_feature[0]+1)+ coeff[2]*log(contig_feature[1])+ coeff[3]*log(contig_feature[2]+1)))
        scal = np.max(mu_vector)
        _norm_contact = []
        
        for i in _map_coor:
            x = i[0]
            y = i[1]
            d = i[2]
            
            d_norm = scal*d/sqrt(mu_vector[x]*mu_vector[y])
            _norm_contact.append(d_norm)
            
            self.seq_map[x , y] = d_norm
            self.seq_map[y , x] = d_norm
            
        logger.info('Eliminating systematic biases finished')
        
        ########Remove spurious contacts###########
        cutoffs = np.percentile(_norm_contact , self.thres*100)
        count = 0
        for j in range(len(_norm_contact)):
            x = _map_row[j]
            y = _map_col[j]
            if _norm_contact[j] < cutoffs:
                self.seq_map[x , y] = 0
                self.seq_map[y , x] = 0
                count += 1
        logger.debug('{}% contacts have been removed with the cutoff {}'.format(round(100*count/len(_norm_contact)) , cutoffs))
        logger.info('Spurious contact detection finished')
        
        del _map_row, _map_col, _map_data, _map_coor, _norm_contact, count
        
        
        
        
class NormCCMap_LC:
    def __init__(self, path , contig_info , seq_map , norm_result , thres):
        '''
        perc: threshold of spurious contacts
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.thres = thres
        self.name = []
        self.len = []
        self.covcc = []

        for i in range(len(contig_info)):
            temp = contig_info[i]
            self.name.append(temp.name)
            self.len.append(temp.length)
            self.covcc.append(temp.covcc)

        del contig_info
        
        ####transfer the list to numpy array to do slicing#####
        self.name = np.array(self.name)
        self.len = np.array(self.len)
        self.covcc = np.array(self.covcc)
        
        #####Normalize raw contacts######
        self.norm()

        

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _index = _map_row<_map_col
        _map_row = _map_row[_index]
        _map_col = _map_col[_index]
        _map_data = _map_data[_index]
        
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result
        
        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(np.float)
        
        mu_vector = []
        for contig_feature in zip(self.len, self.covcc):
            mu_vector.append(exp(coeff[0] + coeff[1]*log(contig_feature[0])+ coeff[2]*log(contig_feature[1]+1)))
        scal = np.max(mu_vector)
        _norm_contact = []
        
        for i in _map_coor:
            x = i[0]
            y = i[1]
            d = i[2]
            
            d_norm = scal*d/sqrt(mu_vector[x]*mu_vector[y])
            _norm_contact.append(d_norm)
            
            self.seq_map[x , y] = d_norm
            self.seq_map[y , x] = d_norm
            
        logger.info('Eliminating systematic biases finished')
        
        ########Remove spurious contacts###########
        cutoffs = np.percentile(_norm_contact , self.thres*100)
        count = 0
        for j in range(len(_norm_contact)):
            x = _map_row[j]
            y = _map_col[j]
            if _norm_contact[j] < cutoffs:
                self.seq_map[x , y] = 0
                self.seq_map[y , x] = 0
                count += 1
        logger.debug('{}% contacts have been removed with the cutoff {}'.format(round(100*count/len(_norm_contact)) , cutoffs))
        logger.info('Spurious contact detection finished')
        
        del _map_row, _map_col, _map_data, _map_coor, _norm_contact, count



