import numpy as np
import pandas as pd
import logging
import os
from math import log, exp, sqrt
from scipy import sparse as scisp

logger = logging.getLogger(__name__)

class NormCCMap:
    def __init__(self, path, contig_info_df, seq_map, norm_result, thres=0.05):
        """
        Initializes the NormCCMap object.

        Parameters:
        - path (str): Path to the metacc folder.
        - contig_info_df (pd.DataFrame): DataFrame containing contig information.
        - seq_map (scipy.sparse matrix): Raw contact matrix.
        - norm_result (list): List of normalization coefficients.
        - thres (float): Threshold for spurious contact detection.
        """
        self.contig_info = contig_info_df
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.thres = thres

        # Log the columns present
        logging.debug(f"Columns in contig_info_metacc.csv: {self.contig_info.columns.tolist()}")

        try:
            self.name = self.contig_info['name'].tolist()
            self.site = self.contig_info['sites'].tolist()
            self.len = self.contig_info['length'].tolist()
            self.covcc = self.contig_info['covcc'].tolist()
        except KeyError as e:
            logging.error(f"Expected column {e} not found in contig_info_metacc.csv")
            raise e

        # Convert lists to numpy arrays for efficient processing
        self.name = np.array(self.name)
        self.site = np.array(self.site)
        self.len = np.array(self.len)
        self.covcc = np.array(self.covcc)

        # Perform normalization
        self.norm()

    def norm(self):
        coeff = self.norm_result
        if coeff is None:
            logger.error("Normalization coefficients are None.")
            raise ValueError("Normalization coefficients are None.")

        logger.debug(f"Normalization coefficients: {coeff}")

        if len(coeff) == 4:
            # Using 'sites', 'length', 'covcc'
            mu_vector = [
                exp(
                    coeff[0] + 
                    coeff[1] * log(site + 1e-9) + 
                    coeff[2] * log(length + 1e-9) + 
                    coeff[3] * log(covcc + 1e-9)
                ) 
                for site, length, covcc in zip(self.site, self.len, self.covcc)
            ]
            logger.debug(f"mu_vector (with 'sites'): min={np.min(mu_vector)}, max={np.max(mu_vector)}, mean={np.mean(mu_vector)}")
        elif len(coeff) == 3:
            # 'sites' is constant; using 'length', 'covcc'
            mu_vector = [
                exp(
                    coeff[0] + 
                    coeff[1] * log(length + 1e-9) + 
                    coeff[2] * log(covcc + 1e-9)
                ) 
                for length, covcc in zip(self.len, self.covcc)
            ]
            logger.debug(f"mu_vector (without 'sites'): min={np.min(mu_vector)}, max={np.max(mu_vector)}, mean={np.mean(mu_vector)}")
        else:
            logger.error(f"Unexpected number of coefficients: {len(coeff)}")
            raise ValueError("Normalization coefficients count mismatch.")

        # Ensure seq_map is in COO format for iteration
        seq_map_coo = self.seq_map.tocoo()
        _map_row = seq_map_coo.row
        _map_col = seq_map_coo.col
        _map_data = seq_map_coo.data

        # Ensure seq_map is in LIL format for efficient modifications
        self.seq_map = self.seq_map.tolil().astype(float)

        # Calculate scaling factor
        scal = np.max(mu_vector)
        logger.debug(f"Scaling factor (scal): {scal}")

        _norm_contact = []

        # Iterate over the non-zero elements of the matrix
        for i, j, d in zip(_map_row, _map_col, _map_data):
            if i < j:  # Ensure each pair is processed once
                try:
                    d_norm = scal * d / sqrt(mu_vector[i] * mu_vector[j])
                except IndexError:
                    logger.error(f"IndexError for i={i}, j={j}. Check mu_vector length.")
                    raise
                except ZeroDivisionError:
                    logger.error(f"ZeroDivisionError for mu_vector[i]={mu_vector[i]}, mu_vector[j]={mu_vector[j]}.")
                    d_norm = 0
                _norm_contact.append(d_norm)
                self.seq_map[i, j] = d_norm
                self.seq_map[j, i] = d_norm

        logger.info('Eliminating systematic biases finished')

        # Remove spurious contacts based on threshold
        cutoffs = np.percentile(_norm_contact, self.thres * 100)
        logger.debug(f"Cutoff for spurious contacts (threshold={self.thres}): {cutoffs}")

        excluded_contacts = 0
        for idx, val in enumerate(_norm_contact):
            if val < cutoffs:
                i, j = _map_row[idx], _map_col[idx]
                self.seq_map[i, j] = 0
                self.seq_map[j, i] = 0
                excluded_contacts += 1

        logger.debug(f"Excluded {excluded_contacts} contacts below cutoff {cutoffs}")
        logger.info('Spurious contact detection finished')

        # Convert back to CSR format for efficient arithmetic operations
        self.seq_map = self.seq_map.tocsr()

        # Log final statistics
        total_contacts = len(_norm_contact)
        remaining_contacts = total_contacts - excluded_contacts
        logger.debug(f"Total contacts: {total_contacts}, Remaining contacts: {remaining_contacts}")

        if remaining_contacts == 0:
            logger.warning("All contacts have been excluded after normalization and thresholding.")

        
class NormCCMap_LC:
    def __init__(self, path, contig_info_df, seq_map, norm_result, thres=0.05):
        """
        Initializes the NormCCMap_LC object.

        Parameters:
        - path (str): Path to the metacc folder.
        - contig_info_df (pd.DataFrame): DataFrame containing contig information.
        - seq_map (scipy.sparse matrix): Raw contact matrix.
        - norm_result (list): List of normalization coefficients.
        - thres (float): Threshold for spurious contact detection.
        """
        self.contig_info = contig_info_df
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.thres = thres

        # Log the columns present
        logging.debug(f"Columns in contig_info_metacc.csv: {self.contig_info.columns.tolist()}")

        try:
            self.name = self.contig_info['name'].tolist()
            self.len = self.contig_info['length'].tolist()
            self.covcc = self.contig_info['covcc'].tolist()
        except KeyError as e:
            logging.error(f"Expected column {e} not found in contig_info_metacc.csv")
            raise e

        # Convert lists to numpy arrays for efficient processing
        self.name = np.array(self.name)
        self.len = np.array(self.len)
        self.covcc = np.array(self.covcc)

        # Perform normalization
        self.norm()

    def norm(self):
        coeff = self.norm_result
        if coeff is None:
            logger.error("Normalization coefficients are None.")
            raise ValueError("Normalization coefficients are None.")

        logger.debug(f"Normalization coefficients: {coeff}")


        if len(coeff) == 4:
            # Using 'sites', 'length', 'covcc'
            mu_vector = [
                exp(
                    coeff[0] + 
                    coeff[1] * log(site + 1e-9) + 
                    coeff[2] * log(length + 1e-9) + 
                    coeff[3] * log(covcc + 1e-9)
                ) 
                for site, length, covcc in zip(self.site, self.len, self.covcc)
            ]
            logger.debug(f"mu_vector (with 'sites'): min={np.min(mu_vector)}, max={np.max(mu_vector)}, mean={np.mean(mu_vector)}")
        elif len(coeff) == 3:
            # 'sites' is constant; using 'length', 'covcc'
            mu_vector = [
                exp(
                    coeff[0] + 
                    coeff[1] * log(length + 1e-9) + 
                    coeff[2] * log(covcc + 1e-9)
                ) 
                for length, covcc in zip(self.len, self.covcc)
            ]
            logger.debug(f"mu_vector (without 'sites'): min={np.min(mu_vector)}, max={np.max(mu_vector)}, mean={np.mean(mu_vector)}")
        else:
            logger.error(f"Unexpected number of coefficients: {len(coeff)}")
            raise ValueError("Normalization coefficients count mismatch.")

        # Ensure seq_map is in COO format for iteration
        seq_map_coo = self.seq_map.tocoo()
        _map_row = seq_map_coo.row
        _map_col = seq_map_coo.col
        _map_data = seq_map_coo.data

        # Ensure seq_map is in LIL format for efficient modifications
        self.seq_map = self.seq_map.tolil().astype(float)

        # Calculate scaling factor
        scal = np.max(mu_vector)
        logger.debug(f"Scaling factor (scal): {scal}")

        _norm_contact = []

        # Iterate over the non-zero elements of the matrix
        for i, j, d in zip(_map_row, _map_col, _map_data):
            if i < j:  # Ensure each pair is processed once
                try:
                    d_norm = scal * d / sqrt(mu_vector[i] * mu_vector[j])
                except IndexError:
                    logger.error(f"IndexError for i={i}, j={j}. Check mu_vector length.")
                    raise
                except ZeroDivisionError:
                    logger.error(f"ZeroDivisionError for mu_vector[i]={mu_vector[i]}, mu_vector[j]={mu_vector[j]}.")
                    d_norm = 0
                _norm_contact.append(d_norm)
                self.seq_map[i, j] = d_norm
                self.seq_map[j, i] = d_norm

        logger.info('Eliminating systematic biases finished')

        # Remove spurious contacts based on threshold
        cutoffs = np.percentile(_norm_contact, self.thres * 100)
        logger.debug(f"Cutoff for spurious contacts (threshold={self.thres}): {cutoffs}")

        excluded_contacts = 0
        for idx, val in enumerate(_norm_contact):
            if val < cutoffs:
                i, j = _map_row[idx], _map_col[idx]
                self.seq_map[i, j] = 0
                self.seq_map[j, i] = 0
                excluded_contacts += 1

        logger.debug(f"Excluded {excluded_contacts} contacts below cutoff {cutoffs}")
        logger.info('Spurious contact detection finished')

        # Convert back to CSR format for efficient arithmetic operations
        self.seq_map = self.seq_map.tocsr()

        # Log final statistics
        total_contacts = len(_norm_contact)
        remaining_contacts = total_contacts - excluded_contacts
        logger.debug(f"Total contacts: {total_contacts}, Remaining contacts: {remaining_contacts}")

        if remaining_contacts == 0:
            logger.warning("All contacts have been excluded after normalization and thresholding.")