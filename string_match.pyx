# cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False, cdivision=True
import numpy as np
import pandas as pd
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from cython.operator cimport dereference as deref
from libcpp.pair cimport pair

# Define C++ types for faster operations
ctypedef vector[string] string_vector
ctypedef unordered_map[string, int] string_int_map
ctypedef pair[string, int] string_int_pair

cpdef tuple string_match_multi(str seq_data, list ref_seqs):
    """
    Improved string matching function that returns (ID, count) directly
    
    Args:
        seq_data: The sequence data to search in
        ref_seqs: List of [ID, sequence] pairs to search for
        
    Returns:
        tuple: (hit_count, list of matched IDs)
    """
    cdef list matched_ids = []
    cdef int hit_count = 0
    
    # Convert Python string to C++ string for faster operations
    cdef string cpp_seq_data = seq_data.encode('utf8')
    
    # Search for each reference sequence
    for i in ref_seqs:
        if i[1] in seq_data:  # Keep the Python string check for compatibility
            matched_ids.append(i[0])
            hit_count += 1
            
    return (hit_count, matched_ids)

def build_kmer_index(ref_seqs, k=10):
    """
    Build a k-mer index for faster sequence matching (optional enhancement)
    
    Args:
        ref_seqs: List of [ID, sequence] pairs
        k: k-mer size
        
    Returns:
        dict: Mapping from k-mer to list of sequence IDs containing that k-mer
    """
    kmer_index = {}
    for id_seq in ref_seqs:
        seq_id, seq = id_seq
        seq_len = len(seq)
        
        if seq_len < k:
            # Handle sequences shorter than k
            if seq not in kmer_index:
                kmer_index[seq] = []
            kmer_index[seq].append(seq_id)
        else:
            # Generate k-mers and add to index
            for i in range(seq_len - k + 1):
                kmer = seq[i:i+k]
                if kmer not in kmer_index:
                    kmer_index[kmer] = []
                if seq_id not in kmer_index[kmer]:
                    kmer_index[kmer].append(seq_id)
                    
    return kmer_index

cpdef dict process_batch(list batch_data, list ref_seqs):
    """
    Process a batch of sequences and return counts in a single dictionary
    
    Args:
        batch_data: List of sequence data
        ref_seqs: List of [ID, sequence] pairs
        
    Returns:
        dict: Mapping from sequence ID to count
    """
    cdef dict count_dict = {}
    cdef list matched_ids
    cdef int hit_count
    
    # Initialize count dictionary with zeros
    for ref_seq in ref_seqs:
        count_dict[ref_seq[0]] = 0
    
    # Process each sequence in the batch
    for seq in batch_data:
        hit_count, matched_ids = string_match_multi(seq, ref_seqs)
        
        # Update counts for matched IDs
        for seq_id in matched_ids:
            count_dict[seq_id] += 1
            
    return count_dict

def read_file(filename, chunksize=1000):
    """
    Read a file in chunks to reduce memory usage
    
    Args:
        filename: File to read
        chunksize: Number of lines to read at once
        
    Yields:
        list: Chunk of lines from the file
    """
    with open(filename, 'r') as f:
        chunk = []
        for i, line in enumerate(f):
            chunk.append(line.strip())
            if (i + 1) % chunksize == 0:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

def merge_count_dicts(dict_list):
    """
    Merge multiple count dictionaries
    
    Args:
        dict_list: List of count dictionaries
        
    Returns:
        dict: Merged count dictionary
    """
    if not dict_list:
        return {}
        
    # Start with the first dictionary
    result = dict_list[0].copy()
    
    # Merge the rest
    for d in dict_list[1:]:
        for key, value in d.items():
            if key in result:
                result[key] += value
            else:
                result[key] = value
                
    return result

def parse_multi(results_iter):
    """
    Parse results from parallel processing, directly returning matched IDs
    
    Args:
        results_iter: Iterator of (hit_count, matched_ids) tuples
        
    Returns:
        list: Flattened list of matched IDs
    """
    matched_ids = []
    
    for hit_count, ids in results_iter:
        matched_ids.extend(ids)
        
    return matched_ids
