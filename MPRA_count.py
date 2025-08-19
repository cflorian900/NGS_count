import sys
import os
import pandas as pd
import numpy as np
import time
import dask.dataframe as dd
import concurrent.futures
from string_match import string_match_multi, parse_multi, process_batch
from itertools import repeat
import gc

def main():
    start_time = time.time()
    print(f"Starting processing at {time.strftime('%H:%M:%S')}")
    
    # Parse command line arguments
    fastq_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_prefix = sys.argv[3]
    partition_index = int(sys.argv[4])
    num_partitions = int(sys.argv[5])
    
    # Create counts directory if it doesn't exist
    os.makedirs('./counts', exist_ok=True)
    
    # Use optimized loading of metadata with reduced memory usage
    print("Loading metadata...")
    metadata = dd.read_csv(metadata_file, sep='\t', usecols=['ID', 'seq'])
    ref_seqs_df = metadata.compute()
    
    # Convert to list format for faster processing 
    # (avoiding dataframe overhead during matching)
    ref_seqs = ref_seqs_df.values.tolist()
    del ref_seqs_df  # Free memory
    gc.collect()
    
    print(f"Loaded {len(ref_seqs)} reference sequences in {time.time() - start_time:.2f} seconds")
    
    # Get basename of input file
    base = os.path.basename(fastq_file).split('.')[0]
    print(f"Processing file: {base}")
    
    # Read in data as dask dataframe
    print(f"Reading partition {partition_index} of {num_partitions}...")
    data_loading_start = time.time()
    data = dd.read_csv(fastq_file, sep=' ', header=None)
    data = data.repartition(npartitions=num_partitions)
    
    # Get only the needed partition
    partition_data = data.partitions[partition_index].compute()
    print(f"Time to read partition: {time.time() - data_loading_start:.2f} seconds")
    print(f"Partition size: {len(partition_data)} rows")
    
    # Process data with improved parallelization
    processing_start = time.time()
    
    # Determine optimal chunk size based on system
    chunk_size = min(1000, max(100, len(partition_data) // (os.cpu_count() * 2)))
    
    # Group reads into batches for more efficient processing
    batches = [partition_data[0][i:i+chunk_size].tolist() for i in range(0, len(partition_data), chunk_size)]
    
    print(f"Processing {len(batches)} batches with chunk size {chunk_size}")
    
    # Initialize count dictionary for all reference sequences
    count_dict = {ref_seq[0]: 0 for ref_seq in ref_seqs}
    
    # Process batches in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Map function processes entire batches at once to reduce overhead
        batch_results = executor.map(
            process_batch,
            batches,
            repeat(ref_seqs, len(batches))
        )
        
        # Merge results from all batches
        for batch_counts in batch_results:
            for seq_id, count in batch_counts.items():
                count_dict[seq_id] += count
    
    print(f"Time for parallel processing: {time.time() - processing_start:.2f} seconds")
    
    # Convert results to dataframe and save
    summary_start = time.time()
    summary_df = pd.DataFrame(list(count_dict.items()), columns=['ID', 'count'])
    summary_df.sort_values('ID', inplace=True)
    
    # Save results
    output_path = f'./counts/{output_prefix}.txt'
    summary_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"Time to create and save summary: {time.time() - summary_start:.2f} seconds")
    
    # Print summary statistics
    total_counts = summary_df['count'].sum()
    total_reads = len(partition_data)
    aligned_percentage = total_counts / (total_reads / 4) * 100  # Assuming FASTQ format (4 lines per read)
    
    print(f"Summary for {output_prefix}:")
    print(f"  Total counts: {total_counts}")
    print(f"  Total reads: {total_reads / 4:.0f}")  # Divide by 4 for FASTQ format
    print(f"  Alignment percentage: {aligned_percentage:.2f}%")
    print(f"Total execution time: {time.time() - start_time:.2f} seconds")

if __name__ == '__main__':
    main()
