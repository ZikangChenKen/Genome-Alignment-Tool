import os
import sys
import time
import re
import psutil
import tracemalloc
from concurrent.futures import ProcessPoolExecutor

# Add the src_path to the system path to import modules from the ../src directory
src_path = os.path.join(os.getcwd(), '../src')
sys.path.append(os.path.abspath(src_path))

# Import necessary modules for substring index, reading FASTA and FASTQ files
from substring_index import Index
from fasta_reader import FastaReader
from fastq_reader import FastqReader

def load_ground_truth(file_path):
    """
    Load the ground truth data from a file.

    Args:
    - file_path (str): Path to the ground truth file containing true mapping positions.
    
    Returns:
    - ground_truth (dict): A dictionary where keys are sequence IDs and values are the true starting positions.
    """
    ground_truth = {}
    with open(file_path, 'r') as f:
        for line in f:
            seqID, start, end = line.strip().split()
            ground_truth[seqID] = [int(start), int(end)]
    return ground_truth

def parse_cigar(cigar):
    """
    Parse the CIGAR string to calculate the length of the alignment on the reference genome.

    Args:
    - cigar (str): The CIGAR string from the alignment.

    Returns:
    - alignment_length (int): The total span of the alignment on the reference genome.
    """
    # Match all CIGAR operations with their respective counts (e.g., "10M2D3I")
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    alignment_length = 0

    # Sum up contributions to alignment length on the reference genome
    for length, op in operations:
        if op in {'M', 'D'}:  
            alignment_length += int(length)
    
    return alignment_length

def process_chunk(read_chunk, ref_sequence, ground_truth, k):
    """
    Process a chunk of reads by aligning them to the reference sequence and comparing
    the results to the ground truth.

    Args:
    - read_chunk (dict): A chunk of read sequences (ID -> sequence).
    - ref_sequence (str): The reference genome sequence.
    - ground_truth (dict): The ground truth positions for each read.
    - k (int): Length of k-mers used in the Index.

    Returns:
    - correct_count (int): Number of correctly predicted read positions.
    """
    # Create an Index object for querying the reference genome
    algo = Index(text=ref_sequence, k=k)

    # Initialize metrics
    correct_count = 0
    TP, FP, FN, TN = 0, 0, 0, 0
    
    # Initialize count of correct predictions
    correct_count = 0
    
    # Process each read in the chunk
    for readID, read in read_chunk.items():
        start_cigar_tuple = algo.query(read)
        predicted_start = start_cigar_tuple[0]  # Query the Index for the predicted start
        cigar_string = start_cigar_tuple[1]
        alignment_length = parse_cigar(cigar_string)
        predicted_end = predicted_start + alignment_length
        if predicted_start == -1: # When there is no alignment
            predicted_end = -1
        true_tuple = ground_truth.get(readID, [-1, -1])  # Get the true start from ground truth (-1 if not found)
        true_start = true_tuple[0]
        true_end = true_tuple[1]
        # print(f'Predicted Start: {predicted_start}, Predicted End: {predicted_end}, True Start: {true_start}, True End: {true_end}')
        # Count as correct if the prediction is within 5 bases of the true start
        if (abs(predicted_start - true_start) <= 5) and (abs(predicted_end - true_end) <= 5):
            correct_count += 1
        # if (abs(predicted_start - true_start) <= 5) and (abs(predicted_end - true_end) <= 5):
        #     correct_count += 1
        #     tp += 1  # True positive: Correct alignment
        # else:
        #     if predicted_start != -1:
        #         fp += 1  # False positive: Incorrect alignment
        #     if true_start != -1:
        #         fn += 1  # False negative: Missed correct alignment
        if true_start != -1:  # Ground truth says this read has a valid mapping
            if predicted_start != -1:  # The algorithm made a prediction
                # Check if the prediction is within 5 positions of the true start
                if abs(predicted_start - true_start) <= 5 and abs(predicted_end - true_end) <= 5:
                    TP += 1  # Correct prediction (True Positive)
                else:
                    FP += 1  # Incorrect prediction (False Positive)
            else:
                FN += 1  # No prediction but there was a valid ground truth (False Negative)
        else:  # Ground truth says this read has no valid mapping
            if predicted_start != -1:
                FP += 1  # Incorrect prediction (False Positive)
            else:
                TN += 1  # Correctly unmapped (True Negative)



        # if (abs(predicted_start - true_start) <= 5):
        #     correct_count += 1
    return correct_count, TP, FP, FN, TN

if __name__ == "__main__":

    start = time.time()
    # Track initial memory usage
    process = psutil.Process()
    initial_memory = process.memory_info().rss  # Memory in bytes

    # File paths for input data
    # fasta_path = "../data/short_reads_reference_genome.fasta"
    # fastq_path = "../data/challenging_dataset_2.fastq"
    # ground_truth_path = "../data/challenging_ground_truth.txt"

    # Midterm Dataset
    fasta_path = "../data/final_data/short_reads_reference_genome.fasta"
    fastq_path = "../data/final_data/challenging_dataset_1.fastq"
    ground_truth_path = "../data/final_data/challenging_ground_truth.txt"

    # Final Dataset
    # fasta_path = "../data/final_data/short_reads_reference_genome.fasta"
    # fastq_path = "../data/final_data/challenging_dataset_1.fastq"
    # # fastq_path = "../data/challenging_dataset_2.fastq"
    # ground_truth_path = "../data/final_data/challenging_ground_truth.txt"

    # Load the FASTA, FASTQ, and ground truth data
    fasta_reader = FastaReader(fasta_path)
    fastq_reader = FastqReader(fastq_path)
    ground_truth = load_ground_truth(ground_truth_path)

    # Retrieve the reference genome sequence and read sequences
    ref_sequence = fasta_reader.get_genome()
    read_sequences = fastq_reader.get_reads()

    k = 20  # Length of k-mers for substring indexing

    mid = time.time()  # Time checkpoint after loading and preprocessing

    # Set the number of workers (processes) to match the number of CPU cores
    num_workers = 48
    print(f"Number of workers: {num_workers}")

    # Convert the read_sequences dict to a list of tuples (readID, read) and split into chunks
    reads = list(read_sequences.items())
    total_reads = len(reads)
    chunk_size = total_reads // num_workers  # Determine chunk size for each worker

    # Create chunks of reads for parallel processing
    chunks = [dict(reads[i:i + chunk_size]) for i in range(0, total_reads, chunk_size)]

    correct_predictions = 0  # Track the total number of correct predictions
    total_tp, total_fp, total_fn, total_tn = 0, 0, 0, 0

    # Use ProcessPoolExecutor to parallelize the query and comparison process
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks to the pool for each chunk of reads
        futures = [executor.submit(process_chunk, chunk, ref_sequence, ground_truth, k) for chunk in chunks]
        
        # Gather the results from each future (processed chunk) and aggregate the correct predictions
        for future in futures:
            correct_count, tp, fp, fn, tn = future.result()
            correct_predictions += correct_count
            total_tp += tp
            total_fp += fp
            total_fn += fn
            total_tn += tn

    end = time.time()  # End time after testing
    # Get final memory usage
    final_memory = process.memory_info().rss
    memory_used_mb = (final_memory - initial_memory) / (1024 ** 2)  # Convert to MB

    # Calculate the accuracy as a percentage
    accuracy = (correct_predictions / total_reads) * 100
    precision = (total_tp / (total_tp + total_fp)) if (total_tp + total_fp) > 0 else 0
    recall = (total_tp / (total_tp + total_fn)) if (total_tp + total_fn) > 0 else 0
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    print(f'TP: {total_tp}')
    print(f'FP: {total_fp}')
    print(f'FN: {total_fn}')
    print(f'TN: {total_tn}')

    print(f"Correctness: {accuracy:.2f}%")
    print(f"Precision: {precision:.2f}")
    print(f"Recall: {recall:.2f}")
    print(f"F1 Score: {f1_score:.2f}")

    # Print timing information
    print("Preprocessing Time:", mid - start)
    print("Test Time:", end - mid)
    print("Total Time:", end - start)
    print(f"Memory Used: {memory_used_mb:.2f} MB")
    print(f"Reads Per minute: {60 / (end - start) * total_reads}")