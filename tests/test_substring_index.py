import os
import sys
import time

# Add the src_path to the system path to allow for importing from the parent directory
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
            ground_truth[seqID] = int(start)
    return ground_truth

def test_read_positions(algo, read_sequences, ground_truth):
    """
    Compare the algorithm's predicted positions with the ground truth.
    
    Args:
    - algo (Index): An instance of the Index class to query for read alignments.
    - read_sequences (dict): A dictionary of read IDs and their corresponding sequences.
    - ground_truth (dict): A dictionary with true starting positions for the sequences.
    
    The function prints out statistics for True Positives (TP), False Positives (FP),
    False Negatives (FN), and True Negatives (TN), and calculates Precision and Recall.
    """
    # Initialize counts for True Positive (TP), False Positive (FP), False Negative (FN), and True Negative (TN)
    TP = 0
    FP = 0
    FN = 0
    TN = 0
    total = 0

    # Test the algorithm's query results against the ground truth
    for readID, read in read_sequences.items():
        predicted_start = algo.query(read)
        true_start = ground_truth.get(readID, -1)  # Get the true start from ground truth (-1 if not found)
        total += 1

        if true_start != -1:  # Ground truth says this read has a valid mapping
            if predicted_start != -1:  # The algorithm made a prediction
                # Check if the prediction is within 5 positions of the true start
                if abs(predicted_start - true_start) <= 5:
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

        # Print comparison for each sequence
        print(f"SeqID: {readID}, Match Index: {predicted_start}, Expected Start: {true_start}, Correct: {abs(predicted_start - true_start) < 5}")

    # Calculate Precision and Recall
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0

    # Print the result summary
    print(f"TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    print(f"Correct: {TP + TN}, Total: {total}")
    print(f"Precision: {precision:.2f}")
    print(f"Recall: {recall:.2f}")

if __name__ == "__main__":
    start = time.time()

    # File paths for input data
    fasta_path = "../data/short_reads_ref_genome.fasta"
    fastq_path = "../data/short_reads_2_1000_subset.fastq"
    ground_truth_path = "../data/short_reads_1000_subset_ground_truth.txt"

    # Load the FASTA, FASTQ, and ground truth data
    fasta_reader = FastaReader(fasta_path)
    fastq_reader = FastqReader(fastq_path)
    ground_truth = load_ground_truth(ground_truth_path)

    # Retrieve the reference genome sequence and read sequences
    ref_sequence = fasta_reader.get_genome()
    read_sequences = fastq_reader.get_reads()

    # Create an Index object for querying
    algo = Index(text=ref_sequence, k=30)

    mid = time.time()  # Time checkpoint after loading and preprocessing

    # Run the test by comparing predicted positions with the ground truth
    test_read_positions(algo, read_sequences, ground_truth)

    end = time.time()

    # Print out the timing information
    print("Preprocessing Time:", mid - start)
    print("Test Time:", end - mid)
    print("Total Time:", end - start)
    print(f"Reads Per minute: {60 / (end - start) * 1000}\n")
