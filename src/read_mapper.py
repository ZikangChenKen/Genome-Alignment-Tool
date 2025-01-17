import os
import re
import time
import argparse
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue
from threading import Thread

# Import necessary modules for substring index, reading FASTA
# and FASTQ files, and writing SAM output.
from substring_index import Index
from fasta_reader import FastaReader
from fastq_reader import FastqReader
from sam_writer import SamWriter

# Define the data directory relative to this script
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def process_chunk(read_chunk, ref_sequence, ref_name, k):
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

    res = []

    # Process each read in the chunk
    for readID, read in read_chunk.items():
        start_cigar_tuple = algo.query(read)
        predicted_start = start_cigar_tuple[0]  # Query the Index for the predicted start
        cigar_string = start_cigar_tuple[1]
        max_score = start_cigar_tuple[2]
        if predicted_start != -1:
            res.append((readID, ref_name, predicted_start, read, cigar_string, max_score))
    if len(res) == 0:
        return None
    return res


def main(fasta_file, fastq_file, output_sam_file=None, threads=48):
    """
    Main function to execute the read mapping process with parallelism.

    Args:
    - fasta_file: Path to the input FASTA file (reference genome)
    - fastq_file: Path to the input FASTQ file (sequenced reads)
    - output_sam_file: Optional path for the output SAM file (default: output.sam in current dir)
    - threads: Number of parallel threads to use for processing
    """
    start = time.time()

    # Initialize FastaReader and FastqReader for reading the input files
    fasta_reader = FastaReader(fasta_file)
    fastq_reader = FastqReader(fastq_file)

    # Extract the reference genome sequence and read sequences
    ref_sequence = fasta_reader.get_genome()
    read_sequences = fastq_reader.get_reads()

    # Set the output SAM file name if not provided
    if output_sam_file is None:
        output_sam_file = os.path.join(os.getcwd(), "output.sam")

    # Initialize the SAM writer to write the results to the SAM file
    sam_writer = SamWriter(output_sam_file)
    sam_writer.write_header()

    # Extract the reference name from the FASTA file name
    ref_name = os.path.basename(fasta_file).split('.')[0]

    print(f"\nStart processing {fasta_file} and {fastq_file}\n")
    mid = time.time()  # Time checkpoint after loading and preprocessing

    # Set the number of workers (processes) to match the number of CPU cores
    num_workers = threads
    print(f"Number of workers: {num_workers}")

    # Convert the read_sequences dict to a list of tuples (readID, read) and split into chunks
    reads = list(read_sequences.items())
    total_reads = len(reads)
    chunk_size = total_reads // num_workers  # Determine chunk size for each worker

    # Create chunks of reads for parallel processing
    chunks = [dict(reads[i:i + chunk_size]) for i in range(0, total_reads, chunk_size)]

    # set seed length to 20
    k = 20

    # Use ProcessPoolExecutor to parallelize the query and comparison process
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks to the pool for each chunk of reads
        futures = [executor.submit(process_chunk, chunk, ref_sequence, ref_name, k) for chunk in chunks]

        # Gather the results from each future (processed chunk) and aggregate the correct predictions
        for future in futures:
            res = future.result()
            if res:
                for result in res:
                    sam_writer.write_mapping(read_name=result[0],
                                             ref_name=result[1],
                                             pos=result[2],
                                             seq=result[3],
                                             cigar=result[4],
                                             max_score=result[5])

    end = time.time()  # End time after testing
    print(f"Result written to {output_sam_file}\n")
    print("Preprocessing Time:", mid - start)
    print("Test Time:", end - mid)
    print("Total Time:", end - start)
    print(f"Reads Per Minute: {len(reads) / (end - mid) * 60}\n")


if __name__ == "__main__":
    # Set up argument parsing for the command-line tool
    parser = argparse.ArgumentParser(description="Read Mapper: Align reads from a FASTQ file \
                                      to a reference genome in a FASTA file.")

    # Required arguments: FASTA and FASTQ file paths, automatically look under ../data
    parser.add_argument('--fasta', required=True, help="Path to the reference genome file in \
                        FASTA format, located under ../data directory.")
    parser.add_argument('--fastq', required=True, help="Path to the sequenced reads file in  \
                        FASTQ format, located under ../data directory.")

    # Optional argument: Output SAM file path
    parser.add_argument('-o', '--output', help="Output SAM file name. If not specified, \
                        the default is 'output.sam' in the current directory.")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Update the paths to look under ../data directory
    fasta_file_path = os.path.join(data_dir, args.fasta)
    fastq_file_path = os.path.join(data_dir, args.fastq)

    # Call the main function with the parsed arguments
    main(fasta_file=fasta_file_path, fastq_file=fastq_file_path, output_sam_file=args.output)
