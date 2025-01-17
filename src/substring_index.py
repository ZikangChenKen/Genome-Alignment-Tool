import math, time, os, sys

# Add the src_path to the system path so that we can import modules from parent directory
src_path = os.path.join(os.getcwd(), '..')
sys.path.append(os.path.abspath(src_path))

# Import necessary modules for reading FASTA and FASTQ files
from fasta_reader import FastaReader
from fastq_reader import FastqReader


class Index:
    """
    Index class to build a k-mer index of the reference genome and perform 
    sequence alignment using the Smith-Waterman algorithm.
    """
    def __init__(self, text, k):
        """
        Initializes the index for the reference sequence (text) by creating a dictionary
        where k-mers are mapped to their positions in the reference sequence.
        
        Args:
        - text: Reference genome sequence as a string
        - k: Length of the k-mers to index
        """
        self.k = k
        self.text = text
        self.index_map = {}

        # Build the k-mer index for the reference genome
        for i in range(len(text) - k + 1):
            sub = text[i: i + k]
            if sub in self.index_map:
                self.index_map[sub].append(i)  # Append new position for existing k-mer
            else:
                self.index_map[sub] = [i]  # Create a new list for this k-mer

    def compress_cigar(self, cigar_ops):
        """
        Compress a list of CIGAR operations into a valid CIGAR string.
        
        Args:
        - cigar_ops: List of operations (e.g., ['M', 'M', 'I', 'D'])
        
        Returns:
        - A string summarizing the operations (e.g., '2M1I1D')
        """
        if not cigar_ops:
            return ''
        
        compressed = []
        current_op = cigar_ops[0]
        count = 1

        for op in cigar_ops[1:]:
            if op == current_op:
                count += 1
            else:
                compressed.append(f"{count}{current_op}")
                current_op = op
                count = 1

        compressed.append(f"{count}{current_op}")
        return ''.join(compressed)

    # Smith-Waterman algorithm for local sequence alignment
    def smith_waterman(self, ref, read, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
        """
        Smith-Waterman algorithm for local sequence alignment.
        Finds the best local alignment between a reference and a read.
        
        Args:
        - ref: Reference genome sequence
        - read: Read sequence to be aligned
        - match_score: Score for a matching base pair
        - mismatch_penalty: Penalty for mismatching base pairs
        - gap_penalty: Penalty for gaps (insertions/deletions)
        
        Returns:
        - max_score: The highest local alignment score for the read
        """
        n = len(ref)
        m = len(read)
        score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
        traceback_matrix = [[None for _ in range(m + 1)] for _ in range(n + 1)]

        max_score = 0
        max_pos = (-1, -1)

        # Fill the score matrix using dynamic programming
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = score_matrix[i - 1][j - 1] + (match_score if ref[i - 1] == read[j - 1] else mismatch_penalty)
                delete = score_matrix[i - 1][j] + gap_penalty
                insert = score_matrix[i][j - 1] + gap_penalty
                score_matrix[i][j] = max(0, match, delete, insert)

                if score_matrix[i][j] == match:
                    traceback_matrix[i][j] = 'M'  
                elif score_matrix[i][j] == delete:
                    traceback_matrix[i][j] = 'D' 
                elif score_matrix[i][j] == insert:
                    traceback_matrix[i][j] = 'I' 

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        
        cigar = []
        i, j = max_pos

        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            if traceback_matrix[i][j] == 'M':  # Match/Mismatch
                cigar.append('M')
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 'D':  # Deletion
                cigar.append('D')
                i -= 1
            elif traceback_matrix[i][j] == 'I':  # Insertion
                cigar.append('I')
                j -= 1

        cigar.reverse()
        cigar_string = self.compress_cigar(cigar)
        return max_score, cigar_string

    def smith_waterman_banded(self, ref, read, match_score=1, mismatch_penalty=-1, gap_penalty=-1, bandwidth=3):
        """
        Smith-Waterman algorithm for local sequence alignment with banded optimization.
        Finds the best local alignment between a reference and a read within a specified bandwidth.
        
        Args:
        - ref: Reference genome sequence
        - read: Read sequence to be aligned
        - match_score: Score for a matching base pair
        - mismatch_penalty: Penalty for mismatching base pairs
        - gap_penalty: Penalty for gaps (insertions/deletions)
        - bandwidth: Maximum allowed difference from the diagonal

        Returns:
        - max_score: The highest local alignment score for the read
        """
        n = len(ref)
        m = len(read)
        max_score = 0
        
        # Initialize score matrix only within the band
        score_matrix = [[0 for _ in range(2 * bandwidth + 1)] for _ in range(n + 1)]

        # Fill the score matrix using dynamic programming within the band
        for i in range(1, n + 1):
            # Determine the range of columns within the bandwidth
            j_start = max(1, i - bandwidth)
            j_end = min(m, i + bandwidth)
            
            for j in range(j_start, j_end + 1):
                # Calculate the offset for the banded matrix
                band_j = j - (i - bandwidth)
                
                match = score_matrix[i - 1][band_j - 1] + (match_score if ref[i - 1] == read[j - 1] else mismatch_penalty)
                delete = score_matrix[i - 1][band_j] + gap_penalty
                insert = score_matrix[i][band_j - 1] + gap_penalty
                
                score_matrix[i][band_j] = max(0, match, delete, insert)

                if score_matrix[i][band_j] > max_score:
                    max_score = score_matrix[i][band_j]

        return max_score

    def smith_waterman_banded_with_cigar(self, ref, read, match_score=1, mismatch_penalty=-1, gap_penalty=-1, bandwidth=8):
        """
        Smith-Waterman algorithm for local sequence alignment with banded optimization.
        Finds the best local alignment between a reference and a read within a specified bandwidth.
        
        Args:
        - ref: Reference genome sequence
        - read: Read sequence to be aligned
        - match_score: Score for a matching base pair
        - mismatch_penalty: Penalty for mismatching base pairs
        - gap_penalty: Penalty for gaps (insertions/deletions)
        - bandwidth: Maximum allowed difference from the diagonal

        Returns:
        - max_score: The highest local alignment score for the read
        - cigar_string: CIGAR string representing the alignment
        """
        n = len(ref)
        m = len(read)
        max_score = 0
        max_pos = (0, 0)
        
        # Initialize score and traceback matrices
        score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
        traceback_matrix = [[None for _ in range(m + 1)] for _ in range(n + 1)]

        # Fill the score matrix using dynamic programming within the band
        for j in range(1, m + 1):
            for i in range(max(1, j - bandwidth), min(n + 1, j + bandwidth + 1)):
                # update the matrix 
                match_value = match_score if ref[i - 1] == read[j - 1] else mismatch_penalty
                diagonal_score = score_matrix[i - 1][j - 1] + match_value
                vertical_score = score_matrix[i - 1][j] + gap_penalty
                horizontal_score = score_matrix[i][j - 1] + gap_penalty

                score_matrix[i][j] = max(0, diagonal_score, vertical_score, horizontal_score)

                # Traceback
                if score_matrix[i][j] == diagonal_score:
                    traceback_matrix[i][j] = 'M'  # Match/Mismatch
                elif score_matrix[i][j] == vertical_score:
                    traceback_matrix[i][j] = 'D'  # Deletion
                elif score_matrix[i][j] == horizontal_score:
                    traceback_matrix[i][j] = 'I'  # Insertion

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)

        # Perform traceback to generate CIGAR string
        cigar = []
        i, j = max_pos

        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            if traceback_matrix[i][j] == 'M':  # Match/Mismatch
                cigar.append('M')
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 'D':  # Deletion
                cigar.append('D')
                i -= 1
            elif traceback_matrix[i][j] == 'I':  # Insertion
                cigar.append('I')
                j -= 1

        cigar.reverse()
        cigar_string = self.compress_cigar(cigar)
        return max_score, cigar_string

    def reverse_complement(self, sequence):
        """
        Compute the reverse complement of a DNA sequence.
        
        Args:
        - sequence: A string representing a DNA sequence
        
        Returns:
        - Reverse complement of the sequence
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(sequence))

    def query(self, read):
        """
        Query function to align a read against the reference genome using the k-mer index
        and Smith-Waterman local alignment algorithm. Considers both forward and reverse complement reads.
        
        Args:
        - read: A DNA read to align against the reference genome
        
        Returns:
        - The start position of the best alignment in the reference genome, or -1 if no match
        """
        num_positions = 30  # Number of positions to check in the read
        kmer_positions = [i * len(read) // num_positions for i in range(0, num_positions)]
        max_score = -math.inf
        res = -1
        early_exit_threshold = 0.1 * len(read)  # Early exit if score exceeds this threshold

        # Search original read
        for pos in kmer_positions:
            if pos + self.k <= len(read):
                kmer = read[pos:pos + self.k]
                if self.index_map.get(kmer) is not None:
                    for i in self.index_map[kmer]:
                        start = i - pos
                        # return start
                        score, cigar_string = self.smith_waterman_banded_with_cigar(self.text[start:start + len(read)], read)
                        # score = self.smith_waterman_banded(self.text[start:start + len(read)], read)
                        if score > max_score:
                            max_score = score
                            if score >= early_exit_threshold:  # Early exit condition
                                # score, cigar_string = self.smith_waterman_banded_with_cigar(self.text[start:start + len(read)], read)
                                # score, cigar_string = self.smith_waterman(self.text[start:start + len(read)], read)
                                # return start, cigar_string
                                return start, cigar_string, max_score

        # Search reverse complement of the read
        reverse_read = self.reverse_complement(read)
        for pos in kmer_positions:
            if pos + self.k <= len(reverse_read):
                kmer = reverse_read[pos:pos + self.k]
                if self.index_map.get(kmer) is not None:
                    for i in self.index_map[kmer]:
                        start = i - pos
                        # return start
                        score, cigar_string = self.smith_waterman_banded_with_cigar(self.text[start:start + len(read)], reverse_read)
                        # score = self.smith_waterman_banded(self.text[start:start + len(reverse_read)], reverse_read)
                        if score > max_score:
                            max_score = score
                            if score >= early_exit_threshold:  # Early exit condition
                                # score, cigar_string = self.smith_waterman_banded_with_cigar(self.text[start:start + len(read)], reverse_read)
                                # score, cigar_string = self.smith_waterman(self.text[start:start + len(read)], reverse_read)
                                # return start, cigar_string
                                return start, cigar_string, max_score

        return res, '', 255

    def process_read(readID, read, algo):
        """
        Process a single read by aligning it against the reference genome using the Index algorithm.
        
        Args:
        - readID: Identifier for the read
        - read: Sequence of the read
        - algo: The Index object for querying
        
        Returns:
        - readID: The identifier of the read
        - position: The best alignment position in the reference genome
        """
        position = algo.query(read)
        return readID, position


if __name__ == "__main__":
    # Timing the entire process
    start = time.time()

    # Define paths for the input FASTA (reference genome) and FASTQ (reads) files
    fasta_path = "../data/short_reads_ref_genome.fasta"
    fastq_path = "../data/short_reads_2_1000_subset.fastq"

    # Initialize FastaReader and FastqReader to read the input files
    fasta_reader = FastaReader(fasta_path)
    fastq_reader = FastqReader(fastq_path)

    # Get the reference genome sequence and read sequences
    ref_sequence = fasta_reader.get_genome()
    read_sequences = fastq_reader.get_reads()

    # Initialize the Index with k-mer size of 10
    algo = Index(text=ref_sequence, k=10)
    mid = time.time()  # Time taken for preprocessing

    cnt = 0  # Counter for aligned reads

    # Iterate through each read and find its alignment position in the reference genome
    for readID, read in read_sequences.items():
        position_cigar_tuple = algo.query(read)
        position = position_cigar_tuple[0]
        cigar = position_cigar_tuple[1]
        max_score = position_cigar_tuple[2]
        if position != -1:
            cnt += 1  # Count aligned reads
        print(f"Read {readID} aligned at position: {position} with max score: {max_score} and cigar string: {cigar}")

    end = time.time()  # Time taken for aligning reads

    # Print results
    print(f"Accuracy upper bound: {cnt / 1000:.2%}.")  # Assuming 1000 reads
    print("Preprocessing Time:", mid - start)
    print("Read Time:", end - mid)
    print("Total Time:", end - start)
    print("Reads Per Minute:", 60 / (end - mid) * 1000)  # Compute reads per minute

    print('================================')
    algo_test = Index(text='AAAGTCTAGAA', k=1)
    score_cigar_tuple_test = algo_test.smith_waterman_banded_with_cigar('AAAGTCTAGAA', 'GTCGATAG')
    score_test = score_cigar_tuple_test[0]
    cigar_test = score_cigar_tuple_test[1]
    print(f"Read {readID} aligned at position: {score_test} with cigar string: {cigar_test}")

    print('========================')
    score_cigar_tuple_test = algo_test.smith_waterman('AAAGTCTAGAA', 'GTCGATAG')
    score_test = score_cigar_tuple_test[0]
    cigar_test = score_cigar_tuple_test[1]
    print(f"Read {readID} aligned at position: {score_test} with cigar string: {cigar_test}")


