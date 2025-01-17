from Bio import SeqIO

class FastaReader:
    """
    A class for reading and parsing a FASTA file to extract the genome sequence.
    """

    def __init__(self, fasta_file):
        """
        Initializes the FastaReader with the provided FASTA file.
        
        Args:
        - fasta_file (str): Path to the FASTA file containing the reference genome.
        """
        self.fasta_file = fasta_file
        self.parse_fastq()  # Parse the FASTA file upon initialization
    
    def parse_fastq(self):
        """
        Parses the FASTA file using the Biopython SeqIO module to extract the full genome sequence.
        The sequence is stored in the `genome` attribute as a string.
        """
        res = ""  # Initialize an empty string to hold the concatenated genome sequence
        with open(self.fasta_file, "r") as handle:
            # Parse the FASTA file and extract the sequences
            for record in SeqIO.parse(handle, "fasta"):
                res += str(record.seq)  # Concatenate all sequences into a single string
        self.genome = res  # Store the genome sequence
    
    def get_genome(self):
        """
        Returns the genome sequence that was parsed from the FASTA file.
        
        Returns:
        - genome (str): The full genome sequence as a single string.
        """
        return self.genome
