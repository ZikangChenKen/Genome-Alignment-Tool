from Bio import SeqIO

class FastqReader:
    """
    A class for reading and parsing a FASTQ file to extract sequencing reads.
    """

    def __init__(self, fastq_file):
        """
        Initializes the FastqReader with the provided FASTQ file.
        
        Args:
        - fastq_file (str): Path to the FASTQ file containing the sequencing reads.
        """
        self.fastq_file = fastq_file
        self.read_fastq()  # Parse the FASTQ file upon initialization

    def read_fastq(self):
        """
        Parses the FASTQ file using the Biopython SeqIO module and stores the reads in a dictionary.
        The dictionary maps read IDs to their corresponding sequences.
        """
        dic = {}  # Dictionary to store read IDs and their sequences
        with open(self.fastq_file, "r") as handle:
            # Parse the FASTQ file and extract the read ID and sequences
            for record in SeqIO.parse(handle, "fastq"):
                dic[str(record.id)] = str(record.seq)  # Store read ID as key and sequence as value
        self.reads = dic  # Store the dictionary of reads

    def get_reads(self):
        """
        Returns the parsed reads from the FASTQ file.
        
        Returns:
        - reads (dict): A dictionary where the keys are read IDs and the values are the read sequences.
        """
        return self.reads
