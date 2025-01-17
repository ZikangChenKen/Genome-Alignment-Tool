import threading
class SamWriter:
    def __init__(self, output_file):
        """
        Initializes the SamWriter with the output SAM file.
        
        Args:
        - output_file (str): Path to the output SAM file where mappings will be written.
        """
        self.output_file = output_file
        self.file = open(output_file, 'w')  # Open the file once for the lifetime of the object
        self.lock = threading.Lock()

    def write_header(self):
        """
        Writes the SAM file header to the output file.
        The header contains metadata about the SAM file, such as version and sorting order.
        In this case, the header indicates version 1.0 and that the file is unsorted.
        """
        # with open(self.output_file, 'w') as file:
        #     file.write("@HD\tVN:1.0\tSO:unsorted\n")  # Example SAM header

        with self.lock:
            self.file.write("@HD\tVN:1.0\tSO:unsorted\n")

    def write_mapping(self, read_name, ref_name, pos, seq, cigar, max_score):
        """
        Writes a single read mapping to the SAM file, including a dynamic CIGAR string.

        Args:
        - read_name (str): The name or identifier of the read.
        - ref_name (str): The name of the reference sequence (usually the reference genome).
        - pos (int): The 1-based position where the read aligns to the reference.
        - seq (str): The sequence of the read.
        - cigar (str): The CIGAR string representing the alignment.

        The SAM format fields in each line include:
        - read_name: The identifier for the read.
        - 0: A flag indicating that the read is mapped with a positive strand (simplified here).
        - ref_name: The name of the reference where the read maps.
        - pos: The alignment position in the reference.
        - 255: The mapping quality (255 means unavailable).
        - cigar: The CIGAR string (e.g., "76M" for a full match of length 76).
        - *: Placeholder for the mate reference name (not used here).
        - 0: Placeholder for the mate's alignment position.
        - 0: Placeholder for the template length.
        - seq: The sequence of the read.
        - *: Placeholder for the read's quality scores (not used here).
        """
        # with open(self.output_file, 'a') as file:
        #     # file.write(f"{read_name}\t0\t{ref_name}\t{pos}\t255\t{cigar}\t*\t0\t0\t{seq}\t*\n")
        #     file.write(f"{read_name}\t0\t{ref_name}\t{pos}\t{max_score}\t{cigar}\t{seq}\n")
        with self.lock:
            self.file.write(
                f"{read_name}\t0\t{ref_name}\t{pos}\t255\t{cigar}\t*\t0\t0\t{seq}\t*\tAS:i:{max_score}\n"
            )


    def close(self):
        """
        Closes the SAM file.
        """
        with self.lock:
            self.file.close()
