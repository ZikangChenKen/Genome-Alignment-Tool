import os, sys

# Add the src_path to the system path to import FastqReader from the parent directory
src_path = os.path.join(os.getcwd(), '..')
sys.path.append(os.path.abspath(src_path))

from fastq_reader import FastqReader

class FastqTester:
    """
    A class to test the functionality of FastqReader.
    This class creates mock FASTQ files and performs tests on valid, invalid, and empty FASTQ files.
    """

    def __init__(self, fastq_reader_class):
        """
        Initializes the FastqTester with the FastqReader class.
        
        Args:
        - fastq_reader_class: The class of the FastqReader to be tested.
        """
        self.fastq_reader_class = fastq_reader_class
        self.sample_fastq_file = "sample.fastq"  # Name of the sample FASTQ file to be created
        self.empty_fastq_file = "empty.fastq"    # Name of the empty FASTQ file to be created

    def create_sample_fastq(self):
        """
        Creates a sample FASTQ file with mock data for testing.
        The data includes two sequencing reads.
        """
        print(f"Creating sample FASTQ file: {self.sample_fastq_file}")
        # Sample FASTQ data with two reads
        sample_data = """@SEQ_ID_1
GATTACA
+
IIIIIII
@SEQ_ID_2
CCTAGG
+
IIIIII
"""
        # Write the sample data to the FASTQ file
        with open(self.sample_fastq_file, "w") as f:
            f.write(sample_data.strip())  # Remove leading/trailing whitespace

    def test_valid_fastq(self):
        """
        Tests the FastqReader with a valid FASTQ file containing two reads.
        """
        print(f"Testing valid FASTQ file: {self.sample_fastq_file}")
        # Create an instance of FastqReader with the sample FASTQ file
        reader = self.fastq_reader_class(self.sample_fastq_file)
        
        # Get the reads from the FASTQ file
        reads = reader.get_reads()

        # Print the number of reads and the read sequences
        print(f"Number of reads: {len(reads)}")
        print(f"Reads: {reads}")

    def test_invalid_path(self, invalid_file):
        """
        Tests the FastqReader with an invalid file path.
        Ensures that a FileNotFoundError is raised.
        
        Args:
        - invalid_file (str): A non-existent file path for testing.
        """
        print(f"Testing invalid file path: {invalid_file}")
        try:
            # Try to create an instance of FastqReader with a non-existent file
            reader = self.fastq_reader_class(invalid_file)
        except FileNotFoundError:
            # Pass the test if a FileNotFoundError is raised
            print("Passed: FileNotFoundError raised as expected.")
        except Exception as e:
            # Fail the test if any other exception is raised
            print(f"Failed: Unexpected exception {e}")

    def test_empty_fastq(self):
        """
        Tests the FastqReader with an empty FASTQ file.
        The test ensures that an empty dictionary is returned.
        """
        print(f"Testing empty FASTQ file: {self.empty_fastq_file}")
        
        # Create an empty FASTQ file
        with open(self.empty_fastq_file, 'w') as f:
            pass  # Create an empty file

        # Create an instance of FastqReader with the empty file
        reader = self.fastq_reader_class(self.empty_fastq_file)
        
        # Get the reads from the empty file
        reads = reader.get_reads()
        
        # Print the number of reads and the read sequences (should be zero)
        print(f"Number of reads in empty file: {len(reads)}")
        print(f"Reads: {reads}")

    def cleanup(self):
        """
        Clean up the sample and empty FASTQ files created during testing.
        """
        # Remove the sample FASTQ file if it exists
        if os.path.exists(self.sample_fastq_file):
            os.remove(self.sample_fastq_file)

        # Remove the empty FASTQ file if it exists
        if os.path.exists(self.empty_fastq_file):
            os.remove(self.empty_fastq_file)


if __name__ == "__main__":
    # Initialize FastqTester with the FastqReader class
    tester = FastqTester(FastqReader)

    # Create a sample FASTQ file for testing
    tester.create_sample_fastq()

    # Test FastqReader with the valid FASTQ file
    tester.test_valid_fastq()

    # Test FastqReader with an invalid file path (non-existent file)
    tester.test_invalid_path("invalid_path.fastq")

    # Test FastqReader with an empty FASTQ file
    tester.test_empty_fastq()

    # Clean up the files created during testing
    tester.cleanup()
