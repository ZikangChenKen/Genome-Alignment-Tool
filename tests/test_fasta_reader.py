import os
import unittest
import sys

# Add the src_path to the system path to import FastaReader from the parent directory
src_path = os.path.join(os.getcwd(), '..')
sys.path.append(os.path.abspath(src_path))

from fasta_reader import FastaReader

class TestFastaReader(unittest.TestCase):
    """
    Unit tests for the FastaReader class.
    """

    def setUp(self):
        """
        Set up a mock FASTA file before each test.
        This method creates a small mock FASTA file that will be used in the tests.
        """
        # Mock content for a small FASTA file
        self.fasta_content = """>seq1
ATGCGTACGTAGCTAG
"""
        self.fasta_file = "test.fasta"

        # Write the mock FASTA content to a file
        with open(self.fasta_file, "w") as f:
            f.write(self.fasta_content)

    def test_parse_fasta(self):
        """
        Test the FastaReader's ability to parse a FASTA file and retrieve the genome.
        """
        print("Testing parse_fasta method.")

        # Create an instance of FastaReader
        reader = FastaReader(self.fasta_file)

        # Get the genome sequence parsed by FastaReader
        genome = reader.get_genome()

        # Expected concatenated genome from the mock FASTA file
        expected_genome = "ATGCGTACGTAGCTAG"

        # Assert that the parsed genome matches the expected sequence
        self.assertEqual(expected_genome, genome)
        print(f"Expected genome: {expected_genome}")
        print(f"Actual genome: {genome}\n")

    def test_incorrect_file_path(self):
        """
        Test how FastaReader handles a non-existent FASTA file.
        The test expects a FileNotFoundError to be raised.
        """
        print("Testing incorrect file path handling.")

        # Use a non-existent file path
        incorrect_file_path = "non_existent_file.fasta"

        # Ensure that FastaReader raises FileNotFoundError for a bad file path
        with self.assertRaises(FileNotFoundError):
            reader = FastaReader(incorrect_file_path)
            reader.get_genome()  # This should raise an error
        print("Passed: FileNotFoundError raised as expected.\n")

    def test_empty_fasta_file(self):
        """
        Test how FastaReader handles an empty FASTA file.
        The genome should be an empty string for an empty file.
        """
        print("Testing empty FASTA file handling.")

        # Create an empty FASTA file
        empty_fasta_file = "empty.fasta"
        with open(empty_fasta_file, "w") as f:
            pass  # Just create an empty file

        # Create an instance of FastaReader for the empty file
        reader = FastaReader(empty_fasta_file)

        # The genome sequence should be empty
        genome = reader.get_genome()
        self.assertEqual(genome, "")
        print(f"Expected genome for empty file: ''")
        print(f"Actual genome for empty file: '{genome}'\n")

        # Clean up the empty file
        if os.path.exists(empty_fasta_file):
            os.remove(empty_fasta_file)

    # Uncomment this test if you want to handle invalid DNA sequences.
    # def test_invalid_dna_sequence(self):
    #     """
    #     Test how FastaReader handles a FASTA file with invalid DNA sequence characters.
    #     This test is relevant if the reader is expected to handle such cases.
    #     """
    #     print("Testing invalid DNA sequence handling.")
    #     # Create a FASTA file with invalid DNA sequence characters
    #     invalid_fasta_file = "invalid.fasta"
    #     invalid_fasta_content = """>seq1
    #     ATGCBTDACGTAGCTAG
    #     """
    #
    #     with open(invalid_fasta_file, "w") as f:
    #         f.write(invalid_fasta_content)
    #
    #     # If FastaReader handles invalid characters, test for expected behavior:
    #     reader = FastaReader(invalid_fasta_file)
    #
    #     # Expected behavior: invalid characters (B, D, etc.) are filtered out
    #     expected_genome = "ATGCACGTAGCTAG"
    #     genome = reader.get_genome()
    #
    #     self.assertEqual(expected_genome, genome)
    #     print(f"Expected genome with filtered characters: {expected_genome}")
    #     print(f"Actual genome with filtered characters: {genome}\n")
    #
    #     # Clean up
    #     if os.path.exists(invalid_fasta_file):
    #         os.remove(invalid_fasta_file)

    def tearDown(self):
        """
        Tear down the mock FASTA file after each test.
        This method removes the mock FASTA file created during setUp.
        """
        # Remove the mock FASTA file after each test
        if os.path.exists(self.fasta_file):
            os.remove(self.fasta_file)

if __name__ == "__main__":
    unittest.main()
