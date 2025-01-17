import os, sys
import unittest
from unittest.mock import mock_open, patch, call

# Add the src_path to the system path to import SAMWriter from the parent directory
src_path = os.path.join(os.getcwd(), '..')
sys.path.append(os.path.abspath(src_path))

from sam_writer import SAMWriter  # Adjust this import based on your file structure


class TestSAMWriter(unittest.TestCase):
    """
    Unit tests for the SAMWriter class. The tests use mocking to avoid actual file I/O.
    """

    def setUp(self):
        """
        Set up the test environment by initializing a SAMWriter instance
        with a mock output file name.
        """
        self.output_file = "test.sam"  # Name of the test output SAM file
        self.sam_writer = SAMWriter(self.output_file)

    @patch("builtins.open", new_callable=mock_open)
    def test_write_header(self, mock_file):
        """
        Test the write_header method of SAMWriter. This method writes the SAM file header.
        The test mocks the open function to check if the correct header is written.
        """
        print("Testing write_header method.")

        # Call the write_header method of SAMWriter
        self.sam_writer.write_header()

        # Assert that the file was opened in 'write' mode
        mock_file.assert_called_once_with(self.output_file, 'w')

        # Check if the header was written correctly to the file
        header_written = mock_file().write.call_args[0][0]
        expected_header = "@HD\tVN:1.0\tSO:unsorted\n"
        self.assertEqual(header_written, expected_header)

        # Print expected and actual results for the header
        print(f"Expected header: {expected_header.strip()}")
        print(f"Actual header: {header_written.strip()}\n")

    @patch("builtins.open", new_callable=mock_open)
    def test_write_mapping(self, mock_file):
        """
        Test the write_mapping method of SAMWriter. This method writes a read mapping to the SAM file.
        The test mocks the open function to check if the correct mapping is written.
        """
        print("Testing write_mapping method.")

        # Write the header first to simulate the actual use case
        self.sam_writer.write_header()

        # Define test values for the mapping
        read_name = "read1"
        ref_name = "ref1"
        pos = 100
        seq = "ACGT"

        # Call the write_mapping method with test data
        self.sam_writer.write_mapping(read_name, ref_name, pos, seq)

        # Check if the file was opened in 'append' mode for the mapping
        mock_file.assert_called_with(self.output_file, 'a')

        # Verify the expected calls to the mock file for writing the header and mapping
        expected_calls = [
            call("@HD\tVN:1.0\tSO:unsorted\n"),  # Call for writing the header
            call(f"{read_name}\t0\t{ref_name}\t{pos}\t255\t{len(seq)}M\t*\t0\t0\t{seq}\t*\n")  # Call for the mapping
        ]
        mock_file().write.assert_has_calls(expected_calls)

        # Retrieve the actual written mapping line and verify its correctness
        mapping_written = mock_file().write.call_args_list[-1][0][0]
        expected_mapping = f"{read_name}\t0\t{ref_name}\t{pos}\t255\t{len(seq)}M\t*\t0\t0\t{seq}\t*\n"
        self.assertEqual(mapping_written, expected_mapping)

        # Print expected and actual results for the mapping
        print(f"Expected mapping: {expected_mapping.strip()}")
        print(f"Actual mapping: {mapping_written.strip()}\n")

    def tearDown(self):
        """
        Clean up the test environment by removing the test SAM file if it exists.
        """
        # Remove the test SAM file if it exists
        if os.path.exists(self.output_file):
            os.remove(self.output_file)


if __name__ == "__main__":
    # Run the unit tests
    unittest.main()
