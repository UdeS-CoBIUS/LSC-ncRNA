import unittest
import csv
import tempfile
import os
from parameterized import parameterized

from experiments.scripts.paper_exeperiments import write_results_to_csv  # Make sure this import matches your actual filename

class TestWriteResultsToCSV(unittest.TestCase):
    EXPECTED_HEADER = [
        'Motif Type', 'Dataset Size', 'CM Min', 'CM Max', 'Number motifs',
        'Time cpp', 'File size GB', 'Model', 'Accuracy', 'Precision', 'Recall',
        'F1 Score', 'Training Time (sec)', 'Testing Time (sec)'
    ]

    def setUp(self):
        self.sample_results = {
            'fixed_length': {
                50: {
                    1: {
                        'EXT': {'accuracy': 0.85, 'precision': 0.80, 'recall': 0.75, 'f1_score': 0.77, 'training_time_sec': 120, 'testing_time_sec': 30},
                        'MLP': {'accuracy': 0.88, 'precision': 0.83, 'recall': 0.78, 'f1_score': 0.81, 'training_time_sec': 100, 'testing_time_sec': 25},
                        'num_motifs': 100, 'cpp_execution_time_sec': 10.5, 'file_size_gb': 2.5
                    }
                }
            },
            'combined_length': {
                350: {
                    (2, 5): {
                        'EXT': {'accuracy': 0.89, 'precision': 0.85, 'recall': 0.80, 'f1_score': 0.83, 'training_time_sec': 130, 'testing_time_sec': 40},
                        'MLP': {'accuracy': 0.87, 'precision': 0.82, 'recall': 0.79, 'f1_score': 0.80, 'training_time_sec': 110, 'testing_time_sec': 35},
                        'num_motifs': 90, 'cpp_execution_time_sec': 9.7, 'file_size_gb': 2.3
                    }
                }
            }
        }

    def test_csv_header(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as temp_file:
            temp_file_name = temp_file.name
        
        write_results_to_csv(self.sample_results, temp_file_name)
        
        with open(temp_file_name, 'r') as file:
            csv_reader = csv.reader(file)
            header = next(csv_reader)
        
        self.assertEqual(header, self.EXPECTED_HEADER)
        
        os.remove(temp_file_name)

    @parameterized.expand([
        ('fixed_length', '50', '1', '1', '100', '10.5', '2.5', 'EXT', '0.85', '0.80', '0.75', '0.77', '120', '30'),
        ('fixed_length', '50', '1', '1', '100', '10.5', '2.5', 'MLP', '0.88', '0.83', '0.78', '0.81', '100', '25'),
        ('combined_length', '350', '2', '5', '90', '9.7', '2.3', 'EXT', '0.89', '0.85', '0.80', '0.83', '130', '40'),
        ('combined_length', '350', '2', '5', '90', '9.7', '2.3', 'MLP', '0.87', '0.82', '0.79', '0.80', '110', '35')
    ])
    def test_csv_data_rows(self, *expected_row):
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as temp_file:
            temp_file_name = temp_file.name
        
        write_results_to_csv(self.sample_results, temp_file_name)
        
        with open(temp_file_name, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)  # Skip header
            for row in csv_reader:
                if row[0] == expected_row[0] and row[7] == expected_row[7]:  # Check motif type and model
                    self.assert_rows_equal(row, list(expected_row))
                    break
            else:
                self.fail(f"Expected row not found: {expected_row}")
        
        os.remove(temp_file_name)

    def assert_rows_equal(self, actual_row, expected_row):
        for actual, expected in zip(actual_row, expected_row):
            try:
                actual_float = float(actual)
                expected_float = float(expected)
                self.assertAlmostEqual(actual_float, expected_float, places=2)
            except ValueError:
                self.assertEqual(actual, expected)

    def test_empty_input(self):
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as temp_file:
            temp_file_name = temp_file.name
        
        write_results_to_csv({}, temp_file_name)
        
        with open(temp_file_name, 'r') as file:
            csv_reader = csv.reader(file)
            rows = list(csv_reader)
        
        self.assertEqual(len(rows), 1)  # Only header should be present
        self.assertEqual(rows[0], self.EXPECTED_HEADER)
        
        os.remove(temp_file_name)

if __name__ == '__main__':
    unittest.main()