import sys
import os
import unittest
import pandas as pd
import numpy as np

# Dynamically add the bin directory to the python path so we can import app.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin')))

# Provide a dummy path in sys.argv so app.py doesn't crash on import
if len(sys.argv) < 2:
    sys.argv.append("dummy_data_dir")

import app

class TestAppHelpers(unittest.TestCase):

    def test_most_common(self):
        # Test typical case
        series = pd.Series(["apple", "banana", "apple", "orange"])
        self.assertEqual(app.most_common(series), "apple")
        
        # Test empty series
        empty_series = pd.Series([], dtype=str)
        self.assertEqual(app.most_common(empty_series, default="Unknown"), "Unknown")
        
        # Test series with only NaNs
        nan_series = pd.Series([np.nan, np.nan])
        self.assertEqual(app.most_common(nan_series, default="NA"), "NA")

    def test_clean_taxon_name(self):
        # Test basic string
        self.assertEqual(app.clean_taxon_name("Panthera leo"), "Panthera leo")
        
        # Test semicolon separation
        self.assertEqual(app.clean_taxon_name("Panthera leo;Panthera tigris"), "Panthera leo")
        
        # Test genus rank (keeps only the first word)
        self.assertEqual(app.clean_taxon_name("Panthera leo", rank="genus"), "Panthera")
        
        # Test species rank (keeps exactly two words)
        self.assertEqual(app.clean_taxon_name("Panthera leo subspecies", rank="species"), "Panthera leo")
        self.assertEqual(app.clean_taxon_name("Panthera", rank="species"), "Panthera")
        
        # Test NaN/missing data
        self.assertEqual(app.clean_taxon_name(np.nan), "Unknown")
        self.assertEqual(app.clean_taxon_name(None), "Unknown")

    def test_shorten(self):
        # Test text shorter than max_len
        self.assertEqual(app.shorten("short string", max_len=20), "short string")
        
        # Test text exactly equal to max_len
        self.assertEqual(app.shorten("1234567890", max_len=10), "1234567890")
        
        # Test text longer than max_len
        self.assertEqual(app.shorten("this is a very long string", max_len=10), "this is a ...")
        
        # Test non-string input
        self.assertEqual(app.shorten(12345, max_len=3), "123...")

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'])
