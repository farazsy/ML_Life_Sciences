import numpy as np


def calculate_mad(values):
    """Calculate the Median Absolute Deviation."""
    median = np.median(values)
    mad = np.median(np.abs(values - median))

    return median, mad
