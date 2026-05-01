import polars as pl
import scipy.stats as stats
import numpy as np

if __name__ == "__main__":

    # Read data in chunk using Polars;
    # Perform statistical tests;
    # and calcutae chunk-level mean depth.
    batch_size = 1000
    for chunk in pl.read_csv_batched("am_extreme_depth.csv", batch_size=batch_size):
        depths = chunk["depth"].to_numpy()
        # Perform the Kolmogorov-Smirnov test against a uniform distribution
        res = stats.chisquare(depths)
        if res.p_value < 0.05:
            raise ValueError("Reject the null hypothesis: Depths are not uniformly distributed.")
        # Calculate and print the mean depth for the chunk
        mean_depth = np.mean(depths)
        print(f"Mean Depth: {mean_depth}")
