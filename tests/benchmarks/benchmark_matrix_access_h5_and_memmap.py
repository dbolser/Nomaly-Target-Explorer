import os
import time

import numpy as np
from tqdm import tqdm

from data_services.genotype import GenotypesHDF5

from config import Config


def benchmark_matrix_access(matrix, n_queries=10):
    """
    Benchmark different approaches to matrix row access
    """

    indexes = sorted(np.random.randint(0, matrix.shape[0], size=n_queries))

    # Test 1: Current approach - individual row access
    start = time.time()
    for idx in tqdm(indexes):
        _ = matrix[idx, :]
    current_time = time.time() - start

    # Test 2: Batched approach - get multiple rows at once
    start = time.time()
    _ = matrix[indexes, :]
    batched_time = time.time() - start

    return {
        "current_approach_time": current_time,
        "batched_approach_time": batched_time,
        "speedup_factor": current_time / batched_time,
    }


if __name__ == "__main__":
    geno = GenotypesHDF5(Config.GENOTYPES_HDF)
    matrix_h5 = geno.genotype_matrix_h5
    print(matrix_h5.shape)

    results = benchmark_matrix_access(matrix_h5)
    print(f"Current approach total time: {results['current_approach_time']:.3f}s")
    print(f"Batched approach total time: {results['batched_approach_time']:.3f}s")
    print(f"Speedup factor: {results['speedup_factor']:.1f}x")

    matrix_path = matrix_h5.file.filename
    matrix_np_path = f"{matrix_path}.npy"

    if not os.path.exists(matrix_np_path):
        # This takes some time
        # 1) Convert entire HDF5 dataset to numpy array first
        np_matrix = matrix_h5[:]  # This creates a numpy array in memory
        # 2) Save numpy array to disk
        np.save(matrix_np_path, np_matrix)
        # 3) Free up memory
        del np_matrix

    time_start = time.time()
    matrix_mm = np.load(matrix_np_path, mmap_mode="r")
    time_end = time.time()
    print(f"Time to memmap: {time_end - time_start:.3f}s")

    assert matrix_h5[0, 0] == matrix_mm[0, 0]
    assert np.all(matrix_h5[0:1000, 0] == matrix_mm[0:1000, 0])
    assert np.all(matrix_h5[0:1000, 1000] == matrix_mm[0:1000, 1000])

    results = benchmark_matrix_access(matrix_mm, n_queries=5000)

    print(f"Current approach total time: {results['current_approach_time']:.3f}s")
    print(f"Batched approach total time: {results['batched_approach_time']:.3f}s")
    print(f"Speedup factor: {results['speedup_factor']:.1f}x")
