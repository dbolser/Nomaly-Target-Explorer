import os
import time

import h5py
import numpy as np
from tqdm import tqdm

from blueprints.nomaly import nomaly_genotype


def benchmark_matrix_access(matrix, n_queries=500):
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
    matrix = nomaly_genotype.genotype_matrix
    print(matrix.shape)
    # results = benchmark_matrix_access(matrix)
    # print(f"Current approach total time: {results['current_approach_time']:.3f}s")
    # print(f"Batched approach total time: {results['batched_approach_time']:.3f}s")
    # print(f"Speedup factor: {results['speedup_factor']:.1f}x")

    assert isinstance(matrix, h5py.Dataset)
    matrix_path = nomaly_genotype.genotype_matrix.file.filename
    matrix_np_path = f"{matrix_path}.npy"

    if not os.path.exists(matrix_np_path):
        np.save(matrix_np_path, matrix)

    memmap_matrix = np.memmap(
        matrix_np_path, mode="r", shape=matrix.shape, dtype=matrix.dtype
    )
    results = benchmark_matrix_access(memmap_matrix, n_queries=5000)

    print(f"Current approach total time: {results['current_approach_time']:.3f}s")
    print(f"Batched approach total time: {results['batched_approach_time']:.3f}s")
    print(f"Speedup factor: {results['speedup_factor']:.1f}x")
