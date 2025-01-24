import numpy as np
import time


def benchmark_matrix_access(matrix_shape=(10000, 100000), n_queries=50):
    """
    Benchmark different approaches to matrix row access
    """
    # Create a test matrix of similar size
    matrix = np.random.randint(0, 3, size=matrix_shape)

    # Test 1: Current approach - individual row access
    start = time.time()
    for _ in range(n_queries):
        idx = np.random.randint(0, matrix_shape[0])
        _ = matrix[idx, :]
    current_time = time.time() - start

    # Test 2: Batched approach - get multiple rows at once
    start = time.time()
    indices = np.random.randint(0, matrix_shape[0], size=n_queries)
    _ = matrix[indices, :]
    batched_time = time.time() - start

    return {
        "current_approach_time": current_time,
        "batched_approach_time": batched_time,
        "speedup_factor": current_time / batched_time,
    }


if __name__ == "__main__":
    results = benchmark_matrix_access()
    print(f"Current approach total time: {results['current_approach_time']:.3f}s")
    print(f"Batched approach total time: {results['batched_approach_time']:.3f}s")
    print(f"Speedup factor: {results['speedup_factor']:.1f}x")
