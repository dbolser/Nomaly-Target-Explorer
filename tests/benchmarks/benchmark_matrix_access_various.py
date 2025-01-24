import numpy as np
import h5py
import zarr
import time
import os
import psutil
from tqdm import tqdm
import resource


def get_memory_usage():
    """Get current memory usage in MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def benchmark_format(data, n_queries=500, batch_size=100):
    """
    Benchmark different formats with memory tracking

    Args:
        data: Original data array/memmap/etc
        n_queries: Number of random queries to perform
        batch_size: Size of batches for batch queries
    """
    shape = data.shape
    results = {}

    # Memory before operations
    initial_mem = get_memory_usage()

    # Test 1: Single row random access
    start = time.time()
    for _ in tqdm(range(n_queries), desc="Single queries"):
        idx = np.random.randint(0, shape[0])
        _ = data[idx : idx + 1, :]
    single_time = time.time() - start

    # Memory after single queries
    single_mem = get_memory_usage() - initial_mem

    # Test 2: Batched random access
    start = time.time()
    n_batches = n_queries // batch_size
    for _ in tqdm(range(n_batches), desc="Batch queries"):
        indices = np.random.randint(0, shape[0], size=batch_size)
        _ = data[indices, :]
    batch_time = time.time() - start

    # Memory after batch queries
    batch_mem = get_memory_usage() - initial_mem

    # Test 3: Sequential access
    start = time.time()
    for i in tqdm(range(0, min(shape[0], n_queries)), desc="Sequential"):
        _ = data[i, :]
    seq_time = time.time() - start

    # Final memory usage
    final_mem = get_memory_usage() - initial_mem

    return {
        "single_query_time": single_time,
        "batch_query_time": batch_time,
        "sequential_time": seq_time,
        "single_query_memory": single_mem,
        "batch_query_memory": batch_mem,
        "final_memory": final_mem,
    }


def run_all_benchmarks(shape, original_data):
    """Run benchmarks for all formats"""
    results = {}

    # 1. Memmap
    print("\nTesting memmap...")
    memmap_file = "test_memmap.npy"
    np.save(memmap_file, original_data)
    memmap_data = np.load(memmap_file, mmap_mode="r")
    results["memmap"] = benchmark_format(memmap_data)

    # 2. Zarr
    print("\nTesting Zarr...")
    zarr_store = zarr.zeros(shape, chunks=(1000, 1000), dtype=np.int8)
    zarr_store[:] = original_data[:]
    results["zarr"] = benchmark_format(zarr_store)

    # 3. HDF5 (current)
    print("\nTesting HDF5...")
    with h5py.File("test.h5", "w") as f:
        dset = f.create_dataset("genotypes", data=original_data)
        results["hdf5"] = benchmark_format(dset)

    return results


def print_results(results):
    """Print formatted results"""
    print("\nBenchmark Results:")
    print("-" * 80)
    formats = list(results.keys())
    metrics = [
        "single_query_time",
        "batch_query_time",
        "sequential_time",
        "single_query_memory",
        "batch_query_memory",
        "final_memory",
    ]

    # Print header
    print(f"{'Format':<10} " + " ".join(f"{m[:15]:>15}" for m in metrics))
    print("-" * 80)

    # Print results
    for fmt in formats:
        line = f"{fmt:<10} "
        for metric in metrics:
            if "time" in metric:
                line += f"{results[fmt][metric]:>15.3f}s "
            else:
                line += f"{results[fmt][metric]:>15.1f}MB "
        print(line)


if __name__ == "__main__":
    # Your matrix dimensions
    SHAPE = (83011, 488377)

    # Run benchmarks
    with h5py.File("your_current.h5", "r") as f:
        original_data = f["genotype_matrix"]
        results = run_all_benchmarks(SHAPE, original_data)

    print_results(results)
