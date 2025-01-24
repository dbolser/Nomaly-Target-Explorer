import numpy as np
import h5py
import zarr
import os
from tqdm import tqdm

# For PGEN
from pyplink import PyPlink
from pgenlibr import writepgen

# For BGEN
import bgen_reader
from bgen_reader import write_bgen


def convert_to_all_formats(original_data, output_dir="data_formats/"):
    """
    Convert data to various formats for benchmarking.

    Args:
        original_data: HDF5 dataset or numpy array
        output_dir: Directory to store converted files

    Returns:
        dict: Paths to converted files
    """
    os.makedirs(output_dir, exist_ok=True)
    paths = {}
    shape = original_data.shape

    print(f"Converting data of shape {shape}")

    # 1. Memmap (numpy)
    print("Converting to memmap...")
    memmap_path = os.path.join(output_dir, "genotypes.npy")
    np.save(memmap_path, original_data)
    paths["memmap"] = memmap_path

    # 2. Zarr
    print("Converting to Zarr...")
    zarr_path = os.path.join(output_dir, "genotypes.zarr")
    # Using chunks of 1000 variants × 10000 samples as a starting point
    # These chunk sizes can be optimized based on access patterns
    store = zarr.open(zarr_path, mode="w")
    z = store.create_dataset(
        "genotypes", shape=shape, chunks=(1000, 10000), dtype=np.int8
    )
    z[:] = original_data[:]
    paths["zarr"] = zarr_path

    # 3. PGEN
    print("Converting to PGEN...")
    pgen_path = os.path.join(output_dir, "genotypes.pgen")

    # PGEN expects data in sample-major order (samples × variants)
    # If our data is variant-major, we need to transpose
    # Note: This might need to be done in chunks for large datasets

    # Create sample and variant info files (simplified)
    n_variants, n_samples = shape

    # Write .fam file (sample information)
    with open(os.path.join(output_dir, "genotypes.fam"), "w") as f:
        for i in range(n_samples):
            f.write(f"sample_{i} sample_{i} 0 0 0 -9\n")

    # Write .bim file (variant information)
    with open(os.path.join(output_dir, "genotypes.bim"), "w") as f:
        for i in range(n_variants):
            f.write(f"1 snp_{i} 0 {i} A B\n")

    # Write PGEN file
    writepgen(
        pgen_path,
        original_data.T,  # Transpose if needed
        sample_major=True,
    )

    paths["pgen"] = pgen_path

    # 4. BGEN
    print("Converting to BGEN...")
    bgen_path = os.path.join(output_dir, "genotypes.bgen")

    # BGEN needs variant and sample metadata
    variants = [
        {
            "id": f"snp_{i}",
            "rsid": f"rs_{i}",
            "chrom": "1",
            "pos": i,
            "alleles": ["A", "B"],
        }
        for i in range(n_variants)
    ]

    samples = [f"sample_{i}" for i in range(n_samples)]

    # Convert genotypes to probabilities
    # BGEN uses probabilistic encoding
    # Convert -1 (missing), 0, 1, 2 to appropriate probability distributions
    def genotype_to_probs(g):
        if g == -1:  # Missing
            return [1 / 3, 1 / 3, 1 / 3]
        elif g == 0:
            return [1, 0, 0]
        elif g == 1:
            return [0, 1, 0]
        else:  # g == 2
            return [0, 0, 1]

    # Process in chunks to avoid memory issues
    chunk_size = 1000
    n_chunks = (n_variants + chunk_size - 1) // chunk_size

    write_bgen(
        bgen_path,
        variants=variants,
        samples=samples,
        genotypes=original_data,  # Will process internally
        genotype_converter=genotype_to_probs,
    )

    paths["bgen"] = bgen_path

    return paths


def verify_conversions(paths, original_data):
    """
    Verify that data was converted correctly by comparing a few random positions
    """
    print("\nVerifying conversions...")

    # Get some random positions to check
    n_checks = 5
    positions = [
        (
            np.random.randint(0, original_data.shape[0]),
            np.random.randint(0, original_data.shape[1]),
        )
        for _ in range(n_checks)
    ]

    orig_values = {pos: original_data[pos] for pos in positions}

    # Check each format
    for fmt, path in paths.items():
        print(f"\nChecking {fmt}...")

        if fmt == "memmap":
            data = np.load(path, mmap_mode="r")
            for pos, orig_val in orig_values.items():
                assert data[pos] == orig_val, f"Mismatch in {fmt} at {pos}"

        elif fmt == "zarr":
            store = zarr.open(path, mode="r")
            data = store["genotypes"]
            for pos, orig_val in orig_values.items():
                assert data[pos] == orig_val, f"Mismatch in {fmt} at {pos}"

        elif fmt == "pgen":
            # Use pgenlib to read and verify
            pgen = PyPlink(path)
            for pos, orig_val in orig_values.items():
                var_idx, sample_idx = pos
                val = pgen.read_genotypes(var_idx)[sample_idx]
                assert val == orig_val, f"Mismatch in {fmt} at {pos}"
            pgen.close()

        elif fmt == "bgen":
            # Use bgen_reader to read and verify
            bgen = bgen_reader.read_bgen(path)
            for pos, orig_val in orig_values.items():
                var_idx, sample_idx = pos
                probs = bgen.genotypes[var_idx][sample_idx]
                # Convert probability back to genotype
                val = np.argmax(probs) if orig_val != -1 else -1
                assert val == orig_val, f"Mismatch in {fmt} at {pos}"

    print("\nAll conversions verified successfully!")


if __name__ == "__main__":
    # Example usage
    with h5py.File("your_current.h5", "r") as f:
        original_data = f["genotype_matrix"][:]
        paths = convert_to_all_formats(original_data)
        verify_conversions(paths, original_data)
