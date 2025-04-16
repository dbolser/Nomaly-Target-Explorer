import numpy as np
import pytest

TEST_VALUES = [
    (100, 100),  # Small
    (10, 100000),  # Tall (more variants)
    (100000, 10),  # Wide (more individuals)
]

BENCHMARK_VALUES = [
    (10, 1000),
    (100, 1000),
    (1000, 1000),
    (10000, 1000),
    (100000, 1000),
    (1000, 10),
    (1000, 100),
    (1000, 1000),
    (1000, 10000),
    (1000, 100000),
]


def calculate_scores_direct_indexing(genotypes_T, sel_vs):
    """Calculates variant scores using direct NumPy indexing, handling -9."""
    num_individuals, num_variants = genotypes_T.shape

    # Create a mask for missing values
    missing_mask = genotypes_T == -9

    # Create temporary genotypes with -9 replaced by 0 (valid index)
    temp_genotypes = np.where(missing_mask, 0, genotypes_T)

    # Create indices for the variants dimension
    variant_indices = np.arange(num_variants)
    # Use advanced integer indexing with the temporary genotypes
    individual_scores = sel_vs[variant_indices, temp_genotypes]

    # Set scores for missing genotypes to 0
    individual_scores[missing_mask] = 0.0

    return individual_scores


def calculate_scores_take_along(genotypes, sel_vs):
    """Calculates scores using np.take_along_axis without broadcasting."""

    # Create a mask for missing values
    missing_mask = genotypes == -9

    # Create temporary genotypes with -9 replaced by 0 (valid index)
    temp_genotypes = np.where(missing_mask, 0, genotypes)

    # Take along axis=1 (the variant dimension)
    individual_scores = np.take_along_axis(sel_vs, temp_genotypes, axis=1)

    # Set scores for missing genotypes to 0
    individual_scores[missing_mask] = 0.0

    return individual_scores


def calculate_scores_take_along_broadcast(genotypes_T, sel_vs):
    """Calculates scores using np.take_along_axis with broadcasting, handling -9."""
    num_individuals, num_variants = genotypes_T.shape
    if num_variants == 0 or sel_vs.shape[0] == 0:  # Handle edge case
        return np.zeros_like(genotypes_T, dtype=float)

    # Create a mask for missing values
    missing_mask = genotypes_T == -9
    # Replace -9 with 0 for indexing purposes
    temp_genotypes = np.where(missing_mask, 0, genotypes_T)

    # Broadcast sel_vs: (V, 3) -> (I, V, 3)
    # This might be memory-intensive for very large I*V
    b_sel_vs = np.broadcast_to(sel_vs, (num_individuals, num_variants, 3))

    # Reshape temporary genotypes: (I, V) -> (I, V, 1)
    genotype_indices = temp_genotypes[..., np.newaxis]

    # Take along axis=2 (the score dimension)
    # Result shape will be (I, V, 1)
    individual_scores_expanded = np.take_along_axis(b_sel_vs, genotype_indices, axis=2)

    # Squeeze the last dimension: (I, V, 1) -> (I, V)
    individual_scores = np.squeeze(individual_scores_expanded, axis=2)

    # Set scores for missing genotypes to 0
    individual_scores[missing_mask] = 0.0

    return individual_scores


def calculate_scores_masking(genotypes_T, sel_vs):
    """Calculates scores using boolean masking, handling -9."""

    # Create an array to hold scores: (num_individuals, num_variants)
    individual_scores = np.zeros_like(genotypes_T, dtype=float)

    # Create masks for each genotype value (0, 1, 2)
    # -9 values will not match any mask, so their scores remain 0.
    mask_0 = genotypes_T == 0
    mask_1 = genotypes_T == 1
    mask_2 = genotypes_T == 2

    # Apply scores based on masks
    # np.where(mask_0) returns a tuple of arrays (row_indices, col_indices)
    # We only need the column indices (variants) to index sel_vs
    individual_scores[mask_0] = sel_vs[np.where(mask_0)[1], 0]
    individual_scores[mask_1] = sel_vs[np.where(mask_1)[1], 1]
    individual_scores[mask_2] = sel_vs[np.where(mask_2)[1], 2]

    return individual_scores


def validate_methods(genotypes_T, sel_vs, tol=1e-12):
    """Ensures all four methods produce the same result."""

    print("\nCalculating scores (direct indexing)...")
    scores_direct = calculate_scores_direct_indexing(genotypes_T, sel_vs)

    print("Calculating scores (masking)...")
    scores_masking = calculate_scores_masking(genotypes_T, sel_vs)

    print("Calculating scores (take_along_axis)...")
    scores_take_along_broadcast = calculate_scores_take_along_broadcast(
        genotypes_T, sel_vs
    )

    print("Calculating scores (take_along_axis_dan)...")
    # NOTE: I DON'T WANT THIS TO BE PART OF THE BENCHMARK, SO IT'S A LITTLE
    # AWKWARD HERE...
    genotypes = genotypes_T.T
    scores_take_along = calculate_scores_take_along(genotypes, sel_vs).T

    print("Comparing direct vs masking...")
    np.testing.assert_allclose(scores_direct, scores_masking, rtol=tol, atol=tol)

    print("Comparing direct vs take_along...")
    np.testing.assert_allclose(
        scores_direct, scores_take_along_broadcast, rtol=tol, atol=tol
    )

    print("Comparing direct vs take_along_dan...")
    np.testing.assert_allclose(scores_direct, scores_take_along, rtol=tol, atol=tol)

    print("Validation successful: All four methods produce identical results.")


@pytest.fixture(params=TEST_VALUES)
def genotype_data(request):
    num_individuals, num_variants = request.param
    return get_random_genotype_and_score_data(num_individuals, num_variants)


def get_random_genotype_and_score_data(num_individuals, num_variants):
    print(f"\nSetting up data: Individuals={num_individuals}, Variants={num_variants}")

    # Use a fixed seed for reproducibility
    rng = np.random.default_rng(42)

    # Generate random genotypes
    genotypes_T = rng.integers(
        0, 3, size=(num_individuals, num_variants), dtype=np.int8
    )
    # Introduce ~5% missing data (-9)
    missing_indices = rng.choice([True, False], size=genotypes_T.shape, p=[0.05, 0.95])
    genotypes_T[missing_indices] = -9

    # Generate random scores
    sel_vs = rng.random((num_variants, 3), dtype=float)

    return genotypes_T, sel_vs


def test_benchmark_direct_indexing(benchmark, genotype_data):
    genotypes_T, sel_vs = genotype_data
    # group="direct_indexing", disable_gc=True are options for benchmark
    result = benchmark.pedantic(
        calculate_scores_direct_indexing,
        args=(genotypes_T, sel_vs),
        rounds=5,
        iterations=1,
    )
    assert result is not None  # Basic check


def test_benchmark_take_along_broadcast(benchmark, genotype_data):
    genotypes_T, sel_vs = genotype_data
    result = benchmark.pedantic(
        calculate_scores_take_along_broadcast,
        args=(genotypes_T, sel_vs),
        rounds=5,
        iterations=1,
    )
    assert result is not None  # Basic check


def test_benchmark_take_along(benchmark, genotype_data):
    genotypes_T, sel_vs = genotype_data

    # NOTE: I DON'T WANT THIS TO BE PART OF THE BENCHMARK, SO IT'S A LITTLE
    # AWKWARD HERE...
    genotypes = genotypes_T.T

    result = benchmark.pedantic(
        calculate_scores_take_along,
        args=(genotypes, sel_vs),
        rounds=5,
        iterations=1,
    )
    assert result is not None  # Basic check


def test_benchmark_masking(benchmark, genotype_data):
    genotypes_T, sel_vs = genotype_data
    result = benchmark.pedantic(
        calculate_scores_masking, args=(genotypes_T, sel_vs), rounds=5, iterations=1
    )
    assert result is not None  # Basic check


def run_vaidations(extended=False):
    """Run standalone validation."""
    print("Running standalone validation...")
    for num_individuals, num_variants in TEST_VALUES:
        val_genotypes_T, val_sel_vs = get_random_genotype_and_score_data(
            num_individuals, num_variants
        )
        validate_methods(val_genotypes_T, val_sel_vs)
    print("Standalone validation finished.")

    if extended:
        print("Running extended standalone validation...")
        for num_individuals, num_variants in BENCHMARK_VALUES:
            val_genotypes_T, val_sel_vs = get_random_genotype_and_score_data(
                num_individuals, num_variants
            )
            validate_methods(val_genotypes_T, val_sel_vs)
        print("Extended standalone validation finished.")


def main():
    run_vaidations()
    run_benchmarks()


if __name__ == "__main__":
    main()
