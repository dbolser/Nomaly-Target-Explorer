"""Library for causal analysis of genes to disease."""

import logging
import re
import sys
from pathlib import Path
from typing import Any, Tuple

import bnlearn as bn
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats
from flask import Blueprint, Response, current_app, request, session
from sklearn.impute import SimpleImputer

from config import Config
from data_services import GenotypeService, NomalyScoreService, PhenotypeService
from db import get_term_variants


def read_files(
    term: str,
    phecode: str,
    phenotype_service: PhenotypeService,
    nomaly_score_service: NomalyScoreService,
    genotype_service: GenotypeService,
    ancestry: str = "EUR",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Get data from the 'Nomaly Data Services'."""
    phenotypes = phenotype_service.get_cases_for_phecode(phecode, ancestry=ancestry)
    # phenotypes.groupby(["sex", "phenotype"]).count()

    # Affected only?
    # phenotypes = phenotypes[phenotypes["phenotype"] == 1]
    # phenotypes.groupby(["sex", "phenotype"]).count()

    assert isinstance(phenotypes["eid"].values, np.ndarray)

    scores = nomaly_score_service._hdf.get_scores_by_eids_unsorted(
        phenotypes["eid"].values, np.array([term])
    )
    assert len(scores) == len(phenotypes)

    # merge phenotypes and scores
    phenotypes["score"] = scores

    # read variants for term
    term_variants = get_term_variants(term)

    # read genotypes for variants
    variants_genotypes = genotype_service.get_genotypes(
        eids=phenotypes.eid, vids=term_variants.variant_id, nomaly_ids=True
    )

    # Go consistency!
    assert variants_genotypes.shape == (term_variants.shape[0], phenotypes.shape[0])

    # And make everything messy again so that we can return to business as usual...

    scores = phenotypes[["eid", "sex", "score", "phenotype"]].set_index("eid")
    scores.rename(columns={"phenotype": "pheno"}, inplace=True)
    variants_info = term_variants[["variant_id", "gene", "aa", "hmm_score"]].set_index(
        "variant_id"
    )
    variants_genotypes = pd.DataFrame(
        data=variants_genotypes.T,
        index=phenotypes.eid,
        columns=term_variants.variant_id,
    )
    return scores, variants_info, variants_genotypes


def process_data(
    scores: pd.DataFrame, variants_info: pd.DataFrame, variants_genotypes: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series]:
    """Process the data."""

    # change -9 to np.nan
    variants_genotypes = variants_genotypes.replace({-9: np.nan})

    # TODO: Duplicate variants for the same term should be carefully removed at
    # source, e.g. by picking the row with the highest hmm_score.

    # remove duplicated columns
    variants_genotypes = variants_genotypes.loc[
        :, ~variants_genotypes.columns.duplicated()
    ]

    # TODO: Regarding the below, I don't understand what this is doing or why...
    # I guess it's counting, looking for cases where alt,alt (2) is more frequent?

    # correcting for the fact that the reference allele is not always the minor allele
    for col in variants_genotypes.columns:
        if 2 * variants_genotypes[col].value_counts().get(0, 0) + variants_genotypes[
            col
        ].value_counts().get(1, 0) < 2 * variants_genotypes[col].value_counts().get(
            2, 0
        ) + variants_genotypes[col].value_counts().get(1, 0):
            variants_genotypes[col] = variants_genotypes[col].replace({0: 2, 2: 0})

    genotypes_weighted2 = variants_genotypes.copy()
    # # sort score and genotypes by index
    # scores.sort_index(inplace=True)
    # genotypes_weighted2.sort_index(inplace=True)
    # # index as integers
    # scores.index = scores.index.astype(int)
    # genotypes_weighted2.index = genotypes_weighted2.index.astype(int)

    # verify if there are duplicates in the index and print if they are
    if scores.index.duplicated().any():
        logging.warning("Duplicate index in scores")
    if genotypes_weighted2.index.duplicated().any():
        logging.warning("Duplicate index in genotypes")
    else:
        logging.info("No duplicates in the index")

    # TODO: We could do this back where we get the data...
    # merge by index
    genotypes_weighted = pd.concat([scores, genotypes_weighted2], axis=1)

    # open phenotype file
    pheno_subset = genotypes_weighted[genotypes_weighted["pheno"] == 1]
    sex_counts = pheno_subset["sex"].value_counts(normalize=True)
    logging.info(
        f"Proportions of M and F within pheno = 1 subset: {sex_counts.to_dict()}"
    )
    # impute missing values
    return genotypes_weighted, scores, genotypes_weighted2, sex_counts


def FATHMM_multiply(
    genotypes_weighted: pd.DataFrame, variants_info: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.Index, pd.DataFrame]:
    variants_info = variants_info.T
    # for those index names which are in the columns of genotypes_weighted,
    # multiply the values in the columns of genotypes_weighted by the values in
    # the hmm_score row of variants_info remove duplicate columns
    genotypes_weighted = genotypes_weighted.loc[
        :, ~genotypes_weighted.columns.duplicated()
    ]
    variants_info = variants_info.T.loc[~variants_info.T.index.duplicated()].T
    genotypes_weighted = genotypes_weighted.loc[~genotypes_weighted.index.duplicated()]
    variants_info = variants_info.loc[~variants_info.index.duplicated()]
    # sort columns in genotypes_weighted and variants_info alphabetically
    genotypes_weighted = genotypes_weighted.reindex(
        sorted(genotypes_weighted.columns), axis=1
    )
    variants_info = variants_info.reindex(sorted(variants_info.columns), axis=1)
    common_cols = genotypes_weighted.columns.intersection(variants_info.columns)
    if "hmm_score" in variants_info.index:
        genotypes_weighted[common_cols] = genotypes_weighted[common_cols].multiply(
            variants_info.loc["hmm_score"], axis=1
        )
    else:
        logging.error("hmm_score row not found in variants_info")
        sys.exit(1)
    # rename columns adding gene names (variantid_GENENAME) from gene row name in variants_info for common_cols
    genotypes_weighted.rename(
        columns={col: f"{col}_{variants_info.loc['gene', col]}" for col in common_cols},
        inplace=True,
    )
    return variants_info, common_cols, genotypes_weighted


def impute_missing_values(genotypes_weighted: pd.DataFrame) -> pd.DataFrame:
    # impute missing values
    genotypes_weighted = genotypes_weighted.dropna(
        thresh=int(genotypes_weighted.shape[0] * 0.1),
        axis=1,
    )
    imputer = SimpleImputer(strategy="most_frequent")
    return pd.DataFrame(
        imputer.fit_transform(genotypes_weighted),
        columns=genotypes_weighted.columns,
        index=genotypes_weighted.index,
    )


def new_columns(
    genotypes_weighted: pd.DataFrame, variants_info: pd.DataFrame
) -> pd.DataFrame:
    # new column which has scores only for affected patients
    genotypes_weighted["score_disease"] = np.where(
        genotypes_weighted["pheno"] == 1, genotypes_weighted["score"], 0
    )
    # remove all columns which have only 0s in them
    genotypes_weighted = genotypes_weighted.loc[
        :, (genotypes_weighted != 0).any(axis=0)
    ]
    return genotypes_weighted


def discretise_genotypes(genotypes_weighted: pd.DataFrame) -> pd.DataFrame:
    genotypes_weighted_discretised = genotypes_weighted.copy()
    genotypes_weighted_discretised["score_disease"] = pd.to_numeric(
        genotypes_weighted_discretised["score_disease"], errors="coerce"
    )
    genotypes_weighted_discretised = genotypes_weighted_discretised.dropna(
        subset=["score_disease"]
    )
    # discretise only non-zero values from score_disease, staring from >0
    genotypes_weighted_discretised.loc[
        genotypes_weighted_discretised["pheno"] == 0, "score_disease"
    ] = np.nan
    # to ensure that 0s are only non-disease individuals, added 1 to the bins
    # Print out binned ranges (not score_disease value counts, because those are only integers)
    binned_ranges = (
        pd.cut(genotypes_weighted_discretised["score_disease"], bins=5)
        .value_counts()
        .sort_index()
    )
    logging.info(f"Binned ranges for score_disease: {binned_ranges}")
    genotypes_weighted_discretised["score_disease"] = (
        pd.cut(
            genotypes_weighted_discretised["score_disease"],
            bins=5,
            labels=False,
            include_lowest=True,
        )
        + 1
    )
    genotypes_weighted_discretised["score_disease"] = genotypes_weighted_discretised[
        "score_disease"
    ].replace(np.nan, 0)
    genotypes_weighted_discretised["score_disease"] = genotypes_weighted_discretised[
        "score_disease"
    ].astype(int)
    discretise_bins = [-np.inf, 0, 1, 3, 5, np.inf]
    for col in genotypes_weighted_discretised.columns:
        if col not in ["score", "score_disease", "pheno", "sex", "suppop"]:
            genotypes_weighted_discretised[col] = pd.to_numeric(
                genotypes_weighted_discretised[col], errors="coerce"
            )
            genotypes_weighted_discretised[col] = pd.cut(
                genotypes_weighted_discretised[col], bins=discretise_bins, labels=False
            )
            genotypes_weighted_discretised[col] = (
                genotypes_weighted_discretised[col].fillna(0).astype(int)
            )
    return genotypes_weighted_discretised


def undersample(
    genotypes_weighted_discretised: pd.DataFrame, sex_counts: pd.Series
) -> pd.DataFrame:
    # if the sex_counts does not equal 0 or 1 (i.e. there are both M and F in the pheno = 1 subset), undersample the pheno = 0 subset to match the sex proportions in the pheno = 1 subset
    pheno_1_subset = genotypes_weighted_discretised[
        genotypes_weighted_discretised["pheno"] == 1
    ]

    pheno_0_subset_undersampled = None

    if not (sex_counts == 0).any() and not (sex_counts == 1).any():
        desired_0_count = max(10000, pheno_1_subset.shape[0])
        pheno_0_subset = genotypes_weighted_discretised[
            genotypes_weighted_discretised["pheno"] == 0
        ]
        desired_0_female = desired_0_count * sex_counts.get("F", 0)
        desired_0_male = desired_0_count * sex_counts.get("M", 0)
        pheno_0_subset = genotypes_weighted_discretised[
            genotypes_weighted_discretised["pheno"] == 0
        ]
        pheno_1_subset = genotypes_weighted_discretised[
            genotypes_weighted_discretised["pheno"] == 1
        ]
        female_0_undersampled = pheno_0_subset[pheno_0_subset["sex"] == "F"].sample(
            n=int(desired_0_female), random_state=42, replace=True
        )
        male_0_undersampled = pheno_0_subset[pheno_0_subset["sex"] == "M"].sample(
            n=int(desired_0_male), random_state=42, replace=True
        )
        pheno_0_subset_undersampled = pd.concat(
            [female_0_undersampled, male_0_undersampled]
        )
    if (sex_counts == 0).any() or (sex_counts == 1).any():
        # undersample the data, so that the number of pheno = 0 samples is equal to desired_0_count
        desired_0_count = max(10000, pheno_1_subset.shape[0])
        pheno_0_subset = genotypes_weighted_discretised[
            genotypes_weighted_discretised["pheno"] == 0
        ]
        pheno_1_subset = genotypes_weighted_discretised[
            genotypes_weighted_discretised["pheno"] == 1
        ]
        pheno_0_subset_undersampled = pheno_0_subset.sample(
            n=desired_0_count, random_state=42, replace=True
        )

    if pheno_0_subset_undersampled is None:
        logging.error("No pheno = 0 samples found.")
        sys.exit(1)

    genotypes_weighted_discretised = pd.concat(
        [pheno_0_subset_undersampled, pheno_1_subset]
    )
    return genotypes_weighted_discretised


def check_missing_values(genotypes_weighted_discretised: pd.DataFrame) -> None:
    if genotypes_weighted_discretised.isnull().values.any():
        logging.error("NaN values found after imputation.")
        sys.exit(1)
    else:
        logging.info("No missing values found after imputation.")


def root(df: pd.DataFrame, class_node: str) -> str:
    """Determine the root node for the Bayesian network."""
    # Create a graph from the correlation matrix.
    # If it will occur to be inefficient in the future consider pivoting the table before running correlation matrix
    nx_graph = nx.from_pandas_adjacency(df.corr().abs())
    eligible_nodes = set(nx_graph.nodes()) - {class_node}
    # Compute Centrality Measures
    degree_centrality = nx.katz_centrality_numpy(nx_graph)
    # Select the Root Node
    root_node = max(eligible_nodes, key=lambda node: degree_centrality[node])
    print(f"Selected Root Node: {root_node}")
    return root_node


def fit_and_plot_bayesian_network(
    df: pd.DataFrame,
    GO_term,
    disease,
    class_node: str,
    root_node: str,
    output_dir: Path,
    genotypes_weighted_discretised: pd.DataFrame,
) -> Tuple[Any, Any]:
    """Fit and plot the Bayesian network model."""
    if df.isnull().values.any():
        # print those rows which have NaN values
        print(df[df.isnull().any(axis=1)].loc[:, df.isnull().any()])
        raise ValueError("Input data contains NaN values")
    else:
        logging.info("No missing values found in the input data.")
    model = bn.structure_learning.fit(
        df, methodtype="tan", class_node=class_node, root_node=root_node
    )
    try:
        model_pruned = bn.independence_test(
            model, df, alpha=0.05, test="g_sq", prune=False
        )
    except ValueError as e:
        logging.error(f"Error during independence test: {e}")
        sys.exit(1)
    # adjmat = model_pruned["adjmat"]
    # adjmat.to_csv(output_dir / "adjmat.csv")
    # bn.adjmat2vec(adjmat, min_weight=1).to_csv(output_dir / "adjmat_vec.csv")
    # adjmat_dict = bn.adjmat2dict(adjmat)
    model_stats = model_pruned["independence_test"]
    model_stats.to_csv(output_dir / "model_stats.csv")

    score_disease_rows = model_stats[
        model_stats.apply(
            lambda row: row.astype(str).str.contains("score_disease").any(), axis=1
        )
    ]
    if score_disease_rows.empty:
        with open(output_dir / "score_disease_rows.txt", "w") as f:
            f.write("No variants with causative relationships with affected outliers\n")
        raise ValueError(
            "No variants with causative relationships with affected outliers"
        )

    score_disease_rows = score_disease_rows.sort_values(by="p_value", ascending=True)
    score_disease_rows["GO_term"] = GO_term
    # header = f"variants for {GO_term} term and {disease} Phecode which have
    # causal relationships with affected outliers, sorted by p values of gsquare
    # independence test\n"
    with open(output_dir / "score_disease_rows.tsv", "w") as f:
        # f.write(header)
        score_disease_rows.to_csv(f, sep="\t", index=False)

    scores = bn.structure_scores(model_pruned, df, scoring_method=["bic"], verbose=3)
    dot_graph = bn.plot_graphviz(model_pruned, edge_labels="pvalue")

    # TODO: Why is this dot file not being written to the output directory?
    if dot_graph is not None:
        with open("/tmp/output.dot", "w") as f:
            f.write(str(dot_graph.source))

    nx_graph = nx.drawing.nx_agraph.read_dot("/tmp/output.dot")
    nodes_with_path = [
        node for node in nx_graph.nodes() if nx.has_path(nx_graph, class_node, node)
    ]
    print(f"Nodes with a path to '{class_node}': {nodes_with_path}")

    # Remove nodes without any path to the class node from the graph (they are not informative)
    nodes_to_remove = [
        node for node in nx_graph.nodes() if not nx.has_path(nx_graph, class_node, node)
    ]
    print(f"Nodes to remove: {nodes_to_remove}")
    nx_graph.remove_nodes_from(nodes_to_remove)

    # Save all paths from the class node to the nodes in the graph for each node
    all_paths = {}
    for node in nx_graph.nodes():
        all_paths[node] = list(
            nx.all_simple_paths(nx_graph, source=class_node, target=node)
        )

    genotypes_weighted_discretised = genotypes_weighted_discretised.round(0).astype(int)

    # For each path (for each node), count how many combinations exist in genotypes_weighted_discretised with non-zero values for all nodes in the path
    all_paths_counts_nonzero_score_disease = {}
    for node in nx_graph.nodes():
        all_paths_counts_nonzero_score_disease[node] = {}
        for path in all_paths[node]:
            all_paths_counts_nonzero_score_disease[node][str(path)] = (
                genotypes_weighted_discretised.loc[
                    (genotypes_weighted_discretised[list(path)] != 0).all(axis=1)
                ].shape[0]
            )

    # Additionally, create new json which will count the same paths as above, but instead of counting paths where all nodes are non-zero, count paths where all nodes are non-zero BUT except from score_disease which is zero
    all_paths_counts_no_score_disease = {}
    for node in nx_graph.nodes():
        all_paths_counts_no_score_disease[node] = {}
        for path in all_paths[node]:
            nodes_except_score_disease = [
                node for node in list(path) if node != "score_disease"
            ]
            all_paths_counts_no_score_disease[node][str(path)] = (
                genotypes_weighted_discretised.loc[
                    (
                        genotypes_weighted_discretised[nodes_except_score_disease] != 0
                    ).all(axis=1)
                    & (genotypes_weighted_discretised["score_disease"] == 0)
                ].shape[0]
            )

    # Remove nodes where all paths to class node have value 0 in all paths to class node (from all_paths_counts_nonzero_score_disease)
    nodes_to_remove = [
        node
        for node in nx_graph.nodes()
        if all(
            all_paths_counts_nonzero_score_disease[node][str(path)] == 0
            for path in all_paths[node]
        )
    ]
    print(f"Nodes to remove: {nodes_to_remove}")
    nx_graph.remove_nodes_from(nodes_to_remove)

    # Create contingency table for each path in all_paths_counts_no_score_disease and all_paths_counts_nonzero_score_disease, including overall counts of score_disease 0 vs non-zero and calculate the p-value for each path
    contingency_table = {}
    for node in nx_graph.nodes():
        contingency_table[node] = {}
        for path in all_paths[node]:
            try:
                # Contingency table with score_disease 0 vs non-zero for paths, compared to score_disease 0 vs non-zero overall
                contingency_table[node][str(path)] = {
                    "table": [
                        [
                            all_paths_counts_nonzero_score_disease[node][str(path)],
                            all_paths_counts_no_score_disease[node][str(path)],
                        ],
                        [
                            genotypes_weighted_discretised.loc[
                                genotypes_weighted_discretised["score_disease"] != 0
                            ].shape[0]
                            - all_paths_counts_nonzero_score_disease[node][str(path)],
                            genotypes_weighted_discretised.loc[
                                genotypes_weighted_discretised["score_disease"] == 0
                            ].shape[0]
                            - all_paths_counts_no_score_disease[node][str(path)],
                        ],
                    ]
                }
                odds_ratio, p_value = scipy.stats.fisher_exact(
                    contingency_table[node][str(path)]["table"]
                )
                contingency_table[node][str(path)]["odds_ratio"] = odds_ratio
                contingency_table[node][str(path)]["p_value"] = p_value
            except ValueError as e:
                logging.error(
                    f"Error creating contingency table for node '{node}' and path '{path}': {e}"
                )
    # with open(output_dir / "contingency_table.json", "w") as f:
    #    json.dump(contingency_table, f)
    # save also as a tsv dataframe
    contingency_table_df = pd.DataFrame(
        [
            (
                node,
                path,
                contingency_table[node][str(path)]["table"][0][0],
                contingency_table[node][str(path)]["table"][0][1],
                contingency_table[node][str(path)]["table"][1][0],
                contingency_table[node][str(path)]["table"][1][1],
                contingency_table[node][str(path)]["odds_ratio"],
                contingency_table[node][str(path)]["p_value"],
            )
            for node in nx_graph.nodes()
            for path in all_paths[node]
        ],
        columns=[
            "node",
            "path",
            "combination_present_disease_present",
            "combination_present_disease_absent",
            "combination_absent_disease_present",
            "combination_absent_disease_absent",
            "odds_ratio",
            "p_value",
        ],
    )
    contingency_table_df.sort_values(by="p_value", ascending=True, inplace=True)
    # remove rows which have only two elements in path
    # contingency_table_df = contingency_table_df[contingency_table_df['path'].apply(lambda x: len(x) > 2)]
    contingency_table_df.to_csv(
        output_dir / "contingency_table.tsv", sep="\t", index=False
    )

    dot_graph = nx.drawing.nx_agraph.to_agraph(nx_graph)

    # show edge labels (pvalues) only for those nodes which are directly connected to the class node through one edge
    for edge in dot_graph.edges():
        if edge[0] == class_node or edge[1] == class_node:
            continue
        else:
            dot_graph.get_edge(edge[0], edge[1]).attr["label"] = ""

    # colour those nodes which have p-value < 0.05 in red
    # Assign colors to nodes based on the number of significant paths
    for node in nx_graph.nodes():
        significant_paths = [
            path
            for path in all_paths[node]
            if contingency_table[node][str(path)]["p_value"] is not None
            and contingency_table[node][str(path)]["p_value"] < 0.05
            and len(path) > 2
        ]
        if significant_paths:
            for path in significant_paths:
                for i in range(len(path) - 1):
                    dot_graph.get_edge(path[i], path[i + 1]).attr["color"] = "red"
                    dot_graph.get_edge(path[i], path[i + 1]).attr["penwidth"] = 2.0

    # colour nodes which have p-value < 0.05 in lightpink
    for node in nx_graph.nodes():
        if any(
            score_disease_rows[
                (score_disease_rows["source"] == node)
                | (score_disease_rows["target"] == node)
            ]["p_value"]
            < 0.05
        ):
            dot_graph.get_node(node).attr["fillcolor"] = "lightpink"
            dot_graph.get_node(node).attr["style"] = "filled"

    if dot_graph is not None:
        dot_graph.draw(path=str(output_dir / "bayesian_network_new.svg"), prog="dot")
        # colour score_disease node in the svg file in lightgrey
        with open(output_dir / "bayesian_network_new.svg", "r") as f:
            svg = f.read()
        svg = re.sub(
            r'(<title>score_disease</title>\s*<ellipse[^>]*fill=")[^"]*',
            r"\1lightgrey",
            svg,
        )
        with open(output_dir / "bayesian_network_new.svg", "w") as f:
            f.write(svg)
    else:
        logging.error("Failed to create dot_graph, skipping rendering.")
    print("BIC Score:", scores["bic"])

    # save to dot file
    dot_graph.write(str(output_dir / "bayesian_network_new.dot"))

    return model_stats, dot_graph


network_analysis_bp = Blueprint("network_analysis", __name__)


logger = logging.getLogger(__name__)


@network_analysis_bp.route("/network_analysis/<disease_code>/<term>")
def do_a_thing(disease_code: str, term: str) -> Response:
    """Perform network analysis for a given disease and GO term."""

    GO_term = term
    disease = disease_code

    flush = bool(request.args.get("flush", False))
    logger.info(f"Flush: {flush}")

    services = current_app.extensions["nomaly_services"]
    phenotype_service = services.phenotype
    nomaly_score_service = services.nomaly_score
    genotype_service = services.genotype

    # Get run version and ancestry from the session
    # TODO: Implement run version!
    run_version = session.get("run_version", "Run-v1")
    ancestry = session.get("ancestry", "EUR")

    output_dir = Config.NETWORK_ANALYSIS_DIR / f"{ancestry}-{disease}-{GO_term}"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Processing GO term: {GO_term}, Disease: {disease}")
    logger.info(f"Output is written to: {output_dir}")

    scores, variants_info, variants_genotypes = read_files(
        GO_term,
        disease,
        phenotype_service,
        nomaly_score_service,
        genotype_service,
        ancestry,
    )

    logger.info(f"Variants info shape: {variants_info.shape}")
    genotypes_weighted, scores, genotypes_weighted2, sex_counts = process_data(
        scores, variants_info, variants_genotypes
    )

    variants_info, common_cols, genotypes_weighted = FATHMM_multiply(
        genotypes_weighted, variants_info
    )

    genotypes_weighted = impute_missing_values(genotypes_weighted)
    if genotypes_weighted.isnull().values.any():
        logging.error("NaN values found after imputation.")
        sys.exit(1)
    genotypes_weighted = new_columns(genotypes_weighted, variants_info)
    # genotypes_weighted = genotypes_weighted.loc[:, genotypes_weighted.columns[genotypes_weighted.loc[(genotypes_weighted['pheno'] == 1), (genotypes_weighted != 0).any()].any()]]
    genotypes_weighted_discretised = discretise_genotypes(genotypes_weighted)
    genotypes_weighted_discretised = genotypes_weighted_discretised.drop(
        ["score"], axis=1
    )
    # all columns as categorical
    genotypes_weighted_discretised = genotypes_weighted_discretised.astype("category")
    # genotypes_weighted_discretised = genotypes_weighted_discretised.drop(['pheno'], axis=1)
    genotypes_weighted_discretised.index.name = None
    genotypes_weighted_discretised = undersample(
        genotypes_weighted_discretised, sex_counts
    )
    genotypes_weighted_discretised = genotypes_weighted_discretised.loc[
        :, (genotypes_weighted_discretised != 0).any(axis=0)
    ]

    logger.info(
        f"Genotypes weighted discretised shape: {genotypes_weighted_discretised.shape}"
    )

    logger.info(genotypes_weighted_discretised["sex"].value_counts())
    logger.info(genotypes_weighted_discretised["sex"].isnull().sum())
    logger.info(genotypes_weighted_discretised.shape)
    logger.info("Proportion of M and F in the pheno = 1 subset")
    logger.info(
        genotypes_weighted_discretised[genotypes_weighted_discretised["pheno"] == 1][
            "sex"
        ].value_counts()
    )
    logger.info("Proportion of M and F in the pheno = 0 subset")
    logger.info(
        genotypes_weighted_discretised[genotypes_weighted_discretised["pheno"] == 0][
            "sex"
        ].value_counts()
    )
    check_missing_values(genotypes_weighted_discretised)

    genotypes_weighted_discretised = genotypes_weighted_discretised.drop(
        ["sex", "pheno"], axis=1
    )

    genotypes_weighted_discretised = genotypes_weighted_discretised.apply(
        pd.to_numeric, errors="ignore"
    )

    logger.info(
        f"Genotypes weighted discretised shape: {genotypes_weighted_discretised.shape}"
    )

    root_node = root(genotypes_weighted_discretised, "score_disease")
    logger.info(f"Root node: {root_node}")
    model_stats, dot_graph = fit_and_plot_bayesian_network(
        genotypes_weighted_discretised,
        GO_term=GO_term,
        disease=disease,
        class_node="score_disease",
        root_node=root_node,
        output_dir=output_dir,
        genotypes_weighted_discretised=genotypes_weighted_discretised,
    )

    return Response(
        "Network analysis completed successfully",
        status=200,
        mimetype="text/plain",
    )


def main():
    do_a_thing("100001", "100001")


if __name__ == "__main__":
    main()
