import pandas as pd


class NomalyDataService:
    """Service for accessing Nomaly variant mapping data."""

    def __init__(self, variants_file_path):
        """Initialize the service with the path to the variants file."""
        self.variants_file_path = variants_file_path
        self._df = None
        self._variant_map_nomaly = None
        self._variant_map_plink = None

    @property
    def df(self):
        """Lazy load the variant mapping dataframe."""
        if self._df is None:
            self._df = pd.read_csv(self.variants_file_path, sep="\t")
        return self._df

    def get_variant_info_nomaly(self, nomaly_variant):
        """Get information about a specific variant by its nomaly_variant ID."""
        if self._variant_map_nomaly is None:
            self._variant_map_nomaly = dict(
                zip(self.df["nomaly_variant"], self.df.to_dict("records"))
            )
        return self._variant_map_nomaly.get(nomaly_variant)

    def get_variant_info_plinky(self, plink_variant):
        """Get information about a specific variant by its plink_variant ID."""
        if self._variant_map_plink is None:
            self._variant_map_plink = dict(
                zip(self.df["CHR_BP_A1_A2"], self.df.to_dict("records"))
            )
        return self._variant_map_plink.get(plink_variant)

    def get_variants_by_gene(self, gene_id):
        """Get all variants for a specific gene."""
        return self.df[self.df["gene_id"] == gene_id]

    def map_variant_to_rsid(self, nomaly_variant):
        """Map a nomaly_variant ID to its RSID."""
        info = self.get_variant_info_nomaly(nomaly_variant)
        return info["RSID"] if info else None
