
genotypes.hdf5
* 488377 * 83011, 2.2GB
* genotypes
* individuals
* variants (identifier by chrom:pos:ref:alt from bim file)

float16_scores.hdf5
* 486145 * 12936, 3.5GB
* scores
* individuals
* ontterms

variants 
* chrom:pos:ref:alt (bim)
* rsid (bim)
+------------+--------------+------+-----+---------+----------------+
| variant_id | varchar(64)  | YES  | MUL | NULL    |                |
| gene       | varchar(64)  | YES  |     | NULL    |                |
| seqid      | varchar(128) | YES  |     | NULL    |                |
| aa         | varchar(10)  | YES  |     | NULL    |                |
| hmm        | varchar(30)  | YES  |     | NULL    |                |
| hmm_pos    | int          | YES  |     | NULL    |                |
| wild       | varchar(10)  | YES  |     | NULL    |                |
| mutant     | varchar(10)  | YES  |     | NULL    |                |
| sf         | int          | YES  |     | NULL    |                |
| sfe        | double       | YES  |     | NULL    |                |
| fa         | int          | YES  |     | NULL    |                |
| fae        | double       | YES  |     | NULL    |                |
+------------+--------------+------+-----+---------+----------------+

individuals
* population
* sex

icd10
* "ICD10"
* "meaning"
* "node_id"
* "parent_id"
* "selectable"

individuals_icd10
* individual_id
* icd10_id

phecode
* "phecode"
* "description"
* "sex"
* "rollup"
* "leaf"
* "groupnum"
* "group"
* "color"
* "phecode_exclude_range"
* "phecode_exclude_phenotypes"

icd10_phecode
* icd10_id
* phecode_id

phenotypes
* "phenotype_id"
* description
* sex
* affected_individuals
* unaffected_individuals

stats
* ontterms
* phenotype_id
* pvalues

ontterms
* "term_id"
* name

term2snps
* "term_id"
* "variant_id"

