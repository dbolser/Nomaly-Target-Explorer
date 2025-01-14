# Nomaly Disease Browser (UKBB edition)

## Installation

```
git clone https://github.com/danbolser/nomaly-disease-browser.git
cd nomaly-disease-browser
```

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

And don't forget...

```bash
cp .env.example .env
```

## Configuration

Edit the .env file with your database credentials (see .env.example)

Look at (config.py)[config.py] for data configuration.

## Authentication

See (schema.sql)[schema.sql] for database schema and auth.py for authentication implementation details.

```sql
-- Create admin user
INSERT INTO users2 (username, password, email, is_active)
VALUES ('admin', 'temporary_password', 'admin@example.com', TRUE);

-- Get the admin user's ID
SET @admin_id = LAST_INSERT_ID();

-- Grant admin full access
INSERT INTO user_permissions (user_id, allowed_paths)
VALUES (@admin_id, '*');

-- Create a limited user with specific phecode access
INSERT INTO users2 (username, password, email, is_active)
VALUES ('phecode_user', 'temp_password', 'phecode@example.com', TRUE);

-- Grant access to specific phecodes
INSERT INTO user_permissions (user_id, allowed_paths)
VALUES (LAST_INSERT_ID(), '705,695,756,256,615');
```


## Running the app

```bash
export FLASK_APP=app.py
export FLASK_ENV=development
flask run
```


## Notes on the input data

# Nomaly Disease Browser (UKBB edition)


## Installation

```
git clone https://github.com/danbolser/nomaly-disease-browser.git
cd nomaly-disease-browser
pip install -r requirements.txt
```


## Configuration

Edit the .env file with your database credentials (see .env.example)

Look at (config.py)[config.py] for data configuration.


## Authentication

See (schema.sql)[schema.sql] for database schema and auth.py for authentication implementation details.

```sql
-- Create admin user
INSERT INTO users2 (username, password, email, is_active)
VALUES ('admin', 'temporary_password', 'admin@example.com', TRUE);

-- Get the admin user's ID
SET @admin_id = LAST_INSERT_ID();

-- Grant admin full access
INSERT INTO user_permissions (user_id, allowed_paths)
VALUES (@admin_id, '*');

-- Create a limited user with specific phecode access
INSERT INTO users2 (username, password, email, is_active)
VALUES ('phecode_user', 'temp_password', 'phecode@example.com', TRUE);

-- Grant access to specific phecodes
INSERT INTO user_permissions (user_id, allowed_paths)
VALUES (LAST_INSERT_ID(), '705,695,756,256,615');
```



## Data information:

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

