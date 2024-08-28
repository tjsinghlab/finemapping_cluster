# RSID SQLite Database Construction

Start with dbSNP files:
[SNP Files](https://ftp.ncbi.nlm.nih.gov/snp/organisms/)

## Download files
```{bash}
wget -O hg19_00-All.vcf.gz \
    https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz  # GRCh37/hg19
wget -O hg38_00-All.vcf.gz \
    https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz   # GRCh38/hg38
```

## Trim and show first few rows
```{bash}
$ zcat 00-All.vcf.gz | tail -n +57 | cut -f1-7 > GRCh38.tsv

$ head GRCh37.tsv 
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
1	10019	rs775809821	TA	T	.	.
1	10039	rs978760828	A	C	.	.
1	10043	rs1008829651	T	A	.	.
1	10051	rs1052373574	A	G	.	.
1	10051	rs1326880612	A	AC	.	.
```

# Setup SQLite DB
```
-- Start SQLite
$ sqlite3 rsid.db
```

### GRCh38/hg38 genome assembly

Replace `hg38` with `hg19` for GRCh37/hg19 build.
```
-- Set the mode to tabs
.mode tabs

-- Import the data
.import /gpfs/commons/home/sfriedman/dbsnp/GRCh38.tsv hg38

-- Update table names
ALTER TABLE hg38 RENAME COLUMN "#CHROM" TO "CHROM";

-- Add integer RSID column
ALTER TABLE hg38 ADD COLUMN RSID INTEGER;
UPDATE hg38 SET RSID = CAST(SUBSTR(ID, 3) AS INTEGER) WHERE ID LIKE 'rs%';

-- Check schema
.schema hg38

-- Create table index
CREATE INDEX idx_hg38_rsid ON hg38 (RSID);

-- Try command of interest
SELECT ID, chromosome, POS 
FROM hg38 
WHERE rsid >= 200000 AND rsid <= 210000;
```
