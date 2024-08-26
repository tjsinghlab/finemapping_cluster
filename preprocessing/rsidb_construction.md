# RSID manipulations & lookup

Starting with file `GCF_000001405.25_slim.tsv` from Max-Drabkin:
```
CHROM POS ID REF ALT
NC_000001.10 10001 rs1570391677 T A,C
NC_000001.10 10002 rs1570391692 A C
```

## Separate CHROM into various components

[HGVS Reference Sequence](https://hgvs-nomenclature.org/stable/background/refseq/)

e.g., **NC_000001.10**: prefix (**NC**), chromosome number (**1**), version (**10**)
```
$ LC_ALL=C awk 'BEGIN {OFS="\t"}
NR==1 {print $0, "prefix", "chromosome", "version"; next}
{
  split($1, a, "[_.]");
  $0 = $0 OFS a[1] OFS a[2]+0 OFS a[3];
  print
}' GCF_000001405.25_slim.tsv > GCF_000001405.25_slim_expanded.tsv
```

```
$ ls -lh GCF_000001405.25_slim_expanded.tsv
-rw-rw-r--  1 sfriedman nslab  51G Aug 24 22:33 GCF_000001405.25_slim_expanded.tsv
$ head GCF_000001405.25_slim_expanded.tsv
CHROM POS ID REF ALT prefix chromosome version
NC_000001.10 10001 rs1570391677 T A,C NC 1 10
NC_000001.10 10002 rs1570391692 A C NC 1 10
```

## Figure out which chr builds are represented in this file
```
$ awk '{print $NF}' GCF_000001405.25_slim_expanded.tsv | sort | uniq -c | sort -rn > version_counts.txt
$ cat version_counts.txt 
395532324 11  # GRCh38/hg38 genome assembly
336060558 10  # GRCh37/hg19 genome assembly
235738966 9
59858844 13
47427291 8
32981289 1
6496723 2
4314669 3
```

## Break up big file into smaller ones per build
```
$ awk -F'\t' 'NR==1{header=$0; next} {print > "version_"$8".txt"} > END {for(f in FILENAME) print header > f}' GCF_000001405.25_slim_expanded.tsv
```

```
$ ls -lh version_*.txt
-rw-rw-r--  1 sfriedman nslab  16G Aug 24 23:41 version_10.txt
-rw-rw-r--  1 sfriedman nslab  18G Aug 24 23:41 version_11.txt
-rw-rw-r--  1 sfriedman nslab 2.7G Aug 24 23:41 version_13.txt
-rw-rw-r--  1 sfriedman nslab 1.6G Aug 24 23:41 version_1.txt
-rw-rw-r--  1 sfriedman nslab 325M Aug 24 23:41 version_2.txt
-rw-rw-r--  1 sfriedman nslab 221M Aug 24 23:41 version_3.txt
-rw-rw-r--  1 sfriedman nslab 2.1G Aug 24 23:41 version_8.txt
-rw-rw-r--  1 sfriedman nslab  11G Aug 24 23:41 version_9.txt
```

## Sort the desired build files by RSID prior to writing to SQLite DB

### GRCh37/hg19 genome assembly
```
$ head -n 1 GCF_000001405.25_slim_expanded.tsv > version_10_sorted.tsv
$ LC_ALL=C sort -k3,3V <(tail -n +2 version_10.txt) >> version_10_sorted.tsv
```

### GRCh38/hg38 genome assembly
```
$ head -n 1 GCF_000001405.25_slim_expanded.tsv > version_11_sorted.tsv
$ LC_ALL=C sort -k3,3V <(tail -n +2 version_11.txt) >> version_11_sorted.tsv
```

# Setup SQLite DB
```
-- Start SQLite
$ sqlite3 rsid.db
```

### GRCh37/hg19 genome assembly
```
-- Set the mode to tabs
.mode tabs

-- Import the data
.import /gpfs/commons/home/sfriedman/version_10_sorted.tsv hg19

-- Add integer rsid column
ALTER TABLE hg19 ADD COLUMN rsid INTEGER;

UPDATE hg19 SET rsid = CAST(SUBSTR(ID, 3) AS INTEGER) WHERE ID LIKE 'rs%';

-- Check schema
.schema hg19

-- Create table index
CREATE INDEX idx_hg19_rsid ON hg19 (rsid);

-- Try command of interest
SELECT ID, chromosome, POS 
FROM hg19 
WHERE rsid >= 200000 AND rsid <= 210000;
```

### GRCh38/hg38 genome assembly
```
-- Set the mode to tabs
.mode tabs

-- Import the data
.import /gpfs/commons/home/sfriedman/version_11_sorted.tsv hg38

-- Add integer rsid column
ALTER TABLE hg38 ADD COLUMN rsid INTEGER;
UPDATE hg38 SET rsid = CAST(SUBSTR(ID, 3) AS INTEGER) WHERE ID LIKE 'rs%';

-- Check schema
.schema hg38

-- Create table index
CREATE INDEX idx_hg38_rsid ON hg38 (rsid);

-- Try command of interest
SELECT ID, chromosome, POS 
FROM hg38
WHERE rsid >= 200000 AND rsid <= 200100;
```

