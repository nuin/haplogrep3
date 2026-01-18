-- MitoMaster DuckDB Schema
-- Portable database for mtDNA variant frequencies

-- Core genome records (minimal, no raw sequence for distribution)
CREATE TABLE IF NOT EXISTS genomes (
    accession VARCHAR PRIMARY KEY,
    organism VARCHAR,
    taxonomy VARCHAR,
    length INTEGER,
    haplogroup VARCHAR,
    haplogroup_quality DOUBLE,
    is_human BOOLEAN,
    country VARCHAR,
    collection_date VARCHAR,
    processed_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Variants table (main query target)
CREATE TABLE IF NOT EXISTS variants (
    id INTEGER PRIMARY KEY,
    accession VARCHAR NOT NULL,
    position INTEGER NOT NULL,
    ref VARCHAR NOT NULL,
    alt VARCHAR NOT NULL,
    mutation_type VARCHAR,
    locus VARCHAR,
    gene_type VARCHAR,
    amino_acid_change VARCHAR,
    FOREIGN KEY (accession) REFERENCES genomes(accession)
);

-- Pre-computed frequency stats (for fast queries)
CREATE TABLE IF NOT EXISTS variant_frequencies (
    position INTEGER NOT NULL,
    ref VARCHAR NOT NULL,
    alt VARCHAR NOT NULL,
    total_count INTEGER DEFAULT 0,
    total_frequency DOUBLE DEFAULT 0.0,
    human_count INTEGER DEFAULT 0,
    human_frequency DOUBLE DEFAULT 0.0,
    PRIMARY KEY (position, ref, alt)
);

-- Haplogroup-specific frequencies
CREATE TABLE IF NOT EXISTS variant_haplogroup_freq (
    position INTEGER NOT NULL,
    ref VARCHAR NOT NULL,
    alt VARCHAR NOT NULL,
    haplogroup VARCHAR NOT NULL,
    count INTEGER DEFAULT 0,
    frequency DOUBLE DEFAULT 0.0,
    PRIMARY KEY (position, ref, alt, haplogroup)
);

-- Gene map (static reference for mtDNA)
CREATE TABLE IF NOT EXISTS gene_map (
    name VARCHAR PRIMARY KEY,
    start_pos INTEGER NOT NULL,
    end_pos INTEGER NOT NULL,
    gene_type VARCHAR NOT NULL,
    strand CHAR(1),
    product VARCHAR
);

-- Create indexes for fast queries
CREATE INDEX IF NOT EXISTS idx_variants_position ON variants(position);
CREATE INDEX IF NOT EXISTS idx_variants_accession ON variants(accession);
CREATE INDEX IF NOT EXISTS idx_variants_locus ON variants(locus);
CREATE INDEX IF NOT EXISTS idx_genomes_haplogroup ON genomes(haplogroup);
CREATE INDEX IF NOT EXISTS idx_genomes_organism ON genomes(organism);
CREATE INDEX IF NOT EXISTS idx_freq_position ON variant_frequencies(position);
CREATE INDEX IF NOT EXISTS idx_hg_freq_position ON variant_haplogroup_freq(position);
CREATE INDEX IF NOT EXISTS idx_hg_freq_haplogroup ON variant_haplogroup_freq(haplogroup);

-- Sequence for variant IDs
CREATE SEQUENCE IF NOT EXISTS variant_id_seq START 1;
