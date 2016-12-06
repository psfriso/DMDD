drop database if exists dmdd_annotations;
create database dmdd_annotations;

use dmdd_annotations;

drop table if exists impc_release;
create table impc_release(
    id int unsigned auto_increment,
    gene_marker_id varchar(128) default null,
    gene_id varchar(128) default null,
    center varchar(16) default null,
    algo_1 varchar(128) default null,
    sex varchar(32) default null,
    zygosity varchar(32) default null,
    allele_marker varchar(128) default null,
    gene_line varchar(64) default null,
    mutation_type varchar(512) default null,
    strain_marker varchar(128) default null,
    strain varchar(128) default null,
    project varchar(128) default null,
    project_name varchar(512) default null,
    pipeline varchar(256) default null,
    pipeline_id varchar(256) default null,
    procedure_id varchar(64) default null,
    procedure_name varchar(256) default null,
    parameter_id varchar(64) default null,
    parameter varchar(128) default null,
    mp_term1 varchar(64) default null,
    mp_term1_desc varchar(512) default null,
    mp_term2 varchar(64) default null,
    mp_term2_desc varchar(512) default null,
    p_value double default null,
    e_size double default null,
    e_size2 double default null,
    statistical_model varchar(512) default null,
    algo_2 varchar(128) default null,
    primary key(id),
    index(gene_marker_id),
    index(gene_id),
    index(allele_marker)
) engine = innodb;


drop table if exists mgi_release;
create table mgi_release(
    id int unsigned auto_increment,
    allele_accession_id varchar(64) default null,
    allele_symbol  varchar(128) default null,
    allele_name varchar(128) default null,
    allele_type varchar(128) default null,
    allele_attribute varchar(128) default null,
    pubmed_id int unsigned default null,
    marker_accession_id varchar(64) default null,
    marker_symbol varchar(64) default null,
    marker_refseq_id varchar(64) default null,
    marker_ensembl_id varchar(64) default null,
    high_leve_mp varchar(128) default null,
    synonyms_mp varchar(2056) default null,
    primary key (id),
    index( allele_accession_id ),
    index( marker_accession_id ),
    index( marker_symbol )
) engine = innodb;

-- permisions for users

grant all on dmdd_annotations.* to 'psfriso'@'localhost' identified by 'PASSWORD';

-- loading data

truncate impc_release;

LOAD DATA LOCAL INFILE '/home/psfriso/Desktop/sid_data/ALL_genotype_phenotype_dr5.0_beta.txt'
INTO TABLE dmdd_annotations.impc_release
FIELDS TERMINATED BY '\t'
ENCLOSED BY ''
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(gene_marker_id, gene_id, center, algo_1, sex, zygosity, gene_line, allele_marker,mutation_type,
strain_marker , strain, project, project_name, pipeline, pipeline_id, procedure_id, procedure_name,
parameter_id, parameter, mp_term1, mp_term1_desc, mp_term2, mp_term2_desc, p_value, e_size, e_size2,
statistical_model,algo_2)
set id = NULL;

show warnings;

truncate mgi_release;
LOAD DATA LOCAL INFILE '/home/psfriso/Desktop/sid_data/MGI_PhenotypicAllele.rpt'
INTO TABLE dmdd_annotations.mgi_release
FIELDS TERMINATED BY '\t'
ENCLOSED BY ''
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(allele_accession_id, allele_symbol, allele_name, allele_type, allele_attribute,
pubmed_id, marker_accession_id, marker_symbol,  marker_refseq_id,
marker_ensembl_id, high_leve_mp, synonyms_mp )
set id = NULL;
