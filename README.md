
Genotyping project
-----

Aims:

Capture the genotype calls for donors from ATRC-000078 in the Atreca production database
The script/program to upload any future donor allele calls to the database

Description: 

1. Upload the reference allele sequences to the database:
 - The germline reference sequences are stored in the binf_germline_seq table.

 - We need to upload the allele sequences (e.g. the germline reference sequences that were used to genotype ATRC-000078) as a separate release to binf_germline_seq table. There is an "allele" field which we are not actively using right now - however, for this, we'd want to use that field.

2. Design a table to store the genotype calls:
       - The genotype information is associated with each donor
       -  The table should have at minimum the following columns:
    + primary_key
    + donor_id: foreign key pointing to study_sample_donor id (The donor_identifier and field in binf_pairs_rpt that you see (e.g. DID-95007) is the "identifier" column in this table)
    + allele_id:  foreign key pointing to binf_germline_seq id
     
      - In addition, these columns could be useful:
    + method (e.g. TIgGER / Partis)
    + number of supporting contigs
    + confidence score (for Partis)
    + date_created
    + ...

3. Create the table
4. Populating the table with the genotype calls that you generated for ATRC-000078
5. The script/program: 
  The idea is that in the future, if we'd like to genotype any additional donors, we'd be able to relatively straightforwardly do it.
