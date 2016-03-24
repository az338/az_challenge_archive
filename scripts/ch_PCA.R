
library(ade4)

DATA= '~/ucc-fileserver/AZ_challenge_data/'
CMP_DATA = file.path(DATA,'dan_descriptors/')
CELL_DATA = file.path('krishna_descriptors/')

## cell line descriptors
cell_disease_binary = read.csv(file.path(CELL_DATA,'cellline_disease_binary_fingerprints.csv'))
cell_disease_inference = read.csv(file.path(CELL_DATA,'cellline_disease_inference_fingerprints.csv'))
cell_fingerprints = read.csv(file.path(CELL_DATA,'cell_line_fingerprints.csv'))
cell_gene_dist = read.csv(file.path(CELL_DATA,'cell_line_gene_distance_fingerprints.csv'))

## compound descriptors

# single compound descriptors
SC_DATA = file.path(CMP_DATA,'compounds/')
cmp_annoTargets = read.csv(file.path(SC_DATA,'dan_compound_annotated_targets.csv'))
cmp_predTargets= read.csv(file.path(SC_DATA,'dan_compound_predicted_targets.csv'))
cmp_pathways = read.csv(file.path(SC_DATA,'dan_compound_signalink_pathways.csv'))
cmp_STN = read.csv(file.path(SC_DATA,'dan_compound_signalink_targets_and_neighbours.csv'))
cmp_fingerprints = read.csv(file.path(SC_DATA,'dan_compound_structural_fingerprints.csv'))





    