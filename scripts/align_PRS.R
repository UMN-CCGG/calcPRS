# load packages
library(tidyverse)

# read in arguments
args = commandArgs(trailing=TRUE)
prs_model = args[1]
bim_file = args[2]
output = args[3]

chr <- str_remove(string = bim_file, pattern = "^.*_") %>% str_remove(pattern= ".bim$")

# select model -- this has all variants listed
model <- read_table(prs_model, comment = "#")
bim <- read_table(bim_file, col_names = c('chr', 'id', 'cm', 'bp', 'A1', 'A2'), col_types = list('n', 'c', 'n', 'n', 'c', 'c'))

joined_info <- bim %>%
	left_join(model, by = c("chr" = "chr_name", "bp" = "chr_position")) %>%
	filter(!is.na(effect_allele))


# take complements in BIM

joined_info_comp <- joined_info %>%
	mutate(A1_comp = case_when(A1 == "C" ~ "G",
							   A1 == "T" ~ "A",
							   A1 == "A" ~ "T",
							   A1 == "G" ~ "C"),
		   A2_comp = case_when(A2 == "C" ~ "G",
		   	                   A2 == "T" ~ "A",
		   	                   A2 == "A" ~ "T",
		   	                   A2 == "G" ~ "C"))

# Get sites that match alleles between genotypes and model
matching_dat <- joined_info_comp %>%
	filter(A1 == effect_allele & A2 == other_allele)

# Get sites that have swapped alleles between genotypes and model
swapped_allele_dat <- joined_info_comp %>%
	filter(A1 == other_allele & A2 == effect_allele)

# Get sites that need to be strand flipped in the model file 
flipped_allele_dat <- joined_info_comp %>%
	filter(A1_comp == effect_allele & A2_comp == other_allele)

# get sites that need to be flipped and swapped
flipped_swapped_dat <- joined_info_comp %>%
	filter(A1_comp == other_allele & A2_comp == effect_allele)

# Create new model file with flipped/swapped 
# matching dat is good to go with no further editing.  Just grab columns I need though
matching_dat_to_merge <- matching_dat %>%
	select(chr, id, cm, bp, A1, A2, effect_allele, other_allele, effect_weight)

#swapped_allele_dat need to have effect allele and other allele swapped and effect size to be multiplied by negative 1
swapped_allele_dat_to_merge <- swapped_allele_dat %>%
	mutate(effect_allele2 = other_allele,
		   other_allele2 = effect_allele,
		   effect_weight2 = effect_weight * -1) %>%
	select(chr, id, cm, bp, A1, A2, "effect_allele" = effect_allele2, "other_allele" = other_allele2, "effect_weight" = effect_weight2)

# flipped_allele_dat need to have effect_alelle set to A1 and other_allele set to A2 (they're complements to the bim genotypes)
flipped_allele_dat_to_merge <- flipped_allele_dat %>%
	mutate(effect_allele = A1,
		   other_allele = A2) %>%
	select(chr, id, cm, bp, A1, A2, effect_allele, other_allele, effect_weight)

# flipped_swapped_dat need to have effect_allele  set to A1, other_allele set to A2, and effect_weight multiplied by -1
flipped_swapped_dat_to_merge <- flipped_swapped_dat %>%
	mutate(effect_allele2 = A1,
		   other_allele2 = A2,
		   effect_weight2 = effect_weight * -1) %>%
	select(chr, id, cm, bp, A1, A2, "effect_allele" = effect_allele2, "other_allele" = other_allele2, "effect_weight" = effect_weight2)


# merge all harmonized data

harmonized_dat <- matching_dat_to_merge %>%
	bind_rows(swapped_allele_dat_to_merge) %>%
	bind_rows(flipped_allele_dat_to_merge) %>%
	bind_rows(flipped_swapped_dat_to_merge) %>%
	arrange(bp)

missing_dat <- joined_info %>%
	anti_join(harmonized_dat, by = "id")

write_delim(x = harmonized_dat, file = output, delim = "\t", col_names = F)
write_delim(x = missing_dat, file = paste("missing_dat_from_model", chr, sep = "_"), delim = "\t", col_names = F)

