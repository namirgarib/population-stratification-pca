#ifndef ANALYSIS_H
#define ANALYSIS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief perform_full_analysis:
 *        1) gather variant calls in a sparse manner
 *        2) run partial PCA
 *        3) write results to an output folder
 *
 * @param ref_file       path to reference genome
 * @param num_individuals
 * @param individuals_files array of individual genome paths
 */
void perform_full_analysis(const char *ref_file, int num_individuals,
                           const char **individuals_files);

#ifdef __cplusplus
}
#endif

#endif
