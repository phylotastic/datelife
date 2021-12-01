test_that("lonely functions work", {
  relevant_age_quantiles(ages = threebirds_median$edge.length)
  numeric_vector_to_html_row(x = threebirds_median$edge.length)
  xx <- patristic_matrix_list_to_array(threebirds_result)
  patristic_matrix_sample(patristic_matrix_array = xx)
  patristic_matrix_subset(
    patristic_matrix =
      birds_yellowstone_result[[which.max(sapply(birds_yellowstone_result, nrow))]],
    taxa = rownames(birds_yellowstone_result[[which.max(sapply(birds_yellowstone_result, nrow))]])[1:10]
  )
  ri <- datelife_result_study_index(datelife_result = threebirds_result)
  datelife_authors_tabulate(results.index = ri)
  relevant_curators_tabulate(results.index = ri)
})
