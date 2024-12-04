library(quarto)

old_file_names <- list.files("delete_candidates/")
quarto_render(
  input = "delete_candidates/minimal_quarto_test.qmd",
  output_file = "minimal_quarto_test.html",
  execute_params = map(c("hub_01", "hub_02"), ~list(dataset_name=.))
)

new_file_names <- setdiff(list.files("scripts/"),old_file_names)
