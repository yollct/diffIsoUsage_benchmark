library(iCOBRA)

path <- "/nfs/home/students/chit/is_benchmark"

cobra <- COBRAData_from_text(truth_file = paste0(path, "/results/truthtable_gene.csv"),
                    result_files = paste0(path, "/results/res_gene.txt"),
                    feature_id = "feature_id")

cobratx <- COBRAData_from_text(truth_file = paste0(path, "/results/truthtable_gene.csv"),
                    result_files = paste0(path, "/results/res_gene.txt"),
                    feature_id = "feature_id")

cobraperf <- calculate_performance(cobra,
binary_truth = "status", aspects = "fprtpr")
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2",
incltruth = TRUE)



COBRAapp(cobra)

