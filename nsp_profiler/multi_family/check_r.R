pkgs <- c("curatedMetagenomicData", "ExperimentHub", "SummarizedExperiment",
          "TreeSummarizedExperiment", "MicrobiomeBenchmarkData")
cat("R:", R.version.string, "\n")
for (p in pkgs) cat(sprintf("%-28s %s\n", p, requireNamespace(p, quietly = TRUE)))
cat("libs:\n")
cat(paste(.libPaths(), collapse = "\n"), "\n")
