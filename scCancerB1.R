library(scCancer)
dataPath <- "/my_path/B1"     # The path of cell ranger processed data
savePath <- "/my_path/scCancerQC/B1/scStatistics"  # A path to save the results
sampleName <- "B1"          # The sample name
authorName <- "Lyu@SCU"           # The author name to mark the report

# Run scStatistics
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName
)

# Run scAnnotation
savePath <- "/my_path/scCancerQC/B1/scAnnotation"  # A path to save the results
statPath <- "/my_path/scCancerQC/B1/scStatistics"
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  nCell.min = 5,
  bool.runDiffExpr = F,
  bool.runGeneSets = F,
  bool.runExprProgram = F,
  bool.runInteraction = F
)