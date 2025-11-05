CellCycleScoring_edited <- function (object, s.features, g2m.features, ctrl = NULL, set.ident = FALSE, 
                                     ...) 
{
  name <- "Cell.Cycle"
  features <- list(S.Score = s.features, G2M.Score = g2m.features)
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  object.cc <- AddModuleScore(object = object, features = features, 
                              name = name, ctrl = ctrl, ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), 
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "S", second = "G2M", null = "G1") {
    # Remover NA de scores
    valid_scores <- scores[!is.na(scores)]
    
    # Si todos los valores son NA, devuelve 'Undecided' o algún valor por defecto
    if (length(valid_scores) == 0) {
      return("Undecided")
    }
    
    if (all(valid_scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = valid_scores == max(valid_scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second)[which(x = valid_scores == 
                                        max(valid_scores))])
      }
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score", 
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Phase"
  }
  return(object)
}
