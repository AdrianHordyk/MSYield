

#' Import Parameters from xlsx file
#'
#' @param scenario Tab name or number
#' @param nspecies Number of species (rows) to use. Default is all
#' @param filename Name of the xlsx file
#'
#' @return
#' @export
#' @importFrom readxl read_excel

importData <- function(scenario=1,  nspecies=NA, filename='parameters.xlsx') {

  DF <- readxl::read_excel(filename, sheet=scenario)
  Names <- c("Species", "M_k", "Lm_Linf", "M", "K", "Linf", "L50", "L95", "SL50", "SL95",
             "B0", "Mpow", "h", "Weight")

  if (!any(Names %in% names(DF))) stop("invalid names in input file. Must be: ", paste(Names, " ", sep=" "))

  if (sum(DF$SL50 > DF$SL95)) {
    ind <- DF$SL50 > DF$SL95
    warning("SL95 < SL50 for some species: \n", paste(DF[ind,1], "\n"),
            "\nDropping species from analysis", call.=FALSE)
    DF <- DF[!ind,]
  }

  if (sum(DF$SL50 > DF$Linf)) {
    ind <- DF$SL50 > DF$Linf
    warning("SL95 > Linf for some species: \n", paste(DF[ind,1], "\n"),
            "\nDropping species from analysis", call.=FALSE)
    DF <- DF[!ind,]
  }

  if (sum(!apply(DF[,2:ncol(DF)], 2, is.finite))) {
    ind <- !apply(DF[,2:ncol(DF)], 2, is.finite)
    ind <- apply(ind, 1, sum)
    warning("Non-finite values found in: ", paste(DF[ind,1], "\n"),
            "\nDropping species from analysis", call.=FALSE)
    DF <- DF[!ind,]
  }

  if (!is.na(nspecies)) {
    return(DF[1:nspecies,])
  }
  DF

}



findR0 <- function(x, B0type=c("SSB0", "B0"), DF, Ftot, LMax, Lint) {
  exp(optimise(optR0, interval=log(c(0.01, 1E7)),
               B0type=B0type, DF=DF, Ftot=Ftot,
               LMax=LMax, Lint=Lint, x)$minimum)
}

optR0 <- function(logR0, B0type=c("SSB0", "B0"), DF, Ftot, LMax, Lint, x) {
  sp <- new("LB_pars", verbose=FALSE)
  sp@MK <- DF$M_k[x]
  sp@Linf <- DF$Linf[x]
  sp@L50 <- DF$L50[x]
  sp@L95 <- DF$L95[x]
  sp@SL50 <- DF$SL50[x]
  sp@SL95 <- DF$SL95[x]
  sp@Mpow <- DF$Mpow[x]
  sp@FM <- Ftot/DF$M[x]
  sp@BinMin <- 0
  sp@BinMax <- LMax
  sp@BinWidth <- Lint
  sp@R0 <- exp(logR0)
  Sim <- LBSPRsim(sp, Control=list(ngtg=31))
  b0 <- slot(Sim, B0type)
  (DF$B0[x] - b0)^2
}

findMLL <- function(x, DF, Ftot, LMax, Lint, fDisc, sdLegal=1, control=1) {
  if (control ==1 )
    return(exp(optimise(optMLL, interval=log(c(1, DF$Linf[x])),
               DF=DF,  Ftot=Ftot, LMax=LMax,
               Lint=Lint, fDisc=fDisc, sdLegal=1, control=control, x=x)$minimum))
  if (control == 2) {
    return(optMLL(logMLL=log(DF$mll[x]),  DF=DF,  Ftot=Ftot, LMax=LMax,
                        Lint=Lint, fDisc=fDisc, sdLegal=1, control=control, x=x))
  }
}

findGroupMLL <- function(group, DF, Ftot, LMax, Lint, fDisc, sdLegal=1, control=1) {
  tempDF <- DF %>% filter(groups==group)
  opt <- optimize(optgroup, interval=log(range(tempDF$optMLL)*c(0.7, 1.3)), tempDF, Ftot, LMax, Lint, fDisc, sdLegal=1, control=1)
  exp(opt$minimum)
}

optgroup <- function(logMLL, tempDF, Ftot, LMax, Lint, fDisc, sdLegal=1, control=1) {
  yield <- NULL
  for (x in 1:nrow(tempDF)) {
    yield[x] <- - optMLL(logMLL=logMLL, tempDF, Ftot, LMax, Lint, fDisc, sdLegal=1, control=1, x)
  }
  -sum(yield * tempDF$Weight)
}

optMLL <- function(logMLL, DF, Ftot, LMax, Lint, fDisc, sdLegal=1, control=1, x) {
  MLL <- exp(logMLL)

  sp <- new("LB_pars", verbose=FALSE)
  sp@Walpha <- 0.01
  sp@MK <- DF$M_k[x]
  sp@R0 <- DF$R0[x]
  sp@Linf <- DF$Linf[x]
  sp@L50 <- DF$L50[x]
  sp@L95 <- DF$L95[x]
  sp@SL50 <- DF$SL50[x]
  sp@SL95 <-  DF$SL95[x]
  sp@Mpow <- DF$Mpow[x]
  sp@FM <- Ftot/DF$M[x]
  sp@BinMin <- 0
  sp@BinMax <- LMax
  sp@BinWidth <- Lint
  sp@fDisc <- fDisc
  sp@MLL <- MLL
  sp@sdLegal <- sdLegal
  Sim <- LBSPRsim(sp, Control = list(ngtg=31))

  if (control==1) return(-Sim@Yield)
  if (control==2) return(cbind(Sim@Yield, Sim@SPR, Sim@SSB))
}


doCalcs <- function(DF2, nSp, FVec, LMax, Lint, fDisc, useParallel=FALSE) {
  saveout <- list()
  message("Calculating yield, SPR, and spawning biomass by F")
  for (XX in seq_along(FVec)) {
    # Yield <- -sapply(1:nSp, function(X)  optMLL(logMLL=log(DF2$mll[X]), DF=DF2, X,
    #                                             Ftot=FVec[XX], LMax=LMax, Lint=Lint,
    #                                             fDisc=fDisc, sdLegal=1))
    if (!useParallel) {
      tt <- sapply(1:nSp, findMLL, DF2, Ftot=FVec[XX], LMax=LMax, Lint=Lint,
                   fDisc=fDisc, sdLegal=1, control=2)
    } else {
      tt <- snowfall::sfSapply(1:nSp, findMLL, DF2, Ftot=FVec[XX], LMax=LMax, Lint=Lint,
                               fDisc=fDisc, sdLegal=1, control=2)
    }

    Yield <- tt[1,]
    SPR <- tt[2,]
    SB <- tt[3,]
    saveout[[XX]] <- data.frame(DF2, Yield=Yield, F=FVec[XX], SPR, SB)
  }
  outDF <- do.call("rbind", saveout)
  return(outDF)
}


calcR0 <- function(scenario=1, plot=TRUE, useParallel=TRUE, Lint=2, size=2.5) {
  suppressMessages(library(ggplot2))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggrepel))
  snowfall::sfInit(TRUE, 3)
  snowfall::sfLibrary(LBSPR)
  snowfall::sfLibrary(dplyr)

  message("Calculating R0 ...")
  DF <- importData(scenario, nspecies=NA)
  LMax <- ceiling(max(DF$Linf)/Lint) * Lint
  nSp <- nrow(DF)
  # --- Calc R0 that gives B0 ----
  if (!useParallel) {
    DF$R0 <-sapply(1:nSp, findR0, B0type="B0", DF, Ftot=0, LMax, Lint)
  } else {
    DF$R0 <-snowfall::sfSapply(1:nSp, findR0, B0type="B0", DF, Ftot=0, LMax, Lint)
  }

  p1 <-  ggplot2::ggplot(DF, ggplot2::aes(x=Linf, y=M_k, size=R0, label=Species)) + ggplot2::geom_point() +
    ggrepel::geom_text_repel(size=size) +ggplot2::theme_classic() + ggplot2::labs(x="Linf", y="M/K")

  p2 <-  ggplot2::ggplot(DF, ggplot2::aes(x=Linf, y=M_k, size=B0, label=Species)) + ggplot2::geom_point() +
    ggrepel::geom_text_repel(size=size) +ggplot2::theme_classic() + ggplot2::labs(x="Linf", y="M/K")

  pout <- gridExtra::grid.arrange(p1, p2, ncol=2)
  if (plot) pout

  df <- DF %>% dplyr::select(Species, M_k, Lm_Linf, R0, B0)
  if (!dir.exists('output')) dir.create('output')
  name <- paste0('output/R0_B0_scen_', scenario, ".csv")
  write.csv(df, file=name)

  return(invisible(DF))

}



runMod <- function(DF, F1, F2, fDisc, Lint=2, plot=TRUE, useParallel=TRUE, optMLL=TRUE) {
  LMax <- ceiling(max(DF$Linf)/Lint) * Lint
  nSp <- nrow(DF)

  # --- Optimise best MLL for each Species (at F = Ftot) ----
  message("Calculating best MLL for each species at F = ", F2)
  if (!useParallel) {
    optMLLs <- sapply(1:nSp, findMLL, DF, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=1)
  } else {
    optMLLs <- snowfall::sfSapply(1:nSp, findMLL, DF, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=1)
  }
  optMLLs <- round(optMLLs,0)
  optMLLs[optMLLs < DF$SL50] <- DF$SL50[optMLLs < DF$SL50]
  DF$optMLL <- optMLLs

  temp <- DF
  temp$mll <- temp$optMLL
  if (!useParallel) {
    temp2  <- sapply(1:nSp, findMLL, temp, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=2)
  } else{
    temp2  <- snowfall::sfSapply(1:nSp, findMLL, temp, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=2)
  }

  DF$optYield <- temp2[1,]
  DF$optSB <- temp2[3,]

  FVec <- c(F1, F2) # seq(0, to=HighF, by=0.05) # vector if fishing mortality

  NLimits <- 1:nSp
  out <- list()
  count <- 0
  doneNLimits <- NULL

  for (NLimit in NLimits) {
    message("Calculating ", NLimit, " Size Limit(s)")
    if (NLimit > nSp) NLimit <- nSp
    breaks <- unique(BAMMtools::getJenksBreaks(DF$optMLL, NLimit+1))
    breaks[which.min(breaks)] <- floor(breaks[which.min(breaks)]) - 2
    breaks[which.max(breaks)] <- ceiling(breaks[which.max(breaks)]) + 2
    DF$groups <- as.numeric(cut(DF$optMLL, breaks))
    actLimit <- length(unique(DF$groups))
    if (actLimit != NLimit)   message("Note: Some size limits are equal. Grouping into ", actLimit,  " instead")
    if (actLimit %in% doneNLimits) {
      message("Already calculated ", actLimit,  ". Skipping")
    } else {
      NLimit <- actLimit
      # DF2 <- DF %>% group_by(groups) %>% mutate(mll=weighted.mean(optMLL, optYield * (1 + 4*optMLL/100)))
      DF2 <- DF %>% group_by(groups) %>% mutate(weighted.mll=weighted.mean(optMLL, Weight))

      # DF2 %>% filter(groups == 1) %>% select(optMLL)

      # optimize for best MLL for group - weighted yield
      if (optMLL) {
        message("Optimizing for ", actLimit, ' size limits')
        if (!useParallel) {
          tt <- sapply(sort(unique(DF$groups)), findGroupMLL, DF2, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=1)
        } else {
          tt <- snowfall::sfSapply(sort(unique(DF$groups)), findGroupMLL, DF2, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=1)
        }

        tempDF <- data.frame(groups=sort(unique(DF$groups)), mll=tt)
        DF3 <- left_join(DF2, tempDF, by='groups')
        outDF <- doCalcs(DF3, nSp, FVec, LMax, Lint, fDisc, useParallel)
      } else {
        DF2$mll <- DF2$weighted.mll
        outDF <- doCalcs(DF2, nSp, FVec, LMax, Lint, fDisc, useParallel)
      }

      outDF$NLimits <- NLimit
      count <- count + 1
      out[[count]] <- outDF
      doneNLimits[count] <- NLimit
    }
  }

  outDF <- do.call("rbind", out)
  outDF$F1 <- F1
  outDF$F2 <- F2

  # Plots

  # Plot total yield by number of limits @ F = F1
  df <- outDF %>% filter(F==F1) %>% group_by(NLimits) %>% summarise(Yield=sum(Yield*Weight))
  p1 <- ggplot2::ggplot(df, aes(x = NLimits, y = Yield/max(Yield))) + geom_line() +
    geom_point() +ggplot2::ylim(0, 1)  + ggplot2::theme_classic() +
    ggplot2::labs(x="# Size Limits", y="Weighted Yield", title=paste("F = ", F1))

  # Plot total yield by number of limits @ F = F2
  df <- outDF %>% filter(F==F2) %>% group_by(NLimits) %>% summarise(Yield=sum(Yield*Weight))
  p2 <- ggplot2::ggplot(df, aes(x = NLimits, y = Yield/max(Yield))) + geom_line() + geom_point() + ggplot2::ylim(0, 1) + ggplot2::theme_classic() +
    ggplot2::labs(x="# Size Limits", y="Weighted Yield", title=paste("F = ", F2))

  # Plot number extinct by number of limits @ F = F1
  df <- outDF %>% filter(F==F1) %>% group_by(NLimits) %>% summarise(n=sum(SB==0))
  p5 <- ggplot2::ggplot(df, aes(x = NLimits, y = n/nSp)) + geom_line() + geom_point() + ggplot2::ylim(0, 1) +
    ggplot2::theme_classic() +
    ggplot2::labs(x="# Size Limits", y="Fraction Extinct", title=paste("F = ", F1))

  # Plot number extinct by number of limits @ F = F2
  df <- outDF %>% filter(F==F2) %>% group_by(NLimits) %>% summarise(n=sum(SB==0))
  p6 <- ggplot2::ggplot(df, aes(x = NLimits, y = n/nSp)) + geom_line() + geom_point() + ggplot2::ylim(0, 1) +
    ggplot2::theme_classic() +
    ggplot2::labs(x="# Size Limits", y="Fraction Extinct", title=paste("F = ", F2))

  pout <- gridExtra::grid.arrange(p1, p2, p5, p6, ncol=2)

  # output file
  Y1 <- outDF %>% filter(F==F1) %>% group_by(NLimits) %>% summarise(Rel.Yield.F1=sum(Yield*Weight),
                                                                    n.extinct.F1=sum(SB==0))
  Y1$Rel.Yield.F1 <- round(Y1$Rel.Yield.F1,0)
  Y2 <- outDF %>% filter(F==F2) %>% group_by(NLimits) %>% summarise(Rel.Yield.F2=sum(Yield*Weight),
                                                                    n.extinct.F2=sum(SB==0))
  Y2$Rel.Yield.F2 <- round(Y2$Rel.Yield.F2,0)

  F2_extinct <- outDF %>% filter(F==F2) %>% group_by(NLimits) %>% filter(SB==0) %>% select(NLimits, Species)
  tab <-  table(F2_extinct)
  ncol <- ncol(tab)
  nrow <- nrow(tab)
  mat <- matrix(NA, ncol=ncol, nrow=nrow)
  for (r in 1:nrow) {
    temp <-  names(tab[r,tab[r,] > 0])
    mat[r,] <- c(temp, rep('', ncol-length(temp)))
  }
  mat <- data.frame(mat, stringsAsFactors = FALSE)
  mat$NLimits <- unique(F2_extinct$NLimits)

  df <- left_join(Y1, Y2, by="NLimits")
  df2 <- left_join(df, mat, by="NLimits")
  df2[is.na(df2)] <- ''

  if (!dir.exists('output')) dir.create('output')
  name <- paste0('output/Yield_extinct_scen_', scenario, ".csv")
  write.csv(df2, file=name)


  return(outDF)
}




getSL <- function(nLimits=5, round=5, DF, fDisc, useParallel=TRUE, useWeight=TRUE, font.size=10) {

  tt <- DF %>% filter(NLimits==nLimits) %>% group_by(Species) %>% summarise(sl=unique(mll))

  sls <- sort(unique(tt$sl))
  n <- max(table(tt$sl))
  mat <- matrix(NA, nrow=n, ncol=length(sls))
  for (x in seq_along(sls)) {
    temp <- tt$Species[tt$sl == sls[x]]
    n2 <- n - length(temp)
    temp2 <- append(temp, rep(NA, n2))
    mat[,x] <- temp2

  }
  mat[is.na(mat)] <- ''
  mat <- as.data.frame(mat)


  F1 <- unique(DF$F1)

  # Plot % contribution by species @ F = F1
  if (!useWeight) {
    df <- DF %>% filter(NLimits==nLimits, F==F1) %>% group_by(Species) %>% summarise(YieldS=sum(Yield))
    df$RelYield <- df$YieldS/sum(df$YieldS)
    df <- df %>% mutate(Species=forcats::fct_reorder(Species, RelYield))

    p1 <- ggplot2::ggplot(df, aes(x=Species, y=RelYield))  + ggplot2::geom_bar(stat='identity') +
      ggplot2::theme_classic() +
      theme(text = element_text(size=font.size)) +
      ggplot2::labs(x="Species", y="Relative yield ", title=paste("F = ", F1)) + coord_flip()

  } else{
    df <- DF %>% filter(NLimits==nLimits, F==F1) %>% group_by(Species) %>% summarise(YieldS=sum(Yield*Weight))
    df$RelYield <- df$YieldS/sum(df$YieldS)
    df <- df %>% mutate(Species=forcats::fct_reorder(Species, RelYield))

    p1 <- ggplot2::ggplot(df, aes(x=Species, y=RelYield))  + ggplot2::geom_bar(stat='identity') +
      ggplot2::theme_classic() +
      theme(text = element_text(size=font.size)) +
      ggplot2::labs(x="Species", y="Weighted yield ", title=paste("F = ", F1, "  # Limits = ", nLimits)) +

      coord_flip()
  }

  print(p1)


  t2 <- DF %>% filter(F==F1, NLimits==nLimits)
  nSp <- nrow(t2)
  Lint <- 2
  LMax <- ceiling(max(t2$Linf)/Lint) * Lint
  t2$mll <- round(t2$mll/round,0) * round
  t2$Yield <- t2$SPR <- t2$F <- t2$SB <- NULL
  outDF <- doCalcs(t2, nSp, FVec=F1, LMax, Lint, fDisc, useParallel)

  chng <- (round(sum(outDF$Yield) / DF %>% filter(F==F1, nLimits==NLimits) %>% summarise(sum(Yield)),2)-1) * 100

  message("\nRounding to nearest ", round)
  message("\nExpected change in yield = ",  chng, "%")
  colnames(mat) <- round(sls/round,0) * round

  print(mat)

  if (!useWeight) {
    saveDF <- DF %>% filter(NLimits==nLimits) %>% group_by(Species, F) %>% summarise(YieldS=sum(Yield))
  } else {
    saveDF <- DF %>% filter(NLimits==nLimits) %>% group_by(Species, F) %>% summarise(YieldS=sum(Yield*Weight))
  }
  saveDF2 <- saveDF %>% select(Species, F, YieldS) %>% tidyr::spread(F, YieldS)

  if (!dir.exists('output')) dir.create('output')
  name <- paste0('output/NLimits_', nLimits, "_scen_",  scenario, ".csv")
  write.csv(saveDF2, file=name)

}





