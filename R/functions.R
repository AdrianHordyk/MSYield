

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
  Names <- c("Species", "M_k", "Lm_Linf", "M", "K", "Linf", "L50", "L95", "SL50", "SL95","B0", "Mpow", "h")

  if (any(!names(DF) %in% Names)) stop("invalid names in input file. Must be: ", paste(Names, " ", sep=" "))

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
#
# M <- DF$M[x]
# K <- DF$K[x]
# Linf <- DF$Linf[x]
# L50 <- DF$L50[x]
# L95 <- DF$L95[x]
# SL50 <- DF$SL50[x]
# SL95 <- DF$SL95[x]
# Mpow <- DF$Mpow[x]
# Steepness <- DF$h[x]
#
# F <- 0.2
# CVLinf <- 0.1
# SDLinf <- CVLinf * Linf
#
# maxsd=2
# gtgby=1
# Wbeta=3
# Walpha <- 0.01
#
# MinL <- 0
# MaxL <- LMax
# BinWidth <- 1
# fDisc <- 0.2
# MLL <- 12
# R0 <- 1000
#
#
# GTGsim <- function(Linf, SDLinf, L50, L95, SL50, SL95, fDisc, MLL, Steepness,
#                    MinL, MaxL, BinWidth, Walpha,  Wbeta=3,
#                    maxsd=2, ngtg=13) {
#   gtgLinfs <- seq(from=Linf-maxsd*SDLinf, to=Linf+maxsd*SDLinf, length.out=ngtg)
#
#   recP <- dnorm(gtgLinfs, Linf, sd=SDLinf) / sum(dnorm(gtgLinfs, Linf, sd=SDLinf))
#
#   LBins <- seq(from=MinL, by=BinWidth, to=MaxL)
#   LMids <- seq(from=LBins[1] + 0.5*BinWidth, by=BinWidth, length.out=length(LBins)-1)
#
#   Weight <- Walpha * LMids^Wbeta
#
#   # Maturity and Fecundity for each GTG
#   L50GTG <- L50/Linf * gtgLinfs # Maturity at same relative size
#   L95GTG <- L95/Linf * gtgLinfs # Assumes maturity age-dependant
#   DeltaGTG <- L95GTG - L50GTG
#   MatLengtg <- sapply(seq_along(gtgLinfs), function (X)
#     1.0/(1+exp(-log(19)*(LMids-L50GTG[X])/DeltaGTG[X])))
#   FecLengtg <- MatLengtg * LMids^Wbeta # Fecundity across GTGs
#
#
#   # Selectivity - asymptotic only at this stage - by careful with knife-edge
#   SelLen <- 1.0/(1+exp(-log(19)*(LBins-(SL50+0.5*BinWidth))/ ((SL95+0.5*BinWidth)-(SL50+0.5*BinWidth))))
#
#   MK <- M/K
#   FM <- F/M
#   # Life-History Ratios
#   MKL <- MK * (Linf/(LBins+0.5*BinWidth))^Mpow # M/K ratio for each length class
#   # Matrix of MK for each GTG
#   MKMat <- matrix(rep(MKL, ngtg), nrow=length(MKL), byrow=FALSE)
#   FK <- FM * MK # F/K ratio
#   FKL <- FK * SelLen # F/K ratio for each length class
#
#   # Minimum legal length
#   plegal2 <- rep(1, length(LMids))
#   if (length(MLL)>0) {
#     plegal <- 1/(1+exp(-(LBins-MLL)/sdLegal))
#     plegal2 <- 1/(1+exp(-(LMids-MLL)/sdLegal))
#     FKL <- FKL * (plegal + (1-plegal) * fDisc)
#   }
#   ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG
#
#
#   # Set Up Empty Matrices
#   # number-per-recruit at length
#   NPRFished <- NPRUnfished <- matrix(0, nrow=length(LBins), ncol=ngtg)
#   NatLUF <- matrix(0, nrow=length(LMids), ncol=ngtg) # N at L unfished
#   NatLF <- matrix(0, nrow=length(LMids), ncol=ngtg) # N at L fished
#   FecGTG <- matrix(0, nrow=length(LMids), ncol=ngtg) # fecundity of GTG
#
#   # Distribute Recruits into first length class
#   NPRFished[1, ] <- NPRUnfished[1, ] <- recP * R0
#   for (L in 2:length(LBins)) { # Calc number at each size class
#     NPRUnfished[L, ] <- NPRUnfished[L-1, ] * ((gtgLinfs-LBins[L])/(gtgLinfs-LBins[L-1]))^MKMat[L-1, ]
#     NPRFished[L, ] <- NPRFished[L-1, ] *  ((gtgLinfs-LBins[L])/(gtgLinfs-LBins[L-1]))^ZKLMat[L-1, ]
#     ind <- gtgLinfs  < LBins[L]
#     NPRFished[L, ind] <- 0
#     NPRUnfished[L, ind] <- 0
#   }
#   NPRUnfished[is.nan(NPRUnfished)] <- 0
#   NPRFished[is.nan(NPRFished)] <- 0
#   NPRUnfished[NPRUnfished < 0] <- 0
#   NPRFished[NPRFished < 0] <- 0
#
#   for (L in 1:length(LMids)) { # integrate over time in each size class
#     NatLUF[L, ] <- (NPRUnfished[L,] - NPRUnfished[L+1,])/MKMat[L, ]
#     NatLF[L, ] <- (NPRFished[L,] - NPRFished[L+1,])/ZKLMat[L, ]
#     FecGTG[L, ] <- NatLUF[L, ] * FecLengtg[L, ]
#   }
#
#   SelLen2 <- 1.0/(1+exp(-log(19)*(LMids-SL50)/(SL95-SL50))) # Selectivity-at-Length
#   NatLV <- NatLUF * SelLen2 # Unfished Vul Pop
#   NatLC <- NatLF * SelLen2 # Catch Vul Pop
#
#
#   # Aggregate across GTGs
#   Nc <- apply(NatLC, 1, sum)/sum(apply(NatLC, 1, sum))
#   VulnUF <- apply(NatLV, 1, sum)/sum(apply(NatLV, 1, sum))
#   PopUF <- apply(NatLUF, 1, sum)/sum(apply(NatLUF, 1, sum))
#   PopF <- apply(NatLF, 1, sum)/sum(apply(NatLF, 1, sum))
#
#   # Calc SPR
#   B0 <- sum(NatLUF * Weight)
#   SSB0 <- sum(NatLUF * Weight * MatLengtg)
#   EPR0 <- sum(NatLUF * FecLengtg) # Eggs-per-recruit Unfished
#   EPRf <- sum(NatLF * FecLengtg) # Eggs-per-recruit Fished
#   SPR <- EPRf/EPR0
#
#   # Equilibrium Relative Recruitment
#   recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio
#   reca <- recK/EPR0
#   recb <- (reca * EPR0 - 1)/(R0*EPR0)
#   RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
#   if (!is.finite(RelRec)) RelRec <- 0
#
#   # RelRec/R0 - relative recruitment
#
#
#
#   YPR <- sum(NatLC  * Weight * SelLen2 * plegal2) * FM
#
#
#   Yield <- YPR * RelRec
#
#   # Spawning Stock Biomass
#   SSB <- sum(NatLF  * RelRec * Weight * MatLengtg)
#
#
# }

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


calcR0 <- function(scenario=1, plot=TRUE, useParallel=TRUE, Lint=2) {
  suppressMessages(library(ggplot2))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggrepel))
  snowfall::sfInit(TRUE, 3)
  snowfall::sfLibrary(LBSPR)

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
    ggrepel::geom_text_repel(size=4) +ggplot2::theme_classic() + ggplot2::labs(x="Linf", y="M/K")

  p2 <-  ggplot2::ggplot(DF, ggplot2::aes(x=Linf, y=M_k, size=B0, label=Species)) + ggplot2::geom_point() +
    ggrepel::geom_text_repel(size=4) +ggplot2::theme_classic() + ggplot2::labs(x="Linf", y="M/K")

  pout <- gridExtra::grid.arrange(p1, p2, ncol=2)
  if (plot) pout

  return(invisible(DF))

}



runMod <- function(DF, F1, F2, fDisc, Lint=2, plot=TRUE, useParallel=TRUE) {
  LMax <- ceiling(max(DF$Linf)/Lint) * Lint
  nSp <- nrow(DF)

  # --- Optimise best MLL for each Species (at F = Ftot) ----
  message("Calculating best MLL for each species at F = ", F2)
  if (!useParallel) {
    mll <- sapply(1:nSp, findMLL, DF, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=1)
  } else {
    mll <- snowfall::sfSapply(1:nSp, findMLL, DF, Ftot=F2, LMax, Lint, fDisc, sdLegal=1, control=1)
  }

  optMLLs <- round(mll,0)
  optMLLs[optMLLs < DF$SL50] <- DF$SL50[optMLLs < DF$SL50]
  DF$optMLL <- optMLLs


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
      DF2 <- DF %>% group_by(groups) %>% mutate(mll=mean(optMLL))
      DF2$mll <- round(DF2$mll,0)
      outDF <- doCalcs(DF2, nSp, FVec, LMax, Lint, fDisc, useParallel)
      outDF$NLimits <- NLimit
      count <- count + 1
      out[[count]] <- outDF
      doneNLimits[count] <- NLimit
    }
  }

  outDF <- do.call("rbind", out)

  # Plots

  # Plot total yield by number of limits @ F = F1
  df <- outDF %>% filter(F==F1) %>% group_by(NLimits) %>% summarise(Yield=sum(Yield))
  p1 <- ggplot2::ggplot(df, aes(x = NLimits, y = Yield/max(Yield))) + ggplot2::geom_line(size=2) + ggplot2::theme_classic() +
    ggplot2::labs(x="# Size Limits", y="Relative Yield", title=paste("F = ", F1))


  # Plot number extinct by number of limits @ F = F2
  df <- outDF %>% filter(F==F2) %>% group_by(NLimits) %>% summarise(n=sum(SB==0))
  p2 <- ggplot2::ggplot(df, aes(x = NLimits, y = n/nSp)) + ggplot2::geom_line(size=2) + ggplot2::ylim(0, 1) +
    ggplot2::theme_classic() +
    ggplot2::labs(x="# Size Limits", y="Fraction Extinct", title=paste("F = ", F2))


  # Plot % contribution by species @ F = F1
  df <- outDF %>% filter(F==F1) %>% group_by(Species)  %>% summarise(YieldS=sum(Yield))
  df$RelYield <- df$YieldS/sum(df$YieldS)

  df <- df %>% arrange(RelYield)
  df$RelYield2 <- cumsum(df$RelYield)
  df$n <- 1:nrow(df)

  p3 <- ggplot2::ggplot(df, aes(x=n, y=RelYield2)) + ggplot2::geom_line(size=2) + ggplot2::ylim(0, 1) +
    ggplot2::theme_classic() +
    ggplot2::labs(x="# Species", y="Cumulative contribution to total yield ", title=paste("F = ", F1))

  pout <- gridExtra::grid.arrange(p1, p2, p3, ncol=2)

  tt <- df %>% filter(RelYield > 0.05)
  if(dim(tt)[1] > 0) {
    message("\nSpecies contributing more than 5% of total catch at F = ", F1)
    df$RelYield <- round(df$RelYield,2)
    print(df %>% filter(RelYield > 0.05) %>% select(Species, RelYield) %>% as.data.frame())
  }

  return(outDF)
}


getSL <- function(DF, nLimits=5) {

  tt <- DF %>% filter(nLimits==NLimits) %>% group_by(Species) %>% summarise(sl=unique(mll))

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
  colnames(mat) <- sls
  mat
}


