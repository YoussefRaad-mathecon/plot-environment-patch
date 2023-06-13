plot2 = function (x, animals = NULL, covs = NULL, ask = TRUE, breaks = "Sturges",
		  hist.ylim = NULL, sepAnimals = FALSE, sepStates = FALSE,
		  col = NULL, cumul = TRUE, plotTracks = TRUE, plotCI = FALSE,
		  alpha = 0.95, plotStationary = FALSE, ...)
{
	m <- x
	m <- delta_bc(m)
	nbAnimals <- length(unique(m$data$ID))
	stateNames <- m$stateNames
	nbStates <- length(stateNames)
	distnames <- names(m$conditions$dist)
	if (is.null(hist.ylim)) {
		hist.ylim <- vector("list", length(distnames))
		names(hist.ylim) <- distnames
	}
	for (i in distnames) {
		if (!is.null(hist.ylim[[i]]) & length(hist.ylim[[i]]) !=
			2)
			stop("hist.ylim$", i, " needs to be a vector of two values (ymin,ymax)")
	}
	if (!is.null(col) & length(col) >= nbStates)
		col <- col[1:nbStates]
	if (!is.null(col) & length(col) < nbStates) {
		warning("Length of 'col' should be at least number of states - argument ignored")
		if (nbStates < 8) {
			pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
					 "#0072B2", "#D55E00", "#CC79A7")
			col <- pal[1:nbStates]
		}
		else {
			hues <- seq(15, 375, length = nbStates + 1)
			col <- hcl(h = hues, l = 65, c = 100)[1:nbStates]
		}
	}
	if (is.null(col) & nbStates < 8) {
		pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
				 "#0072B2", "#D55E00", "#CC79A7")
		col <- pal[1:nbStates]
	}
	if (is.null(col) & nbStates >= 8) {
		hues <- seq(15, 375, length = nbStates + 1)
		col <- hcl(h = hues, l = 65, c = 100)[1:nbStates]
	}
	if (sepStates | nbStates < 2)
		cumul <- FALSE
	if (inherits(x, "miSum"))
		plotEllipse <- m$plotEllipse
	else plotEllipse <- FALSE
	muffWarn <- function(w) {
		if (any(grepl("zero-length arrow is of indeterminate angle and so skipped",
					  w)))
			invokeRestart("muffleWarning")
	}
	coordNames <- attr(m$data, "coords")
	if (!is.null(m$conditions$mvnCoords)) {
		coordNames <- c("x", "y")
		if (m$conditions$dist[[m$conditions$mvnCoords]] %in%
			c("mvnorm3", "rw_mvnorm3"))
			coordNames <- c("x", "y", "z")
		coordNames <- paste0(m$conditions$mvnCoords, ".", coordNames)
	}
	else if (is.null(coordNames))
		coordNames <- c("x", "y")
	if (is.null(animals))
		animalsInd <- 1:nbAnimals
	else {
		if (is.character(animals)) {
			animalsInd <- NULL
			for (zoo in 1:length(animals)) {
				if (length(which(unique(m$data$ID) == animals[zoo])) ==
					0)
					stop("Check 'animals' argument, ID not found")
				animalsInd <- c(animalsInd, which(unique(m$data$ID) ==
												  	animals[zoo]))
			}
		}
		if (is.numeric(animals)) {
			if (length(which(animals < 1)) > 0 | length(which(animals >
															  nbAnimals)) > 0)
				stop("Check 'animals' argument, index out of bounds")
			animalsInd <- animals
		}
	}
	nbAnimals <- length(animalsInd)
	ID <- unique(m$data$ID)[animalsInd]
	if (nbStates > 1) {
		if (inherits(x, "miSum"))
			states <- m$Par$states
		else {
			states <- viterbi(m)
			if (inherits(m, "hierarchical")) {
				hStates <- viterbi(m, hierarchical = TRUE)
			}
		}
	}
	else states <- rep(1, nrow(m$data))
	w <- iStates <- list()
	for (i in distnames) {
		if (!inherits(x, "hierarchical")) {
			if (sepStates | nbStates == 1)
				w[[i]] <- rep(1, nbStates)
			else {
				w[[i]] <- rep(NA, nbStates)
				for (state in 1:nbStates) w[[i]][state] <- length(which(states ==
																			state))/length(states)
			}
			iStates[[i]] <- 1:nbStates
			names(iStates[[i]]) <- stateNames
		}
		else {
			installDataTree()
			w[[i]] <- rep(0, nbStates)
			iLev <- gsub(paste0(".", i), "", names(m$conditions$hierDist$leaves)[grepl(i,
																					   names(m$conditions$hierDist$leaves))])
			iStates[[i]] <- m$conditions$hierStates$Get(function(x) data.tree::Aggregate(x,
																						 "state", min), filterFun = function(x) x$level ==
															(as.numeric(gsub("level", "", iLev)) + 1))
			if (sepStates) {
				w[[i]][iStates[[i]]] <- 1
			}
			else {
				denom <- length(hStates[[iLev]])
				for (state in 1:length(iStates[[i]])) {
					w[[i]][iStates[[i]][state]] <- length(which(hStates[[iLev]] ==
																	names(iStates[[i]][state])))/denom
				}
			}
		}
	}
	if (all(coordNames %in% names(m$data))) {
		x <- list()
		y <- list()
		z <- list()
		if (plotEllipse)
			errorEllipse <- list()
		for (zoo in 1:nbAnimals) {
			ind <- which(m$data$ID == ID[zoo])
			x[[zoo]] <- m$data[[coordNames[1]]][ind]
			y[[zoo]] <- m$data[[coordNames[2]]][ind]
			if (!is.null(m$conditions$mvnCoords)) {
				if (m$conditions$dist[[m$conditions$mvnCoords]] %in%
					c("mvnorm3", "rw_mvnorm3")) {
					z[[zoo]] <- m$data[[coordNames[3]]][ind]
					plotEllipse <- FALSE
					errorEllipse <- NULL
				}
			}
			if (plotEllipse)
				errorEllipse[[zoo]] <- m$errorEllipse[ind]
		}
	}
	covs <- getCovs(m, covs, ID)
	reForm <- formatRecharge(nbStates, m$conditions$formula,
							 m$conditions$betaRef, m$data, par = m$mle)
	recharge <- reForm$recharge
	hierRecharge <- reForm$hierRecharge
	newformula <- reForm$newformula
	nbCovs <- reForm$nbCovs
	aInd <- reForm$aInd
	nbG0covs <- reForm$nbG0covs
	nbRecovs <- reForm$nbRecovs
	g0covs <- reForm$g0covs
	recovs <- reForm$recovs
	if (!is.null(recharge)) {
		rechargeNames <- colnames(reForm$newdata)
		m$data[rechargeNames] <- reForm$newdata
		g0covs <- stats::model.matrix(recharge$g0, covs)
		g0 <- m$mle$g0 %*% t(g0covs)
		recovs <- stats::model.matrix(recharge$theta, covs)
		for (j in rechargeNames) {
			if (is.null(covs[[j]]))
				covs[[j]] <- mean(m$data[[j]])
		}
		covsCol <- cbind(get_all_vars(newformula, m$data), get_all_vars(recharge$theta,
																		m$data))
		if (!all(names(covsCol) %in% names(m$data))) {
			covsCol <- covsCol[, names(covsCol) %in% names(m$data),
							   drop = FALSE]
		}
		rawCovs <- covsCol[which(m$data$ID %in% ID), c(unique(colnames(covsCol))),
						   drop = FALSE]
	}
	else {
		rawCovs <- m$rawCovs[which(m$data$ID %in% ID), , drop = FALSE]
	}
	Par <- m$mle[distnames]
	ncmean <- get_ncmean(distnames, m$conditions$fullDM, m$conditions$circularAngleMean,
						 nbStates)
	nc <- ncmean$nc
	meanind <- ncmean$meanind
	tmPar <- lapply(Par, function(x) c(t(x)))
	parCount <- lapply(m$conditions$fullDM, ncol)
	for (i in distnames[!unlist(lapply(m$conditions$circularAngleMean,
									   isFALSE))]) {
		parCount[[i]] <- length(unique(gsub("cos", "", gsub("sin",
															"", colnames(m$conditions$fullDM[[i]])))))
	}
	parindex <- c(0, cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
	names(parindex) <- distnames
	for (i in distnames) {
		if (!is.null(m$conditions$DM[[i]])) {
			Par[[i]] <- m$mod$estimate[parindex[[i]] + 1:parCount[[i]]]
			if (!isFALSE(m$conditions$circularAngleMean[[i]])) {
				names(Par[[i]]) <- unique(gsub("cos", "", gsub("sin",
															   "", colnames(m$conditions$fullDM[[i]]))))
			}
			else names(Par[[i]]) <- colnames(m$conditions$fullDM[[i]])
		}
	}
	Par <- lapply(Par, function(x) c(t(x)))
	beta <- list(beta = m$mle$beta)
	if (!is.null(m$conditions$recharge)) {
		beta$g0 <- m$mle$g0
		beta$theta <- m$mle$theta
	}
	mixtures <- m$conditions$mixtures
	if (!m$conditions$stationary) {
		nbCovsDelta <- ncol(m$covsDelta) - 1
		foo <- length(m$mod$estimate) - ifelse(nbRecovs, (nbRecovs +
														  	1) + (nbG0covs + 1), 0) - (nbCovsDelta + 1) * (nbStates -
														  												   	1) * mixtures + 1
		delta <- m$mod$estimate[foo:(length(m$mod$estimate) -
									 	ifelse(nbRecovs, (nbRecovs + 1) + (nbG0covs + 1),
									 		   0))]
	}
	else {
		delta <- NULL
	}
	if (mixtures > 1) {
		if (!m$conditions$stationary)
			beta[["pi"]] <- m$mod$estimate[length(m$mod$estimate) -
										   	ncol(m$covsPi) * (mixtures - 1) - ifelse(nbRecovs,
										   											 (nbRecovs + 1) + (nbG0covs + 1), 0) - (nbCovsDelta +
										   											 									   	1) * (nbStates - 1) * mixtures + 1:(ncol(m$covsPi) *
										   											 									   											(mixtures - 1))]
		else beta[["pi"]] <- m$mod$estimate[length(m$mod$estimate) -
												ncol(m$covsPi) * (mixtures - 1) - ifelse(nbRecovs,
																						 (nbRecovs + 1) + (nbG0covs + 1), 0) + 1:(ncol(m$covsPi) *
																						 										 	(mixtures - 1))]
	}
	else beta[["pi"]] <- NULL
	tmpPar <- Par
	tmpConditions <- m$conditions
	for (i in distnames[which(m$conditions$dist %in% angledists)]) {
		if (!m$conditions$estAngleMean[[i]]) {
			tmpConditions$estAngleMean[[i]] <- TRUE
			tmpConditions$userBounds[[i]] <- rbind(matrix(rep(c(-pi,
																pi), nbStates), nbStates, 2, byrow = TRUE),
												   m$conditions$bounds[[i]])
			tmpConditions$workBounds[[i]] <- rbind(matrix(rep(c(-Inf,
																Inf), nbStates), nbStates, 2, byrow = TRUE),
												   m$conditions$workBounds[[i]])
			if (!is.null(m$conditions$DM[[i]])) {
				tmpPar[[i]] <- c(rep(0, nbStates), Par[[i]])
				if (is.list(m$conditions$DM[[i]])) {
					tmpConditions$DM[[i]]$mean <- ~1
				}
				else {
					tmpDM <- matrix(0, nrow(tmpConditions$DM[[i]]) +
										nbStates, ncol(tmpConditions$DM[[i]]) +
										nbStates)
					tmpDM[nbStates + 1:nrow(tmpConditions$DM[[i]]),
						  nbStates + 1:ncol(tmpConditions$DM[[i]])] <- tmpConditions$DM[[i]]
					diag(tmpDM)[1:nbStates] <- 1
					tmpConditions$DM[[i]] <- tmpDM
				}
			}
			else {
				Par[[i]] <- Par[[i]][-(1:nbStates)]
			}
		}
	}
	tmpInputs <- checkInputs(nbStates, tmpConditions$dist, tmpPar,
							 tmpConditions$estAngleMean, tmpConditions$circularAngleMean,
							 tmpConditions$zeroInflation, tmpConditions$oneInflation,
							 tmpConditions$DM, tmpConditions$userBounds, stateNames)
	tmpp <- tmpInputs$p
	splineInputs <- getSplineDM(distnames, tmpInputs$DM, m,
								covs)
	covs <- splineInputs$covs
	DMinputs <- getDM(covs, splineInputs$DM, tmpInputs$dist,
					  nbStates, tmpp$parNames, tmpp$bounds, tmpPar, tmpConditions$zeroInflation,
					  tmpConditions$oneInflation, tmpConditions$circularAngleMean)
	fullDM <- DMinputs$fullDM
	DMind <- DMinputs$DMind
	wpar <- n2w(tmpPar, tmpp$bounds, beta, delta, nbStates,
				tmpInputs$estAngleMean, tmpInputs$DM, tmpp$Bndind, tmpInputs$dist)
	ncmean <- get_ncmean(distnames, fullDM, tmpInputs$circularAngleMean,
						 nbStates)
	nc <- ncmean$nc
	meanind <- ncmean$meanind
	par <- w2n(wpar, tmpp$bounds, tmpp$parSize, nbStates, nbCovs,
			   tmpInputs$estAngleMean, tmpInputs$circularAngleMean,
			   tmpInputs$consensus, stationary = m$conditions$stationary,
			   fullDM, DMind, 1, tmpInputs$dist, tmpp$Bndind, nc, meanind,
			   m$covsDelta, tmpConditions$workBounds, m$covsPi)
	inputs <- checkInputs(nbStates, m$conditions$dist, Par,
						  m$conditions$estAngleMean, m$conditions$circularAngleMean,
						  m$conditions$zeroInflation, m$conditions$oneInflation,
						  m$conditions$DM, m$conditions$userBounds, stateNames)
	p <- inputs$p
	Fun <- lapply(inputs$dist, function(x) paste("d", x, sep = ""))
	for (i in names(Fun)) {
		if (Fun[[i]] == "dcat") {
			if (!requireNamespace("extraDistr", quietly = TRUE))
				stop("Package \"extraDistr\" needed for categorical distribution. Please install it.",
					 call. = FALSE)
			dcat <- extraDistr::dcat
		}
	}
	zeroMass <- oneMass <- vector("list", length(inputs$dist))
	names(zeroMass) <- names(oneMass) <- distnames
	legText <- stateNames
	tmpcovs <- covs
	for (i in which(mapply(is.numeric, covs))) {
		tmpcovs[i] <- round(covs[i], 2)
	}
	for (i in which(mapply(is.factor, covs))) {
		tmpcovs[i] <- as.character(covs[[i]])
	}
	if (inherits(m, "miSum")) {
		if (length(m$conditions$optInd)) {
			Sigma <- matrix(0, length(m$mod$estimate), length(m$mod$estimate))
			Sigma[(1:length(m$mod$estimate))[-m$conditions$optInd],
				  (1:length(m$mod$estimate))[-m$conditions$optInd]] <- m$MIcombine$variance
		}
		else {
			Sigma <- m$MIcombine$variance
		}
	}
	else if (!is.null(m$mod$hessian) && !inherits(m$mod$Sigma,
												  "error")) {
		Sigma <- m$mod$Sigma
	}
	else {
		Sigma <- NULL
		plotCI <- FALSE
	}
	par(mfrow = c(1, 1))
	par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1))
	par(ask = ask)
	if (!missing(...)) {
		arg <- list(...)
		if (any(!(names(arg) %in% plotArgs)))
			stop("additional graphical parameters are currently limited to: ",
				 paste0(plotArgs, collapse = ", "))
		if (!is.null(arg$cex.main))
			cex.main <- arg$cex.main
		else cex.main <- 1
		arg$cex.main <- NULL
		if (!is.null(arg$cex.legend))
			cex.legend <- arg$cex.legend
		else cex.legend <- 1
		arg$cex.legend <- NULL
		if (!is.null(arg[["cex"]]))
			cex <- arg[["cex"]]
		else cex <- 0.6
		arg$cex <- NULL
		if (!is.null(arg$asp))
			asp <- arg$asp
		else asp <- 1
		arg$asp <- NULL
		if (!is.null(arg$lwd))
			lwd <- arg$lwd
		else lwd <- 1.3
		arg$lwd <- NULL
		if (!is.null(arg$legend.pos)) {
			if (!(arg$legend.pos %in% c("bottomright", "bottom",
										"bottomleft", "left", "topleft", "top", "topright",
										"right", "center")))
				stop("legend.pos must be a single keyword from the list \"bottomright\", \"bottom\", \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\", \"right\" and \"center\"")
			legend.pos <- arg$legend.pos
		}
		else legend.pos <- NULL
		arg$legend.pos <- NULL
	}
	else {
		cex <- 0.6
		asp <- 1
		lwd <- 1.3
		cex.main <- 1
		cex.legend <- 1
		legend.pos <- NULL
		arg <- NULL
	}
	marg <- arg
	marg$cex <- NULL
	for (i in distnames) {
		if (m$conditions$dist[[i]] %in% mvndists) {
			if (m$conditions$dist[[i]] == "mvnorm2" || m$conditions$dist[[i]] ==
				"rw_mvnorm2") {
				tmpData <- c(m$data[[paste0(i, ".x")]], m$data[[paste0(i,
																	   ".y")]])
				if (m$conditions$dist[[i]] == "mvnorm2")
					ndim <- as.numeric(gsub("mvnorm", "", m$conditions$dist[[i]]))
				else ndim <- as.numeric(gsub("rw_mvnorm", "",
											 m$conditions$dist[[i]]))
			}
			else if (m$conditions$dist[[i]] == "mvnorm3" ||
					 m$conditions$dist[[i]] == "rw_mvnorm3") {
				tmpData <- c(m$data[[paste0(i, ".x")]], m$data[[paste0(i,
																	   ".y")]], m$data[[paste0(i, ".z")]])
				if (m$conditions$dist[[i]] == "mvnorm3")
					ndim <- as.numeric(gsub("mvnorm", "", m$conditions$dist[[i]]))
				else ndim <- as.numeric(gsub("rw_mvnorm", "",
											 m$conditions$dist[[i]]))
			}
		}
		else {
			tmpData <- m$data[[i]]
		}
		if (sepAnimals) {
			genData <- list()
			for (zoo in 1:nbAnimals) {
				ind <- which(m$data$ID == ID[zoo])
				if (m$conditions$dist[[i]] %in% mvndists) {
					if (m$conditions$dist[[i]] %in% c("mvnorm2",
													  "rw_mvnorm2"))
						genData[[zoo]] <- c(tmpData[ind], tmpData[(-(1:nrow(m$data)))][ind])
					else if (m$conditions$dist[[i]] %in% c("mvnorm3",
														   "rw_mvnorm3"))
						genData[[zoo]] <- c(tmpData[ind], tmpData[-(1:nrow(m$data))][ind],
											tmpData[-(1:(2 * nrow(m$data)))][ind])
				}
				else genData[[zoo]] <- tmpData[ind]
			}
		}
		else {
			ind <- which(m$data$ID %in% ID)
			if (m$conditions$dist[[i]] %in% mvndists) {
				if (m$conditions$dist[[i]] %in% c("mvnorm2",
												  "rw_mvnorm2"))
					genData <- tmpData[c(ind, ind + nrow(m$data))]
				else if (m$conditions$dist[[i]] %in% c("mvnorm3",
													   "rw_mvnorm3"))
					genData <- tmpData[c(ind, ind + nrow(m$data),
										 ind + 2 * nrow(m$data))]
			}
			else genData <- tmpData[ind]
		}
		zeroMass[[i]] <- rep(0, nbStates)
		oneMass[[i]] <- rep(0, nbStates)
		if (m$conditions$zeroInflation[[i]] | m$conditions$oneInflation[[i]]) {
			if (m$conditions$zeroInflation[[i]])
				zeroMass[[i]] <- par[[i]][nrow(par[[i]]) - nbStates *
										  	m$conditions$oneInflation[[i]] - (nbStates -
										  									  	1):0, ]
			if (m$conditions$oneInflation[[i]])
				oneMass[[i]] <- par[[i]][nrow(par[[i]]) - (nbStates -
														   	1):0, ]
			par[[i]] <- par[[i]][-(nrow(par[[i]]) - (nbStates *
													 	m$conditions$oneInflation[[i]] - nbStates *
													 	m$conditions$zeroInflation[[i]] - 1):0), , drop = FALSE]
		}
		infInd <- FALSE
		if (inputs$dist[[i]] %in% angledists)
			if (i == "angle" & ("step" %in% distnames))
				if (inputs$dist$step %in% stepdists & m$conditions$zeroInflation$step)
					if (all(coordNames %in% names(m$data)))
						infInd <- TRUE
		covNames <- getCovNames(m, p, i)
		DMterms <- covNames$DMterms
		DMparterms <- covNames$DMparterms
		if (inputs$consensus[[i]]) {
			for (jj in 1:nbStates) {
				if (!is.null(DMparterms$mean[[jj]]))
					DMparterms$kappa[[jj]] <- c(DMparterms$mean[[jj]],
												DMparterms$kappa[[jj]])
			}
		}
		factorterms <- names(m$data)[unlist(lapply(m$data, is.factor))]
		factorcovs <- paste0(rep(factorterms, times = unlist(lapply(m$data[factorterms],
																	nlevels))), unlist(lapply(m$data[factorterms], levels)))
		if (length(DMterms)) {
			for (jj in 1:length(DMterms)) {
				cov <- DMterms[jj]
				form <- stats::formula(paste("~", cov))
				varform <- all.vars(form)
				if (any(varform %in% factorcovs) && !all(varform %in%
														 factorterms)) {
					factorvar <- factorcovs %in% (varform[!(varform %in%
																factorterms)])
					DMterms[jj] <- rep(factorterms, times = unlist(lapply(m$data[factorterms],
																		  nlevels)))[which(factorvar)]
				}
			}
		}
		DMterms <- unique(DMterms)
		if (length(DMparterms)) {
			for (ii in 1:length(DMparterms)) {
				for (state in 1:nbStates) {
					if (length(DMparterms[[ii]][[state]])) {
						for (jj in 1:length(DMparterms[[ii]][[state]])) {
							cov <- DMparterms[[ii]][[state]][jj]
							form <- stats::formula(paste("~", cov))
							varform <- all.vars(form)
							if (any(varform %in% factorcovs) && !all(varform %in%
																	 factorterms)) {
								factorvar <- factorcovs %in% (varform[!(varform %in%
																			factorterms)])
								DMparterms[[ii]][[state]][jj] <- rep(factorterms,
																	 times = unlist(lapply(m$data[factorterms],
																	 					  nlevels)))[which(factorvar)]
							}
						}
						DMparterms[[ii]][[state]] <- unique(DMparterms[[ii]][[state]])
					}
				}
			}
		}
		covmess <- ifelse(!m$conditions$DMind[[i]], paste0(": ",
														   paste0(DMterms, " = ", tmpcovs[DMterms], collapse = ", ")),
						  "")
		genDensities <- list()
		genFun <- Fun[[i]]
		if (inputs$dist[[i]] %in% angledists) {
			grid <- seq(-pi, pi, length = 1000)
		}
		else {
			if (inputs$dist[[i]] %in% integerdists) {
				if (all(is.na(m$data[[i]])) || !is.finite(max(m$data[[i]],
															  na.rm = TRUE)))
					next
				if (inputs$dist[[i]] == "cat") {
					dimCat <- as.numeric(gsub("cat", "", m$conditions$dist[[i]]))
					grid <- seq(1, dimCat)
				}
				else grid <- seq(0, max(m$data[[i]], na.rm = TRUE))
			}
			else if (inputs$dist[[i]] %in% stepdists) {
				if (all(is.na(m$data[[i]])) || !is.finite(max(m$data[[i]],
															  na.rm = TRUE)))
					next
				grid <- seq(0, max(m$data[[i]], na.rm = TRUE),
							length = 10000)
			}
			else if (inputs$dist[[i]] %in% mvndists) {
				if (inputs$dist[[i]] == "mvnorm2" || inputs$dist[[i]] ==
					"rw_mvnorm2") {
					if (all(is.na(m$data[paste0(i, c(".x", ".y"))])) ||
						!is.finite(max(m$data[paste0(i, c(".x",
														  ".y"))], na.rm = TRUE)))
						next
					grid <- c(seq(min(m$data[[paste0(i, ".x")]],
									  na.rm = TRUE), max(m$data[[paste0(i, ".x")]],
									  				   na.rm = TRUE), length = 100), seq(min(m$data[[paste0(i,
									  				   													 ".y")]], na.rm = TRUE), max(m$data[[paste0(i,
									  				   													 										   ".y")]], na.rm = TRUE), length = 100))
				}
				else if (all(is.na(m$data[paste0(i, c(".x",
													  ".y", ".z"))])) || !is.finite(max(m$data[paste0(i,
													  												c(".x", ".y", ".z"))], na.rm = TRUE)))
					next
			}
			else {
				if (all(is.na(m$data[[i]])) || !is.finite(max(m$data[[i]],
															  na.rm = TRUE)))
					next
				grid <- seq(min(m$data[[i]], na.rm = TRUE),
							max(m$data[[i]], na.rm = TRUE), length = 10000)
			}
		}
		for (state in iStates[[i]]) {
			genArgs <- list(grid)
			if (m$conditions$dist[[i]] %in% mvndists) {
				genArgs[[2]] <- matrix(par[[i]][seq(state, nbStates *
														ndim, nbStates)], ndim, 1)
				sig <- matrix(0, ndim, ndim)
				lowertri <- par[[i]][nbStates * ndim + seq(state,
														   sum(lower.tri(matrix(0, ndim, ndim), diag = TRUE)) *
														   	nbStates, nbStates)]
				sig[lower.tri(sig, diag = TRUE)] <- lowertri
				sig <- t(sig)
				sig[lower.tri(sig, diag = TRUE)] <- lowertri
				genArgs[[3]] <- sig
			}
			else if (grepl("cat", m$conditions$dist[[i]])) {
				genArgs[[2]] <- t(par[[i]][seq(state, dimCat *
											   	nbStates, nbStates), ])
			}
			else {
				for (j in 1:(nrow(par[[i]])/nbStates)) genArgs[[j +
																	1]] <- par[[i]][(j - 1) * nbStates + state,
																	]
			}
			if (inputs$dist[[i]] == "gamma") {
				shape <- genArgs[[2]]^2/genArgs[[3]]^2
				scale <- genArgs[[3]]^2/genArgs[[2]]
				genArgs[[2]] <- shape
				genArgs[[3]] <- 1/scale
			}
			if (m$conditions$zeroInflation[[i]] | m$conditions$oneInflation[[i]]) {
				genDensities[[state]] <- cbind(grid, (1 - zeroMass[[i]][state] -
													  	oneMass[[i]][state]) * w[[i]][state] * do.call(genFun,
													  												   genArgs))
			}
			else if (infInd) {
				genDensities[[state]] <- cbind(grid, (1 - zeroMass$step[state]) *
											   	w[[i]][state] * do.call(genFun, genArgs))
			}
			else if (inputs$dist[[i]] %in% mvndists) {
				if (inputs$dist[[i]] == "mvnorm2" || inputs$dist[[i]] ==
					"rw_mvnorm2") {
					dens <- outer(genArgs[[1]][1:100], genArgs[[1]][101:200],
								  function(x, y) dmvnorm2(c(x, y), matrix(rep(genArgs[[2]],
								  											10000), 2), matrix(rep(genArgs[[3]], 10000),
								  															   2 * 2)))
					genDensities[[state]] <- list(x = genArgs[[1]][1:100],
												  y = genArgs[[1]][101:200], z = w[[i]][state] *
												  	dens)
				}
			}
			else {
				genDensities[[state]] <- cbind(grid, w[[i]][state] *
											   	do.call(genFun, genArgs))
			}
			for (j in p$parNames[[i]]) {
				for (jj in DMparterms[[j]][[state]]) {
					if (!is.factor(m$data[, jj])) {
						gridLength <- 101
						inf <- min(m$data[, jj], na.rm = T)
						sup <- max(m$data[, jj], na.rm = T)
						tempCovs <- data.frame(matrix(covs[jj][[1]],
													  nrow = gridLength, ncol = 1))
						if (length(DMterms) > 1)
							for (ii in DMterms[which(!(DMterms %in%
													   jj))]) tempCovs <- cbind(tempCovs, rep(covs[[ii]],
													   									   gridLength))
						names(tempCovs) <- c(jj, DMterms[which(!(DMterms %in%
																 	jj))])
						tempCovs[, jj] <- seq(inf, sup, length = gridLength)
					}
					else {
						gridLength <- nlevels(m$data[, jj])
						tempCovs <- data.frame(matrix(covs[jj][[1]],
													  nrow = gridLength, ncol = 1))
						if (length(DMterms) > 1)
							for (ii in DMterms[which(!(DMterms %in%
													   jj))]) tempCovs <- cbind(tempCovs, rep(covs[[ii]],
													   									   gridLength))
						names(tempCovs) <- c(jj, DMterms[which(!(DMterms %in%
																 	jj))])
						tempCovs[, jj] <- as.factor(levels(m$data[,
																  jj]))
					}
					for (ii in DMterms[which(unlist(lapply(m$data[DMterms],
														   is.factor)))]) tempCovs[[ii]] <- factor(tempCovs[[ii]],
														   										levels = levels(m$data[[ii]]))
					tmpSplineInputs <- getSplineDM(i, inputs$DM,
												   m, tempCovs)
					tempCovs <- tmpSplineInputs$covs
					DMinputs <- getDM(tempCovs, tmpSplineInputs$DM[i],
									  inputs$dist[i], nbStates, p$parNames[i],
									  p$bounds[i], Par[i], m$conditions$zeroInflation[i],
									  m$conditions$oneInflation[i], m$conditions$circularAngleMean[i])
					fullDM <- DMinputs$fullDM
					DMind <- DMinputs$DMind
					nc[[i]] <- apply(fullDM[[i]], 1:2, function(x) !all(unlist(x) ==
																			0))
					if (!isFALSE(inputs$circularAngleMean[[i]])) {
						meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,
																 , drop = FALSE], 1, function(x) !all(unlist(x) ==
																 									 	0))))
						if (length(meanind[[i]])) {
							angInd <- which(is.na(match(gsub("cos",
															 "", gsub("sin", "", colnames(nc[[i]]))),
														colnames(nc[[i]]), nomatch = NA)))
							sinInd <- colnames(nc[[i]])[which(grepl("sin",
																	colnames(nc[[i]])[angInd]))]
							nc[[i]][meanind[[i]], sinInd] <- ifelse(nc[[i]][meanind[[i]],
																			sinInd], nc[[i]][meanind[[i]], sinInd],
																	nc[[i]][meanind[[i]], gsub("sin", "cos",
																							   sinInd)])
							nc[[i]][meanind[[i]], gsub("sin", "cos",
													   sinInd)] <- ifelse(nc[[i]][meanind[[i]],
													   						   gsub("sin", "cos", sinInd)], nc[[i]][meanind[[i]],
													   						   									 gsub("sin", "cos", sinInd)], nc[[i]][meanind[[i]],
													   						   									 									 sinInd])
						}
					}
					gradfun <- function(wpar, k) {
						w2n(wpar, p$bounds[i], p$parSize[i], nbStates,
							nbCovs, inputs$estAngleMean[i], inputs$circularAngleMean[i],
							inputs$consensus[i], stationary = TRUE,
							fullDM, DMind, gridLength, inputs$dist[i],
							p$Bndind[i], nc[i], meanind[i], m$covsDelta,
							m$conditions$workBounds[c(i, "beta")],
							m$covsPi)[[i]][(which(tmpp$parNames[[i]] ==
												  	j) - 1) * nbStates + state, k]
					}
					est <- w2n(c(m$mod$estimate[parindex[[i]] +
													1:parCount[[i]]], beta$beta, beta[["pi"]]),
							   p$bounds[i], p$parSize[i], nbStates, nbCovs,
							   inputs$estAngleMean[i], inputs$circularAngleMean[i],
							   inputs$consensus[i], stationary = TRUE,
							   fullDM, DMind, gridLength, inputs$dist[i],
							   p$Bndind[i], nc[i], meanind[i], m$covsDelta,
							   m$conditions$workBounds[c(i, "beta")], m$covsPi)[[i]][(which(tmpp$parNames[[i]] ==
							   															 	j) - 1) * nbStates + state, ]
					if (plotCI) {
						dN <- t(mapply(function(x) tryCatch(numDeriv::grad(gradfun,
																		   c(m$mod$estimate[parindex[[i]] + 1:parCount[[i]]],
																		     beta$beta, beta[["pi"]]), k = x), error = function(e) NA),
									   1:gridLength))
						se <- t(apply(dN[, 1:parCount[[i]]], 1,
									  function(x) tryCatch(suppressWarnings(sqrt(x %*%
									  										   	Sigma[parindex[[i]] + 1:parCount[[i]],
									  										   		  parindex[[i]] + 1:parCount[[i]]] %*%
									  										   	x)), error = function(e) NA)))
						uci <- est + qnorm(1 - (1 - alpha)/2) *
							se
						lci <- est - qnorm(1 - (1 - alpha)/2) *
							se
						do.call(plot, c(list(tempCovs[, jj], est,
											 ylim = range(c(lci, est, uci), na.rm = TRUE),
											 xaxt = "n", xlab = jj, ylab = paste(i,
											 									ifelse(j == "kappa", "concentration",
											 										   j), "parameter"), main = paste0(names(iStates[[i]])[match(state,
											 										   														  iStates[[i]])], ifelse(length(tempCovs[,
											 										   														  									   DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]] ==
											 										   														  									   									jj)]]), paste0(": ", paste(DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]] ==
											 										   														  									   																							   	jj)], "=", tmpcovs[, DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]] ==
											 										   														  									   																							   														 	jj)]], collapse = ", ")), "")), type = "l",
											 lwd = lwd, cex.main = cex.main), arg))
						if (!all(is.na(se))) {
							ciInd <- which(!is.na(se))
							withCallingHandlers(do.call(arrows, c(list(as.numeric(tempCovs[ciInd,
																						   jj]), lci[ciInd], as.numeric(tempCovs[ciInd,
																						   									  jj]), uci[ciInd], length = 0.025, angle = 90,
																	   code = 3, col = gray(0.5), lwd = lwd),
																  arg)), warning = muffWarn)
						}
					}
					else do.call(plot, c(list(tempCovs[, jj],
											  est, xaxt = "n", xlab = jj, ylab = paste(i,
											  										 ifelse(j == "kappa", "concentration",
											  										 	   j), "parameter"), main = paste0(names(iStates[[i]])[match(state,
											  										 	   														  iStates[[i]])], ifelse(length(tempCovs[,
											  										 	   														  									   DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]] ==
											  										 	   														  									   									jj)]]), paste0(": ", paste(DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]] ==
											  										 	   														  									   																							   	jj)], "=", tmpcovs[, DMparterms[[j]][[state]][-which(DMparterms[[j]][[state]] ==
											  										 	   														  									   																							   														 	jj)]], collapse = ", ")), "")), type = "l",
											  lwd = lwd, cex.main = cex.main), arg))
					if (is.factor(tempCovs[, jj]))
						do.call(axis, c(list(1, at = tempCovs[,
															  jj], labels = tempCovs[, jj]), arg))
					else do.call(axis, c(list(1), arg))
				}
			}
		}
		if (!(inputs$dist[[i]] %in% mvndists)) {
			if (sepAnimals) {
				for (zoo in 1:nbAnimals) {
					if (sepStates) {
						for (state in iStates[[i]]) {
							gen <- genData[[zoo]][which(states[which(m$data$ID ==
																	 	ID[zoo])] == state)]
							message <- paste0("ID ", ID[zoo], " - ",
											  names(iStates[[i]])[match(state, iStates[[i]])],
											  covmess)
							plotHist(gen, genDensities, inputs$dist[i],
									 message, sepStates, breaks, state, hist.ylim[[i]],
									 col, names(iStates[[i]]), cumul = cumul,
									 iStates[[i]], ...)
						}
					}
					else {
						gen <- genData[[zoo]]
						message <- paste0("ID ", ID[zoo], covmess)
						plotHist(gen, genDensities, inputs$dist[i],
								 message, sepStates, breaks, NULL, hist.ylim[[i]],
								 col, names(iStates[[i]]), cumul = cumul,
								 iStates[[i]], ...)
					}
				}
			}
			else {
				if (sepStates) {
					for (state in iStates[[i]]) {
						gen <- genData[which(states == state)]
						if (nbAnimals > 1)
							message <- paste0("All animals - ", names(iStates[[i]])[match(state,
																						  iStates[[i]])], covmess)
						else message <- paste0("ID ", ID, " - ",
											   names(iStates[[i]])[match(state, iStates[[i]])],
											   covmess)
						plotHist(gen, genDensities, inputs$dist[i],
								 message, sepStates, breaks, state, hist.ylim[[i]],
								 col, names(iStates[[i]]), cumul = cumul,
								 iStates[[i]], ...)
					}
				}
				else {
					gen <- genData
					if (nbAnimals > 1)
						message <- paste0("All animals", covmess)
					else message <- paste0("ID ", ID, covmess)
					plotHist(gen, genDensities, inputs$dist[i],
							 message, sepStates, breaks, NULL, hist.ylim[[i]],
							 col, names(iStates[[i]]), cumul = cumul,
							 iStates[[i]], ...)
				}
			}
		}
		else if (inputs$dist[[i]] == "mvnorm2" || inputs$dist[[i]] ==
				 "rw_mvnorm2") {
			datNames <- paste0(i, c(".x", ".y"))
			if (sepAnimals) {
				for (zoo in 1:nbAnimals) {
					if (sepStates) {
						for (state in iStates[[i]]) {
							gen <- m$data[which(states == state &
													m$data$ID == ID[zoo]), datNames]
							message <- paste0("ID ", ID[zoo], " - ",
											  names(iStates[[i]])[match(state, iStates[[i]])],
											  covmess)
							plotHistMVN(gen, genDensities, inputs$dist[i],
										message, sepStates, breaks, state, col,
										names(iStates[[i]]), cumul = cumul,
										iStates[[i]], ...)
						}
					}
					else {
						gen <- m$data[which(m$data$ID == ID[zoo]),
									  datNames]
						message <- paste0("ID ", ID[zoo], covmess)
						plotHistMVN(gen, genDensities, inputs$dist[i],
									message, sepStates, breaks, NULL, col,
									names(iStates[[i]]), cumul = cumul, iStates[[i]],
									...)
					}
				}
			}
			else {
				if (sepStates) {
					for (state in iStates[[i]]) {
						gen <- m$data[which(states == state), datNames]
						if (nbAnimals > 1)
							message <- paste0("All animals - ", names(iStates[[i]])[match(state,
																						  iStates[[i]])], covmess)
						else message <- paste0("ID ", ID, " - ",
											   names(iStates[[i]])[match(state, iStates[[i]])],
											   covmess)
						plotHistMVN(gen, genDensities, inputs$dist[i],
									message, sepStates, breaks, state, col,
									names(iStates[[i]]), cumul = cumul, iStates[[i]],
									...)
					}
				}
				else {
					gen <- m$data[, datNames]
					if (nbAnimals > 1)
						message <- paste0("All animals", covmess)
					else message <- paste0("ID ", ID, covmess)
					plotHistMVN(gen, genDensities, inputs$dist[i],
								message, sepStates, breaks, NULL, col, names(iStates[[i]]),
								cumul = cumul, iStates[[i]], ...)
				}
			}
			par(mfrow = c(1, 1))
			par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1))
			par(ask = ask)
		}
		else if (inputs$dist[[i]] == "mvnorm3" || inputs$dist[[i]] ==
				 "rw_mvnorm3") {
			datNames <- paste0(i, c(".x", ".y", ".z"))
			tmpdist <- list()
			for (j in 1:length(datNames)) {
				tmpdist[[datNames[j]]] <- "norm"
				genArgs[[1]] <- seq(min(m$data[[datNames[j]]],
										na.rm = TRUE), max(m$data[[datNames[j]]],
														   na.rm = TRUE), length = 10000)
				genArgs[[2]] <- j
				for (state in iStates[[i]]) {
					genArgs[[3]] <- par[[i]][seq(state, nbStates *
												 	ndim, nbStates)]
					sig <- matrix(0, ndim, ndim)
					lowertri <- par[[i]][nbStates * ndim + seq(state,
															   sum(lower.tri(matrix(0, ndim, ndim), diag = TRUE)) *
															   	nbStates, nbStates)]
					sig[lower.tri(sig, diag = TRUE)] <- lowertri
					sig <- t(sig)
					sig[lower.tri(sig, diag = TRUE)] <- lowertri
					genArgs[[4]] <- sig
					genDensities[[state]] <- cbind(genArgs[[1]],
												   w[[i]][state] * do.call("dtmvnorm.marginal",
												   						genArgs))
				}
				if (sepAnimals) {
					for (zoo in 1:nbAnimals) {
						if (sepStates) {
							for (state in iStates[[i]]) {
								gen <- m$data[which(states == state &
														m$data$ID == ID[zoo]), datNames[j]]
								message <- paste0("ID ", ID[zoo], " - ",
												  names(iStates[[i]])[match(state, iStates[[i]])],
												  covmess)
								plotHist(gen, genDensities, tmpdist[datNames[j]],
										 message, sepStates, breaks, state,
										 hist.ylim[[i]], col, names(iStates[[i]]),
										 cumul = cumul, iStates[[i]], ...)
							}
						}
						else {
							gen <- m$data[which(m$data$ID == ID[zoo]),
										  datNames[j]]
							message <- paste0("ID ", ID[zoo], covmess)
							plotHist(gen, genDensities, tmpdist[datNames[j]],
									 message, sepStates, breaks, NULL, hist.ylim[[i]],
									 col, names(iStates[[i]]), cumul = cumul,
									 iStates[[i]], ...)
						}
					}
				}
				else {
					if (sepStates) {
						for (state in iStates[[i]]) {
							gen <- m$data[which(states == state),
										  datNames[j]]
							if (nbAnimals > 1)
								message <- paste0("All animals - ",
												  names(iStates[[i]])[match(state, iStates[[i]])],
												  covmess)
							else message <- paste0("ID ", ID, " - ",
												   names(iStates[[i]])[match(state, iStates[[i]])],
												   covmess)
							plotHist(gen, genDensities, tmpdist[datNames[j]],
									 message, sepStates, breaks, state, hist.ylim[[i]],
									 col, names(iStates[[i]]), cumul = cumul,
									 iStates[[i]], ...)
						}
					}
					else {
						gen <- m$data[, datNames[j]]
						if (nbAnimals > 1)
							message <- paste0("All animals", covmess)
						else message <- paste0("ID ", ID, covmess)
						plotHist(gen, genDensities, tmpdist[datNames[j]],
								 message, sepStates, breaks, NULL, hist.ylim[[i]],
								 col, names(iStates[[i]]), cumul = cumul,
								 iStates[[i]], ...)
					}
				}
			}
		}
	}
	par(mfrow = c(1, 1))
	par(mar = c(5, 4, 4, 2))
	par(ask = FALSE)
}