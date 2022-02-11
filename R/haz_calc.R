#' The hazard calculator
#'
#' This function calculates hazard for the given IMs at the user-defined return periods or IM levels at
#' the given sites using UCERF3 3.1 and 3.2 branches average source model.
#'
#' @param slat The latitude of the interest site
#' @param slon The longitude of the interest site
#' @param Vs30 The associated Vs30 value (in m/s) to the interest site
#' @param z1 The associated z1 value (in km) to the interest site. -999 for unknown.
#' @param max_dist The maximum cut-off distance for considered sources. Default is 300 km.
#' @param periods The list of interest IMs. Use -1 for PGV, 0 for PGA, and other NGA-West2 periods for PSA
#' @param IM_mat A matrix that specifies the IM exceedance levels for hazard curve calculation.
#' If \code{periods} length is 1, then \code{IM_mat} is a vector;
#' otherwise, \code{IM_mat} is a matrix with row number equals to \code{periods} length.
#' Each row in \code{IM_mat} is the specified IM exceedance levels for each of IM in \code{periods}.
#' The default is NA, which uses the default IM exceedance levels.
#' @param site_name The name of the site
#' @param output_dir The directory for the outputs
#' @examples Haz_calc(slat = 36.1, slon = -118.1, Vs30 = 350,
#' z1 = -999, periods = c(-1, 0.01, 5), output_dir = '~/Desktop/RPSHA')
#'
#' Haz_calc(slat = 36.5, slon = -119.5, Vs30 = 567,
#' z1 = 1.5, periods = c(0.01, 1, 5),
#' IM_mat = rbind(c(0.001, 0.0016), c(0.001, 0.0016), c(0.001, 0.0016)), site_name = 'Site2',
#' output_dir = '~/Desktop/RPSHA')
#'
#' @export
Haz_calc <- function(slat, slon, Vs30, z1 = -999, max_dist = 300, periods, IM_mat = NA, site_name = 'Site1', output_dir) {

  if (is.na(IM_mat[1])) {

    IMs_pga_psa <- c(0.0025, 0.0045, 0.0075, 0.0113, 0.0169, 0.0253, 0.038, 0.057, 0.0854,
                     0.128, 0.192, 0.288, 0.432, 0.649, 0.973, 1.46, 2.19, 3.28, 4.92, 7.38)

    IMs_pgv <- c(0.01, 0.0177, 0.0312, 0.0552, 0.0976, 0.173, 0.305, 0.539, 0.953, 1.68,
                 2.98, 5.26, 9.3, 16.4, 29.1, 51.3, 90.8, 160.0, 284.0, 501.0)

    IM_mat <- matrix(data = NA, nrow = length(periods), ncol = 20)

    for (i in 1:nrow(IM_mat)) {

      if (periods[i] == -1) {
        IM_mat[i, ] <- IMs_pgv
      } else {
        IM_mat[i, ] <- IMs_pga_psa
      }
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  haz_curves <- matrix(data = 0, nrow = nrow(IM_mat), ncol = ncol(IM_mat))

  # hazard calculation for Branch 3.1 Flts
  Flts <- flts_filter(slat = slat, slon = slon, branch = 1)

  for (i in 1:length(periods)) {

    tmp_Haz_mat <- Haz_calc_sub(df_source = Flts, source_type = 0, period = periods[i],
                                IM = IM_mat[i, ], z1 = z1, Vs30 = Vs30,
                                weight = 0.5)   # weight of branch 3.1 is 0.5

    haz_curves[i, ] <- colSums(tmp_Haz_mat)
  }

  # hazard calculation for Branch 3.2 Flts
  Flts <- flts_filter(slat = slat, slon = slon, branch = 2)

  for (i in 1:length(periods)) {

    tmp_Haz_mat <- Haz_calc_sub(df_source = Flts, source_type = 0, period = periods[i],
                                IM = IM_mat[i, ], z1 = z1, Vs30 = Vs30,
                                weight = 0.5)   # weight of branch 3.2 is 0.5

    haz_curves[i, ] <- colSums(tmp_Haz_mat) + haz_curves[i, ]
  }

  # hazard calculation for Branch 3.1 Pts
  Pts <- pts_filter(slat = slat, slon = slon, branch = 1)

  for (i in 1:length(periods)) {

    tmp_Haz_mat <- Haz_calc_sub(df_source = Pts, source_type = 1, period = periods[i],
                                IM = IM_mat[i, ], z1 = z1, Vs30 = Vs30,
                                weight = 0.5)   # weight of branch 3.1 is 0.5

    haz_curves[i, ] <- colSums(tmp_Haz_mat) + haz_curves[i, ]
  }

  # hazard calculation for Branch 3.2 Pts
  Pts <- pts_filter(slat = slat, slon = slon, branch = 2)

  for (i in 1:length(periods)) {

    tmp_Haz_mat <- Haz_calc_sub(df_source = Pts, source_type = 1, period = periods[i],
                                IM = IM_mat[i, ], z1 = z1, Vs30 = Vs30,
                                weight = 0.5)   # weight of branch 3.2 is 0.5

    haz_curves[i, ] <- colSums(tmp_Haz_mat) + haz_curves[i, ]
  }

  # loop for each IM in periods and write out hazard curves
  for (i_T in 1:length(periods)) {

    # create sub-folders for each IMs
    if (periods[i_T] == -1) {

      dir.create(paste0(output_dir, '/PGV/'))

      working_dir <- paste0(output_dir, '/PGV/')

      Tmp_res <- cbind(IM_mat[i_T,], haz_curves[i_T,])

      colnames(Tmp_res) <- c('PGV (cm/s)', 'Annual excandance rate')

    } else if (periods[i_T] == 0) {

      dir.create(paste0(output_dir, 'PGA/'))

      working_dir <- paste0(output_dir, '/PGA/')

      Tmp_res <- cbind(IM_mat[i_T,], haz_curves[i_T,])

      colnames(Tmp_res) <- c('PGA (g)', 'Annual excandance rate')

    } else {

      dir.create(paste0(output_dir, '/PSA', gsub('[.]', 'P', as.character(periods[i_T])), '/'))

      working_dir <- paste0(output_dir, '/PSA', gsub('[.]', 'P', as.character(periods[i_T])), '/')

      Tmp_res <- cbind(IM_mat[i_T,], haz_curves[i_T,])

      colnames(Tmp_res) <- c(paste0('PSA', gsub('[.]', 'P', as.character(periods[i_T])), ' (g)'), 'Annual excandance rate')
    }

    # write out hazard curves
    utils::write.csv(Tmp_res, file = paste0(working_dir, site_name, '_haz_curve.csv'), row.names = F)
  }
}



#' The subroutine of hazard curve calculation
#'
#' This is a subroutine of \code{haz_calc} function
#' @param df_source The returned dataframe by \code{flts_filter} and \code{pts_filter} functions
#' @param source_type The indicator of source type: 0 for flts; 1 for pts
#' @param period The interest period
#' @param IM IM level array
#' @param z1 Basin depth (km): depth from the ground surface to the 1km/s shear-wave horizon.
#' -999 if unknown. A numeric value.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s). A numeric value.
#' @param weight The weight that is supposed to be considered for the hazard
#' @return A hazard matrix at the given IM levels for each event in \code{df_sources}
#' @export
Haz_calc_sub <- function(df_source, source_type, period, IM, z1, Vs30, weight = 1) {

  if (source_type == 0) {

    Tmp_AR <- as.numeric(df_source$Rate)

    Tmp_M <- as.numeric(df_source$Mag)

    Tmp_Rjb <- df_source$Rjb

    Tmp_Frake <- as.numeric(df_source$Rake)

    Tmp_FT <- ifelse( Tmp_Frake > 30 & Tmp_Frake < 150, 3, # reverse
                      ifelse( Tmp_Frake > -150 & Tmp_Frake < -30, 2, # normal
                              1 ) ) # strike-slip

    Tmp_region <- rep(1, nrow(df_source))

    Tmp_GMM <- bssa_2014_nga(M = Tmp_M, period = period, Rjb = Tmp_Rjb, Fault_Type = Tmp_FT, region = Tmp_region,
                             z1 = z1, Vs30 = Vs30, coeffs = as.matrix(bssa_2014_coeffs))

    Pexd_mat <- apply(Tmp_GMM, 1, function(x){
      stats::pnorm( q = log( IM ), mean = log( x[1] ), sd = x[2], lower.tail = F )
    })

    Haz_mat <- t(Pexd_mat) * Tmp_AR * weight

  } else if (source_type == 1) {

    Tmp_AR <- as.numeric(df_source$Rate)

    Tmp_M <- as.numeric(df_source$Mag)

    Tmp_Rjb <- df_source$Rjb

    Tmp_FT <- df_source$Fault_Type

    Tmp_FT_weight <- df_source$Weight

    Tmp_region <- rep(1, nrow(df_source))

    Tmp_GMM <- bssa_2014_nga(M = Tmp_M, period = period, Rjb = Tmp_Rjb, Fault_Type = Tmp_FT, region = Tmp_region,
                             z1 = z1, Vs30 = Vs30, coeffs = as.matrix(bssa_2014_coeffs))

    Pexd_mat <- apply(Tmp_GMM, 1, function(x){
      stats::pnorm( q = log( IM ), mean = log( x[1] ), sd = x[2], lower.tail = F )
    })

    Haz_mat <- t(Pexd_mat) * Tmp_AR * Tmp_FT_weight * weight
  }

  return(Haz_mat)
}



#' The subroutine function for hazard curve calculation at a site from a given event
#'
#' This is a subroutine function to calculate the hazard curve (annual exceedance rate curve)
#' at a site from a given event
#'
#' @param Mag The magnitude of the event
#' @param AORate The annual occurrence rate of the event
#' @param Fault_type An indicator for fault type : 0 for unspecified fault;
#' 1 for strike-slip fault; 2 for normal fault; 3 for reverse fault.
#' @param region An indicator for region: 0 for global (incl. Taiwan); 1 for California;
#' 2 for Japan; 3 for China or Turkey; 4 for Italy.
#' @param periods An interest period or period arrays (0 for PGA, -1 for PGV)
#' @param IM_mat A vector or a matrix that specifies the IM exceedance levels for
#' hazard curve calculation.
#' If \code{periods} length is 1, then \code{IM_mat} is a vector;
#' otherwise, \code{IM_mat} is a matrix with row number equals to \code{periods} length.
#' Each row in \code{IM_mat} is the specified IM exceedance levels for each of IM in \code{periods}.
#' @param Rjb The Rjb distance from the event to the interest site
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s). A numeric value.
#' @param z1 Basin depth (km): depth from the ground surface to the 1km/s shear-wave horizon.
#' -999 if unknown. A numeric value.
#' @return A hazard curve matrix at the given IM levels for the event
#' @export
event_haz_calc <- function (Mag, AORate, Fault_type, region, periods, IM_mat, Rjb, Vs30, z1) {

  Haz_mat <- matrix(data = NA, nrow = nrow(IM_mat), ncol = ncol(IM_mat))

  for (i in 1:length(periods)) {

    Tmp_GMM <- bssa_2014_nga(M = Mag, period = periods[i], Rjb = Rjb, Fault_Type = Fault_type,
                             region = region, z1 = z1, Vs30 = Vs30,
                             coeffs = as.matrix(bssa_2014_coeffs))

    Pexd_mat <- stats::pnorm(q = log(IM_mat[i, ]), mean = log(Tmp_GMM$med),
                             sd = Tmp_GMM$sigma, lower.tail = F )

    Haz_mat[i, ] <- Pexd_mat * AORate

  }
  return(Haz_mat)
}




#' The function to construct the hazard matrix given event set and site list
#'
#' This is a function to calculate the hazard curves (annual exceedance rate curve)
#' at all sites from a event set
#'
#' @param Y_haz The target hazard. The column names and data format should be consistent with
#' the example on this page:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Haz.csv}
#' @param EventSet The event set. The column names and data format should be consistent with
#' the example on this page:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/EventSet.csv}
#' @param SiteTable The site list file set. The column names and data format should be consistent with
#' the example on this page:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/SiteTable.csv}
#' @return A dataframe with hazard curves produced by each event in \code{EventSet} is returned.
#' An example result can be found on this shared Google file:
#' \url{https://drive.google.com/file/d/1b_vxrXVSW4h0u_hzcrLif6CvQthj2ZO1/view?usp=sharing}
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom stats quantile
#' @export
events_hazmat_calc <- function (Y_haz, EventSet, SiteTable) {

  # extract IM types and levels in target hazard curves
  haz_IM_type_level <- Y_haz %>% select(IM_type, IM_level) %>% distinct() %>% arrange()
  IM_type <- unique(haz_IM_type_level$IM_type)
  tmp_len <- max(sapply(IM_type, function (x) sum(haz_IM_type_level$IM_type == x)))
  IM_mat <- matrix(data = NA, nrow = length(IM_type), ncol = tmp_len)
  for (i in 1:length(IM_type)) {
    tmp_level <- haz_IM_type_level$IM_level[haz_IM_type_level$IM_type == IM_type[i]]
    IM_mat[i, 1:length(tmp_level)] <- tmp_level
  }

  # develop X_haz matrix from each event in EventSet
  X <- matrix(data = 0, nrow = nrow(Y_haz), ncol = nrow(EventSet))
  periods <- ifelse(IM_type == 'PGA', 0,
                    ifelse(IM_type == 'PGV', -1, as.numeric(IM_type)))
  running_percentages <- round(quantile(seq(1, nrow(EventSet)), probs = seq(0, 1, 0.01)))
  for (i in 1:nrow(EventSet)) {
    Fault_type <- ifelse(EventSet$Fault_type[i] == 'R', 3,
                         ifelse(EventSet$Fault_type[i] == 'N', 2,
                                ifelse(EventSet$Fault_type[i] == 'SS', 1, 0)))
    Tmp_Rjb <- as.numeric(EventSet[i, -c(1:4)])
    Tmp_X <- matrix(data = NA, nrow = length(IM_mat) * nrow(SiteTable), ncol = 4)
    for (j in 1:nrow(SiteTable)) {
      tmp_haz <- event_haz_calc(Mag = EventSet$Magnitude[i], AORate = EventSet$AnnualOccurrenceRate[i],
                                Fault_type = Fault_type, region = 1, periods = periods,
                                IM_mat = IM_mat, Rjb = Tmp_Rjb[j],
                                Vs30 = SiteTable$Vs30[j], z1 = SiteTable$z1[j])
      Tmp_X[((j-1)*length(IM_mat) + 1):(j*length(IM_mat)), 1] <- SiteTable$SiteName[j]
      Tmp_X[((j-1)*length(IM_mat) + 1):(j*length(IM_mat)), 2] <- rep(IM_type, each = ncol(IM_mat))
      Tmp_X[((j-1)*length(IM_mat) + 1):(j*length(IM_mat)), 3] <- c(t(IM_mat))
      Tmp_X[((j-1)*length(IM_mat) + 1):(j*length(IM_mat)), 4] <- c(t(tmp_haz))
    }
    Tmp_X <- as.data.frame(Tmp_X)
    colnames(Tmp_X) <- c('SiteName', 'IM_type', 'IM_level', 'Haz_hat')
    Tmp_X$IM_level <- as.numeric(Tmp_X$IM_level)
    Tmp_X$Haz_hat <- as.numeric(Tmp_X$Haz_hat)
    Tmp_X <- left_join(Y_haz, Tmp_X, by = c('SiteName', 'IM_type', 'IM_level'))
    X[, i] <- Tmp_X$Haz_hat
    if (i %in% running_percentages) {
      print(paste0('The calculation is completed by ',
                   names(running_percentages)[which(running_percentages == i)]))
    }
  }

  X_haz <- Y_haz[c('IM_type', 'IM_level', 'SiteName')]
  X_haz <- cbind(X_haz, X)
  colnames(X_haz)[-c(1:3)] <- paste0('EventID_', EventSet$EventID)

  return(X_haz)
}



