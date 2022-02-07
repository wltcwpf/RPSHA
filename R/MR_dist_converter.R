#' A function to convert deaggregation to magnitude and/or distance hazard distribution
#'
#' This function takes the deaggregation results (could be multiple sites, multiple intensity measures,
#' and multiple return period or intensity measure level) and converts them to magntiude and/or
#' distance hazard distribution
#'
#' @param Deagg_df A dataframe contains the deaggregation results. It must contain
#' at least six or seven (if both magnitude and distance hazard distributions are interested)
#' columns: intensity measure type, intensity measure level, site name,
#' magnitude and/or distance bins for deaggregation, hazard contribution from each magnitude or
#' distance bin, and the hazard annual exceedance rate. Note for the intensity measure type,
#' it could be PGA, PGV or PSA at any oscillator periods. If PSA is interested, please
#' input the numeric oscillator periods (e.g., 0.01 and 0.5 for PSA at 0.01 sec and 0.5 sec,
#' respectively)
#' @param Deagg_level The variable (column) name in \code{Deagg_df} that is corresponding to
#' the variable at which the deaggreation was conducted.
#' Usually deaggregation is conducted based on return periods (e.g., 2475 years),
#' then \code{Deagg_level} should be set by the associated
#' column name in \code{Deagg_df}. For example, the default is "ReturnPeriod", which is
#' based on an example deaggregation result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' Sometimes the deaggregation is conducted based on different intensity measure levels
#' (e.g., PGA = 0.01 g), then \code{Deagg_level} should be set as "IM_level" instead based on the
#' example deaggregation result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' @param IM_type The variable (column) name in \code{Deagg_df} that is corresponding to intensity
#' measure types. The default is "IM_type" based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' @param IM_level The variable (column) name in \code{Deagg_df} that is corresponding to intensity
#' measure level. The default is "IM_level" based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' @param SiteName The variable (column) name in \code{Deagg_df} that is corresponding to
#' site IDs. The default is "SiteName" based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' @param Magnitude The variable (column) name in \code{Deagg_df} that is corresponding to magnitude.
#' The default is "Magnitude" based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' If \code{R_distance} is not \code{NULL}, then \code{Magnitude} could be \code{NULL} if
#' only magnitude hazard distribution is not interested.
#' @param R_distance The variable (column) name in \code{Deagg_df} that is corresponding to
#' distance. The default is "R_distance" (for distance hazard distribution)
#' based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' If \code{Magnitude} is not \code{NULL}, then \code{R_distance} could be \code{NULL} if
#' only distance hazard distribution is not interested.
#' @param HazdContrition The variable (column) name in \code{Deagg_df} that is corresponding to
#' hazard contribution. The default is "HazardContribution"
#' based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' @param AERate The variable (column) name in \code{Deagg_df} that is corresponding to
#' annual exceedance rate. The default is "AnnualExceedanceRate" based on the example result from this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/Y_Deagg.csv}.
#' @param distribution The distribution representation. "cdf" stands for cumulative hazard
#' and "pdf" stands for discrete bin hazard.
#'
#' @return A dataframe with six or seven (if both magnitude and distance hazard distributions are interested)
#' columns is returned, the columns inlcude IM_type, Deagg_level,
#' SiteName, IM_level, Hazd, Magnitude/R_distance (or MagRType and MagR)
#' @importFrom dplyr %>%
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom stats pnorm
#' @export
Deagg_to_MR_Hazd <- function (Deagg_df, Deagg_level = 'ReturnPeriod', IM_type = 'IM_type',
                              IM_level = 'IM_level',
                              SiteName = 'SiteName', Magnitude = 'Magnitude',
                              R_distance = 'R_distance',
                              HazdContrition = 'HazardContribution',
                              AERate = 'AnnualExceedanceRate', distribution = 'cdf') {

  if (is.null(Magnitude) & is.null(R_distance))
    stop('Magnitude and R_distance can not be both NULL!')

  # change colnames to the internally used names
  idx <- which(colnames(Deagg_df) == Deagg_level)
  ifelse (length(idx) != 1,
          stop(paste0('Error: Check the column name of ', Deagg_level, ' in your input dataframe!')),
          colnames(Deagg_df)[idx] <- 'Deagg_level')
  idx <- which(colnames(Deagg_df) == IM_type)
  ifelse (length(idx) != 1,
          stop(paste0('Error: Check the column name of ', IM_type, ' in your input dataframe!')),
          colnames(Deagg_df)[idx] <- 'IM_type')
  idx <- which(colnames(Deagg_df) == IM_level)
  ifelse (length(idx) != 1,
          stop(paste0('Error: Check the column name of ', IM_level, ' in your input dataframe!')),
          colnames(Deagg_df)[idx] <- 'IM_level')
  idx <- which(colnames(Deagg_df) == SiteName)
  ifelse (length(idx) != 1,
          stop(paste0('Error: Check the column name of ', SiteName, ' in your input dataframe!')),
          colnames(Deagg_df)[idx] <- 'SiteName')

  MR_flag <- 0
  if (!is.null(Magnitude)) {
    idx <- which(colnames(Deagg_df) == Magnitude)
    ifelse (length(idx) != 1,
            stop(paste0('Error: Check the column name of ', Magnitude, ' in your input dataframe!')),
            colnames(Deagg_df)[idx] <- 'Magnitude')
    MR_flag <- 1
  }

  if (!is.null(R_distance)) {
    idx <- which(colnames(Deagg_df) == R_distance)
    ifelse (length(idx) != 1,
            stop(paste0('Error: Check the column name of ', R_distance, ' in your input dataframe!')),
            colnames(Deagg_df)[idx] <- 'R_distance')
    MR_flag <- MR_flag + 2
  }

  idx <- which(colnames(Deagg_df) == HazdContrition)
  ifelse (length(idx) != 1,
          stop(paste0('Error: Check the column name of ', HazdContrition, ' in your input dataframe!')),
          colnames(Deagg_df)[idx] <- 'HazdContrition')
  idx <- which(colnames(Deagg_df) == AERate)
  ifelse (length(idx) != 1,
          stop(paste0('Error: Check the column name of ', AERate, ' in your input dataframe!')),
          colnames(Deagg_df)[idx] <- 'AERate')

  if (MR_flag == 1) {

    # only magnitude hazard
    Deagg_df_MR <- Deagg_df %>% group_by(IM_type, Deagg_level, SiteName, Magnitude, IM_level) %>%
      summarise(Hazd_sum = sum(HazdContrition), AERate = mean(AERate))

    if (distribution == 'cdf') {

      Deagg_df_val <- Deagg_df_MR %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = rev(cumsum(rev(Hazd_sum))) / sum(Hazd_sum) * AERate)

      Deagg_df_val$Magnitude <- Deagg_df_MR$Magnitude

      return(Deagg_df_val)

    } else if (distribution == 'pdf') {

      Deagg_df_val <- Deagg_df_MR %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = Hazd_sum / sum(Hazd_sum) * AERate)

      Deagg_df_val$Magnitude <- Deagg_df_MR$Magnitude

      return(Deagg_df_val)

    } else {
      stop('Error: Check your distribution input!')
    }

  } else if (MR_flag == 2) {

    # only distance hazard
    Deagg_df_MR <- Deagg_df %>% group_by(IM_type, Deagg_level, SiteName, R_distance, IM_level) %>%
      summarise(Hazd_sum = sum(HazdContrition), AERate = mean(AERate))

    if (distribution == 'cdf') {

      Deagg_df_val <- Deagg_df_MR %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = rev(cumsum(rev(Hazd_sum))) / sum(Hazd_sum) * AERate)

      Deagg_df_val$R_distance <- Deagg_df_MR$R_distance

      return(Deagg_df_val)

    } else if (distribution == 'pdf') {

      Deagg_df_val <- Deagg_df_MR %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = Hazd_sum / sum(Hazd_sum) * AERate)

      Deagg_df_val$R_distance <- Deagg_df_MR$R_distance

      return(Deagg_df_val)

    } else {
      stop('Error: Check your distribution input!')
    }

  } else if (MR_flag == 3) {

    # both magnitude and distance
    Deagg_df_M <- Deagg_df %>% group_by(IM_type, Deagg_level, SiteName, Magnitude, IM_level) %>%
      summarise(Hazd_sum = sum(HazdContrition), AERate = mean(AERate))

    Deagg_df_R <- Deagg_df %>% group_by(IM_type, Deagg_level, SiteName, R_distance, IM_level) %>%
      summarise(Hazd_sum = sum(HazdContrition), AERate = mean(AERate))

    if (distribution == 'cdf') {

      Deagg_df_val_M <- Deagg_df_M %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = rev(cumsum(rev(Hazd_sum))) / sum(Hazd_sum) * AERate)

      Deagg_df_val_M$Magnitude <- Deagg_df_M$Magnitude

      Deagg_df_val_R <- Deagg_df_R %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = rev(cumsum(rev(Hazd_sum))) / sum(Hazd_sum) * AERate)

      Deagg_df_val_R$R_distance <- Deagg_df_R$R_distance

      Deagg_df_val <- rbind(Deagg_df_val_M[, -6], Deagg_df_val_R[, -6])

      Deagg_df_val$MagRType <- c(rep('M', nrow(Deagg_df_val_M)),
                                 rep('R', nrow(Deagg_df_val_R)))

      Deagg_df_val$MagR <- c(Deagg_df_val_M$Magnitude, Deagg_df_val_R$R_distance)

      return(Deagg_df_val)

    } else if (distribution == 'pdf') {

      Deagg_df_val_M <- Deagg_df_M %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = Hazd_sum / sum(Hazd_sum) * AERate)

      Deagg_df_val_M$Magnitude <- Deagg_df_M$Magnitude

      Deagg_df_val_R <- Deagg_df_R %>% group_by(IM_type, Deagg_level, SiteName, IM_level) %>%
        summarise(Hazd = Hazd_sum / sum(Hazd_sum) * AERate)

      Deagg_df_val_R$R_distance <- Deagg_df_R$R_distance

      Deagg_df_val <- rbind(Deagg_df_val_M[, -6], Deagg_df_val_R[, -6])

      Deagg_df_val$MagRType <- c(rep('M', nrow(Deagg_df_val_M)),
                                 rep('R', nrow(Deagg_df_val_R)))

      Deagg_df_val$MagR <- c(Deagg_df_val_M$Magnitude, Deagg_df_val_R$R_distance)

      return(Deagg_df_val)

    } else {
      stop('Error: Check your distribution input!')
    }

  } else {
    stop('Error: no magnitude or distance hazards are calculated!')
  }
}


#' A function to convert an event to the target magnitude and/or distance hazard distribution
#'
#' This function takes a target magnitude and/or distance hazard distribution and an events to
#' calculate the hazard of the event at each of the corresponding magnitude/distance/
#' intensity measure type/intensity measure level or return period.
#'
#' @param MR_Hazd A dataframe with six or seven (if both magnitude and distance hazard distributions are interested)
#' columns: "IM_type", "Deagg_level", "SiteName", "IM_level", "Hazd", "Magnitude"/"R_distance"
#' (or "MagRType" and "MagR"). Note the column names should exactly match these variable names.
#' \code{MR_Hazd} can be directly provided by the function of \code{Deagg_to_MR_Hazd}.
#' @param Mag The magnitude of the event
#' @param AnnualOccurrenceRate The annual occurrence rate of the event
#' @param Fault_type The fault type of the event: 0 for unspecified fault;
#' 1 for strike-slip fault; 2 for normal fault; 3 for reverse fault
#' @param region The region of the event. Default is 1 for California
#' @param Sitelist A dataframe of sites with four columns: SiteName (should be consistent with
#' the \code{SiteName} in \code{MR_Hazd}), Vs30, z1, and Rjb. The example info for
#' SiteName, Vs30, and z1 can be found on this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/SiteTable.csv}.
#' The example of Rjb can be found on this webpage:
#' \url{https://raw.githubusercontent.com/wltcwpf/Dataset/main/RPSHA/EventSet.csv}.
#' @param distribution The distribution representation. "cdf" stands for cumulative hazard
#' and "pdf" stands for discrete bin hazard.
#'
#' @return A dataframe with six or seven (if both magnitude and distance hazard distributions are interested)
#' columns is returned, the columns inlcude IM_type, Deagg_level,
#' SiteName, IM_level, Hazd_event, Magnitude/R_distance (or MagRType and MagR)
#' @importFrom dplyr %>%
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom stats pnorm
#' @export
Event_to_MR_Hazd <- function (MR_Hazd, Mag, AnnualOccurrenceRate, Fault_type, region = 1,
                              Sitelist, distribution = 'cdf') {

  if (sum(c("IM_type", "Deagg_level", "SiteName", "IM_level") %in% colnames(MR_Hazd)) == 4) {

    if ("Magnitude" %in% colnames(MR_Hazd)) {

      # only magnitude hazard distribution
      Pexd <- rep(NA, nrow(MR_Hazd))
      MR_Hazd_M <- MR_Hazd %>% ungroup() %>% select(IM_type, SiteName) %>% distinct()

      for (i in 1:nrow(MR_Hazd_M)) {
        tmp_period <- MR_Hazd_M$IM_type[i]
        if (tmp_period == 'PGA') {
          tmp_period <- 0
        } else if (tmp_period == 'PGV') {
          tmp_period <- -1
        } else if (!is.na(as.numeric(as.character(tmp_period))) &
                   (as.numeric(as.character(tmp_period)) > 0)) {
          tmp_period <- as.numeric(as.character(tmp_period))
        } else {
          stop('Error: please check if the IM_type is PGA, PGV, or a numeric number
               for oscillator periods for PSA!')
        }
        tmp_Rjb <- Sitelist$Rjb[Sitelist$SiteName == MR_Hazd_M$SiteName[i]]
        tmp_z1 <- Sitelist$z1[Sitelist$SiteName == MR_Hazd_M$SiteName[i]]
        tmp_Vs30 <- Sitelist$Vs30[Sitelist$SiteName == MR_Hazd_M$SiteName[i]]
        tmp_GMM <- bssa_2014_nga(M = Mag, period = tmp_period,
                                 Rjb = tmp_Rjb, Fault_Type = Fault_type,
                                 region = region, z1 = tmp_z1,
                                 Vs30 = tmp_Vs30, coeffs = as.matrix(bssa_2014_coeffs))
        idx <- which((MR_Hazd$IM_type == MR_Hazd_M$IM_type[i]) &
                       (MR_Hazd$SiteName == MR_Hazd_M$SiteName[i]))
        Pexd[idx] <- stats::pnorm(q = log(MR_Hazd$IM_level[idx]),
                                  mean = log(tmp_GMM$med),
                                  sd = tmp_GMM$sigma,
                                  lower.tail = F)
      }

      MR_Hazd$Hazd_hat <- Pexd * AnnualOccurrenceRate
      MR_Hazd$Mag <- Mag

      if (distribution == 'cdf') {
        Tmp_MR_Hazd <- MR_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(Mag >= Magnitude, Hazd_hat, 0))
        # use left_join to ensure the same order as MR_Hazd
        MR_Hazd_Event <- left_join(MR_Hazd, Tmp_MR_Hazd, by = c('IM_type',
                                                                'Deagg_level',
                                                                'SiteName',
                                                                'IM_level',
                                                                'Magnitude'))
        return(MR_Hazd_Event[c('IM_type', 'Deagg_level', 'SiteName', 'IM_level',
                               'Hazd_event', 'Magnitude')])

      } else if (distribution == 'pdf') {
        Tmp_MR_Hazd <- MR_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(abs(Mag - Magnitude) == min(abs(Mag - Magnitude)),
                                   Hazd_hat, 0))
        # use left_join to ensure the same order as MR_Hazd
        MR_Hazd_Event <- left_join(MR_Hazd, Tmp_MR_Hazd, by = c('IM_type',
                                                                'Deagg_level',
                                                                'SiteName',
                                                                'IM_level',
                                                                'Magnitude'))
        return(MR_Hazd_Event[c('IM_type', 'Deagg_level', 'SiteName', 'IM_level',
                               'Hazd_event', 'Magnitude')])
      } else {
        stop('Error: Check your distribution input!')
      }

    } else if ("R_distance" %in% colnames(MR_Hazd)) {

      # only distance hazard distribution
      Pexd <- rep(NA, nrow(MR_Hazd))
      MR_Hazd_R <- MR_Hazd %>% ungroup() %>% select(IM_type, SiteName) %>% distinct()

      for (i in 1:nrow(MR_Hazd_R)) {
        tmp_period <- MR_Hazd_R$IM_type[i]
        if (tmp_period == 'PGA') {
          tmp_period <- 0
        } else if (tmp_period == 'PGV') {
          tmp_period <- -1
        } else if (!is.na(as.numeric(as.character(tmp_period))) &
                   (as.numeric(as.character(tmp_period)) > 0)) {
          tmp_period <- as.numeric(as.character(tmp_period))
        } else {
          stop('Error: please check if the IM_type is PGA, PGV, or a numeric number
               for oscillator periods for PSA!')
        }
        tmp_Rjb <- Sitelist$Rjb[Sitelist$SiteName == MR_Hazd_R$SiteName[i]]
        tmp_z1 <- Sitelist$z1[Sitelist$SiteName == MR_Hazd_R$SiteName[i]]
        tmp_Vs30 <- Sitelist$Vs30[Sitelist$SiteName == MR_Hazd_R$SiteName[i]]
        tmp_GMM <- bssa_2014_nga(M = Mag, period = tmp_period,
                                 Rjb = tmp_Rjb, Fault_Type = Fault_type,
                                 region = region, z1 = tmp_z1,
                                 Vs30 = tmp_Vs30, coeffs = as.matrix(bssa_2014_coeffs))
        idx <- which((MR_Hazd$IM_type == MR_Hazd_R$IM_type[i]) &
                       (MR_Hazd$SiteName == MR_Hazd_R$SiteName[i]))
        Pexd[idx] <- stats::pnorm(q = log(MR_Hazd$IM_level[idx]),
                                  mean = log(tmp_GMM$med),
                                  sd = tmp_GMM$sigma,
                                  lower.tail = F)
      }

      MR_Hazd$Hazd_hat <- Pexd * AnnualOccurrenceRate
      MR_Hazd <- left_join(MR_Hazd, Sitelist[c('SiteName', 'Rjb')], by = c('SiteName'))

      if (distribution == 'cdf') {
        Tmp_MR_Hazd <- MR_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(Rjb >= R_distance, Hazd_hat, 0))
        # use left_join to ensure the same order as MR_Hazd
        MR_Hazd_Event <- left_join(MR_Hazd, Tmp_MR_Hazd, by = c('IM_type',
                                                                'Deagg_level',
                                                                'SiteName',
                                                                'IM_level',
                                                                'R_distance'))
        return(MR_Hazd_Event[c('IM_type', 'Deagg_level', 'SiteName', 'IM_level',
                               'Hazd_event', 'R_distance')])

      } else if (distribution == 'pdf') {
        Tmp_MR_Hazd <- MR_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(abs(Rjb - R_distance) == min(abs(Rjb - R_distance)),
                                     Hazd_hat, 0))
        # use left_join to ensure the same order as MR_Hazd
        MR_Hazd_Event <- left_join(MR_Hazd, Tmp_MR_Hazd, by = c('IM_type',
                                                                'Deagg_level',
                                                                'SiteName',
                                                                'IM_level',
                                                                'R_distance'))
        return(MR_Hazd_Event[c('IM_type', 'Deagg_level', 'SiteName', 'IM_level',
                               'Hazd_event', 'R_distance')])
      } else {
        stop('Error: Check your distribution input!')
      }

    } else if (("MagRType" %in% colnames(MR_Hazd)) & ("MagR" %in% colnames(MR_Hazd))) {

      # both magnitude and distance hazard distributions

      # process magnitude hazard distributions
      M_Hazd <- MR_Hazd[MR_Hazd$MagRType == 'M', ]

      Pexd <- rep(NA, nrow(M_Hazd))
      MR_Hazd_M <- M_Hazd %>% ungroup() %>% select(IM_type, SiteName) %>% distinct()

      for (i in 1:nrow(MR_Hazd_M)) {
        tmp_period <- MR_Hazd_M$IM_type[i]
        if (tmp_period == 'PGA') {
          tmp_period <- 0
        } else if (tmp_period == 'PGV') {
          tmp_period <- -1
        } else if (!is.na(as.numeric(as.character(tmp_period))) &
                   (as.numeric(as.character(tmp_period)) > 0)) {
          tmp_period <- as.numeric(as.character(tmp_period))
        } else {
          stop('Error: please check if the IM_type is PGA, PGV, or a numeric number
               for oscillator periods for PSA!')
        }
        tmp_Rjb <- Sitelist$Rjb[Sitelist$SiteName == MR_Hazd_M$SiteName[i]]
        tmp_z1 <- Sitelist$z1[Sitelist$SiteName == MR_Hazd_M$SiteName[i]]
        tmp_Vs30 <- Sitelist$Vs30[Sitelist$SiteName == MR_Hazd_M$SiteName[i]]
        tmp_GMM <- bssa_2014_nga(M = Mag, period = tmp_period,
                                 Rjb = tmp_Rjb, Fault_Type = Fault_type,
                                 region = region, z1 = tmp_z1,
                                 Vs30 = tmp_Vs30, coeffs = as.matrix(bssa_2014_coeffs))
        idx <- which((M_Hazd$IM_type == MR_Hazd_M$IM_type[i]) &
                       (M_Hazd$SiteName == MR_Hazd_M$SiteName[i]))
        Pexd[idx] <- stats::pnorm(q = log(M_Hazd$IM_level[idx]),
                                  mean = log(tmp_GMM$med),
                                  sd = tmp_GMM$sigma,
                                  lower.tail = F)
      }

      M_Hazd$Hazd_hat <- Pexd * AnnualOccurrenceRate
      M_Hazd$Mag <- Mag

      if (distribution == 'cdf') {
        Tmp_M_Hazd <- M_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(Mag >= MagR, Hazd_hat, 0))
      } else if (distribution == 'pdf') {
        Tmp_M_Hazd <- M_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(abs(Mag - MagR) == min(abs(Mag - MagR)),
                                     Hazd_hat, 0))
      } else {
        stop('Error: Check your distribution input!')
      }


      # process distance hazard distributions
      R_Hazd <- MR_Hazd[MR_Hazd$MagRType == 'R', ]

      Pexd <- rep(NA, nrow(R_Hazd))
      MR_Hazd_R <- R_Hazd %>% ungroup() %>% select(IM_type, SiteName) %>% distinct()

      for (i in 1:nrow(MR_Hazd_R)) {
        tmp_period <- MR_Hazd_R$IM_type[i]
        if (tmp_period == 'PGA') {
          tmp_period <- 0
        } else if (tmp_period == 'PGV') {
          tmp_period <- -1
        } else if (!is.na(as.numeric(as.character(tmp_period))) &
                   (as.numeric(as.character(tmp_period)) > 0)) {
          tmp_period <- as.numeric(as.character(tmp_period))
        } else {
          stop('Error: please check if the IM_type is PGA, PGV, or a numeric number
               for oscillator periods for PSA!')
        }
        tmp_Rjb <- Sitelist$Rjb[Sitelist$SiteName == MR_Hazd_R$SiteName[i]]
        tmp_z1 <- Sitelist$z1[Sitelist$SiteName == MR_Hazd_R$SiteName[i]]
        tmp_Vs30 <- Sitelist$Vs30[Sitelist$SiteName == MR_Hazd_R$SiteName[i]]
        tmp_GMM <- bssa_2014_nga(M = Mag, period = tmp_period,
                                 Rjb = tmp_Rjb, Fault_Type = Fault_type,
                                 region = region, z1 = tmp_z1,
                                 Vs30 = tmp_Vs30, coeffs = as.matrix(bssa_2014_coeffs))
        idx <- which((R_Hazd$IM_type == MR_Hazd_R$IM_type[i]) &
                       (R_Hazd$SiteName == MR_Hazd_R$SiteName[i]))
        Pexd[idx] <- stats::pnorm(q = log(R_Hazd$IM_level[idx]),
                                  mean = log(tmp_GMM$med),
                                  sd = tmp_GMM$sigma,
                                  lower.tail = F)
      }

      R_Hazd$Hazd_hat <- Pexd * AnnualOccurrenceRate
      R_Hazd <- left_join(R_Hazd, Sitelist[c('SiteName', 'Rjb')], by = c('SiteName'))

      if (distribution == 'cdf') {
        Tmp_R_Hazd <- R_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(Rjb >= MagR, Hazd_hat, 0))

        # merge M and R and then use left_join to ensure the same order as MR_Hazd
        Tmp_MR_Hazd <- rbind(Tmp_M_Hazd, Tmp_R_Hazd)
        MR_Hazd_Event <- left_join(MR_Hazd, Tmp_MR_Hazd, by = c('IM_type',
                                                                'Deagg_level',
                                                                'SiteName',
                                                                'IM_level',
                                                                'MagRType',
                                                                'MagR'))
        return(MR_Hazd_Event[c('IM_type', 'Deagg_level', 'SiteName', 'IM_level',
                               'Hazd_event', 'MagRType', 'MagR')])

      } else if (distribution == 'pdf') {
        Tmp_R_Hazd <- R_Hazd %>% group_by(IM_type, SiteName, IM_level) %>%
          mutate(Hazd_event = ifelse(abs(Rjb - MagR) == min(abs(Rjb - MagR)),
                                     Hazd_hat, 0))

        # merge M and R and then use left_join to ensure the same order as MR_Hazd
        Tmp_MR_Hazd <- rbind(Tmp_M_Hazd, Tmp_R_Hazd)
        MR_Hazd_Event <- left_join(MR_Hazd, Tmp_MR_Hazd, by = c('IM_type',
                                                                'Deagg_level',
                                                                'SiteName',
                                                                'IM_level',
                                                                'MagRType',
                                                                'MagR'))
        return(MR_Hazd_Event[c('IM_type', 'Deagg_level', 'SiteName', 'IM_level',
                               'Hazd_event', 'MagRType', 'MagR')])
      } else {
        stop('Error: Check your distribution input!')
      }

    } else {
      stop('Error: The column names of MR_Hazd do not match the required names!')
    }
  } else {
    stop('Error: The column names of MR_Hazd do not match the required names!')
  }
}




