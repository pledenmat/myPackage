#' chi_square_optim_DDM
#'
#' @param params Parameters of the model to be fitted + fixed parameters
#' @param observations Data to be fitted
#' @param returnFit if 0, then returns predicted data. Else, fit the parameters
#'
#' @return predicted data if returnFit == 0, else returns the prediction error
#'
#' @details This is the function used to simultaneously fit reaction time, accuracy and confidence judgments via quantile optimization.
#'
#' @export
chi_square_optim_DDM <- function(params, observations, returnFit){
  #First, generate predictions:
  drift <- params[9:length(params)]
  params <- params[1:8]
  names(params) <- c('a','ter','z','ntrials','sigma','dt','t2time','vratio')

  trial = data.frame(matrix(NA,nrow=0,ncol=7));names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in drift) {
    predictions <- data.frame(myPackage::DDM_with_confidence_slow(v=d,a=params['a'],ter=params['ter'],z=params['z'],ntrials=params['ntrials']/2,s=params['sigma'],dt=params['dt'],t2time=params['t2time'],postdriftmod=params['vratio']))
    predictionsneg <- data.frame(myPackage::DDM_with_confidence_slow(v=-d,a=params['a'],ter=params['ter'],z=params['z'],ntrials=params['ntrials']/2,s=params['sigma'],dt=params['dt'],t2time=params['t2time'],postdriftmod=params['vratio']))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- myPackage::fastmerge(predictions,predictionsneg)
    predictions['drift'] <- d
    trial <- myPackage::fastmerge(trial,predictions)
  }
  predictions <- trial



  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)

    #If we are fitting the model, now compare these predictions to the observations
  }else{
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]

    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]

    coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET

    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (d in 1:length(drift)) {
      # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
      c_quantiles <- quantile(c_observed[c_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      e_quantiles <- quantile(e_observed[e_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      if (any(is.na(e_quantiles))) {
        e_quantiles <- rep(0,5)
      }
      if (any(is.na(c_quantiles))) {
        c_quantiles <- rep(0,5)
      }
      # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
      prop_obs_c <- dim(c_observed[c_observed$coh == coherences[d],])[1] / dim(observations)[1]
      prop_obs_e <- dim(e_observed[e_observed$coh == coherences[d],])[1] / dim(observations)[1]

      c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
      e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
      obs_props <- c(obs_props,c_obs_proportion,e_obs_proportion)

      c_predicted_rt <- sort(c_predicted[c_predicted$drift == drift[d],]$rt)
      e_predicted_rt <- sort(e_predicted[e_predicted$drift == drift[d],]$rt)
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      c_pred_proportion <- c(
        sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[2]) - sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[3]) - sum(c_predicted_rt <= c_quantiles[2]),
        sum(c_predicted_rt <= c_quantiles[4]) - sum(c_predicted_rt <= c_quantiles[3]),
        sum(c_predicted_rt <= c_quantiles[5]) - sum(c_predicted_rt <= c_quantiles[4]),
        sum(c_predicted_rt > c_quantiles[5])
      ) / dim(predictions)[1]

      e_pred_proportion <- c(
        sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[2]) - sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[3]) - sum(e_predicted_rt <= e_quantiles[2]),
        sum(e_predicted_rt <= e_quantiles[4]) - sum(e_predicted_rt <= e_quantiles[3]),
        sum(e_predicted_rt <= e_quantiles[5]) - sum(e_predicted_rt <= e_quantiles[4]),
        sum(e_predicted_rt > e_quantiles[5])
      ) / dim(predictions)[1]
      pred_props <- c(pred_props,c_pred_proportion,e_pred_proportion)

      # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
      pred_props[pred_props==0] <- .0000001

    }


    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}

#' chi_square_optim_DDM_fullconfRT
#'
#' @param params Parameters of the model to be fitted + fixed parameters
#' @param observations Data to be fitted
#' @param returnFit if 0, then returns predicted data. Else, fit the parameters
#' @param confRT_name Name of the column in observations containing the confidence reaction time
#'
#' @return predicted data if returnFit == 0, else returns the prediction error
#'
#' @details This is the function used to simultaneously fit reaction time, accuracy and confidence judgments via quantile optimization.
#'
#' @export
chi_square_optim_DDM_fullconfRT <- function(params, observations, returnFit,confRT_name = "RTconf"){
  #First, generate predictions:
  drift <- params[8:length(params)]
  params <- params[1:7]
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio')

  coherences <- sort(unique(observations$coh)) #/!\ CHANGE ACCORDING TO DATASET

  trial = data.frame(matrix(NA,nrow=0,ncol=7))
  names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in 1:length(drift)) {
    predictions <- data.frame(myPackage::DDM_with_confidence_slow_fullconfRT(
      v=drift[d],a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift),s=params['sigma'],
      dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d],confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    predictionsneg <- data.frame(myPackage::DDM_with_confidence_slow_fullconfRT(
      v=-drift[d],a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift),s=params['sigma'],
      dt=params['dt'],t2distribution=rep(observations[observations$coh==coherences[d],confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- myPackage::fastmerge(predictions,predictionsneg)
    predictions['drift'] <- drift[d]
    trial <- myPackage::fastmerge(trial,predictions)
  }
  predictions <- trial



  #if we're only simulating data, return the predictions
  if(returnFit==0){
    return(predictions)

    #If we are fitting the model, now compare these predictions to the observations
  }else{
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]

    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]

    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (d in 1:length(drift)) {
      # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
      c_quantiles <- quantile(c_observed[c_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      e_quantiles <- quantile(e_observed[e_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      if (any(is.na(e_quantiles))) {
        e_quantiles <- rep(0,5)
      }
      if (any(is.na(c_quantiles))) {
        c_quantiles <- rep(0,5)
      }
      # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
      prop_obs_c <- dim(c_observed[c_observed$coh == coherences[d],])[1] / dim(observations)[1]
      prop_obs_e <- dim(e_observed[e_observed$coh == coherences[d],])[1] / dim(observations)[1]

      c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
      e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
      obs_props <- c(obs_props,c_obs_proportion,e_obs_proportion)

      c_predicted_rt <- sort(c_predicted[c_predicted$drift == drift[d],]$rt)
      e_predicted_rt <- sort(e_predicted[e_predicted$drift == drift[d],]$rt)
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      c_pred_proportion <- c(
        sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[2]) - sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[3]) - sum(c_predicted_rt <= c_quantiles[2]),
        sum(c_predicted_rt <= c_quantiles[4]) - sum(c_predicted_rt <= c_quantiles[3]),
        sum(c_predicted_rt <= c_quantiles[5]) - sum(c_predicted_rt <= c_quantiles[4]),
        sum(c_predicted_rt > c_quantiles[5])
      ) / dim(predictions)[1]

      e_pred_proportion <- c(
        sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[2]) - sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[3]) - sum(e_predicted_rt <= e_quantiles[2]),
        sum(e_predicted_rt <= e_quantiles[4]) - sum(e_predicted_rt <= e_quantiles[3]),
        sum(e_predicted_rt <= e_quantiles[5]) - sum(e_predicted_rt <= e_quantiles[4]),
        sum(e_predicted_rt > e_quantiles[5])
      ) / dim(predictions)[1]
      pred_props <- c(pred_props,c_pred_proportion,e_pred_proportion)

      # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
      pred_props[pred_props==0] <- .0000001

    }


    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    return(chiSquare)
  }
}
