#' Iterative solution to calculate end of ponding time
#'
#' @param tstart    Start of event
#' @param tend      End of event
#' @param tp        Ponding time
#' @param tpp       Equaivalent time for initially ponded condition
#' @param soil.para Soil hydraulic parameters
#' @param tol       Tolerance threshold for abslolute error in iterative solution
#'
#' @export
#'
#'
newtnp = function(tstart, tend, tp, tpp, rfi, amtinf, soil.para, tol = 0.00001, ...){

      Sav = soil.para$Sav
      Ks  = soil.para$Ksat
      M = soil.para$ths - soil.para$thi

  tnp = tend

  bftry = (tnp-tstart)*rfi + amtinf

  arglog = 1.0 + bftry/(M*Sav)

  hnp = bftry - M*Sav*log(arglog) - Ks*(tnp - tp + tpp)

  dhnp= rfi - rfi*(1.0/arglog) - Ks;

   tnp = tnp - hnp/dhnp

 error = abs(hnp)
   # tnp_old = tnp
  # error = abs(tnp_old - tnp)
  itermax = 300
  # error = 1
    # error = abs(tnp - tend)
 iter = 0
  while(error > tol & iter < itermax){

  iter = iter + 1
    # tnpold = tnp

    dhnp = rfi - rfi*(1/arglog) - Ks

    tnp = tnp - hnp/dhnp

   if(tnp < 0){
          tnp = 0
    }

    bftry = (tnp - tstart)*rfi + amtinf
    arglog = 1 + bftry/(M*Sav)

    hnp = bftry - M*Sav*log(arglog) - Ks*(tnp - tp + tpp)
    error = abs(hnp)


 # if(error < tol |iter == itermax) {
 #      break;
 #     tnpnew = tnp
 # }
 #

    print(c(iter, error, tnp))

    if(is.na(hnp)) {
          break;

    } else if(error < tol |iter == itermax) {
      break;

    }

  }

 tnpnew = tnp

 # else if(tnp < 0) {
 #          break;
 #          tnpnew = 0
 # }
 #


  return(tnpnew)


}

