#' Chu (1988)' iterative soilution
#'
#' @param Fpi       Initial value for cumulative infiltration
#' @param Time      Time
#' @param tp        Ponding time
#' @param tpp       Equaivalent time for initially ponded condition
#' @param soil.para Soil hydraulic parameter
#' @param tol       Tolerance threshold for abslolute error in iterative solution
#'
#' @export
#'
#'
sschu<- function (Fpi, Time, tp, tpp, soil.para, tol = 0.00001, ...){

      M = soil.para$ths - soil.para$thi
      Sav = soil.para$Sav
      Ks = soil.para$Ksat

     Fpt = Ks*Time + M*Sav*log(1+ Fpi/(M*Sav))
       iter<- 0
       itermax <- 300
       error = abs(Fpt  -Fpi)

       while(error > tol & (iter< itermax)) {
         Fpi = Fpt
             Fpt = Ks*Time + M*Sav*log(1+ Fpi/(M*Sav))
      # print(c(iter = iter, res = (Fpt- Fpi)))
       iter<- iter+1
       error = abs(Fpi- Fpt)
      print(c(iter, error, Fpt))
      }
      if(error > tol & (iter == itermax)){
          stop("maximum iteration reached!")
      }
 return(Fpt)
}
####
# sschu<- function(Fpi, Time, tp, tpp, soil.para, tol = 0.00001, ...){
#
#       bm = soil.para$ths - soil.para$thi
#       Sav = soil.para$Sav
#       Ks = soil.para$Ksat
#
#       iter = 0
#       hh = Fpi - bm*Sav*log ( 1 + (Fpi/(bm*Sav))) - Ks*(Time - tp + tpp)
#
#       itermax  = 300
#       error = abs(hh)
#
#       while(error > tol & iter < itermax){
#             dhdf = 1.0 - ((bm*Sav) /(bm*Sav + Fpi))
#             Fpt = Fpi - hh/dhdf
#             Fpi = Fpt
#             hh = Fpi - bm*Sav*log ( 1 + (Fpi/(bm*Sav))) - Ks*(Time - tp + tpp)
#
#             error = abs(hh)
#
#             iter = iter + 1
#             print(c(iter, error, Fpt))
#       }
#       return(Fpt)
# }
