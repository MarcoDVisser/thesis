// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List runModelcpp(IntegerVector R,
			    int N,
			    int P,
			    int D,
			    double E,
			    double S,
			    NumericVector uprob = NumericVector::create()
			    ) {

  // patches and sp abundance counts
  IntegerVector patches = RcppArmadillo::sample(R, P, TRUE, uprob) ;
  IntegerVector oldpatches = clone(patches);
  IntegerVector patchid = seq(0, P-1);
  IntegerVector spcounts(R.length());

  // emigration adjusted counts and natal establishment prob
  NumericVector ecounts(R.length());
  NumericVector natalspprob(R.length());
      
  for(int i = 0; i < N; i++) {
    // those doomed to die, and those deemed the succesors
    IntegerVector Deaths = RcppArmadillo::sample(patchid, D, FALSE, uprob);
    IntegerVector Recruits(D);
    IntegerVector Rec(1);

    // fill with zeros
    spcounts.fill(0);
    natalspprob.fill(0);
    Recruits.fill(0);
    ecounts.fill(0);
    
    // calculate the previous status quo    
    for (int j = 0; j < patches.length(); j++) {
      spcounts[patches[j]] += 1;
    }

   
     // for each death d institute succession 
     for (int d = 0; d < Deaths.length(); d++) {

       // Emmigration adjusted sp counts
       for (int r = 0; r < R.length(); r++) {

        // natal site survival
	 ecounts[r] = (E*(double)spcounts[r])/P;

	 if(patches[Deaths[d]] == r){
	 ecounts[r] = S*(1-E)+(S*ecounts[r]);
        } 
	 
       // caluclate  spprob for each species for natal site d
       	 natalspprob[r] = ecounts[r]/sum(ecounts);
	 
	 if(isnan(natalspprob[r])) {
	   natalspprob[r] = 0;
	 }

     }

     Rec = RcppArmadillo::sample(R, 1, FALSE, natalspprob);

     Recruits[d] = mean(Rec);
	 
     }
     
     patches[Deaths] = Recruits;
   }

  //Rcpp::Named("Start Patches") = oldpatche,
return Rcpp::List::create(Rcpp::Named("Patches") = patches,
			  Rcpp::Named("Start species") = R.length(),
			  Rcpp::Named("Species abundance") = spcounts);
			  //	  			  Rcpp::Named("p") = natalspprob);

  //  return(spprob);
}


// [[Rcpp::export]]
Rcpp::List runModelNDDecpp(IntegerVector R,
			    int N,
			    int P,
			    int D,
			    double E,
			    double S,
			    NumericVector uprob = NumericVector::create()
			    ) {

  // patches and sp abundance counts
  IntegerVector patches = RcppArmadillo::sample(R, P, TRUE, uprob) ;
  IntegerVector oldpatches = clone(patches);
  IntegerVector patchid = seq(0, P-1);
  IntegerVector spcounts(R.length());

  // emigration adjusted counts and natal establishment prob
  NumericVector ecounts(R.length());
  NumericVector natalspprob(R.length());

  // emigration and freqeuncy per species
  NumericVector Esp(R.length());
  NumericVector Fsp(R.length());
  
  for(int i = 0; i < N; i++) {
    // those doomed to die, and those deemed the succesors
    IntegerVector Deaths = RcppArmadillo::sample(patchid, D, FALSE, uprob);
    IntegerVector Recruits(D);
    IntegerVector Rec(1);

    // fill with zeros
    spcounts.fill(0);
    ecounts.fill(0);
    Esp.fill(0);
    Fsp.fill(0);
    natalspprob.fill(0);
    Recruits.fill(0);
    
    // calculate the previous status quo    
    for (int j = 0; j < patches.length(); j++) {
          spcounts[patches[j]] += 1;
    }

    for(int sp=0; sp<R.length(); sp++){
      Fsp[sp]=(double)spcounts[sp]/(double)sum(spcounts);
      Esp[sp]= E/(1.0+exp(-(0.15-Fsp[sp])/(E/10)));
    }

    
     // for each death d institute succession 
     for (int d = 0; d < Deaths.length(); d++) {

       // Emmigration adjusted sp counts
       for (int r = 0; r < R.length(); r++) {

       ecounts[r] = (Esp[r]*(double)spcounts[r])/P;
       
	 // natal site survival
       if(patches[Deaths[d]] == r){
	  ecounts[r] = S*(1-Esp[r])+(S*(ecounts[r]/P));
       } 

       // caluclate  spprob for each species for natal site d
       	 natalspprob[r] = ecounts[r]/sum(ecounts);

	 if(isnan(natalspprob[r])) {
	   natalspprob[r] = 0;
	 }
     }

         Rec = RcppArmadillo::sample(R, 1, FALSE, natalspprob);
         Recruits[d] = mean(Rec);
	 
     }
     
         patches[Deaths] = Recruits;
   }

  //Rcpp::Named("Start Patches") = oldpatche,
    return Rcpp::List::create(Rcpp::Named("Patches") = patches,
			  Rcpp::Named("Start species") = R.length(),
			  Rcpp::Named("Species abundance") = spcounts);

//return Rcpp::List::create(Rcpp::Named("Fsp") = Fsp,
//			  Rcpp::Named("Esp") = Esp,
//			  Rcpp::Named("natalspprob") = natalspprob);
//  
}
