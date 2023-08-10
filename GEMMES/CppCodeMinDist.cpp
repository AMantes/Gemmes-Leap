// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include<math.h>
#include <Rcpp.h>
using namespace Rcpp; 
using namespace arma; 
using namespace std;

#define dim 146
#define dimIv 834
#define dimOut 1126



/******************************/
/*  Optimization using CMAES  */
/******************************/

template<typename T>
class CMAES
{
public:
  //Constructor
  CMAES(T modelIn,
        const unsigned int lambdaIn,
        double sigmaIn, 
        const unsigned int nIterMaxIn, 
        const double tolIn,
        const bool useParallelIn) : 
  model(modelIn),
  lambda(lambdaIn), sigmaMem(sigmaIn), nIterMax(nIterMaxIn),
  tol(tolIn), tolFitness(1e-3), useParallel(useParallelIn) {}
  
  // Minimizes (minus) the log-likelihood using CMA-ES
  double Optimize(arma::vec& m) {
    const unsigned int n = m.n_elem;
    double sigma = sigmaMem;
    double mu = floor(lambda/2);
    arma::vec w = createw(lambda, n);
    double muEff = pow(accu(w(span(0, mu-1))), 2)/accu(pow(w(span(0, mu-1)), 2)); 
    double EN01 = std::pow(n, 0.5)*(1-1/(4*n)+1/(21*std::pow(n, 2)));
    double csig = (muEff+2)/(n + muEff+5);
    double dsig = 1 + 2*std::max(0.0, std::pow((muEff-1)/(n+1), 0.5)-1) + csig;
    double cc = (4 + muEff/n)/(n + 4 + 2*muEff/n);
    double c1 = 2/(std::pow(n + 1.3, 2) + muEff);
    double cb = 1;
    double cmu = std::min(1-c1, 2*(muEff- 2 + 1/muEff)/(std::pow(2+n, 2) + 2*muEff/2));
    double hsig = 0;
    unsigned int indexPastVal = 0;
    unsigned int maxPastVal = 95 + floor(30*n/lambda);
    unsigned int indexNewestVal = 0;
    arma::vec NewestVal(25);
    arma::vec OldestVal(25);
    
    unsigned int nTrial;
    arma::vec yw(n);
    arma::mat C = eye(n, n);
    arma::mat zlambda(n, lambda);
    arma::mat ylambda(n, lambda);
    arma::mat xlambda(n, lambda);
    arma::vec flambda(lambda);
    arma::vec psig = zeros(n);
    arma::vec pc = zeros(n);
    arma::Col<unsigned int> flambdaIndex(lambda);
    arma::vec D(n);
    arma::mat B(n, n);
    arma::vec tempM;
    arma::vec pastVal(maxPastVal);
    arma::Col<unsigned int> pastAvgResults(n, fill::zeros);
    arma::vec wo = w;
    unsigned int indexPastAvgResults = 0;
    
    // MAIN LOOP //
    for (unsigned int i=0; i<nIterMax; ++i) {
      // GENERATE NEW GENERATION OF CANDIDATES //
      eig_sym(D, B, C);
      D = pow(D, 0.5);
      // We can't have a parallel construct nested in a single construct, need to start a new parallel regin at each iteration of the for loop
#pragma omp parallel for default(shared) private(nTrial, tempM) if(useParallel)
      for (unsigned int j=0; j<lambda; ++j) { //Generate offsprings and compute their fitness
        nTrial = 0;
        do {
          zlambda.col(j) = randn(n);
          ylambda.col(j) = B*diagmat(D)*zlambda.col(j);
          xlambda.col(j) = m + sigma*ylambda.col(j);
          tempM = vec(xlambda.colptr(j), n, false, false);
          flambda(j) = model.Evaluate(tempM);
          nTrial++;
        } while (flambda(j)>=1e50 && nTrial<100); // Re-generate points untill the LLK is successfully computed
      }
      flambdaIndex = arma::sort_index(flambda, "ascend"); // sort points by fitness
      yw = ylambda.cols(flambdaIndex(arma::span(0, mu-1)))*w(span(0, mu-1));
      // SELECTION AND RECOMBINATION
      m = m + cb*sigma*yw; //update m taking the mean of the selected points
      // STEP-SIZE CONTROL
      psig = (1 - csig)*psig + std::pow(csig*(2-csig)*muEff, 0.5)*B*zlambda.cols(flambdaIndex(span(0, mu-1)))*w(span(0, mu-1));
      sigma = sigma*exp(csig/dsig*(norm(psig)/EN01-1));
      //COVARIANCE MATRIX ADAPTATION
      for (unsigned int j=mu;j<lambda; ++j) {
        wo(j)= w(j)*n/pow(norm(B*zlambda.col(flambdaIndex(j)), 2), 2);
      }
      hsig = ((norm(psig)/std::pow(1 - std::pow(1.0-csig, 2.0*(n+1.0)), 0.5)) < ((1.4 + 2/(n+1))*EN01)) ? 1.0 : 0.0;
      pc = (1-cc)*pc + hsig*std::pow(cc*(2-cc)*muEff, 0.5)*ylambda.cols(flambdaIndex(arma::span(0, mu-1)))*w(span(0, mu-1));
      C = (1 + c1*(1-hsig)*cc*(2-cc) - c1 - cmu*accu(wo))*C + c1*pc*pc.t() + cmu*ylambda.cols(flambdaIndex)*diagmat(wo)*ylambda.cols(flambdaIndex).t();
      C = trimatu(C) + trimatl(C, -1); //Enforce symmetry of C
      //if (std::abs((flambda(flambdaIndex(0))-flambda(flambdaIndex(mu-1)))/flambda(flambdaIndex(0)))<tolFitness) { //check for flatness
      //  sigma = sigma*exp(0.2+csig/dsig);
      //  Rcpp::Rcout<<"Warning: Flat fitness function"<<std::endl;
      //}
      //TERMINATION CRITERIA
      if (max(pc)<tol && max(sigma*C.diag())<tol) {
        if (useParallel==true) Rcpp::Rcout<<"Terminating successfully: pc and C's values are below threshold, at iteration: "<<i<<std::endl;
        break;
      }
      if (max(D)>pow(10, 14)*min(D)) {
        if (useParallel==true) Rcpp::Rcout<<"Terminating successfully: Excessive condition number of the covariance matrix, at iteration: "<<i<<std::endl;
        break;
      }
      (indexPastVal<maxPastVal-1) ? indexPastVal++ : indexPastVal = 0;
      (indexNewestVal<25-1) ? indexNewestVal++ : indexNewestVal=0;
      OldestVal(indexNewestVal) = pastVal(indexPastVal);
      pastVal(indexPastVal) = model.Evaluate(m);
      NewestVal(indexNewestVal) = pastVal(indexPastVal);
      if (i>maxPastVal+25 && std::abs(mean(OldestVal) - mean(NewestVal))<tol & std::abs(median(OldestVal) - median(NewestVal))<tol) {
        if (useParallel==true) Rcpp::Rcout<<"Terminating successfully: Stagnation of function value at m, at iteration: "<<i<<std::endl;
        break;
      }
      (indexPastAvgResults<n-1) ? indexPastAvgResults++ : indexPastAvgResults = 0;
      flambda(flambdaIndex(mu-1)) - flambda(flambdaIndex(0))<tol ? pastAvgResults(indexPastAvgResults)=1 : pastAvgResults(indexPastAvgResults) = 0;
      if (accu(pastAvgResults)>n/3) {
        if (useParallel==true) Rcpp::Rcout<<"Terminating successfully: Constant value for all best mu points, at iteration: "<<i<<std::endl;
        break;
      }
      if (norm(0.1*sigma*D(n-1)*B.col(n-1)/m)<tol) {
        if (useParallel==true) Rcpp::Rcout<<"Terminating successfully: A shock on the principal axis of C does not significantly change m, at iteration: "<<i<<std::endl;
        break;
      }
      //Rcpp::Rcout<<"generation: "<<i<<"flambda: "<<flambda(flambdaIndex(0))<<" "<<flambda(flambdaIndex(mu-1))<<" "<<pastVal(indexPastVal)<<" " <<max(abs(m))<<std::endl;
      if (flambda(flambdaIndex(mu-1))>1e45) break;
    }
    
    // generate and return output (best solution m and value at the best solution)
    return pastVal(indexPastVal);
  }
  
  void editModel(T& newModel) {
    model = newModel;
  }
  
private:
  T model; // a class containing a member function double Evaluate(arma::vec) to minimize
  const unsigned int lambda; 
  double sigmaMem; // used to reset sigma between two calls of the function with he genetic algorithm
  const unsigned int nIterMax;
  const double tol;
  const double tolFitness;
  const bool useParallel;
  const bool destroyPointers = false;
  arma::vec createw(const unsigned int LAMBDA, const unsigned int N) { // used to initialize w
    const unsigned int MU = floor(LAMBDA/2);
    arma::vec out(LAMBDA);
    for (unsigned int i=0; i<LAMBDA; ++i) {
      out(i) = log((LAMBDA+1)/2) - log(i+1);
    }
    double MUEFF = pow(accu(out(span(0, MU-1))), 2)/accu(pow(out(span(0, MU-1)), 2));
    double MUEFFminus = pow(accu(out(span(MU, LAMBDA-1))), 2)/accu(pow(out(span(MU, LAMBDA - 1)), 2));
    double C1 = 2/(std::pow(N + 1.3, 2) + MUEFF);
    double CC = (4 + MUEFF/N)/(N + 4 + 2*MUEFF/N);
    double CMU = std::min(1-C1, 2*(MUEFF- 2 + 1/MUEFF)/(std::pow(2+N, 2) + 2*MUEFF/2));
    double ALPHAMUminus = 1 + C1/CMU;
    double ALPHAMUEFFminus = 1 + 2*MUEFFminus/(MUEFF+2);
    double ALPHAPOSDEFminus = (1 - C1 - CMU)/(N*CMU);
    double SUMTEMPplus = accu(out(span(0, MU-1)));
    double SUMTEMPminus = std::abs(accu(out(span(MU, LAMBDA-1))));
    for (unsigned int i = 0; i<MU;++i) {
      out(i) *= 1/SUMTEMPplus;
    }
    for (unsigned int i = MU; i<LAMBDA;++i) {
      out(i) *= std::min(std::min(ALPHAMUminus, ALPHAMUEFFminus), ALPHAPOSDEFminus)/SUMTEMPminus;
    }
    return out;
  }
};


// SYSTEM OF DIFFERENTIAL EQUATION
void Func(double t, double* y, double* parms, double* ydot, double* x, double** dataExogVar, double** exogSamplingTime, int nExogVar, int* comptExogVar) {

		for (unsigned int it=0;it<nExogVar; it++) {
			while (t>exogSamplingTime[it][comptExogVar[it]]) comptExogVar[it]++;
			comptExogVar[it]--;
		}
ydot[55] = parms[428] * y[55];
ydot[56] = parms[427] * y[56];
x[0] = dataExogVar[0][comptExogVar[0]] + (dataExogVar[0][comptExogVar[0] + 1] - dataExogVar[0][comptExogVar[0]]) * (t - exogSamplingTime[0][comptExogVar[0]])/(exogSamplingTime[0][comptExogVar[0] + 
    1] - exogSamplingTime[0][comptExogVar[0]]);
x[1] = dataExogVar[1][comptExogVar[1]] + (dataExogVar[1][comptExogVar[1] + 1] - dataExogVar[1][comptExogVar[1]]) * (t - exogSamplingTime[1][comptExogVar[1]])/(exogSamplingTime[1][comptExogVar[1] + 
    1] - exogSamplingTime[1][comptExogVar[1]]);
x[2] = dataExogVar[2][comptExogVar[2]] + (dataExogVar[2][comptExogVar[2] + 1] - dataExogVar[2][comptExogVar[2]]) * (t - exogSamplingTime[2][comptExogVar[2]])/(exogSamplingTime[2][comptExogVar[2] + 
    1] - exogSamplingTime[2][comptExogVar[2]]);
x[3] = dataExogVar[3][comptExogVar[3]] + (dataExogVar[3][comptExogVar[3] + 1] - dataExogVar[3][comptExogVar[3]]) * (t - exogSamplingTime[3][comptExogVar[3]])/(exogSamplingTime[3][comptExogVar[3] + 
    1] - exogSamplingTime[3][comptExogVar[3]]);
x[4] = dataExogVar[4][comptExogVar[4]] + (dataExogVar[4][comptExogVar[4] + 1] - dataExogVar[4][comptExogVar[4]]) * (t - exogSamplingTime[4][comptExogVar[4]])/(exogSamplingTime[4][comptExogVar[4] + 
    1] - exogSamplingTime[4][comptExogVar[4]]);
x[5] = dataExogVar[5][comptExogVar[5]] + (dataExogVar[5][comptExogVar[5] + 1] - dataExogVar[5][comptExogVar[5]]) * (t - exogSamplingTime[5][comptExogVar[5]])/(exogSamplingTime[5][comptExogVar[5] + 
    1] - exogSamplingTime[5][comptExogVar[5]]);
x[6] = dataExogVar[6][comptExogVar[6]] + (dataExogVar[6][comptExogVar[6] + 1] - dataExogVar[6][comptExogVar[6]]) * (t - exogSamplingTime[6][comptExogVar[6]])/(exogSamplingTime[6][comptExogVar[6] + 
    1] - exogSamplingTime[6][comptExogVar[6]]);
x[7] = dataExogVar[7][comptExogVar[7]] + (dataExogVar[7][comptExogVar[7] + 1] - dataExogVar[7][comptExogVar[7]]) * (t - exogSamplingTime[7][comptExogVar[7]])/(exogSamplingTime[7][comptExogVar[7] + 
    1] - exogSamplingTime[7][comptExogVar[7]]);
x[8] = dataExogVar[8][comptExogVar[8]] + (dataExogVar[8][comptExogVar[8] + 1] - dataExogVar[8][comptExogVar[8]]) * (t - exogSamplingTime[8][comptExogVar[8]])/(exogSamplingTime[8][comptExogVar[8] + 
    1] - exogSamplingTime[8][comptExogVar[8]]);
x[9] = dataExogVar[9][comptExogVar[9]] + (dataExogVar[9][comptExogVar[9] + 1] - dataExogVar[9][comptExogVar[9]]) * (t - exogSamplingTime[9][comptExogVar[9]])/(exogSamplingTime[9][comptExogVar[9] + 
    1] - exogSamplingTime[9][comptExogVar[9]]);
x[10] = dataExogVar[10][comptExogVar[10]] + (dataExogVar[10][comptExogVar[10] + 1] - dataExogVar[10][comptExogVar[10]]) * (t - exogSamplingTime[10][comptExogVar[10]])/(exogSamplingTime[10][comptExogVar[10] + 
    1] - exogSamplingTime[10][comptExogVar[10]]);
x[11] = dataExogVar[11][comptExogVar[11]] + (dataExogVar[11][comptExogVar[11] + 1] - dataExogVar[11][comptExogVar[11]]) * (t - exogSamplingTime[11][comptExogVar[11]])/(exogSamplingTime[11][comptExogVar[11] + 
    1] - exogSamplingTime[11][comptExogVar[11]]);
x[12] = dataExogVar[12][comptExogVar[12]] + (dataExogVar[12][comptExogVar[12] + 1] - dataExogVar[12][comptExogVar[12]]) * (t - exogSamplingTime[12][comptExogVar[12]])/(exogSamplingTime[12][comptExogVar[12] + 
    1] - exogSamplingTime[12][comptExogVar[12]]);
x[13] = dataExogVar[13][comptExogVar[13]] + (dataExogVar[13][comptExogVar[13] + 1] - dataExogVar[13][comptExogVar[13]]) * (t - exogSamplingTime[13][comptExogVar[13]])/(exogSamplingTime[13][comptExogVar[13] + 
    1] - exogSamplingTime[13][comptExogVar[13]]);
x[14] = dataExogVar[14][comptExogVar[14]] + (dataExogVar[14][comptExogVar[14] + 1] - dataExogVar[14][comptExogVar[14]]) * (t - exogSamplingTime[14][comptExogVar[14]])/(exogSamplingTime[14][comptExogVar[14] + 
    1] - exogSamplingTime[14][comptExogVar[14]]);
x[15] = dataExogVar[15][comptExogVar[15]] + (dataExogVar[15][comptExogVar[15] + 1] - dataExogVar[15][comptExogVar[15]]) * (t - exogSamplingTime[15][comptExogVar[15]])/(exogSamplingTime[15][comptExogVar[15] + 
    1] - exogSamplingTime[15][comptExogVar[15]]);
x[16] = dataExogVar[16][comptExogVar[16]] + (dataExogVar[16][comptExogVar[16] + 1] - dataExogVar[16][comptExogVar[16]]) * (t - exogSamplingTime[16][comptExogVar[16]])/(exogSamplingTime[16][comptExogVar[16] + 
    1] - exogSamplingTime[16][comptExogVar[16]]);
x[17] = dataExogVar[17][comptExogVar[17]] + (dataExogVar[17][comptExogVar[17] + 1] - dataExogVar[17][comptExogVar[17]]) * (t - exogSamplingTime[17][comptExogVar[17]])/(exogSamplingTime[17][comptExogVar[17] + 
    1] - exogSamplingTime[17][comptExogVar[17]]);
x[18] = dataExogVar[18][comptExogVar[18]] + (dataExogVar[18][comptExogVar[18] + 1] - dataExogVar[18][comptExogVar[18]]) * (t - exogSamplingTime[18][comptExogVar[18]])/(exogSamplingTime[18][comptExogVar[18] + 
    1] - exogSamplingTime[18][comptExogVar[18]]);
x[19] = dataExogVar[19][comptExogVar[19]] + (dataExogVar[19][comptExogVar[19] + 1] - dataExogVar[19][comptExogVar[19]]) * (t - exogSamplingTime[19][comptExogVar[19]])/(exogSamplingTime[19][comptExogVar[19] + 
    1] - exogSamplingTime[19][comptExogVar[19]]);
x[20] = dataExogVar[20][comptExogVar[20]] + (dataExogVar[20][comptExogVar[20] + 1] - dataExogVar[20][comptExogVar[20]]) * (t - exogSamplingTime[20][comptExogVar[20]])/(exogSamplingTime[20][comptExogVar[20] + 
    1] - exogSamplingTime[20][comptExogVar[20]]);
x[21] = dataExogVar[21][comptExogVar[21]] + (dataExogVar[21][comptExogVar[21] + 1] - dataExogVar[21][comptExogVar[21]]) * (t - exogSamplingTime[21][comptExogVar[21]])/(exogSamplingTime[21][comptExogVar[21] + 
    1] - exogSamplingTime[21][comptExogVar[21]]);
x[22] = dataExogVar[22][comptExogVar[22]] + (dataExogVar[22][comptExogVar[22] + 1] - dataExogVar[22][comptExogVar[22]]) * (t - exogSamplingTime[22][comptExogVar[22]])/(exogSamplingTime[22][comptExogVar[22] + 
    1] - exogSamplingTime[22][comptExogVar[22]]);
x[23] = dataExogVar[23][comptExogVar[23]] + (dataExogVar[23][comptExogVar[23] + 1] - dataExogVar[23][comptExogVar[23]]) * (t - exogSamplingTime[23][comptExogVar[23]])/(exogSamplingTime[23][comptExogVar[23] + 
    1] - exogSamplingTime[23][comptExogVar[23]]);
x[24] = dataExogVar[24][comptExogVar[24]] + (dataExogVar[24][comptExogVar[24] + 1] - dataExogVar[24][comptExogVar[24]]) * (t - exogSamplingTime[24][comptExogVar[24]])/(exogSamplingTime[24][comptExogVar[24] + 
    1] - exogSamplingTime[24][comptExogVar[24]]);
x[25] = dataExogVar[25][comptExogVar[25]] + (dataExogVar[25][comptExogVar[25] + 1] - dataExogVar[25][comptExogVar[25]]) * (t - exogSamplingTime[25][comptExogVar[25]])/(exogSamplingTime[25][comptExogVar[25] + 
    1] - exogSamplingTime[25][comptExogVar[25]]);
x[26] = dataExogVar[26][comptExogVar[26]] + (dataExogVar[26][comptExogVar[26] + 1] - dataExogVar[26][comptExogVar[26]]) * (t - exogSamplingTime[26][comptExogVar[26]])/(exogSamplingTime[26][comptExogVar[26] + 
    1] - exogSamplingTime[26][comptExogVar[26]]);
x[27] = dataExogVar[27][comptExogVar[27]] + (dataExogVar[27][comptExogVar[27] + 1] - dataExogVar[27][comptExogVar[27]]) * (t - exogSamplingTime[27][comptExogVar[27]])/(exogSamplingTime[27][comptExogVar[27] + 
    1] - exogSamplingTime[27][comptExogVar[27]]);
x[28] = dataExogVar[28][comptExogVar[28]] + (dataExogVar[28][comptExogVar[28] + 1] - dataExogVar[28][comptExogVar[28]]) * (t - exogSamplingTime[28][comptExogVar[28]])/(exogSamplingTime[28][comptExogVar[28] + 
    1] - exogSamplingTime[28][comptExogVar[28]]);
x[29] = dataExogVar[29][comptExogVar[29]] + (dataExogVar[29][comptExogVar[29] + 1] - dataExogVar[29][comptExogVar[29]]) * (t - exogSamplingTime[29][comptExogVar[29]])/(exogSamplingTime[29][comptExogVar[29] + 
    1] - exogSamplingTime[29][comptExogVar[29]]);
x[30] = dataExogVar[30][comptExogVar[30]] + (dataExogVar[30][comptExogVar[30] + 1] - dataExogVar[30][comptExogVar[30]]) * (t - exogSamplingTime[30][comptExogVar[30]])/(exogSamplingTime[30][comptExogVar[30] + 
    1] - exogSamplingTime[30][comptExogVar[30]]);
x[31] = dataExogVar[31][comptExogVar[31]] + (dataExogVar[31][comptExogVar[31] + 1] - dataExogVar[31][comptExogVar[31]]) * (t - exogSamplingTime[31][comptExogVar[31]])/(exogSamplingTime[31][comptExogVar[31] + 
    1] - exogSamplingTime[31][comptExogVar[31]]);
x[32] = dataExogVar[32][comptExogVar[32]] + (dataExogVar[32][comptExogVar[32] + 1] - dataExogVar[32][comptExogVar[32]]) * (t - exogSamplingTime[32][comptExogVar[32]])/(exogSamplingTime[32][comptExogVar[32] + 
    1] - exogSamplingTime[32][comptExogVar[32]]);
x[33] = dataExogVar[33][comptExogVar[33]] + (dataExogVar[33][comptExogVar[33] + 1] - dataExogVar[33][comptExogVar[33]]) * (t - exogSamplingTime[33][comptExogVar[33]])/(exogSamplingTime[33][comptExogVar[33] + 
    1] - exogSamplingTime[33][comptExogVar[33]]);
x[34] = dataExogVar[34][comptExogVar[34]] + (dataExogVar[34][comptExogVar[34] + 1] - dataExogVar[34][comptExogVar[34]]) * (t - exogSamplingTime[34][comptExogVar[34]])/(exogSamplingTime[34][comptExogVar[34] + 
    1] - exogSamplingTime[34][comptExogVar[34]]);
x[35] = dataExogVar[35][comptExogVar[35]] + (dataExogVar[35][comptExogVar[35] + 1] - dataExogVar[35][comptExogVar[35]]) * (t - exogSamplingTime[35][comptExogVar[35]])/(exogSamplingTime[35][comptExogVar[35] + 
    1] - exogSamplingTime[35][comptExogVar[35]]);
x[36] = dataExogVar[36][comptExogVar[36]] + (dataExogVar[36][comptExogVar[36] + 1] - dataExogVar[36][comptExogVar[36]]) * (t - exogSamplingTime[36][comptExogVar[36]])/(exogSamplingTime[36][comptExogVar[36] + 
    1] - exogSamplingTime[36][comptExogVar[36]]);
x[37] = dataExogVar[37][comptExogVar[37]] + (dataExogVar[37][comptExogVar[37] + 1] - dataExogVar[37][comptExogVar[37]]) * (t - exogSamplingTime[37][comptExogVar[37]])/(exogSamplingTime[37][comptExogVar[37] + 
    1] - exogSamplingTime[37][comptExogVar[37]]);
x[38] = dataExogVar[38][comptExogVar[38]] + (dataExogVar[38][comptExogVar[38] + 1] - dataExogVar[38][comptExogVar[38]]) * (t - exogSamplingTime[38][comptExogVar[38]])/(exogSamplingTime[38][comptExogVar[38] + 
    1] - exogSamplingTime[38][comptExogVar[38]]);
x[39] = dataExogVar[39][comptExogVar[39]] + (dataExogVar[39][comptExogVar[39] + 1] - dataExogVar[39][comptExogVar[39]]) * (t - exogSamplingTime[39][comptExogVar[39]])/(exogSamplingTime[39][comptExogVar[39] + 
    1] - exogSamplingTime[39][comptExogVar[39]]);
x[40] = dataExogVar[40][comptExogVar[40]] + (dataExogVar[40][comptExogVar[40] + 1] - dataExogVar[40][comptExogVar[40]]) * (t - exogSamplingTime[40][comptExogVar[40]])/(exogSamplingTime[40][comptExogVar[40] + 
    1] - exogSamplingTime[40][comptExogVar[40]]);
x[41] = dataExogVar[41][comptExogVar[41]] + (dataExogVar[41][comptExogVar[41] + 1] - dataExogVar[41][comptExogVar[41]]) * (t - exogSamplingTime[41][comptExogVar[41]])/(exogSamplingTime[41][comptExogVar[41] + 
    1] - exogSamplingTime[41][comptExogVar[41]]);
x[42] = dataExogVar[42][comptExogVar[42]] + (dataExogVar[42][comptExogVar[42] + 1] - dataExogVar[42][comptExogVar[42]]) * (t - exogSamplingTime[42][comptExogVar[42]])/(exogSamplingTime[42][comptExogVar[42] + 
    1] - exogSamplingTime[42][comptExogVar[42]]);
x[43] = dataExogVar[43][comptExogVar[43]] + (dataExogVar[43][comptExogVar[43] + 1] - dataExogVar[43][comptExogVar[43]]) * (t - exogSamplingTime[43][comptExogVar[43]])/(exogSamplingTime[43][comptExogVar[43] + 
    1] - exogSamplingTime[43][comptExogVar[43]]);
x[44] = dataExogVar[44][comptExogVar[44]] + (dataExogVar[44][comptExogVar[44] + 1] - dataExogVar[44][comptExogVar[44]]) * (t - exogSamplingTime[44][comptExogVar[44]])/(exogSamplingTime[44][comptExogVar[44] + 
    1] - exogSamplingTime[44][comptExogVar[44]]);
x[45] = dataExogVar[45][comptExogVar[45]] + (dataExogVar[45][comptExogVar[45] + 1] - dataExogVar[45][comptExogVar[45]]) * (t - exogSamplingTime[45][comptExogVar[45]])/(exogSamplingTime[45][comptExogVar[45] + 
    1] - exogSamplingTime[45][comptExogVar[45]]);
x[46] = dataExogVar[46][comptExogVar[46]] + (dataExogVar[46][comptExogVar[46] + 1] - dataExogVar[46][comptExogVar[46]]) * (t - exogSamplingTime[46][comptExogVar[46]])/(exogSamplingTime[46][comptExogVar[46] + 
    1] - exogSamplingTime[46][comptExogVar[46]]);
x[47] = dataExogVar[47][comptExogVar[47]] + (dataExogVar[47][comptExogVar[47] + 1] - dataExogVar[47][comptExogVar[47]]) * (t - exogSamplingTime[47][comptExogVar[47]])/(exogSamplingTime[47][comptExogVar[47] + 
    1] - exogSamplingTime[47][comptExogVar[47]]);
x[48] = dataExogVar[48][comptExogVar[48]] + (dataExogVar[48][comptExogVar[48] + 1] - dataExogVar[48][comptExogVar[48]]) * (t - exogSamplingTime[48][comptExogVar[48]])/(exogSamplingTime[48][comptExogVar[48] + 
    1] - exogSamplingTime[48][comptExogVar[48]]);
x[49] = dataExogVar[49][comptExogVar[49]] + (dataExogVar[49][comptExogVar[49] + 1] - dataExogVar[49][comptExogVar[49]]) * (t - exogSamplingTime[49][comptExogVar[49]])/(exogSamplingTime[49][comptExogVar[49] + 
    1] - exogSamplingTime[49][comptExogVar[49]]);
x[50] = dataExogVar[50][comptExogVar[50]] + (dataExogVar[50][comptExogVar[50] + 1] - dataExogVar[50][comptExogVar[50]]) * (t - exogSamplingTime[50][comptExogVar[50]])/(exogSamplingTime[50][comptExogVar[50] + 
    1] - exogSamplingTime[50][comptExogVar[50]]);
x[51] = dataExogVar[51][comptExogVar[51]] + (dataExogVar[51][comptExogVar[51] + 1] - dataExogVar[51][comptExogVar[51]]) * (t - exogSamplingTime[51][comptExogVar[51]])/(exogSamplingTime[51][comptExogVar[51] + 
    1] - exogSamplingTime[51][comptExogVar[51]]);
x[52] = dataExogVar[52][comptExogVar[52]] + (dataExogVar[52][comptExogVar[52] + 1] - dataExogVar[52][comptExogVar[52]]) * (t - exogSamplingTime[52][comptExogVar[52]])/(exogSamplingTime[52][comptExogVar[52] + 
    1] - exogSamplingTime[52][comptExogVar[52]]);
x[53] = dataExogVar[53][comptExogVar[53]] + (dataExogVar[53][comptExogVar[53] + 1] - dataExogVar[53][comptExogVar[53]]) * (t - exogSamplingTime[53][comptExogVar[53]])/(exogSamplingTime[53][comptExogVar[53] + 
    1] - exogSamplingTime[53][comptExogVar[53]]);
x[54] = dataExogVar[54][comptExogVar[54]] + (dataExogVar[54][comptExogVar[54] + 1] - dataExogVar[54][comptExogVar[54]]) * (t - exogSamplingTime[54][comptExogVar[54]])/(exogSamplingTime[54][comptExogVar[54] + 
    1] - exogSamplingTime[54][comptExogVar[54]]);
x[55] = dataExogVar[55][comptExogVar[55]] + (dataExogVar[55][comptExogVar[55] + 1] - dataExogVar[55][comptExogVar[55]]) * (t - exogSamplingTime[55][comptExogVar[55]])/(exogSamplingTime[55][comptExogVar[55] + 
    1] - exogSamplingTime[55][comptExogVar[55]]);
x[56] = dataExogVar[56][comptExogVar[56]] + (dataExogVar[56][comptExogVar[56] + 1] - dataExogVar[56][comptExogVar[56]]) * (t - exogSamplingTime[56][comptExogVar[56]])/(exogSamplingTime[56][comptExogVar[56] + 
    1] - exogSamplingTime[56][comptExogVar[56]]);
x[57] = dataExogVar[57][comptExogVar[57]] + (dataExogVar[57][comptExogVar[57] + 1] - dataExogVar[57][comptExogVar[57]]) * (t - exogSamplingTime[57][comptExogVar[57]])/(exogSamplingTime[57][comptExogVar[57] + 
    1] - exogSamplingTime[57][comptExogVar[57]]);
x[58] = dataExogVar[58][comptExogVar[58]] + (dataExogVar[58][comptExogVar[58] + 1] - dataExogVar[58][comptExogVar[58]]) * (t - exogSamplingTime[58][comptExogVar[58]])/(exogSamplingTime[58][comptExogVar[58] + 
    1] - exogSamplingTime[58][comptExogVar[58]]);
x[59] = dataExogVar[59][comptExogVar[59]] + (dataExogVar[59][comptExogVar[59] + 1] - dataExogVar[59][comptExogVar[59]]) * (t - exogSamplingTime[59][comptExogVar[59]])/(exogSamplingTime[59][comptExogVar[59] + 
    1] - exogSamplingTime[59][comptExogVar[59]]);
x[60] = dataExogVar[60][comptExogVar[60]] + (dataExogVar[60][comptExogVar[60] + 1] - dataExogVar[60][comptExogVar[60]]) * (t - exogSamplingTime[60][comptExogVar[60]])/(exogSamplingTime[60][comptExogVar[60] + 
    1] - exogSamplingTime[60][comptExogVar[60]]);
x[61] = dataExogVar[61][comptExogVar[61]] + (dataExogVar[61][comptExogVar[61] + 1] - dataExogVar[61][comptExogVar[61]]) * (t - exogSamplingTime[61][comptExogVar[61]])/(exogSamplingTime[61][comptExogVar[61] + 
    1] - exogSamplingTime[61][comptExogVar[61]]);
x[62] = dataExogVar[62][comptExogVar[62]] + (dataExogVar[62][comptExogVar[62] + 1] - dataExogVar[62][comptExogVar[62]]) * (t - exogSamplingTime[62][comptExogVar[62]])/(exogSamplingTime[62][comptExogVar[62] + 
    1] - exogSamplingTime[62][comptExogVar[62]]);
x[63] = dataExogVar[63][comptExogVar[63]] + (dataExogVar[63][comptExogVar[63] + 1] - dataExogVar[63][comptExogVar[63]]) * (t - exogSamplingTime[63][comptExogVar[63]])/(exogSamplingTime[63][comptExogVar[63] + 
    1] - exogSamplingTime[63][comptExogVar[63]]);
x[64] = dataExogVar[64][comptExogVar[64]] + (dataExogVar[64][comptExogVar[64] + 1] - dataExogVar[64][comptExogVar[64]]) * (t - exogSamplingTime[64][comptExogVar[64]])/(exogSamplingTime[64][comptExogVar[64] + 
    1] - exogSamplingTime[64][comptExogVar[64]]);
x[65] = dataExogVar[65][comptExogVar[65]] + (dataExogVar[65][comptExogVar[65] + 1] - dataExogVar[65][comptExogVar[65]]) * (t - exogSamplingTime[65][comptExogVar[65]])/(exogSamplingTime[65][comptExogVar[65] + 
    1] - exogSamplingTime[65][comptExogVar[65]]);
x[66] = dataExogVar[66][comptExogVar[66]] + (dataExogVar[66][comptExogVar[66] + 1] - dataExogVar[66][comptExogVar[66]]) * (t - exogSamplingTime[66][comptExogVar[66]])/(exogSamplingTime[66][comptExogVar[66] + 
    1] - exogSamplingTime[66][comptExogVar[66]]);
x[67] = dataExogVar[67][comptExogVar[67]] + (dataExogVar[67][comptExogVar[67] + 1] - dataExogVar[67][comptExogVar[67]]) * (t - exogSamplingTime[67][comptExogVar[67]])/(exogSamplingTime[67][comptExogVar[67] + 
    1] - exogSamplingTime[67][comptExogVar[67]]);
x[68] = dataExogVar[68][comptExogVar[68]] + (dataExogVar[68][comptExogVar[68] + 1] - dataExogVar[68][comptExogVar[68]]) * (t - exogSamplingTime[68][comptExogVar[68]])/(exogSamplingTime[68][comptExogVar[68] + 
    1] - exogSamplingTime[68][comptExogVar[68]]);
x[69] = dataExogVar[69][comptExogVar[69]] + (dataExogVar[69][comptExogVar[69] + 1] - dataExogVar[69][comptExogVar[69]]) * (t - exogSamplingTime[69][comptExogVar[69]])/(exogSamplingTime[69][comptExogVar[69] + 
    1] - exogSamplingTime[69][comptExogVar[69]]);
x[70] = dataExogVar[70][comptExogVar[70]] + (dataExogVar[70][comptExogVar[70] + 1] - dataExogVar[70][comptExogVar[70]]) * (t - exogSamplingTime[70][comptExogVar[70]])/(exogSamplingTime[70][comptExogVar[70] + 
    1] - exogSamplingTime[70][comptExogVar[70]]);
x[71] = dataExogVar[71][comptExogVar[71]] + (dataExogVar[71][comptExogVar[71] + 1] - dataExogVar[71][comptExogVar[71]]) * (t - exogSamplingTime[71][comptExogVar[71]])/(exogSamplingTime[71][comptExogVar[71] + 
    1] - exogSamplingTime[71][comptExogVar[71]]);
x[72] = dataExogVar[72][comptExogVar[72]] + (dataExogVar[72][comptExogVar[72] + 1] - dataExogVar[72][comptExogVar[72]]) * (t - exogSamplingTime[72][comptExogVar[72]])/(exogSamplingTime[72][comptExogVar[72] + 
    1] - exogSamplingTime[72][comptExogVar[72]]);
x[73] = dataExogVar[73][comptExogVar[73]] + (dataExogVar[73][comptExogVar[73] + 1] - dataExogVar[73][comptExogVar[73]]) * (t - exogSamplingTime[73][comptExogVar[73]])/(exogSamplingTime[73][comptExogVar[73] + 
    1] - exogSamplingTime[73][comptExogVar[73]]);
x[74] = dataExogVar[74][comptExogVar[74]] + (dataExogVar[74][comptExogVar[74] + 1] - dataExogVar[74][comptExogVar[74]]) * (t - exogSamplingTime[74][comptExogVar[74]])/(exogSamplingTime[74][comptExogVar[74] + 
    1] - exogSamplingTime[74][comptExogVar[74]]);
x[75] = dataExogVar[75][comptExogVar[75]] + (dataExogVar[75][comptExogVar[75] + 1] - dataExogVar[75][comptExogVar[75]]) * (t - exogSamplingTime[75][comptExogVar[75]])/(exogSamplingTime[75][comptExogVar[75] + 
    1] - exogSamplingTime[75][comptExogVar[75]]);
x[76] = dataExogVar[76][comptExogVar[76]] + (dataExogVar[76][comptExogVar[76] + 1] - dataExogVar[76][comptExogVar[76]]) * (t - exogSamplingTime[76][comptExogVar[76]])/(exogSamplingTime[76][comptExogVar[76] + 
    1] - exogSamplingTime[76][comptExogVar[76]]);
x[77] = dataExogVar[77][comptExogVar[77]] + (dataExogVar[77][comptExogVar[77] + 1] - dataExogVar[77][comptExogVar[77]]) * (t - exogSamplingTime[77][comptExogVar[77]])/(exogSamplingTime[77][comptExogVar[77] + 
    1] - exogSamplingTime[77][comptExogVar[77]]);
x[78] = dataExogVar[78][comptExogVar[78]] + (dataExogVar[78][comptExogVar[78] + 1] - dataExogVar[78][comptExogVar[78]]) * (t - exogSamplingTime[78][comptExogVar[78]])/(exogSamplingTime[78][comptExogVar[78] + 
    1] - exogSamplingTime[78][comptExogVar[78]]);
x[79] = dataExogVar[79][comptExogVar[79]] + (dataExogVar[79][comptExogVar[79] + 1] - dataExogVar[79][comptExogVar[79]]) * (t - exogSamplingTime[79][comptExogVar[79]])/(exogSamplingTime[79][comptExogVar[79] + 
    1] - exogSamplingTime[79][comptExogVar[79]]);
x[80] = dataExogVar[80][comptExogVar[80]] + (dataExogVar[80][comptExogVar[80] + 1] - dataExogVar[80][comptExogVar[80]]) * (t - exogSamplingTime[80][comptExogVar[80]])/(exogSamplingTime[80][comptExogVar[80] + 
    1] - exogSamplingTime[80][comptExogVar[80]]);
x[81] = dataExogVar[81][comptExogVar[81]] + (dataExogVar[81][comptExogVar[81] + 1] - dataExogVar[81][comptExogVar[81]]) * (t - exogSamplingTime[81][comptExogVar[81]])/(exogSamplingTime[81][comptExogVar[81] + 
    1] - exogSamplingTime[81][comptExogVar[81]]);
x[82] = dataExogVar[82][comptExogVar[82]] + (dataExogVar[82][comptExogVar[82] + 1] - dataExogVar[82][comptExogVar[82]]) * (t - exogSamplingTime[82][comptExogVar[82]])/(exogSamplingTime[82][comptExogVar[82] + 
    1] - exogSamplingTime[82][comptExogVar[82]]);
x[83] = dataExogVar[83][comptExogVar[83]] + (dataExogVar[83][comptExogVar[83] + 1] - dataExogVar[83][comptExogVar[83]]) * (t - exogSamplingTime[83][comptExogVar[83]])/(exogSamplingTime[83][comptExogVar[83] + 
    1] - exogSamplingTime[83][comptExogVar[83]]);
x[84] = dataExogVar[84][comptExogVar[84]] + (dataExogVar[84][comptExogVar[84] + 1] - dataExogVar[84][comptExogVar[84]]) * (t - exogSamplingTime[84][comptExogVar[84]])/(exogSamplingTime[84][comptExogVar[84] + 
    1] - exogSamplingTime[84][comptExogVar[84]]);
x[85] = dataExogVar[85][comptExogVar[85]] + (dataExogVar[85][comptExogVar[85] + 1] - dataExogVar[85][comptExogVar[85]]) * (t - exogSamplingTime[85][comptExogVar[85]])/(exogSamplingTime[85][comptExogVar[85] + 
    1] - exogSamplingTime[85][comptExogVar[85]]);
x[86] = dataExogVar[86][comptExogVar[86]] + (dataExogVar[86][comptExogVar[86] + 1] - dataExogVar[86][comptExogVar[86]]) * (t - exogSamplingTime[86][comptExogVar[86]])/(exogSamplingTime[86][comptExogVar[86] + 
    1] - exogSamplingTime[86][comptExogVar[86]]);
x[87] = dataExogVar[87][comptExogVar[87]] + (dataExogVar[87][comptExogVar[87] + 1] - dataExogVar[87][comptExogVar[87]]) * (t - exogSamplingTime[87][comptExogVar[87]])/(exogSamplingTime[87][comptExogVar[87] + 
    1] - exogSamplingTime[87][comptExogVar[87]]);
x[88] = dataExogVar[88][comptExogVar[88]] + (dataExogVar[88][comptExogVar[88] + 1] - dataExogVar[88][comptExogVar[88]]) * (t - exogSamplingTime[88][comptExogVar[88]])/(exogSamplingTime[88][comptExogVar[88] + 
    1] - exogSamplingTime[88][comptExogVar[88]]);
x[89] = dataExogVar[89][comptExogVar[89]] + (dataExogVar[89][comptExogVar[89] + 1] - dataExogVar[89][comptExogVar[89]]) * (t - exogSamplingTime[89][comptExogVar[89]])/(exogSamplingTime[89][comptExogVar[89] + 
    1] - exogSamplingTime[89][comptExogVar[89]]);
x[90] = dataExogVar[90][comptExogVar[90]] + (dataExogVar[90][comptExogVar[90] + 1] - dataExogVar[90][comptExogVar[90]]) * (t - exogSamplingTime[90][comptExogVar[90]])/(exogSamplingTime[90][comptExogVar[90] + 
    1] - exogSamplingTime[90][comptExogVar[90]]);
x[91] = dataExogVar[91][comptExogVar[91]] + (dataExogVar[91][comptExogVar[91] + 1] - dataExogVar[91][comptExogVar[91]]) * (t - exogSamplingTime[91][comptExogVar[91]])/(exogSamplingTime[91][comptExogVar[91] + 
    1] - exogSamplingTime[91][comptExogVar[91]]);
x[92] = dataExogVar[92][comptExogVar[92]] + (dataExogVar[92][comptExogVar[92] + 1] - dataExogVar[92][comptExogVar[92]]) * (t - exogSamplingTime[92][comptExogVar[92]])/(exogSamplingTime[92][comptExogVar[92] + 
    1] - exogSamplingTime[92][comptExogVar[92]]);
x[93] = dataExogVar[93][comptExogVar[93]] + (dataExogVar[93][comptExogVar[93] + 1] - dataExogVar[93][comptExogVar[93]]) * (t - exogSamplingTime[93][comptExogVar[93]])/(exogSamplingTime[93][comptExogVar[93] + 
    1] - exogSamplingTime[93][comptExogVar[93]]);
x[94] = dataExogVar[94][comptExogVar[94]] + (dataExogVar[94][comptExogVar[94] + 1] - dataExogVar[94][comptExogVar[94]]) * (t - exogSamplingTime[94][comptExogVar[94]])/(exogSamplingTime[94][comptExogVar[94] + 
    1] - exogSamplingTime[94][comptExogVar[94]]);
x[95] = dataExogVar[95][comptExogVar[95]] + (dataExogVar[95][comptExogVar[95] + 1] - dataExogVar[95][comptExogVar[95]]) * (t - exogSamplingTime[95][comptExogVar[95]])/(exogSamplingTime[95][comptExogVar[95] + 
    1] - exogSamplingTime[95][comptExogVar[95]]);
x[96] = dataExogVar[96][comptExogVar[96]] + (dataExogVar[96][comptExogVar[96] + 1] - dataExogVar[96][comptExogVar[96]]) * (t - exogSamplingTime[96][comptExogVar[96]])/(exogSamplingTime[96][comptExogVar[96] + 
    1] - exogSamplingTime[96][comptExogVar[96]]);
x[97] = dataExogVar[97][comptExogVar[97]] + (dataExogVar[97][comptExogVar[97] + 1] - dataExogVar[97][comptExogVar[97]]) * (t - exogSamplingTime[97][comptExogVar[97]])/(exogSamplingTime[97][comptExogVar[97] + 
    1] - exogSamplingTime[97][comptExogVar[97]]);
x[98] = dataExogVar[98][comptExogVar[98]] + (dataExogVar[98][comptExogVar[98] + 1] - dataExogVar[98][comptExogVar[98]]) * (t - exogSamplingTime[98][comptExogVar[98]])/(exogSamplingTime[98][comptExogVar[98] + 
    1] - exogSamplingTime[98][comptExogVar[98]]);
x[99] = dataExogVar[99][comptExogVar[99]] + (dataExogVar[99][comptExogVar[99] + 1] - dataExogVar[99][comptExogVar[99]]) * (t - exogSamplingTime[99][comptExogVar[99]])/(exogSamplingTime[99][comptExogVar[99] + 
    1] - exogSamplingTime[99][comptExogVar[99]]);
x[100] = dataExogVar[100][comptExogVar[100]] + (dataExogVar[100][comptExogVar[100] + 1] - dataExogVar[100][comptExogVar[100]]) * (t - exogSamplingTime[100][comptExogVar[100]])/(exogSamplingTime[100][comptExogVar[100] + 
    1] - exogSamplingTime[100][comptExogVar[100]]);
x[101] = dataExogVar[101][comptExogVar[101]] + (dataExogVar[101][comptExogVar[101] + 1] - dataExogVar[101][comptExogVar[101]]) * (t - exogSamplingTime[101][comptExogVar[101]])/(exogSamplingTime[101][comptExogVar[101] + 
    1] - exogSamplingTime[101][comptExogVar[101]]);
x[102] = dataExogVar[102][comptExogVar[102]] + (dataExogVar[102][comptExogVar[102] + 1] - dataExogVar[102][comptExogVar[102]]) * (t - exogSamplingTime[102][comptExogVar[102]])/(exogSamplingTime[102][comptExogVar[102] + 
    1] - exogSamplingTime[102][comptExogVar[102]]);
x[103] = dataExogVar[103][comptExogVar[103]] + (dataExogVar[103][comptExogVar[103] + 1] - dataExogVar[103][comptExogVar[103]]) * (t - exogSamplingTime[103][comptExogVar[103]])/(exogSamplingTime[103][comptExogVar[103] + 
    1] - exogSamplingTime[103][comptExogVar[103]]);
x[104] = dataExogVar[104][comptExogVar[104]] + (dataExogVar[104][comptExogVar[104] + 1] - dataExogVar[104][comptExogVar[104]]) * (t - exogSamplingTime[104][comptExogVar[104]])/(exogSamplingTime[104][comptExogVar[104] + 
    1] - exogSamplingTime[104][comptExogVar[104]]);
x[105] = dataExogVar[105][comptExogVar[105]] + (dataExogVar[105][comptExogVar[105] + 1] - dataExogVar[105][comptExogVar[105]]) * (t - exogSamplingTime[105][comptExogVar[105]])/(exogSamplingTime[105][comptExogVar[105] + 
    1] - exogSamplingTime[105][comptExogVar[105]]);
x[106] = dataExogVar[106][comptExogVar[106]] + (dataExogVar[106][comptExogVar[106] + 1] - dataExogVar[106][comptExogVar[106]]) * (t - exogSamplingTime[106][comptExogVar[106]])/(exogSamplingTime[106][comptExogVar[106] + 
    1] - exogSamplingTime[106][comptExogVar[106]]);
x[107] = dataExogVar[107][comptExogVar[107]] + (dataExogVar[107][comptExogVar[107] + 1] - dataExogVar[107][comptExogVar[107]]) * (t - exogSamplingTime[107][comptExogVar[107]])/(exogSamplingTime[107][comptExogVar[107] + 
    1] - exogSamplingTime[107][comptExogVar[107]]);
x[108] = dataExogVar[108][comptExogVar[108]] + (dataExogVar[108][comptExogVar[108] + 1] - dataExogVar[108][comptExogVar[108]]) * (t - exogSamplingTime[108][comptExogVar[108]])/(exogSamplingTime[108][comptExogVar[108] + 
    1] - exogSamplingTime[108][comptExogVar[108]]);
x[109] = dataExogVar[109][comptExogVar[109]] + (dataExogVar[109][comptExogVar[109] + 1] - dataExogVar[109][comptExogVar[109]]) * (t - exogSamplingTime[109][comptExogVar[109]])/(exogSamplingTime[109][comptExogVar[109] + 
    1] - exogSamplingTime[109][comptExogVar[109]]);
x[110] = dataExogVar[110][comptExogVar[110]] + (dataExogVar[110][comptExogVar[110] + 1] - dataExogVar[110][comptExogVar[110]]) * (t - exogSamplingTime[110][comptExogVar[110]])/(exogSamplingTime[110][comptExogVar[110] + 
    1] - exogSamplingTime[110][comptExogVar[110]]);
x[111] = dataExogVar[111][comptExogVar[111]] + (dataExogVar[111][comptExogVar[111] + 1] - dataExogVar[111][comptExogVar[111]]) * (t - exogSamplingTime[111][comptExogVar[111]])/(exogSamplingTime[111][comptExogVar[111] + 
    1] - exogSamplingTime[111][comptExogVar[111]]);
x[112] = dataExogVar[112][comptExogVar[112]] + (dataExogVar[112][comptExogVar[112] + 1] - dataExogVar[112][comptExogVar[112]]) * (t - exogSamplingTime[112][comptExogVar[112]])/(exogSamplingTime[112][comptExogVar[112] + 
    1] - exogSamplingTime[112][comptExogVar[112]]);
x[113] = dataExogVar[113][comptExogVar[113]] + (dataExogVar[113][comptExogVar[113] + 1] - dataExogVar[113][comptExogVar[113]]) * (t - exogSamplingTime[113][comptExogVar[113]])/(exogSamplingTime[113][comptExogVar[113] + 
    1] - exogSamplingTime[113][comptExogVar[113]]);
x[114] = dataExogVar[114][comptExogVar[114]] + (dataExogVar[114][comptExogVar[114] + 1] - dataExogVar[114][comptExogVar[114]]) * (t - exogSamplingTime[114][comptExogVar[114]])/(exogSamplingTime[114][comptExogVar[114] + 
    1] - exogSamplingTime[114][comptExogVar[114]]);
x[115] = dataExogVar[115][comptExogVar[115]] + (dataExogVar[115][comptExogVar[115] + 1] - dataExogVar[115][comptExogVar[115]]) * (t - exogSamplingTime[115][comptExogVar[115]])/(exogSamplingTime[115][comptExogVar[115] + 
    1] - exogSamplingTime[115][comptExogVar[115]]);
x[116] = dataExogVar[116][comptExogVar[116]] + (dataExogVar[116][comptExogVar[116] + 1] - dataExogVar[116][comptExogVar[116]]) * (t - exogSamplingTime[116][comptExogVar[116]])/(exogSamplingTime[116][comptExogVar[116] + 
    1] - exogSamplingTime[116][comptExogVar[116]]);
x[117] = dataExogVar[117][comptExogVar[117]] + (dataExogVar[117][comptExogVar[117] + 1] - dataExogVar[117][comptExogVar[117]]) * (t - exogSamplingTime[117][comptExogVar[117]])/(exogSamplingTime[117][comptExogVar[117] + 
    1] - exogSamplingTime[117][comptExogVar[117]]);
x[118] = dataExogVar[118][comptExogVar[118]] + (dataExogVar[118][comptExogVar[118] + 1] - dataExogVar[118][comptExogVar[118]]) * (t - exogSamplingTime[118][comptExogVar[118]])/(exogSamplingTime[118][comptExogVar[118] + 
    1] - exogSamplingTime[118][comptExogVar[118]]);
x[119] = dataExogVar[119][comptExogVar[119]] + (dataExogVar[119][comptExogVar[119] + 1] - dataExogVar[119][comptExogVar[119]]) * (t - exogSamplingTime[119][comptExogVar[119]])/(exogSamplingTime[119][comptExogVar[119] + 
    1] - exogSamplingTime[119][comptExogVar[119]]);
x[120] = dataExogVar[120][comptExogVar[120]] + (dataExogVar[120][comptExogVar[120] + 1] - dataExogVar[120][comptExogVar[120]]) * (t - exogSamplingTime[120][comptExogVar[120]])/(exogSamplingTime[120][comptExogVar[120] + 
    1] - exogSamplingTime[120][comptExogVar[120]]);
x[121] = dataExogVar[121][comptExogVar[121]] + (dataExogVar[121][comptExogVar[121] + 1] - dataExogVar[121][comptExogVar[121]]) * (t - exogSamplingTime[121][comptExogVar[121]])/(exogSamplingTime[121][comptExogVar[121] + 
    1] - exogSamplingTime[121][comptExogVar[121]]);
x[122] = dataExogVar[122][comptExogVar[122]] + (dataExogVar[122][comptExogVar[122] + 1] - dataExogVar[122][comptExogVar[122]]) * (t - exogSamplingTime[122][comptExogVar[122]])/(exogSamplingTime[122][comptExogVar[122] + 
    1] - exogSamplingTime[122][comptExogVar[122]]);
x[123] = dataExogVar[123][comptExogVar[123]] + (dataExogVar[123][comptExogVar[123] + 1] - dataExogVar[123][comptExogVar[123]]) * (t - exogSamplingTime[123][comptExogVar[123]])/(exogSamplingTime[123][comptExogVar[123] + 
    1] - exogSamplingTime[123][comptExogVar[123]]);
x[124] = dataExogVar[124][comptExogVar[124]] + (dataExogVar[124][comptExogVar[124] + 1] - dataExogVar[124][comptExogVar[124]]) * (t - exogSamplingTime[124][comptExogVar[124]])/(exogSamplingTime[124][comptExogVar[124] + 
    1] - exogSamplingTime[124][comptExogVar[124]]);
x[125] = dataExogVar[125][comptExogVar[125]] + (dataExogVar[125][comptExogVar[125] + 1] - dataExogVar[125][comptExogVar[125]]) * (t - exogSamplingTime[125][comptExogVar[125]])/(exogSamplingTime[125][comptExogVar[125] + 
    1] - exogSamplingTime[125][comptExogVar[125]]);
x[126] = dataExogVar[126][comptExogVar[126]] + (dataExogVar[126][comptExogVar[126] + 1] - dataExogVar[126][comptExogVar[126]]) * (t - exogSamplingTime[126][comptExogVar[126]])/(exogSamplingTime[126][comptExogVar[126] + 
    1] - exogSamplingTime[126][comptExogVar[126]]);
x[127] = dataExogVar[127][comptExogVar[127]] + (dataExogVar[127][comptExogVar[127] + 1] - dataExogVar[127][comptExogVar[127]]) * (t - exogSamplingTime[127][comptExogVar[127]])/(exogSamplingTime[127][comptExogVar[127] + 
    1] - exogSamplingTime[127][comptExogVar[127]]);
x[128] = dataExogVar[128][comptExogVar[128]] + (dataExogVar[128][comptExogVar[128] + 1] - dataExogVar[128][comptExogVar[128]]) * (t - exogSamplingTime[128][comptExogVar[128]])/(exogSamplingTime[128][comptExogVar[128] + 
    1] - exogSamplingTime[128][comptExogVar[128]]);
x[129] = dataExogVar[129][comptExogVar[129]] + (dataExogVar[129][comptExogVar[129] + 1] - dataExogVar[129][comptExogVar[129]]) * (t - exogSamplingTime[129][comptExogVar[129]])/(exogSamplingTime[129][comptExogVar[129] + 
    1] - exogSamplingTime[129][comptExogVar[129]]);
x[130] = dataExogVar[130][comptExogVar[130]] + (dataExogVar[130][comptExogVar[130] + 1] - dataExogVar[130][comptExogVar[130]]) * (t - exogSamplingTime[130][comptExogVar[130]])/(exogSamplingTime[130][comptExogVar[130] + 
    1] - exogSamplingTime[130][comptExogVar[130]]);
x[131] = dataExogVar[131][comptExogVar[131]] + (dataExogVar[131][comptExogVar[131] + 1] - dataExogVar[131][comptExogVar[131]]) * (t - exogSamplingTime[131][comptExogVar[131]])/(exogSamplingTime[131][comptExogVar[131] + 
    1] - exogSamplingTime[131][comptExogVar[131]]);
x[132] = dataExogVar[132][comptExogVar[132]] + (dataExogVar[132][comptExogVar[132] + 1] - dataExogVar[132][comptExogVar[132]]) * (t - exogSamplingTime[132][comptExogVar[132]])/(exogSamplingTime[132][comptExogVar[132] + 
    1] - exogSamplingTime[132][comptExogVar[132]]);
x[133] = dataExogVar[133][comptExogVar[133]] + (dataExogVar[133][comptExogVar[133] + 1] - dataExogVar[133][comptExogVar[133]]) * (t - exogSamplingTime[133][comptExogVar[133]])/(exogSamplingTime[133][comptExogVar[133] + 
    1] - exogSamplingTime[133][comptExogVar[133]]);
x[134] = dataExogVar[134][comptExogVar[134]] + (dataExogVar[134][comptExogVar[134] + 1] - dataExogVar[134][comptExogVar[134]]) * (t - exogSamplingTime[134][comptExogVar[134]])/(exogSamplingTime[134][comptExogVar[134] + 
    1] - exogSamplingTime[134][comptExogVar[134]]);
x[135] = dataExogVar[135][comptExogVar[135]] + (dataExogVar[135][comptExogVar[135] + 1] - dataExogVar[135][comptExogVar[135]]) * (t - exogSamplingTime[135][comptExogVar[135]])/(exogSamplingTime[135][comptExogVar[135] + 
    1] - exogSamplingTime[135][comptExogVar[135]]);
x[136] = dataExogVar[136][comptExogVar[136]] + (dataExogVar[136][comptExogVar[136] + 1] - dataExogVar[136][comptExogVar[136]]) * (t - exogSamplingTime[136][comptExogVar[136]])/(exogSamplingTime[136][comptExogVar[136] + 
    1] - exogSamplingTime[136][comptExogVar[136]]);
x[137] = dataExogVar[137][comptExogVar[137]] + (dataExogVar[137][comptExogVar[137] + 1] - dataExogVar[137][comptExogVar[137]]) * (t - exogSamplingTime[137][comptExogVar[137]])/(exogSamplingTime[137][comptExogVar[137] + 
    1] - exogSamplingTime[137][comptExogVar[137]]);
x[138] = dataExogVar[138][comptExogVar[138]] + (dataExogVar[138][comptExogVar[138] + 1] - dataExogVar[138][comptExogVar[138]]) * (t - exogSamplingTime[138][comptExogVar[138]])/(exogSamplingTime[138][comptExogVar[138] + 
    1] - exogSamplingTime[138][comptExogVar[138]]);
x[139] = dataExogVar[139][comptExogVar[139]] + (dataExogVar[139][comptExogVar[139] + 1] - dataExogVar[139][comptExogVar[139]]) * (t - exogSamplingTime[139][comptExogVar[139]])/(exogSamplingTime[139][comptExogVar[139] + 
    1] - exogSamplingTime[139][comptExogVar[139]]);
x[140] = dataExogVar[140][comptExogVar[140]] + (dataExogVar[140][comptExogVar[140] + 1] - dataExogVar[140][comptExogVar[140]]) * (t - exogSamplingTime[140][comptExogVar[140]])/(exogSamplingTime[140][comptExogVar[140] + 
    1] - exogSamplingTime[140][comptExogVar[140]]);
x[141] = dataExogVar[141][comptExogVar[141]] + (dataExogVar[141][comptExogVar[141] + 1] - dataExogVar[141][comptExogVar[141]]) * (t - exogSamplingTime[141][comptExogVar[141]])/(exogSamplingTime[141][comptExogVar[141] + 
    1] - exogSamplingTime[141][comptExogVar[141]]);
x[142] = dataExogVar[142][comptExogVar[142]] + (dataExogVar[142][comptExogVar[142] + 1] - dataExogVar[142][comptExogVar[142]]) * (t - exogSamplingTime[142][comptExogVar[142]])/(exogSamplingTime[142][comptExogVar[142] + 
    1] - exogSamplingTime[142][comptExogVar[142]]);
x[143] = dataExogVar[143][comptExogVar[143]] + (dataExogVar[143][comptExogVar[143] + 1] - dataExogVar[143][comptExogVar[143]]) * (t - exogSamplingTime[143][comptExogVar[143]])/(exogSamplingTime[143][comptExogVar[143] + 
    1] - exogSamplingTime[143][comptExogVar[143]]);
x[144] = dataExogVar[144][comptExogVar[144]] + (dataExogVar[144][comptExogVar[144] + 1] - dataExogVar[144][comptExogVar[144]]) * (t - exogSamplingTime[144][comptExogVar[144]])/(exogSamplingTime[144][comptExogVar[144] + 
    1] - exogSamplingTime[144][comptExogVar[144]]);
x[145] = dataExogVar[145][comptExogVar[145]] + (dataExogVar[145][comptExogVar[145] + 1] - dataExogVar[145][comptExogVar[145]]) * (t - exogSamplingTime[145][comptExogVar[145]])/(exogSamplingTime[145][comptExogVar[145] + 
    1] - exogSamplingTime[145][comptExogVar[145]]);
x[146] = dataExogVar[146][comptExogVar[146]] + (dataExogVar[146][comptExogVar[146] + 1] - dataExogVar[146][comptExogVar[146]]) * (t - exogSamplingTime[146][comptExogVar[146]])/(exogSamplingTime[146][comptExogVar[146] + 
    1] - exogSamplingTime[146][comptExogVar[146]]);
x[147] = dataExogVar[147][comptExogVar[147]] + (dataExogVar[147][comptExogVar[147] + 1] - dataExogVar[147][comptExogVar[147]]) * (t - exogSamplingTime[147][comptExogVar[147]])/(exogSamplingTime[147][comptExogVar[147] + 
    1] - exogSamplingTime[147][comptExogVar[147]]);
x[148] = dataExogVar[148][comptExogVar[148]] + (dataExogVar[148][comptExogVar[148] + 1] - dataExogVar[148][comptExogVar[148]]) * (t - exogSamplingTime[148][comptExogVar[148]])/(exogSamplingTime[148][comptExogVar[148] + 
    1] - exogSamplingTime[148][comptExogVar[148]]);
x[149] = dataExogVar[149][comptExogVar[149]] + (dataExogVar[149][comptExogVar[149] + 1] - dataExogVar[149][comptExogVar[149]]) * (t - exogSamplingTime[149][comptExogVar[149]])/(exogSamplingTime[149][comptExogVar[149] + 
    1] - exogSamplingTime[149][comptExogVar[149]]);
x[150] = dataExogVar[150][comptExogVar[150]] + (dataExogVar[150][comptExogVar[150] + 1] - dataExogVar[150][comptExogVar[150]]) * (t - exogSamplingTime[150][comptExogVar[150]])/(exogSamplingTime[150][comptExogVar[150] + 
    1] - exogSamplingTime[150][comptExogVar[150]]);
x[151] = dataExogVar[151][comptExogVar[151]] + (dataExogVar[151][comptExogVar[151] + 1] - dataExogVar[151][comptExogVar[151]]) * (t - exogSamplingTime[151][comptExogVar[151]])/(exogSamplingTime[151][comptExogVar[151] + 
    1] - exogSamplingTime[151][comptExogVar[151]]);
x[152] = dataExogVar[152][comptExogVar[152]] + (dataExogVar[152][comptExogVar[152] + 1] - dataExogVar[152][comptExogVar[152]]) * (t - exogSamplingTime[152][comptExogVar[152]])/(exogSamplingTime[152][comptExogVar[152] + 
    1] - exogSamplingTime[152][comptExogVar[152]]);
x[153] = dataExogVar[153][comptExogVar[153]] + (dataExogVar[153][comptExogVar[153] + 1] - dataExogVar[153][comptExogVar[153]]) * (t - exogSamplingTime[153][comptExogVar[153]])/(exogSamplingTime[153][comptExogVar[153] + 
    1] - exogSamplingTime[153][comptExogVar[153]]);
x[154] = dataExogVar[154][comptExogVar[154]] + (dataExogVar[154][comptExogVar[154] + 1] - dataExogVar[154][comptExogVar[154]]) * (t - exogSamplingTime[154][comptExogVar[154]])/(exogSamplingTime[154][comptExogVar[154] + 
    1] - exogSamplingTime[154][comptExogVar[154]]);
x[155] = dataExogVar[155][comptExogVar[155]] + (dataExogVar[155][comptExogVar[155] + 1] - dataExogVar[155][comptExogVar[155]]) * (t - exogSamplingTime[155][comptExogVar[155]])/(exogSamplingTime[155][comptExogVar[155] + 
    1] - exogSamplingTime[155][comptExogVar[155]]);
x[156] = dataExogVar[156][comptExogVar[156]] + (dataExogVar[156][comptExogVar[156] + 1] - dataExogVar[156][comptExogVar[156]]) * (t - exogSamplingTime[156][comptExogVar[156]])/(exogSamplingTime[156][comptExogVar[156] + 
    1] - exogSamplingTime[156][comptExogVar[156]]);
x[157] = dataExogVar[157][comptExogVar[157]] + (dataExogVar[157][comptExogVar[157] + 1] - dataExogVar[157][comptExogVar[157]]) * (t - exogSamplingTime[157][comptExogVar[157]])/(exogSamplingTime[157][comptExogVar[157] + 
    1] - exogSamplingTime[157][comptExogVar[157]]);
x[158] = dataExogVar[158][comptExogVar[158]] + (dataExogVar[158][comptExogVar[158] + 1] - dataExogVar[158][comptExogVar[158]]) * (t - exogSamplingTime[158][comptExogVar[158]])/(exogSamplingTime[158][comptExogVar[158] + 
    1] - exogSamplingTime[158][comptExogVar[158]]);
x[159] = dataExogVar[159][comptExogVar[159]] + (dataExogVar[159][comptExogVar[159] + 1] - dataExogVar[159][comptExogVar[159]]) * (t - exogSamplingTime[159][comptExogVar[159]])/(exogSamplingTime[159][comptExogVar[159] + 
    1] - exogSamplingTime[159][comptExogVar[159]]);
x[160] = dataExogVar[160][comptExogVar[160]] + (dataExogVar[160][comptExogVar[160] + 1] - dataExogVar[160][comptExogVar[160]]) * (t - exogSamplingTime[160][comptExogVar[160]])/(exogSamplingTime[160][comptExogVar[160] + 
    1] - exogSamplingTime[160][comptExogVar[160]]);
x[161] = dataExogVar[161][comptExogVar[161]] + (dataExogVar[161][comptExogVar[161] + 1] - dataExogVar[161][comptExogVar[161]]) * (t - exogSamplingTime[161][comptExogVar[161]])/(exogSamplingTime[161][comptExogVar[161] + 
    1] - exogSamplingTime[161][comptExogVar[161]]);
x[162] = dataExogVar[162][comptExogVar[162]] + (dataExogVar[162][comptExogVar[162] + 1] - dataExogVar[162][comptExogVar[162]]) * (t - exogSamplingTime[162][comptExogVar[162]])/(exogSamplingTime[162][comptExogVar[162] + 
    1] - exogSamplingTime[162][comptExogVar[162]]);
x[163] = dataExogVar[163][comptExogVar[163]] + (dataExogVar[163][comptExogVar[163] + 1] - dataExogVar[163][comptExogVar[163]]) * (t - exogSamplingTime[163][comptExogVar[163]])/(exogSamplingTime[163][comptExogVar[163] + 
    1] - exogSamplingTime[163][comptExogVar[163]]);
x[164] = dataExogVar[164][comptExogVar[164]] + (dataExogVar[164][comptExogVar[164] + 1] - dataExogVar[164][comptExogVar[164]]) * (t - exogSamplingTime[164][comptExogVar[164]])/(exogSamplingTime[164][comptExogVar[164] + 
    1] - exogSamplingTime[164][comptExogVar[164]]);
x[165] = dataExogVar[165][comptExogVar[165]] + (dataExogVar[165][comptExogVar[165] + 1] - dataExogVar[165][comptExogVar[165]]) * (t - exogSamplingTime[165][comptExogVar[165]])/(exogSamplingTime[165][comptExogVar[165] + 
    1] - exogSamplingTime[165][comptExogVar[165]]);
x[166] = dataExogVar[166][comptExogVar[166]] + (dataExogVar[166][comptExogVar[166] + 1] - dataExogVar[166][comptExogVar[166]]) * (t - exogSamplingTime[166][comptExogVar[166]])/(exogSamplingTime[166][comptExogVar[166] + 
    1] - exogSamplingTime[166][comptExogVar[166]]);
x[167] = dataExogVar[167][comptExogVar[167]] + (dataExogVar[167][comptExogVar[167] + 1] - dataExogVar[167][comptExogVar[167]]) * (t - exogSamplingTime[167][comptExogVar[167]])/(exogSamplingTime[167][comptExogVar[167] + 
    1] - exogSamplingTime[167][comptExogVar[167]]);
x[168] = dataExogVar[168][comptExogVar[168]] + (dataExogVar[168][comptExogVar[168] + 1] - dataExogVar[168][comptExogVar[168]]) * (t - exogSamplingTime[168][comptExogVar[168]])/(exogSamplingTime[168][comptExogVar[168] + 
    1] - exogSamplingTime[168][comptExogVar[168]]);
x[169] = dataExogVar[169][comptExogVar[169]] + (dataExogVar[169][comptExogVar[169] + 1] - dataExogVar[169][comptExogVar[169]]) * (t - exogSamplingTime[169][comptExogVar[169]])/(exogSamplingTime[169][comptExogVar[169] + 
    1] - exogSamplingTime[169][comptExogVar[169]]);
x[170] = dataExogVar[170][comptExogVar[170]] + (dataExogVar[170][comptExogVar[170] + 1] - dataExogVar[170][comptExogVar[170]]) * (t - exogSamplingTime[170][comptExogVar[170]])/(exogSamplingTime[170][comptExogVar[170] + 
    1] - exogSamplingTime[170][comptExogVar[170]]);
x[171] = dataExogVar[171][comptExogVar[171]] + (dataExogVar[171][comptExogVar[171] + 1] - dataExogVar[171][comptExogVar[171]]) * (t - exogSamplingTime[171][comptExogVar[171]])/(exogSamplingTime[171][comptExogVar[171] + 
    1] - exogSamplingTime[171][comptExogVar[171]]);
x[172] = dataExogVar[172][comptExogVar[172]] + (dataExogVar[172][comptExogVar[172] + 1] - dataExogVar[172][comptExogVar[172]]) * (t - exogSamplingTime[172][comptExogVar[172]])/(exogSamplingTime[172][comptExogVar[172] + 
    1] - exogSamplingTime[172][comptExogVar[172]]);
x[173] = dataExogVar[173][comptExogVar[173]] + (dataExogVar[173][comptExogVar[173] + 1] - dataExogVar[173][comptExogVar[173]]) * (t - exogSamplingTime[173][comptExogVar[173]])/(exogSamplingTime[173][comptExogVar[173] + 
    1] - exogSamplingTime[173][comptExogVar[173]]);
x[174] = dataExogVar[174][comptExogVar[174]] + (dataExogVar[174][comptExogVar[174] + 1] - dataExogVar[174][comptExogVar[174]]) * (t - exogSamplingTime[174][comptExogVar[174]])/(exogSamplingTime[174][comptExogVar[174] + 
    1] - exogSamplingTime[174][comptExogVar[174]]);
x[175] = dataExogVar[175][comptExogVar[175]] + (dataExogVar[175][comptExogVar[175] + 1] - dataExogVar[175][comptExogVar[175]]) * (t - exogSamplingTime[175][comptExogVar[175]])/(exogSamplingTime[175][comptExogVar[175] + 
    1] - exogSamplingTime[175][comptExogVar[175]]);
x[176] = dataExogVar[176][comptExogVar[176]] + (dataExogVar[176][comptExogVar[176] + 1] - dataExogVar[176][comptExogVar[176]]) * (t - exogSamplingTime[176][comptExogVar[176]])/(exogSamplingTime[176][comptExogVar[176] + 
    1] - exogSamplingTime[176][comptExogVar[176]]);
x[177] = dataExogVar[177][comptExogVar[177]] + (dataExogVar[177][comptExogVar[177] + 1] - dataExogVar[177][comptExogVar[177]]) * (t - exogSamplingTime[177][comptExogVar[177]])/(exogSamplingTime[177][comptExogVar[177] + 
    1] - exogSamplingTime[177][comptExogVar[177]]);
x[178] = dataExogVar[178][comptExogVar[178]] + (dataExogVar[178][comptExogVar[178] + 1] - dataExogVar[178][comptExogVar[178]]) * (t - exogSamplingTime[178][comptExogVar[178]])/(exogSamplingTime[178][comptExogVar[178] + 
    1] - exogSamplingTime[178][comptExogVar[178]]);
x[179] = dataExogVar[179][comptExogVar[179]] + (dataExogVar[179][comptExogVar[179] + 1] - dataExogVar[179][comptExogVar[179]]) * (t - exogSamplingTime[179][comptExogVar[179]])/(exogSamplingTime[179][comptExogVar[179] + 
    1] - exogSamplingTime[179][comptExogVar[179]]);
x[180] = dataExogVar[180][comptExogVar[180]] + (dataExogVar[180][comptExogVar[180] + 1] - dataExogVar[180][comptExogVar[180]]) * (t - exogSamplingTime[180][comptExogVar[180]])/(exogSamplingTime[180][comptExogVar[180] + 
    1] - exogSamplingTime[180][comptExogVar[180]]);
x[181] = dataExogVar[181][comptExogVar[181]] + (dataExogVar[181][comptExogVar[181] + 1] - dataExogVar[181][comptExogVar[181]]) * (t - exogSamplingTime[181][comptExogVar[181]])/(exogSamplingTime[181][comptExogVar[181] + 
    1] - exogSamplingTime[181][comptExogVar[181]]);
x[182] = dataExogVar[182][comptExogVar[182]] + (dataExogVar[182][comptExogVar[182] + 1] - dataExogVar[182][comptExogVar[182]]) * (t - exogSamplingTime[182][comptExogVar[182]])/(exogSamplingTime[182][comptExogVar[182] + 
    1] - exogSamplingTime[182][comptExogVar[182]]);
x[183] = dataExogVar[183][comptExogVar[183]] + (dataExogVar[183][comptExogVar[183] + 1] - dataExogVar[183][comptExogVar[183]]) * (t - exogSamplingTime[183][comptExogVar[183]])/(exogSamplingTime[183][comptExogVar[183] + 
    1] - exogSamplingTime[183][comptExogVar[183]]);
x[184] = dataExogVar[184][comptExogVar[184]] + (dataExogVar[184][comptExogVar[184] + 1] - dataExogVar[184][comptExogVar[184]]) * (t - exogSamplingTime[184][comptExogVar[184]])/(exogSamplingTime[184][comptExogVar[184] + 
    1] - exogSamplingTime[184][comptExogVar[184]]);
x[185] = dataExogVar[185][comptExogVar[185]] + (dataExogVar[185][comptExogVar[185] + 1] - dataExogVar[185][comptExogVar[185]]) * (t - exogSamplingTime[185][comptExogVar[185]])/(exogSamplingTime[185][comptExogVar[185] + 
    1] - exogSamplingTime[185][comptExogVar[185]]);
x[186] = dataExogVar[186][comptExogVar[186]] + (dataExogVar[186][comptExogVar[186] + 1] - dataExogVar[186][comptExogVar[186]]) * (t - exogSamplingTime[186][comptExogVar[186]])/(exogSamplingTime[186][comptExogVar[186] + 
    1] - exogSamplingTime[186][comptExogVar[186]]);
x[187] = dataExogVar[187][comptExogVar[187]] + (dataExogVar[187][comptExogVar[187] + 1] - dataExogVar[187][comptExogVar[187]]) * (t - exogSamplingTime[187][comptExogVar[187]])/(exogSamplingTime[187][comptExogVar[187] + 
    1] - exogSamplingTime[187][comptExogVar[187]]);
x[188] = dataExogVar[188][comptExogVar[188]] + (dataExogVar[188][comptExogVar[188] + 1] - dataExogVar[188][comptExogVar[188]]) * (t - exogSamplingTime[188][comptExogVar[188]])/(exogSamplingTime[188][comptExogVar[188] + 
    1] - exogSamplingTime[188][comptExogVar[188]]);
x[189] = dataExogVar[189][comptExogVar[189]] + (dataExogVar[189][comptExogVar[189] + 1] - dataExogVar[189][comptExogVar[189]]) * (t - exogSamplingTime[189][comptExogVar[189]])/(exogSamplingTime[189][comptExogVar[189] + 
    1] - exogSamplingTime[189][comptExogVar[189]]);
x[190] = dataExogVar[190][comptExogVar[190]] + (dataExogVar[190][comptExogVar[190] + 1] - dataExogVar[190][comptExogVar[190]]) * (t - exogSamplingTime[190][comptExogVar[190]])/(exogSamplingTime[190][comptExogVar[190] + 
    1] - exogSamplingTime[190][comptExogVar[190]]);
x[191] = dataExogVar[191][comptExogVar[191]] + (dataExogVar[191][comptExogVar[191] + 1] - dataExogVar[191][comptExogVar[191]]) * (t - exogSamplingTime[191][comptExogVar[191]])/(exogSamplingTime[191][comptExogVar[191] + 
    1] - exogSamplingTime[191][comptExogVar[191]]);
x[192] = dataExogVar[192][comptExogVar[192]] + (dataExogVar[192][comptExogVar[192] + 1] - dataExogVar[192][comptExogVar[192]]) * (t - exogSamplingTime[192][comptExogVar[192]])/(exogSamplingTime[192][comptExogVar[192] + 
    1] - exogSamplingTime[192][comptExogVar[192]]);
x[193] = dataExogVar[193][comptExogVar[193]] + (dataExogVar[193][comptExogVar[193] + 1] - dataExogVar[193][comptExogVar[193]]) * (t - exogSamplingTime[193][comptExogVar[193]])/(exogSamplingTime[193][comptExogVar[193] + 
    1] - exogSamplingTime[193][comptExogVar[193]]);
x[194] = dataExogVar[194][comptExogVar[194]] + (dataExogVar[194][comptExogVar[194] + 1] - dataExogVar[194][comptExogVar[194]]) * (t - exogSamplingTime[194][comptExogVar[194]])/(exogSamplingTime[194][comptExogVar[194] + 
    1] - exogSamplingTime[194][comptExogVar[194]]);
x[195] = dataExogVar[195][comptExogVar[195]] + (dataExogVar[195][comptExogVar[195] + 1] - dataExogVar[195][comptExogVar[195]]) * (t - exogSamplingTime[195][comptExogVar[195]])/(exogSamplingTime[195][comptExogVar[195] + 
    1] - exogSamplingTime[195][comptExogVar[195]]);
x[196] = dataExogVar[196][comptExogVar[196]] + (dataExogVar[196][comptExogVar[196] + 1] - dataExogVar[196][comptExogVar[196]]) * (t - exogSamplingTime[196][comptExogVar[196]])/(exogSamplingTime[196][comptExogVar[196] + 
    1] - exogSamplingTime[196][comptExogVar[196]]);
x[197] = dataExogVar[197][comptExogVar[197]] + (dataExogVar[197][comptExogVar[197] + 1] - dataExogVar[197][comptExogVar[197]]) * (t - exogSamplingTime[197][comptExogVar[197]])/(exogSamplingTime[197][comptExogVar[197] + 
    1] - exogSamplingTime[197][comptExogVar[197]]);
x[198] = dataExogVar[198][comptExogVar[198]] + (dataExogVar[198][comptExogVar[198] + 1] - dataExogVar[198][comptExogVar[198]]) * (t - exogSamplingTime[198][comptExogVar[198]])/(exogSamplingTime[198][comptExogVar[198] + 
    1] - exogSamplingTime[198][comptExogVar[198]]);
x[199] = dataExogVar[199][comptExogVar[199]] + (dataExogVar[199][comptExogVar[199] + 1] - dataExogVar[199][comptExogVar[199]]) * (t - exogSamplingTime[199][comptExogVar[199]])/(exogSamplingTime[199][comptExogVar[199] + 
    1] - exogSamplingTime[199][comptExogVar[199]]);
x[200] = dataExogVar[200][comptExogVar[200]] + (dataExogVar[200][comptExogVar[200] + 1] - dataExogVar[200][comptExogVar[200]]) * (t - exogSamplingTime[200][comptExogVar[200]])/(exogSamplingTime[200][comptExogVar[200] + 
    1] - exogSamplingTime[200][comptExogVar[200]]);
x[201] = dataExogVar[201][comptExogVar[201]] + (dataExogVar[201][comptExogVar[201] + 1] - dataExogVar[201][comptExogVar[201]]) * (t - exogSamplingTime[201][comptExogVar[201]])/(exogSamplingTime[201][comptExogVar[201] + 
    1] - exogSamplingTime[201][comptExogVar[201]]);
x[202] = dataExogVar[202][comptExogVar[202]] + (dataExogVar[202][comptExogVar[202] + 1] - dataExogVar[202][comptExogVar[202]]) * (t - exogSamplingTime[202][comptExogVar[202]])/(exogSamplingTime[202][comptExogVar[202] + 
    1] - exogSamplingTime[202][comptExogVar[202]]);
x[203] = dataExogVar[203][comptExogVar[203]] + (dataExogVar[203][comptExogVar[203] + 1] - dataExogVar[203][comptExogVar[203]]) * (t - exogSamplingTime[203][comptExogVar[203]])/(exogSamplingTime[203][comptExogVar[203] + 
    1] - exogSamplingTime[203][comptExogVar[203]]);
x[204] = dataExogVar[204][comptExogVar[204]] + (dataExogVar[204][comptExogVar[204] + 1] - dataExogVar[204][comptExogVar[204]]) * (t - exogSamplingTime[204][comptExogVar[204]])/(exogSamplingTime[204][comptExogVar[204] + 
    1] - exogSamplingTime[204][comptExogVar[204]]);
x[205] = dataExogVar[205][comptExogVar[205]] + (dataExogVar[205][comptExogVar[205] + 1] - dataExogVar[205][comptExogVar[205]]) * (t - exogSamplingTime[205][comptExogVar[205]])/(exogSamplingTime[205][comptExogVar[205] + 
    1] - exogSamplingTime[205][comptExogVar[205]]);
x[206] = dataExogVar[206][comptExogVar[206]] + (dataExogVar[206][comptExogVar[206] + 1] - dataExogVar[206][comptExogVar[206]]) * (t - exogSamplingTime[206][comptExogVar[206]])/(exogSamplingTime[206][comptExogVar[206] + 
    1] - exogSamplingTime[206][comptExogVar[206]]);
x[207] = dataExogVar[207][comptExogVar[207]] + (dataExogVar[207][comptExogVar[207] + 1] - dataExogVar[207][comptExogVar[207]]) * (t - exogSamplingTime[207][comptExogVar[207]])/(exogSamplingTime[207][comptExogVar[207] + 
    1] - exogSamplingTime[207][comptExogVar[207]]);
x[208] = dataExogVar[208][comptExogVar[208]] + (dataExogVar[208][comptExogVar[208] + 1] - dataExogVar[208][comptExogVar[208]]) * (t - exogSamplingTime[208][comptExogVar[208]])/(exogSamplingTime[208][comptExogVar[208] + 
    1] - exogSamplingTime[208][comptExogVar[208]]);
x[209] = dataExogVar[209][comptExogVar[209]] + (dataExogVar[209][comptExogVar[209] + 1] - dataExogVar[209][comptExogVar[209]]) * (t - exogSamplingTime[209][comptExogVar[209]])/(exogSamplingTime[209][comptExogVar[209] + 
    1] - exogSamplingTime[209][comptExogVar[209]]);
x[210] = dataExogVar[210][comptExogVar[210]] + (dataExogVar[210][comptExogVar[210] + 1] - dataExogVar[210][comptExogVar[210]]) * (t - exogSamplingTime[210][comptExogVar[210]])/(exogSamplingTime[210][comptExogVar[210] + 
    1] - exogSamplingTime[210][comptExogVar[210]]);
x[211] = dataExogVar[211][comptExogVar[211]] + (dataExogVar[211][comptExogVar[211] + 1] - dataExogVar[211][comptExogVar[211]]) * (t - exogSamplingTime[211][comptExogVar[211]])/(exogSamplingTime[211][comptExogVar[211] + 
    1] - exogSamplingTime[211][comptExogVar[211]]);
x[212] = dataExogVar[212][comptExogVar[212]] + (dataExogVar[212][comptExogVar[212] + 1] - dataExogVar[212][comptExogVar[212]]) * (t - exogSamplingTime[212][comptExogVar[212]])/(exogSamplingTime[212][comptExogVar[212] + 
    1] - exogSamplingTime[212][comptExogVar[212]]);
x[213] = dataExogVar[213][comptExogVar[213]] + (dataExogVar[213][comptExogVar[213] + 1] - dataExogVar[213][comptExogVar[213]]) * (t - exogSamplingTime[213][comptExogVar[213]])/(exogSamplingTime[213][comptExogVar[213] + 
    1] - exogSamplingTime[213][comptExogVar[213]]);
x[214] = dataExogVar[214][comptExogVar[214]] + (dataExogVar[214][comptExogVar[214] + 1] - dataExogVar[214][comptExogVar[214]]) * (t - exogSamplingTime[214][comptExogVar[214]])/(exogSamplingTime[214][comptExogVar[214] + 
    1] - exogSamplingTime[214][comptExogVar[214]]);
x[215] = dataExogVar[215][comptExogVar[215]] + (dataExogVar[215][comptExogVar[215] + 1] - dataExogVar[215][comptExogVar[215]]) * (t - exogSamplingTime[215][comptExogVar[215]])/(exogSamplingTime[215][comptExogVar[215] + 
    1] - exogSamplingTime[215][comptExogVar[215]]);
x[216] = dataExogVar[216][comptExogVar[216]] + (dataExogVar[216][comptExogVar[216] + 1] - dataExogVar[216][comptExogVar[216]]) * (t - exogSamplingTime[216][comptExogVar[216]])/(exogSamplingTime[216][comptExogVar[216] + 
    1] - exogSamplingTime[216][comptExogVar[216]]);
x[217] = dataExogVar[217][comptExogVar[217]] + (dataExogVar[217][comptExogVar[217] + 1] - dataExogVar[217][comptExogVar[217]]) * (t - exogSamplingTime[217][comptExogVar[217]])/(exogSamplingTime[217][comptExogVar[217] + 
    1] - exogSamplingTime[217][comptExogVar[217]]);
x[218] = dataExogVar[218][comptExogVar[218]] + (dataExogVar[218][comptExogVar[218] + 1] - dataExogVar[218][comptExogVar[218]]) * (t - exogSamplingTime[218][comptExogVar[218]])/(exogSamplingTime[218][comptExogVar[218] + 
    1] - exogSamplingTime[218][comptExogVar[218]]);
x[219] = dataExogVar[219][comptExogVar[219]] + (dataExogVar[219][comptExogVar[219] + 1] - dataExogVar[219][comptExogVar[219]]) * (t - exogSamplingTime[219][comptExogVar[219]])/(exogSamplingTime[219][comptExogVar[219] + 
    1] - exogSamplingTime[219][comptExogVar[219]]);
x[220] = dataExogVar[220][comptExogVar[220]] + (dataExogVar[220][comptExogVar[220] + 1] - dataExogVar[220][comptExogVar[220]]) * (t - exogSamplingTime[220][comptExogVar[220]])/(exogSamplingTime[220][comptExogVar[220] + 
    1] - exogSamplingTime[220][comptExogVar[220]]);
x[221] = dataExogVar[221][comptExogVar[221]] + (dataExogVar[221][comptExogVar[221] + 1] - dataExogVar[221][comptExogVar[221]]) * (t - exogSamplingTime[221][comptExogVar[221]])/(exogSamplingTime[221][comptExogVar[221] + 
    1] - exogSamplingTime[221][comptExogVar[221]]);
x[222] = dataExogVar[222][comptExogVar[222]] + (dataExogVar[222][comptExogVar[222] + 1] - dataExogVar[222][comptExogVar[222]]) * (t - exogSamplingTime[222][comptExogVar[222]])/(exogSamplingTime[222][comptExogVar[222] + 
    1] - exogSamplingTime[222][comptExogVar[222]]);
x[223] = dataExogVar[223][comptExogVar[223]] + (dataExogVar[223][comptExogVar[223] + 1] - dataExogVar[223][comptExogVar[223]]) * (t - exogSamplingTime[223][comptExogVar[223]])/(exogSamplingTime[223][comptExogVar[223] + 
    1] - exogSamplingTime[223][comptExogVar[223]]);
x[224] = dataExogVar[224][comptExogVar[224]] + (dataExogVar[224][comptExogVar[224] + 1] - dataExogVar[224][comptExogVar[224]]) * (t - exogSamplingTime[224][comptExogVar[224]])/(exogSamplingTime[224][comptExogVar[224] + 
    1] - exogSamplingTime[224][comptExogVar[224]]);
x[225] = dataExogVar[225][comptExogVar[225]] + (dataExogVar[225][comptExogVar[225] + 1] - dataExogVar[225][comptExogVar[225]]) * (t - exogSamplingTime[225][comptExogVar[225]])/(exogSamplingTime[225][comptExogVar[225] + 
    1] - exogSamplingTime[225][comptExogVar[225]]);
x[226] = dataExogVar[226][comptExogVar[226]] + (dataExogVar[226][comptExogVar[226] + 1] - dataExogVar[226][comptExogVar[226]]) * (t - exogSamplingTime[226][comptExogVar[226]])/(exogSamplingTime[226][comptExogVar[226] + 
    1] - exogSamplingTime[226][comptExogVar[226]]);
x[227] = dataExogVar[227][comptExogVar[227]] + (dataExogVar[227][comptExogVar[227] + 1] - dataExogVar[227][comptExogVar[227]]) * (t - exogSamplingTime[227][comptExogVar[227]])/(exogSamplingTime[227][comptExogVar[227] + 
    1] - exogSamplingTime[227][comptExogVar[227]]);
x[228] = dataExogVar[228][comptExogVar[228]] + (dataExogVar[228][comptExogVar[228] + 1] - dataExogVar[228][comptExogVar[228]]) * (t - exogSamplingTime[228][comptExogVar[228]])/(exogSamplingTime[228][comptExogVar[228] + 
    1] - exogSamplingTime[228][comptExogVar[228]]);
x[229] = dataExogVar[229][comptExogVar[229]] + (dataExogVar[229][comptExogVar[229] + 1] - dataExogVar[229][comptExogVar[229]]) * (t - exogSamplingTime[229][comptExogVar[229]])/(exogSamplingTime[229][comptExogVar[229] + 
    1] - exogSamplingTime[229][comptExogVar[229]]);
x[230] = dataExogVar[230][comptExogVar[230]] + (dataExogVar[230][comptExogVar[230] + 1] - dataExogVar[230][comptExogVar[230]]) * (t - exogSamplingTime[230][comptExogVar[230]])/(exogSamplingTime[230][comptExogVar[230] + 
    1] - exogSamplingTime[230][comptExogVar[230]]);
x[231] = dataExogVar[231][comptExogVar[231]] + (dataExogVar[231][comptExogVar[231] + 1] - dataExogVar[231][comptExogVar[231]]) * (t - exogSamplingTime[231][comptExogVar[231]])/(exogSamplingTime[231][comptExogVar[231] + 
    1] - exogSamplingTime[231][comptExogVar[231]]);
x[232] = dataExogVar[232][comptExogVar[232]] + (dataExogVar[232][comptExogVar[232] + 1] - dataExogVar[232][comptExogVar[232]]) * (t - exogSamplingTime[232][comptExogVar[232]])/(exogSamplingTime[232][comptExogVar[232] + 
    1] - exogSamplingTime[232][comptExogVar[232]]);
x[233] = dataExogVar[233][comptExogVar[233]] + (dataExogVar[233][comptExogVar[233] + 1] - dataExogVar[233][comptExogVar[233]]) * (t - exogSamplingTime[233][comptExogVar[233]])/(exogSamplingTime[233][comptExogVar[233] + 
    1] - exogSamplingTime[233][comptExogVar[233]]);
x[234] = dataExogVar[234][comptExogVar[234]] + (dataExogVar[234][comptExogVar[234] + 1] - dataExogVar[234][comptExogVar[234]]) * (t - exogSamplingTime[234][comptExogVar[234]])/(exogSamplingTime[234][comptExogVar[234] + 
    1] - exogSamplingTime[234][comptExogVar[234]]);
x[235] = dataExogVar[235][comptExogVar[235]] + (dataExogVar[235][comptExogVar[235] + 1] - dataExogVar[235][comptExogVar[235]]) * (t - exogSamplingTime[235][comptExogVar[235]])/(exogSamplingTime[235][comptExogVar[235] + 
    1] - exogSamplingTime[235][comptExogVar[235]]);
x[236] = dataExogVar[236][comptExogVar[236]] + (dataExogVar[236][comptExogVar[236] + 1] - dataExogVar[236][comptExogVar[236]]) * (t - exogSamplingTime[236][comptExogVar[236]])/(exogSamplingTime[236][comptExogVar[236] + 
    1] - exogSamplingTime[236][comptExogVar[236]]);
x[237] = dataExogVar[237][comptExogVar[237]] + (dataExogVar[237][comptExogVar[237] + 1] - dataExogVar[237][comptExogVar[237]]) * (t - exogSamplingTime[237][comptExogVar[237]])/(exogSamplingTime[237][comptExogVar[237] + 
    1] - exogSamplingTime[237][comptExogVar[237]]);
x[238] = dataExogVar[238][comptExogVar[238]] + (dataExogVar[238][comptExogVar[238] + 1] - dataExogVar[238][comptExogVar[238]]) * (t - exogSamplingTime[238][comptExogVar[238]])/(exogSamplingTime[238][comptExogVar[238] + 
    1] - exogSamplingTime[238][comptExogVar[238]]);
x[239] = dataExogVar[239][comptExogVar[239]] + (dataExogVar[239][comptExogVar[239] + 1] - dataExogVar[239][comptExogVar[239]]) * (t - exogSamplingTime[239][comptExogVar[239]])/(exogSamplingTime[239][comptExogVar[239] + 
    1] - exogSamplingTime[239][comptExogVar[239]]);
x[240] = dataExogVar[240][comptExogVar[240]] + (dataExogVar[240][comptExogVar[240] + 1] - dataExogVar[240][comptExogVar[240]]) * (t - exogSamplingTime[240][comptExogVar[240]])/(exogSamplingTime[240][comptExogVar[240] + 
    1] - exogSamplingTime[240][comptExogVar[240]]);
x[241] = dataExogVar[241][comptExogVar[241]] + (dataExogVar[241][comptExogVar[241] + 1] - dataExogVar[241][comptExogVar[241]]) * (t - exogSamplingTime[241][comptExogVar[241]])/(exogSamplingTime[241][comptExogVar[241] + 
    1] - exogSamplingTime[241][comptExogVar[241]]);
x[242] = dataExogVar[242][comptExogVar[242]] + (dataExogVar[242][comptExogVar[242] + 1] - dataExogVar[242][comptExogVar[242]]) * (t - exogSamplingTime[242][comptExogVar[242]])/(exogSamplingTime[242][comptExogVar[242] + 
    1] - exogSamplingTime[242][comptExogVar[242]]);
x[243] = dataExogVar[243][comptExogVar[243]] + (dataExogVar[243][comptExogVar[243] + 1] - dataExogVar[243][comptExogVar[243]]) * (t - exogSamplingTime[243][comptExogVar[243]])/(exogSamplingTime[243][comptExogVar[243] + 
    1] - exogSamplingTime[243][comptExogVar[243]]);
x[244] = dataExogVar[244][comptExogVar[244]] + (dataExogVar[244][comptExogVar[244] + 1] - dataExogVar[244][comptExogVar[244]]) * (t - exogSamplingTime[244][comptExogVar[244]])/(exogSamplingTime[244][comptExogVar[244] + 
    1] - exogSamplingTime[244][comptExogVar[244]]);
x[245] = dataExogVar[245][comptExogVar[245]] + (dataExogVar[245][comptExogVar[245] + 1] - dataExogVar[245][comptExogVar[245]]) * (t - exogSamplingTime[245][comptExogVar[245]])/(exogSamplingTime[245][comptExogVar[245] + 
    1] - exogSamplingTime[245][comptExogVar[245]]);
x[246] = dataExogVar[246][comptExogVar[246]] + (dataExogVar[246][comptExogVar[246] + 1] - dataExogVar[246][comptExogVar[246]]) * (t - exogSamplingTime[246][comptExogVar[246]])/(exogSamplingTime[246][comptExogVar[246] + 
    1] - exogSamplingTime[246][comptExogVar[246]]);
x[247] = dataExogVar[247][comptExogVar[247]] + (dataExogVar[247][comptExogVar[247] + 1] - dataExogVar[247][comptExogVar[247]]) * (t - exogSamplingTime[247][comptExogVar[247]])/(exogSamplingTime[247][comptExogVar[247] + 
    1] - exogSamplingTime[247][comptExogVar[247]]);
x[248] = dataExogVar[248][comptExogVar[248]] + (dataExogVar[248][comptExogVar[248] + 1] - dataExogVar[248][comptExogVar[248]]) * (t - exogSamplingTime[248][comptExogVar[248]])/(exogSamplingTime[248][comptExogVar[248] + 
    1] - exogSamplingTime[248][comptExogVar[248]]);
x[249] = dataExogVar[249][comptExogVar[249]] + (dataExogVar[249][comptExogVar[249] + 1] - dataExogVar[249][comptExogVar[249]]) * (t - exogSamplingTime[249][comptExogVar[249]])/(exogSamplingTime[249][comptExogVar[249] + 
    1] - exogSamplingTime[249][comptExogVar[249]]);
x[250] = dataExogVar[250][comptExogVar[250]] + (dataExogVar[250][comptExogVar[250] + 1] - dataExogVar[250][comptExogVar[250]]) * (t - exogSamplingTime[250][comptExogVar[250]])/(exogSamplingTime[250][comptExogVar[250] + 
    1] - exogSamplingTime[250][comptExogVar[250]]);
x[251] = dataExogVar[251][comptExogVar[251]] + (dataExogVar[251][comptExogVar[251] + 1] - dataExogVar[251][comptExogVar[251]]) * (t - exogSamplingTime[251][comptExogVar[251]])/(exogSamplingTime[251][comptExogVar[251] + 
    1] - exogSamplingTime[251][comptExogVar[251]]);
x[252] = dataExogVar[252][comptExogVar[252]] + (dataExogVar[252][comptExogVar[252] + 1] - dataExogVar[252][comptExogVar[252]]) * (t - exogSamplingTime[252][comptExogVar[252]])/(exogSamplingTime[252][comptExogVar[252] + 
    1] - exogSamplingTime[252][comptExogVar[252]]);
x[253] = dataExogVar[253][comptExogVar[253]] + (dataExogVar[253][comptExogVar[253] + 1] - dataExogVar[253][comptExogVar[253]]) * (t - exogSamplingTime[253][comptExogVar[253]])/(exogSamplingTime[253][comptExogVar[253] + 
    1] - exogSamplingTime[253][comptExogVar[253]]);
x[254] = dataExogVar[254][comptExogVar[254]] + (dataExogVar[254][comptExogVar[254] + 1] - dataExogVar[254][comptExogVar[254]]) * (t - exogSamplingTime[254][comptExogVar[254]])/(exogSamplingTime[254][comptExogVar[254] + 
    1] - exogSamplingTime[254][comptExogVar[254]]);
x[255] = dataExogVar[255][comptExogVar[255]] + (dataExogVar[255][comptExogVar[255] + 1] - dataExogVar[255][comptExogVar[255]]) * (t - exogSamplingTime[255][comptExogVar[255]])/(exogSamplingTime[255][comptExogVar[255] + 
    1] - exogSamplingTime[255][comptExogVar[255]]);
x[256] = dataExogVar[256][comptExogVar[256]] + (dataExogVar[256][comptExogVar[256] + 1] - dataExogVar[256][comptExogVar[256]]) * (t - exogSamplingTime[256][comptExogVar[256]])/(exogSamplingTime[256][comptExogVar[256] + 
    1] - exogSamplingTime[256][comptExogVar[256]]);
x[257] = dataExogVar[257][comptExogVar[257]] + (dataExogVar[257][comptExogVar[257] + 1] - dataExogVar[257][comptExogVar[257]]) * (t - exogSamplingTime[257][comptExogVar[257]])/(exogSamplingTime[257][comptExogVar[257] + 
    1] - exogSamplingTime[257][comptExogVar[257]]);
x[258] = dataExogVar[258][comptExogVar[258]] + (dataExogVar[258][comptExogVar[258] + 1] - dataExogVar[258][comptExogVar[258]]) * (t - exogSamplingTime[258][comptExogVar[258]])/(exogSamplingTime[258][comptExogVar[258] + 
    1] - exogSamplingTime[258][comptExogVar[258]]);
x[259] = dataExogVar[259][comptExogVar[259]] + (dataExogVar[259][comptExogVar[259] + 1] - dataExogVar[259][comptExogVar[259]]) * (t - exogSamplingTime[259][comptExogVar[259]])/(exogSamplingTime[259][comptExogVar[259] + 
    1] - exogSamplingTime[259][comptExogVar[259]]);
x[260] = dataExogVar[260][comptExogVar[260]] + (dataExogVar[260][comptExogVar[260] + 1] - dataExogVar[260][comptExogVar[260]]) * (t - exogSamplingTime[260][comptExogVar[260]])/(exogSamplingTime[260][comptExogVar[260] + 
    1] - exogSamplingTime[260][comptExogVar[260]]);
x[261] = dataExogVar[261][comptExogVar[261]] + (dataExogVar[261][comptExogVar[261] + 1] - dataExogVar[261][comptExogVar[261]]) * (t - exogSamplingTime[261][comptExogVar[261]])/(exogSamplingTime[261][comptExogVar[261] + 
    1] - exogSamplingTime[261][comptExogVar[261]]);
x[262] = dataExogVar[262][comptExogVar[262]] + (dataExogVar[262][comptExogVar[262] + 1] - dataExogVar[262][comptExogVar[262]]) * (t - exogSamplingTime[262][comptExogVar[262]])/(exogSamplingTime[262][comptExogVar[262] + 
    1] - exogSamplingTime[262][comptExogVar[262]]);
x[263] = dataExogVar[263][comptExogVar[263]] + (dataExogVar[263][comptExogVar[263] + 1] - dataExogVar[263][comptExogVar[263]]) * (t - exogSamplingTime[263][comptExogVar[263]])/(exogSamplingTime[263][comptExogVar[263] + 
    1] - exogSamplingTime[263][comptExogVar[263]]);
x[264] = dataExogVar[264][comptExogVar[264]] + (dataExogVar[264][comptExogVar[264] + 1] - dataExogVar[264][comptExogVar[264]]) * (t - exogSamplingTime[264][comptExogVar[264]])/(exogSamplingTime[264][comptExogVar[264] + 
    1] - exogSamplingTime[264][comptExogVar[264]]);
x[265] = dataExogVar[265][comptExogVar[265]] + (dataExogVar[265][comptExogVar[265] + 1] - dataExogVar[265][comptExogVar[265]]) * (t - exogSamplingTime[265][comptExogVar[265]])/(exogSamplingTime[265][comptExogVar[265] + 
    1] - exogSamplingTime[265][comptExogVar[265]]);
x[266] = dataExogVar[266][comptExogVar[266]] + (dataExogVar[266][comptExogVar[266] + 1] - dataExogVar[266][comptExogVar[266]]) * (t - exogSamplingTime[266][comptExogVar[266]])/(exogSamplingTime[266][comptExogVar[266] + 
    1] - exogSamplingTime[266][comptExogVar[266]]);
x[267] = dataExogVar[267][comptExogVar[267]] + (dataExogVar[267][comptExogVar[267] + 1] - dataExogVar[267][comptExogVar[267]]) * (t - exogSamplingTime[267][comptExogVar[267]])/(exogSamplingTime[267][comptExogVar[267] + 
    1] - exogSamplingTime[267][comptExogVar[267]]);
x[268] = dataExogVar[268][comptExogVar[268]] + (dataExogVar[268][comptExogVar[268] + 1] - dataExogVar[268][comptExogVar[268]]) * (t - exogSamplingTime[268][comptExogVar[268]])/(exogSamplingTime[268][comptExogVar[268] + 
    1] - exogSamplingTime[268][comptExogVar[268]]);
x[269] = dataExogVar[269][comptExogVar[269]] + (dataExogVar[269][comptExogVar[269] + 1] - dataExogVar[269][comptExogVar[269]]) * (t - exogSamplingTime[269][comptExogVar[269]])/(exogSamplingTime[269][comptExogVar[269] + 
    1] - exogSamplingTime[269][comptExogVar[269]]);
x[270] = dataExogVar[270][comptExogVar[270]] + (dataExogVar[270][comptExogVar[270] + 1] - dataExogVar[270][comptExogVar[270]]) * (t - exogSamplingTime[270][comptExogVar[270]])/(exogSamplingTime[270][comptExogVar[270] + 
    1] - exogSamplingTime[270][comptExogVar[270]]);
x[271] = dataExogVar[271][comptExogVar[271]] + (dataExogVar[271][comptExogVar[271] + 1] - dataExogVar[271][comptExogVar[271]]) * (t - exogSamplingTime[271][comptExogVar[271]])/(exogSamplingTime[271][comptExogVar[271] + 
    1] - exogSamplingTime[271][comptExogVar[271]]);
x[272] = dataExogVar[272][comptExogVar[272]] + (dataExogVar[272][comptExogVar[272] + 1] - dataExogVar[272][comptExogVar[272]]) * (t - exogSamplingTime[272][comptExogVar[272]])/(exogSamplingTime[272][comptExogVar[272] + 
    1] - exogSamplingTime[272][comptExogVar[272]]);
x[273] = dataExogVar[273][comptExogVar[273]] + (dataExogVar[273][comptExogVar[273] + 1] - dataExogVar[273][comptExogVar[273]]) * (t - exogSamplingTime[273][comptExogVar[273]])/(exogSamplingTime[273][comptExogVar[273] + 
    1] - exogSamplingTime[273][comptExogVar[273]]);
x[274] = dataExogVar[274][comptExogVar[274]] + (dataExogVar[274][comptExogVar[274] + 1] - dataExogVar[274][comptExogVar[274]]) * (t - exogSamplingTime[274][comptExogVar[274]])/(exogSamplingTime[274][comptExogVar[274] + 
    1] - exogSamplingTime[274][comptExogVar[274]]);
x[275] = dataExogVar[275][comptExogVar[275]] + (dataExogVar[275][comptExogVar[275] + 1] - dataExogVar[275][comptExogVar[275]]) * (t - exogSamplingTime[275][comptExogVar[275]])/(exogSamplingTime[275][comptExogVar[275] + 
    1] - exogSamplingTime[275][comptExogVar[275]]);
x[276] = dataExogVar[276][comptExogVar[276]] + (dataExogVar[276][comptExogVar[276] + 1] - dataExogVar[276][comptExogVar[276]]) * (t - exogSamplingTime[276][comptExogVar[276]])/(exogSamplingTime[276][comptExogVar[276] + 
    1] - exogSamplingTime[276][comptExogVar[276]]);
x[277] = dataExogVar[277][comptExogVar[277]] + (dataExogVar[277][comptExogVar[277] + 1] - dataExogVar[277][comptExogVar[277]]) * (t - exogSamplingTime[277][comptExogVar[277]])/(exogSamplingTime[277][comptExogVar[277] + 
    1] - exogSamplingTime[277][comptExogVar[277]]);
x[278] = dataExogVar[278][comptExogVar[278]] + (dataExogVar[278][comptExogVar[278] + 1] - dataExogVar[278][comptExogVar[278]]) * (t - exogSamplingTime[278][comptExogVar[278]])/(exogSamplingTime[278][comptExogVar[278] + 
    1] - exogSamplingTime[278][comptExogVar[278]]);
x[279] = dataExogVar[279][comptExogVar[279]] + (dataExogVar[279][comptExogVar[279] + 1] - dataExogVar[279][comptExogVar[279]]) * (t - exogSamplingTime[279][comptExogVar[279]])/(exogSamplingTime[279][comptExogVar[279] + 
    1] - exogSamplingTime[279][comptExogVar[279]]);
x[280] = dataExogVar[280][comptExogVar[280]] + (dataExogVar[280][comptExogVar[280] + 1] - dataExogVar[280][comptExogVar[280]]) * (t - exogSamplingTime[280][comptExogVar[280]])/(exogSamplingTime[280][comptExogVar[280] + 
    1] - exogSamplingTime[280][comptExogVar[280]]);
x[281] = dataExogVar[281][comptExogVar[281]] + (dataExogVar[281][comptExogVar[281] + 1] - dataExogVar[281][comptExogVar[281]]) * (t - exogSamplingTime[281][comptExogVar[281]])/(exogSamplingTime[281][comptExogVar[281] + 
    1] - exogSamplingTime[281][comptExogVar[281]]);
x[282] = dataExogVar[282][comptExogVar[282]] + (dataExogVar[282][comptExogVar[282] + 1] - dataExogVar[282][comptExogVar[282]]) * (t - exogSamplingTime[282][comptExogVar[282]])/(exogSamplingTime[282][comptExogVar[282] + 
    1] - exogSamplingTime[282][comptExogVar[282]]);
x[283] = dataExogVar[283][comptExogVar[283]] + (dataExogVar[283][comptExogVar[283] + 1] - dataExogVar[283][comptExogVar[283]]) * (t - exogSamplingTime[283][comptExogVar[283]])/(exogSamplingTime[283][comptExogVar[283] + 
    1] - exogSamplingTime[283][comptExogVar[283]]);
x[284] = dataExogVar[284][comptExogVar[284]] + (dataExogVar[284][comptExogVar[284] + 1] - dataExogVar[284][comptExogVar[284]]) * (t - exogSamplingTime[284][comptExogVar[284]])/(exogSamplingTime[284][comptExogVar[284] + 
    1] - exogSamplingTime[284][comptExogVar[284]]);
x[285] = dataExogVar[285][comptExogVar[285]] + (dataExogVar[285][comptExogVar[285] + 1] - dataExogVar[285][comptExogVar[285]]) * (t - exogSamplingTime[285][comptExogVar[285]])/(exogSamplingTime[285][comptExogVar[285] + 
    1] - exogSamplingTime[285][comptExogVar[285]]);
x[330] = (parms[11] - parms[12]) * pow((1.0/y[21]), parms[13]);
x[338] = 1.0/(1.0 + exp(-0.2 * (t - 1.0))) * (parms[22] - parms[22] * 1.8) + parms[22] * 1.8;
x[340] = y[28];
x[346] = parms[291] * (parms[292] * (y[0]) - y[33]);
x[351] = y[34] * y[46];
x[362] = 1.03 * y[45];
x[365] = ((1.0/(1.0 + exp(-parms[423] * (t - parms[424])))) * (parms[426] - parms[425]) + parms[425]);
x[367] = 1.03 * y[48];
x[378] = y[68];
x[384] = y[61];
x[385] = y[62];
x[387] = y[64];
x[390] = y[65];
x[392] = y[67];
x[394] = y[66];
x[397] = y[69];
x[398] = 1.015 * y[68];
x[399] = y[70];
x[410] = y[75];
x[416] = 1.015 * y[76];
x[418] = 1.015 * y[77];
x[419] = 1.015 * y[78];
x[420] = 1.015 * y[79];
x[421] = 1.015 * y[80];
x[630] = 0.0;
x[672] = y[76] * y[45];
x[677] = y[48] * y[78];
x[678] = y[47] * y[79];
x[740] = y[88] * y[82];
x[748] = y[94] * y[82];
x[755] = parms[208] * (0.35 * y[118] - y[97]);
x[772] = parms[225] * (0.58 * y[118] - y[109]);
x[786] = parms[246] * (0.34 * (y[113] + y[112] + y[111] + y[110]) - y[119]);
x[802] = parms[263] * (parms[262] * y[127] - y[129]);
x[806] = y[99] + y[133] + y[132];
x[811] = y[136]/(y[125] + y[126] + y[128] + y[127]);
x[820] = y[139];
ydot[12] = parms[314] * (x[384]/x[144] - y[12]);
ydot[13] = parms[315] * (x[385]/x[144] - y[13]);
ydot[15] = parms[317] * (x[384]/x[378] - y[15]);
ydot[16] = parms[318] * (x[385]/x[378] - y[16]);
ydot[45] = parms[275] * (x[362] - y[45]);
ydot[48] = parms[288] * (x[367] - y[48]);
ydot[68] = parms[399] * (x[398] - y[68]);
ydot[69] = 0.015 * x[397];
ydot[76] = parms[396] * (x[416] - y[76]);
ydot[77] = parms[154] * (x[418] - y[77]);
ydot[78] = parms[155] * (x[419] - y[78]);
ydot[79] = parms[156] * (x[420] - y[79]);
ydot[80] = parms[157] * (x[421] - y[80]);
ydot[97] = x[755];
ydot[99] = x[159];
ydot[109] = x[772];
ydot[114] = x[128];
ydot[117] = x[133];
ydot[119] = x[786];
ydot[129] = x[802];
ydot[135] = x[127];
x[287] = (y[0] + x[346])/(1.0 - x[65]);
x[288] = x[274] * (1.0 - parms[401]) + x[284] * parms[401];
x[298] = y[56] * (1.0 - parms[52] * x[1]) * (1.0 - parms[54] * x[25]);
x[332] = x[129]/x[397];
x[333] = x[131]/x[397];
x[334] = x[130]/y[68];
x[336] = (y[27]) * (1.0 + parms[26] * x[1]) * (1.0 + parms[27] * x[10]);
x[345] = x[273] + parms[400] * x[282];
x[363] = (1.0 + parms[281] * x[2]) * y[44] * pow((parms[278] * x[272]), parms[279]);
x[371] = y[58] * y[67] * (1.0 + parms[58] * x[18]) * (1.0 + parms[59] * x[9]) * (1.0 + parms[60] * x[5]);
x[375] = y[59] + parms[70] * x[1] + parms[71] * x[3] + parms[72] * x[10] + parms[73] * x[12];
x[386] = y[63] * (1.0 + parms[98] * x[26]);
x[391] = 1.015 * x[390];
x[393] = x[392] * (1.0 + x[268]) * (1.0 - x[251]) * (1.0 + x[201]) * (1.0 + x[255])/((1.0 + x[267]) * (1.0 - x[250]) * (1.0 + x[200])) - x[254];
x[395] = x[394] * (1.0 + x[265]) * (1.0 - x[249]) * (1.0 + x[199]) * (1.0 + x[253])/((1.0 + x[264]) * (1.0 - x[248]) * (1.0 + x[198])) - x[252];
x[404] = y[72] + parms[114] * x[16] + parms[115] * x[7] + parms[116] * x[29] + parms[117] * x[3];
x[406] = y[73] + parms[120] * x[1] + parms[121] * x[11];
x[408] = y[74] + parms[123] * x[36] + parms[124] * x[1];
x[411] = 1.015 * x[410];
x[417] = parms[128] * (x[182] * x[208]/(x[183] * x[209])) + parms[129] * (x[180] * x[208]/(x[181] * x[209])) + parms[130] * (x[178] * x[208]/(x[179] * x[209]));
x[422] = (parms[143] + parms[144] * x[141] + parms[145] * x[142]);
x[423] = (parms[146] + parms[147] * (x[211]/x[210]) + parms[148] * (x[209]/x[208])) * (1.0 + parms[149] * x[39]) * (1.0 - parms[402]);
x[427] = x[130];
x[428] = x[172] * (1.0 - parms[401]) + x[285] * parms[401];
x[432] = x[173] + x[165];
x[480] = x[244] * x[154];
x[488] = x[157] * x[189];
x[514] = x[231] * y[93]/100.0;
x[526] = x[231] * y[117]/100.0;
x[536] = x[227] * y[133]/100.0;
x[538] = x[227] * y[132]/100.0;
x[539] = x[227] * y[99]/100.0;
x[544] = x[231] * y[129]/100.0;
x[547] = x[227] * y[123]/100.0;
x[551] = x[227] * y[126]/100.0;
x[553] = x[227] * y[100]/100.0;
x[554] = x[227] * y[103]/100.0;
x[558] = x[231] * y[106]/100.0;
x[581] = x[54] * x[195];
x[582] = x[55] * x[194];
x[642] = parms[14] * x[167] * x[390] + parms[4] * x[172];
x[643] = parms[16] * x[167] * x[392];
x[644] = (parms[15] + x[3] * parms[6]) * x[167] * x[394] + parms[5] * x[172];
x[646] = x[131] + x[130];
x[652] = x[340] * x[399];
x[676] = x[278] * x[104]/1000000.0;
x[685] = x[242] + x[190];
x[712] = -x[138] - x[136] - x[139] - x[137] - x[140];
x[725] = x[60] * x[192];
x[738] = y[87] + x[740];
x[746] = x[740] - y[93];
x[747] = parms[196] * (x[748] - y[93]);
x[750] = x[225] * y[88];
x[753] = parms[204] * (parms[203] * y[32] * x[399] - y[137]);
x[754] = parms[206] * (parms[205] * y[32] * x[399] - y[96]);
x[764] = y[98] * parms[211] * (1.0 + parms[213] * x[0]) * (1.0 + parms[214] * x[5]) * (1.0 + parms[215] * x[37]) * pow(((1.0 + x[229]/100.0)/(1.0 + x[227]/100.0)), parms[212]);
x[765] = parms[217] * std::max(0.0, y[115]) * (1.0 + parms[218] * x[0]) * (1.0 + parms[219] * x[28]) * (1.0 + parms[220] * x[9]);
x[771] = parms[223] * (parms[222] * y[32] * x[399] - y[108]);
x[785] = 1.0 * (0.22 * x[399] * y[32] - y[118]);
x[789] = parms[248] * y[122] * (1.0 + parms[249] * x[49]) * (1.0 + parms[250] * x[43]);
x[791] = parms[252] * y[122] * (1.0 + parms[253] * x[3]) * (1.0 + parms[254] * x[1]) * (1.0 + parms[255] * x[5]) * (1.0 + parms[256] * x[8]);
x[794] = parms[258] * (y[124] + y[127]) * (parms[259] * pow(((1.0 + x[229]/100.0)/(1.0 + x[227]/100.0)), parms[260]));
x[798] = parms[258] * (1.0 - (parms[259] * pow(((1.0 + x[229]/100.0)/(1.0 + x[227]/100.0)), parms[260]))) * (y[124] + y[127]);
x[812] = y[140] - parms[322]/(1.0 + exp(-parms[323] * (1.0 + parms[326] * x[38]) * x[811]));
x[813] = x[228] + parms[338]/(1.0 + exp(-parms[339] * x[811]));
x[819] = y[138] + parms[324] * x[1] + parms[325] * x[3];
x[821] = y[140] + parms[335] * x[3] + parms[336] * x[10];
x[822] = y[141] + parms[340] * x[1] + parms[341] * x[8];
x[827] = y[142] + parms[347] * x[2];
x[829] = y[143] + parms[353] * x[8] + parms[354] * x[10];
x[831] = y[144] + parms[360] * x[16] + parms[361] * x[5];
x[833] = y[145] + parms[367] * x[2];
ydot[9] = parms[75] * (x[384]/x[375] - y[9]);
ydot[10] = parms[76] * (x[385]/x[375] - y[10]);
ydot[65] = parms[379] * (x[391] - y[65]);
ydot[75] = parms[153] * (x[411] - y[75]);
ydot[81] = parms[158] * (x[422] - y[81]);
ydot[82] = parms[159] * (x[423] - y[82]);
ydot[93] = x[747];
ydot[96] = x[754];
ydot[108] = x[771];
ydot[118] = x[785];
ydot[137] = x[753];
ydot[138] = parms[327] * (x[812] - y[138]);
ydot[141] = parms[342] * (x[813] - y[141]);
x[296] = x[288]/(y[55] + parms[47] * x[40] + parms[48] * x[2]);
x[297] = x[287]/x[298];
x[299] = (parms[49] + parms[50] * (x[336]/y[32] - x[219] - parms[51])) * (1.0 - parms[53] * x[21]);
x[310] = x[109] * x[288];
x[311] = x[111] * x[288];
x[312] = x[112] * x[288];
x[313] = x[110] * x[288];
x[314] = x[124] * x[287];
x[315] = x[123] * x[287];
x[316] = x[119] * x[287];
x[317] = x[121] * x[287];
x[318] = x[120] * x[287];
x[326] = (x[644] + x[330] * (y[83] - x[642] - x[644] - x[643]))/x[394];
x[331] = y[26] * x[333] + y[26] * x[332] + y[25] * x[334];
x[335] = x[336] + x[340] + x[345];
x[344] = (parms[37] * x[288]) * (1.0 + parms[38] * x[22]);
x[352] = (parms[266] + parms[267] * (y[46] - x[288])/(x[288] * pow((1.0 + parms[268] * x[404]/x[144]), parms[269])));
x[354] = 1.0/pow((1.0 + parms[160] * x[408]/x[392]), parms[161]) * ((1.0 + parms[168] * x[39]) * (1.0 + parms[169] * x[16]));
x[355] = 1.0/pow((1.0 + parms[162] * x[408]/x[399]), (parms[163] * (1.0 + parms[170] * x[39])));
x[356] = 1.0/pow((1.0 + parms[164] * x[408]/y[77]), parms[165]);
x[357] = 1.0/pow((1.0 + parms[166] * x[408]/x[385]), (parms[167] * (1.0 + parms[171] * x[1])));
x[359] = 1.0/pow((1.0 + parms[182] * x[406]/x[375]), parms[183]);
x[364] = (1.0 + x[365] * 19.0) * parms[276] * pow((y[81]/(x[375] * y[82])), parms[277]);
x[366] = parms[283] * x[363] + parms[284] * x[104];
x[368] = x[676]/y[80];
x[370] = (1.0 + x[158]) * x[371];
x[377] = parms[62] - parms[63] * y[33]/x[287];
x[400] = (parms[384] + parms[385] * x[144] + parms[386] * x[404] + parms[387] * x[406] + parms[388] * (1.0 - parms[385] - parms[386] - parms[387]) * y[71]) * (1.0 + parms[389] * x[4]) * (1.0 + parms[390] * x[8]);
x[405] = 1.015 * x[404];
x[407] = 1.015 * x[406];
x[409] = 1.015 * x[408];
x[424] = y[0] * x[375];
x[426] = x[375] * x[287];
x[433] = x[63] * x[428];
x[434] = x[66] * x[427];
x[444] = x[432] - x[270];
x[489] = x[156] * x[480];
x[504] = x[232]/(1.0 + x[232]) * x[428];
x[512] = x[822] * y[90]/100.0;
x[513] = x[822] * y[95]/100.0;
x[515] = x[821] * y[91]/100.0;
x[516] = x[827] * y[92]/100.0;
x[518] = x[822] * y[89]/100.0;
x[520] = x[822] * y[86]/100.0;
x[522] = x[822] * y[85]/100.0;
x[524] = std::max(0.0, x[822] * y[115]/100.0);
x[525] = x[822] * y[116]/100.0;
x[527] = x[829] * y[114]/100.0;
x[529] = x[822] * y[110]/100.0;
x[530] = x[821] * y[113]/100.0;
x[531] = x[822] * y[111]/100.0;
x[532] = x[822] * y[112]/100.0;
x[534] = x[822] * y[120]/100.0;
x[537] = x[821] * y[135]/100.0;
x[541] = x[822] * y[122]/100.0;
x[542] = x[822] * y[124]/100.0;
x[543] = x[822] * y[136]/100.0;
x[546] = x[822] * y[121]/100.0;
x[548] = x[822] * y[128]/100.0;
x[549] = x[822] * y[127]/100.0;
x[550] = x[821] * y[125]/100.0;
x[555] = x[822] * y[105]/100.0;
x[557] = x[822] * y[101]/100.0;
x[559] = x[822] * y[104]/100.0;
x[580] = x[581] - x[582];
x[623] = (x[820] - x[822]) * y[86]/100.0;
x[624] = (x[822] - x[819]) * y[85]/100.0;
x[625] = (x[820] - x[822]) * (std::max(0.0, y[115]) + y[116])/100.0;
x[626] = (x[822] - x[819]) * (y[110] + y[112])/100.0;
x[627] = (x[820] - x[822]) * y[95]/100.0;
x[628] = (x[822] - x[819]) * y[89]/100.0;
x[629] = (x[820] - x[822]) * y[105]/100.0;
x[645] = x[129] + x[646];
x[651] = x[336] * x[399];
x[656] = x[399] * x[345];
x[662] = x[404] * x[351];
x[670] = ((x[831] - x[231]) * y[93] + (x[833] - x[231]) * (y[117] + y[129]) + (x[822] - x[819]) * y[101])/100.0;
x[674] = x[363] * y[77];
x[675] = x[277] * x[417];
x[684] = x[685] - x[191];
x[745] = parms[195] * (x[746] - y[92]);
x[749] = parms[197] * (x[750] - y[94]);
x[763] = parms[210] * (x[764] - y[103]);
x[766] = parms[216] * (x[765] - y[104]);
x[769] = x[133] + x[802] + x[747];
x[770] = x[772] + x[771];
x[788] = parms[247] * (x[789] - y[120]);
x[790] = parms[251] * (x[791] - y[121]);
x[793] = parms[257] * (x[794] - y[123]);
x[797] = parms[261] * (x[798] - y[126]);
x[803] = x[785] - x[755] - x[772];
x[804] = x[786] + x[754] + x[771];
x[818] = (x[228] * y[136] + x[819] * y[122])/(y[136] + y[122]);
ydot[31] = parms[40] * (x[400]/x[399] - y[31]);
ydot[34] = parms[273] * (x[352] - y[34]);
ydot[35] = parms[172] * (x[354] - y[35]);
ydot[36] = parms[173] * (x[355] - y[36]);
ydot[37] = parms[174] * (x[356] - y[37]);
ydot[38] = parms[175] * (x[357] - y[38]);
ydot[39] = parms[187] * (x[359] - y[39]);
ydot[44] = parms[282] * (x[364] - y[44]);
ydot[47] = parms[285] * (x[366] - y[47]);
ydot[58] = (parms[55] + parms[56] * parms[427] + parms[57] * x[297]/x[167]) * y[58];
ydot[60] = parms[69] * (x[377] - y[60]);
ydot[92] = x[745];
ydot[94] = x[749];
ydot[103] = x[763];
ydot[104] = x[766];
ydot[106] = x[769];
ydot[107] = x[770];
ydot[120] = x[788];
ydot[121] = x[790];
ydot[123] = x[793];
ydot[126] = x[797];
ydot[130] = x[803];
ydot[131] = x[804];
x[291] = x[287]/y[8] - (y[9] * x[316] + y[10] * (x[318] + x[317]))/y[8];
x[292] = x[288]/y[11] - (y[12] * x[310] + y[13] * (x[313] + x[312]))/y[11];
x[295] = x[296] + x[297] + x[154];
x[301] = x[238] * x[296];
x[302] = x[239] * x[297];
x[303] = x[237] * x[297];
x[304] = x[236] * x[297];
x[319] = x[311] + x[314];
x[329] = x[326]/((1.0 + x[265]) * (1.0 - x[249]) * (1.0 + x[253]) * (1.0 + x[199]));
x[343] = x[335]/y[30] - x[344] * y[31]/y[30];
x[361] = y[49] * y[45] + y[50] * x[363] + y[51] * x[277] + y[52] * y[48] + y[53] * y[47] + y[54] * x[368];
x[372] = x[256] * x[370];
x[376] = (1.0 + y[60]) * (x[370]/(parms[64] * x[298]) + (x[384] * (1.0 + parms[65] * x[7]) * (1.0 + parms[66] * x[16]) * (1.0 + parms[67] * x[36])) * x[119] + x[385] * x[120] + x[385] * x[121]);
x[396] = x[645]/x[331];
x[401] = 1.015 * x[400];
x[412] = parms[151] * (x[407] - y[73]);
x[413] = parms[150] * (x[405] - y[72]);
x[414] = parms[152] * (x[409] - y[74]);
x[429] = x[427] + x[434];
x[435] = -x[434] - x[433];
x[456] = x[216] * x[662];
x[498] = x[444] - x[189] - x[488] - x[143];
x[511] = x[512] + x[513] + x[514] + x[515] + x[516];
x[517] = x[518];
x[519] = x[520];
x[521] = x[522];
x[523] = x[524] + x[525] + x[526] + x[527];
x[528] = x[529] + x[530] + x[531] + x[532];
x[533] = x[534];
x[535] = x[536] + x[537] + x[538] + x[543] + x[539];
x[540] = x[541] + x[542] + x[543] + x[544];
x[545] = x[534] + x[546] + x[547] + x[548] + x[549] + x[550] + x[551];
x[552] = x[553] + x[554] + x[555] + x[539];
x[560] = x[527] + x[516];
x[569] = x[101] * x[426];
x[570] = x[103] * x[426];
x[594] = x[384] * x[310];
x[595] = x[384] * x[316];
x[598] = x[386] * x[311];
x[599] = x[387] * x[314];
x[601] = x[385] * x[313];
x[602] = x[385] * x[318];
x[606] = x[385] * x[312];
x[608] = x[385] * x[317];
x[614] = x[625] + x[626];
x[622] = x[627] + x[628];
x[634] = x[326] * x[394];
x[636] = x[623] + x[624];
x[650] = x[651] + x[652] + x[656];
x[655] = x[344] * x[400];
x[673] = x[674] + x[675] + x[677] + x[678] + x[676];
x[681] = x[629] + x[630];
x[728] = x[58] * x[656];
x[762] = x[745] + x[128];
x[807] = x[763] - x[797];
ydot[1] = parms[290] * (x[343]/y[32] - x[219] - y[1]);
ydot[25] = parms[309] * (y[68]/x[396] - y[25]);
ydot[26] = parms[310] * (x[397]/x[396] - y[26]);
ydot[32] = x[343] - x[219] * y[32];
ydot[61] = parms[85] * x[145] * (1.0 + parms[88] * x[24]) + parms[86] * x[412] + parms[87] * (1.0 - parms[85] - parms[86]) * x[413] * (1.0 + parms[89] * x[24]);
ydot[71] = parms[391] * (x[401] - y[71]);
ydot[72] = x[413];
ydot[73] = x[412];
ydot[74] = x[414];
ydot[102] = x[762];
ydot[132] = x[807];
x[289] = x[429]/x[378];
x[300] = x[302] + x[301] + x[154];
x[305] = (1.0 - x[295]/x[161]) * 100.0;
x[374] = parms[68] * (x[376] - x[375]);
x[452] = x[456];
x[461] = x[252] * x[329];
x[464] = x[264] * x[634]/(1.0 + x[264]);
x[469] = x[198] * x[634]/((1.0 + x[198]) * (1.0 - x[248]) * (1.0 + x[264]));
x[475] = x[248] * x[634]/((1.0 - x[248]) * (1.0 + x[264]));
x[478] = y[57] * x[301];
x[479] = x[371] * x[302];
x[507] = x[233]/(1.0 + x[233]) * x[429];
x[556] = x[557] + x[558] + x[560] + x[559];
x[596] = x[114] * x[429];
x[597] = x[598] + x[599];
x[603] = x[115] * x[429];
x[610] = x[116] * x[429];
x[612] = x[160] * x[614];
x[615] = x[164] * x[622];
x[637] = x[52] * x[636];
x[648] = (1.0 - x[164]) * x[622];
x[654] = x[651] - x[655];
x[671] = x[672] + x[673];
x[682] = x[53] * x[681];
x[702] = (x[498] + x[535] - x[533] + x[136]) * (1.0 - x[260]) - x[134] - x[243];
x[722] = x[260] * (x[498] + x[535] - x[533] + x[136]);
x[729] = x[59] * (x[511] + x[622]);
x[752] = parms[199] * x[429] * (1.0 + parms[200] * x[1]) * (1.0 + parms[201] * x[14]);
x[779] = (y[110] + y[111]) * parms[235] * (1.0 + parms[237] * x[46]) * (1.0 + parms[238] * x[14]) * pow((x[361]/x[287]), parms[236]) + +std::max(0.0, 0.1 * (y[32] * x[399]) - y[115]);
ydot[59] = x[374];
ydot[62] = parms[108] * (1.0 + parms[111] * x[23]) * x[374] + parms[109] * x[412] + parms[110] * (1.0 - parms[108] - parms[109]) * x[414];
ydot[63] = parms[94] * x[374] * (1.0 + parms[97] * x[23]) + parms[95] * x[412] + parms[96] * (1.0 - parms[94] - parms[95]) * x[414];
ydot[64] = parms[99] * (1.0 + parms[102] * x[23]) * x[374] + parms[100] * x[412] + parms[101] * (1.0 - parms[99] - parms[100]) * x[414];
ydot[66] = parms[77] * x[374] * (1.0 + parms[80] * x[23]) + parms[78] * x[412] + parms[79] * (1.0 - parms[77] - parms[78]) * x[414];
ydot[67] = parms[103] * x[374] * (1.0 + parms[106] * x[23]) * (1.0 + parms[107] * x[19]) + parms[104] * x[412] + parms[105] * (1.0 - parms[103] - parms[104]) * x[414];
ydot[70] = parms[81] * (1.0 + parms[84] * x[23]) * x[374] + parms[82] * x[412] + parms[83] * (1.0 - parms[81] - parms[82]) * x[414];
x[320] = x[596]/x[384];
x[321] = x[603]/x[385];
x[415] = x[671]/x[361];
x[467] = x[269]/(1.0 + x[269]) * (x[636] + x[637]);
x[471] = x[202]/(1.0 + x[202]) * (x[636] + x[637]);
x[477] = x[479] + x[478] + x[480];
x[482] = x[105] * x[479] - x[189];
x[484] = x[155] * x[478];
x[485] = x[158] * x[479];
x[583] = x[241] * (x[300] * 1000.0 + x[187])/1000000.0;
x[593] = x[594] + x[595] + x[596];
x[600] = x[601] + x[602] + x[603];
x[613] = x[614] - x[612];
x[617] = x[113] * x[612];
x[618] = x[118] * x[615];
x[649] = x[125] * x[648];
x[653] = x[652] + x[656] + x[654];
x[679] = x[671] - x[681] - x[682];
x[680] = x[673] - x[681] - x[682];
x[703] = x[702];
x[751] = parms[198] * (x[752] - y[95]) * (1.0 + parms[202] * x[14]);
x[778] = parms[234] * (x[779] - y[112]);
ydot[49] = parms[303] * (y[76]/x[415] - y[49]);
ydot[50] = parms[304] * (y[77]/x[415] - y[50]);
ydot[51] = parms[305] * (x[417]/x[415] - y[51]);
ydot[52] = parms[306] * (y[78]/x[415] - y[52]);
ydot[53] = parms[307] * (y[79]/x[415] - y[53]);
ydot[54] = parms[308] * (y[80]/x[415] - y[54]);
ydot[95] = x[751];
ydot[112] = x[778];
x[322] = (x[610] + x[615] + x[618])/x[385];
x[402] = x[653]/x[343];
x[445] = x[429] - x[596] - x[603] - x[610] - x[615] - x[618];
x[455] = x[467] + x[471];
x[466] = x[266] * x[653]/(1.0 + x[266]);
x[481] = x[479] - x[482] - x[189];
x[483] = x[484] + x[485] + x[489];
x[487] = x[157] * x[482];
x[584] = x[176] * x[583]/100.0;
x[585] = x[175] * x[583]/100.0;
x[586] = x[174] * x[583]/100.0;
x[587] = x[177] * x[583]/100.0;
x[605] = x[606] - x[612] - x[617];
x[611] = x[612] + x[613] + x[615];
x[619] = x[122] * x[613];
x[647] = x[131] - x[648] - x[649];
x[704] = x[703];
ydot[30] = parms[39] * (x[402]/x[399] - y[30]);
x[293] = x[289]/y[14] - (y[15] * x[320] + y[16] * (x[321] + x[322]))/y[14];
x[440] = x[428] - x[594] - x[601] - x[605] - x[612] - x[617];
x[486] = x[485] - x[487] - x[488];
x[499] = x[445] - x[480] - x[507] - x[489];
x[607] = x[608] - x[271] - x[270] - x[613] - x[619];
x[616] = x[617] + x[619] + x[618];
x[620] = x[611] + x[636] + x[648] + x[681];
x[726] = x[56] * (x[129] + x[131] + x[586]);
x[743] = x[149] * x[445];
x[294] = y[17] * x[292] + y[18] * x[291] + y[19] * x[293];
x[381] = x[440]/x[292];
x[382] = x[445]/x[293];
x[501] = x[440] - x[478] - x[484] - x[504];
x[590] = x[440] - x[484] - x[504];
x[604] = x[605] + x[607] + x[271] + x[270] + x[610];
x[621] = x[616] + x[637] + x[649] + x[682];
x[724] = x[725] + x[726];
x[742] = parms[194] * (x[743] - y[89]);
x[775] = parms[227] * (x[481] + x[486]) * (1.0 + parms[228] * x[1]) * (1.0 + parms[229] * x[29]) + std::max(0.0, 0.1 * (y[32] * x[399]) - y[115]);
x[777] = parms[231] * (x[481] + x[486]) * (1.0 + parms[232] * x[47]) * (1.0 + parms[233] * x[14]) + std::max(0.0, 0.1 * (y[32] * x[399]) - y[115]);
ydot[11] = parms[313] * (x[381]/x[144] - y[11]);
ydot[14] = parms[316] * (x[382]/x[378] - y[14]);
ydot[89] = x[742];
x[449] = x[61] * (x[620] + x[621]);
x[640] = parms[1] - parms[3] * x[590]/x[168];
x[774] = parms[226] * (x[775] - y[110]);
x[776] = parms[230] * (x[777] - y[111]);
ydot[110] = x[774];
ydot[111] = x[776];
x[436] = x[64] * ((x[620] + x[621]) - (x[449] + x[455] + x[670]));
x[757] = x[776] + x[790];
ydot[90] = x[757];
x[431] = (x[620] + x[621]) - (x[449] + x[455] + x[670]) + x[436];
x[437] = x[435] - x[436] - x[165];
x[430] = x[426] - x[431];
x[443] = x[431] - x[271] - x[173];
x[506] = x[234]/(1.0 + x[234]) * x[431];
x[510] = x[245] * x[431]/((1.0 + x[234]) * (1.0 - x[245]));
x[609] = x[117] * x[431];
x[690] = x[257] * x[431];
x[692] = x[204] * x[431];
x[497] = x[443] - x[482] - x[487] - x[506] + x[510];
x[502] = 0.15 * x[430];
x[505] = x[235]/(1.0 + x[235]) * x[430];
x[509] = x[246] * x[430]/((1.0 + x[235]) * (1.0 - x[246]));
x[563] = 0.023 * x[430];
x[568] = x[102] * x[430];
x[574] = x[71] * x[430];
x[575] = x[70] * x[430];
x[689] = x[258] * x[430];
x[691] = x[206] * x[430];
x[490] = x[203] * (x[477] + x[502]);
x[500] = x[501] + x[502];
x[503] = x[504] + x[505] + x[506] + x[143] + x[507];
x[508] = x[509] + x[510];
x[567] = x[568] + x[569] + x[570];
x[576] = x[68] * x[574];
x[577] = x[69] * x[574];
x[578] = x[67] * x[574];
x[579] = x[72] * x[574];
x[687] = x[108] * (x[689] + x[690] + x[243] + x[259]);
x[688] = (1.0 - x[108]) * (x[689] + x[690] + x[243] + x[259]);
x[491] = x[212] * x[490];
x[492] = x[215] * x[490];
x[493] = x[214] * x[490];
x[494] = x[213] * x[490];
x[572] = x[100] * x[567];
x[573] = x[99] * x[567];
x[562] = x[221] * std::max(0.0, (x[497] + x[545] - x[540] + x[494] - x[585] - x[690] - x[185]) - x[134]);
x[571] = x[572] - x[569];
x[561] = x[563] + x[562] + x[134];
x[699] = (x[497] + x[545] - x[540] + x[138] - x[562] + x[576] + x[494] - x[585] - x[690]) * (1.0 - x[261]);
x[564] = x[222] * x[561];
x[565] = x[220] * x[561];
x[566] = x[224] * x[561];
x[700] = x[699];
x[589] = (x[479] + x[480] + x[502] + x[521] - x[519] + x[139] + x[564] + x[242] + x[687] + x[583] + x[571] + x[577] - x[490]) * (1.0 - x[262]);
x[701] = x[700] + x[692];
x[721] = x[261] * x[700]/(1.0 - x[261]);
x[342] = (parms[32] + parms[33] * y[86]/x[589]) * (1.0 + parms[34] * x[1]) * (1.0 + parms[35] * x[3]);
x[588] = x[589] + x[590];
x[639] = parms[0] - parms[2] * x[589]/x[169];
x[641] = std::max(0.0, (parms[8] - parms[9] * (x[589]/x[167])) * pow((1.0/y[20]), parms[10]));
x[717] = x[262] * x[589]/(1.0 - x[262]);
x[736] = (parms[190] + parms[191] * (parms[192] - 4.0 * x[820]/100.0)) * x[589];
ydot[29] = parms[29] * (x[342] - y[29]);
x[324] = (x[642] + x[641] * (y[83] - x[642] - x[644] - x[643]))/x[390];
x[638] = (x[639] * x[589] + x[640] * x[590]) * (1.0 + parms[17] * x[1] + parms[18] * x[2]);
x[705] = x[588] - y[83];
x[716] = x[97] * x[717];
x[733] = x[283] * x[588];
x[735] = parms[189] * (x[736] - y[86]);
ydot[83] = parms[7] * (x[638] - y[83]);
ydot[86] = x[735];
x[309] = parms[270] * (x[310] + x[316] + x[320] + x[332] + x[324] + x[344] + y[45]);
x[327] = (y[83] - x[326] * x[394] - x[324] * x[390])/x[392];
x[631] = x[324] * x[390];
x[732] = parms[188] * (x[733] - y[84]);
ydot[46] = parms[274] * (x[309] - y[46]);
ydot[84] = x[732];
x[323] = x[324] * y[20] + x[326] * y[21] + x[327] * y[22];
x[325] = x[326] * y[23] + x[327] * y[24];
x[328] = x[327]/((1.0 + x[268]) * (1.0 - x[251]) * (1.0 + x[255]) * (1.0 + x[201]));
x[473] = x[247] * x[631]/(1.0 - x[247]);
x[591] = x[129] + x[672] + x[655] + x[631] + x[593];
x[633] = x[327] * x[392];
x[809] = x[788] + x[732];
ydot[134] = x[809];
x[353] = (y[35] * x[323] + y[36] * (x[335] - x[344]) + y[37] * x[363] + y[38] * x[319]);
x[358] = (y[39] * pow((x[323] + x[335] - x[344] + x[361] + x[333]), parms[181])) * (1.0 + parms[184] * x[1]) * (1.0 + parms[185] * x[42]) * (1.0 + parms[186] * x[31]);
x[388] = y[83]/x[323];
x[462] = x[254] * x[328];
x[465] = x[267] * x[633]/(1.0 + x[267]);
x[470] = x[200] * x[633]/((1.0 + x[200]) * (1.0 - x[250]) * (1.0 + x[267]));
x[476] = x[250] * x[633]/((1.0 - x[250]) * (1.0 + x[267]));
x[632] = x[633] + x[634];
x[635] = x[633] - x[636] - x[637];
ydot[20] = parms[19] * (x[390]/x[388] - y[20]);
ydot[21] = parms[20] * (x[394]/x[388] - y[21]);
ydot[22] = parms[21] * (x[392]/x[388] - y[22]);
x[306] = x[588]/x[388];
x[307] = x[590]/x[388];
x[308] = x[589]/x[388];
x[360] = (parms[176] * x[353] + parms[177] * x[358] + parms[178] * x[351]) * (1.0 + parms[179] * x[15]) * (1.0 + parms[180] * x[14]);
x[369] = x[276] * x[388];
x[373] = x[275] * x[388];
x[389] = x[632]/x[325];
x[460] = x[461] + x[462];
x[463] = x[464] + x[465] + x[466];
x[468] = x[469] + x[470];
x[474] = x[475] + x[476];
x[592] = x[597] + x[632] + x[131] + x[653] + x[673];
x[664] = x[353] * x[408];
x[665] = x[358] * x[406];
ydot[23] = parms[311] * (x[394]/x[389] - y[23]);
ydot[24] = parms[312] * (x[392]/x[389] - y[24]);
ydot[57] = parms[61] * (x[369] - y[57]);
x[341] = x[308] * (parms[30] - parms[31] * x[342]);
x[350] = y[40] * x[351] + y[41] * x[353] + y[42] * x[358] + y[43] * x[360];
x[447] = x[62] * x[592];
x[472] = x[473] + x[474];
x[666] = x[360] * x[410];
x[667] = x[664] * x[240];
x[668] = x[664] * (1.0 - x[240]);
x[715] = x[106] * x[463];
x[727] = x[57] * (x[473] + x[474]);
ydot[28] = parms[36] * (x[341] - x[340]);
x[446] = -x[447];
x[448] = x[447] - x[449];
x[458] = x[217] * x[667];
x[459] = x[218] * (x[668] + x[665] + x[666]);
x[659] = x[428] - x[447] + x[452] - x[473] + x[662] - x[591] - x[433];
x[663] = x[664] + x[665] + x[666];
x[723] = x[724] + x[727] + x[729] + x[728];
x[348] = x[659]/x[147];
x[457] = x[458] + x[459];
x[661] = x[662] + x[663];
x[669] = x[663] - x[670];
ydot[33] = (x[426] - x[435] + x[447] + x[457] + x[468] - x[474] + x[460] + x[463] + x[663] - x[592])/x[148];
x[403] = x[661]/x[350];
x[453] = x[457] + x[460] + x[463] + x[468];
x[660] = x[661] - x[671];
x[683] = x[671] - x[661] + x[135] + x[197] + x[684];
x[741] = std::max(x[226] * y[87]/y[82], y[88] + (x[166] * x[661] - x[806])/y[82]);
x[761] = x[279] * (x[671] + x[661]);
x[768] = x[280] * (x[671] + x[661]);
x[784] = parms[241] * (x[671] + x[661]) * (1.0 + parms[242] * x[45]) * (1.0 + parms[243] * x[50]) * (1.0 + parms[244] * x[11]);
ydot[40] = parms[299] * (x[404]/x[403] - y[40]);
ydot[41] = parms[302] * (x[408]/x[403] - y[41]);
ydot[42] = parms[300] * (x[406]/x[403] - y[42]);
ydot[43] = parms[301] * (x[410]/x[403] - y[43]);
ydot[88] = parms[193] * (x[741] - y[88]);
x[425] = x[592] - x[447] - x[453] + x[474] - x[663];
x[451] = x[452] + x[453];
x[454] = x[453] - x[455];
x[658] = x[426] + x[447] + x[453] - x[474] - x[435] + x[663] - x[592];
x[710] = x[660] + x[556] - x[552] + x[140] + x[566] - x[570] + x[579] - x[580] - x[242] + x[492] - x[587] - x[259];
x[760] = parms[209] * (x[761] - y[101]);
x[767] = parms[221] * (x[768] - y[105]);
x[783] = parms[240] * (x[784] - y[116]);
ydot[101] = x[760];
ydot[105] = x[767];
ydot[116] = x[783];
x[286] = x[425]/x[375];
x[349] = x[658]/x[148];
x[442] = (x[600] + x[604] + x[634] + x[635] + x[647] + x[653] + x[658] + x[680]) - (x[448] + x[454] - x[474]) + x[437] - x[669] - (x[595] + x[602] + x[607] + x[613] + x[619]);
x[450] = x[451] - x[472];
x[657] = x[659] + x[658];
x[795] = x[778] + x[760];
x[801] = x[783] + x[767];
ydot[0] = parms[289] * (x[286] - y[0]) + y[0] * y[1];
ydot[124] = x[795];
ydot[127] = x[801];
x[347] = x[657]/x[146];
x[438] = y[83] + x[653] + x[655] + x[645] + x[671] - x[661] + x[657];
x[441] = x[442] + x[443] + x[444];
x[496] = x[442] - x[481] - x[486] - x[505] + x[509];
x[781] = x[281] * x[442] + std::max(0.0, 0.1 * (y[32] * x[399]) - y[115]);
x[290] = y[2] * x[323] + y[3] * x[335] + y[4] * x[331] + y[5] * x[361] - y[6] * x[350] + y[7] * x[347];
x[339] = x[496] - x[502] + x[528] - x[523] + x[712] - x[568] - x[574] + x[575];
x[380] = x[441]/x[291];
x[439] = x[440] + x[441] + x[445];
x[495] = x[496] + x[497] + x[498];
x[686] = x[683]/x[438] * 100.0;
x[693] = x[205] * x[438];
x[694] = x[207] * x[438];
x[696] = (x[496] - x[502] + x[528] - x[523] + x[712] - x[563] - x[568] - x[574] + x[575] + x[580] + x[493] - x[584] - x[689]) * (1.0 - x[263]);
x[737] = ((1.0/(1.0 + exp(-parms[411] * (t - parms[412])))) * (parms[414] - parms[413]) + parms[413]) * x[438];
x[756] = ((1.0/(1.0 + exp(-parms[415] * (t - parms[416])))) * (parms[418] - parms[417]) + parms[417]) * x[438];
x[773] = ((1.0/(1.0 + exp(-parms[419] * (t - parms[420])))) * (parms[422] - parms[421]) + parms[421]) * x[438];
x[780] = parms[239] * (x[781] - y[113]);
x[787] = ((1.0/(1.0 + exp(-parms[403] * (t - parms[404])))) * (parms[406] - parms[405]) + parms[405]) * x[438];
x[805] = ((1.0/(1.0 + exp(-parms[407] * (t - parms[408])))) * (parms[410] - parms[409]) + parms[409]) * x[438];
x[814] = (y[114] + y[115] + y[116] + y[117])/x[496];
x[816] = (y[91] + y[92] + y[95] + y[93])/x[438];
x[823] = (y[114])/x[438];
x[824] = (y[91] + y[92] + y[95] + y[93])/x[438];
x[825] = (-y[100] + y[101] + y[92] - y[103] + y[104] - y[105] + y[106] + y[108] + y[109])/x[438];
ydot[8] = parms[74] * (x[380]/x[375] - y[8]);
ydot[113] = x[780];
x[337] = y[32] * (x[338] + parms[23] * (x[291]/y[32]) + parms[24] * (x[339]/(y[32] * x[399])));
x[379] = x[438]/x[290];
x[383] = x[439]/x[294];
x[695] = -x[691] - x[692] - x[693] - x[694];
x[697] = x[696];
x[706] = x[705] + x[693] - x[652] - x[655];
x[711] = x[710] + x[694];
x[810] = -(x[787] + x[805] + x[737] + x[756] + x[773]);
x[815] = x[818] + parms[328] * (1.0 + parms[331] * x[21]) + parms[329]/(1.0 + exp(-parms[330] * x[814]));
x[817] = x[228] + parms[333]/(1.0 + exp(-parms[334] * x[816]));
x[828] = x[229] + parms[349]/(1.0 + exp(-parms[350] * x[825])) + parms[351]/(1.0 + exp(-parms[352] * x[823]));
x[832] = x[231] + parms[363]/(1.0 + exp(-parms[364] * x[825])) + parms[365]/(1.0 + exp(-parms[366] * x[823]));
ydot[2] = parms[293] * (x[388]/x[379] - y[2]);
ydot[3] = parms[294] * (x[399]/x[379] - y[3]);
ydot[4] = parms[295] * (x[396]/x[379] - y[4]);
ydot[5] = parms[297] * (x[415]/x[379] - y[5]);
ydot[6] = parms[296] * (x[403]/x[379] - y[6]);
ydot[7] = parms[298] * (x[146]/x[379] - y[7]);
ydot[17] = parms[319] * (x[381]/x[383] - y[17]);
ydot[18] = parms[320] * (x[380]/x[383] - y[18]);
ydot[19] = parms[321] * (x[382]/x[383] - y[19]);
ydot[27] = parms[28] * (x[337] - x[336]) * (1.0 + parms[25] * x[23]);
ydot[139] = parms[332] * (x[815] - x[820]);
ydot[140] = parms[337] * (x[817] - y[140]);
ydot[143] = parms[355] * (x[828] - y[143]);
ydot[145] = parms[368] * (x[832] - y[145]);
x[698] = x[697] + x[691] - x[654] - x[657];
x[720] = x[263] * x[697]/(1.0 - x[263]);
x[734] = x[706] - x[732] + x[735] + x[737];
x[758] = -x[711] + x[760] + x[762] + x[766] - x[767] + x[769] + x[770] - x[773];
ydot[85] = x[734];
ydot[98] = x[758];
x[719] = x[720] + x[721] + x[722];
x[759] = x[758] - x[159] - x[763];
x[782] = -x[698] + x[774] + x[776] + x[778] + x[780] - x[128] - x[783] - x[133] + x[786] - x[785] - x[787];
x[792] = x[774] + x[734] + x[742];
ydot[100] = x[759];
ydot[115] = x[782];
ydot[122] = x[792];
x[707] = x[499] + x[451] - x[472] + x[503] - x[508] + x[517] - x[511] + x[137] + x[565] + x[573] + x[578] + x[719] + x[717] + x[491] + x[483] - x[586] + x[688];
x[718] = x[98] * x[719];
x[800] = x[782] + x[735] + x[751] - x[766];
x[808] = x[759] - x[793];
ydot[128] = x[800];
ydot[133] = x[808];
x[708] = x[707] - x[645];
x[714] = x[716] + x[718] + x[184] + x[456] + x[457] + x[715] + x[193] + x[126] + x[188] + x[186];
x[709] = x[708] + x[695] - x[656];
x[713] = x[714] - x[715] * 0.3;
x[730] = x[713] - x[723] + x[196];
x[744] = -x[709] + x[742] - x[757] - x[745] - x[751] - x[747] + x[753] + x[754] + x[755] - x[756];
ydot[91] = x[744];
x[731] = x[730]/x[438] * 100.0;
x[739] = x[744] + x[751];
x[796] = x[744] - x[127] - x[780];
ydot[87] = x[739];
ydot[125] = x[796];
x[799] = -x[701] + x[788] - x[792] + x[790] + x[793] - x[795] + x[796] + x[797] + x[800] + x[801] - x[802] - x[804] + x[803] - x[805];
x[826] = x[229] + parms[343]/(1.0 + exp(-parms[344] * x[825])) + parms[345]/(1.0 + exp(parms[346] * x[731]));
x[830] = x[231] + parms[356]/(1.0 + exp(-parms[357] * x[825])) + parms[358]/(1.0 + exp(-parms[359] * x[731]));
ydot[136] = x[799];
ydot[142] = parms[348] * (x[826] - y[142]);
ydot[144] = parms[362] * (x[830] - y[144]);
}
	

class minDistRK4
{
  
public:
  // CONSTRUCTOR
  minDistRK4(double* yInitIn,
             unsigned int ntIn,
             double bytIn,
             double** dataExogVarIn, 
             double** exogSamplingTimeIn, 
             int nExogVarIn,
             unsigned int* varMinDistIn,
             unsigned int nVarMinDistIn,
             double** dataMinDistIn,
             double** samplingTimeIn,
             double** pointsWeightIn) :
  yInit(yInitIn),
  nt(ntIn),
  byt(bytIn),
  dataExogVar(dataExogVarIn),
  exogSamplingTime(exogSamplingTimeIn),
  nExogVar(nExogVarIn),
  varMinDist(varMinDistIn), 
  nVarMinDist(nVarMinDistIn),
	dataMinDist(dataMinDistIn),
  samplingTime(samplingTimeIn),
  pointsWeight(pointsWeightIn) {}
	
	// COMPUTES INTERMEDIATE VARIABLES AT T=0	
	void getFullInitialPosition(double t, double* y, double* parms, double* yx, double** dataExogVar, double** exogSamplingTime, int nExogVar) {
		double ydot0[dim];
		for (unsigned int it=0; it<dim; it++) yx[it] = y[it];
		double* yptr = &yx[0];
		double* xptr = &yx[dim];
		int comptExogVar[nExogVar];
		for (unsigned int it=0; it<nExogVar; it++) comptExogVar[it]=1;
		Func(0, yptr, parms, ydot0, xptr, dataExogVar, exogSamplingTime, nExogVar, comptExogVar); // because yptr and xptr are pointers to elements of yx, elements in yx are directly filled by Func without having to do any copy
	}
	
	double EvaluateWithoutTryCatch(double* parms) {
		// INIT
		int it, it1, indexData;
		
		// To store number of points already explored in dataMinDist
		unsigned int arraycomptVarMinDist[nVarMinDist];
		for (it=0; it<nVarMinDist; it++) {
			arraycomptVarMinDist[it] = 1;
		}
		unsigned int* comptVarMinDist = &arraycomptVarMinDist[0];
		
		// Initialize pointer y
		double yy[dim];
		for (it=0; it<dim; ++it) {
			yy[it] = yInit[it];
		}
		double *y = &yy[0];
		double yx[dim+dimIv];
		double y1[dim], y2[dim], y3[dim], ydot0[dim], ydot1[dim], ydot2[dim], ydot3[dim], ydots[dim], x0[dimIv], x1[dimIv], x2[dimIv], x3[dimIv];
		
		double fit = 0;
		
		int comptExogVar[nExogVar];
		for (it=0; it<nExogVar; it++) {
		  comptExogVar[it] = 1;
		}
		// get intermediateVar and compute distance at t=0 //
		getFullInitialPosition(0, y, parms, yx, dataExogVar, exogSamplingTime, nExogVar);
		for (it1=0; it1<nVarMinDist; it1++) { // update fit
			if ( 0 >= samplingTime[it1][0]) {
				fit+=pointsWeight[it1][0]*std::abs((yx[varMinDist[it1]] - dataMinDist[it1][0])/dataMinDist[it1][0]);
			}
		}
		
		// MAIN LOOP (RK4 and compute distance)
		for (it=0; it<(nt-1); it++) {
			
			
			
			Func(it*byt, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++)
				y1[it1] = y[it1] + ydot0[it1]*0.5*byt;
			Func((it + 0.5)*byt, y1, parms, ydot1, x1, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++)
				y2[it1] = y[it1] + ydot1[it1]*0.5*byt;
			Func((it + 0.5)*byt, y2, parms, ydot2, x2, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++)
				y3[it1] = y[it1] + ydot2[it1]*byt;
			Func((it+1)*byt, y3, parms, ydot3, x3, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++) {
				ydots[it1] = (ydot0[it1] + 2.0*ydot1[it1] + 2.0*ydot2[it1] + ydot3[it1])/6.0;
				y[it1] = y[it1] + byt*ydots[it1];
				yx[it1] = y[it1];
			}
			
			for (it1=0; it1<dimIv; it1++) {
				yx[dim + it1] = (x0[it1] + 2.0*x1[it1] + 2.0*x2[it1] + x3[it1])/6.0;
			}
			for (it1=0; it1<nVarMinDist; it1++) { // update fit
				if ( (double) (it + 1)*byt>= samplingTime[it1][comptVarMinDist[it1]]) {
					fit+=pointsWeight[it1][comptVarMinDist[it1]]*std::abs((yx[varMinDist[it1]] - dataMinDist[it1][comptVarMinDist[it1]])/dataMinDist[it1][comptVarMinDist[it1]]);
					comptVarMinDist[it1]++;
				}
			}
		}
		return fit;
	}
	
	double Evaluate(double* parms) {
		double out = 0;
		try {
			out = EvaluateWithoutTryCatch(parms);
		}
		catch(...) {
			out = 1e50;
		}
		if (std::isnan(out)) out = 1e50;
		
		return out;
	}
	
private:
	double* yInit;
	unsigned int nt;
	double byt;
	double** dataExogVar; 
	double** exogSamplingTime; 
	int nExogVar;
	unsigned int* varMinDist;
	unsigned int nVarMinDist;
	double** dataMinDist;
	double** samplingTime;
	double** pointsWeight;
};

template<typename RK4>
class fixedMinDistRK4 {
public:
	fixedMinDistRK4(RK4* myRK4In, 
                 const bool* parmsFixedIn, 
                 const double* parmsFullInitIn,
                 const double* parmsLowerIn, 
                 const double* parmsUpperIn, 
                 const unsigned int nEltsParmsFullIn, 
                 const bool standardizeParmsIn) :
			myRK4(myRK4In), 
			parmsFixed(parmsFixedIn),
			parmsFullInit(parmsFullInitIn),
			parmsLower(parmsLowerIn),
			parmsUpper(parmsUpperIn), 
			nEltsParmsFull(nEltsParmsFullIn), 
			standardizeParms(standardizeParmsIn) {}
	
	void completeParms(arma::vec parmsPart, double* parmsFull) {
		unsigned int compt=0;
		for (unsigned int it=0; it<nEltsParmsFull; it++) {
			if (parmsFixed[it]==true) {
				parmsFull[it] = parmsPart[compt];
				compt++;
			} else {
				parmsFull[it] = parmsFullInit[it];
			}
		}
	}
	
	void standardizeParmsPart(arma::vec& parmsPart) {
		unsigned int compt=0;
		for (unsigned int it=0; it<nEltsParmsFull; it++) {
			if (parmsFixed[it]==true) {
				parmsPart[compt]-=parmsLower[it];
				parmsPart[compt]/=(parmsUpper[it]-parmsLower[it]);
				compt++;
			}
		}
	}
	
	void unstandardizeParmsPart(arma::vec& parmsPart) {
		unsigned int compt=0;
		for (unsigned int it=0; it<nEltsParmsFull; it++) {
			if (parmsFixed[it]==true) {
				parmsPart[compt]*=(parmsUpper[it]-parmsLower[it]);
				parmsPart[compt]+=parmsLower[it];
				compt++;
			}
		}
	}
	
	double Evaluate(arma::vec armaParmsPart) {
		if (standardizeParms==true) unstandardizeParmsPart(armaParmsPart);
		double parmsFullArray[nEltsParmsFull];
		double* parmsFull = &parmsFullArray[0];
		completeParms(armaParmsPart, parmsFull);
		double out = myRK4->Evaluate(parmsFull);
		// add penalty
		bool outOfBounds = false;
		for (unsigned int it=0; it<nEltsParmsFull; it++) {
			if (parmsFull[it]<parmsLower[it] || parmsFull[it]>parmsUpper[it]) {
				outOfBounds = true;
			}
		}
		if (outOfBounds==true) { //arbitrary large value added as penalty if outOfBounds
			out+=1e50; 
		}
		return out;
	}

private:
RK4* myRK4;
const bool* parmsFixed;
const double* parmsFullInit;
const double* parmsLower;
const double* parmsUpper;
const unsigned int nEltsParmsFull;
const bool standardizeParms;
};

template<typename T>
void RcppListToPptr(Rcpp::List L, T**& pptr) {
	for (unsigned int it=0; it<L.size(); it++) {
		std::vector<double> tempVec = L[it];
		pptr[it] = (double*) malloc(sizeof(*pptr[it]) * tempVec.size());
		for (unsigned int it2=0; it2<tempVec.size(); it2++) {
			pptr[it][it2] = tempVec[it2];
		}
	}
}

// [[Rcpp::export]]
Rcpp::List minDistOptimForR(std::vector<double> Rparms,
                            std::vector<bool> RparmsFixed,
                            std::vector<double> RparmsLower,
                            std::vector<double> RparmsUpper,
                            std::vector<double> RyInit, 
                            unsigned int nt, 
                            double byt,
                            Rcpp::List RdataExogVar,
                            Rcpp::List RexogSamplingTime,
                            std::vector<unsigned int> RvarMinDist,
                            Rcpp::List RdataMinDist,
                            Rcpp::List RsamplingTime,
                            Rcpp::List RpointsWeight,
                            const unsigned int lambda,
                            double sigma, 
                            const unsigned int nIterMax, 
                            const double tol,
                            const bool standardizeParms,
                            const bool useParallel) {
	// INITIALIZATION //
	arma::vec RparmsInit(Rparms);
	double* parmsInit = &RparmsInit[0];
	double* parms = &Rparms[0];
	bool parmsFixedArray[RparmsFixed.size()];
	for (unsigned int it=0; it<RparmsFixed.size(); it++) {
		parmsFixedArray[it] = RparmsFixed[it];
	}
	bool* parmsFixed = &parmsFixedArray[0];
	double* parmsLower = &RparmsLower[0];
	double* parmsUpper = &RparmsUpper[0];
	unsigned int nEltsParmsPart = 0;
	arma::vec parmsPart(Rparms.size());
	for (unsigned int it=0; it<Rparms.size(); it++) {
		if (parmsFixed[it]==true) {
			parmsPart[nEltsParmsPart] = parms[it];
			nEltsParmsPart++;
		}
	}
	parmsPart.resize(nEltsParmsPart);
	
	double* yInit = &RyInit[0];
	unsigned int* varMinDist = &RvarMinDist[0];
	
	int nExogVar = RdataExogVar.size();
	double** dataExogVar = (double**) malloc(sizeof(double*)*RdataExogVar.size());
	double** exogSamplingTime = (double**) malloc(sizeof(double*)*RexogSamplingTime.size());
	RcppListToPptr(RdataExogVar, dataExogVar);
	RcppListToPptr(RexogSamplingTime, exogSamplingTime);

	double** dataMinDist = (double**) malloc(sizeof(double*)*RdataMinDist.size());
	double** samplingTime = (double**) malloc(sizeof(double*)*RsamplingTime.size());
	double** pointsWeight = (double**) malloc(sizeof(double*)*RpointsWeight.size());
	RcppListToPptr(RdataMinDist, dataMinDist);
	RcppListToPptr(RsamplingTime, samplingTime);
	RcppListToPptr(RpointsWeight, pointsWeight);
	// return Rcpp::List::create();
	minDistRK4 myMinDistRK4(yInit, nt, byt,
                         dataExogVar,
                         exogSamplingTime,
                         nExogVar,
                         varMinDist,
                         RvarMinDist.size(),
                         dataMinDist, 
                         samplingTime,
                         pointsWeight);
	
	fixedMinDistRK4<minDistRK4> myFixedMinDistRK4(&myMinDistRK4, 
										                            parmsFixed, 
										                            parms,
										                            parmsLower, 
										                            parmsUpper, 
										                            Rparms.size(), 
										                            standardizeParms);
	CMAES<fixedMinDistRK4<minDistRK4>> myCMAES(myFixedMinDistRK4,
                                           lambda,
                                           sigma, 
                                           nIterMax, 
                                           tol,
                                           useParallel);
	
	if (standardizeParms==true) myFixedMinDistRK4.standardizeParmsPart(parmsPart);
	double crit = myCMAES.Optimize(parmsPart);
	if (standardizeParms==true) myFixedMinDistRK4.unstandardizeParmsPart(parmsPart);
	// double crit = myFixedMinDistRK4.Evaluate(parmsPart);
	// double crit = 1;
	Rcpp::Rcout<<crit<<std::endl;
	unsigned int compt=0;
	for (unsigned int it=0; it<Rparms.size(); it++) {
		if (parmsFixed[it]==true) {
			Rparms[it]=parmsPart[compt];
			compt++;
		}
	}
	
	// FREE MEMORY
	for (unsigned int it=0; it<RdataMinDist.size(); it++) {
		free(dataMinDist[it]);
		free(samplingTime[it]);
		free(pointsWeight[it]);
	}
	free(dataMinDist);
	free(samplingTime);
	free(pointsWeight);
	
	return Rcpp::List::create(
		Named("criterion") = crit, 
		Named("parms") = Rparms);
}

// [[Rcpp::export]]
double ComputeDistanceForR(std::vector<double> Rparms, 
                          	std::vector<double> RyInit, 
                          	unsigned int nt, 
                          	double byt,
                          	Rcpp::List RdataExogVar,
                          	Rcpp::List RexogSamplingTime,
                          	std::vector<unsigned int> RvarMinDist,
                          	Rcpp::List RdataMinDist,
                          	Rcpp::List RsamplingTime,
                          	Rcpp::List RpointsWeight) {
	// INITIALIZATION //
	arma::vec RparmsInit(Rparms);
	double* parmsInit = &RparmsInit[0];
	double* parms = &Rparms[0];
	double* yInit = &RyInit[0];
	unsigned int* varMinDist = &RvarMinDist[0];
	
	double** dataExogVar = (double**) malloc(sizeof(double*)*RdataExogVar.size());
	double** exogSamplingTime = (double**) malloc(sizeof(double*)*RexogSamplingTime.size());
	RcppListToPptr(RdataExogVar, dataExogVar);
	RcppListToPptr(RexogSamplingTime, exogSamplingTime);
	int nExogVar = RdataExogVar.size();
	
	double** dataMinDist = (double**) malloc(sizeof(double*)*RdataMinDist.size());
	double** samplingTime = (double**) malloc(sizeof(double*)*RsamplingTime.size());
	double** pointsWeight = (double**) malloc(sizeof(double*)*RpointsWeight.size());
	RcppListToPptr(RdataMinDist, dataMinDist);
	RcppListToPptr(RsamplingTime, samplingTime);
	RcppListToPptr(RpointsWeight, pointsWeight);
	
	
	minDistRK4 myMinDistRK4(yInit, nt, byt,
                         dataExogVar,
                         exogSamplingTime,
                         nExogVar,
                         varMinDist,
                         RvarMinDist.size(),
                         dataMinDist, 
                         samplingTime,
                         pointsWeight);
	double out = myMinDistRK4.Evaluate(parms);
	for (unsigned int it=0; it<RdataMinDist.size(); it++) {
		free(dataMinDist[it]);
		free(samplingTime[it]);
		free(pointsWeight[it]);
	}
	free(dataMinDist);
	free(samplingTime);
	free(pointsWeight);
	return out;
}

// [[Rcpp::export]]
std::vector<double> getFulInitialPositionForR(std::vector<double> Rparms, 
                  							        			std::vector<double> Ry, 
                  							        			Rcpp::List RdataExogVar,
                  							        			Rcpp::List RexogSamplingTime) {
	// INITIALIZATION //
	double* parms = &Rparms[0];
	double* y = &Ry[0];
	
	double** dataExogVar = (double**) malloc(sizeof(double*)*RdataExogVar.size());
	double** exogSamplingTime = (double**) malloc(sizeof(double*)*RexogSamplingTime.size());
	RcppListToPptr(RdataExogVar, dataExogVar);
	RcppListToPptr(RexogSamplingTime, exogSamplingTime);
	int nExogVar = RdataExogVar.size();
	
	
	// init empty objects to initialize class (ugly but it works...)
	unsigned int nt;
	double byt;
	double** dataMinDist;
	double** samplingTime;
	double** pointsWeight;
	unsigned int* varMinDist;
	// END OF DATA CONVERSION
	
	minDistRK4 myMinDistRK4(y, nt, byt,
                         dataExogVar,
                         exogSamplingTime,
                         nExogVar,
                         varMinDist,
                         0,
                         dataMinDist, 
                         samplingTime,
                         pointsWeight);
	
	std::vector<double> yx(dim+dimIv);
	myMinDistRK4.getFullInitialPosition(0, y, parms, &yx[0], dataExogVar, exogSamplingTime, nExogVar);
	return yx;
}	

Rcpp::NumericMatrix RK4(int nt, 
                        double byT,
                        std::vector<double> Ry0,
                        std::vector<double> Rparms, 
                        double** dataExogVar,
                        double** exogSamplingTime, 
                        int nExogVar) {
	int it, it1;
	double *y = &Ry0[0];
	double *parms = &Rparms[0];
	double y1[dim], y2[dim], y3[dim], ydot0[dim], ydot1[dim], ydot2[dim], ydot3[dim], ydots[dim], x0[dimIv], x1[dimIv], x2[dimIv], x3[dimIv];
	Rcpp::NumericMatrix out(nt, dimOut);
	
	for (it=0; it<dim;it++) { //init out vector
		out(0, it)=y[it];
	}
	int comptExogVar[nExogVar];
	for (it=0; it<nExogVar; it++) comptExogVar[it]=1;
	// get intermediateVar and compute distance at t=0 //
	Func(0, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar); 
						for (it1=0; it1<dim; it1++) {
							out(0, dim+it1) = ydot0[it1];
						}
						for (it1=0; it1<dimIv; it1++) {
							out(0, 2*dim+it1) = x0[it1];
						}
						 
		
		for (it=0; it<(nt-1); it++) {
			
			
			
			
			Func(it*byT, y, parms, ydot0, x0, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++)
				y1[it1] = y[it1] + ydot0[it1]*0.5*byT;
			Func((it + 0.5)*byT, y1, parms, ydot1, x1, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++)
				y2[it1] = y[it1] + ydot1[it1]*0.5*byT;
			Func((it + 0.5)*byT, y2, parms, ydot2, x2, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++)
				y3[it1] = y[it1] + ydot2[it1]*byT;
			Func((it+1)*byT, y3, parms, ydot3, x3, dataExogVar, exogSamplingTime, nExogVar, comptExogVar);
			
			for (it1=0; it1<dim; it1++) {
				ydots[it1] = (ydot0[it1] + 2.0*ydot1[it1] + 2.0*ydot2[it1] + ydot3[it1])/6.0;
				y[it1] = y[it1] + byT*ydots[it1];
				out(it+1, it1) = y[it1];
			}
			
			for(it1=0;it1<dim;it1++){
							out(it+1, dim+it1) = ydots[it1];
						}
						for(it1=0;it1<dimIv;it1++){
							out(it+1, 2*dim+it1) = x0[it1];
						}
				
		}
		return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RK4(int nt, 
                        double byT,
                        std::vector<double> Ry0,
                        std::vector<double> Rparms, 
                        Rcpp::List RdataExogVar,
                        Rcpp::List RexogSamplingTime) {
	double** dataExogVar = (double**) malloc(sizeof(double*)*RdataExogVar.size());
	RcppListToPptr(RdataExogVar, dataExogVar);
	double** exogSamplingTime = (double**) malloc(sizeof(double*)*RexogSamplingTime.size());
	RcppListToPptr(RexogSamplingTime, exogSamplingTime);
	int nExogVar = RdataExogVar.size();
	Rcpp::NumericMatrix out = RK4(nt, byT, Ry0, Rparms, dataExogVar, exogSamplingTime, nExogVar);
	for (unsigned int it=0; it<RdataExogVar.size(); it++) {
		free(dataExogVar[it]);
		free(exogSamplingTime[it]);
	}
	free(dataExogVar);
	free(exogSamplingTime);
	
	return out;
}