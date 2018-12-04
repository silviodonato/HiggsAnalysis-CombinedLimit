//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooDijet6ParamSilvioPdf_h
#define HiggsAnalysis_CombinedLimit_RooDijet6ParamSilvioPdf_h
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooDijet6ParamSilvioPdf : public RooAbsPdf
{
public:
   RooDijet6ParamSilvioPdf() {} ;
   RooDijet6ParamSilvioPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		    RooAbsReal& _p2, RooAbsReal& _p3, RooAbsReal& _p4, RooAbsReal& _p5,
		    RooAbsReal& _sqrts, RooAbsReal& _meff, RooAbsReal& _seff);
   RooDijet6ParamSilvioPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		    RooAbsReal& _p2, RooAbsReal& _p3, RooAbsReal& _p4, RooAbsReal& _p5,
		    RooAbsReal& _sqrts);
   RooDijet6ParamSilvioPdf(const RooDijet6ParamSilvioPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooDijet6ParamSilvioPdf(*this,newname); }
   inline virtual ~RooDijet6ParamSilvioPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
   RooRealProxy p4;        // p4
   RooRealProxy p5;        // p5
   RooRealProxy sqrts;        // sqrts
   RooRealProxy meff;        // meff
   RooRealProxy seff;        // seff
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooDijet6ParamSilvioPdf,1) // RazorDijet6ParamSilvioPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class Dijet6ParamSilvioFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
//     double pdf = pow(1-x/p[0],p[1])/pow(x/p[0],p[2]+p[3]*log(x/p[0])+p[4]*log(x/p[0])*log(x/p[0])+p[5]*log(x/p[0])*log(x/p[0])*log(x/p[0]));
//     double pdf = 1 / pow(x + p[2] , 5 + p[3] * log(x/p[0]) ) / ( exp(p[1]/(x+p[3])) - 1) * (1+TMath::Erf((x - p[4])/p[5]));


//     double pdf = exp(-p[2]*(x/p[0])-p[3]*pow(x/p[0],2))/pow(x/p[0],p[1])  * (1. - p[5]/TMath::Max((x-p[4]),1.));
     double pdf = exp(-p[2]*(x/p[0])-p[3]*pow(x/p[0],2))/pow(x/p[0],p[1])  * (1. - p[5]/(pow(x-p[4],2)));
//     double pdf = exp(-p[2]*(x/p[0])-p[3]*pow(x/p[0],2))/pow(x/p[0],p[1] + p[5] * log(x/p[0]))  / (exp(1/TMath::Max(x-p[4],1e-9))-1);

     double eff = 1.;
     //if (p[6]>0) eff = 0.5 * (1.0 + TMath::Erf((x - p[5])/p[6])) ; // Error function
     if (p[7]>0) eff = 1.0/(1.0 + exp(-2.4*(x - p[6])/p[7])) ; // Sigmoid function
     return pdf*eff;
   }
   
   double DoEval(double x) const
   {
//     double pdf = pow(1-x/pars[0],pars[1])/pow(x/pars[0],pars[2]+pars[3]*log(x/pars[0])+pars[4]*log(x/pars[0])*log(x/pars[0])+pars[5]*log(x/pars[0])*log(x/pars[0])*log(x/pars[0]));
//     double pdf = 1 / pow(x + pars[2] , 5 + pars[3] * log(x/pars[0]) ) / ( exp(pars[1]/(x+pars[3])) - 1) * (1+TMath::Erf((x - pars[4])/pars[5]));

//     double pdf = exp(-pars[2]*(x/pars[0])-pars[3]*pow(x/pars[0],2))/pow(x/pars[0],pars[1])  * (1. - pars[5]/TMath::Max((x-pars[4]),1.));
     double pdf = exp(-pars[2]*(x/pars[0])-pars[3]*pow(x/pars[0],2))/pow(x/pars[0],pars[1])  * (1. - pars[5]/(pow(x-pars[4],2)));
//     double pdf = exp(-pars[2]*(x/pars[0])-pars[3]*pow(x/pars[0],2))/pow(x/pars[0],pars[1] + pars[5] * log(x/pars[0]) )  / (exp(1/TMath::Max(x-pars[4],1e-9))-1);

     double eff = 1.;
     //if (pars[6]>0) eff = 0.5 * (1.0 + TMath::Erf((x - pars[5])/pars[6])); // Error function     
     if (pars[7]>0) eff = 1.0/(1.0 + exp(-2.4*(x - pars[6])/pars[7])); // Sigmoid function
     return pdf*eff;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new Dijet6ParamSilvioFunction();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 8;
   }
};
