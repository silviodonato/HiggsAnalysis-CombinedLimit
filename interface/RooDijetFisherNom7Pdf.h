
#ifndef HiggsAnalysis_CombinedLimit_RooDijetFisherNom7Pdf_h
#define HiggsAnalysis_CombinedLimit_RooDijetFisherNom7Pdf_h
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
class RooDijetFisherNom7Pdf : public RooAbsPdf
{
public:
   RooDijetFisherNom7Pdf() {} ;
   RooDijetFisherNom7Pdf(const char *name, const char *title,
        RooAbsReal& _th1x, 
        RooAbsReal& _p1,
        RooAbsReal& _p2, 
        RooAbsReal& _p3, 
        RooAbsReal& _p4, 
        RooAbsReal& _p5,
        RooAbsReal& _p6,
        RooAbsReal& _sqrts);
   RooDijetFisherNom7Pdf(const RooDijetFisherNom7Pdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooDijetFisherNom7Pdf(*this,newname); }
   inline virtual ~RooDijetFisherNom7Pdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
   RooRealProxy p4;        // p4
   RooRealProxy p5;        // p5
   RooRealProxy p6;        // p6
   RooRealProxy sqrts;        // sqrts
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooDijetFisherNom7Pdf,1) // RazorDijetFisherNom7Pdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class DijetFisherNom7Function: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = exp(log(p[2] * x/p[0] - 1) + p[3]*x/p[0] + p[4]*pow(x/p[0],2) + p[5]*pow(x/p[0],3) + p[6]*pow(x/p[0],4))/pow(x/p[0],p[1]);
     return pdf;
   }
   double DoEval(double x) const
   {
     double pdf = exp(log(pars[2] * x/pars[0] - 1) + pars[3]*x/pars[0] + pars[4]*pow(x/pars[0],2) + pars[5]*pow(x/pars[0],3) + pars[6]*pow(x/pars[0],4))/pow(x/pars[0],pars[1]);
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new DijetFisherNom7Function();
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
      return 7;
   }
};
