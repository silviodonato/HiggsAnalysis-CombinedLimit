
#ifndef HiggsAnalysis_CombinedLimit_RooDijetTestPdf_h
#define HiggsAnalysis_CombinedLimit_RooDijetTestPdf_h
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
class RooDijetTestPdf : public RooAbsPdf
{
public:
   RooDijetTestPdf() {} ;
   RooDijetTestPdf(const char *name, const char *title,
        RooAbsReal& _th1x, 
        RooAbsReal& _p1,
        RooAbsReal& _p2, 
        RooAbsReal& _p3, 
//        RooAbsReal& _p4, 
//        RooAbsReal& _p5,
//        RooAbsReal& _p6,
        RooAbsReal& _sqrts);
   RooDijetTestPdf(const RooDijetTestPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooDijetTestPdf(*this,newname); }
   inline virtual ~RooDijetTestPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
//   RooRealProxy p4;        // p4
//   RooRealProxy p5;        // p5
//   RooRealProxy p6;        // p6
   RooRealProxy sqrts;        // sqrts
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooDijetTestPdf,1) // RazorDijetTestPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class DijetTestFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = exp(p[2]*x/p[0]) * (- 1 + x/p[0]*p[3]) /pow(x/p[0],p[1]);
     return pdf;
   }
   double DoEval(double x) const
   {
     double pdf = exp(pars[2]*x/pars[0]) * (- 1 + x/pars[0]*pars[3]) /pow(x/pars[0],pars[1]);
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new DijetTestFunction();
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
      return 4;
   }
};
