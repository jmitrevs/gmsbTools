#ifndef USERANALYSISUTILS_EMJESFIX_HPP
#define USERANALYSISUTILS_EMJESFIX_HPP

///******************************************
///******************************************
///
/// A temporary routine for fixing null EMJES moment.
/// 
/// This is *temporary* and should *NOT* be applied on analysis files
/// built with athena > 15.6.8.9
///
/// Usage in C++
///   1) include this header in your analysis code
///   2) instantiate a EMJESFixer object once :
///  EMJESFixer jetEMJESfixer;
///
///   3) use it if the EMJES moment is null (this example for AntiKt4H1TopoJets)   :
///
///  double emjes = ...; // retrieve EMJES from jet or ntuple
///  if( emjes == 0 ){
///    emjes = jetEMJESfixer.fixAntiKt4H1Topo(jet_pt,jet_eta);
///  }
///  // use emjes as you wish
///
/// Usage in python
///  1) Load the macro and instantiate an object :
///
/// ROOT.gSystem.CompileMacro( 'EMJESfix.hpp')
/// jetEMJESfixer = ROOT.EMJESFixer()
///
/// 2) use it if the EMJES moment is null (this example for AntiKt4H1TopoJets)  
///
/// # retrieve EMJES from jet or ntuple in emjes
/// if emjes == 0:
///     emjes = jetEMJESfixer.fixAntiKt4H1Topo(jet_pt,jet_eta)
/// # use emjes as you wish
///
///******************************************

#include <cmath>

class EMJESFixer {

public:
  EMJESFixer(){
    fillConstants();
  }



  double fix(double& jet_pt, double& jet_eta, const double calibConstants[][4]) const;

  double fixAntiKt6H1Tower(double jet_pt, double jet_eta) const {
    return fix(jet_pt , jet_eta , m_AntiKt6H1Tower);
  }

  double fixAntiKt6H1Topo(double jet_pt, double jet_eta) const {
    return fix(jet_pt , jet_eta, m_AntiKt6H1Topo );
  }

  double fixAntiKt4H1Tower(double jet_pt, double jet_eta) const {
    return fix(jet_pt , jet_eta , m_AntiKt4H1Tower);
  }

  double fixAntiKt4H1Topo(double jet_pt, double jet_eta) const {
    return fix(jet_pt , jet_eta, m_AntiKt4H1Topo );
  }



private:
  double m_AntiKt6H1Tower[45][4];
  double m_AntiKt6H1Topo[45][4];
  double m_AntiKt4H1Tower[45][4];
  double m_AntiKt4H1Topo[45][4];




  void fillConstants();


};
#endif
