#ifndef LikelihoodSwitches_H
#define LikelihoodSwitches_H

struct LikelihoodSwitches {

  LikelihoodSwitches () :
    m_useEoverP (false) ,
    m_useOneOverEMinusOneOverP (false) ,
    m_useDeltaEta (true) ,
    m_useDeltaPhi (true) ,
    m_useHoverE (true) ,
    m_useFBrem (true) ,
    m_useSigmaEtaEta (true) ,
    m_useSigmaPhiPhi (false) {} ;

  bool m_useEoverP ;
  bool m_useOneOverEMinusOneOverP ;
  bool m_useDeltaEta ;
  bool m_useDeltaPhi ;
  bool m_useHoverE ;
  bool m_useFBrem ;
  bool m_useSigmaEtaEta ;
  bool m_useSigmaPhiPhi ;

};

#endif // LikelihoodSwitches_H
