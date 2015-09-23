#ifndef Pythia8_RopeUserHooks_H
#define Pythia8_RopeUserHooks_H

// Includes
#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "ParameterHandler.h"
#include <exception>

namespace Pythia8 {


class RopeUserHooks : public UserHooks {

// Convenient typedefs
typedef map<string,double> PytPars;

 public:
   // Exception class
  class RopeException: public exception {
  public:
    RopeException(string msg) : _msg(msg) { }
    virtual ~RopeException() throw() {}
    virtual const char* what() const throw(){
      return _msg.c_str();
    }

  private:
    string _msg;
  };


  public:


 RopeUserHooks() : _ph(NULL), _h(-1.0), _deltaY(0.5), _alpha(0.10), _veto(false), _event(NULL) {

  }

  ~RopeUserHooks() {
  }

  virtual bool canChangeFragPar() {
    return true;
  }

  virtual bool doChangeFragPar(StringFlav* flavPtr, StringZ* zPtr,
   StringPT * pTPtr, int endFlavour, double m2Had, vector<int> iParton) {

    // Get new parameters
    PytPars newPar;
    try{
      newPar = fetchParameters(endFlavour, m2Had, iParton); 
      if(newPar.find("null") != newPar.end())
        throw RopeException("Failed to fetch parameters! No Ropes.");
      // Change settings to new settings
      for(PytPars::iterator itr = newPar.begin(); itr!=newPar.end(); ++itr) 
        settingsPtr->parm(itr->first,itr->second);
      // Re-initialize flavour, z, and pT selection with new settings
      flavPtr->init(*settingsPtr,rndmPtr);
      zPtr->init(*settingsPtr,*particleDataPtr,rndmPtr);
      pTPtr->init(*settingsPtr,*particleDataPtr,rndmPtr);
    }
    catch(RopeException& e){
      cerr << e.what() << endl;
      return false;
    }
    return true;
  }

  virtual bool doVetoFragmentation(Particle p){
    // This is just an example of a veto.
    // Really this should not be implemented,
    // as it introduces a horrible bias against
    // high-pT particles!

    if(p.pT() > 5.0 && _enh > 1.0){
      _veto = true;
      return true;
    }
    return false;

  }

  void setParameterHandler(ParameterHandler * ph){
    _ph = ph;
  }

  void setEvent(Event& event){
    _event = &event;
  }

  // Update event, needs to be done from main file
  void newEvent(){
    _hadronized.clear();
  }

  // Set enhancement manually from main file
  // Will overrule others if set.
  void setEnhancement(double h){
    _h = h;
  }

  // Set the rap. span around a break to look for 
  // overlapping gluons
  void setRapiditySpan(double dy){
    _deltaY = dy;
  }

  // Set the ratio of string radius to
  // radius of a pp collision 
  void setRadiusRatio(double alpha){
    _alpha = alpha;
  }

 private:

  PytPars fetchParameters(int endFlavour, double m2Had, vector<int> iParton){
    if(!_event)
      throw RopeException("Can't find event!");

    if(!(_ph))
      throw RopeException("Missing parameter handler!");

    if(_veto){
      _veto = false;
      return _ph->GetEffectiveParameters(1.0);
    }

    if(find(_hadronized.begin(),_hadronized.end(),*iParton.begin()) == _hadronized.end()){
      _hadronized.reserve(_hadronized.size() + iParton.size());
      _hadronized.insert(_hadronized.end(),iParton.begin(),iParton.end());
    }

    // If we have manually set enhancement, we should just use that, and avoid complicated behavior.
    if(_h > 0) return _ph->GetEffectiveParameters(_h); 

   // Quark string ends, default mode
   if (endFlavour != 21){
    // Test consistency
    if(_event->at(*(iParton.begin())).id() != endFlavour && 
        _event->at(*(iParton.end() - 1)).id() != endFlavour)
      throw RopeException("Quark end inconsistency!");
   
      // First we must let the string vector point in the right direction
      if(_event->at(*(iParton.begin())).id() != endFlavour)
       reverse(iParton.begin(),iParton.end());

      // Initialize a bit
      Vec4 hadronic4Momentum(0,0,0,0);
      _enh = 1.0;
      double dipFrac;
      vector<int>::iterator dipItr;
      // Find out when invariant mass exceeds m2Had
      for(dipItr = iParton.begin(); dipItr != iParton.end(); ++dipItr){
       double m2Big = hadronic4Momentum.m2Calc();
        if( m2Had <= m2Big){
          // Approximate the fraction we are in on the dipole, this goes in three cases
          // We are at the beginning
          if(m2Had == 0){
            dipFrac = 0;
          }
          // We are somewhere in the first dipole
          else if(dipItr - 1  == iParton.begin()){
            dipFrac = sqrt(m2Had/m2Big);
          }
          else{
            if(_event->at(*(dipItr - 1)).id() != 21)
              throw RopeException("Connecting partons should always be gluons");
    
            hadronic4Momentum -= 0.5*_event->at(*(dipItr -1)).p();
            double m2Small = hadronic4Momentum.m2Calc();  
      
            dipFrac = (sqrt(m2Had) - sqrt(m2Small)) / (sqrt(m2Big) - sqrt(m2Small));
          }
          break;
        }
        hadronic4Momentum += _event->at(*dipItr).id() == 21 ? 0.5*_event->at(*dipItr).p() : _event->at(*dipItr).p();
      }
      // If we reached the end
      // we are in a small string that should just be collapsed
      if(dipItr == iParton.end())
        return _ph->GetEffectiveParameters(1.0);
      // Sanity check
      if(dipFrac < 0 || dipFrac > 1)
        throw RopeException("We can never be less than 0 or more than 100 perc. in on a dipole.");
      // We now figure out at what rapidity value,
      // in the lab system, the string is breaking
      double yBreak;
      // Trivial case, just inherit
      if(dipFrac == 0) 
        yBreak = _event->at(*dipItr).y();
      else{
        // Sanity check
        if(dipItr == iParton.begin())
          throw RopeException("We are somehow before the first dipole.");
        double dy = _event->at(*dipItr).y() - _event->at(*(dipItr - 1)).y();
        yBreak = _event->at(*(dipItr - 1)).y() + dipFrac*dy;
      }
      // Count the number of partons in the whole
      // event within deltay of breaking point
      double p = 1;
      double q = 0;

      for(int i = 0; i < _event->size(); ++i){
        // Don't double count partons from this
        // string (no self-overlap)
        if(find(iParton.begin(),iParton.end(), i) != iParton.end())
          continue;
        // Don't count strings that are already hadronized
        if(find(_hadronized.begin(),_hadronized.end(),i) != _hadronized.end())
          continue;
        double pRap = _event->at(i).y();
        if(pRap > yBreak - _deltaY && pRap < yBreak + _deltaY ){
          // Do a "Buffon" selection to decide whether
          // two strings overlap in impact parameter space
          // given ratio of string diameter to collision diameter (_alpha)
          double r1 = rndmPtr->flat();
          double r2 = rndmPtr->flat();
          double theta1 = 2*M_PI*rndmPtr->flat();
          double theta2 = 2*M_PI*rndmPtr->flat();
          // Overlap?
          if(4*_alpha*_alpha > pow2(sqrt(r1)*cos(theta1) - sqrt(r2)*cos(theta2)) +
            pow2(sqrt(r1)*sin(theta1) - sqrt(r2)*sin(theta2))){
            if(rndmPtr->flat() < 0.5) p += 0.5;
            else q += 0.5;
          } 
        }
      }
      _enh = 0.25*(2.0*p+q+2.0);

      return _ph->GetEffectiveParameters(_enh);
    }
   // For closed gluon loops we cannot distinguish the ends.
   // Do an average instead
   else{
      return _ph->GetEffectiveParameters(1.0);
   }
  
    return _ph->GetEffectiveParameters(1.0);
}
  ParameterHandler* _ph;
  double _h, _deltaY, _alpha, _enh;
  bool _veto;
  Event*  _event;
  vector<int> _hadronized;
};

}
#endif
