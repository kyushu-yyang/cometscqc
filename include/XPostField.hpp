/**
 *  @file   XPostField.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   3 Aug 2016
 */

#ifndef XPostField_HH
#define XPostField_HH

#include <vector>
#include <string>

#ifndef XFieldHandle_HH
#include "XFieldHandle.hpp"
#endif

class TH2F;

/// class to handle the plot of magnetic field
//
class XPostField
{
  public:
    /*! constructor */
    XPostField();

    /*! constructor */
    XPostField(Quench::XFieldHandle* fld);

    /*! deconstructor */
    ~XPostField();

    /*! @brief setup field handler */
    void SetFieldHandler(Quench::XFieldHandle* fld);

    /*! @brief plot */
    void Plot();

  protected:
    /*! @brief plot the 2d hist */
    void plot2d(TH2F* hist);

  private:
    std::string fMagnet;
    Quench::XMagnetInfoContainer* fInfo;
    std::vector<Quench::XFieldContainer*> fCollect;
};

#endif
