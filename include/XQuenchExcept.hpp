/**
 *  @file   XQuenchLogger.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   1st Aug 2016
 */

#ifndef XQuenchExcept_HH
#define XQuenchExcept_HH

#include <string>
#include <exception>

/// class description:
/// exception class to handle the exception
///
class XQuenchExcept : public std::exception
{
  public:
    /*! constructor */
    XQuenchExcept(const std::string& aMessage);

    /*! deconstructor */
    virtual ~XQuenchExcept() throw();

    /*! send a message to exception */
    virtual void SetMessage(const std::string& message) { fMessage = message; }

    /*! return a exception message */
    virtual const char* what() const throw();    

  private:
    std::string fMessage;
};

#endif
