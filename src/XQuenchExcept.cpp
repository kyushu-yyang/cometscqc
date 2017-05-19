#include "XQuenchExcept.hpp"

XQuenchExcept::XQuenchExcept(const std::string& aMessage)
    : fMessage(aMessage)
{}


XQuenchExcept::~XQuenchExcept() throw()
{}


const char* XQuenchExcept::what() const throw()
{ return fMessage.c_str(); }

