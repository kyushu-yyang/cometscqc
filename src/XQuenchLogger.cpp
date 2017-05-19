#include <iostream>
#include "XQuenchLogger.hpp"

using Quench::XQuenchLogger;

XQuenchLogger* XQuenchLogger::instance = NULL;

XQuenchLogger::XQuenchLogger()
    : fActive(false)
{}

XQuenchLogger* XQuenchLogger :: GetInstance()
{
  if (instance==NULL)
    instance = new XQuenchLogger();

  return instance;
}

void XQuenchLogger :: Start(Level minpriority, const std::string& filename)
{
  instance->fActive = true;
  instance->fPriority = minpriority;

  if ( filename!=" " ) 
    instance->fFile.open( filename.c_str() );
}

void XQuenchLogger :: Stop()
{
  instance->fActive = false;
  
  if ( instance->fFile.is_open() )
    instance->fFile.close();
}

void XQuenchLogger :: Write(Level priority, const std::string& message)
{
  if ( instance->fActive && priority>=instance->fPriority ) {
    std::ostream& stream = instance->fFile.is_open() ? instance->fFile : std::cout;

    stream << "["
           << instance->GetPriorityName(priority)
           << "]"
           << ": "
           << message;
  }
}

std::ostream& XQuenchLogger :: GetLogStream(Level priority)
{
  if (!instance->fFile.is_open())
    return std::cout;

  instance->fFile << "["
                  << instance->GetPriorityName(priority)
                  << "]"
                  << ": ";

  return instance->fFile;
}

const std::string XQuenchLogger :: GetPriorityName(Level level)
{
  const std::string priority[5] = {
      "DEBUG",
      "CONFIG",
      "INFO",
      "WARNING",
      "ERROR"
  };

  switch (level) {
    case   DEBUG: return priority[0]; break;
    case  CONFIG: return priority[1]; break;
    case    INFO: return priority[2]; break;
    case WARNING: return priority[3]; break;
    case   ERROR: return priority[4]; break;
    default: return ""; break;
  }
}
