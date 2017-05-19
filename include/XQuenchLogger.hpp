/**
 *  @file   XQuenchLogger.hpp
 *  @author Y.Yang (Kyushu University)
 *  @date   1st. Aug. 2016
 */

#ifndef XQuenchLogger_HH
#define XQuenchLogger_HH

#include <iomanip>
#include <string>
#include <fstream>

namespace Quench
{ class XQuenchLogger; }

/// class description:
/// logger class is a singleton class to log the error and warning
/// 
class Quench::XQuenchLogger
{
  public:
    /** @enum Level
     *  logging level
     */
    enum Level {
      DEBUG,
      CONFIG,
      INFO,
      WARNING,
      ERROR
    };

    /// @fn get instance 
    /// @return static logging class
    static XQuenchLogger* GetInstance();

    /// @brief start loging and write message to the file.
    static void Start(Level minpriority, const std::string& filename);

    /*! stop the logging */
    static void Stop();

    /*! write message */
    static void Write(Level priority, const std::string& message);

    /*! return ostream of file */
    static std::ostream& GetLogStream(Level priority);


  private:
    /*! constructor */
    XQuenchLogger();

    /*! constructor */
    XQuenchLogger(const XQuenchLogger& log) {}

    /*! = operator */
    //XQuenchLogger& operator = (const XQuenchLogger& log) {}

    /*! static instance of this class */
    static XQuenchLogger* instance;

    /*! get the priority names */
    static const std::string GetPriorityName(Level level);

  private:
    bool fActive;
    std::ofstream fFile;
    Level fPriority;

};


/// @def a macro to handle the output of error message
/*
#ifndef _QUENCH_OUTPUT_ERROR
#define _QUENCH_OUTPUT_ERROR(level, message)                               \
    do {                                                                   \
      XQuenchLogger* log = XQuenchLogger::GetInstance();                   \
      log->GetLogStream(level) << message << "  "                          \
          << " ||    file: " << __FILE__                                   \
          << " - line: " << __LINE__                                       \
          << std::setprecision(6) << std::setw(0)                          \
          << std::setfill(' ')    << std::endl;                            \
    } while(0)
#endif
*/
#ifndef _QUENCH_OUTPUT_ERROR
#define _QUENCH_OUTPUT_ERROR(level, message)                               \
    do {                                                                   \
      XQuenchLogger* log = XQuenchLogger::GetInstance();                   \
      log->GetLogStream(level) << message << "  "                          \
          << std::setprecision(6) << std::setw(0)                          \
          << std::setfill(' ')    << std::endl;                            \
    } while(0)
#endif

/// @def turn on or off the error output
#ifndef QUENCH_ERROR_OUTPUT
#define QUENCH_ERROR_OUTPUT true
#endif

/// @def error output
#ifndef QuenchError
#define QuenchError(level, message)            \
    do {                                       \
      if (QUENCH_ERROR_OUTPUT) {               \
        _QUENCH_OUTPUT_ERROR(level, message);  \
      }                                        \
    } while(0)
#endif

/// @def logger for debug level
#ifndef QuenchDebug
#define QuenhcDebug(message)                        \
    do {                                            \
      QuenchError( XQuenchLogger::DEBUG, message);  \
    } while(0)
#endif

/// @def logger for fatal error level
#ifndef QuenchFatal
#define QuenchFatal(message)                       \
    do {                                           \
      QuenchError( XQuenchLogger::ERROR, message); \
    } while(0)
#endif

/// @def logger for warning level
#ifndef QuenchWarning
#define QuenchWarning(message)                        \
    do {                                              \
      QuenchError( XQuenchLogger::WARNING, message);  \
    } while(0)
#endif

/// @def logger for info level
#ifndef QuenchInfo
#define QuenchInfo(message)                        \
    do {                                           \
      QuenchError( XQuenchLogger::INFO, message);  \
    } while(0)
#endif


#endif
