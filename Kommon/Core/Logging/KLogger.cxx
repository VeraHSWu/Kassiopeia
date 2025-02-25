/*
 * KLogger.cxx
 *
 *  Created on: 18.11.2011
 *      Author: Marco Kleesiek <marco.kleesiek@kit.edu>
 */

#include "KLogger.h"

#include "KException.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace katrin;

static const char* skEndColor = COLOR_PREFIX COLOR_NORMAL COLOR_SUFFIX;
static const char* skFatalColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_PURPLE COLOR_SUFFIX;
static const char* skErrorColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_RED COLOR_SUFFIX;
static const char* skWarnColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_YELLOW COLOR_SUFFIX;
static const char* skInfoColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_GREEN COLOR_SUFFIX;
static const char* skDebugColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_CYAN COLOR_SUFFIX;
static const char* skTraceColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_BLUE COLOR_SUFFIX;
static const char* skOtherColor = COLOR_PREFIX COLOR_BRIGHT COLOR_SEPARATOR COLOR_FOREGROUND_WHITE COLOR_SUFFIX;

namespace
{

inline const char* level2Color(KLogger::ELevel level)
{
    switch (level) {
        case KLogger::eFatal:
            return skFatalColor;
        case KLogger::eError:
            return skErrorColor;
        case KLogger::eWarn:
            return skWarnColor;
        case KLogger::eInfo:
            return skInfoColor;
        case KLogger::eDebug:
            return skDebugColor;
        case KLogger::eTrace:
            return skTraceColor;
        default:
            return skOtherColor;
    }
}

}  // namespace

#ifdef LOG4CXX

/*
 * Default implementation for systems with the 'log4cxx' library installed.
 */

#include <log4cxx/basicconfigurator.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/level.h>
#include <log4cxx/logger.h>
#include <log4cxx/logmanager.h>
#include <log4cxx/logstring.h>
#include <log4cxx/patternlayout.h>
#include <log4cxx/propertyconfigurator.h>

using namespace log4cxx;

namespace
{

inline LevelPtr level2Ptr(KLogger::ELevel level)
{
    switch (level) {
        case KLogger::eTrace:
            return Level::getTrace();
        case KLogger::eDebug:
            return Level::getDebug();
        case KLogger::eInfo:
            return Level::getInfo();
        case KLogger::eWarn:
            return Level::getWarn();
        case KLogger::eError:
            return Level::getError();
        case KLogger::eFatal:
            return Level::getFatal();
        default:
            return Level::getOff();
    }
}

inline KLogger::ELevel ptr2Level(LevelPtr ptr)
{
    if (!ptr)
        return KLogger::eUndefined;
    switch (ptr->toInt()) {
        case Level::TRACE_INT:
            return KLogger::eTrace;
        case Level::DEBUG_INT:
            return KLogger::eDebug;
        case Level::INFO_INT:
            return KLogger::eInfo;
        case Level::WARN_INT:
            return KLogger::eWarn;
        case Level::ERROR_INT:
            return KLogger::eError;
        case Level::FATAL_INT:
            return KLogger::eFatal;
        default:
            return KLogger::eUndefined;
    }
}

}  // namespace

#ifndef _LOG4CXX_COLORED_PATTERN_LAYOUT_H
#define _LOG4CXX_COLORED_PATTERN_LAYOUT_H

namespace log4cxx
{

class LOG4CXX_EXPORT ColoredPatternLayout : public PatternLayout
{
  public:
    DECLARE_LOG4CXX_OBJECT(ColoredPatternLayout)
    BEGIN_LOG4CXX_CAST_MAP()
    LOG4CXX_CAST_ENTRY(ColoredPatternLayout)
    LOG4CXX_CAST_ENTRY_CHAIN(Layout)  // NOLINT
    END_LOG4CXX_CAST_MAP()

    ColoredPatternLayout() : PatternLayout() {}
    ColoredPatternLayout(const LogString& pattern) : PatternLayout(pattern) {}
    ~ColoredPatternLayout() override = default;

  protected:
    void format(LogString& output, const spi::LoggingEventPtr& event, helpers::Pool& pool) const override;
};

LOG4CXX_PTR_DEF(ColoredPatternLayout);

}  // namespace log4cxx

#endif /* _LOG4CXX_COLORED_PATTERN_LAYOUT_H */

IMPLEMENT_LOG4CXX_OBJECT(ColoredPatternLayout)

void ColoredPatternLayout::format(LogString& output, const spi::LoggingEventPtr& event, helpers::Pool& pool) const
{
    PatternLayout::format(output, event, pool);
    output = level2Color(ptr2Level(event->getLevel())) + output + skEndColor;
    return;
}


struct StaticInitializer
{
    StaticInitializer()
    {

        if (LogManager::getLoggerRepository()->isConfigured())
            return;

        char* envLoggerConfig;
        envLoggerConfig = getenv("LOGGER_CONFIGURATION");
        if (envLoggerConfig != nullptr) {
            PropertyConfigurator::configure(envLoggerConfig);
            //            #ifndef NDEBUG
            //                std::cout << "Logger configuration: " << envLoggerConfig << std::endl;
            //            #endif
        }
        else {
#ifdef LOGGER_CONFIGURATION
            PropertyConfigurator::configure(TOSTRING(LOGGER_CONFIGURATION));
#else
            LogManager::getLoggerRepository()->setConfigured(true);
            LoggerPtr root = Logger::getRootLogger();
#ifdef NDEBUG
            Logger::getRootLogger()->setLevel(Level::getInfo());
#endif
            const LogString TTCC_CONVERSION_PATTERN(LOG4CXX_STR("%r [%-5p] %16c: %m%n"));
            LayoutPtr layout(new PatternLayout(TTCC_CONVERSION_PATTERN));
            AppenderPtr appender(new ConsoleAppender(layout));
            root->addAppender(appender);
#endif
        }
    }
} static sLoggerInitializer;

struct KLogger::Private
{
    void log(const LevelPtr& level, const string& message, const Location& loc)
    {
        fLogger->forcedLog(level,
                           message,
                           ::log4cxx::spi::LocationInfo(loc.fFileName, loc.fFunctionName, loc.fLineNumber));
    }

    LoggerPtr fLogger;
    string fLoggerName;
};

KLogger::KLogger(const char* name) : fPrivate(new Private())
{
    fPrivate->fLoggerName = (name == nullptr) ? "root" : name;
    fPrivate->fLogger = (name == nullptr) ? Logger::getRootLogger() : Logger::getLogger(name);

    KLoggerTable::GetInstance().Add(this);
}

KLogger::KLogger(const std::string& name) : fPrivate(new Private())
{
    fPrivate->fLoggerName = name;
    fPrivate->fLogger = Logger::getLogger(name);

    KLoggerTable::GetInstance().Add(this);
}

KLogger::~KLogger()
{
    if (KLoggerTable::IsInitialized())
        KLoggerTable::GetInstance().Remove(this);

    if (fPrivate)
        //delete fPrivate;  // BUG: should be deleted, but causes segfault in log4cxx:~Logger()
        fPrivate = nullptr;
}

const string& KLogger::GetName() const
{
    return fPrivate->fLoggerName;
}

bool KLogger::IsLevelEnabled(ELevel level) const
{
    ELevel displayedLevel = KLoggerTable::GetInstance().CorrectedLevel(level);
    return fPrivate->fLogger->isEnabledFor(level2Ptr(displayedLevel)) || level >= eFatal;
}

KLogger::ELevel KLogger::GetLevel() const
{
    return ptr2Level(fPrivate->fLogger->getLevel());
}

void KLogger::SetLevel(ELevel level)
{
    fPrivate->fLogger->setLevel(level2Ptr(level));
}

void KLogger::Log(ELevel level, const string& message, const Location& loc)
{
    fPrivate->log(level2Ptr(level), message, loc);

    if (level >= eFatal) {
        throw KException() << "fatal error in " << loc.fFunctionName;
        exit(EXIT_FAILURE);
    }
#ifdef KLOGGER_THROW_EXCEPTIONS
    else if (level >= ELevel::eError)
        throw KException() << "error in " << loc.fFunctionName;
#endif
}


#else  // LOG4CXX

/**
 * Fallback solution for systems without log4cxx.
 */

#include <iomanip>
#include <map>
#include <utility>

static const bool skColored = true;

namespace
{

const char* level2Str(KLogger::ELevel level)
{
    switch (level) {
        case KLogger::eTrace:
            return "TRACE";
        case KLogger::eDebug:
            return "DEBUG";
        case KLogger::eInfo:
            return "INFO";
        case KLogger::eWarn:
            return "WARN";
        case KLogger::eError:
            return "ERROR";
        case KLogger::eFatal:
            return "FATAL";
        default:
            return "XXX";
    }
}

void printTime(ostream& strm)
{
    using namespace std::chrono;

    const auto now = system_clock::now();
    const time_t cTime = system_clock::to_time_t(now);

    auto duration = now.time_since_epoch();
    duration -= duration_cast<seconds>(duration);

    /*
     * Unfortunately, g++ < 5.0 does not implement std::put_time, so I have to
     * resort to strftime at this point:
     */
    char dateTimeStr[24];
    strftime(dateTimeStr, sizeof(dateTimeStr), "%F %T", localtime(&cTime));
    strm << dateTimeStr;
    strm << "." << setfill('0') << setw(3) << duration_cast<milliseconds>(duration).count();
    strm << setfill(' ');
}

}  // namespace

struct KLogger::Private
{
    Private(string name) :
        fLoggerName(name),
#ifdef NDEBUG
        fLogLevel(eInfo)
#else
        fLogLevel(eDebug)
#endif
    {}

    std::string fLoggerName;
    KLogger::ELevel fLogLevel;

    void logCout(const char* level, const string& message, const Location& /*loc*/, const char* color = skOtherColor)
    {
        if (skColored) {
            cout << color;
            printTime(cout);
            cout << " [" << setw(5) << level << "] " << setw(10) << fLoggerName << ": " << message;
            cout << skEndColor;
            cout << endl;
        }
        else {
            printTime(cout);
            cout << " [" << setw(5) << level << "] " << setw(10) << fLoggerName << ": " << message;
            cout << endl;
        }
    }

    void logCerr(const char* level, const string& message, const Location& /*loc*/, const char* color = skOtherColor)
    {
        if (skColored) {
            cerr << color;
            printTime(cerr);
            cerr << " [" << setw(5) << level << "] " << setw(10) << fLoggerName << ": " << message;
            cerr << skEndColor;
            cerr << endl;
        }
        else {
            printTime(cerr);
            cerr << " [" << setw(5) << level << "] " << setw(10) << fLoggerName << ": " << message;
            cerr << endl;
        }
    }
};

KLogger::KLogger(const char* name) : fPrivate(new Private(name == nullptr ? "root" : name))
{
    KLoggerTable::GetInstance().Add(this);
}

KLogger::KLogger(const std::string& name) : fPrivate(new Private(name.empty() ? "root" : name))
{
    KLoggerTable::GetInstance().Add(this);
}

KLogger::~KLogger()
{
    if (KLoggerTable::IsInitialized())
        KLoggerTable::GetInstance().Remove(this);

    delete fPrivate;
}

const string& KLogger::GetName() const
{
    return fPrivate->fLoggerName;
}

bool KLogger::IsLevelEnabled(ELevel level) const
{
    ELevel displayedLevel = KLoggerTable::GetInstance().CorrectedLevel(level);
    return fPrivate->fLogLevel <= displayedLevel || level >= eFatal;
}

KLogger::ELevel KLogger::GetLevel() const
{
    return fPrivate->fLogLevel;
}

void KLogger::SetLevel(ELevel level)
{
    fPrivate->fLogLevel = level;
}

void KLogger::Log(ELevel level, const string& message, const Location& loc)
{
    const char* levelStr = level2Str(level);
    const char* color = level2Color(level);

#if 0
    if (level >= eError)
        fPrivate->logCerr(levelStr, message, loc, color);
    else
        fPrivate->logCout(levelStr, message, loc, color);
#else
    fPrivate->logCerr(levelStr, message, loc, color);
#endif

    if (level >= eFatal) {
        throw KException() << "fatal error in " << loc.fFunctionName;
        exit(EXIT_FAILURE);
    }
#ifdef KLOGGER_THROW_EXCEPTIONS
    else if (level >= eError)
        throw KException() << "error in " << loc.fFunctionName;
#endif
}
#endif  // LOG4CXX

KLoggerTable::KLoggerTable() : fLoggerMap(), fVerbosityLevel(0) {}

KLoggerTable::~KLoggerTable() = default;

void KLoggerTable::Add(KLogger* logger)
{
    const string& name = logger->GetName();
    fLoggerMap[name].push_front(logger);
}

list<KLogger*> KLoggerTable::Get(const string& name)
{
    auto mapIt = fLoggerMap.find(name);
    if (mapIt != fLoggerMap.end())
        return mapIt->second;
    return {};
}

void KLoggerTable::Remove(const string& name)
{
    fLoggerMap.erase(name);
}

void KLoggerTable::Remove(KLogger* logger)
{
    for (auto & mapIt : fLoggerMap) {
        mapIt.second.remove(logger);
    }
}

void KLoggerTable::SetLevel(const KLogger::ELevel& level)
{
    for (auto& mapIt : fLoggerMap) {
        for (auto& logger : mapIt.second)
            logger->SetLevel(level);
    }
}

void KLoggerTable::SetLevel(const KLogger::ELevel& level, const string& name)
{
    auto mapIt = fLoggerMap.find(name);
    if (mapIt != fLoggerMap.end()) {
        for (auto& logger : mapIt->second)
            logger->SetLevel(level);
    }
}

void KLoggerTable::SetVerbosityLevel(int level)
{
    fVerbosityLevel = level;
}

KLogger::ELevel KLoggerTable::CorrectedLevel(KLogger::ELevel level) const
{
    // larger value means higher message severity
    auto result = static_cast<KLogger::ELevel>(level + fVerbosityLevel);
    if (result <= KLogger::eTrace)
        return KLogger::eTrace;
    else if (result >= KLogger::eFatal)
        return KLogger::eFatal;
    return result;
}
