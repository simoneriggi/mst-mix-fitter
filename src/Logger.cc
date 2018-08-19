#include <Logger.h>

ClassImp(MSTMixFitter_ns::Logger)
ClassImp(MSTMixFitter_ns::ConsoleLogger)
ClassImp(MSTMixFitter_ns::LoggerManager)

namespace MSTMixFitter_ns {

int LoggerManager::m_target;
Logger* LoggerManager::m_logger= 0;

}//close namespace
