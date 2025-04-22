#include "cmdLine.h"
#include <iomanip>
#include <sstream>

void CmdLine::addArgument(const std::vector<std::string>& flags, const std::string& help, Value value) {
  arguments.emplace_back(flags, help, value);
}
void CmdLine::addCategory(const std::string& category) {
  categories.emplace_back((int)arguments.size(), category);
}

void CmdLine::printHelp(std::ostream& os) const {
  os << description << std::endl;

  int maxFlagLength = 0;
  for (const Argument& argument : arguments) {
    int flagLength = 0;
    for (const std::string& flag : argument.flags) {
      flagLength += (int)(flag.size() + 2);
    }
    maxFlagLength = std::max(maxFlagLength, flagLength);
  }

  size_t catIx = 0;
  for (size_t i = 0; i < arguments.size(); ++i) {
    if (catIx < categories.size() && categories[catIx].first == i) {
      std::cout << std::endl << categories[catIx].second << std::endl;
      catIx += 1;
    }

    std::string flags;
    const Argument& argument = arguments[i];
    for (const std::string& flag : argument.flags) {
      flags += flag + ", ";
    }

    std::stringstream sstr;
    sstr << std::left << std::setw(maxFlagLength) << flags.substr(0, flags.size() - 2);

    size_t spacePos  = 0;
    size_t lineWidth = 0;
    while (spacePos != std::string::npos) {
      size_t nextspacePos = argument.help.find_first_of(' ', spacePos + 1);
      sstr << argument.help.substr(spacePos, nextspacePos - spacePos);
      lineWidth += nextspacePos - spacePos;
      spacePos = nextspacePos;
      if (lineWidth > 60) {
        os << sstr.str() << std::endl;
        sstr = std::stringstream();
        sstr << std::left << std::setw(maxFlagLength - 1) << " ";
        lineWidth = 0;
      }
    }
  }
}

bool CmdLine::parse(int argc, const char* argv[]) const {
  int i = 1;
  while (i < argc) {
    std::string flag(argv[i]);
    std::string value;
    bool valueIsSeparate = false;

    size_t equalPos = flag.find('=');
    if (equalPos != std::string::npos) {
      value = flag.substr(equalPos + 1);
      flag = flag.substr(0, equalPos);
    } else if (i + 1 < argc && argv[i + 1][0] != '-') {
      value = argv[i + 1];
      valueIsSeparate = true;
    }

    bool foundArgument = false;
    for (const Argument& argument : arguments) {
      if (std::find(argument.flags.begin(), argument.flags.end(), flag) != std::end(argument.flags)) {
        foundArgument = true;
        if (std::holds_alternative<bool*>(argument.value)) {
          if (!value.empty() && value != "true" && value != "false") {
            valueIsSeparate = false;
          }
          *std::get<bool*>(argument.value) = (value != "false");
        } else if (value.empty()) {
          std::cout << "ERROR: Missing value for argument \"" + flag + "\"." << std::endl;
          return false;
        } else if (std::holds_alternative<std::string*>(argument.value)) {
          *std::get<std::string*>(argument.value) = value;
        } else {
          std::visit([&value](auto&& arg) {
            std::stringstream sstr(value);
            sstr >> *arg;
          }, argument.value);
        }
        break;
      }
    }

    if (!foundArgument) {
      std::cout << "ERROR: Unknown command line argument \"" << flag << "\"." << std::endl;
      return false;
    }

    ++i;
    if (foundArgument && valueIsSeparate) {
      ++i;
    }
  }
  return true;
}
