#pragma once
#include <iostream>
#include <string>
#include <variant>
#include <vector>

class CmdLine {
public:
  typedef std::variant<int32_t*, float*, bool*, std::string*> Value;

  explicit CmdLine(const std::string& _description) : description(_description) {}

  void addArgument(const std::vector<std::string>& flags, const std::string& help, Value value);
  void addCategory(const std::string& category);
  void printHelp(std::ostream& os = std::cout) const;
  bool parse(int argc, const char* argv[]) const;

private:
  struct Argument {
    Argument(const std::vector<std::string>& _flags, const std::string& _help, Value _value) :
      flags(_flags), help(_help), value(_value) {}
    std::vector<std::string> flags;
    std::string help;
    Value value;
  };

  std::string description;
  std::vector<std::pair<int,std::string>> categories;
  std::vector<Argument> arguments;
};
