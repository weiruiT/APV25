#ifndef CRCONFIG_H
#define CRCONFIG_H 1

#include <string>

class CRConfig
{
public:
  bool Load(const std::string &filename = "config/configCR.json");
  bool GetBool(const std::string &key, bool defaultValue) const;
  double GetDouble(const std::string &key, double defaultValue) const;
  bool IsLoaded() const { return loaded; }

private:
  std::string FindValue(const std::string &key) const;

  std::string text;
  bool loaded = false;
};

#endif
