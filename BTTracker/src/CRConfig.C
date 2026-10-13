#include "CRConfig.h"

#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace
{
std::string Trim(const std::string &value)
{
  const std::string whitespace = " \t\r\n";
  const size_t first = value.find_first_not_of(whitespace);
  if (first == std::string::npos)
    return "";

  const size_t last = value.find_last_not_of(whitespace);
  return value.substr(first, last - first + 1);
}

std::string Lower(std::string value)
{
  for (char &character : value)
    character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
  return value;
}
}

bool CRConfig::Load(const std::string &filename)
{
  std::ifstream input(filename);
  if (!input.is_open())
  {
    std::cerr << "WARNING: cannot open " << filename
              << "; all configured stages will use their defaults." << std::endl;
    text.clear();
    loaded = false;
    return false;
  }

  text.assign(std::istreambuf_iterator<char>(input),
              std::istreambuf_iterator<char>());
  loaded = !text.empty();
  if (!loaded)
    std::cerr << "WARNING: configuration file is empty: " << filename << std::endl;
  return loaded;
}

std::string CRConfig::FindValue(const std::string &key) const
{
  const std::string quotedKey = "\"" + key + "\"";
  size_t position = text.find(quotedKey);
  while (position != std::string::npos)
  {
    const size_t colon = text.find(':', position + quotedKey.size());
    const size_t nextKey = text.find('"', position + quotedKey.size());
    if (colon != std::string::npos &&
        (nextKey == std::string::npos || colon < nextKey))
    {
      size_t valueBegin = colon + 1;
      while (valueBegin < text.size() &&
             std::isspace(static_cast<unsigned char>(text[valueBegin])))
        ++valueBegin;

      if (valueBegin >= text.size())
        return "";

      if (text[valueBegin] == '"')
      {
        const size_t valueEnd = text.find('"', valueBegin + 1);
        if (valueEnd == std::string::npos)
          return "";
        return text.substr(valueBegin + 1, valueEnd - valueBegin - 1);
      }

      size_t valueEnd = text.find_first_of(",}\n\r", valueBegin);
      if (valueEnd == std::string::npos)
        valueEnd = text.size();
      return Trim(text.substr(valueBegin, valueEnd - valueBegin));
    }

    position = text.find(quotedKey, position + quotedKey.size());
  }
  return "";
}

bool CRConfig::GetBool(const std::string &key, bool defaultValue) const
{
  const std::string value = Lower(Trim(FindValue(key)));
  if (value == "true" || value == "1" || value == "yes" || value == "on")
    return true;
  if (value == "false" || value == "0" || value == "no" || value == "off")
    return false;
  return defaultValue;
}

double CRConfig::GetDouble(const std::string &key, double defaultValue) const
{
  const std::string value = Trim(FindValue(key));
  if (value.empty())
    return defaultValue;
  char *end = nullptr;
  errno = 0;
  const double parsed = std::strtod(value.c_str(), &end);
  if (errno != 0 || end == value.c_str() || *end != '\0' ||
      !std::isfinite(parsed))
    return defaultValue;
  return parsed;
}
