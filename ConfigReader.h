#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <string>
#include <unordered_map>

class ConfigReader {
public:
    ConfigReader(const std::string& filename);

    int getInt(const std::string& key) const;
    double getDouble(const std::string& key) const;
    std::string getString(const std::string& key) const;

private:
    std::unordered_map<std::string, std::string> configMap;

    void loadConfig(const std::string& filename);
};

#endif // CONFIG_READER_H

