#include "ConfigReader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

ConfigReader::ConfigReader(const std::string& filename) {
    loadConfig(filename);
}

void ConfigReader::loadConfig(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open config file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string key, value;
        if (std::getline(ss, key, '=') && std::getline(ss, value)) {
            // Remove any leading/trailing whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);

            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            configMap[key] = value;
        }
    }
}

int ConfigReader::getInt(const std::string& key) const {
    if (configMap.find(key) != configMap.end()) {
        return std::stoi(configMap.at(key));
    } else {
        throw std::runtime_error("Key not found: " + key);
    }
}

double ConfigReader::getDouble(const std::string& key) const {
    if (configMap.find(key) != configMap.end()) {
        return std::stod(configMap.at(key));
    } else {
        throw std::runtime_error("Key not found: " + key);
    }
}

std::string ConfigReader::getString(const std::string& key) const {
    if (configMap.find(key) != configMap.end()) {
        return configMap.at(key);
    } else {
        throw std::runtime_error("Key not found: " + key);
    }
}
