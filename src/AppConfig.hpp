#ifndef APPCONFIG_H
#define APPCONFIG_H

#include <vector>
#include <string>
#include <fstream>

// PARAMETERS   ------------------------------------------------------
struct AppConfig{
    std::string dir;
    std::vector<std::string> ifiles;
    std::string outf;
    std::string logf;
    std::string app;
    bool help;
    std::ofstream ofs;
    std::ofstream lfs;
    int kv;
    int method;
    int extend;
    bool only_lcp;

    void printHelp(std::ostream& ots);
    bool validate(std::ostream& ots);
    void write(std::ostream& ots);
    AppConfig(const std::vector<std::string>& files,
              const std::string& of,
              const std::string& lf,
              int kval, int mt, int ext = 0);
    AppConfig(int argc, char** argv);
    ~AppConfig();
};


#endif /* APPCONFIG_H */
