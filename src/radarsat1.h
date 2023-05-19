#pragma once
#include <fstream>
#include <filesystem>

class CEOSFormat
{
public:
	CEOSFormat(std::filesystem::path &pathVolumeFile, 
        std::filesystem::path &pathTrailerFile,
        std::filesystem::path &pathLeaderFile,
        std::filesystem::path &pathDataSetFile,
        std::filesystem::path &pathNullVolumeFile);
    
	~CEOSFormat();
    
    int parseVolumeDirectoryFile();
    int parseSARLeaderFile();
    int parseSARDataFile();
    int parseSARTrailerFile();
    int parseNullVolumeDirectoryFile();

private:
    std::filesystem::path _pathVolumeFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.vol";
    std::filesystem::path _pathTrailerFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.tlr";
    std::filesystem::path _pathLeaderFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.ldr";
    std::filesystem::path _pathDataSetFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.raw";
    std::filesystem::path _pathNullVolumeFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.nul";
};

