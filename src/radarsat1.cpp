#include "radarsat1.h"

CEOSFormat::CEOSFormat(std::filesystem::path& pathVolumeFile,
    std::filesystem::path& pathTrailerFile,
    std::filesystem::path& pathLeaderFile,
    std::filesystem::path& pathDataSetFile,
    std::filesystem::path& pathNullVolumeFile)
{
    _pathVolumeFile = pathVolumeFile;
    _pathTrailerFile = pathTrailerFile;
    _pathLeaderFile = pathLeaderFile;
    _pathDataSetFile = pathDataSetFile;
    _pathNullVolumeFile = pathNullVolumeFile;
}

int CEOSFormat::parseVolumeDirectoryFile() {
    // Volume descriptor record
    std::ifstream vol(_pathVolumeFile);

    uint8_t tmp8 = 0;
    uint8_t tmp16 = 0;
    uint32_t tmp32 = 0;

    // Bytes 1 -4
    vol.read(reinterpret_cast<char*>(&tmp32), sizeof(tmp32));
    tmp32 = _byteswap_ulong(tmp32);
    // Byte 5
    vol.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
    // Byte 6
    vol.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
    // Byte 7
    vol.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
    // Byte 8
    vol.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
    // Byte 9-12
    vol.read(reinterpret_cast<char*>(&tmp32), sizeof(tmp32));
    tmp32 = _byteswap_ulong(tmp32);
    // Byte 13-14
    vol.read(reinterpret_cast<char*>(&tmp16), sizeof(tmp16));
    tmp16 = _byteswap_ulong(tmp16);

    // File pointer records 
    // Text record
    return 0;
}

int CEOSFormat::parseSARLeaderFile() {
    return 0;
}

int CEOSFormat::parseSARDataFile() {
    return 0;
}

int CEOSFormat::parseSARTrailerFile() {
    return 0;
}

int CEOSFormat::parseNullVolumeDirectoryFile() {
    return 0;
}

CEOSFormat::~CEOSFormat() {

}