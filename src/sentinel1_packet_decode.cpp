#include "sentinel1_packet_decode.h"
#include <intrin.h>
#include <bitset>



Sentinel1PacketDecode::Sentinel1PacketDecode()
{

}

Sentinel1PacketDecode::~Sentinel1PacketDecode()
{

}

float Sentinel1PacketDecode::RangeDecimationToSampleRate(RangeDecimation& rangeDecimation) {
    float fRef = 37.53472224;
    float sampleRate = 0.0;
    switch (rangeDecimation) {
    case RangeDecimation::Filter0:
        sampleRate = 3 * fRef;
        break;
    case RangeDecimation::Filter1:
        sampleRate = (8 / 3) * fRef;
        break;
    case RangeDecimation::Filter3:
        sampleRate = (20 / 9) * fRef;
        break;
    case RangeDecimation::Filter4:
        sampleRate = (16 / 9) * fRef;
        break;
    case RangeDecimation::Filter5:
        sampleRate = (3 / 2) * fRef;
        break;
    case RangeDecimation::Filter6:
        sampleRate = (4 / 3) * fRef;
        break;
    case RangeDecimation::Filter7:
        sampleRate = (2 / 3) * fRef;
        break;
    case RangeDecimation::Filter8:
        sampleRate = (12 / 7) * fRef;
        break;
    case RangeDecimation::Filter9:
        sampleRate = (5 / 4) * fRef;
        break;
    case RangeDecimation::Filter10:
        sampleRate = (6 / 13) * fRef;
        break;
    case RangeDecimation::Filter11:
        sampleRate = (16 / 11) * fRef;
        break;
    default:
        break;
    }
    return sampleRate;
}

float Sentinel1PacketDecode::calcRxGain(uint8_t rawRxGain) {
    return -(static_cast<float>(rawRxGain) * 0.5);
}

float Sentinel1PacketDecode::calcTxPulseRampRate(uint16_t rawTxPulseRampRate) {
    rawTxPulseRampRate = _byteswap_ushort(rawTxPulseRampRate);
    float fRef = 37.53472224;
    uint8_t sign = std::bitset<16>(rawTxPulseRampRate)[15];
    uint16_t value = rawTxPulseRampRate & 0x7FFF;
    return std::pow(-1, sign) * static_cast<float>(value) * std::pow(fRef, 2) * std::pow(2, -21);
}

float Sentinel1PacketDecode::calcTxPulseStartFreq(uint16_t rawTxPulseStartFreq, float TxPulseRampRate) {
    rawTxPulseStartFreq = _byteswap_ushort(rawTxPulseStartFreq);
    float fRef = 37.53472224;
    uint8_t sign = std::bitset<16>(rawTxPulseStartFreq)[15];
    uint16_t value = rawTxPulseStartFreq & 0x7FFF;
    return TxPulseRampRate / (4 * fRef) + std::pow(-1, sign) * static_cast<float>(value) * fRef * std::pow(2, -14);
}

uint32_t Sentinel1PacketDecode::calcTxPulseLength(uint32_t rawTxPulseStartFreq) {
    float fRef = 37.53472224;
    return static_cast<float>(rawTxPulseStartFreq) / fRef;
}

float Sentinel1PacketDecode::calcPulseRepetitionInterval(uint32_t rawPulseRepetitionInterval) {
    float fRef = 37.53472224;
    return static_cast<float>(rawPulseRepetitionInterval) / fRef;
}

float Sentinel1PacketDecode::calcSamplingWindowStartTime(uint32_t rawSamplingWindowStartTime) {
    float fRef = 37.53472224;
    return static_cast<float>(rawSamplingWindowStartTime) / fRef;
}

float Sentinel1PacketDecode::calcSamplingWindowLength(uint32_t rawSamplingWindowLength) {
    float fRef = 37.53472224;
    return static_cast<float>(rawSamplingWindowLength) / fRef;
}

float Sentinel1PacketDecode::calcFineTime(uint16_t rawFineTime) {
    rawFineTime = _byteswap_ushort(rawFineTime);
    return (static_cast<float>(rawFineTime) + 0.5) * std::pow(2, -16.0);
}

uint32_t Sentinel1PacketDecode::calcInstrumentConfigID(uint32_t rawInstrumentConfigID) {
    return _byteswap_ulong(rawInstrumentConfigID);
}

uint32_t Sentinel1PacketDecode::read24Bit(std::ifstream& f) {
    uint8_t tmp8 = 0;
    uint16_t tmp16 = 0;
    f.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
    f.read(reinterpret_cast<char*>(&tmp16), sizeof(tmp16));
    tmp16 = _byteswap_ushort(tmp16);
    return (tmp8 << 8) | tmp16;
}

int Sentinel1PacketDecode::ReadSARParam(std::filesystem::path pathToRawData) {
    uint8_t tmp8 = 0;
    uint16_t tmp16 = 0;
    uint32_t tmp32 = 0;

    std::ifstream rawData(pathToRawData, std::ios::binary);
    if (!rawData.is_open()) {
        return -1;
    }

    Sentinel1RawPacket sentinelOneParam;

    while (true) {
        // Octet 0,1
        rawData.read(reinterpret_cast<char*>(&tmp16), 2);
        if (rawData.gcount() == 0) {
            break;
        }
        tmp16 = _byteswap_ushort(tmp16);
        sentinelOneParam.PacketVersionNumber = (tmp16 >> 12);
        sentinelOneParam.PacketType = std::bitset<16>(tmp16)[12];
        sentinelOneParam.SecondaryHeaderFlag = std::bitset<16>(tmp16)[11];
        sentinelOneParam.ProcessID = (tmp16 >> 4) & 0x7F;
        sentinelOneParam.PacketCategory = (tmp16) & 0xF;

        // Octet 2,3
        rawData.read(reinterpret_cast<char*>(&tmp16), 2);
        tmp16 = _byteswap_ushort(tmp16);
        sentinelOneParam.SequenceFlags = (tmp16 >> 14);
        sentinelOneParam.PacketSequenceCount = (tmp16 & 0x3FF);

        // Octet 4,5
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.PacketDataLength), sizeof(sentinelOneParam.PacketDataLength));
        sentinelOneParam.PacketDataLength = _byteswap_ushort(sentinelOneParam.PacketDataLength) + 1;
        if (((sentinelOneParam.PacketDataLength + 6) % 4) != 0) {
            std::cerr << "ERROR: Length not multiple of 4." << std::endl;
        }

        // Octet 6,7,8,9
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.Time), sizeof(sentinelOneParam.Time));
        sentinelOneParam.Time = _byteswap_ulong(sentinelOneParam.Time);

        // Octet 10,11
        rawData.read(reinterpret_cast<char*>(&tmp16), sizeof(uint16_t));
        sentinelOneParam.FineTime = calcFineTime(tmp16);

        // Octet 12,13,14,15
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SyncMarker), sizeof(sentinelOneParam.SyncMarker));
        sentinelOneParam.SyncMarker = _byteswap_ulong(sentinelOneParam.SyncMarker);

        if (sentinelOneParam.SyncMarker != 0x352EF853) {
            std::cerr << "ERROR: Sync marker != 352EF853" << std::endl;
        }

        // Octet 16,17,18,19
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.DataTakeID), sizeof(sentinelOneParam.DataTakeID));
        sentinelOneParam.DataTakeID = _byteswap_ulong(sentinelOneParam.DataTakeID);

        // Octet 20
        rawData.read(reinterpret_cast<char*>(&tmp8), 1);
        sentinelOneParam.ECCNumber = static_cast<ECCNumber>(tmp8);

        // Octet 21
        rawData.read(reinterpret_cast<char*>(&tmp8), 1);
        sentinelOneParam.TestMode = static_cast<TestMode>(0/*tmp8 & 0x0E*/);
        sentinelOneParam.RxChannelId = std::bitset<8>(tmp8)[0];

        // Octet 22,23,24,25
        rawData.read(reinterpret_cast<char*>(&tmp32), sizeof(uint32_t));
        sentinelOneParam.InstrumentConfigID = calcInstrumentConfigID(tmp32);

        // Octet 26
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.WordIndex), sizeof(sentinelOneParam.WordIndex));

        // Octet 27, 28
        rawData.read(reinterpret_cast<char*>(&tmp16), sizeof(uint16_t));
        sentinelOneParam.WordVal = _byteswap_ushort(tmp16);

        // Octet 29, 30, 31, 32
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SpacePacketCount), sizeof(sentinelOneParam.SpacePacketCount));
        sentinelOneParam.SpacePacketCount = _byteswap_ulong(sentinelOneParam.SpacePacketCount);

        // Octet 33, 34, 35, 36
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.PRICount), sizeof(sentinelOneParam.PRICount));
        sentinelOneParam.PRICount = _byteswap_ulong(sentinelOneParam.PRICount);

        // Octet 37
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.BAQMode = static_cast<BAQMode>(tmp8 & 0x1f);
        sentinelOneParam.ErrorFlag = std::bitset<8>(tmp8)[7];

        // Octet 38
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.BAQBlockLength = 8 * (tmp8 + 1);

        // Octet 39
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));

        // Octet 40
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.RangeDecimation = static_cast<RangeDecimation>(tmp8);

        // Octet 41
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.RxGain = calcRxGain(tmp8);

        // Octet 42, 43
        rawData.read(reinterpret_cast<char*>(&tmp16), sizeof(tmp16));
        sentinelOneParam.TxRampRate = calcTxPulseRampRate(tmp16);

        // Octet 44, 45
        rawData.read(reinterpret_cast<char*>(&tmp16), sizeof(tmp16));
        sentinelOneParam.TxPulseStartFreq = calcTxPulseStartFreq(tmp16, sentinelOneParam.TxRampRate);

        // Octet 46, 47, 48
        sentinelOneParam.TxPulseLength = calcTxPulseLength(read24Bit(rawData));

        // Octet 49
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.Rank = tmp8 & 0x1F;

        // Octet 50, 51, 52
        sentinelOneParam.PulseRepetitionInterval = calcPulseRepetitionInterval(read24Bit(rawData));
        // Octet 53, 54, 55
        sentinelOneParam.SamplingWindowStartTime = calcSamplingWindowStartTime(read24Bit(rawData));
        // Octet 56, 57, 58
        sentinelOneParam.SamplingWindowLength = calcSamplingWindowLength(read24Bit(rawData));

        // Octet 59
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.SSBFlag = std::bitset<8>(tmp8)[7];
        sentinelOneParam.Polarisation = static_cast<Polarisation>((tmp8 & 0x70) >> 4);
        sentinelOneParam.TemperatureCompensation = static_cast<TemperatureCompensation>(tmp8 & 0x0C);
        if (sentinelOneParam.SSBFlag == 0) {
            // Octet 60
            rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
            sentinelOneParam.ElevationBeamAddress = (tmp8 & 0xF0) >> 4;
            tmp16 = tmp8 & 0xC0;
            // Octet 61
            rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
            sentinelOneParam.AzimuthBeamAddress = tmp16 & (tmp8 << 2);
        }
        else {
            // Octet 60
            rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
            sentinelOneParam.SASTestMode = std::bitset<8>(tmp8)[0];
            sentinelOneParam.CalType = static_cast<CalType>(tmp8 & 0x0E);
            tmp16 = tmp8 & 0xC0;
            // Octet 61
            rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
            sentinelOneParam.CalibrationBeamAddress = tmp16 & (tmp8 << 2);
        }
        // Octet 62
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.CalMode = (tmp8 & 0xC0) >> 6;
        sentinelOneParam.TxPulseNumber = tmp8 & 0x1F;
        // Octet 63
        rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        sentinelOneParam.SignalType = static_cast<SignalType>((tmp8 & 0xF0) >> 4);
        sentinelOneParam.SwapFlag = std::bitset<8>(tmp8)[0];
        // Octet 64
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SwathNumber), sizeof(sentinelOneParam.SwathNumber));

        // Octet 65, 66
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.NumberOfQuads), sizeof(sentinelOneParam.NumberOfQuads));
        sentinelOneParam.NumberOfQuads = _byteswap_ushort(sentinelOneParam.NumberOfQuads);

        if (static_cast<int>(sentinelOneParam.SignalType) > 7) {
            if (static_cast<int>(sentinelOneParam.BAQMode) != 0) {
                std::cerr << "Calibration data(SIGTYPcode > 7, all CALTYPcode) are only with BAQMODcode = 0";
            }
        }

        if (static_cast<int>(sentinelOneParam.SignalType) == 1) {
            if (static_cast<int>(sentinelOneParam.BAQMode) != 0 && static_cast<int>(sentinelOneParam.BAQMode) != 3
                && static_cast<int>(sentinelOneParam.BAQMode) != 4 && static_cast<int>(sentinelOneParam.BAQMode) != 5) {
                std::cerr << "Noise data (SIGTYPcode =1) are only with BAQMODcode =0 or 3 or 4 or 5 ";
            }
        }

        for (size_t i = 0; i < sentinelOneParam.PacketDataLength - 62 + 1; i++) {
            rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
        }

        if (sentinelOneParam.BAQMode == BAQMode::BypassMode && (static_cast<int>(sentinelOneParam.CalType) > 7)) {

        }

        if ((sentinelOneParam.BAQMode == BAQMode::FDBAQMode0) && (sentinelOneParam.CalType == CalType::TxCal)) {

        }

        header.emplace_back(sentinelOneParam);
    };

    PositionVelocityTime pvt;
    Attitude att;

    uint16_t tmp[4] = { 0, 0, 0, 0 };

    for (size_t i = 0; i < header.size(); i++) {
        if (header.at(i).WordIndex == 1) {
            tmp[3] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 2) {
            tmp[2] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 3) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 4) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&pvt.XAxisPositionECEF, tmp, sizeof(double));
        }
        else if (header.at(i).WordIndex == 5) {
            tmp[3] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 6) {
            tmp[2] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 7) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 8) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&pvt.YAxisPositionECEF, tmp, sizeof(double));
        }
        if (header.at(i).WordIndex == 9) {
            tmp[3] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 10) {
            tmp[2] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 11) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 12) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&pvt.ZAxisPositionECEF, tmp, sizeof(double));
        }
        else if (header.at(i).WordIndex == 13) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 14) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&pvt.XvelocityECEF, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 15) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 16) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&pvt.YvelocityECEF, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 17) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 18) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&pvt.ZvelocityECEF, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 19) {
            pvt.DataStampOne = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 20) {
            pvt.DataStampTwo = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 21) {
            pvt.DataStampThree = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 22) {
            pvt.DataStampFour = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 23) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 24) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.Q0AttitudeQuaternion, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 25) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 26) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.Q1AttitudeQuaternion, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 27) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 28) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.Q2AttitudeQuaternion, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 29) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 30) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.Q3AttitudeQuaternion, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 31) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 32) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.XangularRate, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 33) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 34) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.YangularRate, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 35) {
            tmp[1] = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 36) {
            tmp[0] = header.at(i).WordVal;
            std::memcpy(&att.ZangularRate, tmp, sizeof(float));
        }
        else if (header.at(i).WordIndex == 37) {
            att.DataStampOne = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 38) {
            att.DataStampTwo = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 39) {
            att.DataStampThree = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 40) {
            att.DataStampFour = header.at(i).WordVal;
        }
        else if (header.at(i).WordIndex == 41) {
            //header.at(i).WordVal;
            positionVelocityTime.emplace_back(pvt);
            attitude.emplace_back(att);
        }
    }
}