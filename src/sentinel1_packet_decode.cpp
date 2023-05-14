#include "sentinel1_packet_decode.h"
#include <intrin.h>
#include <bitset>

Sentinel1PacketDecode::Sentinel1PacketDecode() {

}

Sentinel1PacketDecode::~Sentinel1PacketDecode() {

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
        int i = 0, j = 0;

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
            sentinelOneParam.SASTestMode = std::bitset<8>(tmp8)[7];
            sentinelOneParam.CalType = static_cast<CalType>(tmp8 >> 4 & 0x7);
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

        if (sentinelOneParam.BAQMode == BAQMode::FDBAQMode0) {
            if (sentinelOneParam.SignalType == SignalType::Echo) {
                j += 1;
            }
            i += 1;
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

void Sentinel1PacketDecode::reconstruction(unsigned char* BRCn, unsigned char* THIDXn, ShCode* hcode, int NQ, float* result)
{
    int hcode_index = 0, h;
    int BRCindex = 0;
    int inc = 128;
    do
    {
        if ((hcode_index + 128) > NQ) {
            inc = (NQ - hcode_index);                      
        }

        for (h = 0; h < inc; h++) {
            switch (BRCn[BRCindex]) {
            case 0:
                if (THIDXn[BRCindex] <= 3) {
                    if (hcode[hcode_index].mcode < 3) {
                        result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
                    }
                    else if (hcode[hcode_index].mcode == 3) {
                        result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC0[THIDXn[BRCindex]];
                    }
                }
                else {
                    result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL0[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
                }
                break;
            case 1:
                if (THIDXn[BRCindex] <= 3)
                {
                    if (hcode[hcode_index].mcode < 4)
                        result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
                    else
                        result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC1[THIDXn[BRCindex]];
                }
                else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL1[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
                break;
            case 2:
                if (THIDXn[BRCindex] <= 5)
                {
                    if (hcode[hcode_index].mcode < 6)
                        result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
                    else
                        result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC2[THIDXn[BRCindex]];
                }
                else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL2[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
                break;
            case 3:
                if (THIDXn[BRCindex] <= 6)
                {
                    if (hcode[hcode_index].mcode < 9)
                        result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
                    else
                        result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC3[THIDXn[BRCindex]];
                }
                else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL3[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
                break;
            case 4:
                if (THIDXn[BRCindex] <= 8)
                {
                    if (hcode[hcode_index].mcode < 15)
                        result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
                    else
                        result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC4[THIDXn[BRCindex]];
                }
                else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL4[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
                break;
            default: 
                printf("UNHANDLED CASE\n"); 
                exit(-1);
            }
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
}

int Sentinel1PacketDecode::next_bit(unsigned char* p, int* cposition, int* bposition)
{
    int bit = ((p[*cposition] >> (*bposition)) & 1);
    (*bposition)--;
    if ((*bposition) < 0) { (*cposition)++; (*bposition) = 7; }
    return bit;
}

unsigned char Sentinel1PacketDecode::get_THIDX(unsigned char* p, int* cposition, int* bposition)
{
    int res = 0;
    int k;
    for (k = 0; k < 8; k++)
    {
        res = res << 1;
        res += next_bit(p, cposition, bposition);
    }
    return(res);
}

ShCode Sentinel1PacketDecode::BRC_4(unsigned char* p, int* cposition, int* bposition)
{
    int hcode, sign;
    ShCode sol;
    int b;
    sign = next_bit(p, cposition, bposition);
    if (sign == 0) sol.sign = 1; else sol.sign = -1;
    hcode = 0;
    do {
        b = next_bit(p, cposition, bposition);
        switch (hcode)
        {
        case 5:
            if (b == 0) // must be BEFORE hcode=5 at bottom
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0) { sol.mcode = 5; return(sol); }  // BRC4,1100
                else { sol.mcode = 6; return(sol); }  // BRC4,1101
            }
            else hcode = 6;
            break;
        case 6:
        case 7:
        case 8:
            if (b == 0) { sol.mcode = hcode + 1; return(sol); }
            else hcode++;
            break;
        case 9:
            if (b == 0)
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0) { sol.mcode = 10; return(sol); } // BRC4,11111100
                else { sol.mcode = 11; return(sol); } // BRC4,11111101
            }
            else hcode++;
            break;
        case 10:
            if (b == 0)
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0) { sol.mcode = 12; return(sol); } // BRC4,111111100
                else { sol.mcode = 13; return(sol); } // BRC4,111111101
            }
            else
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0) { sol.mcode = 14; return(sol); } // BRC4,111111100
                else { sol.mcode = 15; return(sol); } // BRC4,111111101
            }
            break;
        case 0:
            if (b == 0)   // first 0
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0) { sol.mcode = 0; return(sol); }  // BRC4,00 
                else
                {
                    b = next_bit(p, cposition, bposition);
                    if (b == 0) { sol.mcode = 1; return(sol); } // BRC4,010
                    else { sol.mcode = 2; return(sol); } // BRC4,011
                }
            }
            else
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0)  // BRC4,00 
                {
                    b = next_bit(p, cposition, bposition);
                    if (b == 0) { sol.mcode = 3; return(sol); } // BRC4,100
                    else { sol.mcode = 4; return(sol); } // BRC4,101
                }
                else hcode = 5;
            }
            break;
        }
    } while (hcode <= 15);
    sol.mcode = 999;
    exit(-1);    // should never get here
    return(sol);
}

ShCode Sentinel1PacketDecode::BRC(int BRCn, unsigned char* p, int* cposition, int* bposition)
{
    int hcode;
    int sign;
    int b;
    struct ShCode sol;
    switch (BRCn) {
    case 0: 
        BRCn = 3; 
        break; // number of steps to reach the leaves BRC0
    case 1: 
        BRCn = 4; 
        break; // number of steps to reach the leaves BRC1
    case 2: 
        BRCn = 6; 
        break; // number of steps to reach the leaves BRC2
    case 3: 
        BRCn = 9; 
        break; // number of steps to reach the leaves BRC3
    case 4: 
        return(BRC_4(p, cposition, bposition));
        printf("\nCheck if BRC4 output is correct\n"); 
        exit(0); 
        break;
    default: 
        printf("ERROR"); 
        exit(-1);
    }
    sign = next_bit(p, cposition, bposition);
    if (sign == 0) sol.sign = 1; else sol.sign = -1;
    hcode = 0;
    do {
        b = next_bit(p, cposition, bposition);
        if (b == 0)                                   // 0: end of tree
            if ((BRCn == 9) && (hcode == 0))             // first 0 (hcode==0) of BRC3
            {
                b = next_bit(p, cposition, bposition);
                if (b == 0)
                {
                    sol.mcode = hcode; return(sol);
                }     // we reached 00 of BRC3
                else
                {
                    sol.mcode = hcode + 1; return(sol);
                }   // we reached 01 of BRC3
            }
            else
            {
                sol.mcode = hcode; return(sol);
            }         // we reached 0 -> return
        else
        {
            hcode++;                                // 1: continue
            if ((BRCn == 9) && (hcode == 1)) hcode++;
            if (hcode == BRCn)                        // unless last 1 was reached 
            {
                sol.mcode = hcode; return(sol);
            }       // end of tree reached
        }
    } while (hcode < BRCn);
    exit(-1);                                     // ERROR in decoding Huffman
    sol.mcode = 99;
    return(sol);                                  // should never be reached
}

int Sentinel1PacketDecode::packet_decode(unsigned char* p, int NQ, float* IE, float* IO, float* QE, float* QO, char* brc, int* brcpos) {

    std::vector<ShCode> hcodeIE(52378);
    std::vector<ShCode> hcodeIO(52378);
    std::vector<ShCode> hcodeQE(52378);
    std::vector<ShCode> hcodeQO(52378);
    std::vector<uint8_t> BRCn(410);
    std::vector<uint8_t> THIDXn(410);

    int BRCindex;
    int h;
    int hcode_index;
    int cposition = 0;
    int bposition = 7;
    int inc = 128;  
    BRCindex = 0;
    hcode_index = 0;
    do {
        BRCn[BRCindex] = next_bit(p, &cposition, &bposition) * 4;  // MSb=4
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 2; // then 2
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 1; // then 1 ...
        brc[*brcpos] = BRCn[BRCindex];
        (*brcpos)++;
        if ((hcode_index + 128) > NQ) {
            inc = (NQ - hcode_index);      
        }
        for (h = 0; h < inc; h++) {
            hcodeIE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    
    inc = 128;
    if (bposition != 7) { 
        bposition = 7; 
        cposition++; 
    } 
    
    if ((cposition & 1) != 0) { 
        cposition++; 
    }        
    BRCindex = 0;
    hcode_index = 0;
    do
    {
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeIO[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    inc = 128;
    if (bposition != 7) { 
        bposition = 7; 
        cposition++; 
    } 
    if ((cposition & 1) != 0) { 
        cposition++; 
    } 
    BRCindex = 0;
    hcode_index = 0;
    do
    {
        THIDXn[BRCindex] = get_THIDX(p, &cposition, &bposition);
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeQE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    inc = 128;
    if (bposition != 7) { 
        bposition = 7; 
        cposition++; 
    } // start at new position
    if ((cposition & 1) != 0) { 
        cposition++; 
    }     // odd address => +1 
    BRCindex = 0;
    hcode_index = 0;
    do
    {
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) {
            hcodeQO[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    if (bposition != 7) { 
        bposition = 7; 
        cposition++; 
    } // start at new position
    if ((cposition & 1) != 0) { 
        cposition++; 
    }  // odd address => +1 
    reconstruction(BRCn.data(), THIDXn.data(), hcodeIE.data(), NQ, IE);
    reconstruction(BRCn.data(), THIDXn.data(), hcodeIO.data(), NQ, IO);
    reconstruction(BRCn.data(), THIDXn.data(), hcodeQE.data(), NQ, QE);
    reconstruction(BRCn.data(), THIDXn.data(), hcodeQO.data(), NQ, QO);
    return(cposition);
}

