#include <vector>
#include <iostream>
#include "ipp.h"
#include <complex>
#include <filesystem>
#include <fstream>
#include <intrin.h>

template<typename T>
struct Point {
    T x = static_cast <typeid(T).name()>(0);
    T y = static_cast <typeid(T).name()>(0);
    T z = static_cast <typeid(T).name()>(0);
};

struct ComplexMatrix {
    std::vector<std::vector<double>> re;
    std::vector<std::vector<double>> im;
    uint64_t rows = 0;
    uint64_t cols = 0;
};


struct RangeReferenceFunction {
    std::vector<std::complex<float>> refFunc;
};

// Calculate the mean of the raw data.
template<typename T>
double MeanOfRawData(ComplexMatrix &rawData) {
    double sum = 0.0;

    for (size_t i = 0; i < re.size(); i++) {
        sum += static_cast<double>(accumulate(rawData.re[i].begin(), rawData.re[i].end(), 0));
    }

    sum = sum / (static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols));
    return sum;
}

// Calculate the standard deviations of the raw data.
template<typename T>
T StandardDeviationsOfRawData(std::vector<std::vector<T>> rawData) {
    
    std::sqrt();
    return static_cast <typeid(T).name()>(0);
}

// Calculate the IQ gain imbalance
template<typename T>
T IQGainImbalance(std::vector<std::vector<T>> rawData) {
    return static_cast <typeid(T).name()>(0);
}

// IQ quadrature departure
template<typename T>
T IQQuadratureDeparture(std::vector<std::vector<T>> rawData) {
    return static_cast <typeid(T).name()>(0);
}

// Set the statistics significance flags

constexpr double speedOfLight = 299792458.0;

template<typename T>
T MigrationFactor(T carrierFrequency, T freqAzimuth, T effectiveRadarVelocity) {
    return std::sqrt(1 - (std::pow(speedOfLight / carrierFrequency,2) * std::pow(freqAzimuth, 2) / 4 * std::pow(effectiveRadarVelocity,2)));
}

template<typename T>
T AzimuthFMRate(T carrierFrequency, T frequencyDopplerCentroid , T effectiveRadarVelocity, T slantRange) {
    return 2 * std::pow(effectiveRadarVelocity, 2) * MigrationFactor<T>(carrierFrequency, frequencyDopplerCentroid, effectiveRadarVelocity) / ((speedOfLight / carrierFrequency) * slantRange);
}

template<typename T>
T EffectiveRadarVelocity() {
    return static_cast <typeid(T).name()>(0);
}

template<typename T>
int rangeCompression(ComplexMatrix &rawData, bool calcIfft = true) {
    Ipp32fc* data_ref = nullptr;
    Ipp32fc* correl = nullptr;
    
    int sizeDFTSpec = 0;
    int sizeDFTInitBuf = 0;
    int sizeDFTWorkBuf = 0;

    Ipp8u* pDFTSpec = nullptr;
    Ipp8u* pDFTInitBuf = nullptr;
    Ipp8u* pDFTWorkBuf = nullptr;
    
    IppsFFTSpec_C_32fc* pSpec = nullptr;
    
    int fft_len = 0; 
     
    data_ref = ippsMalloc_32fc(fft_len);
    if (data_ref == nullptr)
        return -1;

    correl = ippsMalloc_32fc(fft_len);
    if (correl == nullptr)
        return -1;

    ippsZero_32fc(correl, fft_len);
    ippsZero_32fc(data_ref, fft_len);

    pDFTSpec = ippsMalloc_8u(sizeDFTSpec);
    if (pDFTSpec == nullptr)
        return -1;
    pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
    if (pDFTInitBuf == nullptr)
        return -1;
    pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);
    if (pDFTWorkBuf == nullptr)
        return -1;

    //ippsFFTGetSize_C_32fc(_correlParam.size, IPP_NODIV_BY_ANY, ippAlgHintAccurate, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);
    
    //ippsFFTInit_C_32fc(&pSpec, _correlParam.size, IPP_NODIV_BY_ANY, ippAlgHintAccurate, pDFTSpec, pDFTInitBuf);

    ippsFFTFwd_CToC_32fc(data_ref, data_ref, pSpec, pDFTWorkBuf);

    ippsConj_32fc(data_ref, data_ref, fft_len);

    for (size_t i = 0; i < rawData.rows; i++)
    {
        //ippsFFTFwd_CToC_32fc(, correl, pSpec, pDFTWorkBuf);

        ippsMul_32fc(correl, data_ref, correl, fft_len);
        if (calcIfft) {
            ippsFFTInv_CToC_32fc(correl, correl, pSpec, pDFTWorkBuf);
        }
    }

    ippsFree(pDFTSpec);
    ippsFree(pDFTInitBuf);
    ippsFree(pDFTWorkBuf);
    ippsFree(data_ref);
    ippsFree(correl);

    return 0;
}

template<typename T>
void  NominalImageReplicaGeneration(T  TXPL, T TXPSF, T TXPRR) {
    std::vector<std::complex<float>> refFunc;
    
    for (size_t i = 0; i < length; i++) {
        refFunc.emplace_back(std::exp(std::complex<float>(2 * M_PI * (((TXPSF - TXPRR * (-TXPL / 2))* i +(TXPRR/2) * std::pow(i, 2))))));
    }
}

template<typename T>
T getRangeReferenceFunction(T  TXPL, T TXPSF, T TXPRR) {
    return static_cast <typeid(T).name()>(0);
}

struct Sentinel1RawPacket {
    uint8_t PacketVersionNumber;
    uint8_t PacketType;
    uint8_t SecondaryHeaderFlag;
    uint8_t ProcessID;
    uint8_t PacketCategory;
    uint8_t SequenceFlags;
    uint16_t PacketSequenceCount;
    uint16_t PacketDataLength;
    uint32_t Time;
    uint32_t FineTime;
    uint32_t SyncMarker;
    uint32_t DataTakeID;
    uint8_t ECCNumber;
    uint8_t TestMode;
    bool RxChannelId;
    uint32_t InstrumentConfigID;
    uint8_t WordIndex;
    uint16_t WordVal;
    uint32_t SpacePacketCount;
    uint32_t PRICount;
    uint8_t BAQMode;
    bool ErrorFlag;
    uint8_t BAQBlockLength;
    uint8_t RangeDecimation;
    uint32_t RxGain;
    uint16_t TxRampRate;
    uint16_t TxPulseStartFreq;
    uint32_t TxPulseLength;
    uint8_t Rank;
    uint32_t PulseRepetitionInterval;
    uint32_t SamplingWindowStartTime;
    uint32_t SamplingWindowLength;
    uint8_t Polarisation;
    uint8_t TemperatureCompensation;
    uint8_t ElevationBeamAddress;
    uint16_t AzimuthBeamAddress;
    bool SASTestMode;
    uint8_t CalType;
    uint16_t CalibrationBeamAddress;
    uint8_t CalMode;
    uint8_t TxPulseNumber;
    uint8_t SignalType;
    bool SwapFlag;
    uint8_t SwathNumber;
};

enum TestMode {
    Default = 0,
    GroundTestingOnlyModeOne = 0x4,
    GroundTestingOnlyModeTwo = 0x5,
    Oper = 0x6,
    Bypass = 0x7
};

enum BAQMode {
    BypassMode = 0,
    BAQ3BitMode = 3,
    BAQ4BitMode = 4,
    BAQ5BitMode = 5,
    FDBAQMode0 = 12,
    FDBAQMode1 = 13,
    FDBAQMode2 = 14
};

enum RangeDecimation {
    Filter0,
    Filter1,
    Filter2,
    Filter3,
    Filter4,
    Filter5,
    Filter6,
    Filter7,
    Filter8,
    Filter9,
    Filter10,
    Filter11
};

enum TemperatureCompensation {
    TemperatureCompensationAntennaFE_OFF_TA_OFF,
    TemperatureCompensationAntennaFE_ON_TA_OFF,
    TemperatureCompensationAntennaFE_OFF_TA_ON,
    TemperatureCompensationAntennaFE_ON_TA_ON
};

enum Polarisation {
    TxHorizontalOnly,
    TxHRxH,
    TxHRxV,
    TxHRxVRxH,
    TxVerticalalOnly,
    TxVRxH,
    TxVRxV,
    TxVRxVRxH,
};

enum CalType {
    TxCal = 0,
    RxCal = 1,
    EPDNCal = 2,
    TACal = 3,
    APDNCal = 4,
    TxHcalIso = 7
};

enum SignalType {
    Echo = 0,
    Noise = 1,
    TxCal = 8,
    RxCal = 9,
    EPDNCal = 10,
    TACal = 11,
    APDNCal = 12,
    TxHcalIso = 15
};

enum ECCNumber {
    Stripmap1 = 1,
    Stripmap2 = 2,
    Stripmap3 = 3,
    Stripmap4 = 4,
    Stripmap5 = 5,
    Stripmap6 = 6,
    InterferometricWideSwath = 8,
    WaveMode = 9,
    Stripmap5S = 10,
    Stripmap1WithoutInterlCal = 11,
    Stripmap2WithoutInterlCal = 12,
    Stripmap3WithoutInterlCal = 13,
    Stripmap4WithoutInterlCal = 14,
    RFCmode = 15,
    TestModeOperBypass = 16,
    ElevationNotchS3 = 17,
    AzimuthNotchS1 = 19,
    AzimuthNotchS2 = 20,
    AzimuthNotchS3 = 21,
    AzimuthNotchS4 = 22,
    AzimuthNotchS5N = 23,
    AzimuthNotchS5S = 24,
    AzimuthNotchS6 = 25,
    Stripmap5NWithoutInterlCal = 26,
    Stripmap5SWithoutInterlCal = 27,
    Stripmap6WithoutInterlCal = 28,
    ElevationNotchS3WithoutInterlCal = 31,
    ExtraWideSwath = 32,
    AzimuthNotchS1WithoutInterlCal = 33,
    AzimuthNotchS3WithoutInterlCal = 34,
    AzimuthNotchS6WithoutInterlCal = 35,
    NoiseCharacterisationS1 = 37,
    NoiseCharacterisationS2 = 38,
    NoiseCharacterisationS3 = 39,
    NoiseCharacterisationS4 = 40,
    NoiseCharacterisationS5N = 41,
    NoiseCharacterisationS5S = 42,
    NoiseCharacterisationS6 = 43,
    NoiseCharacterisationEWS = 44,
    NoiseCharacterisationIWS = 45,
    NoiseCharacterisationWave = 46
};

int ReadSARParam(std::filesystem::path pathToRawData) {
    uint8_t tmp8 = 0;
    uint16_t tmp16 = 0;
    uint32_t tmp32 = 0;
    
    std::ifstream rawData(pathToRawData);
    if (!rawData.is_open()) {
        return -1;
    }

    Sentinel1RawPacket sentinelOneParam;

    while (!rawData.eof()) {
        rawData.read(reinterpret_cast<char*>(&tmp16), 2);
        tmp16 = _byteswap_ushort(tmp16);
        sentinelOneParam.PacketVersionNumber = tmp16 & 0x7;
        sentinelOneParam.PacketType = tmp16 & 0x8;
        sentinelOneParam.SecondaryHeaderFlag = (tmp16 >> 11) & 0x01;
        sentinelOneParam.ProcessID = (tmp16 >> 4) & 0x7F;
        sentinelOneParam.PacketCategory = (tmp16) & 0xF;

        rawData.read(reinterpret_cast<char*>(&tmp16), 2);
        tmp16 = _byteswap_ushort(tmp16);
        sentinelOneParam.SequenceFlags = (tmp16 >> 14);
        sentinelOneParam.PacketSequenceCount = (tmp16 & 0x3FF);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.PacketDataLength), sizeof(sentinelOneParam.PacketDataLength));
        sentinelOneParam.PacketDataLength = _byteswap_ushort(sentinelOneParam.PacketDataLength) + 1;
        if (((sentinelOneParam.PacketDataLength + 6) % 4) != 0) {
            printf("\nERROR: Length not multiple of 4\n");
        }
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.Time), sizeof(sentinelOneParam.Time));
        sentinelOneParam.Time = _byteswap_ulong(sentinelOneParam.Time);
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.FineTime), sizeof(sentinelOneParam.FineTime));
        sentinelOneParam.FineTime = _byteswap_ushort(sentinelOneParam.FineTime);
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SyncMarker), sizeof(sentinelOneParam.SyncMarker));
        sentinelOneParam.SyncMarker = _byteswap_ulong(sentinelOneParam.SyncMarker);
        
        if (sentinelOneParam.SyncMarker != 0x352EF853) {
            std::cerr << "ERROR: Sync marker != 352EF853" << std::endl;
        }
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.DataTakeID), sizeof(sentinelOneParam.DataTakeID));

        rawData.read(reinterpret_cast<char*>(&tmp8), 1);
        sentinelOneParam.ECCNumber = tmp8;

        rawData.read(reinterpret_cast<char*>(&tmp8), 1);
        sentinelOneParam.TestMode = tmp8;
        sentinelOneParam.RxChannelId = tmp8;

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.InstrumentConfigID), sizeof(sentinelOneParam.InstrumentConfigID));
        sentinelOneParam.InstrumentConfigID = _byteswap_ulong(sentinelOneParam.InstrumentConfigID);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.WordIndex), sizeof(sentinelOneParam.WordIndex));
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.WordVal), sizeof(sentinelOneParam.WordVal));
        sentinelOneParam.WordVal = _byteswap_ushort(sentinelOneParam.WordVal);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SpacePacketCount), sizeof(sentinelOneParam.SpacePacketCount));
        sentinelOneParam.SpacePacketCount = _byteswap_ulong(sentinelOneParam.SpacePacketCount);
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.PRICount), sizeof(sentinelOneParam.PRICount));
        sentinelOneParam.PRICount = _byteswap_ulong(sentinelOneParam.PRICount);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.PRICount), sizeof(sentinelOneParam.PRICount));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.BAQMode), sizeof(sentinelOneParam.BAQMode));
        sentinelOneParam.BAQMode = sentinelOneParam.BAQMode & 0x1f;
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.ErrorFlag), sizeof(sentinelOneParam.ErrorFlag));
        
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.BAQBlockLength), sizeof(sentinelOneParam.BAQBlockLength));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.RangeDecimation), sizeof(sentinelOneParam.RangeDecimation));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.RxGain), sizeof(sentinelOneParam.RxGain));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.TxRampRate), sizeof(sentinelOneParam.TxRampRate));
        sentinelOneParam.TxRampRate = _byteswap_ushort(sentinelOneParam.TxRampRate);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.TxPulseStartFreq), sizeof(sentinelOneParam.TxPulseStartFreq));
        sentinelOneParam.TxPulseStartFreq = _byteswap_ushort(sentinelOneParam.TxPulseStartFreq);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.TxPulseLength), sizeof(sentinelOneParam.TxPulseLength));
        sentinelOneParam.TxPulseLength = _byteswap_ulong(sentinelOneParam.TxPulseLength);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.Rank), sizeof(sentinelOneParam.Rank));
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.PulseRepetitionInterval), sizeof(sentinelOneParam.PulseRepetitionInterval));
        sentinelOneParam.PulseRepetitionInterval = _byteswap_ulong(sentinelOneParam.PulseRepetitionInterval);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SamplingWindowStartTime), sizeof(sentinelOneParam.SamplingWindowStartTime));
        sentinelOneParam.SamplingWindowStartTime = _byteswap_ulong(sentinelOneParam.SamplingWindowStartTime);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SamplingWindowLength), sizeof(sentinelOneParam.SamplingWindowLength));
        sentinelOneParam.SamplingWindowLength = _byteswap_ulong(sentinelOneParam.SamplingWindowLength);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.Polarisation), sizeof(sentinelOneParam.Polarisation));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.TemperatureCompensation), sizeof(sentinelOneParam.TemperatureCompensation));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.ElevationBeamAddress), sizeof(sentinelOneParam.ElevationBeamAddress));
        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.AzimuthBeamAddress), sizeof(sentinelOneParam.AzimuthBeamAddress));
        sentinelOneParam.AzimuthBeamAddress = _byteswap_ushort(sentinelOneParam.AzimuthBeamAddress);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SASTestMode), sizeof(sentinelOneParam.SASTestMode));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.CalType), sizeof(sentinelOneParam.CalType));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.CalibrationBeamAddress), sizeof(sentinelOneParam.CalibrationBeamAddress));
        sentinelOneParam.CalibrationBeamAddress = _byteswap_ushort(sentinelOneParam.CalibrationBeamAddress);

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.CalMode), sizeof(sentinelOneParam.CalMode));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.TxPulseNumber), sizeof(sentinelOneParam.TxPulseNumber));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SignalType), sizeof(sentinelOneParam.SignalType));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SwapFlag), sizeof(sentinelOneParam.SwapFlag));

        rawData.read(reinterpret_cast<char*>(&sentinelOneParam.SwathNumber), sizeof(sentinelOneParam.SwathNumber));

        if ((sentinelOneParam.BAQMode == 0x00) && (sentinelOneParam.CalType > 7)) {
            //cposition = bypass(user, NQ, IE, IO, QE, QO);
        }
        if ((sentinelOneParam.BAQMode == 0x0c) && (sentinelOneParam.CalType == 0)) {
            //cposition = packet_decode(user, NQ, IE, IO, QE, QO, brc, &brcpos);
        }
    }; 
}

int next_bit(unsigned char* p, int* cposition, int* bposition)
{
    int bit = ((p[*cposition] >> (*bposition)) & 1);
    (*bposition)--;
    if ((*bposition) < 0) { 
        (*cposition)++; 
        (*bposition) = 7; 
    }
    return(bit);
}

unsigned char get_THIDX(unsigned char* p, int* cposition, int* bposition)
{
    int res = 0;
    int k;
    for (k = 0; k < 8; k++) {
        res = res << 1;
        res += next_bit(p, cposition, bposition);
    }
    return(res);
}

struct sh_code { 
    int sign; 
    int mcode; 
};

static struct sh_code BRC4Func(unsigned char* p, int* cposition, int* bposition) // TODO: never tested !
{
    int hcode, sign;
    struct sh_code sol;
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

static struct sh_code BRC(int BRCn, unsigned char* p, int* cposition, int* bposition)
{
    int hcode;
    int sign;
    int b;
    struct sh_code sol;

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
            return(BRC4Func(p, cposition, bposition));
            printf("\nCheck if BRC4 output is correct\n"); 
            exit(0); 
            break;
        default: 
            printf("ERROR"); exit(-1);
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

// p.78
float BRC0[4] = { 3.,3.,3.16,3.53 };
float BRC1[4] = { 4.,4.,4.08,4.37 };
float BRC2[6] = { 6.,6.,6.,6.15, 6.5,6.88 };
float BRC3[7] = { 9.,9.,9.,9.,9.36,9.50, 10.1 };
float BRC4[9] = { 15.,15.,15.,15.,15.,15., 15.22, 15.50, 16.05 };

// p.79
float NRL0[4] = { .3637,1.0915,1.8208,2.6406 };
float NRL1[5] = { .3042,.9127,1.5216,2.1313,2.8426 };
float NRL2[7] = { .2305,.6916,1.1528,1.6140,2.0754,2.5369,3.1191 };
float NRL3[10] = { .1702,.5107,.8511,1.1916,1.5321,1.8726,2.2131,2.5536,2.8942,3.3744 };
float NRL4[16] = { .1130,.3389,.5649,.7908,1.0167,1.2428,1.4687,1.6947,1.9206,2.1466,2.3725,2.5985,2.8244,3.0504,3.2764,3.6623 };

// p.80
float SF[256] = { 0., 0.630, 1.250, 1.880, 2.510, 3.130, 3.760, 4.390, 5.010, 5.640, 6.270, 6.890, 7.520,
            8.150, 8.770, 9.40, 10.030, 10.650, 11.280, 11.910, 12.530, 13.160, 13.790, 14.410, 15.040,
            15.670, 16.290, 16.920, 17.550, 18.170, 18.80, 19.430, 20.050, 20.680, 21.310, 21.930, 22.560,
            23.190, 23.810, 24.440, 25.070, 25.690, 26.320, 26.950, 27.570, 28.20, 28.830, 29.450, 30.080,
            30.710, 31.330, 31.960, 32.590, 33.210, 33.840, 34.470, 35.090, 35.720, 36.350, 36.970, 37.60,
            38.230, 38.850, 39.480, 40.110, 40.730, 41.360, 41.990, 42.610, 43.240, 43.870, 44.490, 45.120,
            45.750, 46.370, 47., 47.630, 48.250, 48.880, 49.510, 50.130, 50.760, 51.390, 52.010, 52.640,
            53.270, 53.890, 54.520, 55.150, 55.770, 56.40, 57.030, 57.650, 58.280, 58.910, 59.530, 60.160,
            60.790, 61.410, 62.040, 62.980, 64.240, 65.490, 66.740, 68., 69.250, 70.50, 71.760, 73.010,
            74.260, 75.520, 76.770, 78.020, 79.280, 80.530, 81.780, 83.040, 84.290, 85.540, 86.80, 88.050,
            89.30, 90.560, 91.810, 93.060, 94.320, 95.570, 96.820, 98.080, 99.330, 100.580, 101.840, 103.090,
            104.340, 105.60, 106.850, 108.10, 109.350, 110.610, 111.860, 113.110, 114.370, 115.620, 116.870,
            118.130, 119.380, 120.630, 121.890, 123.140, 124.390, 125.650, 126.90, 128.150, 129.410, 130.660,
            131.910, 133.170, 134.420, 135.670, 136.930, 138.180,139.430, 140.690, 141.940, 143.190, 144.450,
            145.70, 146.950, 148.210, 149.460, 150.710, 151.970, 153.220, 154.470, 155.730, 156.980, 158.230,
            159.490, 160.740, 161.990, 163.250, 164.50, 165.750, 167.010, 168.260, 169.510, 170.770, 172.020,
            173.270, 174.530, 175.780, 177.030, 178.290, 179.540, 180.790, 182.050, 183.30, 184.550, 185.810,
            187.060, 188.310, 189.570, 190.820, 192.070, 193.330, 194.580, 195.830, 197.090, 198.340, 199.590,
            200.850, 202.10, 203.350, 204.610, 205.860, 207.110, 208.370, 209.620, 210.870, 212.130, 213.380,
            214.630, 215.890, 217.140, 218.390, 219.650, 220.90, 222.150, 223.410, 224.660, 225.910, 227.170,
            228.420, 229.670, 230.930, 232.180, 233.430, 234.690, 235.940, 237.190, 238.450, 239.70, 240.950,
            242.210, 243.460, 244.710, 245.970, 247.220, 248.470, 249.730, 250.980, 252.230, 253.490, 254.740, 255.990, 255.990 };

void reconstruction(unsigned char* BRCn, unsigned char* THIDXn, struct sh_code* hcode, int NQ, float* result)
{
    int hcode_index = 0, h;
    int BRCindex = 0;
    int inc = 128;
    do
    {
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);
        for (h = 0; h < inc; h++)
        {
            switch (BRCn[BRCindex])
            {
            case 0:
                if (THIDXn[BRCindex] <= 3)
                {
                    if (hcode[hcode_index].mcode < 3)
                        result[hcode_index] = (float)(hcode[hcode_index].sign * hcode[hcode_index].mcode);
                    else
                        result[hcode_index] = (float)(hcode[hcode_index].sign) * BRC0[THIDXn[BRCindex]];
                }
                else  result[hcode_index] = (float)(hcode[hcode_index].sign) * NRL0[hcode[hcode_index].mcode] * SF[THIDXn[BRCindex]];
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
            default: printf("UNHANDLED CASE\n"); exit(-1);
            }
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
}

int packet_decode(unsigned char* p, int NQ, float* IE, float* IO, float* QE, float* QO, char* brc, int* brcpos) // FDBAQ: section 4.4 p.67
{   
    sh_code hcodeIE[52378];
    sh_code hcodeIO[52378];
    sh_code hcodeQE[52378];
    sh_code hcodeQO[52378];
    unsigned char BRCn[410];   // max value p.55: 52378/128=409.2
    unsigned char THIDXn[410];
    int BRCindex;
    int h, hcode_index;
    int cposition = 0, bposition = 7;
    int inc = 128;  // 128 samples until NQ is reached
    
    BRCindex = 0;
    hcode_index = 0;
    do // repeat until end of packet: depends on NQ
    {//msg("\npos:%d:",bposition);
        BRCn[BRCindex] = next_bit(p, &cposition, &bposition) * 4;  // MSb=4
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 2; // then 2
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 1; // then 1 ...
        brc[*brcpos] = BRCn[BRCindex];
        (*brcpos)++;
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeIE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
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
    }        // odd address => +1 
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
    } // start at new position
    if ((cposition & 1) != 0) { 
        cposition++; 
    }     // odd address => +1 
    BRCindex = 0;
    hcode_index = 0;
    do
    {
        THIDXn[BRCindex] = get_THIDX(p, &cposition, &bposition);
        //   THIDXn[BRCindex]=p[cposition]; // 8-bit THIDX
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeQE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
      //      msg("%d",hcodeQE[hcode_index]);
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
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeQO[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
      //      msg("%d",hcodeQO[hcode_index]);
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
    }        // odd address => +1 
    reconstruction(BRCn, THIDXn, hcodeIE, NQ, IE);
    reconstruction(BRCn, THIDXn, hcodeIO, NQ, IO);
    reconstruction(BRCn, THIDXn, hcodeQE, NQ, QE);
    reconstruction(BRCn, THIDXn, hcodeQO, NQ, QO);
    return(cposition);
}