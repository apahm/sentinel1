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

/*
1. The reference replica, Rep(t), is calculated from the amplitude and phase
coefficients that are derived from the extracted chirp replica or from the
nominal chirp.
2. Rep(t) is zero padded to length Nfft and transformed 
to the frequency domain using an FFT of this size, resulting in R(f) of size Nfft.
3.

*/

template<typename T>
T RangeReferenceFunction(T  TXPL, T TXPSF, T TXPRR) {
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
    uint16_t FineTime;
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
    uint8_t RxGain;
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

int ReadSARParam(std::filesystem::path pathToRawData) {
    int f, res; 
    uint16_t c, tmp16, NQ = 0;
    uint32_t tmp32, Time;
    int Secondary;
    int Count, DataLen, PID, PCAT, Sequence;
    unsigned char tablo[65536], BAQ, Swath, Typ;
    int cal_p, cposition, brcpos;
    unsigned char* user;
    float IE[52378]; // results
    float IO[52378];
    float QE[52378];
    float QO[52378];
    char brc[52378];
    
    int result;
    int numline = 0;
    int file_swath_number = 0;
    int file_nq = 0; // for new file creation
    std::ifstream rawData(pathToRawData);
    if (!rawData.is_open()) {
        return -1;
    }

    do {
        rawData.read(reinterpret_cast<char*>(&c), 2);
                
        c = _byteswap_ushort(c);
        Secondary = (c >> 11) & 0x01;
        PID = (c >> 4) & 0x7f;
        PCAT = (c) & 0xf;
        printf("%04x: %d(1)\t%d(65)\t%d(12)\t", c, Secondary, PID, PCAT);

        rawData.read(reinterpret_cast<char*>(&c), 2);
        c = _byteswap_ushort(c);
        Sequence = (c >> 14);
        Count = (c & 0x3f);

        rawData.read(reinterpret_cast<char*>(&c), 2);
        DataLen = _byteswap_ushort(c) + 1;
        printf("%04x: %x(3)\tCount=%02d\tLen=%d(61..65533)\n", c, Sequence, Count, DataLen);
        
        if (((DataLen + 6) % 4) != 0) {
            printf("\nERROR: Length not multiple of 4\n");
        }
        
        // Coarse Time
        rawData.read(reinterpret_cast<char*>(&Time), 4);
        Time = _byteswap_ulong(Time);
        printf("\tTime: %d", Time);               
        
        // Fine Time   
        tmp16 = *(uint16_t*)(tablo + 4); 
        tmp16 = _byteswap_ushort(tmp16);
        printf(":%d", tmp16);                     
        
        // Sync Marker 
        tmp32 = *(uint32_t*)(tablo + 6); 
        tmp32 = _byteswap_ulong(tmp32);
        printf("\t%08x(352EF853)", tmp32);       
        
        if (tmp32 != 0x352EF853) printf("\nERROR: Sync marker != 352EF853");
        
        tmp32 = *(uint32_t*)(tablo + 10); 
        tmp32 = _byteswap_ulong(tmp32);
        //printf(" %08x",tmp32);                 // Data Take ID= begin+16 (TBD)
        //printf("\t%hhx",*(uint8_t*)tablo+14); // ECC Number  = begin+20  ?! /!\ ?!
        //printf("\t%hhx",*(uint8_t*)tablo+15); // TestMode/RXID=begin+21
        tmp32 = *(uint32_t*)(tablo + 16); 
        tmp32 = _byteswap_ulong(tmp32);
        //printf("\t%08x",tmp32);                  // Config ID   = begin+22 (TBD)
        // page 22

        // Sub-commutation ancillary data: PVT/attitude will be accumulated as 42 words along headers
        // p.27: Counter+value                    (3 bytes)
        printf("\tWordIndex=%hhx", *(uint8_t*)(tablo + 20)); // Word index will increment from 1 to 0x40 to fill
        tmp16 = *(uint16_t*)(tablo + 21); 
        tmp16 = _byteswap_ushort(tmp16);  // the array described in p.23 with ... vvv
        printf("\tWordVal=%hx", tmp16);                   // Word value

        // Space packet count (4 bytes)
        tmp32 = *(uint32_t*)(tablo + 23); 
        tmp32 = _byteswap_ulong(tmp32);
        printf("\t%08x", tmp32);                  
        
        // PRI count (4 bytes)
        tmp32 = *(uint32_t*)(tablo + 27); 
        tmp32 = _byteswap_ulong(tmp32);
        printf(" %08x", tmp32);                   
        
        // BAQ mode 0x0c=FDBAQ mode0 nominal p.33: 
        BAQ = (*(uint8_t*)(tablo + 31)) & 0x1f;       
        printf(" BAQ=%02x(c)", (*(uint8_t*)(tablo + 31)));
        // p.32: RADAR Configuration Support     (28 bytes) -> total=31+28=59
        printf(" BlockLen=%hhx(1F)\n", (*(uint8_t*)(tablo + 32)));

        printf("\tDecim=%hhx", *(uint8_t*)(tablo + 34));     // RGDEC
        tmp16 = *(uint16_t*)(tablo + 36); 
        tmp16 = _byteswap_ushort(tmp16);  // Tx Pulse Ramp Rate
        if ((tmp16 & 0x8000) == 0) 
            printf("\tTXPRR=v%x=", tmp16); 
        else 
            printf("\tTXPRR=^%x=", tmp16);
        printf("%d", tmp16 & 0x7fff);
        tmp16 = *(uint16_t*)(tablo + 38); 
        tmp16 = _byteswap_ushort(tmp16);  // Tx Pulse Start Freq
        if ((tmp16 & 0x8000) == 0) 
            printf("\tTXPSF=-0x%x", tmp16); 
        else 
            printf("\tTXPSF=+0x%x", tmp16); // 0 negative, 1 positive ?!
        printf("=%d ", tmp16 & 0x7fff);             // Word value
        tmp32 = *(uint32_t*)(tablo + 40); 
        tmp32 = ntohl(tmp32) >> 8; // keep last 24 bits
        printf("TXPL=%08x=%d ", tmp32, tmp32);     // TX Pulse length (3 bytes)
        tmp32 = *(uint32_t*)(tablo + 44); 
        tmp32 = ntohl(tmp32) >> 8; // keep last 24 bits
        printf("PRI=%08x=%d", tmp32, tmp32);       // PRI, needed for Doppler centroid analysis (3 bytes)

        cal_p = ((*(uint8_t*)(tablo + 53)) >> 4) & 0x07;
        printf(" Polar=%x", cal_p);              // SSB Data calibration (p.47) 
        Typ = (*(uint8_t*)(tablo + 57)); Typ = Typ >> 4; // p.52: signal type
        printf(" Typ=%hhx(0)", Typ);
        Swath = (*(uint8_t*)(tablo + 58));          // p.54: swath number
        printf(" Swath=%hhx", Swath);
        // p.54 RADAR Sample Count (2 bytes)
        NQ = *(uint16_t*)(tablo + 59); 
        NQ = _byteswap_ushort(NQ); // number of quads NQ => Nsamples=2xNQ 
        printf(" NQ=%d\n", NQ);
        // if (NQ==0) fprintf(result,"%d\n",2*NQ);
        // tablo+61 = n/a (p.54: index 67)

        // p.56 User Data Field -- DataLen-62 ; 4 sections with IE, IO, QE, QO (NOT interleaved)
        user = (uint8_t*)(tablo + 62);  // User Data Field starts @ end of Secondary Header
        // p.58: format D is nominally used to output radar echo data = Decimation + FDBAQ
        // NO BYTESWAP SINCE WE WORK BIT BY BIT: Keep bytes in read() order
        // calibration (p.52: Typ>7 for cal, and p.33: BAQMOD=0=BYPASS for cal)
        if ((BAQ == 0x00) && (Typ > 7)) {

            cposition = bypass(user, NQ, IE, IO, QE, QO);

            printf(", finished processing %d calibration\n", DataLen - 62);

            if ((DataLen - 62 - cposition) > 2) {
                printf("Not enough data processed: DataLen %d v.s. cposition %d\n", DataLen - 62, cposition); exit(-1);
            }
        }
        // echo data (p.52: Typ=0 for echo data and p.33 BAQ Mode=0x0C => nominal FDBAQ p.67)
        if ((BAQ == 0x0c) && (Typ == 0)) {
            brcpos = 0;
            cposition = packet_decode(user, NQ, IE, IO, QE, QO, brc, &brcpos);

            printf(", finished processing %d echo\n", DataLen - 62);
            if ((DataLen - 62 - cposition) > 2) { 
                printf("Not enough data processed: DataLen %d v.s. cposition %d\n", DataLen - 62, cposition); 
                exit(-1); 
            }
        }
        numline++;
    } while ((res > 0)); // until EOF
}

int next_bit(unsigned char* p, int* cposition, int* bposition)
{
    int bit = ((p[*cposition] >> (*bposition)) & 1);
    (*bposition)--;
    if ((*bposition) < 0) { (*cposition)++; (*bposition) = 7; }
    return(bit);
}

static struct sh_code BRC4(unsigned char* p, int* cposition, int* bposition) // TODO: never tested !
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
    case 0: BRCn = 3; break; // number of steps to reach the leaves BRC0
    case 1: BRCn = 4; break; // number of steps to reach the leaves BRC1
    case 2: BRCn = 6; break; // number of steps to reach the leaves BRC2
    case 3: BRCn = 9; break; // number of steps to reach the leaves BRC3
    case 4: return(BRC4(p, cposition, bposition)); printf("\nCheck if BRC4 output is correct\n"); exit(0); break;
    default: printf("ERROR"); exit(-1);
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

unsigned char get_THIDX(unsigned char* p, int* cposition, int* bposition)
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

int packet_decode(unsigned char* p, int NQ, float* IE, float* IO, float* QE, float* QO, char* brc, int* brcpos) // FDBAQ: section 4.4 p.67
{// IE 1st 3 bits = BRC
 // QE first 8 bits = THIDX
    struct sh_code hcodeIE[52378];
    struct sh_code hcodeIO[52378];
    struct sh_code hcodeQE[52378];
    struct sh_code hcodeQO[52378];
    unsigned char BRCn[410];   // max value p.55: 52378/128=409.2
    unsigned char THIDXn[410];
    int BRCindex;
    int h, hcode_index;
    int cposition = 0, bposition = 7;
    int inc = 128;  // 128 samples until NQ is reached
    msg("\nstarting IE\n");
    BRCindex = 0;
    hcode_index = 0;
    do // repeat until end of packet: depends on NQ
    {//msg("\npos:%d:",bposition);
        BRCn[BRCindex] = next_bit(p, &cposition, &bposition) * 4;  // MSb=4
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 2; // then 2
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 1; // then 1 ...
#ifdef dump_payload
        msg("\n");
#endif
        msg("%d>", BRCn[BRCindex]);
        brc[*brcpos] = BRCn[BRCindex];
        (*brcpos)++;
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeIE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
#ifdef dump_payload
            msg("%d", hcodeIE[hcode_index]);
#endif
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    msg("\nIE finished, starting IO @ %d\n", cposition);
    inc = 128;
    if (bposition != 7) { msg("bposition=%d->7\n", bposition); bposition = 7; cposition++; } // start at new position
    if ((cposition & 1) != 0) { msg("cposition=%d++\n", cposition); cposition++; }        // odd address => +1 
    BRCindex = 0;
    hcode_index = 0;
    do
    {
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeIO[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
      //      msg("%d",hcodeIO[hcode_index]);
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    msg("IO finished, starting QE @ %d\n", cposition);
    inc = 128;
    if (bposition != 7) { msg("bposition=%d\n", bposition); bposition = 7; cposition++; } // start at new position
    if ((cposition & 1) != 0) { msg("cposition=%d++\n", cposition); cposition++; }     // odd address => +1 
    BRCindex = 0;
    hcode_index = 0;
    do
    {
        THIDXn[BRCindex] = get_THIDX(p, &cposition, &bposition);
        //   THIDXn[BRCindex]=p[cposition]; // 8-bit THIDX
        msg("#%d", THIDXn[BRCindex]);
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
        {
            hcodeQE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); // 128 samples with same BRC
      //      msg("%d",hcodeQE[hcode_index]);
            hcode_index++;
        }
        BRCindex++;
    } while (hcode_index < NQ);
    msg("\nQE finished, starting QO @ %d\n", cposition);
    inc = 128;
    if (bposition != 7) { msg("bposition=%d\n", bposition); bposition = 7; cposition++; } // start at new position
    if ((cposition & 1) != 0) { msg("cposition=%d++\n", cposition); cposition++; }     // odd address => +1 
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
    if (bposition != 7) { msg("bposition=%d->7\n", bposition); bposition = 7; cposition++; } // start at new position
    if ((cposition & 1) != 0) { msg("cposition=%d++\n", cposition); cposition++; }        // odd address => +1 
    msg("finished at %d", cposition);
    reconstruction(BRCn, THIDXn, hcodeIE, NQ, IE);
    reconstruction(BRCn, THIDXn, hcodeIO, NQ, IO);
    reconstruction(BRCn, THIDXn, hcodeQE, NQ, QE);
    reconstruction(BRCn, THIDXn, hcodeQO, NQ, QO);
    return(cposition);
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
float SF[256] = { 0., 0.630, 1.250, 1.880, 2.510, 3.130, 3.760, 4.390, 5.010, 5.640, 6.270, 6.890, 7.520, 8.150, 8.770, 9.40, 10.030, 10.650, 11.280, 11.910, 12.530, 13.160, 13.790, 14.410, 15.040, 15.670, 16.290, 16.920, 17.550, 18.170, 18.80, 19.430, 20.050, 20.680, 21.310, 21.930, 22.560, 23.190, 23.810, 24.440, 25.070, 25.690, 26.320, 26.950, 27.570, 28.20, 28.830, 29.450, 30.080, 30.710, 31.330, 31.960, 32.590, 33.210, 33.840, 34.470, 35.090, 35.720, 36.350, 36.970, 37.60, 38.230, 38.850, 39.480, 40.110, 40.730, 41.360, 41.990, 42.610, 43.240, 43.870, 44.490, 45.120, 45.750, 46.370, 47., 47.630, 48.250, 48.880, 49.510, 50.130, 50.760, 51.390, 52.010, 52.640, 53.270, 53.890, 54.520, 55.150, 55.770, 56.40, 57.030, 57.650, 58.280, 58.910, 59.530, 60.160, 60.790, 61.410, 62.040, 62.980, 64.240, 65.490, 66.740, 68., 69.250, 70.50, 71.760, 73.010, 74.260, 75.520, 76.770, 78.020, 79.280, 80.530, 81.780, 83.040, 84.290, 85.540, 86.80, 88.050, 89.30, 90.560, 91.810, 93.060, 94.320, 95.570, 96.820, 98.080, 99.330, 100.580, 101.840, 103.090, 104.340, 105.60, 106.850, 108.10, 109.350, 110.610, 111.860, 113.110, 114.370, 115.620, 116.870, 118.130, 119.380, 120.630, 121.890, 123.140, 124.390, 125.650, 126.90, 128.150, 129.410, 130.660, 131.910, 133.170, 134.420, 135.670, 136.930, 138.180, 139.430, 140.690, 141.940, 143.190, 144.450, 145.70, 146.950, 148.210, 149.460, 150.710, 151.970, 153.220, 154.470, 155.730, 156.980, 158.230, 159.490, 160.740, 161.990, 163.250, 164.50, 165.750, 167.010, 168.260, 169.510, 170.770, 172.020, 173.270, 174.530, 175.780, 177.030, 178.290, 179.540, 180.790, 182.050, 183.30, 184.550, 185.810, 187.060, 188.310, 189.570, 190.820, 192.070, 193.330, 194.580, 195.830, 197.090, 198.340, 199.590, 200.850, 202.10, 203.350, 204.610, 205.860, 207.110, 208.370, 209.620, 210.870, 212.130, 213.380, 214.630, 215.890, 217.140, 218.390, 219.650, 220.90, 222.150, 223.410, 224.660, 225.910, 227.170, 228.420, 229.670, 230.930, 232.180, 233.430, 234.690, 235.940, 237.190, 238.450, 239.70, 240.950, 242.210, 243.460, 244.710, 245.970, 247.220, 248.470, 249.730, 250.980, 252.230, 253.490, 254.740, 255.990, 255.990 };

void reconstruction(unsigned char* BRCn, unsigned char* THIDXn, struct sh_code* hcode, int NQ, float* result)
{
    int hcode_index = 0, h;
    int BRCindex = 0;
    int inc = 128;
    do
    {
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      // smaller increment to match NQ
        for (h = 0; h < inc; h++) // p.68: 128 HCodes
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

int bypass(unsigned char* p, int NQ, float* IE, float* IO, float* QE, float* QO) // bypass
{
    int NW = (10 * NQ) / 16 * 2;
    short res;
    int pos = 0, sign, index = 0;
    printf("\nNQ=%d -> NW=%d\n", NQ, NW);
    while (index < (NQ))  // 8*5=40 : 4 10=bit values for 5 bytes
    {// printf("%02hhx %02hhx %02hhx %02hhx %02hhx = ",p[pos+0],p[pos+1],p[pos+2],p[pos+3],p[pos+4]);
        if (index < NQ)
        {
            res = (short)(p[pos + 0] & 0xff) * 4 + (short)(p[pos + 1] >> 6);       // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IE[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 1] & 0x3f) * 4 * 4 + (short)(p[pos + 2] >> 4);     // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IE[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 2] & 0x0f) * 4 * 4 * 4 + (short)(p[pos + 3] >> 2);   // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IE[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 3] & 0x03) * 4 * 4 * 4 * 4 + (short)(p[pos + 4]);    // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IE[index] = (float)res;
            index++;
        }
        if (index < NQ) pos += 5;
    }
    index = 0;
    pos = NW + 2; // next short
    while (index < (NQ))  // 8*5=40 : 4 10=bit values for 5 bytes
    {// printf("%02hhx %02hhx %02hhx %02hhx %02hhx = ",p[pos+0],p[pos+1],p[pos+2],p[pos+3],p[pos+4]);
        if (index < NQ)
        {
            res = (short)(p[pos + 0] & 0xff) * 4 + (short)(p[pos + 1] >> 6);       // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IO[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 1] & 0x3f) * 4 * 4 + (short)(p[pos + 2] >> 4);     // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IO[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 2] & 0x0f) * 4 * 4 * 4 + (short)(p[pos + 3] >> 2);   // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IO[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 3] & 0x03) * 4 * 4 * 4 * 4 + (short)(p[pos + 4]);    // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            IO[index] = (float)res;
            index++;
        }
        if (index < NQ) pos += 5;
    }
    printf("\npos=%d index=%d nw=%d\n", pos, index, (2 * NW + 2));
    index = 0;
    pos = 2 * (NW + 2); // next short
    while (index < (NQ))  // 8*5=40 : 4 10=bit values for 5 bytes
    {// printf("%02hhx %02hhx %02hhx %02hhx %02hhx = ",p[pos+0],p[pos+1],p[pos+2],p[pos+3],p[pos+4]);
        if (index < NQ)
        {
            res = (short)(p[pos + 0] & 0xff) * 4 + (short)(p[pos + 1] >> 6);       // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QE[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 1] & 0x3f) * 4 * 4 + (short)(p[pos + 2] >> 4);     // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QE[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 2] & 0x0f) * 4 * 4 * 4 + (short)(p[pos + 3] >> 2);   // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QE[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 3] & 0x03) * 4 * 4 * 4 * 4 + (short)(p[pos + 4]);    // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QE[index] = (float)res;
            index++;
        }
        if (index < NQ) pos += 5;
    }
    index = 0;
    pos = 3 * (NW + 2); // next short
    while (index < (NQ))  // 8*5=40 : 4 10=bit values for 5 bytes
    {// printf("%02hhx %02hhx %02hhx %02hhx %02hhx = ",p[pos+0],p[pos+1],p[pos+2],p[pos+3],p[pos+4]);
        if (index < NQ)
        {
            res = (short)(p[pos + 0] & 0xff) * 4 + (short)(p[pos + 1] >> 6);       // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QO[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 1] & 0x3f) * 4 * 4 + (short)(p[pos + 2] >> 4);     // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QO[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 2] & 0x0f) * 4 * 4 * 4 + (short)(p[pos + 3] >> 2);   // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QO[index] = (float)res;
            index++;
        }
        if (index < NQ)
        {
            res = (short)(p[pos + 3] & 0x03) * 4 * 4 * 4 * 4 + (short)(p[pos + 4]);    // printf("%hx=",res);   
            sign = (res & 0x200); res = res & (0x1ff); if (sign > 0) res = -res; // printf("%03d\t",res);
            QO[index] = (float)res;
            index++;
        }
        if (index < NQ) pos += 5;
    }
    printf("\npos=%d index=%d\n", pos, index);
    pos = 4 * (NW + 2); // next short
    return(pos);
}

int main(int argc, char *argv[])
{

    return 0;
}
