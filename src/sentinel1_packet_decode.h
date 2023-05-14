#pragma once
#include <iostream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <filesystem>
#include <fstream>
#include "../Eigen/Dense"

enum class ECCNumber {
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

enum class TestMode {
    Default = 0,
    GroundTestingOnlyModeOne = 0x4,
    GroundTestingOnlyModeTwo = 0x5,
    Oper = 0x6,
    Bypass = 0x7
};

enum class BAQMode {
    BypassMode = 0,
    BAQ3BitMode = 3,
    BAQ4BitMode = 4,
    BAQ5BitMode = 5,
    FDBAQMode0 = 12,
    FDBAQMode1 = 13,
    FDBAQMode2 = 14
};

enum class RangeDecimation {
    Filter0,
    Filter1,
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

enum class TemperatureCompensation {
    TemperatureCompensationAntennaFE_OFF_TA_OFF,
    TemperatureCompensationAntennaFE_ON_TA_OFF,
    TemperatureCompensationAntennaFE_OFF_TA_ON,
    TemperatureCompensationAntennaFE_ON_TA_ON
};

enum class Polarisation {
    TxHorizontalOnly,
    TxHRxH,
    TxHRxV,
    TxHRxVRxH,
    TxVerticalalOnly,
    TxVRxH,
    TxVRxV,
    TxVRxVRxH,
};

enum class CalType {
    TxCal = 0,
    RxCal = 1,
    EPDNCal = 2,
    TACal = 3,
    APDNCal = 4,
    TxHcalIso = 7
};

enum class SignalType {
    Echo = 0,
    Noise = 1,
    TxCal = 8,
    RxCal = 9,
    EPDNCal = 10,
    TACal = 11,
    APDNCal = 12,
    TxHcalIso = 15
};

/*Position, Velocity, Time data from raw Sentinel 1 packet*/

struct PositionVelocityTime {
    double XAxisPositionECEF = 0.0;
    double YAxisPositionECEF = 0.0;
    double ZAxisPositionECEF = 0.0;
    float XvelocityECEF = 0.0;
    float YvelocityECEF = 0.0;
    float ZvelocityECEF = 0.0;
    uint16_t DataStampOne = 0;
    uint16_t DataStampTwo = 0;
    uint16_t DataStampThree = 0;
    uint16_t DataStampFour = 0;
};

/*
* The reference frames for the S / C Attitude Quaternions and the S / C angular rates as given in[NRD 02] are as
follows :
1) the S / C Attitude Quaternions represent the S / C attitude in the ECI � J2000 Reference Frame. Q0 is the real
component and Q1, Q2, Q3 are the vector components of the Attitude Quaternion.
2) the S / C inertial angular rate vector is measured in the Body Fixed Reference Frame.
*/

struct Attitude {
    float Q0AttitudeQuaternion = 0.0;
    float Q1AttitudeQuaternion = 0.0;
    float Q2AttitudeQuaternion = 0.0;
    float Q3AttitudeQuaternion = 0.0;
    float XangularRate = 0.0;
    float YangularRate = 0.0;
    float ZangularRate = 0.0;
    uint16_t DataStampOne = 0;
    uint16_t DataStampTwo = 0;
    uint16_t DataStampThree = 0;
    uint16_t DataStampFour = 0;
    uint8_t AOCSOperationalMode = 0;
    bool RollErrorStatus = false;
    bool PitchErrorStatus = false;
    bool YawErrorStatus = false;
};

/*
Raw data analysis is required in order to perform corrections of the I and Q channels
of the raw signal data.
    I/Q bias removal
    I/Q gain imbalance correction
    I/Q non - orthogonality correction
For Sentinel-1 however, the instrument�s receive module performs the demodulation
in the digital domain, therefore the I/Q gain imbalance and I/Q non-orthogonality
corrections are no longer necessary.
*/

struct RawDataAnalysis {
    double meanOfRawDataI = 0.0;
    double meanOfRawDataQ = 0.0;
    double standardDeviationsOfRawDataI = 0.0;
    double standardDeviationsOfRawDataQ = 0.0;
    double IQGainImbalance = 0.0;
    double IQGainImbalanceLowerBounds = 0.0;
    double IQGainImbalanceUpperBounds = 0.0;
    double IQQuadratureDeparture = 0.0;
    double IQQuadratureDepartureLowerBounds = 0.0;
    double IQQuadratureDepartureUpperBounds = 0.0;
    bool IBiasSignificanceFlag = false;
    bool QBiasSignificanceFlag = false;
    bool IQGainSignificanceFlag = false;
    bool IQQuadratureDepartureSignificantFlag = false;
};

/*
Each Space Packet generated by the Instrument contains the complete SAR data acquired in one PRI.
The standard limits the maximum packet size to (65536+6) octets but, due to decimation and BAQ
compression of SAR data in nominal operation the packet size will stay well below this limit.
The format of the Space Packet is described in section 3. In each packet a SAR ancillary data field is included.
The ancillary data provide the information how to interprete, decode and process the SAR radar data in the packet.
In addition, they provide information about Instrument status and configuration at the moment of data acquisition (digitizing).
General approach here is to provide appropriate information for ground to support SAR image decoding and processing.
*/

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
    float FineTime;
    uint32_t SyncMarker;
    uint32_t DataTakeID;
    ECCNumber ECCNumber;
    TestMode TestMode;
    bool RxChannelId;
    uint32_t InstrumentConfigID;
    uint8_t WordIndex;
    uint16_t WordVal;
    uint32_t SpacePacketCount;
    uint32_t PRICount;
    BAQMode BAQMode;
    bool ErrorFlag;
    uint32_t BAQBlockLength;
    RangeDecimation RangeDecimation;
    float RxGain;
    float TxRampRate;
    float TxPulseStartFreq;
    float TxPulseLength;
    uint8_t Rank;
    float PulseRepetitionInterval;
    float SamplingWindowStartTime;
    float SamplingWindowLength;
    Polarisation Polarisation;
    bool SSBFlag;
    TemperatureCompensation  TemperatureCompensation;
    uint8_t ElevationBeamAddress;
    uint16_t AzimuthBeamAddress;
    bool SASTestMode;
    CalType CalType;
    uint16_t CalibrationBeamAddress;
    uint8_t CalMode;
    uint8_t TxPulseNumber;
    SignalType SignalType;
    bool SwapFlag;
    uint8_t SwathNumber;
    uint16_t NumberOfQuads;
};

struct ShCode { 
    int sign; 
    int mcode; 
};

class Sentinel1PacketDecode
{
public:
    Sentinel1PacketDecode();
	~Sentinel1PacketDecode();
	float RangeDecimationToSampleRate(RangeDecimation& rangeDecimation);
    float calcRxGain(uint8_t rawRxGain);
    float calcTxPulseRampRate(uint16_t rawTxPulseRampRate);
    float calcTxPulseStartFreq(uint16_t rawTxPulseStartFreq, float TxPulseRampRate);
    uint32_t calcTxPulseLength(uint32_t rawTxPulseStartFreq);
    float calcPulseRepetitionInterval(uint32_t rawPulseRepetitionInterval);
    float calcSamplingWindowStartTime(uint32_t rawSamplingWindowStartTime);
    float calcSamplingWindowLength(uint32_t rawSamplingWindowLength);
    float calcFineTime(uint16_t rawFineTime);
    uint32_t calcInstrumentConfigID(uint32_t rawInstrumentConfigID);
    uint32_t read24Bit(std::ifstream& f);
    int ReadSARParam(std::filesystem::path pathToRawData);
    void reconstruction(unsigned char* BRCn, unsigned char* THIDXn, ShCode* hcode, int NQ, float* result);
    int packet_decode(unsigned char* p, int NQ, float* IE, float* IO, float* QE, float* QO, char* brc, int* brcpos);
    unsigned char get_THIDX(unsigned char* p, int* cposition, int* bposition);
    int next_bit(unsigned char* p, int* cposition, int* bposition);
    ShCode BRC(int BRCn, unsigned char* p, int* cposition, int* bposition);
    ShCode BRC_4(unsigned char* p, int* cposition, int* bposition);

    std::vector<Sentinel1RawPacket> header;
    std::vector<PositionVelocityTime> positionVelocityTime;
    std::vector<Attitude> attitude;

    Eigen::MatrixXd matrix;

    float BRC0[4] = { 3.,3.,3.16,3.53 };
    float BRC1[4] = { 4.,4.,4.08,4.37 };
    float BRC2[6] = { 6.,6.,6.,6.15, 6.5,6.88 };
    float BRC3[7] = { 9.,9.,9.,9.,9.36,9.50, 10.1 };
    float BRC4[9] = { 15.,15.,15.,15.,15.,15., 15.22, 15.50, 16.05 };

    float NRL0[4] = { .3637,1.0915,1.8208,2.6406 };
    float NRL1[5] = { .3042,.9127,1.5216,2.1313,2.8426 };
    float NRL2[7] = { .2305,.6916,1.1528,1.6140,2.0754,2.5369,3.1191 };
    float NRL3[10] = { .1702,.5107,.8511,1.1916,1.5321,1.8726,2.2131,2.5536,2.8942,3.3744 };
    float NRL4[16] = { .1130,.3389,.5649,.7908,1.0167,1.2428,1.4687,1.6947,1.9206,2.1466,2.3725,2.5985,2.8244,3.0504,3.2764,3.6623 };

    float SF[256] = { 0., 0.630, 1.250, 1.880, 2.510, 3.130, 3.760, 4.390, 5.010, 5.640, 6.270, 
        6.890, 7.520, 8.150, 8.770, 9.40, 10.030, 10.650, 11.280, 11.910, 12.530, 13.160, 13.790, 
        14.410, 15.040, 15.670, 16.290, 16.920, 17.550, 18.170, 18.80, 19.430, 20.050, 20.680, 
        21.310, 21.930, 22.560, 23.190, 23.810, 24.440, 25.070, 25.690, 26.320, 26.950, 27.570,
        28.20, 28.830, 29.450, 30.080, 30.710, 31.330, 31.960, 32.590, 33.210, 33.840, 34.470, 
        35.090, 35.720, 36.350, 36.970, 37.60, 38.230, 38.850, 39.480, 40.110, 40.730, 41.360, 
        41.990, 42.610, 43.240, 43.870, 44.490, 45.120, 45.750, 46.370, 47., 47.630, 48.250, 
        48.880, 49.510, 50.130, 50.760, 51.390, 52.010, 52.640, 53.270, 53.890, 54.520, 55.150, 
        55.770, 56.40, 57.030, 57.650, 58.280, 58.910, 59.530, 60.160, 60.790, 61.410, 62.040, 
        62.980, 64.240, 65.490, 66.740, 68., 69.250, 70.50, 71.760, 73.010, 74.260, 75.520, 76.770, 
        78.020, 79.280, 80.530, 81.780, 83.040, 84.290, 85.540, 86.80, 88.050, 89.30, 90.560, 91.810, 93.060,
        94.320, 95.570, 96.820, 98.080, 99.330, 100.580, 101.840, 103.090, 104.340, 105.60, 106.850, 
        108.10, 109.350, 110.610, 111.860, 113.110, 114.370, 115.620, 116.870, 118.130, 119.380, 
        120.630, 121.890, 123.140, 124.390, 125.650, 126.90, 128.150, 129.410, 130.660, 131.910, 
        133.170, 134.420, 135.670, 136.930, 138.180, 139.430, 140.690, 141.940, 143.190, 144.450, 
        145.70, 146.950, 148.210, 149.460, 150.710, 151.970, 153.220, 154.470, 155.730, 156.980, 
        158.230, 159.490, 160.740, 161.990, 163.250, 164.50, 165.750, 167.010, 168.260, 169.510, 
        170.770, 172.020, 173.270, 174.530, 175.780, 177.030, 178.290, 179.540, 180.790, 182.050, 
        183.30, 184.550, 185.810, 187.060, 188.310, 189.570, 190.820, 192.070, 193.330, 194.580, 
        195.830, 197.090, 198.340, 199.590, 200.850, 202.10, 203.350, 204.610, 205.860, 207.110, 
        208.370, 209.620, 210.870, 212.130, 213.380, 214.630, 215.890, 217.140, 218.390, 219.650, 
        220.90, 222.150, 223.410, 224.660, 225.910, 227.170, 228.420, 229.670, 230.930, 232.180, 
        233.430, 234.690, 235.940, 237.190, 238.450, 239.70, 240.950, 242.210, 243.460, 244.710, 
        245.970, 247.220, 248.470, 249.730, 250.980, 252.230, 253.490, 254.740, 255.990, 255.990 };
};

