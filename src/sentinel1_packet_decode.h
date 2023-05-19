#pragma once
#include <iostream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <filesystem>
#include <fstream>
#include <array>
#include <complex>
#include "ipp.h"

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
1) the S / C Attitude Quaternions represent the S / C attitude in the ECI – J2000 Reference Frame. Q0 is the real
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
For Sentinel-1 however, the instrument’s receive module performs the demodulation
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

struct State {
    int _start_time;
    int _ancillary_data_index;
    std::vector<size_t> _header_offset;
    std::vector<std::array<uint8_t, 62 + 6>> _header_data;
    size_t _mmap_filesize;
    uint8_t* _mmap_data;
    char const* _filename;
    std::vector<uint8_t> file;
};

enum { 
    MAX_NUMBER_QUADS = 52378 
}; 

struct sequential_bit_t {
    size_t current_bit_count;
    uint8_t* data;
};

typedef struct sequential_bit_t sequential_bit_t;

class Sentinel1PacketDecode
{
public:
    Sentinel1PacketDecode();
	~Sentinel1PacketDecode();
    float calcRxGain(uint8_t rawRxGain);
    float calcTxPulseRampRate(uint16_t rawTxPulseRampRate);
    float calcTxPulseStartFreq(uint16_t rawTxPulseStartFreq, float TxPulseRampRate);
    float calcTxPulseLength(uint32_t rawTxPulseStartFreq);
    float calcPulseRepetitionInterval(uint32_t rawPulseRepetitionInterval);
    float calcSamplingWindowStartTime(uint32_t rawSamplingWindowStartTime);
    float calcSamplingWindowLength(uint32_t rawSamplingWindowLength);
    float calcFineTime(uint16_t rawFineTime);
    uint32_t calcInstrumentConfigID(uint32_t rawInstrumentConfigID);
    uint32_t read24Bit(std::ifstream& f);
    int ReadSARParam(std::filesystem::path pathToRawData);
    int initDecodePacket(Ipp32fc* output, Sentinel1RawPacket& sentinelOneParam);

    std::vector<Sentinel1RawPacket> header;
    std::vector<PositionVelocityTime> positionVelocityTime;
    std::vector<Attitude> attitude;
    std::vector<std::vector<Ipp32fc>> out;
    // table 5.2-1 simple reconstruction parameter values B
    const std::array<const float, 4> table_b0 = { 3, 3, (3.160f), (3.530f) };
    const std::array<const float, 4> table_b1 = { 4, 4, (4.080f), (4.370f) };
    const std::array<const float, 6> table_b2 = { 6,        6,       6,
                                                 (6.150f), (6.50f), (6.880f) };
    const std::array<const float, 7> table_b3 = { 9,        9,       9,       9,
                                                 (9.360f), (9.50f), (10.10f) };
    const std::array<const float, 9> table_b4 = {
        15, 15, 15, 15, 15, 15, (15.220f), (15.50f), (16.050f) };
    // table 5.2-2 normalized reconstruction levels
    const std::array<const float, 4> table_nrl0 = { (0.36370f), (1.09150f),
                                                   (1.82080f), (2.64060f) };
    const std::array<const float, 5> table_nrl1 = {
        (0.30420f), (0.91270f), (1.52160f), (2.13130f), (2.84260f) };
    const std::array<const float, 7> table_nrl2 = {
        (0.23050f), (0.69160f), (1.15280f), (1.6140f),
        (2.07540f), (2.53690f), (3.11910f) };
    const std::array<const float, 10> table_nrl3 = {
        (0.17020f), (0.51070f), (0.85110f), (1.19160f), (1.53210f),
        (1.87260f), (2.21310f), (2.55360f), (2.89420f), (3.37440f) };
    const std::array<const float, 16> table_nrl4 = {
        (0.1130f),  (0.33890f), (0.56490f), (0.79080f), (1.01670f), (1.24280f),
        (1.46870f), (1.69470f), (1.92060f), (2.14660f), (2.37250f), (2.59850f),
        (2.82440f), (3.05040f), (3.27640f), (3.66230f) };
    // table 5.2-3 sigma factors
    const std::array<const float, 256> table_sf = {
        (0.f),      (0.630f),   (1.250f),   (1.880f),   (2.510f),   (3.130f),
        (3.760f),   (4.390f),   (5.010f),   (5.640f),   (6.270f),   (6.890f),
        (7.520f),   (8.150f),   (8.770f),   (9.40f),    (10.030f),  (10.650f),
        (11.280f),  (11.910f),  (12.530f),  (13.160f),  (13.790f),  (14.410f),
        (15.040f),  (15.670f),  (16.290f),  (16.920f),  (17.550f),  (18.170f),
        (18.80f),   (19.430f),  (20.050f),  (20.680f),  (21.310f),  (21.930f),
        (22.560f),  (23.190f),  (23.810f),  (24.440f),  (25.070f),  (25.690f),
        (26.320f),  (26.950f),  (27.570f),  (28.20f),   (28.830f),  (29.450f),
        (30.080f),  (30.710f),  (31.330f),  (31.960f),  (32.590f),  (33.210f),
        (33.840f),  (34.470f),  (35.090f),  (35.720f),  (36.350f),  (36.970f),
        (37.60f),   (38.230f),  (38.850f),  (39.480f),  (40.110f),  (40.730f),
        (41.360f),  (41.990f),  (42.610f),  (43.240f),  (43.870f),  (44.490f),
        (45.120f),  (45.750f),  (46.370f),  (47.f),     (47.630f),  (48.250f),
        (48.880f),  (49.510f),  (50.130f),  (50.760f),  (51.390f),  (52.010f),
        (52.640f),  (53.270f),  (53.890f),  (54.520f),  (55.150f),  (55.770f),
        (56.40f),   (57.030f),  (57.650f),  (58.280f),  (58.910f),  (59.530f),
        (60.160f),  (60.790f),  (61.410f),  (62.040f),  (62.980f),  (64.240f),
        (65.490f),  (66.740f),  (68.f),     (69.250f),  (70.50f),   (71.760f),
        (73.010f),  (74.260f),  (75.520f),  (76.770f),  (78.020f),  (79.280f),
        (80.530f),  (81.780f),  (83.040f),  (84.290f),  (85.540f),  (86.80f),
        (88.050f),  (89.30f),   (90.560f),  (91.810f),  (93.060f),  (94.320f),
        (95.570f),  (96.820f),  (98.080f),  (99.330f),  (100.580f), (101.840f),
        (103.090f), (104.340f), (105.60f),  (106.850f), (108.10f),  (109.350f),
        (110.610f), (111.860f), (113.110f), (114.370f), (115.620f), (116.870f),
        (118.130f), (119.380f), (120.630f), (121.890f), (123.140f), (124.390f),
        (125.650f), (126.90f),  (128.150f), (129.410f), (130.660f), (131.910f),
        (133.170f), (134.420f), (135.670f), (136.930f), (138.180f), (139.430f),
        (140.690f), (141.940f), (143.190f), (144.450f), (145.70f),  (146.950f),
        (148.210f), (149.460f), (150.710f), (151.970f), (153.220f), (154.470f),
        (155.730f), (156.980f), (158.230f), (159.490f), (160.740f), (161.990f),
        (163.250f), (164.50f),  (165.750f), (167.010f), (168.260f), (169.510f),
        (170.770f), (172.020f), (173.270f), (174.530f), (175.780f), (177.030f),
        (178.290f), (179.540f), (180.790f), (182.050f), (183.30f),  (184.550f),
        (185.810f), (187.060f), (188.310f), (189.570f), (190.820f), (192.070f),
        (193.330f), (194.580f), (195.830f), (197.090f), (198.340f), (199.590f),
        (200.850f), (202.10f),  (203.350f), (204.610f), (205.860f), (207.110f),
        (208.370f), (209.620f), (210.870f), (212.130f), (213.380f), (214.630f),
        (215.890f), (217.140f), (218.390f), (219.650f), (220.90f),  (222.150f),
        (223.410f), (224.660f), (225.910f), (227.170f), (228.420f), (229.670f),
        (230.930f), (232.180f), (233.430f), (234.690f), (235.940f), (237.190f),
        (238.450f), (239.70f),  (240.950f), (242.210f), (243.460f), (244.710f),
        (245.970f), (247.220f), (248.470f), (249.730f), (250.980f), (252.230f),
        (253.490f), (254.740f), (255.990f), (255.990f) };
    State state;

    void init_sequential_bit_function(sequential_bit_t* seq_state, size_t byte_pos);
    void consume_padding_bits(sequential_bit_t* s);

    inline int get_bit_rate_code(sequential_bit_t* s);
    inline int decode_huffman_brc0(sequential_bit_t* s);
    inline int decode_huffman_brc1(sequential_bit_t* s);
    inline int decode_huffman_brc2(sequential_bit_t* s);
    inline int decode_huffman_brc3(sequential_bit_t* s);
    inline int decode_huffman_brc4(sequential_bit_t * s);

    inline bool get_sequential_bit(sequential_bit_t* seq_state) {
        auto current_byte = *(seq_state->data);
        auto res = static_cast<bool>((((current_byte) >> (((7) - (seq_state->current_bit_count)))) & (1)));
        (seq_state->current_bit_count)++;
        if ((7) < (seq_state->current_bit_count)) {
            seq_state->current_bit_count = 0;
            (seq_state->data)++;
        }
        return res;
    }
    inline int get_threshold_index(sequential_bit_t* s) {
        return ((((0x80) * (get_sequential_bit(s)))) +
            (((0x40) * (get_sequential_bit(s)))) +
            (((0x20) * (get_sequential_bit(s)))) +
            (((0x10) * (get_sequential_bit(s)))) +
            (((0x8) * (get_sequential_bit(s)))) +
            (((0x4) * (get_sequential_bit(s)))) +
            (((0x2) * (get_sequential_bit(s)))) +
            (((0x1) * (get_sequential_bit(s)))));
    }
};

