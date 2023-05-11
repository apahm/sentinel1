#pragma once
#include <iostream>
#include <cstdint>
#include <cstring>
#include <vector>
#include <filesystem>
#include <fstream>

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

    std::vector<Sentinel1RawPacket> header;
    std::vector<PositionVelocityTime> positionVelocityTime;
    std::vector<Attitude> attitude;
};
