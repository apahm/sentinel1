#include <vector>
#include <iostream>
#include "ipp.h"
#include <complex>
#include <filesystem>
#include <fstream>
#include <intrin.h>
#include <bitset>
#include <algorithm>
#include <numeric>

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

struct ComplexMatrix {
    std::vector<std::vector<double>> re;
    std::vector<std::vector<double>> im;
    uint64_t rows = 0;
    uint64_t cols = 0;
};



void MeanOfRawData(ComplexMatrix &rawData, RawDataAnalysis& rawDataAnalysis) {
    double sumI = 0.0;
    double sumQ = 0.0;

    for (size_t i = 0; i < rawData.re.size(); i++) {
        sumI += static_cast<double>(std::accumulate(rawData.re[i].begin(), rawData.re[i].end(), 0));
    }

    for (size_t i = 0; i < rawData.im.size(); i++) {
        sumQ += static_cast<double>(std::accumulate(rawData.im[i].begin(), rawData.im[i].end(), 0));
    }

    rawDataAnalysis.meanOfRawDataI = sumI;
    rawDataAnalysis.meanOfRawDataQ = sumQ;

    sumI = sumI / (static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols));
    sumQ = sumQ / (static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols));
}

// Calculate the standard deviations of the raw data.
void StandardDeviationsOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);
    
    double sumIMinusMean = 0.0;
    double sumQMinusMean = 0.0;

    for (size_t i = 0; i < rawData.rows; i++) {
        for (size_t j = 0; j < rawData.cols; j++) {
            sumIMinusMean += std::pow(rawData.re[i][j] - rawDataAnalysis.meanOfRawDataI, 2);
            sumQMinusMean += std::pow(rawData.im[i][j] - rawDataAnalysis.meanOfRawDataQ, 2);
        }
    }

    rawDataAnalysis.standardDeviationsOfRawDataI = std::sqrt(k * sumIMinusMean);
    rawDataAnalysis.standardDeviationsOfRawDataQ = std::sqrt(k * sumQMinusMean);
}

// Calculate the IQ gain imbalance
void IQGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);

    rawDataAnalysis.IQGainImbalance = rawDataAnalysis.standardDeviationsOfRawDataI / rawDataAnalysis.standardDeviationsOfRawDataQ;
    rawDataAnalysis.IQGainImbalanceLowerBounds = 1 - 3 / std::sqrt(k);
    rawDataAnalysis.IQGainImbalanceUpperBounds = 1 + 3 / std::sqrt(k);
}

double getRowSumI(ComplexMatrix& rawData, uint64_t number) {
    return std::accumulate(rawData.re[number].begin(), rawData.re[number].end(), 0);
}

double getRow2SumI(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += std::pow(rawData.re[number].at(i),2);
    }
    return sum;
}

double getRowSumQ(ComplexMatrix& rawData, uint64_t number) {
    return std::accumulate(rawData.im[number].begin(), rawData.im[number].end(), 0);
}

double getRow2SumQ(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += std::pow(rawData.im[number].at(i), 2);
    }
    return sum;
}

double getRowSumMultIQ(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += rawData.re[number].at(i) * rawData.im[number].at(i);
    }
    return sum;
}

// IQ quadrature departure
void IQQuadratureDeparture(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    std::vector<double> C;

    for (size_t i = 0; i < rawData.rows; i++) {
        C.push_back((getRowSumMultIQ(rawData, i) - getRow2SumI(rawData, i) * getRow2SumI(rawData, i) / rawData.cols) /
            std::sqrt((getRow2SumI(rawData, i) - getRow2SumI(rawData, i) * getRow2SumI(rawData, i) / rawData.cols) *
                (getRow2SumQ(rawData, i) - getRow2SumQ(rawData, i) * getRow2SumQ(rawData, i) / rawData.cols)));
    }

    for (size_t i = 0; i < rawData.rows; i++)
    {
        C[i] = 0.5 * std::log((1 + C[i])/ (1 - C[i]));
    }

    double sum = std::accumulate(C.begin(), C.end(), 0) / C.size();
    
    rawDataAnalysis.IQQuadratureDeparture = std::asin(sum);
    rawDataAnalysis.IQQuadratureDepartureLowerBounds = 0;
    rawDataAnalysis.IQQuadratureDepartureUpperBounds = 0;
}

// Set the statistics significance flags
void SetStatisticsSignificanceFlags(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);
    
    rawDataAnalysis.IBiasSignificanceFlag = 0.0;
    rawDataAnalysis.IQGainSignificanceFlag = 0.0;
    rawDataAnalysis.QBiasSignificanceFlag = 0.0;
}

// Remove constant biases from the I and Q channels. This correction enforces a zero mean to the
// I and Q components of the signal.
void CorrectBiases(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.re[i].begin(), rawData.re[i].end(), rawData.re[i].begin(), 
            [rawDataAnalysis](double d) -> double { return d + rawDataAnalysis.meanOfRawDataI; });
    }

    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.im[i].begin(), rawData.im[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d + rawDataAnalysis.meanOfRawDataQ; });
    }
}

// Correct for gain imbalance between the I and Q channels.This correction
// imposes the same standard deviation for the I and Q channels.
void CorrectGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.im[i].begin(), rawData.im[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d * rawDataAnalysis.IQGainImbalance; });
    }
}

// Correct the Q channel for non - orthogonality
void CorrectQChannelForNonOrthogonality(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        for (size_t j = 0; j < rawData.cols; j++)
        {
            rawData.im[i][j] = rawData.im[i][j] / cos(rawDataAnalysis.IQQuadratureDeparture) - rawData.re[i][j] * std::tan(rawDataAnalysis.IQQuadratureDeparture);
        }
    }
}

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

int azimuthCompression(ComplexMatrix& rawData) {
    return 0;
}

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
void NominalImageReplicaGeneration(T  TXPL, T TXPSF, T TXPRR) {
    std::vector<std::complex<float>> refFunc;
    
    for (size_t i = 0; i < length; i++) {
        refFunc.emplace_back(std::exp(std::complex<float>(2 * M_PI * (((TXPSF - TXPRR * (-TXPL / 2))* i +(TXPRR/2) * std::pow(i, 2))))));
    }
}

template<typename T>
T getRangeReferenceFunction(T  TXPL, T TXPSF, T TXPRR) {
    return static_cast <typeid(T).name()>(0);
}

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

struct PositionVelocityTime {
    double XAxisPositionECEF = 0.0;
    double YAxisPositionECEF = 0.0;
    double ZAxisPositionECEF = 0.0;
    double XvelocityECEF = 0.0;
    double YvelocityECEF = 0.0;
    double ZvelocityECEF = 0.0;
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
    uint16_t PointingStatus = 0;
};

float calcRxGain(uint8_t rawRxGain) {
    return -(static_cast<float>(rawRxGain) * 0.5);
}

float calcTxPulseRampRate(uint16_t rawTxPulseRampRate) {
    rawTxPulseRampRate = _byteswap_ushort(rawTxPulseRampRate);
    float fRef = 37.53472224;
    uint8_t sign = std::bitset<16>(rawTxPulseRampRate)[15];
    uint16_t value = rawTxPulseRampRate & 0x7FFF;
    return std::pow(-1, sign) * static_cast<float>(value) * std::pow(fRef, 2) * std::pow(2, -21);
}

float calcTxPulseStartFreq(uint16_t rawTxPulseStartFreq, float TxPulseRampRate) {
    rawTxPulseStartFreq = _byteswap_ushort(rawTxPulseStartFreq);
    float fRef = 37.53472224;
    uint8_t sign = std::bitset<16>(rawTxPulseStartFreq)[15];
    uint16_t value = rawTxPulseStartFreq & 0x7FFF;
    return TxPulseRampRate /(4 * fRef) + std::pow(-1, sign) * static_cast<float>(value) * fRef * std::pow(2, -14);
}

uint32_t calcTxPulseLength(uint32_t rawTxPulseStartFreq) {
    float fRef = 37.53472224;
    return static_cast<float>(rawTxPulseStartFreq) / fRef;
}

float calcPulseRepetitionInterval(uint32_t rawPulseRepetitionInterval) {
    float fRef = 37.53472224;
    return static_cast<float>(rawPulseRepetitionInterval) / fRef;
}

float calcSamplingWindowStartTime(uint32_t rawSamplingWindowStartTime) {
    float fRef = 37.53472224;
    return static_cast<float>(rawSamplingWindowStartTime) / fRef;
}

float calcSamplingWindowLength(uint32_t rawSamplingWindowLength) {
    float fRef = 37.53472224;
    return static_cast<float>(rawSamplingWindowLength) / fRef;
}

float calcFineTime(uint16_t rawFineTime) {
    rawFineTime = _byteswap_ushort(rawFineTime);
    return (static_cast<float>(rawFineTime) + 0.5) * std::pow(2, -16.0);
}

uint32_t calcInstrumentConfigID(uint32_t rawInstrumentConfigID) {
    return _byteswap_ulong(rawInstrumentConfigID);
}

uint32_t read24Bit(std::ifstream &f) {
    uint8_t tmp8 = 0;
    uint16_t tmp16 = 0;
    f.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));
    f.read(reinterpret_cast<char*>(&tmp16), sizeof(tmp16));
    tmp16 = _byteswap_ushort(tmp16);
    return (tmp8 << 8) | tmp16;
}

/*packet_version_number	0
packet_type	0
secondary_header_flag	1
application_process_id_process_id	65
application_process_id_packet_category	12
sequence_flags	3
sequence_count	0
data_length	32381
coarse_time	1366374707
fine_time	63095
sync_marker	892270675
data_take_id	194628384
ecc_number	11
ignore_0	0
test_mode	0
rx_channel_id	1
instrument_configuration_id	7
sub_commutated_index	1
sub_commutated_data	49426
space_packet_count	0
pri_count	3775
error_flag	0
ignore_1	0
baq_mode	5
baq_block_length	31
ignore_2	0
range_decimation	1
rx_gain	12
tx_ramp_rate_polarity	1
tx_ramp_rate_magnitude	2869
tx_pulse_start_frequency_polarity	0
tx_pulse_start_frequency_magnitude	19125
tx_pulse_length	1706
ignore_3	0
rank	9
pulse_repetition_interval	20065
sampling_window_start_time	4366
sampling_window_length	9715
sab_ssb_calibration_p	0
sab_ssb_polarisation	7
sab_ssb_temp_comp	0
sab_ssb_ignore_0	0
sab_ssb_elevation_beam_address	0
sab_ssb_ignore_1	0
sab_ssb_azimuth_beam_address	0
ses_ssb_cal_mode	1
ses_ssb_ignore_0	0
ses_ssb_tx_pulse_number	0
ses_ssb_signal_type	1
ses_ssb_ignore_1	0
ses_ssb_swap	0
ses_ssb_swath_number	0
number_of_quads	12886
ignore_4	0
azi	0
baq_n	31
baqmod	5
cal_iter	0
ele_count	0
cal_mode	1
cal_p	0
cal_type	0
data_delay	4406
offset	0
packet_idx	0
pol	7
rgdec	1
rx	1
signal_type	1
swath	0
swl	9715
swst	4366
tstmod	0
txpl	45.4512
txpl_	1706
txprr	-1.92738
txprr_	-2869
txpsf	43.8013
*/

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
        // Octet 0,1
        rawData.read(reinterpret_cast<char*>(&tmp16), 2);
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

        if ((sentinelOneParam.BAQMode == BAQMode::FDBAQMode0) && (sentinelOneParam.CalType == CalType::TxCal)) {
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

struct sh_code BRC4Func(unsigned char* p, int* cposition, int* bposition) // TODO: never tested !
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
    do {
        BRCn[BRCindex] = next_bit(p, &cposition, &bposition) * 4;  // MSb=4
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 2; // then 2
        BRCn[BRCindex] += next_bit(p, &cposition, &bposition) * 1; // then 1 ...
        brc[*brcpos] = BRCn[BRCindex];
        (*brcpos)++;
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);      // smaller increment to match NQ
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
        if ((hcode_index + 128) > NQ) inc = (NQ - hcode_index);                      
        for (h = 0; h < inc; h++) {
            hcodeQE[hcode_index] = BRC(BRCn[BRCindex], p, &cposition, &bposition); 
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

int main() {
    ReadSARParam("C:/Houston/s1a-s1-raw-s-vh-20230424t123129-20230424t123155-048238-05cce5.dat");

    return 0;
}