#include "sentinel1.h"
#include "ipp.h"

Sentinel::Sentinel()
{
    sentinel1PacketDecode.ReadSARParam("C:/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE.SAFE/s1a-s3-raw-s-hh-20220710t213600-20220710t213625-044043-0541db.dat");
    calcParams();
}

Sentinel::~Sentinel()
{

}

void Sentinel::MeanOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
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

float Sentinel::getNormSateliteVelocity(Sentinel1PacketDecode& sentinel1PacketDecode, uint64_t numberOfNavigPacket) {
    return std::sqrt(std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).XvelocityECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).YvelocityECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).ZvelocityECEF, 2));
}

float Sentinel::getNormGroundVelocity(Sentinel1PacketDecode& sentinel1PacketDecode) {
    float sateliteVelocity = getNormSateliteVelocity(sentinel1PacketDecode, 0);
    float satelitePosition = getNormSatelitePosition(sentinel1PacketDecode, 0);

    float groundVelocity = 0.0;
    
    double lat = std::atan(sentinel1PacketDecode.positionVelocityTime.at(0).ZAxisPositionECEF /
        sentinel1PacketDecode.positionVelocityTime.at(0).XAxisPositionECEF);

    return groundVelocity;
}

float Sentinel::getEffectiveVelocity(Sentinel1PacketDecode& sentinel1PacketDecode) {
    return std::sqrt(getNormGroundVelocity(sentinel1PacketDecode) * getNormSateliteVelocity(sentinel1PacketDecode, 0));
}

float Sentinel::getNormSatelitePosition(Sentinel1PacketDecode& sentinel1PacketDecode, uint64_t numberOfNavigPacket) {
    return std::sqrt(std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).XAxisPositionECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).YAxisPositionECEF, 2) +
        std::pow(sentinel1PacketDecode.positionVelocityTime.at(numberOfNavigPacket).ZAxisPositionECEF, 2));
}

int Sentinel::rangeCompression(ComplexMatrix& rawData, bool calcIfft) {
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

void Sentinel::StandardDeviationsOfRawData(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
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

void Sentinel::IQGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);

    rawDataAnalysis.IQGainImbalance = rawDataAnalysis.standardDeviationsOfRawDataI / rawDataAnalysis.standardDeviationsOfRawDataQ;
    rawDataAnalysis.IQGainImbalanceLowerBounds = 1 - 3 / std::sqrt(k);
    rawDataAnalysis.IQGainImbalanceUpperBounds = 1 + 3 / std::sqrt(k);
}
double Sentinel::getRowSumI(ComplexMatrix& rawData, uint64_t number) {
    return std::accumulate(rawData.re[number].begin(), rawData.re[number].end(), 0);
}

double Sentinel::getRow2SumI(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += std::pow(rawData.re[number].at(i), 2);
    }
    return sum;
}

double Sentinel::getRowSumQ(ComplexMatrix& rawData, uint64_t number) {
    return std::accumulate(rawData.im[number].begin(), rawData.im[number].end(), 0);
}

double Sentinel::getRow2SumQ(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += std::pow(rawData.im[number].at(i), 2);
    }
    return sum;
}

double Sentinel::getRowSumMultIQ(ComplexMatrix& rawData, uint64_t number) {
    double sum = 0.0;
    for (size_t i = 0; i < rawData.cols; i++) {
        sum += rawData.re[number].at(i) * rawData.im[number].at(i);
    }
    return sum;
}

void Sentinel::IQQuadratureDeparture(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    std::vector<double> C;

    for (size_t i = 0; i < rawData.rows; i++) {
        C.push_back((getRowSumMultIQ(rawData, i) - getRow2SumI(rawData, i) * getRow2SumI(rawData, i) / rawData.cols) /
            std::sqrt((getRow2SumI(rawData, i) - getRow2SumI(rawData, i) * getRow2SumI(rawData, i) / rawData.cols) *
                (getRow2SumQ(rawData, i) - getRow2SumQ(rawData, i) * getRow2SumQ(rawData, i) / rawData.cols)));
    }

    for (size_t i = 0; i < rawData.rows; i++)
    {
        C[i] = 0.5 * std::log((1 + C[i]) / (1 - C[i]));
    }

    double sum = std::accumulate(C.begin(), C.end(), 0) / C.size();

    rawDataAnalysis.IQQuadratureDeparture = std::asin(sum);
    rawDataAnalysis.IQQuadratureDepartureLowerBounds = 0;
    rawDataAnalysis.IQQuadratureDepartureUpperBounds = 0;
}

void Sentinel::SetStatisticsSignificanceFlags(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    double k = 1 / static_cast<double>(rawData.rows) * static_cast<double>(rawData.cols);

    rawDataAnalysis.IBiasSignificanceFlag = 0.0;
    rawDataAnalysis.IQGainSignificanceFlag = 0.0;
    rawDataAnalysis.QBiasSignificanceFlag = 0.0;
}

void Sentinel::CorrectBiases(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.re[i].begin(), rawData.re[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d + rawDataAnalysis.meanOfRawDataI; });
    }

    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.im[i].begin(), rawData.im[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d + rawDataAnalysis.meanOfRawDataQ; });
    }
}

void Sentinel::CorrectGainImbalance(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        std::transform(rawData.im[i].begin(), rawData.im[i].end(), rawData.re[i].begin(),
            [rawDataAnalysis](double d) -> double { return d * rawDataAnalysis.IQGainImbalance; });
    }
}

void Sentinel::CorrectQChannelForNonOrthogonality(ComplexMatrix& rawData, RawDataAnalysis& rawDataAnalysis) {
    for (size_t i = 0; i < rawData.rows; i++) {
        for (size_t j = 0; j < rawData.cols; j++)
        {
            rawData.im[i][j] = rawData.im[i][j] / cos(rawDataAnalysis.IQQuadratureDeparture) - rawData.re[i][j] * std::tan(rawDataAnalysis.IQQuadratureDeparture);
        }
    }
}

bool Sentinel::checkRefFunc() {
    double TXPRR = sentinel1PacketDecode.header[0].TxRampRate;
    double TXPSF = sentinel1PacketDecode.header[0].TxPulseStartFreq;
    double TXPL = sentinel1PacketDecode.header[0].TxPulseLength;

    bool f = false;

    for (size_t i = 0; i < sentinel1PacketDecode.out.size(); i++) {
        if (TXPRR != sentinel1PacketDecode.header[i].TxRampRate) {
            f = true;
        }
        if (TXPSF != sentinel1PacketDecode.header[i].TxPulseStartFreq) {
            f = true;
        }
        if (TXPL != sentinel1PacketDecode.header[i].TxPulseLength) {
            f = true;
        }
    }

    return f;
}

template<typename T>
std::complex<double> Sentinel::NominalImageReplicaGeneration(T  TXPL, T TXPSF, T TXPRR, double i) {
    return std::exp(std::complex<double>(2 * M_PI * (((TXPSF - TXPRR * (-TXPL / 2)) * i + (TXPRR / 2) * std::pow(i, 2)))));
}

template<typename T>
T Sentinel::MigrationFactor(T carrierFrequency, T freqAzimuth, T effectiveRadarVelocity) {
    return std::sqrt(1 - (std::pow(speedOfLight / carrierFrequency, 2) * std::pow(freqAzimuth, 2) / 4 * std::pow(effectiveRadarVelocity, 2)));
}

template<typename T>
T Sentinel::AzimuthFMRate(T carrierFrequency, T frequencyDopplerCentroid, T effectiveRadarVelocity, T slantRange) {
    return 2 * std::pow(effectiveRadarVelocity, 2) * MigrationFactor<T>(carrierFrequency, frequencyDopplerCentroid, effectiveRadarVelocity) / ((speedOfLight / carrierFrequency) * slantRange);
}

int Sentinel::calcParams() {
    rangeSampleFreq = RangeDecimationToSampleRate(sentinel1PacketDecode.header[0].RangeDecimation);
    rangeSamplePeriod = 1.0 / rangeSampleFreq;
    rangeStartTime = sentinel1PacketDecode.header[0].SamplingWindowStartTime + 320.0 * 1e-6 / (8 * 37.53472224);
    azSampleFreq = 1 / sentinel1PacketDecode.header[0].PulseRepetitionInterval;
    azSamplePeriod = sentinel1PacketDecode.header[0].PulseRepetitionInterval;

    for (size_t i = 0; i < sentinel1PacketDecode.out[0].size(); i++) {
        fastTime.push_back(rangeStartTime + i * rangeSamplePeriod);
    }

    for (size_t i = 0; i < sentinel1PacketDecode.out[0].size(); i++) {
        slantRange.push_back((sentinel1PacketDecode.header[0].Rank * sentinel1PacketDecode.header[0].PulseRepetitionInterval + fastTime.at(i)) * speedOfLight / 2);
    }
    float SWL = sentinel1PacketDecode.header[0].NumberOfQuads * 2 / rangeSampleFreq;

    for (size_t i = 0; i < sentinel1PacketDecode.out.size(); i++) {
        azimuthFreq.push_back((-azSampleFreq / 2) + i / (sentinel1PacketDecode.header[0].PulseRepetitionInterval * sentinel1PacketDecode.out.size()));
    }

    double TXPRR = 1e12 * sentinel1PacketDecode.header[0].TxRampRate;
    double TXPSF = 1e6 * sentinel1PacketDecode.header[0].TxPulseStartFreq;
    double TXPL = 1e-6 * sentinel1PacketDecode.header[0].TxPulseLength;

    for (size_t i = 0; i < static_cast<int>(TXPL * rangeSampleFreq); i++) {
        double samplingT = (double)i / (double)rangeSampleFreq - (double)TXPL / 2.0;
        double sum1 = (TXPSF + TXPRR * (TXPL / 2)) * samplingT;
        double sum2 = TXPRR * std::pow(samplingT, 2) / 2;

        double arg = 2 * M_PI * (sum1 + sum2);
        double re = cos(arg);
        double im = sin(arg);

        refFunc.push_back(std::complex<double>(re, im));
    }



    return 0;
}