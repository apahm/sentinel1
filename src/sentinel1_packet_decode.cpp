#include "sentinel1_packet_decode.h"
#include <intrin.h>
#include <bitset>
#include <cassert>

Sentinel1PacketDecode::Sentinel1PacketDecode() {

}

Sentinel1PacketDecode::~Sentinel1PacketDecode() {

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
    return static_cast<float>(rawTxPulseStartFreq)  / fRef;
}

float Sentinel1PacketDecode::calcPulseRepetitionInterval(uint32_t rawPulseRepetitionInterval) {
    float fRef = 37.53472224;
    return static_cast<float>(rawPulseRepetitionInterval) * 1e-6 / fRef;
}

float Sentinel1PacketDecode::calcSamplingWindowStartTime(uint32_t rawSamplingWindowStartTime) {
    float fRef = 37.53472224;
    return static_cast<float>(rawSamplingWindowStartTime) * 1e-6 / fRef;
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

		// Octet 67
		rawData.read(reinterpret_cast<char*>(&tmp8), sizeof(tmp8));

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

		std::vector<uint8_t> data;
		std::vector<std::complex<float>> output;
		
		output.resize(2 * sentinelOneParam.NumberOfQuads);
		data.resize(sentinelOneParam.PacketDataLength - 62);

		rawData.read(reinterpret_cast<char*>(data.data()), sentinelOneParam.PacketDataLength - 62);

		state._mmap_data = data.data();

        if (sentinelOneParam.BAQMode == BAQMode::BypassMode && (static_cast<int>(sentinelOneParam.CalType) > 7)) {

        }

        if (sentinelOneParam.BAQMode == BAQMode::FDBAQMode0) {
            if (sentinelOneParam.SignalType == SignalType::Echo) {
				if (initDecodePacket(output.data(), sentinelOneParam) == 2 * sentinelOneParam.NumberOfQuads) {
					out.push_back(output);
					header.emplace_back(sentinelOneParam);
					j += 1;
				}
				if (j == 1000) {
					break;
				}
			}
        }
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

void Sentinel1PacketDecode::init_sequential_bit_function(sequential_bit_t* seq_state, size_t byte_pos) {
	seq_state->data = state._mmap_data;
	seq_state->current_bit_count = 0;
}
void Sentinel1PacketDecode::consume_padding_bits(sequential_bit_t* s) {
	auto byte_offset = static_cast<int>(((s->data) - (static_cast<uint8_t*>(state._mmap_data))));
	// make sure we are at first bit of an even byte in the next read
	if ((0) == (byte_offset % 2)) {
		// we are in an even byte
		if ((0) == (s->current_bit_count)) {
			// nothing to be done
		}
		else {
			(s->data) += (2);
			s->current_bit_count = 0;
		}
	}
	else {
		// we are in an odd byte
		(s->data) += (1);
		s->current_bit_count = 0;
	}
}
inline int Sentinel1PacketDecode::get_bit_rate_code(sequential_bit_t* s) {
	// note: evaluation order is crucial
	auto brc = ((((0x4) * (get_sequential_bit(s)))) +
		(((0x2) * (get_sequential_bit(s)))) +
		(((0x1) * (get_sequential_bit(s)))));
	if (!((((0) == (brc)) || ((1) == (brc)) || ((2) == (brc)) || ((3) == (brc)) ||
		((4) == (brc))))) {
		std::cout << (std::setw(8)) << (" brc=") << (brc) << (std::endl);
		throw std::out_of_range("brc");
	}
	return brc;
}
inline int Sentinel1PacketDecode::decode_huffman_brc0(sequential_bit_t* s) {
	if (get_sequential_bit(s)) {
		if (get_sequential_bit(s)) {
			if (get_sequential_bit(s)) {
				return 3;
			}
			else {
				return 2;
			}
		}
		else {
			return 1;
		}
	}
	else {
		return 0;
	}
}
inline int Sentinel1PacketDecode::decode_huffman_brc1(sequential_bit_t* s) {
	if (get_sequential_bit(s)) {
		if (get_sequential_bit(s)) {
			if (get_sequential_bit(s)) {
				if (get_sequential_bit(s)) {
					return 4;
				}
				else {
					return 3;
				}
			}
			else {
				return 2;
			}
		}
		else {
			return 1;
		}
	}
	else {
		return 0;
	}
}
inline int Sentinel1PacketDecode::decode_huffman_brc2(sequential_bit_t* s) {
	if (get_sequential_bit(s)) {
		if (get_sequential_bit(s)) {
			if (get_sequential_bit(s)) {
				if (get_sequential_bit(s)) {
					if (get_sequential_bit(s)) {
						if (get_sequential_bit(s)) {
							return 6;
						}
						else {
							return 5;
						}
					}
					else {
						return 4;
					}
				}
				else {
					return 3;
				}
			}
			else {
				return 2;
			}
		}
		else {
			return 1;
		}
	}
	else {
		return 0;
	}
}
inline int Sentinel1PacketDecode::decode_huffman_brc3(sequential_bit_t* s) {
	if (get_sequential_bit(s)) {
		if (get_sequential_bit(s)) {
			if (get_sequential_bit(s)) {
				if (get_sequential_bit(s)) {
					if (get_sequential_bit(s)) {
						if (get_sequential_bit(s)) {
							if (get_sequential_bit(s)) {
								if (get_sequential_bit(s)) {
									return 9;
								}
								else {
									return 8;
								}
							}
							else {
								return 7;
							}
						}
						else {
							return 6;
						}
					}
					else {
						return 5;
					}
				}
				else {
					return 4;
				}
			}
			else {
				return 3;
			}
		}
		else {
			return 2;
		}
	}
	else {
		if (get_sequential_bit(s)) {
			return 1;
		}
		else {
			return 0;
		}
	}
}
inline int Sentinel1PacketDecode::decode_huffman_brc4(sequential_bit_t* s) {
	if (get_sequential_bit(s)) {
		if (get_sequential_bit(s)) {
			if (get_sequential_bit(s)) {
				if (get_sequential_bit(s)) {
					if (get_sequential_bit(s)) {
						if (get_sequential_bit(s)) {
							if (get_sequential_bit(s)) {
								if (get_sequential_bit(s)) {
									if (get_sequential_bit(s)) {
										return 15;
									}
									else {
										return 14;
									}
								}
								else {
									if (get_sequential_bit(s)) {
										return 13;
									}
									else {
										return 12;
									}
								}
							}
							else {
								if (get_sequential_bit(s)) {
									return 11;
								}
								else {
									return 10;
								}
							}
						}
						else {
							return 9;
						}
					}
					else {
						return 8;
					}
				}
				else {
					return 7;
				}
			}
			else {
				if (get_sequential_bit(s)) {
					return 6;
				}
				else {
					return 5;
				}
			}
		}
		else {
			if (get_sequential_bit(s)) {
				return 4;
			}
			else {
				return 3;
			}
		}
	}
	else {
		if (get_sequential_bit(s)) {
			if (get_sequential_bit(s)) {
				return 2;
			}
			else {
				return 1;
			}
		}
		else {
			return 0;
		}
	}
}

int Sentinel1PacketDecode::initDecodePacket(std::complex<float>* output, Sentinel1RawPacket &sentinelOneParam) {
	std::array<uint8_t, 205> brcs;
	std::array<uint8_t, 205> thidxs;
	float fref = 37.53472f;
	
	auto number_of_quads = sentinelOneParam.NumberOfQuads;
	auto baq_block_length = sentinelOneParam.BAQBlockLength;
	auto number_of_baq_blocks = static_cast<int>(round(ceil((((((2.0f)) * (number_of_quads))) / (256)))));
	auto baq_mode = sentinelOneParam.BAQMode;
	auto swst = sentinelOneParam.SamplingWindowStartTime;
	auto delta_t_suppressed = (((3.20e+2)) / (((8) * (fref))));
	auto data_delay_us = swst + delta_t_suppressed;
	sequential_bit_t s;

	init_sequential_bit_function(&s, 0);
	auto decoded_ie_symbols = 0;
	std::array<float, MAX_NUMBER_QUADS> decoded_ie_symbols_a;
	for (auto i = 0; (i) < (MAX_NUMBER_QUADS); i += 1) {
		decoded_ie_symbols_a[i] = (0.f);
	}
	// parse ie data
	for (int block = 0; decoded_ie_symbols < number_of_quads; block++) {
		auto brc = get_bit_rate_code(&s);
		brcs[block] = brc;
		switch (brc) {
		case 0: {
			{
				// reconstruction law block=ie thidx-choice=thidx-unknown brc=0
				for (int i = 0; (((i) < (128)) && ((decoded_ie_symbols) < (number_of_quads))); i++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc0(&s);
					auto symbol_sign = 1.0f;
					if (sign_bit) {
						symbol_sign = -1.0f;
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_ie_symbols_a[decoded_ie_symbols] = v;
					(decoded_ie_symbols)++;
				}
				break;
			}
			break;
		}
		case 1: {
			{

				// reconstruction law block=ie thidx-choice=thidx-unknown brc=1
				for (int i = 0;
					(((i) < (128)) && ((decoded_ie_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc1(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_ie_symbols_a[decoded_ie_symbols] = v;
					(decoded_ie_symbols)++;
				}
				break;
			}
			break;
		}
		case 2: {
			{

				// reconstruction law block=ie thidx-choice=thidx-unknown brc=2
				for (int i = 0;
					(((i) < (128)) && ((decoded_ie_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc2(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_ie_symbols_a[decoded_ie_symbols] = v;
					(decoded_ie_symbols)++;
				}
				break;
			}
			break;
		}
		case 3: {
			{

				// reconstruction law block=ie thidx-choice=thidx-unknown brc=3
				for (int i = 0;
					(((i) < (128)) && ((decoded_ie_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc3(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_ie_symbols_a[decoded_ie_symbols] = v;
					(decoded_ie_symbols)++;
				}
				break;
			}
			break;
		}
		case 4: {
			{

				// reconstruction law block=ie thidx-choice=thidx-unknown brc=4
				for (int i = 0;
					(((i) < (128)) && ((decoded_ie_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc4(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_ie_symbols_a[decoded_ie_symbols] = v;
					(decoded_ie_symbols)++;
				}
				break;
			}
			break;
		}
		default: {
			{
				assert(0);
				break;
			}
			break;
		}
		}
	}
	consume_padding_bits(&s);
	auto decoded_io_symbols = 0;
	std::array<float, MAX_NUMBER_QUADS> decoded_io_symbols_a;
	for (auto i = 0; (i) < (MAX_NUMBER_QUADS); (i) += (1)) {
		decoded_io_symbols_a[i] = (0.f);
	}
	// parse io data
	for (int block = 0; (decoded_io_symbols) < (number_of_quads); (block)++) {
		auto brc = brcs[block];
		switch (brc) {
		case 0: {
			{

				// reconstruction law block=io thidx-choice=thidx-unknown brc=0
				for (int i = 0;
					(((i) < (128)) && ((decoded_io_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc0(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_io_symbols_a[decoded_io_symbols] = v;
					(decoded_io_symbols)++;
				}
				break;
			}
			break;
		}
		case 1: {
			{

				// reconstruction law block=io thidx-choice=thidx-unknown brc=1
				for (int i = 0;
					(((i) < (128)) && ((decoded_io_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc1(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_io_symbols_a[decoded_io_symbols] = v;
					(decoded_io_symbols)++;
				}
				break;
			}
			break;
		}
		case 2: {
			{

				// reconstruction law block=io thidx-choice=thidx-unknown brc=2
				for (int i = 0;
					(((i) < (128)) && ((decoded_io_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc2(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_io_symbols_a[decoded_io_symbols] = v;
					(decoded_io_symbols)++;
				}
				break;
			}
			break;
		}
		case 3: {
			{

				// reconstruction law block=io thidx-choice=thidx-unknown brc=3
				for (int i = 0;
					(((i) < (128)) && ((decoded_io_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc3(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_io_symbols_a[decoded_io_symbols] = v;
					(decoded_io_symbols)++;
				}
				break;
			}
			break;
		}
		case 4: {
			{

				// reconstruction law block=io thidx-choice=thidx-unknown brc=4
				for (int i = 0;
					(((i) < (128)) && ((decoded_io_symbols) < (number_of_quads)));
					(i)++) {
					auto sign_bit = get_sequential_bit(&s);
					auto mcode = decode_huffman_brc4(&s);
					auto symbol_sign = (1.0f);
					if (sign_bit) {
						symbol_sign = (-1.0f);
					}
					auto v = ((symbol_sign) * (mcode));
					// in ie and io we don't have thidx yet, will be processed later
					decoded_io_symbols_a[decoded_io_symbols] = v;
					(decoded_io_symbols)++;
				}
				break;
			}
			break;
		}
		default: {
			{
				assert(0);
				break;
			}
			break;
		}
		}
	}
	consume_padding_bits(&s);
	auto decoded_qe_symbols = 0;
	std::array<float, MAX_NUMBER_QUADS> decoded_qe_symbols_a;
	for (auto i = 0; (i) < (MAX_NUMBER_QUADS); (i) += (1)) {
		decoded_qe_symbols_a[i] = (0.f);
	}
	// parse qe data
	for (int block = 0; (decoded_qe_symbols) < (number_of_quads); (block)++) {
		auto thidx = get_threshold_index(&s);
		auto brc = brcs[block];
		thidxs[block] = thidx;
		switch (brc) {
		case 0: {
			{

				if ((thidx) <= (3)) {
					// reconstruction law block=qe thidx-choice=simple brc=0
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc0(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							if ((mcode) < (3)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (3)) {
									v = ((symbol_sign) * (table_b0.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				else {
					// reconstruction law block=qe thidx-choice=normal brc=0
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc0(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl0.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 1: {
			{

				if ((thidx) <= (3)) {
					// reconstruction law block=qe thidx-choice=simple brc=1
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc1(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							if ((mcode) < (4)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (4)) {
									v = ((symbol_sign) * (table_b1.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				else {
					// reconstruction law block=qe thidx-choice=normal brc=1
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc1(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl1.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 2: {
			{

				if ((thidx) <= (5)) {
					// reconstruction law block=qe thidx-choice=simple brc=2
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc2(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							if ((mcode) < (6)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (6)) {
									v = ((symbol_sign) * (table_b2.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				else {
					// reconstruction law block=qe thidx-choice=normal brc=2
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc2(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl2.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 3: {
			{

				if ((thidx) <= (6)) {
					// reconstruction law block=qe thidx-choice=simple brc=3
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc3(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							if ((mcode) < (9)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (9)) {
									v = ((symbol_sign) * (table_b3.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				else {
					// reconstruction law block=qe thidx-choice=normal brc=3
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc3(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl3.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 4: {
			{

				if ((thidx) <= (8)) {
					// reconstruction law block=qe thidx-choice=simple brc=4
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc4(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							if ((mcode) < (15)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (15)) {
									v = ((symbol_sign) * (table_b4.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				else {
					// reconstruction law block=qe thidx-choice=normal brc=4
					for (int i = 0;
						(((i) < (128)) && ((decoded_qe_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc4(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qe p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl4.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qe_symbols_a[decoded_qe_symbols] = v;
						(decoded_qe_symbols)++;
					}
				}
				break;
			}
			break;
		}
		default: {
			{

				assert(0);
				break;
			}
			break;
		}
		}
	}
	consume_padding_bits(&s);
	auto decoded_qo_symbols = 0;
	std::array<float, MAX_NUMBER_QUADS> decoded_qo_symbols_a;
	for (auto i = 0; (i) < (MAX_NUMBER_QUADS); (i) += (1)) {
		decoded_qo_symbols_a[i] = (0.f);
	}
	// parse qo data
	for (int block = 0; (decoded_qo_symbols) < (number_of_quads); (block)++) {
		auto brc = brcs[block];
		auto thidx = thidxs[block];
		switch (brc) {
		case 0: {
			{

				if ((thidx) <= (3)) {
					// reconstruction law block=qo thidx-choice=simple brc=0
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc0(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							if ((mcode) < (3)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (3)) {
									v = ((symbol_sign) * (table_b0.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				else {
					// reconstruction law block=qo thidx-choice=normal brc=0
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc0(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl0.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 1: {
			{

				if ((thidx) <= (3)) {
					// reconstruction law block=qo thidx-choice=simple brc=1
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc1(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							if ((mcode) < (4)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (4)) {
									v = ((symbol_sign) * (table_b1.at(thidx)));
								}
								else {
									std::setprecision(3);
									(std::cout) << (std::setw(10))
										<< (((std::chrono::high_resolution_clock::now()
											.time_since_epoch()
											.count()) -
											(state._start_time)))
										<< (" ") << (__FILE__) << (":") << (__LINE__)
										<< (" ") << (__func__) << (" ")
										<< ("mcode too large") << (" ") << (std::setw(8))
										<< (" mcode=") << (mcode) << (std::endl);
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				else {
					// reconstruction law block=qo thidx-choice=normal brc=1
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc1(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl1.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 2: {
			{

				if ((thidx) <= (5)) {
					// reconstruction law block=qo thidx-choice=simple brc=2
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc2(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							if ((mcode) < (6)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (6)) {
									v = ((symbol_sign) * (table_b2.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				else {
					// reconstruction law block=qo thidx-choice=normal brc=2
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc2(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl2.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 3: {
			{

				if ((thidx) <= (6)) {
					// reconstruction law block=qo thidx-choice=simple brc=3
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc3(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							if ((mcode) < (9)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (9)) {
									v = ((symbol_sign) * (table_b3.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				else {
					// reconstruction law block=qo thidx-choice=normal brc=3
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc3(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl3.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				break;
			}
			break;
		}
		case 4: {
			{

				if ((thidx) <= (8)) {
					// reconstruction law block=qo thidx-choice=simple brc=4
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc4(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							if ((mcode) < (15)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (15)) {
									v = ((symbol_sign) * (table_b4.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				else {
					// reconstruction law block=qo thidx-choice=normal brc=4
					for (int i = 0;
						(((i) < (128)) && ((decoded_qo_symbols) < (number_of_quads)));
						(i)++) {
						auto sign_bit = get_sequential_bit(&s);
						auto mcode = decode_huffman_brc4(&s);
						auto symbol_sign = (1.0f);
						if (sign_bit) {
							symbol_sign = (-1.0f);
						}
						// decode qo p.75
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl4.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_qo_symbols_a[decoded_qo_symbols] = v;
						(decoded_qo_symbols)++;
					}
				}
				break;
			}
			break;
		}
		default: {
			{

				assert(0);
				break;
			}
			break;
		}
		}
	}
	consume_padding_bits(&s);
	for (auto block = 0; (block) < (number_of_baq_blocks); (block) += (1)) {
		auto brc = brcs[block];
		auto thidx = thidxs[block];
		switch (brc) {
		case 0: {
			{

				// decode ie p.74 reconstruction law middle choice brc=0
				if ((thidx) <= (3)) {
					// decode ie p.74 reconstruction law simple brc=0
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (3)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (3)) {
									v = ((symbol_sign) * (table_b0.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				else {
					// decode ie p.74 reconstruction law normal brc=0
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl0.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 1: {
			{

				// decode ie p.74 reconstruction law middle choice brc=1
				if ((thidx) <= (3)) {
					// decode ie p.74 reconstruction law simple brc=1
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (4)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (4)) {
									v = ((symbol_sign) * (table_b1.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				else {
					// decode ie p.74 reconstruction law normal brc=1
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl1.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 2: {
			{

				// decode ie p.74 reconstruction law middle choice brc=2
				if ((thidx) <= (5)) {
					// decode ie p.74 reconstruction law simple brc=2
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (6)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (6)) {
									v = ((symbol_sign) * (table_b2.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				else {
					// decode ie p.74 reconstruction law normal brc=2
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl2.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 3: {
			{

				// decode ie p.74 reconstruction law middle choice brc=3
				if ((thidx) <= (6)) {
					// decode ie p.74 reconstruction law simple brc=3
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (9)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (9)) {
									v = ((symbol_sign) * (table_b3.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				else {
					// decode ie p.74 reconstruction law normal brc=3
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl3.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 4: {
			{

				// decode ie p.74 reconstruction law middle choice brc=4
				if ((thidx) <= (8)) {
					// decode ie p.74 reconstruction law simple brc=4
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (15)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (15)) {
									v = ((symbol_sign) * (table_b4.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				else {
					// decode ie p.74 reconstruction law normal brc=4
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_ie_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_ie_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode ie p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl4.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_ie_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		default: {
			{

				assert(0);
				break;
			}
			break;
		}
		}
	}
	for (auto block = 0; (block) < (number_of_baq_blocks); (block) += (1)) {
		auto brc = brcs[block];
		auto thidx = thidxs[block];
		switch (brc) {
		case 0: {
			{

				// decode io p.74 reconstruction law middle choice brc=0
				if ((thidx) <= (3)) {
					// decode io p.74 reconstruction law simple brc=0
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (3)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (3)) {
									v = ((symbol_sign) * (table_b0.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				else {
					// decode io p.74 reconstruction law normal brc=0
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl0.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 1: {
			{

				// decode io p.74 reconstruction law middle choice brc=1
				if ((thidx) <= (3)) {
					// decode io p.74 reconstruction law simple brc=1
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (4)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (4)) {
									v = ((symbol_sign) * (table_b1.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				else {
					// decode io p.74 reconstruction law normal brc=1
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl1.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 2: {
			{

				// decode io p.74 reconstruction law middle choice brc=2
				if ((thidx) <= (5)) {
					// decode io p.74 reconstruction law simple brc=2
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (6)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (6)) {
									v = ((symbol_sign) * (table_b2.at(thidx)));
								}
								else {

									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				else {
					// decode io p.74 reconstruction law normal brc=2
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl2.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {

							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 3: {
			{

				// decode io p.74 reconstruction law middle choice brc=3
				if ((thidx) <= (6)) {
					// decode io p.74 reconstruction law simple brc=3
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (9)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (9)) {
									v = ((symbol_sign) * (table_b3.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				else {
					// decode io p.74 reconstruction law normal brc=3
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl3.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		case 4: {
			{

				// decode io p.74 reconstruction law middle choice brc=4
				if ((thidx) <= (8)) {
					// decode io p.74 reconstruction law simple brc=4
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							if ((mcode) < (15)) {
								v = ((symbol_sign) * (mcode));
							}
							else {
								if ((mcode) == (15)) {
									v = ((symbol_sign) * (table_b4.at(thidx)));
								}
								else {
									assert(0);
								}
							}
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				else {
					// decode io p.74 reconstruction law normal brc=4
					for (int i = 0; (((i) < (128)) && ((((i)+(((128) * (block))))) <
						(decoded_io_symbols)));
						(i)++) {
						auto pos = ((i)+(((128) * (block))));
						auto scode = decoded_io_symbols_a[pos];
						auto mcode = static_cast<int>(fabsf(scode));
						auto symbol_sign = copysignf((1.0f), scode);
						// decode io p.74 reconstruction law right side
						auto v = (0.f);
						try {
							v = ((symbol_sign) * (table_nrl4.at(mcode)) *
								(table_sf.at(thidx)));
						}
						catch (std::out_of_range e) {
							assert(0);
						};
						decoded_io_symbols_a[pos] = v;
					}
				}
				break;
			}
			break;
		}
		default: {
			{
				assert(0);
				break;
			}
			break;
		}
		}
	}
	assert((decoded_ie_symbols) == (decoded_io_symbols));
	assert((decoded_ie_symbols) == (decoded_qe_symbols));
	assert((decoded_qo_symbols) == (decoded_qe_symbols));
	for (auto i = 0; (i) < (decoded_ie_symbols); (i) += (1)) {
		output[((2) * (i))].real(decoded_ie_symbols_a[i]);
		output[((2) * (i))].imag(decoded_qe_symbols_a[i]);
		output[((1) + (((2) * (i))))].real(decoded_io_symbols_a[i]);
		output[((1) + (((2) * (i))))].imag(decoded_qo_symbols_a[i]);
	}
	auto n = ((decoded_ie_symbols)+(decoded_io_symbols));
	return n;
}