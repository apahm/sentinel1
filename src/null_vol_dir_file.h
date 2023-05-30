#pragma once
#include <iostream>
#include <string>

struct null_vol_rec {
	uint32_t rec_seq;
	uint8_t rec_sub1;
	uint8_t rec_type;
	uint8_t rec_sub2;
	uint8_t rec_sub3;
	uint32_t length;
	std::string ascii_flag;
	std::string spare1;
	std::string format_doc;
	std::string format_ver;
	std::string format_rev;
	std::string software_id;
	std::string tape_id;
	std::string logvol_id;
	std::string phyvol_id;
	uint16_t n_phyvol;
	uint16_t first_phyvol;
	uint16_t last_phyvol;
	uint16_t curr_phyvol;
	uint32_t first_file;
	uint32_t volset_log;
	uint32_t logvol_vol;
	std::string spare2;
};