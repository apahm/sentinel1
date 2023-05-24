#pragma once
#include <iosfwd>
#include <string>
// Appendix B-1: Volume Descriptor Record Contents
struct vol_desc_rec {
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
    std::string phyvol_id;
    std::string logvol_id;
    std::string volset_id;
    uint16_t phyvol_cnt;
    uint16_t first_phyvol;
    uint16_t last_phyvol;
    uint16_t curr_phyvol;
    uint32_t first_file;
    uint32_t volset_log;
    uint32_t phyvol_log;
    std::string logvol_date;
    std::string logvol_time;
    std::string logvol_country;
    std::string logvol_agency;
    std::string logvol_facility;
    uint32_t n_filepoint;
    uint32_t n_voldir;
    std::string spare2;
    std::string product_id;
    std::string spare3;
};

// Appendix B-2
struct vol_ptr_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    std::string ascii_flag;
    std::string spare1;
    uint32_t file_num;
    std::string file_name;
    std::string file_class;
    std::string file_code;
    std::string data_type;
    std::string data_code;
    uint64_t nrec;
    uint64_t first_len;
    uint64_t max_len;
    std::string len_type;
    std::string len_code;
    uint16_t first_phyvol;
    uint16_t last_phyvol;
    uint64_t first_rec;
    uint64_t last_rec;
    std::string spare2;
    std::string spare3;
};

struct vol_file_pntr_rec {
    vol_ptr_rec _ldr_file_ptr_rec;
    vol_ptr_rec _img_opt_ptr_rec;
    vol_ptr_rec _trl_file_ptr_rec;
};

// Appendix B-5
struct vol_text_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    std::string ascii_flag;
    std::string cont_flag;
    std::string product_type;
    std::string product_create;
    std::string phyvol_id;
    std::string scene_id;
    std::string scene_loc;
    std::string copyright_info;
    std::string spare2;
};

struct vol_dir_file {
    vol_desc_rec _vol_desc_rec;
    vol_file_pntr_rec _file_pntr_rec;
    vol_text_rec _text_rec;
};