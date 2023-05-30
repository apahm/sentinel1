#pragma once
#include <fstream>
#include <filesystem>
#include <string>

struct sar_desc_rec {
    uint32_t rec_seq;
    uint8_t rec_sub1;
    uint8_t rec_type;
    uint8_t rec_sub2;
    uint8_t rec_sub3;
    uint32_t length;
    std::string ascii_flag;
    std::string spare1;
    std::string format_doc;
    std::string format_rev;
    std::string design_rev;
    std::string software_id;
    uint32_t file_num;
    std::string file_name;
    std::string rec_seq;
    uint64_t seq_loc;
    uint32_t seq_len;
    std::string rec_code;
    uint64_t code_loc;
    uint32_t code_len;
    std::string rec_len;
    uint64_t rlen_loc;
    uint32_t rlen_len;
    std::string spare2;
    std::string spare3;
    uint64_t n_dataset;
    uint64_t l_dataset;
    uint64_t n_map_proj;
    uint64_t l_map_proj;
    uint64_t n_plat_pos;
    uint64_t l_plat_pos;
    uint64_t n_att_data;
    uint64_t l_att_data;
    uint64_t n_radi_data;
    uint64_t l_radi_data;
    uint64_t n_radi_comp;
    uint64_t l_radi_comp;
    uint64_t n_qual_sum;
    uint64_t l_qual_sum;
    uint64_t n_data_hist;
    uint64_t l_data_hist;
    uint64_t n_rang_spec;
    uint64_t l_rang_spec;
    uint64_t n_dem_desc;
    uint64_t l_dem_desc;
    uint64_t n_radar_par;
    uint64_t l_radar_par;
    uint64_t n_anno_data;
    uint64_t l_anno_data;
    uint64_t n_det_proc;
    uint64_t l_det_proc;
    uint64_t n_cal;
    uint64_t l_cal;
    uint64_t n_gcp;
    uint64_t l_gcp;
    std::vector<uint64_t> spare4;
    uint64_t n_fac_data;
    uint64_t l_fac_data;
    std::string spare5;
};

struct dataset_sum_rec {
    rec_seq
        rec_sub1
        rec_type
        rec_sub2
        rec_sub3
        Length
        seq_num
        sar_chn
        scene_id
        scene_des
        inp_sctim
        asc-des
        pro_lat
        pro_long
        pro_head
        ellip_des
        ellip_maj
        ellip_min
        earth_mass
        grav_const
        ellip_j
        spare2
        terrain_h
        sc_lin
        sc_pix
        scene_len
        scene_wid
        spare3
        nchn
        spare5
        mission_id
        sensor_id
        orbit_num
        plat_lat
        plat_long
        plat_head
        clock_ang
        incident_ang
        spare15
        wave_length
        motion_comp
        pulse_code
        ampl_coef
        phas_coef
        chirp_ext_ind
        spare6
        fr
        rng_gate
        rng_length
        baseband_f
        rngcmp_f
        gn_polar
        gn_cross
        chn_bits
        quant_desc
        i_bias
        q_bias
        iq_ratio
        spare7
        spare8
        ele_sight
        mech_sight
        echo_track
        fa
        elev_beam
        azim_beam
        sat_bintim
        sat_clktim
        sat_clkinc
        spare9
        fac_id
        sys_id
        ver_id
        fac_code
        lev_code
        prod_type
        algor_id
        n_azilok
        n_rnglok
        bnd_azilok
        bnd_rnglok
        bnd_azi
        bnd_rng
        azi_weight
        rng_weight
        echo_track
        fa
        elev_beam
        azim_beam
        sat_bintim
        sat_clktim
        sat_clkinc
        spare9
        fac_id
        sys_id
        ver_id
        fac_code
        lev_code
        prod_type
        algor_id
        n_azilok
        n_rnglok
        bnd_azilok
        bnd_rnglok
        bnd_azi
        bnd_rng
        azi_weight
        rng_weight
        data_inpsrc
        rng_res
        azi_res
        radi_stretch
        alt_dopcen
        spare10
        crt_dopcen
        time_dir_pix
        time_dir_lin
        alt_rate
        spare12
        crt_rate
        spare13
        line_cont
        clutter_lock
        auto_focus
        line_spacing
        pix_spacing
        rngcmp_desg
        spare14
};

struct qual_sum_rec {

};

struct sdr_hist_rec {

};

struct pdr16_hist_rec {

};

struct proc_parm_rec {

};

struct map_proj_rec {

};

struct pos_data_rec {

};

struct att_data_rec {

};

struct radi_data_rec {

};

struct radi_comp_rec {

};

struct sar_ldr_file {
    sar_desc_rec _sar_desc_rec;
    dataset_sum_rec _dataset_sum_rec;
    qual_sum_rec _qual_sum_rec;
    sdr_hist_rec _sdr_hist_rec;
    pdr16_hist_rec _pdr16_hist_rec;
    proc_parm_rec _proc_parm_rec;
    map_proj_rec _map_proj_rec;
    pos_data_rec _pos_data_rec;
    att_data_rec _att_data_rec;
    radi_data_rec _radi_data_rec;
    radi_comp_rec _radi_comp_rec;
};


class CEOSFormat
{
public:
	CEOSFormat(std::filesystem::path &pathVolumeFile, 
        std::filesystem::path &pathTrailerFile,
        std::filesystem::path &pathLeaderFile,
        std::filesystem::path &pathDataSetFile,
        std::filesystem::path &pathNullVolumeFile);
    
	~CEOSFormat();
    
    int parseVolumeDirectoryFile();
    int parseSARLeaderFile();
    int parseSARDataFile();
    int parseSARTrailerFile();
    int parseNullVolumeDirectoryFile();

private:
    std::filesystem::path _pathVolumeFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.vol";
    std::filesystem::path _pathTrailerFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.tlr";
    std::filesystem::path _pathLeaderFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.ldr";
    std::filesystem::path _pathDataSetFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.raw";
    std::filesystem::path _pathNullVolumeFile = "C:/R1_33400_FN1_L0_F331/R1_33400_FN1_L0_F331.000.nul";
};

