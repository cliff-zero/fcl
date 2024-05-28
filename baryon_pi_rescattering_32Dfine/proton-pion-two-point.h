#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "two-point-function-init&correct.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{

    inline std::string get_two_point_correlation_function_all_psel_func_psel_path(const std::string &job_tag,
                                                                                  const int &traj)
    {
        return ssprintf("analysis/two_point_correlation_function/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_two_point_correlation_function_all_psel_func_psel_type(const std::string &job_tag, const int &traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();

        const std::string path = get_two_point_correlation_function_all_psel_func_psel_path(job_tag, traj);

        std::string fn_meson_twop_list;
        std::string fn_baryon_twop_list;
        std::string fn_baryon_twop_list_no_projection;

        fn_meson_twop_list = path + ssprintf("/pi_%s_two_point_psel.lat", (is_meson_smear ? "smear" : "point"));
        fn_baryon_twop_list = path + ssprintf("/baryon_%s_two_point_psel.lat", (is_baryon_smear ? "smear" : "point"));
        fn_baryon_twop_list_no_projection = path + ssprintf("/baryon_%s_two_point_no_projection_psel.lat", (is_baryon_smear ? "smear" : "point"));

        TIMER_VERBOSE("compute_two_point_correlation_function_all_psel_func_psel_type");
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const PointsSelection &psel_meson = (!is_meson_smear ? psel : psel_smear);
        const long &psel_meson_num = (!is_meson_smear ? n_points : n_points_smear);
        const PointsSelection &psel_baryon = (!is_baryon_smear ? psel : psel_smear);
        const long &psel_baryon_num = (!is_baryon_smear ? n_points : n_points_smear);

        const FieldSelection &fsel = get_field_selection(job_tag, traj);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();

        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_baryon_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);

        std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> block_meson;
        std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> block_baryon;
        std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> block_baryon_no_projection;

        baryon_pion_two_point_simple(
            job_tag, traj,
            psel_baryon_num_list,
            psel_meson_num_list,
            list_n_from_idx_baryon,
            list_n_from_idx_meson,
            geo, psel_baryon, psel_meson,
            is_baryon_smear, is_meson_smear);

        contract_meson_psel_two_point_block(block_meson, psel_meson_num_list, list_n_from_idx_meson, job_tag, traj, geo, psel_meson, total_site[3] - 1, is_meson_smear);
        contract_baryon_psel_two_point_block(block_baryon, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);
        contract_baryon_psel_two_point_no_projection_block(block_baryon_no_projection, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

        displayln_info("contract_two_point_or_pselpsel_block_finish;");

        
        std::vector<Complex> meson_two_raw(total_site[3]);
        std::vector<Complex> baryon_two_raw(total_site[3]);
        std::vector<SpinMatrix> baryon_two_no_projection_raw(total_site[3]);

        contract_meson_psel_two_point_block_all_psel_gys_omp_optimiztion(
        meson_two_raw,
        block_meson,
        job_tag, traj,
        psel_meson,
        psel_meson_num,
        psel_meson_num_list,
        list_n_from_idx_meson,
        geo,
        is_meson_smear);

        contract_baryon_psel_two_point_block_all_psel_gys_omp_optimiztion(
        baryon_two_raw,
        block_baryon,
        job_tag, traj,
        psel_baryon,
        psel_baryon_num,
        psel_baryon_num_list,
        list_n_from_idx_baryon,
        geo,
        is_baryon_smear);

        contract_baryon_psel_two_point_block_no_projection_all_psel_gys_omp_optimiztion(
        baryon_two_no_projection_raw,
        block_baryon_no_projection,
        job_tag, traj,
        psel_baryon,
        psel_baryon_num,
        psel_baryon_num_list,
        list_n_from_idx_baryon,
        geo,
        is_baryon_smear);

        LatData meson_twop_data = mk_two_point_sum_table(total_site);
        LatData baryon_twop_data = mk_two_point_sum_table(total_site);
        LatData baryon_twop_data_no_projection = mk_two_point_sum_table_no_projection(total_site);

        contract_two_point_function_psel_correction(
        meson_twop_data,
        baryon_twop_data,
        baryon_twop_data_no_projection,
        meson_two_raw,
        baryon_two_raw,
        baryon_two_no_projection_raw,
        psel_meson_num_list,
        psel_baryon_num_list,
        total_site);

        lat_data_save_info(fn_meson_twop_list, meson_twop_data);
        lat_data_save_info(fn_baryon_twop_list, baryon_twop_data);
        lat_data_save_info(fn_baryon_twop_list_no_projection, baryon_twop_data_no_projection);
    }

    inline void compute_two_point_correlation_function_all_psel_func_psel(const std::string &job_tag,
                                                                          const int &traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_two_point_correlation_function_all_psel_func_psel");
        const std::string path = get_two_point_correlation_function_all_psel_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_two_point_no_projection_psel.lat", (is_baryon_smear ? "smear" : "point"))))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_two_pointn_psel.lat", (is_baryon_smear ? "smear" : "point"))))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/meson_%s_two_pointn_psel.lat", (is_meson_smear ? "smear" : "point"))))
        {
            return;
        }
        if (not(check_prop_smear_psel(job_tag, traj, 0)))
        {
            displayln_info("err_check_prop_smear;");
            displayln_info(get_prop_smear_path(job_tag, traj, 0));
            displayln_info(get_psel_prop_smear_path(job_tag, traj, 0));
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-baryon-pi-psel-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/two_point_correlation_function");
        qmkdir_info(ssprintf("analysis/two_point_correlation_function/%s", job_tag.c_str()));
        qmkdir_info(path);

        std::vector<int> tsep_list;
        // tsep_list.push_back(-1);
        tsep_list.push_back(0);
        // tsep_list.push_back(1);

        compute_two_point_correlation_function_all_psel_func_psel_type(job_tag, traj, is_baryon_smear, is_meson_smear);

        release_lock();
    }
}