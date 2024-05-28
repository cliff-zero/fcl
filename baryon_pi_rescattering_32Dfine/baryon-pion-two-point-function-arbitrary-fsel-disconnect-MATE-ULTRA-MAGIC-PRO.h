#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{

    inline std::string get_baryon_pi_rescattering_all_psel_disconnect_arbitrarty_func_psel_path(const std::string &job_tag,
                                                                                     const int &traj)
    {
        return ssprintf("analysis/baryon_pi_rescattering_fsel_mate_ultra_magic_pro_disconnect/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline LatData mk_baryon_pi_rescattering_fsel_for_arbitrary_mate_disconnect(
        const int &num_dtxy, const int &gram_num,
        const std::vector<std::vector<int>> &Momentum_Targets_src,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("txy", 0, num_dtxy - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_src", 0, Momentum_Targets_src.size() - 1));
        ld.info.push_back(lat_dim_number("momentum_mode_curr", 0, Momentum_Targets_curr.size() - 1));
        ld.info.push_back(lat_dim_number("gram", 0, gram_num - 1));
        ld.info.push_back(lat_dim_number("VAmu", 0, 8 - 1));
        ld.info.push_back(lat_dim_number("mu", 0, 3));
        ld.info.push_back(lat_dim_number("nu", 0, 3));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline void contract_disconnect_diagram_topology_fsel_mate_ultra(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>>> &block_meson,
        const std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> &block_baryon,
        const std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>>> &baryon_pi_two_func_coeff,
        const std::vector<std::vector<int>> &Momentum_Targets_src,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<long>> &idx_class_by_time_baryon,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<long>> &idx_class_by_time_meson,
        const std::vector<int> &fsel_num_list,
        const std::vector<int> &list_n_from_idx_fsel,
        const std::vector<std::vector<long>> &idx_class_by_time_fsel,
        const FieldSelection &fsel, const Geometry &geo,
        const int &ti, const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy, const int &t_src_snk,
        std::vector<std::vector<ComplexD>> &disconnect_diagram_result,
        const double &anti_period,
        const int &gram_num)
    {
        TIMER_VERBOSE("contract_disconnect_diagram_topology_fsel_mate_ultra");
        const Coordinate &total_site = geo.total_site();

        const int &mom_num_src = Momentum_Targets_src.size();
        const int &mom_num_curr = Momentum_Targets_curr.size();

#pragma omp parallel for collapse(3)
        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    // if (0 == fsel_num_list[xt])
                    // {
                    //     continue;
                    // }

                    for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; mom_curr_id++)
                    {
                        const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                        const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                        const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                        const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                        const Coordinate &xg_y = psel_meson[xg_y_psel_idx];

                        // long xg_x_fsel_idx = idx_class_by_time_fsel[xt][idx_x];
                        // const long index = fsel.indices[xg_x_fsel_idx];
                        // const Coordinate xg_x_l = geo.coordinate_from_index(index);
                        // const Coordinate xg_x_g = geo.coordinate_g_from_l(xg_x_l);
                        // const Coordinate &xg_x = xg_x_g;
                        // const int &t_x = xg_x[3];
                        // qassert(t_x == xt);
                        // qassert(idx_x == list_n_from_idx_fsel[xg_x_fsel_idx]);

                        if (is_baryon_smear == is_meson_smear && (xg_snk_psel_idx == xg_y_psel_idx || xg_src_psel_idx == xg_y_psel_idx))
                        {
                            continue;
                        }

                        qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                        qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                        qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);

                        const int del_t = mod(t_baryon_snk - t_baryon_src, total_site[3]);
                        const int &dt = del_t;
                        const int n_t = list_n_from_idx_baryon[xg_src_psel_idx];
                        const Baryon_Block_2PF_No_Projection_Psel_Psel &baryon_block = block_baryon[xg_snk_psel_idx][dt][n_t];
                        qassert(baryon_block.is_build);
                        qassert(baryon_block.xg_snk_psel_idx == xg_snk_psel_idx);
                        qassert(baryon_block.xg_src_psel_idx == xg_src_psel_idx);
                        qassert(baryon_block.n_snk == idx_baryon_snk);
                        qassert(baryon_block.n_src == idx_baryon_src);

                        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                        {
                            const Meson_Block_Scalar_Psel_Psel &meson_block = block_meson[xg_y_psel_idx][dt][mom_curr_id][VA_idx];

                            qassert(meson_block.is_build);
                            qassert(meson_block.xg_src_psel_idx == xg_y_psel_idx);
                            qassert(meson_block.n_src == idx_meson_ty);

                            for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
                            {
                                long phase_temp = space_dot(Momentum_Targets_src[mom_src_id], xg_src) - space_dot(Momentum_Targets_src[mom_src_id], xg_y) + space_dot(Momentum_Targets_curr[mom_curr_id], xg_snk);
                                int phase = mod(phase_temp, total_site[0]);

                                const ComplexD normal_coeff = baryon_pi_two_func_coeff[ti][t_src_snk][mom_src_id][mom_curr_id][dtxy][phase] * (ComplexD)anti_period;

                                for (int topology_id = 0; topology_id < 2; topology_id++)
                                {
                                    SpinMatrix spin_result_list = (normal_coeff * meson_block.block_pp) * baryon_block.block_pp[topology_id];
                                    // std::cout<<spin_result_list(0, 0)<<std::endl;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            if (spin_result_list(mu1, mu2) == 0.0)
                                            {
                                                continue;
                                            }
                                            disconnect_diagram_result[omp_get_thread_num()][((((mom_src_id * mom_num_curr + mom_curr_id) * gram_num + topology_id) * 8 + VA_idx) * 4 + mu1) * 4 + mu2] += spin_result_list(mu1, mu2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    inline void acc_data_to_LatData_disconnect_fsel(
        const std::vector<std::vector<ComplexD>> &disconnect_diagram_result,
        LatData &baryon_pi_data_temp,
        const std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>>> &baryon_pi_two_func_coeff,
        const std::vector<std::vector<int>> &Momentum_Targets_src,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const bool &is_baryon_smear,
        const PointsSelection &psel_baryon,
        const long &psel_baryon_num,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<std::vector<long>> &idx_class_by_time_baryon,
        const bool &is_meson_smear,
        const PointsSelection &psel_meson,
        const long &psel_meson_num,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const std::vector<std::vector<long>> &idx_class_by_time_meson,
        const std::vector<int> &fsel_num_list,
        const std::vector<int> &list_n_from_idx_fsel,
        const std::vector<std::vector<long>> &idx_class_by_time_fsel,
        const Geometry &geo, const int &ti,
        const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy, const int &t_src_snk, const int &gram_num)
    {
        TIMER_VERBOSE("acc_data_to_LatData_disconnect_fsel");

        const Coordinate &total_site = geo.total_site();
        const int &mom_num_src = Momentum_Targets_src.size();
        const int &mom_num_curr = Momentum_Targets_curr.size();

#pragma omp parallel for collapse(2)
        for (int mom_src_id = 0; mom_src_id < mom_num_src; ++mom_src_id)
        {
            for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
            {
                Vector<ComplexD> baryon_pi_data_vec_temp = lat_data_complex_get(baryon_pi_data_temp, make_array(dtxy, mom_src_id, mom_curr_id));
                for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
                {
                    for (int topology_id = 0; topology_id < gram_num; topology_id++)
                    {
                        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                        {
                            for (int mu1 = 0; mu1 < 4; mu1++)
                            {
                                for (int mu2 = 0; mu2 < 4; mu2++)
                                {
                                    baryon_pi_data_vec_temp[((topology_id * 8 + VA_idx) * 4 + mu1) * 4 + mu2] += disconnect_diagram_result[thread_id][((((mom_src_id * mom_num_curr + mom_curr_id) * gram_num + topology_id) * 8 + VA_idx) * 4 + mu1) * 4 + mu2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    inline void compute_baryon_pi_rescattering_fsel_disconnect_arbitrarty_func_psel_type(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<int>> &Momentum_Targets_src,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list,
        const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 2;

        qassert(tsep_src_list.size() == tsep_snk_list.size());
        const int tsep_num = (int)tsep_src_list.size();

        TIMER_VERBOSE("compute_baryon_pi_rescattering_fsel_disconnect_arbitrarty_func_psel_type");

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
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);
        const int &dtmax = num_dtxy + (tsep_snk_list.back() > 0 ? tsep_snk_list.back() : 0) + (tsep_src_list.back() > 0 ? tsep_src_list.back() : 0);
        const int &dtmin = +(tsep_snk_list.front() > 0 ? tsep_snk_list.front() : 0) + (tsep_src_list.front() > 0 ? tsep_src_list.front() : 0); // block 所需时间片长度

        const std::string path = get_baryon_pi_rescattering_all_psel_disconnect_arbitrarty_func_psel_path(job_tag, traj);
        std::vector<std::vector<std::string>> fn_gram_list(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];

            fn_gram_list[ti].resize(num_dtxy);
            for (int t_src_snk = 0; t_src_snk < num_dtxy; ++t_src_snk)
            {
                fn_gram_list[ti][t_src_snk] = path + ssprintf("/baryon_%s_pi_%s_rescattering_mate_ultra_magic_fsel_tsep_%d_%d_tss_%d_disconnect.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), tsep_src, tsep_snk, t_src_snk + tsep_src + tsep_snk);
            }
        }

        bool is_complete = true;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int t_src_snk = 0; t_src_snk < num_dtxy; ++t_src_snk)
            {
                if (fn_gram_list[ti][t_src_snk] != "")
                {
                    if (not is_d_field(fn_gram_list[ti][t_src_snk]))
                    {
                        is_complete = false;
                        break;
                    }
                }
            }
        }

        if (is_complete)
        {
            return;
        }

        std::vector<std::vector<LatData>> baryon_pi_data_dc_result(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_data_dc_result[ti].resize(num_dtxy);
            for (int t_src_snk = 0; t_src_snk < num_dtxy; ++t_src_snk)
            {
                baryon_pi_data_dc_result[ti][t_src_snk] = mk_baryon_pi_rescattering_fsel_for_arbitrary_mate_disconnect(t_src_snk, gram_num, Momentum_Targets_src, Momentum_Targets_curr, total_site);
            }
        }

        std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<ComplexD>>>>>> baryon_pi_two_func_coeff;
        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> fsel_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_baryon_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);
        std::vector<int> list_n_from_idx_fsel(fsel.n_elems);
        std::vector<std::vector<long>> idx_class_by_time_baryon(total_site[3]);
        std::vector<std::vector<long>> idx_class_by_time_meson(total_site[3]);
        std::vector<std::vector<long>> idx_class_by_time_fsel(total_site[3]);

        baryon_pion_two_point_arbitrary_coeffience_mate_magic_fsel(
            job_tag, traj,
            baryon_pi_two_func_coeff,
            psel_baryon_num_list,
            psel_meson_num_list,
            list_n_from_idx_baryon,
            list_n_from_idx_meson,
            idx_class_by_time_baryon,
            idx_class_by_time_meson,
            fsel_num_list,
            list_n_from_idx_fsel,
            idx_class_by_time_fsel,
            Momentum_Targets_src,
            Momentum_Targets_curr,
            fsel, geo,
            psel_baryon, psel_meson,
            tsep_src_list, tsep_snk_list,
            is_baryon_smear, is_meson_smear);

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate xg_src = psel_smear[xg_src_psel_idx];
            const PselProp &prop_src = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);
        }

        displayln_info("baryon_pion_two_point_coeffience_finish;");

        qassert(is_baryon_smear == is_meson_smear);
        std::vector<std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>>> block_meson;
        std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> block_baryon;

        contract_meson_psel_two_point_block_mate_ultra(
            block_meson,
            job_tag, traj,
            is_meson_smear,
            psel_meson,
            psel_meson_num,
            psel_meson_num_list,
            list_n_from_idx_meson,
            idx_class_by_time_meson,
            fsel_num_list,
            list_n_from_idx_fsel,
            idx_class_by_time_fsel,
            fsel,
            geo,
            total_site,
            Momentum_Targets_curr);

        contract_baryon_psel_two_point_no_projection_block(block_baryon, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

        displayln_info("contract_sequential_block_mate_coll;");
        Timer::autodisplay();

        const int y_idx_sparsity = 1;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep_src = tsep_src_list[ti];
            const int &tsep_snk = tsep_snk_list[ti];
            for (int t_src_snk = 0; t_src_snk < num_dtxy; ++t_src_snk)
            {
                long iter = 0;
                for (int t_baryon_src = 0; t_baryon_src < total_site[3]; ++t_baryon_src)
                {
                    // 获取通信域中的进程总数（节点数）
                    int MPI_SIZE;
                    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);

                    // 获取当前进程的进程号
                    int MPI_RANK;
                    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);

                    if (t_baryon_src % MPI_SIZE != MPI_RANK)
                    {
                        continue;
                    }
                    for (int dt = 0; dt < t_src_snk; ++dt)
                    {
                        const int &dtxy = dt;
                        int xt, t_meson, t_baryon_snk;
                        double anti_period;
                        three_point_time_slice(dtxy, t_src_snk, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt, anti_period);

                        Timer::autodisplay();
                        TIMER_VERBOSE("compute_baryon_pi_rescattering_all_psel_func_noama_type-iter");

                        std::vector<std::vector<ComplexD>> disconnect_diagram_result;

                        const int &mom_num_src = Momentum_Targets_src.size();
                        const int &mom_num_curr = Momentum_Targets_curr.size();

                        disconnect_diagram_result.resize(omp_get_max_threads());
                        for (int thread_id = 0; thread_id < omp_get_max_threads(); thread_id++)
                        {
                            disconnect_diagram_result[thread_id].resize(mom_num_src * mom_num_curr * gram_num * 8 * 4 * 4);
                            set_zero(disconnect_diagram_result[thread_id]);
                        }

                        contract_disconnect_diagram_topology_fsel_mate_ultra(
                            job_tag, traj,
                            block_meson,
                            block_baryon,
                            baryon_pi_two_func_coeff,
                            Momentum_Targets_src,
                            Momentum_Targets_curr,
                            is_baryon_smear,
                            psel_baryon,
                            psel_baryon_num,
                            psel_baryon_num_list,
                            list_n_from_idx_baryon,
                            idx_class_by_time_baryon,
                            is_meson_smear,
                            psel_meson,
                            psel_meson_num,
                            psel_meson_num_list,
                            list_n_from_idx_meson,
                            idx_class_by_time_meson,
                            fsel_num_list,
                            list_n_from_idx_fsel,
                            idx_class_by_time_fsel,
                            fsel, geo,
                            ti, dtmax, dtmin,
                            t_baryon_snk, xt, t_meson, t_baryon_src, dtxy, t_src_snk,
                            disconnect_diagram_result,
                            anti_period,
                            gram_num);

                        acc_data_to_LatData_disconnect_fsel(
                            disconnect_diagram_result,
                            baryon_pi_data_dc_result[ti][t_src_snk],
                            baryon_pi_two_func_coeff,
                            Momentum_Targets_src,
                            Momentum_Targets_curr,
                            is_baryon_smear,
                            psel_baryon,
                            psel_baryon_num,
                            psel_baryon_num_list,
                            list_n_from_idx_baryon,
                            idx_class_by_time_baryon,
                            is_meson_smear,
                            psel_meson,
                            psel_meson_num,
                            psel_meson_num_list,
                            list_n_from_idx_meson,
                            idx_class_by_time_meson,
                            fsel_num_list,
                            list_n_from_idx_fsel,
                            idx_class_by_time_fsel,
                            geo, ti, dtmax, dtmin,
                            t_baryon_snk, xt, t_meson, t_baryon_src, dtxy, t_src_snk, gram_num);

                        iter += 1;
                        displayln_info(fname + ssprintf(":iter=%ld ti=%ld t_src_snk=%ld tmeson=%ld dt=%ld", iter, ti, t_src_snk, t_meson, dt));
                    }
                }
                // MPI_Barrier(MPI_COMM_WORLD);
                glb_sum(baryon_pi_data_dc_result[ti][t_src_snk]);
                lat_data_save_info(fn_gram_list[ti][t_src_snk], baryon_pi_data_dc_result[ti][t_src_snk]);
                Timer::autodisplay();
            }
        }
    }

    inline void compute_baryon_pi_rescattering_fsel_disconnect_arbitrarty_func_psel(const std::string &job_tag,
                                                                         const int &traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_baryon_pi_rescattering_fsel_disconnect_arbitrarty_func_psel");
        const std::string path = get_baryon_pi_rescattering_all_psel_disconnect_arbitrarty_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_pi_%s_rescattering_mate_ultra_magic_fsel_tsep_%d_%d_tss_%d_disconnect.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), 0, 0, 1)))
        {
            return;
        }
        if (not(check_prop_smear(job_tag, traj, 0)))
        {
            displayln_info("err_check_prop_smear;");
            displayln_info(get_prop_smear_path(job_tag, traj, 0));
            displayln_info(get_psel_prop_smear_path(job_tag, traj, 0));
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-baryon-pi-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/baryon_pi_rescattering_fsel_mate_ultra_magic_pro_disconnect");
        qmkdir_info(ssprintf("analysis/baryon_pi_rescattering_fsel_mate_ultra_magic_pro_disconnect/%s", job_tag.c_str()));
        qmkdir_info(path);

        std::vector<std::vector<int>> Momentum_Targets_src;

        Momentum_Targets_src.push_back({0, 0, 0});

        Momentum_Targets_src.push_back({0, 0, 1});
        Momentum_Targets_src.push_back({0, 0, -1});
        Momentum_Targets_src.push_back({0, 1, 0});
        Momentum_Targets_src.push_back({0, -1, 0});
        Momentum_Targets_src.push_back({1, 0, 0});
        Momentum_Targets_src.push_back({-1, 0, 0});

        Momentum_Targets_src.push_back({0, 1, 1});
        Momentum_Targets_src.push_back({0, 1, -1});
        Momentum_Targets_src.push_back({0, -1, 1});
        Momentum_Targets_src.push_back({0, -1, -1});
        Momentum_Targets_src.push_back({1, 0, 1});
        Momentum_Targets_src.push_back({1, 0, -1});
        Momentum_Targets_src.push_back({1, 1, 0});
        Momentum_Targets_src.push_back({1, -1, 0});
        Momentum_Targets_src.push_back({-1, 0, 1});
        Momentum_Targets_src.push_back({-1, 0, -1});
        Momentum_Targets_src.push_back({-1, 1, 0});
        Momentum_Targets_src.push_back({-1, -1, 0});

        std::vector<std::vector<int>> Momentum_Targets_curr;
        Momentum_Targets_curr.push_back({0, 0, 0});

        Momentum_Targets_curr.push_back({0, 0, 1});
        // Momentum_Targets_curr.push_back({0, 0, -1});
        // Momentum_Targets_curr.push_back({0, 1, 0});
        // Momentum_Targets_curr.push_back({0, -1, 0});
        // Momentum_Targets_curr.push_back({1, 0, 0});
        // Momentum_Targets_curr.push_back({-1, 0, 0});

        // Momentum_Targets_curr.push_back({0, 1, 1});
        // Momentum_Targets_curr.push_back({0, 1, -1});
        // Momentum_Targets_curr.push_back({0, -1, 1});
        // Momentum_Targets_curr.push_back({0, -1, -1});
        // Momentum_Targets_curr.push_back({1, 0, 1});
        // Momentum_Targets_curr.push_back({1, 0, -1});
        // Momentum_Targets_curr.push_back({1, 1, 0});
        // Momentum_Targets_curr.push_back({1, -1, 0});
        // Momentum_Targets_curr.push_back({-1, 0, 1});
        // Momentum_Targets_curr.push_back({-1, 0, -1});
        // Momentum_Targets_curr.push_back({-1, 1, 0});
        // Momentum_Targets_curr.push_back({-1, -1, 0});

        std::vector<int> tsep_src_list;
        tsep_src_list.push_back(0);

        std::vector<int> tsep_snk_list;
        tsep_snk_list.push_back(0);

        compute_baryon_pi_rescattering_fsel_disconnect_arbitrarty_func_psel_type(
            job_tag, traj,
            Momentum_Targets_src,
            Momentum_Targets_curr,
            tsep_src_list, tsep_snk_list,
            is_baryon_smear, is_meson_smear);

        release_lock();
    }
}
