#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "two-point-function-init&correct.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{

    inline std::string get_baryon_pi_two_point_disconnect_all_psel_func_psel_path(const std::string &job_tag,
                                                                                  const int traj)
    {
        return ssprintf("analysis/baryon_pi_two_point_disconnect/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline LatData mk_baryon_pi_two_point_disconnect(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsnk", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("txy", 0, get_baryon_pi_max_sep(total_site[3]) - 1));
        ld.info.push_back(lat_dim_number("gram", 0, 1));
        // ld.info.push_back(lat_dim_number("mu", 0, 8 - 1));
        // ld.info.push_back(lat_dim_number("nu", 0, 8 - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    // inline LatData mk_baryon_pi_two_point_disconnect_for_omp(const Coordinate &total_site)
    // {
    //     LatData ld;
    //     ld.info.push_back(lat_dim_number("omp", 0, omp_get_num_threads() - 1));
    //     ld.info.push_back(lat_dim_number("txy", 0, get_baryon_pi_max_sep(total_site[3]) - 1));
    //     ld.info.push_back(lat_dim_number("gram", 0, 1));
    //     // ld.info.push_back(lat_dim_number("mu", 0, 8 - 1));
    //     // ld.info.push_back(lat_dim_number("nu", 0, 8 - 1));
    //     ld.info.push_back(lat_dim_re_im());
    //     lat_data_alloc(ld);
    //     set_zero(ld);
    //     return ld;
    // }

    inline void
    contract_baryon_pi_two_point_disconnect_all_psel_acc_psel(
        LatData &baryon_pi_disconnect_data,
        LatData &baryon_twop_data_test,
        LatData &meson_twop_data_test,
        LatData &baryon_twop_data_test_corr,
        LatData &meson_twop_data_test_corr,
        const std::string &job_tag, const int traj,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const std::vector<int> &psel_meson_num_list,
        const std::vector<int> &list_n_from_idx_meson,
        const long &xg_y_psel_idx, const int &tsep, const int &ti,
        const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const FieldSelection &fsel,
        const std::vector<std::vector<long>> &baryon_pi_two_func_coeff,
        const int dtmin, const std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> &block_meson,
        const std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> &block_baryon, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pi_unshifted_psel_typeI");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const Coordinate xg_y = psel_meson[xg_y_psel_idx];
        const int t_meson = xg_y[3];
        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        for (int dt = 0; dt < num_dtxy; ++dt)
        {
            const int dtxy = dt - num_dtxy / 2;
            const int total_time_length = abs(dtxy) + 2 * abs(tsep);

            int xt = (tsep >= 0) ? mod(t_meson + dtxy, total_site[3]) : (dtxy >= 0 ? mod(t_meson - 2 * tsep + dtxy, total_site[3]) : mod(t_meson + 2 * tsep + dtxy, total_site[3]));
            int t_baryon_src, t_baryon_snk;
            if (dtxy >= 0)
            {
                t_baryon_src = mod(t_meson - tsep, total_site[3]);
                t_baryon_snk = mod(xt + tsep, total_site[3]);
            }
            else
            {
                t_baryon_src = mod(xt - tsep, total_site[3]);
                t_baryon_snk = mod(t_meson + tsep, total_site[3]);
            }

            qassert(total_time_length == (tsep >= 0 ? mod(t_baryon_snk - t_baryon_src, total_site[3]) : mod(xt - t_meson, total_site[3])));

            const Complex coef = 1.0 / (double)baryon_pi_two_func_coeff[ti][dt];

            const int dt_meson = mod(xt - t_meson, total_site[3]);
            const int dt_baryon = mod(t_baryon_snk - t_baryon_src, total_site[3]);
            

            for (long nsnk = 0; nsnk < psel_baryon_num; ++nsnk)
            {
                const long xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                if (xg_snk[3] != t_baryon_snk)
                {
                    continue;
                }
                for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; ++idx_baryon_src)
                {
                    const Baryon_Block_Scalar_Psel_Psel &block_bar = block_baryon[xg_snk_psel_idx][dt_baryon][idx_baryon_src];

                    qassert(idx_baryon_src == list_n_from_idx_baryon[block_bar.xg_src_psel_idx]);
                    qassert(idx_baryon_src == block_bar.n_src);
                    qassert(psel_baryon[block_bar.xg_src_psel_idx][3] == t_baryon_src);
                    qassert(psel_baryon[block_bar.xg_snk_psel_idx][3] == t_baryon_snk);
                    qassert(block_bar.xg_snk_psel_idx != -1);
                    qassert(block_bar.is_build);

                    if (!(block_bar.is_build))
                    {
                        displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d, xg_snk_psel_idx: %d, block_bar.xg_snk_psel_idx: %d, list_n_from_idx_baryon[xg_snk_psel_idx]: %d, list_n_from_idx_baryon[block_bar.xg_snk_psel_idx]: %d, dt:%d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, xg_snk_psel_idx, block_bar.xg_snk_psel_idx, list_n_from_idx_baryon[xg_snk_psel_idx], list_n_from_idx_baryon[block_bar.xg_snk_psel_idx], dt));
                    }
                    qassert(block_bar.xg_snk_psel_idx == xg_snk_psel_idx);
                    qassert(block_bar.is_build);

                    // std::vector<Complex> gc_ts(omp_get_num_threads() * 2);
                    // std::vector<Complex> gc_ts0(omp_get_num_threads());
                    // std::vector<Complex> gc_ts1(omp_get_num_threads());
                    // set_zero(gc_ts);
                    // set_zero(gc_ts0);
                    // set_zero(gc_ts1);

                        std::vector<Complex> c_ts(2);
                        std::vector<Complex> c_ts0(1);
                        std::vector<Complex> c_ts1(1);
                        std::vector<Complex> c_ts2(1);
                        std::vector<Complex> c_ts3(1);
                        set_zero(c_ts);
                        set_zero(c_ts0);
                        set_zero(c_ts1);
                        set_zero(c_ts2);
                        set_zero(c_ts3);
#pragma omp parallel for reduction(+: c_ts,c_ts0,c_ts1)
                        for (long idx = 0; idx < psel_meson_num; ++idx)
                        {
                            const long xg_x_idx = idx;
                            const Coordinate &xg_x = psel_meson[xg_x_idx];
                            // const long index = fsel.indices[xg_x_idx];
                            // const Coordinate xl = geo.coordinate_from_index(index);
                            // const Coordinate xg = geo.coordinate_g_from_l(xl);
                            // const Coordinate &xg_x = xg;
                            if (xg_x[3] != xt)
                            {
                                continue;
                            }
                            const Meson_Block_Scalar_Psel_Psel &block_m = block_meson[xg_x_idx][dt_meson][list_n_from_idx_meson[xg_y_psel_idx]];

                            qassert(block_m.n_src == list_n_from_idx_meson[xg_y_psel_idx]);
                            qassert(block_m.xg_src_psel_idx == xg_y_psel_idx);

                            if (!(block_m.is_build))
                            {
                                displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d, xg_y_psel_idx: %d, block_m.xg_snk_psel_idx: %d, dt_meson:%d, psel_meson_num_list[t_meson]:%d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, xg_y_psel_idx, block_m.xg_snk_psel_idx, dt_meson, psel_meson_num_list[t_meson]));
                            }
                            qassert(block_m.is_build);

                            // Vector<Complex> baryon_pi_data_temp = lat_data_complex_get(baryon_pi_disconnect_data, make_array(omp_get_thread_num(), dt));

                            if (tsep == 0 && is_baryon_smear == is_meson_smear)
                            {
                                if (dtxy >= 0)
                                {
                                    if ((block_bar.xg_src_psel_idx == xg_y_psel_idx) || (xg_snk_psel_idx == xg_x_idx))
                                    {
                                        continue;
                                    }
                                }
                                else
                                {
                                    if ((block_bar.xg_src_psel_idx == xg_x_idx) || (xg_snk_psel_idx == xg_y_psel_idx))
                                    {
                                        continue;
                                    }
                                }
                            }

                            c_ts[0] += coef * block_m.block_pp * block_bar.block_pp[0];
                            c_ts[1] += coef * block_m.block_pp * block_bar.block_pp[1];

                            c_ts0[0] += coef * block_bar.block_pp[0];
                            c_ts0[0] -= coef * block_bar.block_pp[1];

                            c_ts1[0] += coef * block_m.block_pp;

                            c_ts2[0] += block_bar.block_pp[0];
                            c_ts2[0] -= block_bar.block_pp[1];

                            c_ts3[0] += block_m.block_pp;
                        }

                  
                        // for (int plot_id = 0; plot_id < 2; ++plot_id)
                        // {
                        //     gc_ts[omp_get_thread_num() * 2 + plot_id] = c_ts[plot_id];
                        // }

                        // gc_ts0[omp_get_thread_num()] = c_ts0[0];
                        // gc_ts1[omp_get_thread_num()] = c_ts1[0];
                    

                    // std::vector<Complex> m_ts(2);
                    // std::vector<Complex> m_ts0(1);
                    // std::vector<Complex> m_ts1(1);
                    // set_zero(m_ts);
                    // set_zero(m_ts0);
                    // set_zero(m_ts1);

                    std::vector<Complex> &m_ts = c_ts;
                    std::vector<Complex> &m_ts0 = c_ts0;
                    std::vector<Complex> &m_ts1 = c_ts1;
                    std::vector<Complex> &m_ts2 = c_ts2;
                    std::vector<Complex> &m_ts3 = c_ts3;

                    // for (int thread_id = 0; thread_id < omp_get_num_threads(); ++thread_id)
                    // {
                    //     for (int plot_id = 0; plot_id < 2; ++plot_id)
                    //     {
                    //         m_ts[plot_id] += gc_ts[thread_id * 2 + plot_id];
                    //     }
                    //     m_ts0[0] += gc_ts0[thread_id];
                    //     m_ts1[0] += gc_ts1[thread_id];
                    // }

                    Vector<Complex> m_src_snk = lat_data_complex_get(baryon_pi_disconnect_data, make_array(t_baryon_snk, dt));
                    Vector<Complex> m_src_snk0 = lat_data_complex_get(baryon_twop_data_test, make_array(t_baryon_snk));
                    Vector<Complex> m_src_snk2 = lat_data_complex_get(baryon_twop_data_test_corr, make_array(t_baryon_snk));

                    if (tsep >= 0)
                    {
                        if (t_baryon_src <= t_baryon_snk)
                        {
                            for (int plot_id = 0; plot_id < 2; ++plot_id)
                            {
                                m_src_snk[plot_id] += m_ts[plot_id];
                            }
                            m_src_snk0[dt_baryon] += m_ts0[0];
                            m_src_snk2[dt_baryon] += m_ts2[0];
                        }
                        else
                        {
                            for (int plot_id = 0; plot_id < 2; ++plot_id)
                            {
                                m_src_snk[plot_id] -= m_ts[plot_id];
                            }
                            m_src_snk0[dt_baryon] -= m_ts0[0];
                            m_src_snk2[dt_baryon] -= m_ts2[0];
                        }
                    }
                    else
                    {
                        if (mod(t_baryon_src + tsep, total_site[3]) <= mod(t_baryon_snk - tsep, total_site[3]))
                        {
                            for (int plot_id = 0; plot_id < 2; ++plot_id)
                            {
                                m_src_snk[plot_id] += m_ts[plot_id];
                            }
                            m_src_snk0[dt_baryon] += m_ts0[0];
                            m_src_snk2[dt_baryon] += m_ts2[0];
                        }
                        else
                        {
                            for (int plot_id = 0; plot_id < 2; ++plot_id)
                            {
                                m_src_snk[plot_id] -= m_ts[plot_id];
                            }
                            m_src_snk0[dt_baryon] -= m_ts0[0];
                            m_src_snk2[dt_baryon] -= m_ts2[0];
                        }
                    }

                    Vector<Complex> m_src_snk1 = lat_data_complex_get(meson_twop_data_test, make_array(xt));
                    m_src_snk1[dt_meson] += m_ts1[0];
                    Vector<Complex> m_src_snk3 = lat_data_complex_get(meson_twop_data_test_corr, make_array(xt));
                    m_src_snk3[dt_meson] += m_ts3[0];
                }
            }
        }
    }

    inline void compute_baryon_pi_two_point_disconnect_all_psel_func_psel_type(const std::string &job_tag, const int traj, const std::vector<int> &tsep_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 2;
        const int tsep_num = (int)tsep_list.size();
        const std::string path = get_baryon_pi_two_point_disconnect_all_psel_func_psel_path(job_tag, traj);
        std::vector<std::string> fn_disconnect_gram_list(tsep_num);
        std::vector<std::string> fn_baryon_twop_list(tsep_num);
        std::vector<std::string> fn_meson_twop_list(tsep_num);
        std::vector<std::string> fn_baryon_test_list(tsep_num);
        std::vector<std::string> fn_meson_test_list(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep = tsep_list[ti];
            fn_disconnect_gram_list[ti] = path + ssprintf("/baryon_%s_pi_%s_two_point_dt%d_disconnect.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), tsep);
            fn_baryon_twop_list[ti] = path + ssprintf("/baryon_%s_two_point.field", (is_baryon_smear ? "smear" : "point"));
            fn_meson_twop_list[ti] = path + ssprintf("/pi_%s_two_point.field", (is_meson_smear ? "smear" : "point"));
            fn_baryon_test_list[ti] = path + ssprintf("/baryon_%s_two_point.test", (is_baryon_smear ? "smear" : "point"));
            fn_meson_test_list[ti] = path + ssprintf("/pi_%s_two_point.test", (is_meson_smear ? "smear" : "point"));
        }

        bool is_complete = true;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            if (fn_disconnect_gram_list[ti] != "")
            {
                if (not is_d_field(fn_disconnect_gram_list[ti]))
                {
                    is_complete = false;
                    break;
                }
            }
        }

        if (is_complete)
        {
            return;
        }
        TIMER_VERBOSE("compute_baryon_pi_two_point_disconnect_all_psel_func_psel_type");
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const PointsSelection &psel_meson = (!is_meson_smear ? psel : psel_smear);
        const long &psel_meson_num = (!is_meson_smear ? n_points : n_points_smear);
        const PointsSelection &psel_baryon = (!is_baryon_smear ? psel : psel_smear);
        const long &psel_baryon_num = (!is_baryon_smear ? n_points : n_points_smear);

        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        // std::vector<std::vector<FieldM<Complex, 8 * 8>>> baryon_pi_disconnect_data_omp(tsep_num);
        std::vector<LatData> baryon_pi_disconnect_data(tsep_num);
        std::vector<LatData> baryon_twop_data(tsep_num);
        std::vector<LatData> meson_twop_data(tsep_num);
        std::vector<LatData> baryon_twop_data_test(tsep_num);
        std::vector<LatData> meson_twop_data_test(tsep_num);
        
        std::vector<LatData> baryon_twop_data_test_corr(tsep_num);
        std::vector<LatData> meson_twop_data_test_corr(tsep_num);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + 2 * (tsep_list[tsep_num - 1] > 0 ? tsep_list[tsep_num - 1] : 0);
        const int dtmin = 2 * (tsep_list[0] > 0 ? tsep_list[0] : 0); // block 所需时间片长度
        const int num_dt = dtmax - dtmin + 1;
        qassert(dtmin == 0);

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_disconnect_data[ti] = mk_baryon_pi_two_point_disconnect(total_site);
            baryon_twop_data[ti] = mk_two_point_remain_snk_table(total_site);
            meson_twop_data[ti] = mk_two_point_remain_snk_table(total_site);
            baryon_twop_data_test[ti] = mk_two_point_remain_snk_table(total_site);
            meson_twop_data_test[ti] = mk_two_point_remain_snk_table(total_site);

            baryon_twop_data_test_corr[ti] = mk_two_point_remain_snk_table(total_site);
            meson_twop_data_test_corr[ti] = mk_two_point_remain_snk_table(total_site);
        }

        std::vector<std::vector<long>> baryon_pi_two_func_coeff;
        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_baryon_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);

        std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> block_meson;
        std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> block_baryon;

        baryon_pion_two_point_coeffience(
            job_tag, traj, baryon_pi_two_func_coeff,
            psel_baryon_num_list, psel_meson_num_list,
            list_n_from_idx_baryon, list_n_from_idx_meson,
            geo, psel_baryon, psel_meson, tsep_list, is_baryon_smear, is_meson_smear);

        displayln_info("baryon_pion_two_point_coeffience_finish;");

        contract_meson_psel_two_point_block(block_meson, psel_meson_num_list, list_n_from_idx_meson, job_tag, traj, geo, psel_meson, total_site[3] - 1, is_meson_smear);
        contract_baryon_psel_two_point_block(block_baryon, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

        displayln_info("contract_two_point_or_pselpsel_block_finish;");

        const int y_idx_sparsity = 1;
        long iter = 0;
        // #pragma omp parallel for
        for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
        {
            // load_prop_psrc_baryon_pi(job_tag, traj, tsep_list.back(), num_dtxy, t_meson, psel, geo);
            for (long ny = 0; ny < psel_meson_num; ++ny)
            {
                const long xg_y_psel_idx = ny;
                const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                if (t_meson != xg_y[3])
                {
                    continue;
                }
                if (mod((int)xg_y_psel_idx, y_idx_sparsity) != 0)
                {
                    continue;
                }
                // if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
                // {
                //     continue;
                // }
                Timer::autodisplay();
                TIMER_VERBOSE("compute_baryon_pi_two_point_disconnect_all_psel_func_noama_type-iter");

                iter += 1;
                // const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);
                // const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);

                displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter));

                for (int ti = 0; ti < tsep_num; ++ti)
                {
                    const int tsep = tsep_list[ti];

                    contract_baryon_pi_two_point_disconnect_all_psel_acc_psel(baryon_pi_disconnect_data[ti], baryon_twop_data_test[ti], meson_twop_data_test[ti],   baryon_twop_data_test_corr[ti], meson_twop_data_test_corr[ti], job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon, psel_meson_num_list, list_n_from_idx_meson, xg_y_psel_idx, tsep, ti, psel_baryon, psel_meson, fsel, baryon_pi_two_func_coeff, dtmin, block_meson, block_baryon, is_baryon_smear, is_meson_smear);
                    // displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter) + "all_psel_acc_psel_typeII");
                }
            }
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            lat_data_save_info(fn_disconnect_gram_list[ti], baryon_pi_disconnect_data[ti]);

            lat_data_save_info(fn_baryon_test_list[ti], baryon_twop_data_test[ti]);
            lat_data_save_info(fn_meson_test_list[ti], meson_twop_data_test[ti]);

            lat_data_save_info(fn_baryon_test_list[ti] + "_uncorr", baryon_twop_data_test_corr[ti]);
            lat_data_save_info(fn_meson_test_list[ti] + "_uncorr", meson_twop_data_test_corr[ti]);
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            contract_baryon_psel_two_point_block_all_psel_acc_psel(baryon_twop_data[ti], block_baryon,
                                                                   job_tag, traj,
                                                                   psel_baryon_num,
                                                                   psel_baryon_num_list, list_n_from_idx_baryon,
                                                                   geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

            contract_meson_psel_two_point_block_all_psel_acc_psel(meson_twop_data[ti], block_meson,
                                                                  job_tag, traj,
                                                                  psel_meson_num,
                                                                  psel_meson_num_list, list_n_from_idx_meson,
                                                                  geo, psel_meson,
                                                                  total_site[3] - 1, is_meson_smear);


            lat_data_save_info(fn_baryon_twop_list[ti] + "_snk_remain", baryon_twop_data[ti]);
            lat_data_save_info(fn_meson_twop_list[ti] + "_snk_remain", meson_twop_data[ti]);

            LatData baryon_twop_data_corr = contract_two_point_function_psel_correction(baryon_twop_data[ti], psel_baryon_num_list, geo);
            LatData meson_twop_data_corr = contract_two_point_function_psel_correction(meson_twop_data[ti], psel_meson_num_list, geo);

            lat_data_save_info(fn_baryon_twop_list[ti], baryon_twop_data_corr);
            lat_data_save_info(fn_meson_twop_list[ti], meson_twop_data_corr);
        }
    }

    inline void compute_baryon_pi_two_point_disconnect_all_psel_func_psel(const std::string &job_tag,
                                                                          const int traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_baryon_pi_two_point_disconnect_all_psel_func_psel");
        const std::string path = get_baryon_pi_two_point_disconnect_all_psel_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_pi_%s_two_point_dt%d_disconnect.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), 0)))
        {
            return;
        }
        if (not(check_prop_psrc(job_tag, traj, 0)))
        {
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
        qmkdir_info("analysis/baryon_pi_two_point_disconnect");
        qmkdir_info(ssprintf("analysis/baryon_pi_two_point_disconnect/%s", job_tag.c_str()));
        qmkdir_info(path);

        std::vector<int> tsep_list;
        // tsep_list.push_back(-1);
        tsep_list.push_back(0);
        // tsep_list.push_back(1);

        compute_baryon_pi_two_point_disconnect_all_psel_func_psel_type(job_tag, traj, tsep_list, is_baryon_smear, is_meson_smear);

        release_lock();
    }
}
