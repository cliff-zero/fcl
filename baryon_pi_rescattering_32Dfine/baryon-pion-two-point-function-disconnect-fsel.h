#pragma once

#include "data-load.h"
#include "baryon-pion-block-acc.h"
#include "two-point-function-init&correct.h"
#include "baryon-pion-coefficient.h"
#include "compute-utils.h"

namespace qlat
{

    inline std::string get_baryon_pi_two_point_disconnect_fsel_psel_func_psel_path(const std::string &job_tag,
                                                                                   const int traj)
    {
        return ssprintf("analysis/baryon_pi_two_point_disconnect_fsel/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void contract_baryon_pion_two_point_unshifted_fsel_psel_acc_x_typeI(std::vector<Vector<Complex>> &baryon_pi_data_temp, std::vector<Vector<Complex>> &meson_t2pf_data_temp, const Complex coef,
                                                                               const WilsonMatrix &wm_x_y, const Baryon_Block_Scalar_Psel_Psel &block_bar)
    {
        Complex pion_2pt_temp = pion_twop_block(wm_x_y, wm_x_y);
        meson_t2pf_data_temp[0][0] += coef * pion_2pt_temp;
        for (int gram = 0; gram < 2; gram++)
        {
            baryon_pi_data_temp[gram][0] += coef * block_bar.block_pp[gram] * pion_2pt_temp;
        }
        return;
    }

    inline void contract_baryon_pion_two_point_unshifted_fsel_psel_disconnect(std::vector<SelectedField<Complex>> &baryon_pi_disconnect_data_sfs, SelectedField<Complex> &meson_t2pf_data_sfs,
                                                                              const std::string &job_tag, const int &traj,
                                                                              const std::vector<int> &psel_baryon_num_list,
                                                                              const std::vector<int> &list_n_from_idx_baryon,
                                                                              const SelProp &prop_x_y, const Coordinate &xg_y,
                                                                              const long &xg_y_psel_idx, const int &tsep, const int &ti, const PointsSelection &psel_baryon, const PointsSelection &psel_meson, const FieldSelection &fsel,
                                                                              const std::vector<std::vector<long>> baryon_pi_two_func_coeff, const int dtmin, const std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> &block_baryon, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pi_unshifted_fsel_psel_typeI");
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const int num_dxxy = num_dtxy + 8;
        const int t_meson = xg_y[3];
        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();
        const int multiplicity = 1;
        clear(baryon_pi_disconnect_data_sfs);
        baryon_pi_disconnect_data_sfs.resize(2);

        for (int gram = 0; gram < (int)baryon_pi_disconnect_data_sfs.size(); gram++)
        {
            baryon_pi_disconnect_data_sfs[gram].init(fsel, multiplicity);
            set_zero(baryon_pi_disconnect_data_sfs[gram]);
        }

        meson_t2pf_data_sfs.init(fsel, multiplicity);
        set_zero(meson_t2pf_data_sfs);

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
            // const int dt_y = mod(t_meson - t_baryon_src, total_site[3]);
            // const int dt_meson = mod(xt - t_meson, total_site[3]);
            const int dt_baryon = mod(t_baryon_snk - t_baryon_src, total_site[3]);
            const Complex coef = 1.0 / (double)baryon_pi_two_func_coeff[ti][dt];

            // if (!(total_time_length == (tsep >= 0 ? abs(mod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(mod(xt - t_meson, total_site[3])))))
            // {
            //     displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d,total_time_length: %d,,dtxy: %d,%d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, total_time_length, dtxy, (tsep >= 0 ? std::abs(mod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(mod(xt - t_meson, total_site[3])))));
            // }
            qassert(total_time_length == (tsep >= 0 ? abs(mod(t_baryon_snk - t_baryon_src, total_site[3])) : abs(mod(xt - t_meson, total_site[3]))));

            int idx_snk = 0;
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
                    // const SelProp &prop_x_snk = get_prop(job_tag, traj, xg_snk, is_baryon_smear);
                    // const Baryon_Block_Psel_Psel &block1 = ppblock[dt][idx_snk];
                    const Baryon_Block_Scalar_Psel_Psel &block_bar = block_baryon[xg_snk_psel_idx][dt_baryon][idx_baryon_src];

                    qassert(idx_snk == list_n_from_idx_baryon[xg_snk_psel_idx]);
                    if (block_bar.xg_snk_psel_idx == -1)
                    {
                        continue;
                    }
                    qassert(block_bar.is_build);
                    // if (!(block_bar.xg_snk_psel_idx == xg_snk_psel_idx))
                    // {
                    //     displayln(ssprintf("thread_id: %d/%d, tsep: %d, t_baryon_snk: %d, xt: %d, t_meson: %d, t_baryon_src: %d, xg_snk_psel_idx: %d, block_bar.xg_snk_psel_idx: %d, list_n_from_idx_baryon[xg_snk_psel_idx]: %d, list_n_from_idx_baryon[block_bar.xg_snk_psel_idx]: %d, idx_snk:%d, dt:%d", omp_get_thread_num(), omp_get_num_threads(), tsep, t_baryon_snk, xt, t_meson, t_baryon_src, xg_snk_psel_idx, block_bar.xg_snk_psel_idx, list_n_from_idx_baryon[xg_snk_psel_idx], list_n_from_idx_baryon[block_bar.xg_snk_psel_idx], idx_snk, dt));
                    // }
                    qassert(block_bar.xg_snk_psel_idx == xg_snk_psel_idx);

#pragma omp parallel for
                    for (long idx = 0; idx < fsel.n_elems; ++idx)
                    {
                        const long xg_x_idx = idx;
                        // const Coordinate &xg_x = psel_meson[xg_x_idx];
                        const long index = fsel.indices[xg_x_idx];
                        const Coordinate xl = geo.coordinate_from_index(index);
                        const Coordinate xg = geo.coordinate_g_from_l(xl);
                        const Coordinate &xg_x = xg;
                        if (xg_x[3] != xt)
                        {
                            continue;
                        }
                        // const Coordinate xg_sep = mod(xg_x - xg_y, total_site);
                        // if ((abs(xg_sep[0]) > num_dxxy / 2) || (abs(xg_sep[1]) > num_dxxy / 2) ||
                        //     (abs(xg_sep[2]) > num_dxxy / 2))
                        // {
                        //     continue;
                        // }
                        const WilsonMatrix &wm_x_y = prop_x_y.get_elem(xg_x_idx);

                        std::vector<Vector<Complex>> baryon_pi_data_temp(2);
                        std::vector<Vector<Complex>> meson_t2pf_data_temp(1);

                        for (int gram = 0; gram < (int)baryon_pi_disconnect_data_sfs.size(); gram++)
                        {
                            baryon_pi_data_temp[gram] = baryon_pi_disconnect_data_sfs[gram].get_elems(xg_x_idx);
                        }
                        meson_t2pf_data_temp[0] = meson_t2pf_data_sfs.get_elems(xg_x_idx);
                        // Vector<Complex> baryon_pi_data_temp = lat_data_complex_get(baryon_pi_data, make_array(omp_get_thread_num(), dt));

                        if (tsep >= 0)
                        {
                            if (t_baryon_src <= t_baryon_snk)
                            {
                                contract_baryon_pion_two_point_unshifted_fsel_psel_acc_x_typeI(baryon_pi_data_temp, meson_t2pf_data_temp, coef, wm_x_y, block_bar);
                            }
                            else
                            {
                                const Complex coef1 = -1.0 * coef;
                                contract_baryon_pion_two_point_unshifted_fsel_psel_acc_x_typeI(baryon_pi_data_temp, meson_t2pf_data_temp, coef1, wm_x_y, block_bar);
                            }
                        }
                        else
                        {
                            if (mod(t_baryon_src + tsep, total_site[3]) <= mod(t_baryon_snk - tsep, total_site[3]))
                            {
                                contract_baryon_pion_two_point_unshifted_fsel_psel_acc_x_typeI(baryon_pi_data_temp, meson_t2pf_data_temp, coef, wm_x_y, block_bar);
                                // baryon_pi_data +=baryon_pi_data_temp;
                            }
                            else
                            {
                                const Complex coef1 = -1.0 * coef;
                                contract_baryon_pion_two_point_unshifted_fsel_psel_acc_x_typeI(baryon_pi_data_temp, meson_t2pf_data_temp, coef1, wm_x_y, block_bar);
                            }
                        }
                    }
                }
                idx_snk += 1;
            }
        }
        return;
    }

    inline void
    contract_baryon_pion_two_point_all_fsel_psel_acc_psel_disconnect(
        std::vector<FieldM<Complex, 1>> &baryon_pi_disconnect_data,
        FieldM<Complex, 1> &meson_t2pf_data,
        const std::string &job_tag, const int &traj,
        const std::vector<int> &psel_baryon_num_list,
        const std::vector<int> &list_n_from_idx_baryon,
        const SelProp &prop_x_y, const Coordinate &xg_y,
        const long &xg_y_psel_idx, const int &tsep, const int &ti, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp,
        const std::vector<std::vector<long>> &baryon_pi_two_func_coeff,
        const int dtmin,
        std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> block_baryon, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pion_two_point_all_fsel_psel_acc_psel_disconnect");
        const Geometry &geo = fsel.f_rank.geo();
        qassert(psel_meson[xg_y_psel_idx] == xg_y);
        qassert(baryon_pi_disconnect_data.size() == 2);
        qassert(geo == prop_x_y.geo());
        qassert(fsel.n_elems == prop_x_y.n_elems);
        qassert(is_initialized(prop_x_y));
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);

        std::vector<SelectedField<Complex>> baryon_pi_disconnect_data_sfs;
        SelectedField<Complex> meson_t2pf_data_sfs;
        contract_baryon_pion_two_point_unshifted_fsel_psel_disconnect(baryon_pi_disconnect_data_sfs, meson_t2pf_data_sfs, job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon, prop_x_y, xg_y, xg_y_psel_idx, tsep, ti, psel_baryon, psel_meson,
                                                                      fsel, baryon_pi_two_func_coeff, dtmin, block_baryon, is_baryon_smear, is_meson_smear);
        qassert(baryon_pi_disconnect_data_sfs.size() == 2);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;
        for (int gram = 0; gram < 2; gram++)
        {
            SelectedField<Complex> &s_data = baryon_pi_disconnect_data_sfs[gram];
            acc_field(baryon_pi_disconnect_data[gram], coef, s_data, ssp);
        }

        SelectedField<Complex> &meson_disconnect_s_data = meson_t2pf_data_sfs;
        acc_field(meson_t2pf_data, coef, meson_disconnect_s_data, ssp);
    }

    inline void compute_baryon_pi_two_point_disconnect_fsel_psel_func_psel_type(const std::string &job_tag, const int traj, const std::vector<int> &tsep_list, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 2;
        const int tsep_num = (int)tsep_list.size();
        const std::string path = get_baryon_pi_two_point_disconnect_fsel_psel_func_psel_path(job_tag, traj);
        std::vector<std::string> fn_disconnect_gram_list(tsep_num * gram_num);
        std::vector<std::string> fn_baryon_twop_list(tsep_num);
        std::vector<std::string> fn_meson_twop_list(tsep_num);
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep = tsep_list[ti];

            for (int gram = 0; gram < gram_num; gram++)
            {
                fn_disconnect_gram_list[gram + gram_num * ti] = path + ssprintf("/baryon_%s_pi_%s_two_point_gram_%d_dt%d_disconnect_fsel.field", (is_baryon_smear ? "smear" : "point"), (is_meson_smear ? "smear" : "point"), gram, tsep);
            }
            fn_baryon_twop_list[ti] = path + ssprintf("/baryon_%s_two_point_psel_dt%d.lat", (is_baryon_smear ? "smear" : "point"), tsep);
            fn_meson_twop_list[ti] = path + ssprintf("/meson_%s_two_point_fsel_dt%d.field", (is_meson_smear ? "smear" : "point"), tsep);
        }

        bool is_complete = true;
        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int gram = 0; gram < gram_num; gram++)
            {
                if (fn_disconnect_gram_list[gram + gram_num * ti] != "")
                {
                    if (not is_d_field(fn_disconnect_gram_list[gram + gram_num * ti]))
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
        TIMER_VERBOSE("compute_baryon_pi_two_point_disconnect_fsel_psel_func_psel_type");
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();

        const PointsSelection &psel_meson = (!is_meson_smear ? psel : psel_smear);
        const long &psel_meson_num = (!is_meson_smear ? n_points : n_points_smear);
        const PointsSelection &psel_baryon = (!is_baryon_smear ? psel : psel_smear);
        const long &psel_baryon_num = (!is_baryon_smear ? n_points : n_points_smear);

        const FieldSelection &fsel = get_field_selection(job_tag, traj);

        std::vector<std::vector<FieldM<Complex, 1>>> baryon_pi_disconnect_data(tsep_num);
        std::vector<FieldM<Complex, 1>> meson_t2pf_data(tsep_num);
        std::vector<LatData> baryon_twop_data(tsep_num);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + 2 * (tsep_list[tsep_num - 1] > 0 ? tsep_list[tsep_num - 1] : 0);
        const int dtmin = 2 * (tsep_list[0] > 0 ? tsep_list[0] : 0); // block 所需时间片长度
        const int num_dt = dtmax - dtmin + 1;

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            baryon_pi_disconnect_data[ti].resize(gram_num);
            baryon_twop_data[ti] = mk_two_point_remain_snk_table(total_site);
        }

        std::vector<std::vector<long>> baryon_pi_two_func_coeff;
        std::vector<int> psel_baryon_num_list;
        std::vector<int> psel_meson_num_list;
        std::vector<int> list_n_from_idx_baryon(psel_baryon_num);
        std::vector<int> list_n_from_idx_meson(psel_meson_num);

        std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> block_baryon;

        baryon_pion_two_point_fsel4meson_coeffience(
            job_tag, traj, baryon_pi_two_func_coeff,
            psel_baryon_num_list, psel_meson_num_list,
            list_n_from_idx_baryon, list_n_from_idx_meson,
            geo, psel_baryon, psel_meson, tsep_list, is_baryon_smear, is_meson_smear);

        displayln_info("baryon_pion_two_point_coeffience_finish;");

        contract_baryon_psel_two_point_block(block_baryon, psel_baryon_num_list, list_n_from_idx_baryon, job_tag, traj, geo, psel_baryon, total_site[3] - 1, is_baryon_smear);

        displayln_info("contract_two_point_or_pselpsel_block_finish;");

        const int y_idx_sparsity = 1;
        long iter = 0;
        // #pragma omp parallel for
        for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
        {
            load_prop_baryon_pi(job_tag, traj, tsep_list.back(), num_dtxy, t_meson, psel_meson, is_meson_smear, geo);
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
                Timer::autodisplay();
                TIMER_VERBOSE("compute_baryon_pi_two_point_disconnect_fsel_psel_func_noama_type-iter");

                iter += 1;
                const SelProp &prop_x_y = get_prop(job_tag, traj, xg_y, is_meson_smear);
                const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);

                displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter));

                for (int ti = 0; ti < tsep_num; ++ti)
                {
                    const int tsep = tsep_list[ti];

                    contract_baryon_pion_two_point_all_fsel_psel_acc_psel_disconnect(baryon_pi_disconnect_data[ti], meson_t2pf_data[ti], job_tag, traj, psel_baryon_num_list, list_n_from_idx_baryon,
                                                                                     prop_x_y, xg_y, xg_y_psel_idx, tsep, ti,
                                                                                     psel_baryon, psel_meson, fsel, ssp, baryon_pi_two_func_coeff, dtmin, block_baryon, is_baryon_smear, is_meson_smear);
                    // displayln_info(fname + ssprintf(":n=%ld iter=%ld", ny, iter) + "fsel_psel_acc_psel_typeII");
                }
            }
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            for (int gram = 0; gram < 2; gram++)
            {
                const double coef_wall = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
                const double coef = coef_wall * coef_wall;
                baryon_pi_disconnect_data[ti][gram] *= coef;
            }
        }

        FieldM<Complex, 1> avg;
        for (int gram = 0; gram < gram_num; gram++)
        {
            for (int ti = 0; ti < tsep_num; ++ti)
            {
                qassert(is_initialized(baryon_pi_disconnect_data[ti][gram]));
                qassert(fn_disconnect_gram_list[gram + gram_num * ti] != "");
                avg = baryon_pi_disconnect_data[ti][gram];
                write_field_float_from_double(avg, fn_disconnect_gram_list[gram + gram_num * ti]);
            }
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            qassert(is_initialized(meson_t2pf_data[ti]));
            qassert(fn_meson_twop_list[ti] != "");
            avg = meson_t2pf_data[ti];
            write_field_float_from_double(avg, fn_meson_twop_list[ti]);
        }

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            contract_baryon_psel_two_point_block_all_psel_acc_psel(baryon_twop_data[ti], block_baryon,
                                                                   job_tag, traj,
                                                                   psel_baryon_num,
                                                                   psel_baryon_num_list, list_n_from_idx_baryon,
                                                                   geo, psel_baryon, total_site[3] - 1, is_baryon_smear);
                                                                   
            lat_data_save_info(fn_baryon_twop_list[ti] + "_snk_remain", baryon_twop_data[ti]);

            LatData baryon_twop_data_corr = contract_two_point_function_psel_correction(baryon_twop_data[ti], psel_baryon_num_list, geo);

            lat_data_save_info(fn_baryon_twop_list[ti], baryon_twop_data_corr);
        }
    }

    inline void compute_baryon_pi_two_point_disconnect_fsel_psel_func_psel(const std::string &job_tag,
                                                                           const int traj, const bool is_baryon_smear, const bool is_meson_smear)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_baryon_pi_two_point_disconnect_fsel_psel_func_psel");
        const std::string path = get_baryon_pi_two_point_disconnect_fsel_psel_func_psel_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + ssprintf("/baryon_%s_pi_%s_two_point_gram_%d_dt%d_disconnect_fsel.field", (true ? "smear" : "point"), (false ? "smear" : "point"), 1, 0)))
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
                ssprintf("lock-baryon-pi-fsel-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/baryon_pi_two_point_disconnect_fsel");
        qmkdir_info(ssprintf("analysis/baryon_pi_two_point_disconnect_fsel/%s", job_tag.c_str()));
        qmkdir_info(path);

        std::vector<int> tsep_list;
        // tsep_list.push_back(-1);
        tsep_list.push_back(0);
        // tsep_list.push_back(1);

        compute_baryon_pi_two_point_disconnect_fsel_psel_func_psel_type(job_tag, traj, tsep_list, is_baryon_smear, is_meson_smear);

        release_lock();
    }
}
