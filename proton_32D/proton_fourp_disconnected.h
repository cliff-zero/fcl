#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{

    inline LatData contract_proton_two_point_psrc_pselsnk_function1(const PselProp &prop1, const Geometry &geo,
                                                                    const int tslice,
                                                                    const PointSelection &psel)
    // m_ts[tsep][op_src][op_snk] = trace( (\sum_x prop1(x) gms[op_src] gamma5
    // prop2(x)^\dagger gamma5) gms[op_snk] ) 0 <= tsep < total_site[3]
    {
        TIMER_VERBOSE("contract_proton_two_point_psrc_pselsnk_function");
        const Coordinate total_site = geo.total_site();
        std::vector<Complex> gc_ts(omp_get_max_threads() *
                                   total_site[3]);
        set_zero(gc_ts);
#pragma omp parallel
        {
            std::vector<Complex> c_ts(total_site[3]);
            set_zero(c_ts);
#pragma omp for
            for (long idx = 0; idx < (long)psel.size(); ++idx)
            {
                const Coordinate xg = psel[idx];
                const int tsep = mod(xg[3] - tslice, total_site[3]);
                const WilsonMatrix wm = prop1.get_elem(idx);
                c_ts[tsep] -= proton_twop_block(wm, wm, wm, 0);
                c_ts[tsep] += proton_twop_block(wm, wm, wm, 1);
            }
            for (int t = 0; t < total_site[3]; ++t)
            {
                gc_ts[omp_get_thread_num() * total_site[3] + t] = c_ts[t];
            }
        }
        std::vector<Complex> m_ts(total_site[3]);
        set_zero(m_ts);
        for (int i = 0; i < omp_get_max_threads(); ++i)
        {
            for (int t = 0; t < total_site[3]; ++t)
            {

                m_ts[t] += gc_ts[i * total_site[3] + t];
            }
        }
        LatData ld = mk_proton_two_point_table(total_site);
        set_zero(ld);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {
            Vector<Complex> m_src_snk = lat_data_complex_get(ld, make_array(tslice));
            if (tslice + tsep < total_site[3])
            {
                m_src_snk[tsep] += m_ts[tsep];
            }
            else
            {
                m_src_snk[tsep] -= m_ts[tsep];
            }
        }
        return ld;
    }

    inline void compute_proton_two_point_psrc_pselsnk1(std::vector<std::vector<Complex>> &block,
                                                       const std::string &job_tag, const int traj,
                                                       std::vector<std::vector<long>> &psel_dt_num_list)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_psrc_pselsnk");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        std::vector<int> psel_p_tslice(total_site[3]);
        set_zero(psel_p_tslice);

        LatData ld_two_point_psrc_func = mk_proton_two_point_table(total_site);
        set_zero(ld_two_point_psrc_func);
        for (long n = 0; n < n_points; n++)
        {
            const long xg_src_psel_idx = n;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int &tslice = xg_src[3];
            if (get_point_src_info(job_tag, traj, xg_src, 0).size() == 0)
            {
                continue;
            }
            psel_p_tslice[tslice] += 1;
            const PselProp &psprop = get_psel_prop_psrc(job_tag, traj, xg_src, 0, 0);
            LatData ld = contract_proton_two_point_psrc_pselsnk_function1(psprop, geo, tslice, psel);
            ld_two_point_psrc_func += ld;
        }
        const int num_dt = total_site[3];
        const int dtmin = 0;
        psel_dt_num_list.resize(num_dt);
        for (int dt = 0; dt < num_dt; ++dt)
        {
            const int del_t = dtmin + dt;
            psel_dt_num_list[dt].resize(del_t + 1);
            set_zero(psel_dt_num_list[dt]);
        }
        for (int tsrc = 0; tsrc < total_site[3]; ++tsrc)
        {
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int tsnk = mod(tsrc + del_t, total_site[3]);
                const int num_src_snk = psel_p_tslice[tsrc] * psel_p_tslice[tsnk];
                for (int dt_y = 0; dt_y <= del_t; ++dt_y)
                {
                    const int yt = mod(tsrc + dt_y, total_site[3]);
                    const int num_src_y_snk = num_src_snk * psel_p_tslice[yt];
                    psel_dt_num_list[dt][dt_y] += num_src_y_snk;
                }
            }
        }
        block.resize(total_site[3]);
        for (int tslice = 0; tslice < total_site[3]; ++tslice)
        {
            block[tslice].resize(total_site[3]);
            set_zero(block[tslice]);
            const Vector<Complex> m_src_snk = lat_data_complex_get_const(ld_two_point_psrc_func, make_array(tslice));
            for (int tsep = 0; tsep < total_site[3]; ++tsep)
            {
                block[tslice][tsep] = m_src_snk[tsep];
            }
        }
    }

    inline std::string get_proton_four_point_disc_func_noama_path(const std::string &job_tag,
                                                                  const int traj)
    {
        return ssprintf("analysis/proton_fourp_disconnected/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_four_point_disc_func_type_noama(const std::string &job_tag, const int traj,
                                                               const std::vector<int> &tsep_list)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int tsep_num = (int)tsep_list.size();
        const std::string path = get_proton_four_point_disc_func_noama_path(job_tag, traj);
        std::vector<std::string> fn_gram_list(tsep_num);

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const int &tsep = tsep_list[ti];
            fn_gram_list[ti] = path + ssprintf("/fourp-data-dt%d-disconnected.field", tsep);
        }

        bool is_complete = true;

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            if (fn_gram_list[ti] != "")
            {
                if (not is_d_field(fn_gram_list[ti]))
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
        TIMER_VERBOSE("compute_proton_four_point_disc_func_type_noama");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        std::vector<std::vector<Complex>> block;
        std::vector<std::vector<long>> psel_dt_num_list;
        compute_proton_two_point_psrc_pselsnk1(block, job_tag, traj, psel_dt_num_list);
        const std::vector<std::vector<Complex>> &twop_block = block;

        std::vector<FieldM<Complex, 8 * 8>> fourp_data(tsep_num);

        long iter = 0;
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_y_psel_idx = n;
            const Coordinate &xg_y = psel[xg_y_psel_idx];
            if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
            {
                continue;
            }
            Timer::autodisplay();
            TIMER_VERBOSE("compute_proton_four_point_disc_func_type_noama-iter");
            iter += 1;
            const ShiftShufflePlan ssp = make_shift_shuffle_plan(fsel, -xg_y);

            displayln_info(fname + ssprintf(":n=%ld iter=%ld", n,
                                            iter));
            for (int ti = 0; ti < tsep_num; ++ti)
            {
                const int &tsep = tsep_list[ti];
                contract_proton_four_point_disc_acc_noama(fourp_data[ti], job_tag, traj, xg_y, xg_y_psel_idx,
                                                          psel, fsel, ssp, tsep, twop_block, psel_dt_num_list);
            }
        }
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            const double coef = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
            fourp_data[ti] *= (coef * coef);
        }

        FieldM<Complex, 8 * 8> avg;

        for (int ti = 0; ti < tsep_num; ++ti)
        {
            qassert(is_initialized(fourp_data[ti]));
            qassert(fn_gram_list[ti] != "");
            avg = fourp_data[ti];
            write_field_float_from_double(avg, fn_gram_list[ti]);
        }
    }

    inline void compute_proton_four_point_disconnected_func_noama(const std::string &job_tag,
                                                                  const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_four_point_disconnected_func_noama");
        const std::string path = get_proton_four_point_disc_func_noama_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(path + "/fourp-data-dt2-disconnected.field"))
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
                ssprintf("lock-proton-fourp-disc-noama-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_fourp_disconnected");
        qmkdir_info(ssprintf("analysis/proton_fourp_disconnected/%s", job_tag.c_str()));
        qmkdir_info(path);
        std::vector<int> tsep_list;
        tsep_list.push_back(2);
        compute_proton_four_point_disc_func_type_noama(job_tag, traj, tsep_list);
        release_lock();
    }
}
