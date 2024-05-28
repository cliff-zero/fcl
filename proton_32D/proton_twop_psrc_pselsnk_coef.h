#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{ //
    inline LatData mk_proton_two_point_sum_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }
    inline void save_psel_dt_num(const std::vector<double> &psel_dt_num_list,
                                 const std::string &path)
    {
        TIMER_VERBOSE("save psel_per_tslice");
        FILE *fp = qopen(path + ".partial", "w");
        fprintf(fp, "%ld\n", (long)psel_dt_num_list.size());
        for (long i = 0; i < (long)psel_dt_num_list.size(); ++i)
        {
            const double &num = psel_dt_num_list[i];
            fprintf(fp, "%5ld %.6f\n", i, num);
        }
        qclose(fp);
        qrename(path + ".partial", path);
    }

    inline LatData contract_proton_psel_correction(const LatData &ld0, const std::vector<double> &psel_dt_num_list,
                                                   const Geometry &geo)
    {
        const Coordinate total_site = geo.total_site();
        LatData ld = mk_proton_two_point_sum_table(total_site);
        set_zero(ld);
        Vector<Complex> c_sum = lat_data_cget(ld);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {
            for (int tslice = 0; tslice < total_site[3]; ++tslice)
            {
                const int tsnk = mod(tslice + tsep, total_site[3]);
                const Complex &c_tslice_tsep = lat_data_complex_get_const(ld0, make_array(tslice, tsep))[0];
                c_sum[tsep] += c_tslice_tsep;
            }
            const double coef = -1.0 / psel_dt_num_list[tsep] * (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
            c_sum[tsep] *= coef;
        }
        return ld;
    }

    inline void psel_twop_coef_list(std::vector<double> &psel_dt_num_list, const PointSelection &psel,
                                    const Geometry &geo, const LatData &coef_table)
    {
        TIMER_VERBOSE("psel_dt_list");
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const double dx_max = sqrt(32.0 * 32.0 + 12.0 * 12.0 * 3.0);

        qassert(psel_dt_num_list.size() == total_site[3]);
        set_zero(psel_dt_num_list);

        std::vector<double> gd_ts(omp_get_max_threads() *
                                  total_site[3]);
        set_zero(gd_ts);
#pragma omp parallel
        {
            std::vector<double> d_ts(total_site[3]);
            set_zero(d_ts);
#pragma omp for
            for (long nsrc = 0; nsrc < n_points; ++nsrc)
            {
                const long xg_src_psel_idx = nsrc;
                const Coordinate &xg_src = psel[xg_src_psel_idx];
                const int tsrc = xg_src[3];
                for (long nsnk = 0; nsnk < n_points; ++nsnk)
                {
                    const long xg_snk_psel_idx = nsnk;
                    const Coordinate &xg_snk = psel[xg_snk_psel_idx];
                    const int tsnk = xg_snk[3];
                    const int del_t = mod(tsnk - tsrc, total_site[3]);
                    const Coordinate xg_sep = smod(xg_snk - xg_src, total_site);
                    const int dist = ceil(sqrt(sum(xg_sep * xg_sep)) / dx_max * 20.0);
                    const double coef = real(lat_data_complex_get_const(coef_table, make_array(dist))[0]);
                    d_ts[del_t] += coef;
                }
            }
            for (int t = 0; t < total_site[3]; ++t)
            {
                gd_ts[omp_get_thread_num() * total_site[3] + t] = d_ts[t];
            }
        }
        for (int i = 0; i < omp_get_max_threads(); ++i)
        {
            for (int t = 0; t < total_site[3]; ++t)
            {

                psel_dt_num_list[t] += gd_ts[i * total_site[3] + t];
            }
        }
        return;
    }

    inline LatData contract_proton_two_point_psrc_pselsnk_function(const PselProp &prop1, const Geometry &geo,
                                                                   const int tslice, const Coordinate &xg_src,
                                                                   const PointSelection &psel,
                                                                   const LatData &coef_twop)
    // m_ts[tsep][op_src][op_snk] = trace( (\sum_x prop1(x) gms[op_src] gamma5
    // prop2(x)^\dagger gamma5) gms[op_snk] ) 0 <= tsep < total_site[3]
    {
        TIMER_VERBOSE("contract_proton_two_point_psrc_pselsnk_function");
        const Coordinate total_site = geo.total_site();
        const double dx_max = sqrt(32.0 * 32.0 + 12.0 * 12.0 * 3.0);
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
                const Coordinate xg_sep = smod(xg - xg_src, total_site);
                const int dist = ceil(sqrt(sum(xg_sep * xg_sep)) / dx_max * 20.0);
                const double coef = real(lat_data_complex_get_const(coef_twop, make_array(dist))[0]);
                const int tsep = mod(xg[3] - tslice, total_site[3]);
                const WilsonMatrix wm = prop1.get_elem(idx);
                c_ts[tsep] += proton_twop_block(wm, wm, wm, 0) * coef;
                c_ts[tsep] -= proton_twop_block(wm, wm, wm, 1) * coef;
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

    inline std::string get_proton_two_point_psrc_pselsnk_path(const std::string &job_tag,
                                                              const int traj)
    {
        return ssprintf("analysis/proton_twop_pselsnk_coef/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_two_point_psrc_pselsnk(const std::string &job_tag,
                                                      const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_psrc_pselsnk");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        LatData ld_two_point_psrc_func;
        const std::string path = get_proton_two_point_psrc_pselsnk_path(job_tag, traj);
        const std::string path_two_point_psrc =
            path + ssprintf("/two-point-psrc-pselsnk.lat");
        const std::string path_two_point_psrc_sum =
            path + ssprintf("/two-point-psrc-pselsnk-sum.lat");
        const std::string path_psel_dt_num_list =
            path + ssprintf("/psel-dt_num_list.txt");
        const std::string path_coef_table = ssprintf("coef-list/twop/coef-twop-20slice.lat");
        LatData coef_twop;

        if (does_file_exist_sync_node(path_coef_table))
        {
            coef_twop =
                lat_data_load_info(path_coef_table);
        }
        else
        {
            qassert(false);
        }

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        std::vector<int> psel_p_tslice(total_site[3]);
        // std::vector<int> pselpertsep(total_site[3]);
        set_zero(psel_p_tslice);
        // set_zero(pselpertsep);
        std::vector<double> psel_dt_num_list(total_site[3]);
        psel_twop_coef_list(psel_dt_num_list, psel, geo, coef_twop);
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
            const PselProp &pprop = get_psel_prop_psrc(job_tag, traj, xg_src, 0, 0);
            LatData ld = contract_proton_two_point_psrc_pselsnk_function(pprop, geo, tslice, xg_src, psel, coef_twop);
            ld_two_point_psrc_func += ld;
        }
        lat_data_save_info(path_two_point_psrc, ld_two_point_psrc_func);
        LatData ld_two_point_psrc_func_sum =
            contract_proton_psel_correction(ld_two_point_psrc_func, psel_dt_num_list, geo);
        lat_data_save_info(path_two_point_psrc_sum, ld_two_point_psrc_func_sum);
        save_psel_dt_num(psel_dt_num_list, path_psel_dt_num_list);
    }

    inline void compute_proton_two_point_psrc_pselsnk_func(const std::string &job_tag, const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_psrc_pselsnk_func");
        const std::string path = get_proton_two_point_psrc_pselsnk_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(
                path + "/two-point-psrc-pselsnk.lat"))
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
                ssprintf("lock-two-point-pselsnk-func-coef-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_twop_pselsnk_coef");
        qmkdir_info(ssprintf("analysis/proton_twop_pselsnk_coef/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_proton_two_point_psrc_pselsnk(job_tag, traj);
        release_lock();
    }

} // namespace qlat
