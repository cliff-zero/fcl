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
    inline void save_psel_per_tslice(const std::vector<int> &psel_p_tslice,
                                     const std::string &path)
    {
        TIMER_VERBOSE("save psel_per_tslice");
        FILE *fp = qopen(path + ".partial", "w");
        fprintf(fp, "%ld\n", (long)psel_p_tslice.size());
        for (long i = 0; i < (long)psel_p_tslice.size(); ++i)
        {
            const int &num = psel_p_tslice[i];
            fprintf(fp, "%5ld %3d\n", i, num);
        }
        qclose(fp);
        qrename(path + ".partial", path);
    }

    inline LatData contract_proton_psel_correction(const LatData &ld0, const std::vector<int> &psel_p_tslice,
                                                   const Geometry &geo)
    {
        const Coordinate total_site = geo.total_site();
        LatData ld = mk_proton_two_point_sum_table(total_site);
        set_zero(ld);
        Vector<Complex> c_sum = lat_data_cget(ld);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {
            int pselsum = 0;
            for (int tslice = 0; tslice < total_site[3]; ++tslice)
            {
                const int tsnk = mod(tslice + tsep, total_site[3]);
                pselsum += psel_p_tslice[tslice] * psel_p_tslice[tsnk];
                const Complex &c_tslice_tsep = lat_data_complex_get_const(ld0, make_array(tslice, tsep))[0];
                c_sum[tsep] += c_tslice_tsep;
            }
            const double coef = -1.0 / (double)pselsum * (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
            c_sum[tsep] *= coef;
        }
        return ld;
    }

    inline LatData contract_proton_two_point_psrc_pselsnk_function(const PselProp &prop1, const Geometry &geo,
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
                c_ts[tsep] += proton_twop_block(wm, wm, wm, 0);
                c_ts[tsep] -= proton_twop_block(wm, wm, wm, 1);
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

    inline std::string get_proton_two_point_smear_path(const std::string &job_tag,
                                                       const int traj)
    {
        return ssprintf("analysis/proton_twop_smear/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_two_point_smear(const std::string &job_tag,
                                               const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_smear");
        const PointSelection &psel = get_point_selection_smear(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        LatData ld_two_point_psrc_func;
        const std::string path = get_proton_two_point_smear_path(job_tag, traj);
        const std::string path_two_point_psrc =
            path + ssprintf("/two-point-smear.lat");
        const std::string path_two_point_psrc_sum =
            path + ssprintf("/two-point-smear-sum.lat");
        std::string path_psel_p_tslice =
            path + ssprintf("/psel-per-tslice.txt");

        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        std::vector<int> psel_p_tslice(total_site[3]);
        // std::vector<int> pselpertsep(total_site[3]);
        set_zero(psel_p_tslice);
        // set_zero(pselpertsep);
        for (long n = 0; n < n_points; n++)
        {
            const long xg_src_psel_idx = n;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int &tslice = xg_src[3];
            psel_p_tslice[tslice] += 1;
            const PselProp &sp_prop = get_psel_prop_smear(job_tag, traj, xg_src, 0, 0, 1);
            LatData ld = contract_proton_two_point_psrc_pselsnk_function(sp_prop, geo, tslice, psel);
            ld_two_point_psrc_func += ld;
        }
        lat_data_save_info(path_two_point_psrc, ld_two_point_psrc_func);
        LatData ld_two_point_psrc_func_sum =
            contract_proton_psel_correction(ld_two_point_psrc_func, psel_p_tslice, geo);
        lat_data_save_info(path_two_point_psrc_sum, ld_two_point_psrc_func_sum);
        save_psel_per_tslice(psel_p_tslice, path_psel_p_tslice);
    }

    inline void compute_proton_two_point_smear_func(const std::string &job_tag, const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_two_point_smear_func");
        const std::string path = get_proton_two_point_smear_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(
                path + "/two-point-smear.lat"))
        {
            return;
        }
        if (not(check_prop_smear(job_tag, traj, 0)))
        {
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-two-point-smear-func-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_twop_smear");
        qmkdir_info(ssprintf("analysis/proton_twop_smear/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_proton_two_point_smear(job_tag, traj);
        release_lock();
    }

} // namespace qlat
