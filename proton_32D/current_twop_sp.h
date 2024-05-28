#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{ //

    inline LatData mk_current_two_point_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("isrc", 0, 2));
        ld.info.push_back(lat_dim_number("isnk", 0, 2));
        ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline void pion_twop_block(const WilsonMatrix &m1, std::vector<Complex> &ret1)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const WilsonMatrix m1g = gamma5 * (WilsonMatrix)matrix_adjoint(m1) * gamma5;

        for (int mu = 0; mu < 3; mu++)
        {
            const WilsonMatrix m11 = m1 * va_ms[mu] * m1g;
            for (int nu = 0; nu < 3; nu++)
            {
                const int mu_nu = 3 * nu + mu;
                ret1[mu_nu] += matrix_trace(m11, va_ms[nu]);
            }
        }
        return;
    };

    inline LatData contract_pion_two_point_psrc_function(const SelProp &prop1,
                                                         const int tslice,
                                                         const FieldSelection &fsel)
    // m_ts[tsep][op_src][op_snk] = trace( (\sum_x prop1(x) gms[op_src] gamma5
    // prop2(x)^\dagger gamma5) gms[op_snk] ) 0 <= tsep < total_site[3]
    {
        TIMER_VERBOSE("contract_pion_two_point_psrc_function");
        const Geometry &geo = prop1.geo();
        const Coordinate total_site = geo.total_site();
        std::vector<std::vector<Complex>> gc_ts(omp_get_max_threads() *
                                                total_site[3]);
        for (int i = 0; i < omp_get_max_threads() * total_site[3]; ++i)
        {
            gc_ts[i].resize(9);
            set_zero(gc_ts[i]);
        }
#pragma omp parallel
        {
            std::vector<std::vector<Complex>> c_ts(total_site[3]);
            for (int i = 0; i < total_site[3]; ++i)
            {
                c_ts[i].resize(9);
                set_zero(c_ts[i]);
            }
#pragma omp for
            for (long idx = 0; idx < (long)fsel.indices.size(); ++idx)
            {
                const long index = fsel.indices[idx];
                const Coordinate xl = geo.coordinate_from_index(index);
                const Coordinate xg = geo.coordinate_g_from_l(xl);
                const int tsep = mod(xg[3] - tslice, total_site[3]);
                const WilsonMatrix wm = prop1.get_elem(idx);
                pion_twop_block(wm, c_ts[tsep]);
            }
            for (int t = 0; t < total_site[3]; ++t)
            {
                for (int i = 0; i < 9; ++i)
                {
                    gc_ts[omp_get_thread_num() * total_site[3] + t][i] = c_ts[t][i];
                }
            }
        }
        std::vector<Complex> m_ts(total_site[3] * 9);

        set_zero(m_ts);

        for (int i = 0; i < omp_get_max_threads(); ++i)
        {
            for (int t = 0; t < total_site[3]; ++t)
            {
                for (int iss = 0; iss < 9; ++iss)
                {
                    m_ts[t + iss * total_site[3]] += gc_ts[i * total_site[3] + t][iss];
                }
            }
        }
        glb_sum_double_vec(get_data(m_ts));
        LatData ld = mk_current_two_point_table(total_site);
        set_zero(ld);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {
            for (int isrc = 0; isrc < 3; ++isrc)
            {
                for (int isnk = 0; isnk < 3; ++isnk)
                {
                    Vector<Complex> m_src_snk = lat_data_complex_get(ld, make_array(isrc, isnk, tslice));

                    m_src_snk[tsep] += m_ts[tsep + (isrc + isnk * 3) * total_site[3]];
                }
            }
        }
        ld *= 1.0 / fsel.prob;
        return ld;
    }
    inline std::string get_pion_two_point_psrc_path(const std::string &job_tag,
                                                    const int traj)
    {
        return ssprintf("analysis/current_twop_sp/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_pion_two_point_psrc(const std::string &job_tag,
                                            const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_pion_two_point_psrc");
        const PointSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points = psel_smear.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        LatData ld_two_point_psrc_func;
        const std::string path = get_pion_two_point_psrc_path(job_tag, traj);
        const std::string path_two_point_psrc =
            path + ssprintf("/two-point-psrc.lat");

        long iter = 0;
        for (long n = 0; n < n_points; n++)
        {
            const long xg_src_psel_idx = n;
            const Coordinate &xg_src = psel_smear[xg_src_psel_idx];
            const int &tslice = xg_src[3];
            iter += 1;
            const SelProp &pprop = get_prop_smear(job_tag, traj, xg_src, 0, 0, 0);
            LatData ld = contract_pion_two_point_psrc_function(pprop, tslice, fsel);
            ld_two_point_psrc_func += ld;
        }
        const long n_iter = iter;
        const double coef = 1.0 / (double)n_iter;
        ld_two_point_psrc_func *= coef;
        displayln_info(ssprintf("n_iter=%ld", n_iter));
        lat_data_save_info(path_two_point_psrc, ld_two_point_psrc_func);
    }

    inline void compute_pion_two_point_psrc_func(const std::string &job_tag, const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_pion_two_point_psrc_func");
        const std::string path = get_pion_two_point_psrc_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        if (does_file_exist_sync_node(
                path + "/two-point-psrc.lat"))
        {
            return;
        }
        check_sigterm();
        check_time_limit();
        if (not obtain_lock(
                ssprintf("lock-two-point-fsel-func-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/current_twop_sp");
        qmkdir_info(ssprintf("analysis/current_twop_sp/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_pion_two_point_psrc(job_tag, traj);
        release_lock();
    }

} // namespace qlat