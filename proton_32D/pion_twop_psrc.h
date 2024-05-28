#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{ //
    template <class T>
    qacc Complex pion_twop_block(const WilsonMatrixT<T> &m1)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &g55 = va_ms[7];
        const SpinMatrix &g5 = gamma5;
        Complex ret = 0;
        const WilsonMatrix m1g = gamma5 * (WilsonMatrix)matrix_adjoint(m1) * gamma5;
        const WilsonMatrix m11 = m1 * g5 * m1g;
        ret = matrix_trace(m11 * g5);

        return ret;
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
        std::vector<Complex> gc_ts(omp_get_max_threads() *
                                   total_site[3]);
        set_zero(gc_ts);
#pragma omp parallel
        {
            std::vector<Complex> c_ts(total_site[3]);
            set_zero(c_ts);
#pragma omp for
            for (long idx = 0; idx < (long)fsel.indices.size(); ++idx)
            {
                const long index = fsel.indices[idx];
                const Coordinate xl = geo.coordinate_from_index(index);
                const Coordinate xg = geo.coordinate_g_from_l(xl);
                const int tsep = mod(xg[3] - tslice, total_site[3]);
                const WilsonMatrix wm = prop1.get_elem(idx);
                c_ts[tsep] += pion_twop_block(wm);
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
        glb_sum_double_vec(get_data(m_ts));
        LatData ld = mk_proton_two_point_table(total_site);
        set_zero(ld);
        for (int tsep = 0; tsep < total_site[3]; ++tsep)
        {
            Vector<Complex> m_src_snk = lat_data_complex_get(ld, make_array(tslice));

            m_src_snk[tsep] += m_ts[tsep];
        }
        ld *= 1.0 / fsel.prob;
        return ld;
    }
    inline std::string get_pion_two_point_psrc_path(const std::string &job_tag,
                                                    const int traj)
    {
        return ssprintf("analysis/pion_twop/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_pion_two_point_psrc(const std::string &job_tag,
                                            const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_pion_two_point_psrc");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        LatData ld_two_point_psrc_func;
        const std::string path = get_pion_two_point_psrc_path(job_tag, traj);
        const std::string path_two_point_psrc =
            path + ssprintf("/two-point-psrc.lat");

        long iter = 0;
        for (long n = 0; n < n_points; n++)
        {
            const long xg_src_psel_idx = n;
            const Coordinate &xg_src = psel[xg_src_psel_idx];
            const int &tslice = xg_src[3];
            iter += 1;
            const SelProp &pprop = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
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
	if (not(check_prop_smear(job_tag, traj, 0)))
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
        qmkdir_info("analysis/pion_twop");
        qmkdir_info(ssprintf("analysis/pion_twop/%s", job_tag.c_str()));
        qmkdir_info(path);
        compute_pion_two_point_psrc(job_tag, traj);
        release_lock();
    }

} // namespace qlat
