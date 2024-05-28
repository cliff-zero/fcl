#pragma once

#include "data-load.h"
#include "contract_proton.h"
#include "compute-utils.h"

namespace qlat
{
#ifndef PROTONPPBLOCK
#define PROTONPPBLOCK
    struct ProtonPPBlock
    {
        bool is_build;
        array<WilsonMatrix, 5> block_pp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        ProtonPPBlock()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            for (int i = 0; i < 5; i++)
            {
                set_zero(block_pp[i]);
            }
        }
    };
    struct ProtonSPBlock
    {
        bool is_build;
        array<WilsonMatrix, 8 * 5> block_sp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        ProtonSPBlock()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            for (int i = 0; i < 5; ++i)
            {
                for (int j = 0; j < 8; ++j)
                {
                    set_zero(block_sp[j + 5 * i]);
                }
            }
        }
    };
#endif

    inline LatData mk_proton_four_point_split_table(const int &out_size)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("ipsel", 0, out_size - 1));
        ld.info.push_back(lat_dim_number("gram", 0, 9));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    typedef std::vector<Coordinate> PselSelection;

    inline void load_psel_selection(const std::string &path, PointSelection &pssel)
    {
        TIMER_VERBOSE("load_psel_selection");
        const std::vector<std::string> lines = qgetlines(path);

        qassert(lines.size() > 0);
        const long len = read_long(lines[0]);
        qassert(len + 1 <= (long)lines.size());
        for (long k = 1; k < len + 1; ++k)
        {
            const std::vector<std::string> strs = split_line_with_spaces(lines[k]);
            if (strs.size() >= 5)
            {
                qassert(k - 1 == read_long(strs[0]));
                const Coordinate xg(read_long(strs[1]), read_long(strs[2]),
                                    read_long(strs[3]), read_long(strs[4]));
                pssel.push_back(xg);
            }
            else
            {
                displayln(fname + ssprintf(": line is '%s'.", lines[k].c_str()));
                qassert(false);
            }
        }

        return;
    }

    inline void load_prop_psrc_fourp(const std::string &job_tag, const int traj, const int tsep,
                                     const int num_dtxy, const int yt, const PointSelection &psel,
                                     const Geometry &geo)
    {
        TIMER_VERBOSE("load_prop_psrc_fourp");
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const int num_t = 2 * tsep + num_dtxy;
        // const int num_t = total_site[3];
        const int t_begin = mod(yt - num_t / 2, total_site[3]);
        for (int ti = 0; ti < num_t; ++ti)
        {
            const int tsrc = mod(t_begin + ti, total_site[3]);
            for (long nsrc = 0; nsrc < n_points; ++nsrc)
            {
                const long xg_src_psel_idx = nsrc;
                const Coordinate &xg_src = psel[xg_src_psel_idx];
                if (xg_src[3] != tsrc)
                {
                    continue;
                }
                const SelProp &prop = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
                qassert(prop.initialized);
            }
        }
        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long xg_src_psel_idx = nsrc;
            const Coordinate xg_src = psel[xg_src_psel_idx];
            if (xg_src[3] != t_begin)
            {
                continue;
            }
            const SelProp &prop = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
            qassert(prop.initialized);
            break;
        }
    }

    inline void contract_proton_pp_block(std::vector<int> &psel_num_list,
                                         std::vector<std::vector<double>> &psel_dt_num_list,
                                         const std::string &job_tag, const int traj,
                                         const Geometry &geo, const PointSelection &psel,
                                         const int dtmax, const int dtmin, const LatData &coef_table)
    {
        TIMER_VERBOSE("contract_proton_pp_block");
        const int num_dt = dtmax - dtmin + 1;
        qassert(num_dt > 0);
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const double dx_max = sqrt(32.0 * 32.0 + 12.0 * 12.0 * 3.0);

        psel_num_list.resize(total_site[3]);
        set_zero(psel_num_list);
        psel_dt_num_list.resize(num_dt);
        for (int dt = 0; dt < num_dt; ++dt)
        {
            const int del_t = dtmin + dt;
            psel_dt_num_list[dt].resize(del_t + 1);
            set_zero(psel_dt_num_list[dt]);
        }
        std::vector<int> list_n_from_idx(n_points);
        for (long n = 0; n < n_points; ++n)
        {
            const long xg_psel_idx = n;
            const Coordinate xg = psel[xg_psel_idx];
            const int t = xg[3];
            list_n_from_idx[xg_psel_idx] = psel_num_list[t];
            psel_num_list[t] += 1;
        }

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
                const int dt = del_t - dtmin;
                if ((dt < 0) || (dt >= num_dt))
                {
                    continue;
                }
                const Coordinate xg_sep1 = smod(xg_snk - xg_src, total_site);
                const int dist1 = ceil(sqrt(sum(xg_sep1 * xg_sep1)) / dx_max * 80.0);
                for (long ny = 0; ny < n_points; ++ny)
                {
                    const long xg_y_psel_idx = ny;
                    const Coordinate &xg_y = psel[xg_y_psel_idx];
                    const int ty = xg_y[3];
                    const int dt_y_src = mod(ty - tsrc, total_site[3]);
                    const int dt_snk_y = mod(tsnk - ty, total_site[3]);
                    if ((dt_y_src != (dtmin / 2)) && (dt_snk_y != (dtmin / 2)))
                    {
                        continue;
                    }
                    const Coordinate xg_sep2 = smod(xg_y - xg_src, total_site);
                    const Coordinate xg_sep3 = smod(xg_snk - xg_y, total_site);
                    const int dist2 = ceil(sqrt(sum(xg_sep2 * xg_sep2)) / dx_max * 80.0);
                    const int dist3 = ceil(sqrt(sum(xg_sep3 * xg_sep3)) / dx_max * 80.0);
                    const double coef_sphere = real(lat_data_complex_get_const(coef_table,
                                                                               make_array(dist1, dist2, dist3))[0]);
                    psel_dt_num_list[dt][dt_y_src] += coef_sphere;
                }
            }
        }
        return;
    }

    inline void contract_proton_fourp_unshifted_psel_acc_x_typeI(
        std::vector<Complex> &vv, const Complex coef, const WilsonMatrix &wm_x_snk,
        const WilsonMatrix &wm_x_y, const ProtonPPBlock &block)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;

        for (int gram = 0; gram < 5; gram++)
        {
            const WilsonMatrix wm_pp_y_tsrc_tsnk_x = block.block_pp[gram] * wm_snk_x;

            const WilsonMatrix wm_pp = wm_x_y * va_ms[3] * wm_pp_y_tsrc_tsnk_x;

            vv[gram] += coef * matrix_trace(wm_pp, va_ms[3]);
        }
        return;
    }

    inline void contract_proton_fourp_unshifted_psel_acc_x_typeII(
        std::vector<Complex> &vv, const Complex coef, const WilsonMatrix &wm_x_snk,
        const WilsonMatrix &wm_x_src, const ProtonSPBlock &block)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_snk_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_snk) * gamma5;
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();

        for (int gram = 0; gram < 5; gram++)
        {
            const WilsonMatrix wm_sp = wm_x_src * block.block_sp[3 + gram * 8] * wm_snk_x;

            vv[gram] += coef * matrix_trace(wm_sp, va_ms[3]);
        }
    }

    inline void contract_proton_split_typeI(LatData &ld_fourp_split,
                                            const FieldSelection &fsel, const Geometry &geo,
                                            const int tx, const int tsrc, const int tsnk,
                                            const SelProp &prop_x_snk, const SelProp &prop_x_y,
                                            const double &coef_typeI, const ProtonPPBlock &ppblock,
                                            const std::vector<int> &gram_convert_list, const int ipsel)
    {
        std::vector<Complex> gc_gram(omp_get_max_threads() * 5);
        set_zero(gc_gram);
#pragma omp parallel
        {
            std::vector<Complex> c_gram(5);
            set_zero(c_gram);
#pragma omp for
            for (long idx = 0; idx < fsel.n_elems; ++idx)
            {
                const long xg_x_idx = idx;
                const long index = fsel.indices[xg_x_idx];
                const Coordinate xl = geo.coordinate_from_index(index);
                const Coordinate xg = geo.coordinate_g_from_l(xl);
                const Coordinate &xg_x = xg;
                if (xg_x[3] != tx)
                {
                    continue;
                }
                const WilsonMatrix &wm_x_y = prop_x_y.get_elem(idx);
                const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(idx);
                if (tsrc <= tsnk)
                {
                    contract_proton_fourp_unshifted_psel_acc_x_typeI(c_gram, coef_typeI, wm_x_snk, wm_x_y, ppblock);
                }
                else
                {
                    const Complex coef1 = -1.0 * coef_typeI;
                    contract_proton_fourp_unshifted_psel_acc_x_typeI(c_gram, coef1, wm_x_snk, wm_x_y, ppblock);
                }
            }
            for (int igram = 0; igram < 5; ++igram)
            {
                gc_gram[omp_get_thread_num() * 5 + igram] = c_gram[igram];
            }
        }
        std::vector<Complex> m_gram(5);
        set_zero(m_gram);
        for (int i = 0; i < omp_get_max_threads(); ++i)
        {
            for (int igram = 0; igram < 5; ++igram)
            {

                m_gram[igram] += gc_gram[i * 5 + igram];
            }
        }
        glb_sum_double_vec(get_data(m_gram));
        Vector<Complex> c_ipsel = lat_data_complex_get(ld_fourp_split, make_array(ipsel));
        for (int gram = 0; gram < 5; ++gram)
        {
            int gram_idx = gram_convert_list[gram];
            c_ipsel[gram_idx - 1] = m_gram[gram];
        }
        return;
    }

    inline void contract_proton_split_typeII(LatData &ld_fourp_split, const int &x_sparsity,
                                             const FieldSelection &fsel, const Geometry &geo,
                                             const int tx, const int tsrc, const int tsnk,
                                             const SelProp &prop_x_snk, const SelProp &prop_x_src,
                                             const double &coef_typeII, const ProtonSPBlock &spblock,
                                             const std::vector<int> &gram_convert_list, const int ipsel)
    {
        std::vector<Complex> gc_gram(omp_get_max_threads() * 5);
        set_zero(gc_gram);
#pragma omp parallel
        {
            std::vector<Complex> c_gram(5);
            set_zero(c_gram);
#pragma omp for
            for (long idx = 0; idx < fsel.n_elems; ++idx)
            {
                if (mod(idx, x_sparsity) != 0)
                {
                    continue;
                }
                const long xg_x_idx = idx;
                const long index = fsel.indices[xg_x_idx];
                const Coordinate xl = geo.coordinate_from_index(index);
                const Coordinate xg = geo.coordinate_g_from_l(xl);
                const Coordinate &xg_x = xg;
                if (xg_x[3] != tx)
                {
                    continue;
                }
                const WilsonMatrix &wm_x_snk = prop_x_snk.get_elem(idx);
                const WilsonMatrix &wm_x_src = prop_x_src.get_elem(idx);
                if (tsrc <= tsnk)
                {
                    contract_proton_fourp_unshifted_psel_acc_x_typeII(c_gram, coef_typeII, wm_x_snk, wm_x_src, spblock);
                }
                else
                {
                    const Complex coef1 = -1.0 * coef_typeII;
                    contract_proton_fourp_unshifted_psel_acc_x_typeII(c_gram, coef1, wm_x_snk, wm_x_src, spblock);
                }
            }
            for (int igram = 0; igram < 5; ++igram)
            {
                gc_gram[omp_get_thread_num() * 5 + igram] = c_gram[igram];
            }
        }
        std::vector<Complex> m_gram(5);
        set_zero(m_gram);
        for (int i = 0; i < omp_get_max_threads(); ++i)
        {
            for (int igram = 0; igram < 5; ++igram)
            {

                m_gram[igram] += gc_gram[i * 5 + igram];
            }
        }
        glb_sum_double_vec(get_data(m_gram));
        Vector<Complex> c_ipsel = lat_data_complex_get(ld_fourp_split, make_array(ipsel));
        for (int gram = 0; gram < 5; ++gram)
        {
            int gram_idx = gram_convert_list[gram + 5];
            c_ipsel[gram_idx - 1] = m_gram[gram];
        }
        return;
    }

    inline std::string get_proton_four_point_split_path(const std::string &job_tag,
                                                        const int traj)
    {
        return ssprintf("analysis/proton_fourp_split/%s/results=%d", job_tag.c_str(),
                        traj);
    }

    inline void compute_proton_four_point_split_type(const std::string &job_tag, const int traj,
                                                     const std::vector<int> &gram_convert_list,
                                                     const int tsep)
    {
        check_sigterm();
        check_time_limit();
        Timer::autodisplay();
        const int gram_num = 10;
        qassert(gram_convert_list.size() == gram_num);
        const std::string path = get_proton_four_point_split_path(job_tag, traj);
        const std::string path_four_point_split =
            path + ssprintf("/fourp-split-dt%d.lat", tsep);
        TIMER_VERBOSE("compute_proton_four_point_func_psel_type");
        const PointSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const FieldSelection &fsel = get_field_selection(job_tag, traj);
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int dtmax = num_dtxy / 2 + 2 * tsep;
        const int dtmin = 2 * tsep;
        const double dx_max = sqrt(32.0 * 32.0 + 12.0 * 12.0 * 3.0);

        LatData coef_fourp;
        const std::string path_coef_table = ssprintf("coef-list/fourp/coef-fourp-80slice.lat");
        if (does_file_exist_sync_node(path_coef_table))
        {
            coef_fourp = lat_data_load_info(path_coef_table);
        }
        else
        {
            qassert(false);
        }
        PointSelection pssel;
        const std::string path_psel_sel = ssprintf("/thfs1/home/fengxu/wxh/proton/psel-selection/dt%d/%d.txt", tsep, traj);
        if (does_file_exist_sync_node(path_psel_sel))
        {
            load_psel_selection(path_psel_sel, pssel);
        }
        else
        {
            return;
        }

        std::vector<int> psel_num_list;
        std::vector<std::vector<double>> psel_dt_num_list;
        contract_proton_pp_block(psel_num_list, psel_dt_num_list, job_tag, traj, geo, psel, dtmax, dtmin, coef_fourp);

        const long n_psel = pssel.size();
        displayln_info(fname + ssprintf("n_psel=%ld", n_psel));
        LatData ld_fourp_split = mk_proton_four_point_split_table(300000);
        set_zero(ld_fourp_split);

        for (int yt = 0; yt < total_site[3]; ++yt)
        {
            load_prop_psrc_fourp(job_tag, traj, tsep, num_dtxy, yt, psel, geo);
            for (long ipsel = 0; ipsel < n_psel; ++ipsel)
            {
                const Coordinate psel_data = pssel[ipsel];
                const long xg_src_psel_idx = psel_data[0];
                const long xg_snk_psel_idx = psel_data[1];
                const long xg_y_psel_idx = psel_data[2];
                const Coordinate &xg_y = psel[xg_y_psel_idx];
                if (yt != xg_y[3])
                {
                    continue;
                }
                if (get_point_src_info(job_tag, traj, xg_y, 0).size() == 0)
                {
                    continue;
                }
                Timer::autodisplay();
                TIMER_VERBOSE("compute_proton_four_point_func_noama_type-iter");
                displayln_info(fname + ssprintf("ipsel=%ld, src_idx=%ld, snk_idx=%ld, y_idx=%ld", ipsel,
                                                xg_src_psel_idx, xg_snk_psel_idx, xg_y_psel_idx));
                const Coordinate xg_src = psel[xg_src_psel_idx];
                const Coordinate xg_snk = psel[xg_snk_psel_idx];
                const int tsrc = xg_src[3];
                const int tsnk = xg_snk[3];
                int tx, dtxy;
                if (mod(yt - tsrc, total_site[3]) == tsep)
                {
                    tx = mod(tsnk - tsep, total_site[3]);
                    dtxy = mod(tx - yt, total_site[3]);
                }
                else
                {
                    tx = mod(tsrc + tsep, total_site[3]);
                    dtxy = -mod(yt - tx, total_site[3]);
                }
                const Coordinate xg_sep1 = smod(xg_snk - xg_src, total_site);
                const Coordinate xg_sep2 = smod(xg_y - xg_src, total_site);
                const Coordinate xg_sep3 = smod(xg_snk - xg_y, total_site);
                const int dist1 = ceil(sqrt(sum(xg_sep1 * xg_sep1)) / dx_max * 80.0);
                const int dist2 = ceil(sqrt(sum(xg_sep2 * xg_sep2)) / dx_max * 80.0);
                const int dist3 = ceil(sqrt(sum(xg_sep3 * xg_sep3)) / dx_max * 80.0);
                const double coef_cor = real(lat_data_complex_get_const(coef_fourp, make_array(dist1, dist2, dist3))[0]);

                const double coef_fsel = 1.0 / fsel.prob;
                const int dt_y = mod(yt - tsrc, total_site[3]);
                const int i_dt = abs(dtxy) + 2 * tsep - dtmin;
                const double coef_psel = 1.0 / psel_dt_num_list[i_dt][dt_y];
                const double coef_typeI = coef_cor * coef_fsel * coef_psel;
                const int x_sparsity = 4;
                const double coef_typeII = coef_cor * coef_fsel * coef_psel * (double)x_sparsity;

                const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);
                const SelProp &prop_x_src = get_prop_psrc(job_tag, traj, xg_src, 0, 0);
                const SelProp &prop_x_snk = get_prop_psrc(job_tag, traj, xg_snk, 0, 0);

                const PselProp &prop_src = get_psel_prop_psrc(job_tag, traj, xg_src, 0, 0);
                const PselProp &prop_snk = get_psel_prop_psrc(job_tag, traj, xg_snk, 0, 0);

                const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
                const WilsonMatrix &wm_snk_src = prop_src.get_elem(xg_snk_psel_idx);
                const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                const WilsonMatrix &wm_y_snk = prop_snk.get_elem(xg_y_psel_idx);
                const WilsonMatrix &wm_snk_y = gamma5 * (WilsonMatrix)matrix_adjoint(wm_y_snk) * gamma5;

                ProtonPPBlock ppblock;
                proton_fourp_block_pp(ppblock.block_pp, wm_snk_src, wm_snk_src);
                for (int gram = 0; gram < (int)ppblock.block_pp.size(); ++gram)
                {
                    ppblock.block_pp[gram] = (wm_y_src * ppblock.block_pp[gram]);
                }

                ProtonSPBlock spblock;
                proton_fourp_block_sp(spblock.block_sp, wm_snk_src, wm_y_src, wm_snk_y);
                contract_proton_split_typeI(ld_fourp_split, fsel, geo, tx, tsrc, tsnk,
                                            prop_x_snk, prop_x_y, coef_typeI, ppblock,
                                            gram_convert_list, ipsel);
                contract_proton_split_typeII(ld_fourp_split, x_sparsity, fsel, geo, tx, tsrc, tsnk,
                                             prop_x_snk, prop_x_src, coef_typeII, spblock,
                                             gram_convert_list, ipsel);
            }
        }

        const double coef_wall = (double)total_site[0] * (double)total_site[1] * (double)total_site[2];
        ld_fourp_split *= (coef_wall * coef_wall);
        lat_data_save_info(path_four_point_split, ld_fourp_split);
    }

    inline void compute_proton_four_point_split(const std::string &job_tag,
                                                const int traj)
    {
        check_sigterm();
        Timer::autodisplay();
        TIMER_VERBOSE("compute_proton_four_point_split");
        const std::string path = get_proton_four_point_split_path(job_tag, traj);
        const std::string path_checkpoint = path + "/checkpoint.txt";
        if (does_file_exist_sync_node(path_checkpoint))
        {
            return;
        }
        const int tsep = 1;
        if (does_file_exist_sync_node(path + ssprintf("/fourp-split-dt%d.lat", tsep)))
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
                ssprintf("lock-proton-fourp-split-%s-%d", job_tag.c_str(), traj)))
        {
            return;
        }
        setup(job_tag, traj);
        qmkdir_info("analysis/proton_fourp_split");
        qmkdir_info(ssprintf("analysis/proton_fourp_split/%s", job_tag.c_str()));
        qmkdir_info(path);
        std::vector<int> gram_convert_list;
        gram_convert_list.push_back(3);
        gram_convert_list.push_back(4);
        gram_convert_list.push_back(8);
        gram_convert_list.push_back(9);
        gram_convert_list.push_back(10);
        gram_convert_list.push_back(1);
        gram_convert_list.push_back(2);
        gram_convert_list.push_back(5);
        gram_convert_list.push_back(6);
        gram_convert_list.push_back(7);
        compute_proton_four_point_split_type(job_tag, traj, gram_convert_list, tsep);
        release_lock();
    }
}
