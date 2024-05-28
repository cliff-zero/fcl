#pragma once

#include <qlat/contract-pion.h>
#include "new-matrix.h"

namespace qlat
{ //
    template <class T>
    qacc Complex proton_twop_block(const WilsonMatrixT<T> &m1, const WilsonMatrixT<T> &m2, const WilsonMatrixT<T> &m3, const int type)
    {
        const SpinMatrix &Cg5 = ProtonMatrixConstants::get_Cgamma5();
        const SpinMatrix &proj = ProtonMatrixConstants::get_projector();
        Complex ret = 0;
        if (type == 0)
        {
            const WilsonMatrix wm1 = m1 * Cg5;
            const WilsonMatrix wm2 = Cg5 * m2;
            const WilsonMatrix wm3 = m3 * proj;
            ret = lv_wm3(wm1, wm2, wm3, 0);
        }
        else if (type == 1)
        {
            const WilsonMatrix wm1 = m1 * Cg5;
            const WilsonMatrix wm2 = Cg5 * m2;
            const WilsonMatrix wm3 = m3 * proj;
            ret = lv_wm3(wm1, wm2, wm3, 1);
        }
        else
        {
            qassert(false);
        }
        return ret;
    };

    qacc int get_proton_fourp_max_sep(int t)
    {
        if (t == 16)
        {
            return 3;
        }
        else if (t == 64)
        {
            return 23;
        }
        else
        {
            qassert(false);
            return 0;
        }
    }

    qacc WilsonMatrix proton_fourp_block_type(const WilsonMatrix &m1, const WilsonMatrix &m2, const int &type)
    {
        const SpinMatrix &Cg5 = ProtonMatrixConstants::get_Cgamma5();
        const SpinMatrix &proj = ProtonMatrixConstants::get_projector();
        WilsonMatrix block;
        WilsonMatrix wm1, wm2;
        set_zero(block);

        switch (type)
        {
        case 1:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = m2;
            block = lv_wm2(wm1, wm2, proj, type);
            break;

        case 2:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = m2 * proj;
            block = lv_wm2(wm1, wm2, proj, type);
            break;

        case 3:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = proj * m2;
            block = lv_wm2(wm1, wm2, proj, type);
            break;

        case 4:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = m2 * proj;
            block = lv_wm2(wm1, wm2, proj, type);
            break;

        case 5:
            wm1 = Cg5 * m1 * proj;
            wm2 = m2 * Cg5;
            block = lv_wm2(wm1, wm2, proj, type);
            break;

        default:
            qassert(false);
            break;
        }
        return block;
    };

    qacc void proton_fourp_block_pp(array<WilsonMatrix, 5> &block, const WilsonMatrix &m1, const WilsonMatrix &m2)
    {
        for (int type = 1; type < 6; type++)
        {
            block[type - 1] = proton_fourp_block_type(m1, m2, type);
        }
    };

    qacc void proton_fourp_block_sp(array<WilsonMatrix, 8 * 5> &block, const WilsonMatrix &m1, const WilsonMatrix &m2, const WilsonMatrix &m3)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        for (int mu = 0; mu < 8; mu++)
        {
            const WilsonMatrix wm = m3 * va_ms[mu] * m2;
            block[mu + 0 * 8] = proton_fourp_block_type(m1, wm, 1);
            block[mu + 1 * 8] = proton_fourp_block_type(wm, m1, 2);
            block[mu + 2 * 8] = proton_fourp_block_type(wm, m1, 4);
            block[mu + 3 * 8] = proton_fourp_block_type(wm, m1, 3);
            block[mu + 4 * 8] = proton_fourp_block_type(m1, wm, 4);
        }
    };

    inline LatData mk_proton_two_point_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    //######################################################################################################

    inline LatData mk_proton_three_point_table(const Coordinate &total_site)
    {
        LatData ld;
        ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("top", 0, total_site[3] - 1));
        ld.info.push_back(lat_dim_number("gram", 0, 4));
        ld.info.push_back(lat_dim_number("op", 0, 15));
        ld.info.push_back(lat_dim_re_im());
        lat_data_alloc(ld);
        set_zero(ld);
        return ld;
    }

    inline LatData contract_proton_three_point_function(const SelProp &prop_a,
                                                        const SelProp &prop_b,
                                                        const WilsonMatrix &wm_ab,
                                                        const WilsonMatrix &wm_ba,
                                                        const int ta, const int tb,
                                                        const FieldSelection &fsel)
    // ``wm_ab'' is prop from ``tb'' to ``ta''.
    // prop_a (type1)
    // prop_b (type2)
    // wm_ab (type3)
    {
        TIMER("contract_proton_three_point_function");
        const array<SpinMatrix, 16> &gms = SpinMatrixConstants::get_cps_gms();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        qassert(is_matching_geo_mult(prop_a.geo(), geo));
        qassert(is_matching_geo_mult(prop_b.geo(), geo));
        const int tsep = mod(tb - ta, total_site[3]);

        const WilsonMatrix wmd = gamma5 * (WilsonMatrix)matrix_adjoint(wm_ab) * gamma5;
        array<WilsonMatrix, 5> block;
        proton_fourp_block_pp(block, wm_ba, wmd);

        std::vector<WilsonMatrix> gwm_ts(omp_get_max_threads() * total_site[3] * 5);
        set_zero(gwm_ts);
#pragma omp parallel
        {
            std::vector<WilsonMatrix> wm_ts(total_site[3] * 5);
            set_zero(wm_ts);
#pragma omp for
            for (long idx = 0; idx < (long)fsel.indices.size(); ++idx)
            {
                const long index = fsel.indices[idx];
                const Coordinate xl = geo.coordinate_from_index(index);
                const Coordinate xg = geo.coordinate_g_from_l(xl);
                if ((mod(xg[3] - ta, total_site[3]) < tsep) && (mod(xg[3] - ta, total_site[3]) > 0))
                {
                    for (int gram = 0; gram < 5; gram++)
                    {
                        wm_ts[xg[3] * 5 + gram] +=
                            (WilsonMatrix)(prop_a.get_elem(idx) * block[gram]) *
                            (gamma5 * (WilsonMatrix)matrix_adjoint(prop_b.get_elem(idx)) *
                             gamma5);
                    }
                }
            }
            for (int t = 0; t < total_site[3]; ++t)
            {
                for (int gram = 0; gram < 5; gram++)
                {
                    gwm_ts[omp_get_thread_num() * total_site[3] * 5 + t * 5 + gram] = wm_ts[t * 5 + gram];
                }
            }
        }
        std::vector<WilsonMatrix> wm_ts(total_site[3] * 5);
        set_zero(wm_ts);
        for (int i = 0; i < omp_get_max_threads(); ++i)
        {
            for (int t = 0; t < total_site[3]; ++t)
            {
                for (int gram = 0; gram < 5; gram++)
                {
                    wm_ts[t * 5 + gram] += gwm_ts[i * total_site[3] * 5 + t * 5 + gram];
                }
            }
        }
        glb_sum_double_vec(get_data(wm_ts));
        LatData ld = mk_proton_three_point_table(total_site);
        for (int t = 0; t < total_site[3]; ++t)
        {
            const int top = mod(t - ta, total_site[3]);
            if ((top > 0) && (top < tsep))
            {
                for (int gram = 0; gram < 5; gram++)
                {
                    Vector<Complex> v = lat_data_complex_get(ld, make_array(tsep, top, gram));
                    for (int op = 0; op < 16; ++op)
                    {
                        v[op] = matrix_trace(wm_ts[t * 5 + gram], gms[op]);
                    }
                }
            }
        }
        ld *= 1.0 / fsel.prob;
        return ld;
    }

    //##########################################################################################################

    struct ProtonFourBlock
    {
        bool is_build = false;
        array<WilsonMatrix, 5> block_pp;
        array<WilsonMatrix, 8 * 5> block_sp;
        WilsonMatrix block_twop;

        ProtonFourBlock(){};
    };

    //##############################################################################################
    // contract proton two point point_source function
    inline LatData contract_proton_two_point_psrc_function(const SelProp &prop1,
                                                           const int tslice,
                                                           const FieldSelection &fsel)
    // m_ts[tsep][op_src][op_snk] = trace( (\sum_x prop1(x) gms[op_src] gamma5
    // prop2(x)^\dagger gamma5) gms[op_snk] ) 0 <= tsep < total_site[3]
    {
        TIMER_VERBOSE("contract_proton_two_point_psrc_function");
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

    //###############################################################################################
    // contract va two point

    inline void contract_va_twop_unshifted_acc_x(
        Vector<Complex> vv, const Complex coef, const WilsonMatrix &wm_x_y, const Coordinate &xg_x,
        const long xg_x_idx, const Coordinate &xg_y, const long xg_y_psel_idx)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm_y_x =
            gamma5 * (WilsonMatrix)matrix_adjoint(wm_x_y) * gamma5;
        for (int mu = 0; mu < 8; mu++)
        {
            const WilsonMatrix wm = wm_x_y * va_ms[mu] * wm_y_x;
            for (int nu = 0; nu < 8; nu++)
            {
                const int mu_nu = 8 * nu + mu;
                vv[mu_nu] += coef * matrix_trace(wm, va_ms[nu]);
            }
        }
    }

    inline void contract_va_twop_unshifted(
        SelectedField<Complex> &va_twop,
        const SelProp &prop_x_y, const Coordinate &xg_y,
        const long xg_y_psel_idx, const PointSelection &psel,
        const FieldSelection &fsel)
    // fsel.prob is NOT accounted.
    {
        TIMER_VERBOSE("contract_va_twop_unshifted");
        qassert(psel[xg_y_psel_idx] == xg_y);
        const Geometry &geo = fsel.f_rank.geo();
        const Coordinate total_site = geo.total_site();
        const int multiplicity = 8 * 8;

        va_twop.init(fsel, multiplicity);
        set_zero(va_twop);

        const int num_dtxy = get_proton_fourp_max_sep(total_site[3]);
        const int num_dxxy = num_dtxy + 18;

#pragma omp parallel
        {
#pragma omp for
            for (long idx = 0; idx < fsel.n_elems; ++idx)
            {
                const long xg_x_idx = idx;
                const long index = fsel.indices[idx];
                const Coordinate xl = geo.coordinate_from_index(index);
                const Coordinate xg = geo.coordinate_g_from_l(xl);
                const Coordinate &xg_x = xg;
                const Coordinate xg_sep = smod(xg_x - xg_y, total_site);

                if ((abs(xg_sep[0]) < num_dxxy / 2 + 1) && (abs(xg_sep[1]) < num_dxxy / 2 + 1) &&
                    (abs(xg_sep[2]) < num_dxxy / 2 + 1) && (abs(xg_sep[3]) < num_dtxy / 2 + 1))
                {

                    const WilsonMatrix &wm_x_y = prop_x_y.get_elem(idx);
                    Vector<Complex> vv;

                    vv = va_twop.get_elems(idx);

                    contract_va_twop_unshifted_acc_x(vv, 1.0, wm_x_y, xg_x, xg_x_idx,
                                                     xg_y, xg_y_psel_idx);
                }
            }
        }
        sync_node();
    }

    inline void contract_va_two_point_acc(
        FieldM<Complex, 8 * 8> &fourp_data, const std::string &job_tag, const int traj,
        const Coordinate &xg_y, const long xg_y_psel_idx, const PointSelection &psel,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_va_two_point_acc");

        const std::vector<PointInfo> &pis_xgt =
            get_point_src_info(job_tag, traj, xg_y, 0);
        if (pis_xgt.size() == 0)
        {
            return;
        }
        const int num_acc = pis_xgt.size();
        qassert(num_acc >= 1);
        SelectedField<Complex> sfama;
        const Geometry &geo = fsel.f_rank.geo();
        qassert(psel[xg_y_psel_idx] == xg_y);
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);

        if (num_acc == 1)
        {
            qassert(pis_xgt[0].accuracy == 0);
            const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);
            qassert(geo == prop_x_y.geo());
            qassert(fsel.n_elems == prop_x_y.n_elems);
            qassert(is_initialized(prop_x_y));

            contract_va_twop_unshifted(sfama, prop_x_y, xg_y,
                                       xg_y_psel_idx, psel, fsel);

            qassert(fsel.prob == ssp.fsel.prob);
            const Complex coef = 1.0 / fsel.prob;
            SelectedField<Complex> &s_data = sfama;
            acc_field(fourp_data, coef, s_data, ssp);
            return;
        }

        const TypeAccuracyTable &tat = get_type_accuracy_table(job_tag, traj);
        qassert(get_accuracy_weight(tat, 0, 0) == 1.0);
        qassert(num_acc > 1);
        std::vector<double> coefs(num_acc, 0.0);
        coefs[0] = 1.0;
        for (int acc = 1; acc < num_acc; ++acc)
        {
            const double weight = get_accuracy_weight(tat, 0, acc);
            coefs[acc] += weight;
            coefs[acc - 1] -= weight;
        }
        for (int acc = 0; acc < num_acc; ++acc)
        {
            qassert(pis_xgt[0].accuracy == 0);
            const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, acc);
            qassert(geo == prop_x_y.geo());
            qassert(fsel.n_elems == prop_x_y.n_elems);
            qassert(is_initialized(prop_x_y));
            SelectedField<Complex> sf;
            contract_va_twop_unshifted(sf, prop_x_y, xg_y,
                                       xg_y_psel_idx, psel, fsel);
            sf *= coefs[acc];
            sfama += sf;
        }

        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;
        SelectedField<Complex> &s_data = sfama;
        acc_field(fourp_data, coef, s_data, ssp);
    }
    inline void contract_va_two_point_acc_noama(
        FieldM<Complex, 8 * 8> &fourp_data, const std::string &job_tag, const int traj,
        const Coordinate &xg_y, const long xg_y_psel_idx, const PointSelection &psel,
        const FieldSelection &fsel, const ShiftShufflePlan &ssp)
    // xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
    // ssp = make_shift_shuffle_plan(fsel, -xg_y);
    {
        TIMER_VERBOSE("contract_va_two_point_acc_noama");
        const Geometry &geo = fsel.f_rank.geo();
        const SelProp &prop_x_y = get_prop_psrc(job_tag, traj, xg_y, 0, 0);
        qassert(geo == prop_x_y.geo());
        qassert(fsel.n_elems == prop_x_y.n_elems);
        qassert(is_initialized(prop_x_y));
        qassert(psel[xg_y_psel_idx] == xg_y);
        qassert(ssp.shift == -xg_y);
        qassert(ssp.is_reflect == false);
        SelectedField<Complex> sf;
        contract_va_twop_unshifted(sf, prop_x_y, xg_y,
                                   xg_y_psel_idx, psel, fsel);
        qassert(fsel.prob == ssp.fsel.prob);
        const Complex coef = 1.0 / fsel.prob;

        SelectedField<Complex> &s_data = sf;
        acc_field(fourp_data, coef, s_data, ssp);
    }
    //#################################################################################################
}
