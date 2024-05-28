#pragma once

#include <qlat/selected-points.h>
#include "data-load.h"
#include "contract_block.h"
#include "compute-utils.h"
#include "four-point-time-system.h"

#include <algorithm>
#include <iterator>

// #pragma omp declare reduction(+ : std::vector<qlat::ComplexT> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<qlat::ComplexT>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

// #pragma omp declare reduction(+ : std::vector<long> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long>()))  initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

namespace qlat
{
    inline const SelProp &get_prop(const std::string &job_tag, const int &traj, const Coordinate &xg, bool is_smear)
    {
        if (!is_smear)
        {
            return get_prop_psrc(job_tag, traj, xg, 0, 0);
        }
        else
        {
            return get_prop_smear(job_tag, traj, xg, 0, 0, 0);
        }
    }

    inline const PselProp &get_psel_prop(const std::string &job_tag, const int &traj, const Coordinate &xg, bool is_snk_smear, bool is_src_smear)
    {
        if (!is_src_smear)
        {
            if (!is_snk_smear)
            {
                return get_psel_prop_psrc(job_tag, traj, xg, 0, 0);
            }
            else
            {
                qassert(false);
            }
        }
        else
        {
            if (!is_snk_smear)
            {
                return get_psel_prop_smear(job_tag, traj, xg, 0, 0, 0);
            }
            else
            {
                return get_psel_prop_smear(job_tag, traj, xg, 0, 0, 1);
            }
        }
    }

#ifndef Baryon_Block
#define Baryon_Block
    struct Baryon_Block_Psel_Psel
    {
        bool is_build;
        array<WilsonMatrix, 5> block_pp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        Baryon_Block_Psel_Psel()
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

    struct Baryon_Block_Scalar_Psel_Psel
    {
        bool is_build;
        array<ComplexD, 2> block_pp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        Baryon_Block_Scalar_Psel_Psel()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            block_pp[0] = 0.0;
            block_pp[1] = 0.0;
        }
    };

    struct Baryon_Block_2PF_No_Projection_Psel_Psel
    {
        bool is_build;
        array<SpinMatrix, 2> block_pp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        Baryon_Block_2PF_No_Projection_Psel_Psel()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            set_zero(block_pp[0]);
            set_zero(block_pp[1]);
        }
    };

    struct Meson_Block_Scalar_Psel_Psel
    {
        bool is_build;
        ComplexD block_pp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        Meson_Block_Scalar_Psel_Psel()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            block_pp = 0.0;
        }
    };

    struct Baryon_Block_Sequential_Psel_Psel
    {
        bool is_build;
        array<WilsonMatrix, 8 * 5> block_sp;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        Baryon_Block_Sequential_Psel_Psel()
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

    struct Sequential_Prop_2pt_Block
    {
        bool is_build;
        array<WilsonMatrix, 2> block;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        long xg_x_psel_idx;
        long xg_y_psel_idx;
        int n_snk;
        int n_src;

        Sequential_Prop_2pt_Block()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            xg_x_psel_idx = -1;
            xg_y_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            set_zero(block[0]);
            set_zero(block[1]);
        }
    };

    struct Sequential_Prop_1pt_Block
    {
        bool is_build;
        WilsonMatrix block;
        long xg_par1_psel_idx;
        long xg_par2_psel_idx;
        long xg_par3_psel_idx;
        int n_snk;
        int n_src;

        Sequential_Prop_1pt_Block()
        {
            is_build = false;
            xg_par1_psel_idx = -1;
            xg_par2_psel_idx = -1;
            xg_par3_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            set_zero(block);
        }
    };

    struct Baryon_Block_No_Projection_Psel_Psel
    {
        bool is_build;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        std::vector<WilsonMatrix> block_0_12;
        WilsonMatrix block_0_3;
        std::vector<WilsonMatrix> block_1_1;
        std::vector<WilsonMatrix> block_1_2;
        std::vector<WilsonMatrix> block_1_3;

        Baryon_Block_No_Projection_Psel_Psel()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;

            block_0_12.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_0_12[mu1 * 4 + mu2]);
                }
            }

            set_zero(block_0_3);

            block_1_1.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_1_1[mu1 * 4 + mu2]);
                }
            }

            block_1_2.resize(4);
            block_1_3.resize(4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                set_zero(block_1_2[mu1]);
                set_zero(block_1_3[mu1]);
            }
        }
    };

    struct Leading_Twist_Block_No_Projection_Psel_Psel
    {
        bool is_build;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        long xg_y_psel_idx;
        int n_snk;
        int n_src;
        int n_y;

        std::vector<std::vector<WilsonMatrix>> block_0_12;
        std::vector<std::vector<WilsonMatrix>> block_0_3;
        std::vector<std::vector<WilsonMatrix>> block_1_1;
        std::vector<std::vector<WilsonMatrix>> block_1_2;
        std::vector<std::vector<WilsonMatrix>> block_1_3;

        Leading_Twist_Block_No_Projection_Psel_Psel()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            xg_y_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            n_y = -1;

            block_0_12.resize(2);
            block_0_3.resize(2);
            block_1_1.resize(2);
            block_1_2.resize(2);
            block_1_3.resize(2);

            for (int is_inv = 0; is_inv < 2; is_inv++)
            {
                block_0_12[is_inv].resize(4 * 4);
                for (int mu1 = 0; mu1 < 4; mu1++)
                {
                    for (int mu2 = 0; mu2 < 4; mu2++)
                    {
                        set_zero(block_0_12[is_inv][mu1 * 4 + mu2]);
                    }
                }

                block_0_3[is_inv].resize(4);
                for (int mu1 = 0; mu1 < 4; mu1++)
                {
                    set_zero(block_0_3[is_inv][mu1]);
                }

                block_1_1[is_inv].resize(4 * 4);
                for (int mu1 = 0; mu1 < 4; mu1++)
                {
                    for (int mu2 = 0; mu2 < 4; mu2++)
                    {
                        set_zero(block_1_1[is_inv][mu1 * 4 + mu2]);
                    }
                }
            }

            block_1_2[0].resize(4);
            block_1_3[1].resize(4);
            block_1_2[1].resize(16);
            block_1_3[0].resize(16);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                set_zero(block_1_2[0][mu1]);
                set_zero(block_1_3[1][mu1]);
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_1_2[1][mu1 * 4 + mu2]);
                    set_zero(block_1_3[0][mu1 * 4 + mu2]);
                }
            }
        }
    };

    struct Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra
    {
        bool is_build;
        long xg_snk_psel_idx;
        long xg_y_psel_idx;
        int n_snk;
        int n_src;
        int n_y;

        std::vector<WilsonMatrix> block_0_12;
        std::vector<WilsonMatrix> block_0_3;
        std::vector<WilsonMatrix> block_1_1;
        std::vector<WilsonMatrix> block_1_2;
        std::vector<WilsonMatrix> block_1_3;

        Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_y_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            n_y = -1;

            block_0_12.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_0_12[mu1 * 4 + mu2]);
                }
            }

            block_0_3.resize(4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                set_zero(block_0_3[mu1]);
            }

            block_1_1.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_1_1[mu1 * 4 + mu2]);
                }
            }

            block_1_2.resize(4);
            block_1_3.resize(16);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                set_zero(block_1_2[mu1]);
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_1_3[mu1 * 4 + mu2]);
                }
            }
        }
    };

    struct Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra
    {
        bool is_build;
        long xg_snk_psel_idx;
        long xg_src_psel_idx;
        long xg_y_psel_idx;
        int n_snk;
        int n_src;
        int n_y;

        std::vector<WilsonMatrix> block_0_12;
        std::vector<WilsonMatrix> block_0_3;
        std::vector<WilsonMatrix> block_1_1;
        std::vector<WilsonMatrix> block_1_2;
        std::vector<WilsonMatrix> block_1_3;

        Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_src_psel_idx = -1;
            xg_y_psel_idx = -1;
            n_snk = -1;
            n_src = -1;
            n_y = -1;

            block_0_12.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_0_12[mu1 * 4 + mu2]);
                }
            }

            block_0_3.resize(4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                set_zero(block_0_3[mu1]);
            }

            block_1_1.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_1_1[mu1 * 4 + mu2]);
                }
            }

            block_1_2.resize(16);
            block_1_3.resize(4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                set_zero(block_1_3[mu1]);
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(block_1_2[mu1 * 4 + mu2]);
                }
            }
        }
    };

    struct Baryon_SPBlock_No_Projection_Psel_Psel
    {
        bool is_build;
        long xg_snk_psel_idx;
        long xg_y_psel_idx;
        long xg_src_psel_idx;
        int n_snk;
        int n_src;

        WilsonMatrix spblock_5;
        std::vector<WilsonMatrix> spblock_6;
        std::vector<WilsonMatrix> spblock_7;
        std::vector<WilsonMatrix> spblock_8;
        std::vector<WilsonMatrix> spblock_9;

        std::vector<WilsonMatrix> spblock_ex_5;
        std::vector<WilsonMatrix> spblock_ex_6;
        std::vector<WilsonMatrix> spblock_ex_7;
        std::vector<WilsonMatrix> spblock_ex_8;
        std::vector<WilsonMatrix> spblock_ex_9;

        Baryon_SPBlock_No_Projection_Psel_Psel()
        {
            is_build = false;
            xg_snk_psel_idx = -1;
            xg_y_psel_idx = -1;
            xg_src_psel_idx = -1;
            n_snk = -1;
            n_src = -1;

            set_zero(spblock_5);

            spblock_7.resize(4);
            spblock_8.resize(4);
            spblock_9.resize(4);
            spblock_ex_9.resize(4);
            for (int mu = 0; mu < 4; mu++)
            {
                set_zero(spblock_7[mu]);
                set_zero(spblock_8[mu]);
                set_zero(spblock_9[mu]);
                set_zero(spblock_ex_9[mu]);
            }

            spblock_6.resize(4 * 4);
            spblock_ex_5.resize(4 * 4);
            spblock_ex_6.resize(4 * 4);
            spblock_ex_7.resize(4 * 4);
            spblock_ex_8.resize(4 * 4);
            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    set_zero(spblock_6[mu1 * 4 + mu2]);
                    set_zero(spblock_ex_5[mu1 * 4 + mu2]);
                    set_zero(spblock_ex_6[mu1 * 4 + mu2]);
                    set_zero(spblock_ex_7[mu1 * 4 + mu2]);
                    set_zero(spblock_ex_8[mu1 * 4 + mu2]);
                }
            }
        }

        Baryon_SPBlock_No_Projection_Psel_Psel &operator*=(const ComplexD &coeff)
        {
            spblock_5 *= coeff;

            for (int mu = 0; mu < 4; mu++)
            {
                spblock_7[mu] *= coeff;
                spblock_8[mu] *= coeff;
                spblock_9[mu] *= coeff;
                spblock_ex_9[mu] *= coeff;
            }

            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    spblock_6[mu1 * 4 + mu2] *= coeff;
                    spblock_ex_5[mu1 * 4 + mu2] *= coeff;
                    spblock_ex_6[mu1 * 4 + mu2] *= coeff;
                    spblock_ex_7[mu1 * 4 + mu2] *= coeff;
                    spblock_ex_8[mu1 * 4 + mu2] *= coeff;
                }
            }

            return *this;
        }

        Baryon_SPBlock_No_Projection_Psel_Psel &operator+=(const Baryon_SPBlock_No_Projection_Psel_Psel &other)
        {
            qassert(xg_snk_psel_idx == other.xg_snk_psel_idx);
            qassert(xg_src_psel_idx == other.xg_src_psel_idx);

            spblock_5 += other.spblock_5;

            for (int mu = 0; mu < 4; mu++)
            {
                spblock_7[mu] += other.spblock_7[mu];
                spblock_8[mu] += other.spblock_8[mu];
                spblock_9[mu] += other.spblock_9[mu];
                spblock_ex_9[mu] += other.spblock_ex_9[mu];
            }

            for (int mu1 = 0; mu1 < 4; mu1++)
            {
                for (int mu2 = 0; mu2 < 4; mu2++)
                {
                    spblock_6[mu1 * 4 + mu2] += other.spblock_6[mu1 * 4 + mu2];
                    spblock_ex_5[mu1 * 4 + mu2] += other.spblock_ex_5[mu1 * 4 + mu2];
                    spblock_ex_6[mu1 * 4 + mu2] += other.spblock_ex_6[mu1 * 4 + mu2];
                    spblock_ex_7[mu1 * 4 + mu2] += other.spblock_ex_7[mu1 * 4 + mu2];
                    spblock_ex_8[mu1 * 4 + mu2] += other.spblock_ex_8[mu1 * 4 + mu2];
                }
            }

            return *this;
        }
    };
#endif

    qacc int get_baryon_pi_max_sep(int t)
    {
        if (t == 16)
        {
            return 3;
        }
        else if (t == 64)
        {
            return 27;
        }
        else
        {
            qassert(false);
            return 0;
        }
    }

    qacc int get_baryon_pi_half_sep(int t)
    {
        if (t == 16)
        {
            return 2;
        }
        else if (t == 64)
        {
            return 13;
        }
        else
        {
            qassert(false);
            return 1;
        }
    }

    qacc long space_dot(const std::vector<int> &coor1, const Coordinate &coor2)
    {
        long ret = 0;
        qassert(coor1.size() == 3);
        qassert(4 == coor2.size());
        for (long unsigned int i = 0; i < 3; i++)
        {
            ret += coor1[i] * coor2[i];
        }
        return ret;
    }

    qacc long space_dot(const Coordinate &coor1, const Coordinate &coor2)
    {
        long ret = 0;
        qassert(coor1.size() == coor2.size());
        for (long unsigned int i = 0; i < coor1.size() - 1; i++)
        {
            ret += coor1[i] * coor2[i];
        }
        return ret;
    }

    inline void load_prop_baryon_pi(const std::string &job_tag, const int &traj, const int tsep,
                                    const int num_dtxy, const int t_meson, const PointsSelection &psel, const bool &is_smear,
                                    const Geometry &geo)
    {
        TIMER_VERBOSE("load_prop_baryon_pi");
        const long n_points = psel.size();
        const Coordinate total_site = geo.total_site();
        const int num_t = 2 * tsep + num_dtxy;
        // const int num_t = total_site[3];
        const int t_begin = mod(t_meson - num_t / 2, total_site[3]);
        for (int ti = 0; ti < num_t; ++ti)
        {
            const int t_baryon_src = mod(t_begin + ti, total_site[3]);
            for (long nsrc = 0; nsrc < n_points; ++nsrc)
            {
                const long &xg_src_psel_idx = nsrc;
                const Coordinate &xg_src = psel[xg_src_psel_idx];
                if (xg_src[3] != t_baryon_src)
                {
                    continue;
                }
                const SelProp &prop = get_prop(job_tag, traj, xg_src, is_smear);
                qassert(prop.initialized);
            }
        }
        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate xg_src = psel[xg_src_psel_idx];
            if (xg_src[3] != t_begin)
            {
                continue;
            }
            const SelProp &prop = get_prop(job_tag, traj, xg_src, is_smear);
            qassert(prop.initialized);
            break;
        }
    }

    inline void contract_baryon_pselpsel_block(std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
                                               const std::vector<int> &psel_baryon_num_list, const std::vector<int> &list_n_from_idx_baryon,
                                               const std::string &job_tag, const int &traj,
                                               const Geometry &geo, const PointsSelection &psel_baryon,
                                               const int &dtmax, const int &dtmin, const bool &is_baryon_smear)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_block");
        const int num_dt = dtmax - dtmin + 1;
        qassert(num_dt > 0);
        const long n_points = psel_baryon.size();
        const Coordinate total_site = geo.total_site();
        block.resize(n_points);
        for (long n = 0; n < n_points; ++n)
        {
            block[n].resize(num_dt);
        }

        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long &xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            const int t_baryon_snk = xg_snk[3];
            for (int dt = 0; dt < num_dt; ++dt)
            {
                const int del_t = dtmin + dt;
                const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_baryon_num_list[t_baryon_src]);
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            // const int t_baryon_src = xg_src[3];
            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);
        }

#pragma omp parallel for
        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            const int t_baryon_src = xg_src[3];
            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long &xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                const int t_baryon_snk = xg_snk[3];
                const int del_t = mod(t_baryon_snk - t_baryon_src, total_site[3]);
                const int &dt = del_t - dtmin;
                if ((dt < 0) || (dt >= num_dt))
                {
                    continue;
                }
                Baryon_Block_Psel_Psel ppblock;
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);
                baryon_pi_block_pp(ppblock.block_pp, wm, wm);
                ppblock.xg_snk_psel_idx = xg_snk_psel_idx;
                ppblock.xg_src_psel_idx = xg_src_psel_idx;
                ppblock.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                ppblock.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                ppblock.is_build = true;

                const int n_t = list_n_from_idx_baryon[xg_src_psel_idx];
                block[xg_snk_psel_idx][dt][n_t] = ppblock;
            }
        }
        return;
    }

    inline void contract_baryon_pselpsel_no_projection_block(
        const std::string &job_tag, const int &traj,
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
        const Geometry &geo, const int &ti, const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy,
        std::vector<std::vector<Baryon_Block_No_Projection_Psel_Psel>> &block)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_no_projection_block");

        block.resize(psel_baryon_num_list[t_baryon_snk]);

        // #pragma omp parallel for
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

            block[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
        }

        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);
        }

#pragma omp parallel for collapse(2)
        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
            {
                long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);

                long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                Baryon_Block_No_Projection_Psel_Psel ppblock_np;
                const SpinMatrix &Cg5 = BaryonMatrixConstants::get_Cgamma5();
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);

                const WilsonMatrix &Cg5_wm = Cg5 * wm;
                const WilsonMatrix &wm_Cg5 = wm * Cg5;
                const WilsonMatrix &Cg5_wm1_Cg5 = Cg5_wm * Cg5;

                baryon_pi_no_projection_block_0_12(ppblock_np.block_0_12, Cg5_wm1_Cg5, wm);
                baryon_pi_no_projection_block_0_3(ppblock_np.block_0_3, Cg5_wm1_Cg5, wm);
                baryon_pi_no_projection_block_1_1(ppblock_np.block_1_1, wm_Cg5, Cg5_wm);
                baryon_pi_no_projection_block_1_2(ppblock_np.block_1_2, Cg5_wm1_Cg5, wm);
                baryon_pi_no_projection_block_1_3(ppblock_np.block_1_3, Cg5_wm1_Cg5, wm);

                ppblock_np.xg_snk_psel_idx = xg_snk_psel_idx;
                ppblock_np.xg_src_psel_idx = xg_src_psel_idx;
                ppblock_np.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                ppblock_np.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                ppblock_np.is_build = true;

                const int n_t = list_n_from_idx_baryon[xg_src_psel_idx];
                block[idx_baryon_snk][idx_baryon_src] = ppblock_np;
            }
        }
        return;
    }

    inline void ppblock_np_2_LT_block_np(Leading_Twist_Block_No_Projection_Psel_Psel &LT_block_np, const Baryon_Block_No_Projection_Psel_Psel &ppblock_np, const WilsonMatrix &g5_wm_y_src, const WilsonMatrix &wm_snk_y_g5);

    inline void contract_leading_twist_pselpsel_no_projection_block(
        const std::string &job_tag, const int &traj,
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
        const Geometry &geo, const int &ti, const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy,
        const std::vector<std::vector<Baryon_Block_No_Projection_Psel_Psel>> &block_baryon_no_proj,
        std::vector<std::vector<std::vector<Leading_Twist_Block_No_Projection_Psel_Psel>>> &block)
    {
        TIMER_VERBOSE("contract_leading_twist_pselpsel_no_projection_block");
        block.resize(psel_baryon_num_list[t_baryon_snk]);

        // #pragma omp parallel for
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            qassert(xg_snk[3] == t_baryon_snk);

            block[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                qassert(xg_src[3] == t_baryon_src);

                block[idx_baryon_snk][idx_baryon_src].resize(psel_meson_num_list[t_meson]);
            }
        }

#pragma omp parallel for collapse(2)
        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
            {
                const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                // const PselProp &prop_snk = get_psel_prop(job_tag, traj, xg_snk, is_meson_smear, is_baryon_smear);

                Baryon_Block_No_Projection_Psel_Psel ppblock_np = block_baryon_no_proj[idx_baryon_snk][idx_baryon_src];

                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                    const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                    const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                    Leading_Twist_Block_No_Projection_Psel_Psel LT_block_np;

                    const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

                    const WilsonMatrix g5_wm_y_src = gamma5 * prop_src.get_elem(xg_y_psel_idx);
                    const WilsonMatrix wm_snk_y_g5 = prop_y.get_elem(xg_snk_psel_idx) * gamma5;

                    ppblock_np_2_LT_block_np(LT_block_np, ppblock_np, g5_wm_y_src, wm_snk_y_g5);
                    LT_block_np.xg_snk_psel_idx = xg_snk_psel_idx;
                    LT_block_np.xg_src_psel_idx = xg_src_psel_idx;
                    LT_block_np.xg_y_psel_idx = xg_y_psel_idx;
                    LT_block_np.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                    LT_block_np.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                    LT_block_np.n_y = list_n_from_idx_meson[xg_y_psel_idx];
                    qassert(LT_block_np.is_build);
                    block[idx_baryon_snk][idx_baryon_src][idx_meson_ty] = LT_block_np;
                }
            }
        }
    }

    inline void ppblock_np_2_LT_block_np(
        Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra &LT_block_unex,
        Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra &LT_block_ex,
        const Baryon_Block_No_Projection_Psel_Psel &ppblock_np,
        const WilsonMatrix &g5_wm_y_src, const WilsonMatrix &wm_snk_y_g5);

    inline void contract_leading_twist_pselpsel_no_projection_block_ultra_sum(
        const std::string &job_tag, const int &traj,
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
        const Geometry &geo, const int &ti, const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src, const int &dtxy,
        const std::vector<std::vector<Baryon_Block_No_Projection_Psel_Psel>> &block_baryon_no_proj,
        std::vector<std::vector<Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra>> &block_unex,
        std::vector<std::vector<std::vector<Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra>>> &block_ex)
    {
        TIMER_VERBOSE("contract_leading_twist_pselpsel_no_projection_block_ultra_sum");

        // #pragma omp parallel for
        block_unex.resize(psel_baryon_num_list[t_baryon_snk]);
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            qassert(xg_snk[3] == t_baryon_snk);

            block_unex[idx_baryon_snk].resize(psel_meson_num_list[t_meson]);
        }

        block_ex.resize(psel_meson_num_list[t_meson]);
        for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
        {
            const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
            const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
            qassert(xg_y[3] == t_meson);

            block_ex[idx_meson_ty].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                block_ex[idx_meson_ty][idx_baryon_src].resize(psel_baryon_num_list[t_baryon_snk]);
            }
        }

#pragma omp parallel for collapse(2)
        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
            {
                const long &xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                const long &xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                // const PselProp &prop_snk = get_psel_prop(job_tag, traj, xg_snk, is_meson_smear, is_baryon_smear);

                Baryon_Block_No_Projection_Psel_Psel ppblock_np = block_baryon_no_proj[idx_baryon_snk][idx_baryon_src];

                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    const long &xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                    const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                    const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                    if (is_baryon_smear == is_meson_smear && (xg_snk_psel_idx == xg_y_psel_idx || xg_src_psel_idx == xg_y_psel_idx))
                    {
                        continue;
                    }

                    Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra &LT_block_unex = block_unex[idx_baryon_snk][idx_meson_ty];
                    Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra &LT_block_ex = block_ex[idx_meson_ty][idx_baryon_src][idx_baryon_snk];

                    const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

                    const WilsonMatrix g5_wm_y_src = gamma5 * prop_src.get_elem(xg_y_psel_idx);
                    const WilsonMatrix wm_snk_y_g5 = prop_y.get_elem(xg_snk_psel_idx) * gamma5;

                    ppblock_np_2_LT_block_np(LT_block_unex, LT_block_ex, ppblock_np, g5_wm_y_src, wm_snk_y_g5);
                    LT_block_unex.xg_snk_psel_idx = xg_snk_psel_idx;
                    LT_block_unex.xg_y_psel_idx = xg_y_psel_idx;
                    LT_block_unex.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                    LT_block_unex.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                    LT_block_unex.n_y = list_n_from_idx_meson[xg_y_psel_idx];

                    LT_block_ex.xg_src_psel_idx = xg_src_psel_idx;
                    LT_block_ex.xg_y_psel_idx = xg_y_psel_idx;
                    LT_block_ex.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                    LT_block_ex.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                    LT_block_ex.n_y = list_n_from_idx_meson[xg_y_psel_idx];

                    qassert(LT_block_unex.is_build);
                    qassert(LT_block_ex.is_build);
                }
            }
        }
    }

    inline void contract_baryon_psel_two_point_block(std::vector<std::vector<std::vector<Baryon_Block_Scalar_Psel_Psel>>> &block,
                                                     const std::vector<int> &psel_baryon_num_list, const std::vector<int> &list_n_from_idx_baryon,
                                                     const std::string &job_tag, const int &traj,
                                                     const Geometry &geo, const PointsSelection &psel_baryon,
                                                     const int &dtmax, const bool &is_baryon_smear)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_block");
        const Coordinate total_site = geo.total_site();

        const long n_points = psel_baryon.size();
        block.resize(n_points);

#pragma omp parallel for
        for (long n = 0; n < n_points; ++n)
        {
            block[n].resize(total_site[3]);
        }

#pragma omp parallel for
        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long &xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            const int t_baryon_snk = xg_snk[3];
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int del_t = dt;
                const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_baryon_num_list[t_baryon_src]);
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            const int t_baryon_src = xg_src[3];
            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);

#pragma omp parallel for
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long &xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                const int t_baryon_snk = xg_snk[3];
                const int del_t = mod(t_baryon_snk - t_baryon_src, total_site[3]);
                const int &dt = del_t;
                Baryon_Block_Scalar_Psel_Psel ppblock;
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);
                ppblock.block_pp[0] = proton_twop_block(wm, wm, wm, 0);
                ppblock.block_pp[1] = proton_twop_block(wm, wm, wm, 1);
                ppblock.xg_snk_psel_idx = xg_snk_psel_idx;
                ppblock.xg_src_psel_idx = xg_src_psel_idx;
                ppblock.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                ppblock.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                ppblock.is_build = true;

                const int n_t = list_n_from_idx_baryon[xg_src_psel_idx];
                block[xg_snk_psel_idx][dt][n_t] = ppblock;
            }
        }
        return;
    }

    inline void contract_baryon_psel_two_point_no_projection_block(std::vector<std::vector<std::vector<Baryon_Block_2PF_No_Projection_Psel_Psel>>> &block,
                                                                   const std::vector<int> &psel_baryon_num_list, const std::vector<int> &list_n_from_idx_baryon,
                                                                   const std::string &job_tag, const int &traj,
                                                                   const Geometry &geo, const PointsSelection &psel_baryon,
                                                                   const int &dtmax, const bool &is_baryon_smear)
    {
        TIMER_VERBOSE("contract_baryon_psel_two_point_no_projection_block");
        const Coordinate total_site = geo.total_site();

        const long n_points = psel_baryon.size();
        block.resize(n_points);

#pragma omp parallel for
        for (long n = 0; n < n_points; ++n)
        {
            block[n].resize(total_site[3]);
        }

#pragma omp parallel for
        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long &xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            const int t_baryon_snk = xg_snk[3];
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int del_t = dt;
                const int t_baryon_src = mod(t_baryon_snk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_baryon_num_list[t_baryon_src]);
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            const int t_baryon_src = xg_src[3];
            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);

#pragma omp parallel for
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long &xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                const int t_baryon_snk = xg_snk[3];
                const int del_t = mod(t_baryon_snk - t_baryon_src, total_site[3]);
                const int &dt = del_t;
                Baryon_Block_2PF_No_Projection_Psel_Psel ppblock;
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);
                ppblock.block_pp[0] = proton_twop_no_projection_block(wm, wm, wm, 0);
                ppblock.block_pp[1] = proton_twop_no_projection_block(wm, wm, wm, 1);
                ppblock.xg_snk_psel_idx = xg_snk_psel_idx;
                ppblock.xg_src_psel_idx = xg_src_psel_idx;
                ppblock.n_snk = list_n_from_idx_baryon[xg_snk_psel_idx];
                ppblock.n_src = list_n_from_idx_baryon[xg_src_psel_idx];
                ppblock.is_build = true;

                const int n_t = list_n_from_idx_baryon[xg_src_psel_idx];
                block[xg_snk_psel_idx][dt][n_t] = ppblock;
            }
        }
        return;
    }

    inline void contract_meson_psel_two_point_block(std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> &block,
                                                    const std::vector<int> &psel_meson_num_list, const std::vector<int> &list_n_from_idx_meson,
                                                    const std::string &job_tag, const int &traj,
                                                    const Geometry &geo, const PointsSelection &psel_meson,
                                                    const int &dtmax, const bool &is_meson_smear)
    {
        TIMER_VERBOSE("contract_meson_pselpsel_block");
        const Coordinate total_site = geo.total_site();
        // const int num_dt = dtmax - dtmin + 1;
        // qassert(num_dt > 0);
        const long n_points = psel_meson.size();
        block.resize(n_points);
#pragma omp parallel for
        for (long n = 0; n < n_points; ++n)
        {
            block[n].resize(total_site[3]);
        }

#pragma omp parallel for
        for (long nsnk = 0; nsnk < n_points; ++nsnk)
        {
            const long &xg_snk_psel_idx = nsnk;
            const Coordinate &xg_snk = psel_meson[xg_snk_psel_idx];
            const int t_meson_snk = xg_snk[3];
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int del_t = dt;
                const int t_meson_src = mod(t_meson_snk - del_t, total_site[3]);
                block[xg_snk_psel_idx][dt].resize(psel_meson_num_list[t_meson_src]);
            }
        }

        for (long nsrc = 0; nsrc < n_points; ++nsrc)
        {
            const long &xg_src_psel_idx = nsrc;
            const Coordinate &xg_src = psel_meson[xg_src_psel_idx];
            const int t_meson_src = xg_src[3];
            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_meson_smear);

#pragma omp parallel for
            for (long nsnk = 0; nsnk < n_points; ++nsnk)
            {
                const long &xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_meson[xg_snk_psel_idx];
                const int t_meson_snk = xg_snk[3];
                const int del_t = mod(t_meson_snk - t_meson_src, total_site[3]);
                const int &dt = del_t;

                Meson_Block_Scalar_Psel_Psel ppblock;
                const WilsonMatrix &wm = prop_src.get_elem(xg_snk_psel_idx);
                ppblock.block_pp = pion_twop_block(wm, wm);
                ppblock.xg_snk_psel_idx = xg_snk_psel_idx;
                ppblock.xg_src_psel_idx = xg_src_psel_idx;
                ppblock.n_snk = list_n_from_idx_meson[xg_snk_psel_idx];
                ppblock.n_src = list_n_from_idx_meson[xg_src_psel_idx];
                ppblock.is_build = true;

                const int n_t = list_n_from_idx_meson[xg_src_psel_idx];
                block[xg_snk_psel_idx][dt][n_t] = ppblock;
            }
        }
        return;
    }

    inline void contract_meson_psel_two_point_block_mate_ultra(std::vector<std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>>> &block,
                                                               const std::string &job_tag, const int &traj,
                                                               const bool &is_meson_smear,
                                                               const PointsSelection &psel_meson,
                                                               const long &psel_meson_num,
                                                               const std::vector<int> &psel_meson_num_list,
                                                               const std::vector<int> &list_n_from_idx_meson,
                                                               const std::vector<std::vector<long>> &idx_class_by_time_meson,
                                                               const std::vector<int> &fsel_num_list,
                                                               const std::vector<int> &list_n_from_idx_fsel,
                                                               const std::vector<std::vector<long>> &idx_class_by_time_fsel,
                                                               const FieldSelection &fsel,
                                                               const Geometry &geo,
                                                               const Coordinate &total_site,
                                                               const std::vector<std::vector<int>> &Momentum_Targets_curr)
    {
        TIMER_VERBOSE("contract_meson_psel_two_point_block_mate_ultra");
        const int &mom_num_curr = Momentum_Targets_curr.size();

        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        block.resize(psel_meson_num);
#pragma omp parallel for
        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            block[xg_src_psel_idx].resize(total_site[3]);
        }

#pragma omp parallel for
        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            const Coordinate &xg_meson_src = psel_meson[xg_src_psel_idx];
            const int t_meson_src = xg_meson_src[3];
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int &del_t = dt;
                const int t_meson_snk = mod(t_meson_src + del_t, total_site[3]);
                block[xg_src_psel_idx][dt].resize(mom_num_curr);
                for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                {
                    block[xg_src_psel_idx][dt][mom_curr_id].resize(8);
                }
            }
        }

        array<SpinMatrix, 8> g5_va_ms;
        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
        {
            g5_va_ms[VA_idx] = gamma5 * va_ms[VA_idx];
        }

        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            const Coordinate &xg_meson_src = psel_meson[xg_src_psel_idx];
            const int t_meson_src = xg_meson_src[3];
            const SelProp &prop_src = get_prop(job_tag, traj, xg_meson_src, is_meson_smear);

#pragma omp parallel for
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int &del_t = dt;
                const int t_meson_snk = mod(t_meson_src + del_t, total_site[3]);
                for (int idx_snk = 0; idx_snk < fsel_num_list[t_meson_snk]; idx_snk++)
                {
                    long xg_snk_fsel_id = idx_class_by_time_fsel[t_meson_snk][idx_snk];
                    const long index = fsel.indices[xg_snk_fsel_id];
                    const Coordinate xg_meson_snk_l = geo.coordinate_from_index(index);
                    const Coordinate xg_meson_snk_g = geo.coordinate_g_from_l(xg_meson_snk_l);
                    const int &t_meson_snk = xg_meson_snk_g[3];

                    qassert(list_n_from_idx_fsel[xg_snk_fsel_id] == idx_snk);

                    const WilsonMatrix &wm = prop_src.get_elem(xg_snk_fsel_id);
                    for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                    {
                        const WilsonMatrix g5_ja_wm = g5_va_ms[VA_idx] * wm;
                        ComplexD block_value_temp = pion_2pt_block(g5_ja_wm, wm);
                        for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                        {
                            long phase_temp = -space_dot(Momentum_Targets_curr[mom_curr_id], xg_meson_snk_g);
                            int phase = mod(phase_temp, total_site[0]);
                            Meson_Block_Scalar_Psel_Psel &ppblock = block[xg_src_psel_idx][dt][mom_curr_id][VA_idx];

                            // ppblock.xg_snk_psel_idx = xg_snk_fsel_id;
                            ppblock.xg_src_psel_idx = xg_src_psel_idx;
                            // ppblock.n_snk = list_n_from_idx_fsel[xg_snk_fsel_id];
                            ppblock.n_src = list_n_from_idx_meson[xg_src_psel_idx];
                            ppblock.is_build = true;

                            ppblock.block_pp += block_value_temp * exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]);
                        }
                    }
                }
                if (fsel_num_list[t_meson_snk] == 0)
                {
                    for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                    {
                        for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                        {
                            Meson_Block_Scalar_Psel_Psel &ppblock = block[xg_src_psel_idx][dt][mom_curr_id][VA_idx];

                            ppblock.xg_src_psel_idx = xg_src_psel_idx;
                            ppblock.n_src = list_n_from_idx_meson[xg_src_psel_idx];
                            ppblock.is_build = true;
                        }
                    }
                }
            }
        }

        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            const Coordinate &xg_meson_src = psel_meson[xg_src_psel_idx];
            const int t_meson_src = xg_meson_src[3];
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int &del_t = dt;
                const int t_meson_snk = mod(t_meson_src + del_t, total_site[3]);

                for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                {
                    for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                    {
                        Meson_Block_Scalar_Psel_Psel &ppblock = block[xg_src_psel_idx][dt][mom_curr_id][VA_idx];
                        ComplexD &block_value = ppblock.block_pp;
                        glb_sum(block_value);
                    }
                }
            }
        }
        return;
    }

    inline void contract_meson_psel_two_point_block_ultra(std::vector<std::vector<std::vector<Meson_Block_Scalar_Psel_Psel>>> &block,
                                                          const std::string &job_tag, const int &traj,
                                                          const bool &is_meson_smear,
                                                          const PointsSelection &psel_meson,
                                                          const long &psel_meson_num,
                                                          const std::vector<int> &psel_meson_num_list,
                                                          const std::vector<int> &list_n_from_idx_meson,
                                                          const std::vector<std::vector<long>> &idx_class_by_time_meson,
                                                          const std::vector<int> &fsel_num_list,
                                                          const std::vector<int> &list_n_from_idx_fsel,
                                                          const std::vector<std::vector<long>> &idx_class_by_time_fsel,
                                                          const FieldSelection &fsel,
                                                          const Geometry &geo,
                                                          const Coordinate &total_site)
    {
        TIMER_VERBOSE("contract_meson_psel_two_point_block_ultra");

        block.resize(psel_meson_num);
#pragma omp parallel for
        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            block[xg_src_psel_idx].resize(total_site[3]);
        }

#pragma omp parallel for
        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            const Coordinate &xg_meson_src = psel_meson[xg_src_psel_idx];
            const int t_meson_src = xg_meson_src[3];
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int &del_t = dt;
                const int t_meson_snk = mod(t_meson_src + del_t, total_site[3]);
                block[xg_src_psel_idx][dt].resize(fsel_num_list[t_meson_snk]);
            }
        }

        for (long xg_src_psel_idx = 0; xg_src_psel_idx < psel_meson_num; ++xg_src_psel_idx)
        {
            const Coordinate &xg_meson_src = psel_meson[xg_src_psel_idx];
            const int t_meson_src = xg_meson_src[3];
            const SelProp &prop_src = get_prop(job_tag, traj, xg_meson_src, is_meson_smear);

#pragma omp parallel for
            for (int dt = 0; dt < total_site[3]; ++dt)
            {
                const int &del_t = dt;
                const int t_meson_snk = mod(t_meson_src + del_t, total_site[3]);
                for (int idx_snk = 0; idx_snk < fsel_num_list[t_meson_snk]; idx_snk++)
                {
                    long xg_snk_fsel_id = idx_class_by_time_fsel[t_meson_snk][idx_snk];
                    const long index = fsel.indices[xg_snk_fsel_id];
                    const Coordinate xg_meson_snk_l = geo.coordinate_from_index(index);
                    const Coordinate xg_meson_snk_g = geo.coordinate_g_from_l(xg_meson_snk_l);
                    const int &t_meson_snk = xg_meson_snk_g[3];

                    qassert(list_n_from_idx_fsel[xg_snk_fsel_id] == idx_snk);

                    Meson_Block_Scalar_Psel_Psel ppblock;
                    const WilsonMatrix &wm = prop_src.get_elem(xg_snk_fsel_id);
                    ppblock.block_pp = pion_twop_block(wm, wm);
                    ppblock.xg_snk_psel_idx = xg_snk_fsel_id;
                    ppblock.xg_src_psel_idx = xg_src_psel_idx;
                    ppblock.n_snk = list_n_from_idx_fsel[xg_snk_fsel_id];
                    ppblock.n_src = list_n_from_idx_meson[xg_src_psel_idx];
                    ppblock.is_build = true;

                    block[xg_src_psel_idx][dt][idx_snk] = ppblock;
                }
            }
        }
        return;
    }

    inline void contract_baryon_pselpsel_y_block(const std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
                                                 std::vector<std::vector<Baryon_Block_Psel_Psel>> &ppblock,
                                                 const std::vector<int> &psel_baryon_num_list, const std::vector<int> &list_n_from_idx_baryon,
                                                 const std::string &job_tag, const int &traj,
                                                 const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
                                                 const long &xg_y_psel_idx, const Coordinate &xg_y,
                                                 const int &dtmin, const int tsep, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_y_block");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        qassert(num_dtxy == (int)ppblock.size());
        const int t_meson = xg_y[3];
        const long psel_baryon_num = psel_baryon.size();
        const long psel_meson_num = psel_meson.size();
        for (int dt = 0; dt < num_dtxy; dt++)
        {
            const int &dtxy = dt - num_dtxy / 2;
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
            const int num_psel_src = psel_baryon_num_list[t_baryon_src];
            const int num_psel_snk = psel_baryon_num_list[t_baryon_snk];
            ppblock[dt].resize(num_psel_snk);

            // #pragma omps parallel for
            for (long nsnk = 0; nsnk < psel_baryon_num; ++nsnk)
            {
                const long &xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                if (xg_snk[3] != t_baryon_snk)
                {
                    continue;
                }
                const int &dt_idx = mod(t_baryon_snk - t_baryon_src, total_site[3]) - dtmin;
                for (int n_t = 0; n_t < num_psel_src; ++n_t)
                {
                    const Baryon_Block_Psel_Psel &block1 = block[xg_snk_psel_idx][dt_idx][n_t];

                    // if (!block1.is_build)
                    // {
                    //     displayln(ssprintf("thread_id: %d/%d, tsep:%d, xg_snk_psel_idx:%d, t_baryon_snk:%d, xt:%d, t_meson:%d, t_baryon_src:%d, dt_idx:%d, n_t:%d,dt:%d", omp_get_thread_num(), omp_get_num_threads(), tsep, xg_snk_psel_idx, t_baryon_snk, xt, t_meson, t_baryon_src, dt_idx, n_t, dt));
                    // }

                    qassert(block1.is_build);
                    qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
                    const int i_snk = block1.n_snk;
                    const long &xg_src_psel_idx = block1.xg_src_psel_idx;
                    const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                    if (tsep == 0 && is_baryon_smear == is_meson_smear)
                    {
                        if (dtxy >= 0)
                        {
                            if (xg_src_psel_idx == xg_y_psel_idx)
                            {
                                continue;
                            }
                        }
                        else
                        {
                            if (xg_snk_psel_idx == xg_y_psel_idx)
                            {
                                continue;
                            }
                        }
                    }

                    const PselProp &prop_y_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
                    const WilsonMatrix &wm = prop_y_src.get_elem(xg_y_psel_idx);

                    for (int gram = 0; gram < (int)block1.block_pp.size(); ++gram)
                    {
                        ppblock[dt][i_snk].block_pp[gram] += wm * block1.block_pp[gram];
                    }

                    qassert(block1.n_snk == list_n_from_idx_baryon[xg_snk_psel_idx]);
                    ppblock[dt][i_snk].xg_snk_psel_idx = xg_snk_psel_idx;
                    ppblock[dt][i_snk].n_snk = i_snk;

                    ppblock[dt][i_snk].xg_src_psel_idx = xg_y_psel_idx;
                    ppblock[dt][i_snk].n_src = list_n_from_idx_baryon[xg_y_psel_idx]; // 加不加无所谓因为src位置在调用的时候已在上层确定

                    ppblock[dt][i_snk].is_build = true;
                }
            }
        }
    }

    inline void contract_baryon_pselpsel_y_block_direct_phase(const std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
                                                              std::vector<std::vector<Baryon_Block_Psel_Psel>> &ppblock,
                                                              const std::vector<int> &psel_baryon_num_list, const std::vector<int> &list_n_from_idx_baryon,
                                                              const std::string &job_tag, const int &traj,
                                                              const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
                                                              const long &xg_y_psel_idx, const Coordinate &xg_y,
                                                              const int &dtmin, const int tsep, const bool is_baryon_smear, const bool is_meson_smear, const std::vector<std::vector<int>> &Momentum_Targets_in)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_y_block_direct_phase");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        qassert(num_dtxy == (int)ppblock.size());

        qassert(Momentum_Targets_in.size() == 1);
        const int t_meson = xg_y[3];
        const long psel_baryon_num = psel_baryon.size();
        const long psel_meson_num = psel_meson.size();
        for (int dt = 0; dt < num_dtxy; dt++)
        {
            const int &dtxy = dt - num_dtxy / 2;
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
            const int num_psel_src = psel_baryon_num_list[t_baryon_src];
            const int num_psel_snk = psel_baryon_num_list[t_baryon_snk];
            ppblock[dt].resize(num_psel_snk);

            // #pragma omps parallel for
            for (long nsnk = 0; nsnk < psel_baryon_num; ++nsnk)
            {
                const long &xg_snk_psel_idx = nsnk;
                const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
                if (xg_snk[3] != t_baryon_snk)
                {
                    continue;
                }
                const int &dt_idx = mod(t_baryon_snk - t_baryon_src, total_site[3]) - dtmin;
                for (int n_t = 0; n_t < num_psel_src; ++n_t)
                {
                    const Baryon_Block_Psel_Psel &block1 = block[xg_snk_psel_idx][dt_idx][n_t];

                    // if (!block1.is_build)
                    // {
                    //     displayln(ssprintf("thread_id: %d/%d, tsep:%d, xg_snk_psel_idx:%d, t_baryon_snk:%d, xt:%d, t_meson:%d, t_baryon_src:%d, dt_idx:%d, n_t:%d,dt:%d", omp_get_thread_num(), omp_get_num_threads(), tsep, xg_snk_psel_idx, t_baryon_snk, xt, t_meson, t_baryon_src, dt_idx, n_t, dt));
                    // }

                    qassert(block1.is_build);
                    qassert(block1.xg_snk_psel_idx == xg_snk_psel_idx);
                    const int i_snk = block1.n_snk;
                    const long &xg_src_psel_idx = block1.xg_src_psel_idx;
                    const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];

                    if (tsep == 0 && is_baryon_smear == is_meson_smear)
                    {
                        if (dtxy >= 0)
                        {
                            if (xg_src_psel_idx == xg_y_psel_idx)
                            {
                                continue;
                            }
                        }
                        else
                        {
                            if (xg_snk_psel_idx == xg_y_psel_idx)
                            {
                                continue;
                            }
                        }
                    }

                    const PselProp &prop_y_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
                    const WilsonMatrix &wm = prop_y_src.get_elem(xg_y_psel_idx);

                    long phase_temp = space_dot(Momentum_Targets_in[0], xg_src);
                    int phase = mod(phase_temp, total_site[0]);

                    for (int gram = 0; gram < (int)block1.block_pp.size(); ++gram)
                    {
                        ppblock[dt][i_snk].block_pp[gram] += wm * block1.block_pp[gram] * exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]);
                    }

                    qassert(block1.n_snk == list_n_from_idx_baryon[xg_snk_psel_idx]);
                    ppblock[dt][i_snk].xg_snk_psel_idx = xg_snk_psel_idx;
                    ppblock[dt][i_snk].n_snk = i_snk;

                    ppblock[dt][i_snk].xg_src_psel_idx = xg_y_psel_idx;
                    ppblock[dt][i_snk].n_src = list_n_from_idx_baryon[xg_y_psel_idx]; // 加不加无所谓因为src位置在调用的时候已在上层确定

                    ppblock[dt][i_snk].is_build = true;
                }
            }
        }
    }

    inline void contract_baryon_pselpsel_y_block_coincide_corr(const std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
                                                               Baryon_Block_Psel_Psel &ppblock_cor,
                                                               const std::string &job_tag, const int &traj,
                                                               const Geometry &geo, const PointsSelection &psel,
                                                               const std::vector<int> &psel_num_list, const std::vector<int> &list_n_from_idx,
                                                               const long &xg_snk_psel_idx,
                                                               const long &xg_x_idx,
                                                               const long &xg_y_psel_idx, const Coordinate &xg_y, const int &tsep, const int &dtxy, const int &dtmin, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_y_block_coincide_corr");
        qassert(is_baryon_smear == is_meson_smear);

        const Coordinate &total_site = geo.total_site();
        const int &t_baryon_snk = psel[xg_snk_psel_idx][3];
        const int &xt = psel[xg_x_idx][3];
        const int &t_meson = xg_y[3];
        const int &t_baryon_src = xt;

        const int &dt_idx = mod(t_baryon_snk - t_baryon_src, total_site[3]) - dtmin;

        qassert(t_meson == t_baryon_snk);

        const Baryon_Block_Psel_Psel &block0 = block[xg_snk_psel_idx][dt_idx][list_n_from_idx[xg_x_idx]];

        qassert(block0.is_build);
        qassert(block0.xg_snk_psel_idx == xg_snk_psel_idx);
        qassert(block0.n_src == list_n_from_idx[xg_x_idx]);
        qassert(list_n_from_idx[block0.xg_src_psel_idx] == list_n_from_idx[xg_x_idx]);

        if (!(block0.xg_src_psel_idx == xg_x_idx))
        {
            displayln(ssprintf("%d, %d, %d, %d, %d, %d, %d, %d, %d", block0.xg_src_psel_idx, xg_x_idx, list_n_from_idx[block0.xg_src_psel_idx], list_n_from_idx[xg_x_idx], psel_num_list[xt], t_baryon_snk, xt, t_meson, psel[block0.xg_src_psel_idx][3]));
        }
        qassert(block0.xg_src_psel_idx == xg_x_idx);
        qassert(psel[block0.xg_src_psel_idx][3] == psel[xg_x_idx][3]);
        qassert(block0.n_snk == list_n_from_idx[xg_snk_psel_idx]);

        const int &i_snk = block0.n_snk;
        const long &xg_src_psel_idx = block0.xg_src_psel_idx;
        const Coordinate &xg_src = psel[xg_src_psel_idx];

        const PselProp &prop_y_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        const WilsonMatrix &wm = prop_y_src.get_elem(xg_y_psel_idx);

        for (int gram = 0; gram < (int)block0.block_pp.size(); ++gram)
        {
            ppblock_cor.block_pp[gram] = wm * block0.block_pp[gram];
        }

        ppblock_cor.xg_snk_psel_idx = xg_snk_psel_idx;
        ppblock_cor.n_snk = i_snk;

        ppblock_cor.xg_src_psel_idx = xg_y_psel_idx;
        ppblock_cor.n_src = list_n_from_idx[xg_y_psel_idx]; // 加不加无所谓因为src位置在调用的时候已在上层确定

        ppblock_cor.is_build = true;
    }

    inline void contract_baryon_pselpsel_y_block_coincide_corr_direct_phase(const std::vector<std::vector<std::vector<Baryon_Block_Psel_Psel>>> &block,
                                                                            Baryon_Block_Psel_Psel &ppblock_cor,
                                                                            const std::string &job_tag, const int &traj,
                                                                            const Geometry &geo, const PointsSelection &psel,
                                                                            const std::vector<int> &psel_num_list, const std::vector<int> &list_n_from_idx,
                                                                            const long &xg_snk_psel_idx,
                                                                            const long &xg_x_idx,
                                                                            const long &xg_y_psel_idx, const Coordinate &xg_y, const int &tsep, const int &dtxy, const int &dtmin, const bool is_baryon_smear, const bool is_meson_smear, const std::vector<std::vector<int>> &Momentum_Targets_in)
    {
        TIMER_VERBOSE("contract_baryon_pselpsel_y_block_coincide_corr_direct_phase");
        qassert(is_baryon_smear == is_meson_smear);

        const Coordinate &total_site = geo.total_site();
        const int &t_baryon_snk = psel[xg_snk_psel_idx][3];
        const int &xt = psel[xg_x_idx][3];
        const int &t_meson = xg_y[3];
        const int &t_baryon_src = xt;

        const int &dt_idx = mod(t_baryon_snk - t_baryon_src, total_site[3]) - dtmin;

        qassert(t_meson == t_baryon_snk);

        const Baryon_Block_Psel_Psel &block0 = block[xg_snk_psel_idx][dt_idx][list_n_from_idx[xg_x_idx]];

        qassert(block0.is_build);
        qassert(block0.xg_snk_psel_idx == xg_snk_psel_idx);
        qassert(block0.n_src == list_n_from_idx[xg_x_idx]);
        qassert(list_n_from_idx[block0.xg_src_psel_idx] == list_n_from_idx[xg_x_idx]);

        if (!(block0.xg_src_psel_idx == xg_x_idx))
        {
            displayln(ssprintf("%d, %d, %d, %d, %d, %d, %d, %d, %d", block0.xg_src_psel_idx, xg_x_idx, list_n_from_idx[block0.xg_src_psel_idx], list_n_from_idx[xg_x_idx], psel_num_list[xt], t_baryon_snk, xt, t_meson, psel[block0.xg_src_psel_idx][3]));
        }
        qassert(block0.xg_src_psel_idx == xg_x_idx);
        qassert(psel[block0.xg_src_psel_idx][3] == psel[xg_x_idx][3]);
        qassert(block0.n_snk == list_n_from_idx[xg_snk_psel_idx]);

        const int &i_snk = block0.n_snk;
        const long &xg_src_psel_idx = block0.xg_src_psel_idx;
        const Coordinate &xg_src = psel[xg_src_psel_idx];

        const PselProp &prop_y_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        const WilsonMatrix &wm = prop_y_src.get_elem(xg_y_psel_idx);

        long phase_temp = space_dot(Momentum_Targets_in[0], xg_src);
        int phase = mod(phase_temp, total_site[0]);
        for (int gram = 0; gram < (int)block0.block_pp.size(); ++gram)
        {
            ppblock_cor.block_pp[gram] = wm * block0.block_pp[gram] * exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]);
        }

        ppblock_cor.xg_snk_psel_idx = xg_snk_psel_idx;
        ppblock_cor.n_snk = i_snk;

        ppblock_cor.xg_src_psel_idx = xg_y_psel_idx;
        ppblock_cor.n_src = list_n_from_idx[xg_y_psel_idx]; // 加不加无所谓因为src位置在调用的时候已在上层确定

        ppblock_cor.is_build = true;
    }

    inline void get_all_wilson_matrix(std::vector<WilsonMatrix> &wm_snk_src_all,
                                      std::vector<WilsonMatrix> &wm_y_src_all,
                                      std::vector<WilsonMatrix> &wm_snk_y_all,
                                      const std::string &job_tag, const int &traj,
                                      const long &xg_y_psel_idx,
                                      const int t_baryon_src, const int t_baryon_snk, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
                                      const int num_src, const int num_snk,
                                      std::vector<long> &list_idx_from_isrc,
                                      std::vector<long> &list_idx_from_isnk, const bool is_baryon_smear, const bool is_meson_smear)
    {
        if (num_src == 0)
        {
            return;
        }
        if (num_snk == 0)
        {
            return;
        }

        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const long &psel_meson_num = psel_meson.size();
        const long &psel_baryon_num = psel_baryon.size();

        int i_snk = 0;
        int i_src = 0;
        wm_snk_src_all.resize(num_src * num_snk);
        wm_snk_y_all.resize(num_snk);
        wm_y_src_all.resize(num_src);
        list_idx_from_isnk.resize(num_snk);
        list_idx_from_isrc.resize(num_src);
        for (long n = 0; n < psel_baryon_num; ++n)
        {
            const long &xg_snk_psel_idx = n;
            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];
            if (xg_snk[3] == t_baryon_snk)
            {
                const PselProp &prop_y_snk = get_psel_prop(job_tag, traj, xg_snk, is_meson_smear, is_baryon_smear);
                const WilsonMatrix &wm_y_snk = prop_y_snk.get_elem(xg_y_psel_idx);
                wm_snk_y_all[i_snk] =
                    gamma5 * (WilsonMatrix)matrix_adjoint(wm_y_snk) * gamma5;
                list_idx_from_isnk[i_snk] = xg_snk_psel_idx;
                i_snk += 1;
            }
            const long &xg_src_psel_idx = n;
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            if (xg_src[3] == t_baryon_src)
            {
                const PselProp &prop_y_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
                const PselProp &prop_snk_src = get_psel_prop(job_tag, traj, xg_src, is_baryon_smear, is_baryon_smear);
                const WilsonMatrix &wm_y_src = prop_y_src.get_elem(xg_y_psel_idx);
                wm_y_src_all[i_src] = wm_y_src;
                list_idx_from_isrc[i_src] = xg_src_psel_idx;
                int i_snk1 = 0;
                for (long nsnk = 0; nsnk < psel_baryon_num; ++nsnk)
                {
                    const long &xg_snk_psel_idx1 = nsnk;
                    const Coordinate &xg_snk1 = psel_baryon[xg_snk_psel_idx1];
                    if (xg_snk1[3] == t_baryon_snk)
                    {
                        const WilsonMatrix &wm_snk_src = prop_snk_src.get_elem(xg_snk_psel_idx1);
                        const int i_src_snk = i_src + num_src * i_snk1;
                        wm_snk_src_all[i_src_snk] = wm_snk_src;
                        i_snk1 += 1;
                    }
                }
                i_src += 1;
            }
        }
        return;
    }

    inline void contract_baryon_sequential_meson_y_block(std::vector<std::vector<Baryon_Block_Psel_Psel>> &spblock,
                                                         const std::vector<int> &psel_baryon_num_list,
                                                         const std::vector<int> &list_n_from_idx_baryon,
                                                         const std::string &job_tag, const int &traj,
                                                         const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
                                                         const long &xg_y_psel_idx, const Coordinate &xg_y,
                                                         const int &dtmin, const int tsep, const bool is_baryon_smear, const bool is_meson_smear)
    {
        TIMER_VERBOSE("contract_baryon_sequential_meson_y_block");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        qassert(num_dtxy == (int)spblock.size());
        const int t_meson = xg_y[3];
        for (int dt = 0; dt < num_dtxy; dt++)
        {
            const int &dtxy = dt - num_dtxy / 2;
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
            const int num_psel_src = psel_baryon_num_list[t_baryon_src];
            const int num_psel_snk = psel_baryon_num_list[t_baryon_snk];
            const int num_psel_snk_src = num_psel_src * num_psel_snk;
            spblock[dt].resize(num_psel_snk_src);

            std::vector<WilsonMatrix> wm_snk_src_all;
            std::vector<WilsonMatrix> wm_y_src_all;
            std::vector<WilsonMatrix> wm_snk_y_all;
            std::vector<long> list_idx_from_isrc;
            std::vector<long> list_idx_from_isnk;
            get_all_wilson_matrix(wm_snk_src_all, wm_y_src_all, wm_snk_y_all, job_tag, traj,
                                  xg_y_psel_idx, t_baryon_src, t_baryon_snk, psel_baryon, psel_meson, num_psel_src,
                                  num_psel_snk, list_idx_from_isrc, list_idx_from_isnk, is_baryon_smear, is_meson_smear);

#pragma omp parallel for
            for (int i_src_snk = 0; i_src_snk < num_psel_snk_src; ++i_src_snk)
            {
                const int i_src = mod(i_src_snk, num_psel_src);
                const int i_snk = (i_src_snk - i_src) / num_psel_src;
                const long &xg_src_psel_idx = list_idx_from_isrc[i_src];
                const long &xg_snk_psel_idx = list_idx_from_isnk[i_snk];
                const WilsonMatrix &wm_snk_src = wm_snk_src_all[i_src_snk];
                const WilsonMatrix &wm_y_src = wm_y_src_all[i_src];
                const WilsonMatrix &wm_snk_y = wm_snk_y_all[i_snk];

                if (tsep == 0 && is_baryon_smear == is_meson_smear)
                {
                    if (dtxy >= 0)
                    {
                        if (xg_src_psel_idx == xg_y_psel_idx)
                        {
                            continue;
                        }
                    }
                    else
                    {
                        if (xg_snk_psel_idx == xg_y_psel_idx)
                        {
                            continue;
                        }
                    }
                }

                Baryon_Block_Psel_Psel block1;
                baryon_pi_block_sp(block1.block_pp, wm_snk_src, wm_y_src, wm_snk_y);
                block1.xg_snk_psel_idx = xg_snk_psel_idx;
                block1.xg_src_psel_idx = xg_src_psel_idx;
                block1.n_snk = i_snk;
                block1.n_src = i_src;
                block1.is_build = true;
                spblock[dt][i_src_snk] = block1;
            }
        }
    }

    inline void contract_baryon_sequential_meson_y_block_direct_phase(std::vector<std::vector<Baryon_Block_Psel_Psel>> &spblock,
                                                                      const std::vector<int> &psel_baryon_num_list,
                                                                      const std::vector<int> &list_n_from_idx_baryon,
                                                                      const std::string &job_tag, const int &traj,
                                                                      const Geometry &geo, const PointsSelection &psel_baryon, const PointsSelection &psel_meson,
                                                                      const long &xg_y_psel_idx, const Coordinate &xg_y,
                                                                      const int &dtmin, const int tsep, const bool is_baryon_smear, const bool is_meson_smear, const std::vector<std::vector<int>> &Momentum_Targets_in)
    {
        TIMER_VERBOSE("contract_baryon_sequential_meson_y_block_direct_phase");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_max_sep(total_site[3]);
        qassert(num_dtxy == (int)spblock.size());
        const int t_meson = xg_y[3];
        for (int dt = 0; dt < num_dtxy; dt++)
        {
            const int &dtxy = dt - num_dtxy / 2;
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
            const int num_psel_src = psel_baryon_num_list[t_baryon_src];
            const int num_psel_snk = psel_baryon_num_list[t_baryon_snk];
            const int num_psel_snk_src = num_psel_src * num_psel_snk;
            spblock[dt].resize(num_psel_snk_src);

            std::vector<WilsonMatrix> wm_snk_src_all;
            std::vector<WilsonMatrix> wm_y_src_all;
            std::vector<WilsonMatrix> wm_snk_y_all;
            std::vector<long> list_idx_from_isrc;
            std::vector<long> list_idx_from_isnk;
            get_all_wilson_matrix(wm_snk_src_all, wm_y_src_all, wm_snk_y_all, job_tag, traj,
                                  xg_y_psel_idx, t_baryon_src, t_baryon_snk, psel_baryon, psel_meson, num_psel_src,
                                  num_psel_snk, list_idx_from_isrc, list_idx_from_isnk, is_baryon_smear, is_meson_smear);

#pragma omp parallel for
            for (int i_src_snk = 0; i_src_snk < num_psel_snk_src; ++i_src_snk)
            {
                const int i_src = mod(i_src_snk, num_psel_src);
                const int i_snk = (i_src_snk - i_src) / num_psel_src;
                const long &xg_src_psel_idx = list_idx_from_isrc[i_src];
                const long &xg_snk_psel_idx = list_idx_from_isnk[i_snk];
                const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                const WilsonMatrix &wm_snk_src = wm_snk_src_all[i_src_snk];
                const WilsonMatrix &wm_y_src = wm_y_src_all[i_src];
                const WilsonMatrix &wm_snk_y = wm_snk_y_all[i_snk];

                if (tsep == 0 && is_baryon_smear == is_meson_smear)
                {
                    if (dtxy >= 0)
                    {
                        if (xg_src_psel_idx == xg_y_psel_idx)
                        {
                            continue;
                        }
                    }
                    else
                    {
                        if (xg_snk_psel_idx == xg_y_psel_idx)
                        {
                            continue;
                        }
                    }
                }

                long phase_temp = space_dot(Momentum_Targets_in[0], xg_src);
                int phase = mod(phase_temp, total_site[0]);

                qassert(spblock[dt][i_src_snk].block_pp.size() == 5);
                Baryon_Block_Psel_Psel block1;
                baryon_pi_block_sp(block1.block_pp, wm_snk_src, wm_y_src, wm_snk_y);

                for (size_t i = 0; i < block1.block_pp.size(); i++)
                {
                    block1.block_pp[i] = block1.block_pp[i] * exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]);
                }

                block1.xg_snk_psel_idx = xg_snk_psel_idx;
                block1.xg_src_psel_idx = xg_src_psel_idx;
                block1.n_snk = i_snk;
                block1.n_src = i_src;
                block1.is_build = true;
                spblock[dt][i_src_snk] = block1;
            }
        }
    }

    inline void contract_sequential_block_psel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>> &block,
        const std::vector<std::vector<int>> &Momentum_Targets,
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
        const Geometry &geo,
        const int &tslice_1, const int &tslice_2, const int &tslice_3,
        const bool &is_baryon_1, const bool &is_baryon_2, const bool &is_baryon_3)
    {
        TIMER_VERBOSE("contract_sequential_block_psel");
        const Coordinate &total_site = geo.total_site();

        const unsigned mom_num = Momentum_Targets.size();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const bool &is_par1_smear = is_baryon_1 ? is_baryon_smear : is_meson_smear;
        const bool &is_par2_smear = is_baryon_2 ? is_baryon_smear : is_meson_smear;
        const bool &is_par3_smear = is_baryon_3 ? is_baryon_smear : is_meson_smear;

        const PointsSelection &psel_par1 = is_baryon_1 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par2 = is_baryon_2 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par3 = is_baryon_3 ? psel_baryon : psel_meson;

        const long &psel_par1_num = is_baryon_1 ? psel_baryon_num : psel_meson_num;
        const long &psel_par2_num = is_baryon_2 ? psel_baryon_num : psel_meson_num;
        const long &psel_par3_num = is_baryon_3 ? psel_baryon_num : psel_meson_num;

        const std::vector<int> &psel_par1_num_list = is_baryon_1 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par2_num_list = is_baryon_2 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par3_num_list = is_baryon_3 ? psel_baryon_num_list : psel_meson_num_list;

        const std::vector<int> &list_n_from_idx_par1 = is_baryon_1 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par2 = is_baryon_2 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par3 = is_baryon_3 ? list_n_from_idx_baryon : list_n_from_idx_meson;

        const std::vector<std::vector<long>> &idx_class_by_time_par1 = is_baryon_1 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par2 = is_baryon_2 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par3 = is_baryon_3 ? idx_class_by_time_baryon : idx_class_by_time_meson;

        // for (long idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; ++idx_par1)
        // {
        //     long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
        //     const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
        //     const PselProp &prop_par1 = get_psel_prop(job_tag, traj, xg_par1, is_par2_smear, is_par1_smear);
        // }
        for (long idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; ++idx_par2)
        {
            long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
            const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
            const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
        }
        for (long idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; ++idx_par3)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            const PselProp &prop_par3 = get_psel_prop(job_tag, traj, xg_par3, is_par2_smear, is_par3_smear);
        }

        block.resize(psel_par1_num_list[tslice_1]);
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            block[idx_par1].resize(psel_par3_num_list[tslice_3]);
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                block[idx_par1][idx_par3].resize(psel_par2_num_list[tslice_2]);
            }
        }

#pragma omp parallel for collapse(3)
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                for (int idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; idx_par2++)
                {
                    long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
                    const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
                    qassert(tslice_3 == xg_par3[3]);
                    const PselProp &prop_par3 = get_psel_prop(job_tag, traj, xg_par3, is_par2_smear, is_par3_smear);

                    long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
                    const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
                    const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
                    qassert(tslice_2 == xg_par2[3]);

                    long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
                    const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
                    qassert(tslice_1 == xg_par1[3]);

                    const WilsonMatrix &wm_1_2 = prop_par2.get_elem(xg_par1_psel_idx);
                    const WilsonMatrix &wm_2_3 = prop_par3.get_elem(xg_par2_psel_idx);

                    WilsonMatrix swm_1_2_3 = wm_1_2 * gamma5 * wm_2_3;

                    Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3][idx_par2];

                    qassert(!block_temp.is_build);
                    block_temp.is_build = true;
                    block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                    block_temp.xg_par2_psel_idx = xg_par2_psel_idx;
                    block_temp.xg_par3_psel_idx = xg_par3_psel_idx;

                    block_temp.n_snk = idx_par1;
                    block_temp.n_src = idx_par3;
                    block_temp.block = swm_1_2_3;
                }
            }
        }
    }

       inline void contract_sequential_block_psel_ultra_sum(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<Sequential_Prop_1pt_Block>> &block,
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
        const Geometry &geo,
        const int &tslice_1, const int &tslice_2, const int &tslice_3,
        const bool &is_baryon_1, const bool &is_baryon_2, const bool &is_baryon_3)
    {
        TIMER_VERBOSE("contract_sequential_block_psel");
        const Coordinate &total_site = geo.total_site();

        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const bool &is_par1_smear = is_baryon_1 ? is_baryon_smear : is_meson_smear;
        const bool &is_par2_smear = is_baryon_2 ? is_baryon_smear : is_meson_smear;
        const bool &is_par3_smear = is_baryon_3 ? is_baryon_smear : is_meson_smear;

        const PointsSelection &psel_par1 = is_baryon_1 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par2 = is_baryon_2 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par3 = is_baryon_3 ? psel_baryon : psel_meson;

        const long &psel_par1_num = is_baryon_1 ? psel_baryon_num : psel_meson_num;
        const long &psel_par2_num = is_baryon_2 ? psel_baryon_num : psel_meson_num;
        const long &psel_par3_num = is_baryon_3 ? psel_baryon_num : psel_meson_num;

        const std::vector<int> &psel_par1_num_list = is_baryon_1 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par2_num_list = is_baryon_2 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par3_num_list = is_baryon_3 ? psel_baryon_num_list : psel_meson_num_list;

        const std::vector<int> &list_n_from_idx_par1 = is_baryon_1 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par2 = is_baryon_2 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par3 = is_baryon_3 ? list_n_from_idx_baryon : list_n_from_idx_meson;

        const std::vector<std::vector<long>> &idx_class_by_time_par1 = is_baryon_1 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par2 = is_baryon_2 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par3 = is_baryon_3 ? idx_class_by_time_baryon : idx_class_by_time_meson;

        // for (long idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; ++idx_par1)
        // {
        //     long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
        //     const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
        //     const PselProp &prop_par1 = get_psel_prop(job_tag, traj, xg_par1, is_par2_smear, is_par1_smear);
        // }
        for (long idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; ++idx_par2)
        {
            long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
            const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
            const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
        }
        for (long idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; ++idx_par3)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            const PselProp &prop_par3 = get_psel_prop(job_tag, traj, xg_par3, is_par2_smear, is_par3_smear);
        }

        block.resize(psel_par1_num_list[tslice_1]);
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            block[idx_par1].resize(psel_par3_num_list[tslice_3]);
        }

#pragma omp parallel for collapse(2)
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                for (int idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; idx_par2++)
                {
                    long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
                    const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
                    qassert(tslice_3 == xg_par3[3]);
                    const PselProp &prop_par3 = get_psel_prop(job_tag, traj, xg_par3, is_par2_smear, is_par3_smear);

                    long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
                    const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
                    const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
                    qassert(tslice_2 == xg_par2[3]);

                    long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
                    const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
                    qassert(tslice_1 == xg_par1[3]);

                    const WilsonMatrix &wm_1_2 = prop_par2.get_elem(xg_par1_psel_idx);
                    const WilsonMatrix &wm_2_3 = prop_par3.get_elem(xg_par2_psel_idx);


                    Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3];

                    qassert(!block_temp.is_build);
                    block_temp.is_build = true;
                    block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                    // block_temp.xg_par2_psel_idx = xg_par2_psel_idx;
                    block_temp.xg_par3_psel_idx = xg_par3_psel_idx;

                    block_temp.n_snk = idx_par1;
                    block_temp.n_src = idx_par3;

                    block_temp.block += wm_1_2 * gamma5 * wm_2_3;
                }
            }
        }
    }

    inline void contract_sequential_block_psel_flike(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>> &block,
        const std::vector<std::vector<int>> &Momentum_Targets,
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
        const Geometry &geo,
        const int &tslice_1, const int &tslice_2, const int &tslice_3,
        const bool &is_baryon_1, const bool &is_baryon_2, const bool &is_baryon_3)
    {
        TIMER_VERBOSE("contract_sequential_block_psel_flike");
        const Coordinate &total_site = geo.total_site();

        const unsigned mom_num = Momentum_Targets.size();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const bool &is_par1_smear = is_baryon_1 ? is_baryon_smear : is_meson_smear;
        const bool &is_par2_smear = is_baryon_2 ? is_baryon_smear : is_meson_smear;
        const bool &is_par3_smear = is_baryon_3 ? is_baryon_smear : is_meson_smear;

        const PointsSelection &psel_par1 = is_baryon_1 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par2 = is_baryon_2 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par3 = is_baryon_3 ? psel_baryon : psel_meson;

        const long &psel_par1_num = is_baryon_1 ? psel_baryon_num : psel_meson_num;
        const long &psel_par2_num = is_baryon_2 ? psel_baryon_num : psel_meson_num;
        const long &psel_par3_num = is_baryon_3 ? psel_baryon_num : psel_meson_num;

        const std::vector<int> &psel_par1_num_list = is_baryon_1 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par2_num_list = is_baryon_2 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par3_num_list = is_baryon_3 ? psel_baryon_num_list : psel_meson_num_list;

        const std::vector<int> &list_n_from_idx_par1 = is_baryon_1 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par2 = is_baryon_2 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par3 = is_baryon_3 ? list_n_from_idx_baryon : list_n_from_idx_meson;

        const std::vector<std::vector<long>> &idx_class_by_time_par1 = is_baryon_1 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par2 = is_baryon_2 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par3 = is_baryon_3 ? idx_class_by_time_baryon : idx_class_by_time_meson;

        for (long idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; ++idx_par1)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
            const PselProp &prop_par1 = get_psel_prop(job_tag, traj, xg_par1, is_par2_smear, is_par1_smear);
        }
        // for (long idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; ++idx_par2)
        // {
        //     long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
        //     const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
        //     const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
        // }
        for (long idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; ++idx_par3)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            const PselProp &prop_par3 = get_psel_prop(job_tag, traj, xg_par3, is_par2_smear, is_par3_smear);
        }

        block.resize(psel_par1_num_list[tslice_1]);
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            block[idx_par1].resize(psel_par3_num_list[tslice_3]);
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                block[idx_par1][idx_par3].resize(psel_par2_num_list[tslice_2]);
            }
        }

#pragma omp parallel for collapse(3)
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                for (int idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; idx_par2++)
                {
                    long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
                    const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
                    qassert(tslice_3 == xg_par3[3]);
                    qassert(list_n_from_idx_par3[xg_par3_psel_idx] == idx_par3);
                    const PselProp &prop_par3 = get_psel_prop(job_tag, traj, xg_par3, is_par2_smear, is_par3_smear);

                    long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
                    const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
                    qassert(tslice_1 == xg_par1[3]);
                    qassert(list_n_from_idx_par1[xg_par1_psel_idx] == idx_par1);
                    const PselProp &prop_par1 = get_psel_prop(job_tag, traj, xg_par1, is_par2_smear, is_par1_smear);

                    long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
                    const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
                    qassert(list_n_from_idx_par2[xg_par2_psel_idx] == idx_par2);

                    const WilsonMatrix &wm_2_1 = prop_par1.get_elem(xg_par2_psel_idx);
                    const WilsonMatrix &wm_2_3 = prop_par3.get_elem(xg_par2_psel_idx);

                    WilsonMatrix swm_1_2_3 = gamma5 * (WilsonMatrix)matrix_adjoint(wm_2_1) * wm_2_3;

                    Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3][idx_par2];

                    qassert(!block_temp.is_build);
                    block_temp.is_build = true;
                    block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                    block_temp.xg_par2_psel_idx = xg_par2_psel_idx;
                    block_temp.xg_par3_psel_idx = xg_par3_psel_idx;

                    block_temp.n_snk = idx_par1;
                    block_temp.n_src = idx_par3;
                    block_temp.block = swm_1_2_3;
                }
            }
        }
    }

    inline void contract_sequential_block_x_psel_y_psel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_2pt_Block>>>> &block_x_y,
        const std::vector<std::vector<int>> &Momentum_Targets,
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
        const Geometry &geo,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src)
    {
        TIMER_VERBOSE("contract_sequential_block_x_psel_y_psel");
        const Coordinate &total_site = geo.total_site();

        const unsigned mom_num = Momentum_Targets.size();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            qassert(t_baryon_src == xg_src[3]);

            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        }

        block_x_y.resize(psel_baryon_num_list[t_baryon_snk]);
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            block_x_y[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                block_x_y[idx_baryon_snk][idx_baryon_src].resize(psel_meson_num_list[xt]);
                for (int idx_x = 0; idx_x < psel_meson_num_list[xt]; idx_x++)
                {
                    block_x_y[idx_baryon_snk][idx_baryon_src][idx_x].resize(psel_meson_num_list[t_meson]);
                }
            }
        }

        std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>> block_x4snk_y;
        std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>> block_x4y_src;

        contract_sequential_block_psel_flike(
            job_tag, traj, block_x4snk_y,
            Momentum_Targets,
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
            geo,
            t_baryon_snk, xt, t_meson,
            true, false, false);

        contract_sequential_block_psel_flike(
            job_tag, traj, block_x4y_src,
            Momentum_Targets,
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
            geo,
            t_meson, xt, t_baryon_src,
            false, false, true);

#pragma omp parallel for collapse(2)
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    for (int idx_x = 0; idx_x < psel_meson_num_list[xt]; idx_x++)
                    {
                        long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                        const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                        long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                        long xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                        const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                        const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                        long xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                        const Coordinate &xg_x = psel_meson[xg_x_psel_idx];

                        qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                        qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);
                        qassert(list_n_from_idx_meson[xg_x_psel_idx] == idx_x);

                        Sequential_Prop_1pt_Block &block_x4snk_y_temp = block_x4snk_y[idx_baryon_snk][idx_meson_ty][idx_x];
                        Sequential_Prop_1pt_Block &block_x4_y_src_temp = block_x4y_src[idx_meson_ty][idx_baryon_src][idx_x];

                        qassert(block_x4snk_y_temp.is_build);
                        qassert(block_x4snk_y_temp.xg_par1_psel_idx == xg_snk_psel_idx);
                        qassert(block_x4snk_y_temp.xg_par2_psel_idx == xg_x_psel_idx);
                        qassert(block_x4snk_y_temp.xg_par3_psel_idx == xg_y_psel_idx);

                        qassert(block_x4_y_src_temp.is_build);
                        qassert(block_x4_y_src_temp.xg_par1_psel_idx == xg_y_psel_idx);
                        qassert(block_x4_y_src_temp.xg_par2_psel_idx == xg_x_psel_idx);
                        qassert(block_x4_y_src_temp.xg_par3_psel_idx == xg_src_psel_idx);

                        const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                        const WilsonMatrix &wm_snk_y = prop_y.get_elem(xg_snk_psel_idx);

                        qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                        WilsonMatrix wm_sequence_0 = block_x4snk_y_temp.block * gamma5 * wm_y_src;
                        WilsonMatrix wm_sequence_1 = wm_snk_y * gamma5 * block_x4_y_src_temp.block;

                        Sequential_Prop_2pt_Block &block_x_y_temp = block_x_y[idx_baryon_snk][idx_baryon_src][idx_x][idx_meson_ty];

                        qassert(!block_x_y_temp.is_build);
                        block_x_y_temp.is_build = true;
                        block_x_y_temp.xg_snk_psel_idx = xg_snk_psel_idx;
                        block_x_y_temp.xg_src_psel_idx = xg_src_psel_idx;
                        block_x_y_temp.xg_x_psel_idx = xg_x_psel_idx;
                        block_x_y_temp.xg_y_psel_idx = xg_y_psel_idx;

                        block_x_y_temp.n_snk = idx_baryon_snk;
                        block_x_y_temp.n_src = idx_baryon_src;
                        block_x_y_temp.block[0] = wm_sequence_0;
                        block_x_y_temp.block[1] = wm_sequence_1;
                    }
                }
            }
        }
    }

    inline void contract_sequential_block_fsel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const int &tslice_1, const int &tslice_2, const int &tslice_3,
        const bool &is_baryon_1, const bool &is_baryon_3)
    {
        TIMER_VERBOSE("contract_sequential_block_fsel");

        const bool is_par2_smear = false;

        const Coordinate &total_site = geo.total_site();

        // const unsigned mom_num = Momentum_Targets.size();
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const bool &is_par1_smear = is_baryon_1 ? is_baryon_smear : is_meson_smear;
        const bool &is_par3_smear = is_baryon_3 ? is_baryon_smear : is_meson_smear;

        const PointsSelection &psel_par1 = is_baryon_1 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par3 = is_baryon_3 ? psel_baryon : psel_meson;

        const long &psel_par1_num = is_baryon_1 ? psel_baryon_num : psel_meson_num;
        const long &psel_par3_num = is_baryon_3 ? psel_baryon_num : psel_meson_num;

        const std::vector<int> &psel_par1_num_list = is_baryon_1 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par3_num_list = is_baryon_3 ? psel_baryon_num_list : psel_meson_num_list;

        const std::vector<int> &list_n_from_idx_par1 = is_baryon_1 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par3 = is_baryon_3 ? list_n_from_idx_baryon : list_n_from_idx_meson;

        const std::vector<std::vector<long>> &idx_class_by_time_par1 = is_baryon_1 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par3 = is_baryon_3 ? idx_class_by_time_baryon : idx_class_by_time_meson;

        /* -------------------------------- */
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();
        /* -------------------------------- */

        for (long idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; ++idx_par1)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];

            qassert(tslice_1 == xg_par1[3]);
            // qassert(xg_par1 == psel_smear[xg_par1_psel_idx]);

            const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);
        }
        // for (long idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; ++idx_par2)
        // {
        //     long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
        //     const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
        //     const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
        // }
        for (long idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; ++idx_par3)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);
        }

        block.resize(psel_par1_num_list[tslice_1]);
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            block[idx_par1].resize(psel_par3_num_list[tslice_3]);
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                block[idx_par1][idx_par3].resize(fsel_num_list[tslice_2]);
                for (int idx_par2 = 0; idx_par2 < fsel_num_list[tslice_2]; idx_par2++)
                {
                    block[idx_par1][idx_par3][idx_par2].resize(8);
                }
            }
        }

        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
            qassert(tslice_1 == xg_par1[3]);
            qassert(list_n_from_idx_par1[xg_par1_psel_idx] == idx_par1);

            const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);

            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
                const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
                qassert(tslice_3 == xg_par3[3]);
                qassert(list_n_from_idx_par3[xg_par3_psel_idx] == idx_par3);

                const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);

#pragma omp parallel for
                for (int idx_par2 = 0; idx_par2 < fsel_num_list[tslice_2]; idx_par2++)
                {
                    const long &xg_par2_fsel_idx = idx_class_by_time_fsel[tslice_2][idx_par2];
                    const long index = fsel.indices[xg_par2_fsel_idx];
                    const Coordinate xg_par2_l = geo.coordinate_from_index(index);
                    const Coordinate xg_par2_g = geo.coordinate_g_from_l(xg_par2_l);
                    const Coordinate &xg_par2 = xg_par2_g;
                    const int &t_par2 = xg_par2[3];
                    qassert(t_par2 == tslice_2);
                    qassert(idx_par2 == list_n_from_idx_fsel[xg_par2_fsel_idx]);

                    const WilsonMatrix &wm_2_1 = prop_par1.get_elem(xg_par2_fsel_idx);
                    const WilsonMatrix &wm_2_3 = prop_par3.get_elem(xg_par2_fsel_idx);

                    for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                    {
                        WilsonMatrix swm_1_2_3 = gamma5 * (WilsonMatrix)matrix_adjoint(wm_2_1) * gamma5 * va_ms[VA_idx] * wm_2_3;

                        Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3][idx_par2][VA_idx];

                        qassert(!block_temp.is_build);
                        block_temp.is_build = true;
                        block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                        block_temp.xg_par2_psel_idx = xg_par2_fsel_idx;
                        block_temp.xg_par3_psel_idx = xg_par3_psel_idx;

                        block_temp.n_snk = idx_par1;
                        block_temp.n_src = idx_par3;
                        block_temp.block = swm_1_2_3;
                    }
                }
            }
        }
    }

    inline void contract_sequential_block_fsel_ultra_sum(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const int &tslice_1, const int &tslice_2, const int &tslice_3,
        const bool &is_baryon_1, const bool &is_baryon_3)
    {
        TIMER_VERBOSE("contract_sequential_block_fsel_ultra_sum");

        const bool is_par2_smear = false;

        const Coordinate &total_site = geo.total_site();

        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const bool &is_par1_smear = is_baryon_1 ? is_baryon_smear : is_meson_smear;
        const bool &is_par3_smear = is_baryon_3 ? is_baryon_smear : is_meson_smear;

        const PointsSelection &psel_par1 = is_baryon_1 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par3 = is_baryon_3 ? psel_baryon : psel_meson;

        const long &psel_par1_num = is_baryon_1 ? psel_baryon_num : psel_meson_num;
        const long &psel_par3_num = is_baryon_3 ? psel_baryon_num : psel_meson_num;

        const std::vector<int> &psel_par1_num_list = is_baryon_1 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par3_num_list = is_baryon_3 ? psel_baryon_num_list : psel_meson_num_list;

        const std::vector<int> &list_n_from_idx_par1 = is_baryon_1 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par3 = is_baryon_3 ? list_n_from_idx_baryon : list_n_from_idx_meson;

        const std::vector<std::vector<long>> &idx_class_by_time_par1 = is_baryon_1 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par3 = is_baryon_3 ? idx_class_by_time_baryon : idx_class_by_time_meson;

        /* -------------------------------- */
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();
        /* -------------------------------- */

        for (long idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; ++idx_par1)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];

            qassert(tslice_1 == xg_par1[3]);
            // qassert(xg_par1 == psel_smear[xg_par1_psel_idx]);

            const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);
        }
        // for (long idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; ++idx_par2)
        // {
        //     long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
        //     const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
        //     const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
        // }
        for (long idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; ++idx_par3)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);
        }

        block.resize(psel_par1_num_list[tslice_1]);
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            block[idx_par1].resize(psel_par3_num_list[tslice_3]);
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                block[idx_par1][idx_par3].resize(1);
                for (int idx_par2 = 0; idx_par2 < (0 == fsel_num_list[tslice_2] ? 0 : 1); idx_par2++)
                {
                    block[idx_par1][idx_par3][idx_par2].resize(8);
                }
            }
        }

        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
            qassert(tslice_1 == xg_par1[3]);
            qassert(list_n_from_idx_par1[xg_par1_psel_idx] == idx_par1);

            const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);
        }

        for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            qassert(tslice_3 == xg_par3[3]);
            qassert(list_n_from_idx_par3[xg_par3_psel_idx] == idx_par3);

            const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);
        }

#pragma omp parallel for collapse(2)
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
                const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
                qassert(tslice_1 == xg_par1[3]);
                qassert(list_n_from_idx_par1[xg_par1_psel_idx] == idx_par1);

                const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);

                long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
                const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
                qassert(tslice_3 == xg_par3[3]);
                qassert(list_n_from_idx_par3[xg_par3_psel_idx] == idx_par3);

                const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);

                for (int idx_par2 = 0; idx_par2 < fsel_num_list[tslice_2]; idx_par2++)
                {
                    const long &xg_par2_fsel_idx = idx_class_by_time_fsel[tslice_2][idx_par2];
                    const long index = fsel.indices[xg_par2_fsel_idx];
                    const Coordinate xg_par2_l = geo.coordinate_from_index(index);
                    const Coordinate xg_par2_g = geo.coordinate_g_from_l(xg_par2_l);
                    const Coordinate &xg_par2 = xg_par2_g;
                    const int &t_par2 = xg_par2[3];
                    qassert(t_par2 == tslice_2);
                    qassert(idx_par2 == list_n_from_idx_fsel[xg_par2_fsel_idx]);

                    const WilsonMatrix &wm_2_1 = prop_par1.get_elem(xg_par2_fsel_idx);
                    const WilsonMatrix &wm_2_3 = prop_par3.get_elem(xg_par2_fsel_idx);

                    for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                    {
                        WilsonMatrix swm_1_2_3 = gamma5 * (WilsonMatrix)matrix_adjoint(wm_2_1) * gamma5 * va_ms[VA_idx] * wm_2_3;

                        Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3][0][VA_idx];

                        // qassert(!block_temp.is_build);
                        block_temp.is_build = true;
                        block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                        // block_temp.xg_par2_psel_idx = xg_par2_fsel_idx;
                        block_temp.xg_par3_psel_idx = xg_par3_psel_idx;

                        block_temp.n_snk = idx_par1;
                        block_temp.n_src = idx_par3;
                        block_temp.block += swm_1_2_3;
                    }
                }
            }
        }
    }

    inline void contract_sequential_block_fsel_mate_ultra_sum(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const int &tslice_1, const int &tslice_2, const int &tslice_3,
        const bool &is_baryon_1, const bool &is_baryon_3,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        bool need_glb_sum = true)
    {
        TIMER_VERBOSE("contract_sequential_block_fsel_mate_ultra_sum");

        const bool is_par2_smear = false;
        const int &mom_num_curr = Momentum_Targets_curr.size();

        const Coordinate &total_site = geo.total_site();

        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        const bool &is_par1_smear = is_baryon_1 ? is_baryon_smear : is_meson_smear;
        const bool &is_par3_smear = is_baryon_3 ? is_baryon_smear : is_meson_smear;

        const PointsSelection &psel_par1 = is_baryon_1 ? psel_baryon : psel_meson;
        const PointsSelection &psel_par3 = is_baryon_3 ? psel_baryon : psel_meson;

        const long &psel_par1_num = is_baryon_1 ? psel_baryon_num : psel_meson_num;
        const long &psel_par3_num = is_baryon_3 ? psel_baryon_num : psel_meson_num;

        const std::vector<int> &psel_par1_num_list = is_baryon_1 ? psel_baryon_num_list : psel_meson_num_list;
        const std::vector<int> &psel_par3_num_list = is_baryon_3 ? psel_baryon_num_list : psel_meson_num_list;

        const std::vector<int> &list_n_from_idx_par1 = is_baryon_1 ? list_n_from_idx_baryon : list_n_from_idx_meson;
        const std::vector<int> &list_n_from_idx_par3 = is_baryon_3 ? list_n_from_idx_baryon : list_n_from_idx_meson;

        const std::vector<std::vector<long>> &idx_class_by_time_par1 = is_baryon_1 ? idx_class_by_time_baryon : idx_class_by_time_meson;
        const std::vector<std::vector<long>> &idx_class_by_time_par3 = is_baryon_3 ? idx_class_by_time_baryon : idx_class_by_time_meson;

        /* -------------------------------- */
        const PointsSelection &psel = get_point_selection(job_tag, traj);
        const long n_points = psel.size();
        const PointsSelection &psel_smear = get_point_selection_smear(job_tag, traj);
        const long n_points_smear = psel_smear.size();
        /* -------------------------------- */

        for (long idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; ++idx_par1)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];

            qassert(tslice_1 == xg_par1[3]);
            // qassert(xg_par1 == psel_smear[xg_par1_psel_idx]);

            const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);
        }
        // for (long idx_par2 = 0; idx_par2 < psel_par2_num_list[tslice_2]; ++idx_par2)
        // {
        //     long xg_par2_psel_idx = idx_class_by_time_par2[tslice_2][idx_par2];
        //     const Coordinate &xg_par2 = psel_par2[xg_par2_psel_idx];
        //     const PselProp &prop_par2 = get_psel_prop(job_tag, traj, xg_par2, is_par1_smear, is_par2_smear);
        // }
        for (long idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; ++idx_par3)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);
        }

        block.resize(psel_par1_num_list[tslice_1]);
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            block[idx_par1].resize(psel_par3_num_list[tslice_3]);
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                block[idx_par1][idx_par3].resize(mom_num_curr);
                for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                {
                    block[idx_par1][idx_par3][mom_curr_id].resize(8);
                }
            }
        }

        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
            const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
            qassert(tslice_1 == xg_par1[3]);
            qassert(list_n_from_idx_par1[xg_par1_psel_idx] == idx_par1);

            const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);
        }

        for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
        {
            long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
            const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
            qassert(tslice_3 == xg_par3[3]);
            qassert(list_n_from_idx_par3[xg_par3_psel_idx] == idx_par3);

            const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);
        }

#pragma omp parallel for collapse(2)
        for (int idx_par1 = 0; idx_par1 < psel_par1_num_list[tslice_1]; idx_par1++)
        {
            for (int idx_par3 = 0; idx_par3 < psel_par3_num_list[tslice_3]; idx_par3++)
            {
                long xg_par1_psel_idx = idx_class_by_time_par1[tslice_1][idx_par1];
                const Coordinate &xg_par1 = psel_par1[xg_par1_psel_idx];
                qassert(tslice_1 == xg_par1[3]);
                qassert(list_n_from_idx_par1[xg_par1_psel_idx] == idx_par1);

                const SelProp &prop_par1 = get_prop(job_tag, traj, xg_par1, is_par1_smear);

                long xg_par3_psel_idx = idx_class_by_time_par3[tslice_3][idx_par3];
                const Coordinate &xg_par3 = psel_par3[xg_par3_psel_idx];
                qassert(tslice_3 == xg_par3[3]);
                qassert(list_n_from_idx_par3[xg_par3_psel_idx] == idx_par3);

                const SelProp &prop_par3 = get_prop(job_tag, traj, xg_par3, is_par3_smear);

                // for (int idx_par2 = 0; idx_par2 < (fsel_num_list[tslice_2] == 0 ? 0 : 1); idx_par2++)
                for (int idx_par2 = 0; idx_par2 < fsel_num_list[tslice_2]; idx_par2++)
                {
                    const long &xg_par2_fsel_idx = idx_class_by_time_fsel[tslice_2][idx_par2];
                    const long index = fsel.indices[xg_par2_fsel_idx];
                    const Coordinate xg_par2_l = geo.coordinate_from_index(index);
                    const Coordinate xg_par2_g = geo.coordinate_g_from_l(xg_par2_l);
                    const Coordinate &xg_par2 = xg_par2_g;
                    const int &t_par2 = xg_par2[3];
                    qassert(t_par2 == tslice_2);
                    qassert(idx_par2 == list_n_from_idx_fsel[xg_par2_fsel_idx]);

                    const WilsonMatrix &wm_2_1 = prop_par1.get_elem(xg_par2_fsel_idx);
                    const WilsonMatrix &wm_2_3 = prop_par3.get_elem(xg_par2_fsel_idx);

                    for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                    {
                        WilsonMatrix swm_1_2_3 = gamma5 * (WilsonMatrix)matrix_adjoint(wm_2_1) * gamma5 * va_ms[VA_idx] * wm_2_3;
                        for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                        {

                            Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3][mom_curr_id][VA_idx];

                            // qassert(!block_temp.is_build);
                            block_temp.is_build = true;
                            block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                            // block_temp.xg_par2_psel_idx = xg_par2_fsel_idx;
                            block_temp.xg_par3_psel_idx = xg_par3_psel_idx;
                            block_temp.n_snk = idx_par1;
                            block_temp.n_src = idx_par3;

                            long phase_temp = -space_dot(Momentum_Targets_curr[mom_curr_id], xg_par2);
                            int phase = mod(phase_temp, total_site[0]);
                            block_temp.block += swm_1_2_3 * exp(ii * 2.0 * PI * (double)phase / (double)total_site[0]);
                        }
                    }
                }

                if (need_glb_sum && (fsel_num_list[tslice_2] == 0))
                {
                    for (int mom_curr_id = 0; mom_curr_id < mom_num_curr; ++mom_curr_id)
                    {
                        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                        {
                            Sequential_Prop_1pt_Block &block_temp = block[idx_par1][idx_par3][mom_curr_id][VA_idx];

                            // qassert(!block_temp.is_build);
                            block_temp.is_build = true;
                            block_temp.xg_par1_psel_idx = xg_par1_psel_idx;
                            // block_temp.xg_par2_psel_idx = xg_par2_fsel_idx;
                            block_temp.xg_par3_psel_idx = xg_par3_psel_idx;
                            block_temp.n_snk = idx_par1;
                            block_temp.n_src = idx_par3;
                        }
                    }
                }
            }
        }
    }

    struct Sequential_Prop_1pt_Block_Collection
    {
        bool is_build;
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> block_coll;
        int tslice_1;
        int tslice_2;
        int tslice_3;

        Sequential_Prop_1pt_Block_Collection()
        {
            is_build = false;
            tslice_1 = -1;
            tslice_2 = -1;
            tslice_3 = -1;
        }

        Sequential_Prop_1pt_Block_Collection(
            const std::string &job_tag, const int &traj,
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
            const FieldSelection &fsel,
            const Geometry &geo,
            const int &tslice_1_, const int &tslice_2_, const int &tslice_3_,
            const bool &is_baryon_1, const bool &is_baryon_3,
            const std::vector<std::vector<int>> &Momentum_Targets_curr)
        {
            TIMER_VERBOSE("Sequential_Prop_1pt_Block_Collection");
            is_build = true;
            tslice_1 = tslice_1_;
            tslice_2 = tslice_2_;
            tslice_3 = tslice_3_;
            contract_sequential_block_fsel_mate_ultra_sum(
                job_tag, traj,
                block_coll,
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
                fsel,
                geo,
                tslice_1, tslice_2, tslice_3,
                is_baryon_1, is_baryon_3,
                Momentum_Targets_curr);
        }
    };

    inline void contract_sequential_block_mate_coll(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const bool &is_baryon_1, const bool &is_baryon_3,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const int &dtmin, const int &dtmax,
        const std::vector<int> &tsep_src_list)
    {
        TIMER_VERBOSE("contract_sequential_block_mate_coll");
        const Coordinate &total_site = geo.total_site();
        int num_dt = dtmax - dtmin;

        block_coll.resize(total_site[3]);
        for (int t_1 = 0; t_1 < total_site[3]; ++t_1)
        {
            block_coll[t_1].resize(num_dt);
            for (int dt_1_3 = 0; dt_1_3 < num_dt; ++dt_1_3)
            {
                int t_3 = mod(t_1 - (dt_1_3 + dtmin), total_site[3]);
                int num_dt_2_3 = dt_1_3 + dtmin;
                block_coll[t_1][dt_1_3].resize(num_dt_2_3);
                for (int dt_2_3 = 0; dt_2_3 < num_dt_2_3; ++dt_2_3)
                {
                    displayln_info(fname + ssprintf(": t_1=%ld, dt_1_3=%ld, dt_2_3=%ld", t_1, dt_1_3, dt_2_3));
                    int t_2 = mod(t_3 + dt_2_3, total_site[3]);
                    block_coll[t_1][dt_1_3][dt_2_3] = Sequential_Prop_1pt_Block_Collection(
                        job_tag, traj,
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
                        fsel,
                        geo,
                        t_1, t_2, t_3,
                        is_baryon_1, is_baryon_3,
                        Momentum_Targets_curr);
                }
            }
            // Timer::autodisplay();
        }
    }

    inline void contract_sequential_block_mate_ultra_coll(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll_ex,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const bool &is_baryon_1, const bool &is_baryon_3,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const int &dtmin, const int &dtmax,
        const std::vector<int> &tsep_src_list)
    {
        TIMER_VERBOSE("contract_sequential_block_mate_ultra_coll");
        const Coordinate &total_site = geo.total_site();
        const int num_dtxy = get_baryon_pi_half_sep(total_site[3]);

        std::vector<int> tsep_src_list0 = tsep_src_list;
        std::vector<int>::iterator tsep_src_list_ite = std::unique(tsep_src_list0.begin(), tsep_src_list0.end());
        tsep_src_list0.erase(tsep_src_list_ite, tsep_src_list0.end());

        qassert(tsep_src_list0.size() == 1);
        block_coll_ex.resize(total_site[3]);
        for (int t_meson = 0; t_meson < total_site[3]; ++t_meson)
        {
            block_coll_ex[t_meson].resize(tsep_src_list0.size());
            for (long unsigned int ti0 = 0; ti0 < tsep_src_list0.size(); ++ti0)
            {
                int t_baryon_src = mod(t_meson - tsep_src_list0[ti0], total_site[3]);
                block_coll_ex[t_meson][ti0].resize(num_dtxy);
                for (int dt_xy = 0; dt_xy < num_dtxy; ++dt_xy)
                {
                    int xt = mod(t_meson + dt_xy, total_site[3]);
                    displayln_info(fname + ssprintf(": t_meson=%d, xt=%d, t_baryon_src=%d", t_meson, xt, t_baryon_src));
                    block_coll_ex[t_meson][ti0][dt_xy] = Sequential_Prop_1pt_Block_Collection(
                        job_tag, traj,
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
                        fsel,
                        geo,
                        t_meson, xt, t_baryon_src,
                        is_baryon_1, is_baryon_3,
                        Momentum_Targets_curr);
                }
            }
            // Timer::autodisplay();
        }
    }

    inline void glb_sum(Sequential_Prop_1pt_Block &block)
    {
        TIMER_VERBOSE("glb_sum Sequential_Prop_1pt_Block>");
        bool is_build_check;
        long xg_par1_psel_idx_check;
        long xg_par3_psel_idx_check;
        int n_snk_check;
        int n_src_check;

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank == 0)
        {
            is_build_check = block.is_build;
            xg_par1_psel_idx_check = block.xg_par1_psel_idx;
            xg_par3_psel_idx_check = block.xg_par3_psel_idx;
            n_snk_check = block.n_snk;
            n_src_check = block.n_src;
        }
        MPI_Bcast(&is_build_check, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xg_par1_psel_idx_check, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xg_par3_psel_idx_check, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n_snk_check, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n_src_check, 1, MPI_INT, 0, MPI_COMM_WORLD);

        qassert(is_build_check == block.is_build);
        qassert(xg_par1_psel_idx_check == block.xg_par1_psel_idx);
        qassert(xg_par3_psel_idx_check == block.xg_par3_psel_idx);
        qassert(n_snk_check == block.n_snk);
        qassert(n_src_check == block.n_src);

        // MPI_Barrier(MPI_COMM_WORLD);
        glb_sum(block.block);
    }

    inline void glb_sum(std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block_coll)
    {
        TIMER_VERBOSE("glb_sum std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>>");
        Vector<Sequential_Prop_1pt_Block> block_coll_1d;
        for (auto &block_col : block_coll)
        {
            for (auto &block_co : block_col)
            {
                for (auto &block_c : block_co)
                {
                    for (auto &block_ : block_c)
                    {
                        glb_sum(block_);
                    }
                }
            }
        }
        // for (int i = 0; i < block_coll.size(); i++)
        // {
        //     auto &block_col = block_coll[i];
        //     for (int j = 0; j < block_col.size(); j++)
        //     {
        //         auto &block_co = block_col[j];
        //         for (int k = 0; k < block_co.size(); k++)
        //         {
        //             auto &block_c = block_co[k];
        //             for (int l = 0; l < block_c.size(); l++)
        //             {
        //                 auto &block_ = block_c[l];
        //                 glb_sum(block_);
        //             }
        //         }
        //     }
        // }
    }

    inline void glb_sum(Sequential_Prop_1pt_Block_Collection &block_coll)
    {

        TIMER_VERBOSE("glb_sum Sequential_Prop_1pt_Block_Collection");
        bool is_build_check;
        int tslice_1_check;
        int tslice_2_check;
        int tslice_3_check;

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank == 0)
        {
            is_build_check = block_coll.is_build;
            tslice_1_check = block_coll.tslice_1;
            tslice_2_check = block_coll.tslice_2;
            tslice_3_check = block_coll.tslice_3;
        }
        MPI_Bcast(&is_build_check, sizeof(bool), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tslice_1_check, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tslice_2_check, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tslice_3_check, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

        qassert(is_build_check == block_coll.is_build);
        qassert(tslice_1_check == block_coll.tslice_1);
        qassert(tslice_2_check == block_coll.tslice_2);
        qassert(tslice_3_check == block_coll.tslice_3);

        glb_sum(block_coll.block_coll);

        // Timer::autodisplay();
    }

    inline void glb_sum(std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll)
    {
        TIMER_VERBOSE("glb_sum std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>>");
        for (auto &block_col : block_coll)
        {
            for (auto &block_co : block_col)
            {
                for (auto &block_c : block_co)
                {
                    glb_sum(block_c);
                }
            }
        }
    }

    inline void contract_sequential_block_x_fsel_y_psel(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<std::vector<Sequential_Prop_2pt_Block>>>>> &block_x_y,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src)
    {
        TIMER_VERBOSE("contract_sequential_block_x_fsel_y_psel");

        const Coordinate &total_site = geo.total_site();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            qassert(t_baryon_src == xg_src[3]);

            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        }

        block_x_y.resize(psel_baryon_num_list[t_baryon_snk]);
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            block_x_y[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                block_x_y[idx_baryon_snk][idx_baryon_src].resize(psel_meson_num_list[t_meson]);
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty].resize(fsel_num_list[xt]);
                    for (int idx_x = 0; idx_x < fsel_num_list[xt]; idx_x++)
                    {
                        block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x].resize(8);
                    }
                }
            }
        }

        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> block_x4snk_y;
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> block_x4y_src;

        contract_sequential_block_fsel(
            job_tag, traj, block_x4snk_y,
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
            fsel,
            geo,
            t_baryon_snk, xt, t_meson,
            true, false);
        contract_sequential_block_fsel(
            job_tag, traj, block_x4y_src,
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
            fsel,
            geo,
            t_meson, xt, t_baryon_src,
            false, true);

#pragma omp parallel for collapse(2)
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    for (int idx_x = 0; idx_x < fsel_num_list[xt]; idx_x++)
                    {
                        long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                        const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                        long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                        long xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                        const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                        const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                        // long xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                        // const Coordinate &xg_x = psel_meson[xg_x_psel_idx];

                        long xg_x_fsel_idx = idx_class_by_time_fsel[xt][idx_x];
                        const long index = fsel.indices[xg_x_fsel_idx];
                        const Coordinate xg_x_l = geo.coordinate_from_index(index);
                        const Coordinate xg_x_g = geo.coordinate_g_from_l(xg_x_l);
                        const Coordinate &xg_x = xg_x_g;
                        const int &t_x = xg_x[3];
                        qassert(t_x == xt);
                        qassert(idx_x == list_n_from_idx_fsel[xg_x_fsel_idx]);

                        qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                        qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                        qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);
                        qassert(list_n_from_idx_fsel[xg_x_fsel_idx] == idx_x);

                        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                        {
                            Sequential_Prop_1pt_Block &block_x4snk_y_temp = block_x4snk_y[idx_baryon_snk][idx_meson_ty][idx_x][VA_idx];
                            Sequential_Prop_1pt_Block &block_x4_y_src_temp = block_x4y_src[idx_meson_ty][idx_baryon_src][idx_x][VA_idx];

                            qassert(block_x4snk_y_temp.is_build);
                            qassert(block_x4snk_y_temp.xg_par1_psel_idx == xg_snk_psel_idx);
                            qassert(block_x4snk_y_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                            qassert(block_x4snk_y_temp.xg_par3_psel_idx == xg_y_psel_idx);

                            qassert(block_x4_y_src_temp.is_build);
                            qassert(block_x4_y_src_temp.xg_par1_psel_idx == xg_y_psel_idx);
                            qassert(block_x4_y_src_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                            qassert(block_x4_y_src_temp.xg_par3_psel_idx == xg_src_psel_idx);

                            const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                            const WilsonMatrix &wm_snk_y = prop_y.get_elem(xg_snk_psel_idx);

                            qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                            WilsonMatrix wm_sequence_0 = block_x4snk_y_temp.block * gamma5 * wm_y_src;
                            WilsonMatrix wm_sequence_1 = wm_snk_y * gamma5 * block_x4_y_src_temp.block;

                            Sequential_Prop_2pt_Block &block_x_y_temp = block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x][VA_idx];

                            qassert(!block_x_y_temp.is_build);
                            block_x_y_temp.is_build = true;
                            block_x_y_temp.xg_snk_psel_idx = xg_snk_psel_idx;
                            block_x_y_temp.xg_src_psel_idx = xg_src_psel_idx;
                            block_x_y_temp.xg_x_psel_idx = xg_x_fsel_idx;
                            block_x_y_temp.xg_y_psel_idx = xg_y_psel_idx;

                            block_x_y_temp.n_snk = idx_baryon_snk;
                            block_x_y_temp.n_src = idx_baryon_src;
                            block_x_y_temp.block[0] = wm_sequence_0;
                            block_x_y_temp.block[1] = wm_sequence_1;
                        }
                    }
                }
            }
        }
    }

    inline void contract_sequential_block_x_fsel_y_psel_ultra(
        const std::string &job_tag, const int &traj,
        std::vector<std::vector<std::vector<std::vector<std::vector<Sequential_Prop_2pt_Block>>>>> &block_x_y,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src)
    {
        TIMER_VERBOSE("contract_sequential_block_x_fsel_y_psel");

        const Coordinate &total_site = geo.total_site();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();

        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            qassert(t_baryon_src == xg_src[3]);

            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        }

        block_x_y.resize(psel_baryon_num_list[t_baryon_snk]);
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            block_x_y[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                block_x_y[idx_baryon_snk][idx_baryon_src].resize(psel_meson_num_list[t_meson]);
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty].resize((0 == fsel_num_list[xt] ? 0 : 1));
                    for (int idx_x = 0; idx_x < (0 == fsel_num_list[xt] ? 0 : 1); idx_x++)
                    {
                        block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x].resize(8);
                    }
                }
            }
        }

        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> block_x4snk_y;
        std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> block_x4y_src;

        contract_sequential_block_fsel_ultra_sum(
            job_tag, traj, block_x4snk_y,
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
            fsel,
            geo,
            t_baryon_snk, xt, t_meson,
            true, false);
        contract_sequential_block_fsel_ultra_sum(
            job_tag, traj, block_x4y_src,
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
            fsel,
            geo,
            t_meson, xt, t_baryon_src,
            false, true);

#pragma omp parallel for collapse(2)
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    for (int idx_x = 0; idx_x < (0 == fsel_num_list[xt] ? 0 : 1); idx_x++)
                    {
                        long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                        const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                        long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                        long xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                        const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                        const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                        // long xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                        // const Coordinate &xg_x = psel_meson[xg_x_psel_idx];

                        // long xg_x_fsel_idx = idx_class_by_time_fsel[xt][idx_x];
                        // const long index = fsel.indices[xg_x_fsel_idx];
                        // const Coordinate xg_x_l = geo.coordinate_from_index(index);
                        // const Coordinate xg_x_g = geo.coordinate_g_from_l(xg_x_l);
                        // const Coordinate &xg_x = xg_x_g;
                        // const int &t_x = xg_x[3];
                        // qassert(t_x == xt);
                        // qassert(idx_x == list_n_from_idx_fsel[xg_x_fsel_idx]);

                        qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                        qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                        qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);
                        // qassert(list_n_from_idx_fsel[xg_x_fsel_idx] == idx_x);

                        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                        {
                            Sequential_Prop_1pt_Block &block_x4snk_y_temp = block_x4snk_y[idx_baryon_snk][idx_meson_ty][idx_x][VA_idx];
                            Sequential_Prop_1pt_Block &block_x4_y_src_temp = block_x4y_src[idx_meson_ty][idx_baryon_src][idx_x][VA_idx];

                            qassert(block_x4snk_y_temp.is_build);
                            qassert(block_x4snk_y_temp.xg_par1_psel_idx == xg_snk_psel_idx);
                            // qassert(block_x4snk_y_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                            qassert(block_x4snk_y_temp.xg_par3_psel_idx == xg_y_psel_idx);

                            qassert(block_x4_y_src_temp.is_build);
                            qassert(block_x4_y_src_temp.xg_par1_psel_idx == xg_y_psel_idx);
                            // qassert(block_x4_y_src_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                            qassert(block_x4_y_src_temp.xg_par3_psel_idx == xg_src_psel_idx);

                            const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                            const WilsonMatrix &wm_snk_y = prop_y.get_elem(xg_snk_psel_idx);

                            qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                            WilsonMatrix wm_sequence_0 = block_x4snk_y_temp.block * gamma5 * wm_y_src;
                            WilsonMatrix wm_sequence_1 = wm_snk_y * gamma5 * block_x4_y_src_temp.block;

                            Sequential_Prop_2pt_Block &block_x_y_temp = block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x][VA_idx];

                            qassert(!block_x_y_temp.is_build);
                            block_x_y_temp.is_build = true;
                            block_x_y_temp.xg_snk_psel_idx = xg_snk_psel_idx;
                            block_x_y_temp.xg_src_psel_idx = xg_src_psel_idx;
                            // block_x_y_temp.xg_x_psel_idx = xg_x_fsel_idx;
                            block_x_y_temp.xg_y_psel_idx = xg_y_psel_idx;

                            block_x_y_temp.n_snk = idx_baryon_snk;
                            block_x_y_temp.n_src = idx_baryon_src;
                            block_x_y_temp.block[0] = wm_sequence_0;
                            block_x_y_temp.block[1] = wm_sequence_1;
                        }
                    }
                }
            }
        }
    }

    inline void contract_sequential_block_x_fsel_y_psel_mate_ultra(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll,
        const std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll_ex,
        std::vector<std::vector<std::vector<std::vector<std::vector<Sequential_Prop_2pt_Block>>>>> &block_x_y,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src)
    {
        TIMER_VERBOSE("contract_sequential_block_x_fsel_y_psel_mate_ultra");

        const Coordinate &total_site = geo.total_site();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const int &mom_num_curr = Momentum_Targets_curr.size();

        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            qassert(t_baryon_src == xg_src[3]);

            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        }

        block_x_y.resize(psel_baryon_num_list[t_baryon_snk]);
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            block_x_y[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                block_x_y[idx_baryon_snk][idx_baryon_src].resize(psel_meson_num_list[t_meson]);
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    if (0 != fsel_num_list[xt])
                    {
                        block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty].resize(mom_num_curr);
                        for (int idx_x = 0; idx_x < mom_num_curr; idx_x++)
                        {
                            block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x].resize(8);
                        }
                    }
                }
            }
        }

        const Sequential_Prop_1pt_Block_Collection &block_x4snk_y_ = block_coll[t_baryon_snk][mod(t_baryon_snk - t_meson, total_site[3]) - dtmin][mod(xt - t_meson, total_site[3])];

        // displayln(ssprintf("&block_x4snk_y_ = block_coll[%d][%d][%d]", t_baryon_snk, (mod(t_baryon_snk - t_meson, total_site[3]) - dtmin), (mod(xt - t_meson, total_site[3]))));

        qassert(block_x4snk_y_.is_build);
        qassert(block_x4snk_y_.tslice_1 == t_baryon_snk);
        qassert(block_x4snk_y_.tslice_2 == xt);
        qassert(block_x4snk_y_.tslice_3 == t_meson);

        const Sequential_Prop_1pt_Block_Collection &block_x4y_src_ = block_coll_ex[t_meson][0][mod(xt - t_meson, total_site[3])];

        // displayln_info(ssprintf("&t_meson=%d, xt=%d, t_baryon_src=%d, block_x4y_src_.tslice_2 = %d", t_meson, xt, t_baryon_src, block_x4y_src_.tslice_2));

        qassert(block_x4y_src_.is_build);
        qassert(block_x4y_src_.tslice_1 == t_meson);
        qassert(block_x4y_src_.tslice_2 == xt);
        qassert(block_x4y_src_.tslice_3 == t_baryon_src);

        const std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block_x4snk_y = block_x4snk_y_.block_coll;
        const std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block_x4y_src = block_x4y_src_.block_coll;

#pragma omp parallel for collapse(2)
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    if (0 != fsel_num_list[xt])
                    {
                        for (int idx_x = 0; idx_x < mom_num_curr; idx_x++)
                        {
                            long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                            const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                            long xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                            const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                            const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                            // long xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                            // const Coordinate &xg_x = psel_meson[xg_x_psel_idx];

                            // long xg_x_fsel_idx = idx_class_by_time_fsel[xt][idx_x];
                            // const long index = fsel.indices[xg_x_fsel_idx];
                            // const Coordinate xg_x_l = geo.coordinate_from_index(index);
                            // const Coordinate xg_x_g = geo.coordinate_g_from_l(xg_x_l);
                            // const Coordinate &xg_x = xg_x_g;
                            // const int &t_x = xg_x[3];
                            // qassert(t_x == xt);
                            // qassert(idx_x == list_n_from_idx_fsel[xg_x_fsel_idx]);

                            qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                            qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                            qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);
                            // qassert(list_n_from_idx_fsel[xg_x_fsel_idx] == idx_x);

                            for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                            {
                                const Sequential_Prop_1pt_Block &block_x4snk_y_temp = block_x4snk_y[idx_baryon_snk][idx_meson_ty][idx_x][VA_idx];
                                const Sequential_Prop_1pt_Block &block_x4_y_src_temp = block_x4y_src[idx_meson_ty][idx_baryon_src][idx_x][VA_idx];

                                qassert(block_x4snk_y_temp.is_build);
                                qassert(block_x4snk_y_temp.xg_par1_psel_idx == xg_snk_psel_idx);
                                // qassert(block_x4snk_y_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                                qassert(block_x4snk_y_temp.xg_par3_psel_idx == xg_y_psel_idx);

                                qassert(block_x4_y_src_temp.is_build);
                                qassert(block_x4_y_src_temp.xg_par1_psel_idx == xg_y_psel_idx);
                                // qassert(block_x4_y_src_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                                qassert(block_x4_y_src_temp.xg_par3_psel_idx == xg_src_psel_idx);

                                const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                                const WilsonMatrix &wm_snk_y = prop_y.get_elem(xg_snk_psel_idx);

                                qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                                WilsonMatrix wm_sequence_0 = block_x4snk_y_temp.block * gamma5 * wm_y_src;
                                WilsonMatrix wm_sequence_1 = wm_snk_y * gamma5 * block_x4_y_src_temp.block;

                                Sequential_Prop_2pt_Block &block_x_y_temp = block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x][VA_idx];

                                qassert(!block_x_y_temp.is_build);
                                block_x_y_temp.is_build = true;
                                block_x_y_temp.xg_snk_psel_idx = xg_snk_psel_idx;
                                block_x_y_temp.xg_src_psel_idx = xg_src_psel_idx;
                                // block_x_y_temp.xg_x_psel_idx = xg_x_fsel_idx;
                                block_x_y_temp.xg_y_psel_idx = xg_y_psel_idx;

                                block_x_y_temp.n_snk = idx_baryon_snk;
                                block_x_y_temp.n_src = idx_baryon_src;
                                block_x_y_temp.block[0] = wm_sequence_0;
                                block_x_y_temp.block[1] = wm_sequence_1;
                            }
                        }
                    }
                }
            }
        }
    }

    inline void contract_sequential_block_x_fsel_y_psel_mate_ultra_magic(
        const std::string &job_tag, const int &traj,
        const std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll,
        const std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block_Collection>>> &block_coll_ex,
        std::vector<std::vector<std::vector<std::vector<std::vector<Sequential_Prop_2pt_Block>>>>> &block_x_y,
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
        const FieldSelection &fsel,
        const Geometry &geo,
        const std::vector<std::vector<int>> &Momentum_Targets_curr,
        const int &dtmax, const int &dtmin,
        const int &t_baryon_snk, const int &xt, const int &t_meson, const int &t_baryon_src)
    {
        TIMER_VERBOSE("contract_sequential_block_x_fsel_y_psel_mate_ultra_magic");

        const Coordinate &total_site = geo.total_site();
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const int &mom_num_curr = Momentum_Targets_curr.size();

        for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
        {
            long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
            const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
            qassert(t_baryon_src == xg_src[3]);

            const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);
        }

        block_x_y.resize(psel_baryon_num_list[t_baryon_snk]);
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            block_x_y[idx_baryon_snk].resize(psel_baryon_num_list[t_baryon_src]);
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                block_x_y[idx_baryon_snk][idx_baryon_src].resize(psel_meson_num_list[t_meson]);
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty].resize(mom_num_curr);
                    for (int idx_x = 0; idx_x < mom_num_curr; idx_x++)
                    {
                        block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x].resize(8);
                    }
                }
            }
        }

        const Sequential_Prop_1pt_Block_Collection &block_x4snk_y_ = block_coll[t_baryon_snk][mod(t_baryon_snk - t_meson, total_site[3]) - dtmin][mod(xt - t_meson, total_site[3])];

        // displayln(ssprintf("&block_x4snk_y_ = block_coll[%d][%d][%d]", t_baryon_snk, (mod(t_baryon_snk - t_meson, total_site[3]) - dtmin), (mod(xt - t_meson, total_site[3]))));

        qassert(block_x4snk_y_.is_build);
        qassert(block_x4snk_y_.tslice_1 == t_baryon_snk);
        qassert(block_x4snk_y_.tslice_2 == xt);
        qassert(block_x4snk_y_.tslice_3 == t_meson);

        const Sequential_Prop_1pt_Block_Collection &block_x4y_src_ = block_coll_ex[t_meson][0][mod(xt - t_meson, total_site[3])];

        // displayln_info(ssprintf("&t_meson=%d, xt=%d, t_baryon_src=%d, block_x4y_src_.tslice_2 = %d", t_meson, xt, t_baryon_src, block_x4y_src_.tslice_2));

        qassert(block_x4y_src_.is_build);
        qassert(block_x4y_src_.tslice_1 == t_meson);
        qassert(block_x4y_src_.tslice_2 == xt);
        qassert(block_x4y_src_.tslice_3 == t_baryon_src);

        const std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block_x4snk_y = block_x4snk_y_.block_coll;
        const std::vector<std::vector<std::vector<std::vector<Sequential_Prop_1pt_Block>>>> &block_x4y_src = block_x4y_src_.block_coll;

#pragma omp parallel for collapse(2)
        for (int idx_baryon_snk = 0; idx_baryon_snk < psel_baryon_num_list[t_baryon_snk]; idx_baryon_snk++)
        {
            for (int idx_baryon_src = 0; idx_baryon_src < psel_baryon_num_list[t_baryon_src]; idx_baryon_src++)
            {
                for (int idx_meson_ty = 0; idx_meson_ty < psel_meson_num_list[t_meson]; idx_meson_ty++)
                {
                    for (int idx_x = 0; idx_x < mom_num_curr; idx_x++)
                    {
                        long xg_snk_psel_idx = idx_class_by_time_baryon[t_baryon_snk][idx_baryon_snk];
                        const Coordinate &xg_snk = psel_baryon[xg_snk_psel_idx];

                        long xg_src_psel_idx = idx_class_by_time_baryon[t_baryon_src][idx_baryon_src];
                        const Coordinate &xg_src = psel_baryon[xg_src_psel_idx];
                        const PselProp &prop_src = get_psel_prop(job_tag, traj, xg_src, is_meson_smear, is_baryon_smear);

                        long xg_y_psel_idx = idx_class_by_time_meson[t_meson][idx_meson_ty];
                        const Coordinate &xg_y = psel_meson[xg_y_psel_idx];
                        const PselProp &prop_y = get_psel_prop(job_tag, traj, xg_y, is_baryon_smear, is_meson_smear);

                        // long xg_x_psel_idx = idx_class_by_time_meson[xt][idx_x];
                        // const Coordinate &xg_x = psel_meson[xg_x_psel_idx];

                        // long xg_x_fsel_idx = idx_class_by_time_fsel[xt][idx_x];
                        // const long index = fsel.indices[xg_x_fsel_idx];
                        // const Coordinate xg_x_l = geo.coordinate_from_index(index);
                        // const Coordinate xg_x_g = geo.coordinate_g_from_l(xg_x_l);
                        // const Coordinate &xg_x = xg_x_g;
                        // const int &t_x = xg_x[3];
                        // qassert(t_x == xt);
                        // qassert(idx_x == list_n_from_idx_fsel[xg_x_fsel_idx]);

                        // qassert(list_n_from_idx_baryon[xg_snk_psel_idx] == idx_baryon_snk);
                        // qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                        // qassert(list_n_from_idx_meson[xg_y_psel_idx] == idx_meson_ty);
                        // qassert(list_n_from_idx_fsel[xg_x_fsel_idx] == idx_x);

                        for (int VA_idx = 0; VA_idx < 8; VA_idx++)
                        {
                            const Sequential_Prop_1pt_Block &block_x4snk_y_temp = block_x4snk_y[idx_baryon_snk][idx_meson_ty][idx_x][VA_idx];
                            const Sequential_Prop_1pt_Block &block_x4_y_src_temp = block_x4y_src[idx_meson_ty][idx_baryon_src][idx_x][VA_idx];

                            qassert(block_x4snk_y_temp.is_build);
                            qassert(block_x4snk_y_temp.xg_par1_psel_idx == xg_snk_psel_idx);
                            // qassert(block_x4snk_y_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                            qassert(block_x4snk_y_temp.xg_par3_psel_idx == xg_y_psel_idx);

                            qassert(block_x4_y_src_temp.is_build);
                            qassert(block_x4_y_src_temp.xg_par1_psel_idx == xg_y_psel_idx);
                            // qassert(block_x4_y_src_temp.xg_par2_psel_idx == xg_x_fsel_idx);
                            qassert(block_x4_y_src_temp.xg_par3_psel_idx == xg_src_psel_idx);

                            const WilsonMatrix &wm_y_src = prop_src.get_elem(xg_y_psel_idx);
                            const WilsonMatrix &wm_snk_y = prop_y.get_elem(xg_snk_psel_idx);

                            qassert(list_n_from_idx_baryon[xg_src_psel_idx] == idx_baryon_src);
                            WilsonMatrix wm_sequence_0 = block_x4snk_y_temp.block * gamma5 * wm_y_src;
                            WilsonMatrix wm_sequence_1 = wm_snk_y * gamma5 * block_x4_y_src_temp.block;

                            Sequential_Prop_2pt_Block &block_x_y_temp = block_x_y[idx_baryon_snk][idx_baryon_src][idx_meson_ty][idx_x][VA_idx];

                            qassert(!block_x_y_temp.is_build);
                            block_x_y_temp.is_build = true;
                            block_x_y_temp.xg_snk_psel_idx = xg_snk_psel_idx;
                            block_x_y_temp.xg_src_psel_idx = xg_src_psel_idx;
                            // block_x_y_temp.xg_x_psel_idx = xg_x_fsel_idx;
                            block_x_y_temp.xg_y_psel_idx = xg_y_psel_idx;

                            block_x_y_temp.n_snk = idx_baryon_snk;
                            block_x_y_temp.n_src = idx_baryon_src;
                            block_x_y_temp.block[0] = wm_sequence_0;
                            block_x_y_temp.block[1] = wm_sequence_1;
                        }
                    }
                }
            }
        }
    }

    inline void ppblock_np_2_LT_block_np(Leading_Twist_Block_No_Projection_Psel_Psel &LT_block_np, const Baryon_Block_No_Projection_Psel_Psel &ppblock_np, const WilsonMatrix &g5_wm_y_src, const WilsonMatrix &wm_snk_y_g5)
    {
        TIMER_VERBOSE("ppblock_np_2_LT_block_np");

        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            LT_block_np.block_1_2[0][mu3] = g5_wm_y_src * ppblock_np.block_1_2[mu3];
            LT_block_np.block_1_3[1][mu3] = ppblock_np.block_1_3[mu3] * wm_snk_y_g5;

            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                LT_block_np.block_0_12[0][mu4 * 4 + mu3] = g5_wm_y_src * ppblock_np.block_0_12[mu4 * 4 + mu3];
                LT_block_np.block_0_12[1][mu4 * 4 + mu3] = ppblock_np.block_0_12[mu4 * 4 + mu3] * wm_snk_y_g5;
                LT_block_np.block_1_1[0][mu4 * 4 + mu3] = g5_wm_y_src * ppblock_np.block_1_1[mu4 * 4 + mu3];
                LT_block_np.block_1_1[1][mu4 * 4 + mu3] = ppblock_np.block_1_1[mu4 * 4 + mu3] * wm_snk_y_g5;

                for (int mu = 0; mu < 4; mu++)
                {
                    for (int c = 0; c < NUM_COLOR; c++)
                    {
                        for (int c1 = 0; c1 < NUM_COLOR; c1++)
                        {
                            for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                            {
                                LT_block_np.block_0_3[0][mu4](mu * NUM_COLOR + c, mu3 * NUM_COLOR + c1) += g5_wm_y_src(mu * NUM_COLOR + c, mu4 * NUM_COLOR + c1p) * ppblock_np.block_0_3(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1);

                                LT_block_np.block_0_3[1][mu3](mu4 * NUM_COLOR + c1p, mu * NUM_COLOR + c) += ppblock_np.block_0_3(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * wm_snk_y_g5(mu3 * NUM_COLOR + c1, mu * NUM_COLOR + c);

                                for (int mup = 0; mup < 4; mup++)
                                {
                                    LT_block_np.block_1_2[1][mu4 * 4 + mu3](mup * NUM_COLOR + c1p, mu * NUM_COLOR + c) += ppblock_np.block_1_2[mu4](mup * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * wm_snk_y_g5(mu3 * NUM_COLOR + c1, mu * NUM_COLOR + c);
                                    LT_block_np.block_1_3[0][mu4 * 4 + mu3](mu * NUM_COLOR + c, mup * NUM_COLOR + c1) += g5_wm_y_src(mu * NUM_COLOR + c, mu4 * NUM_COLOR + c1p) * ppblock_np.block_1_3[mu3](mu4 * NUM_COLOR + c1p, mup * NUM_COLOR + c1);
                                }
                            }
                        }
                    }
                }
            }
        }

        LT_block_np.is_build = true;

        return;
    }

    inline void ppblock_np_2_LT_block_np(
        Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra &LT_block_unex,
        Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra &LT_block_ex,
        const Baryon_Block_No_Projection_Psel_Psel &ppblock_np,
        const WilsonMatrix &g5_wm_y_src, const WilsonMatrix &wm_snk_y_g5)
    {
        TIMER_VERBOSE("ppblock_np_2_LT_block_np_pro");

        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            LT_block_unex.block_1_2[mu3] += g5_wm_y_src * ppblock_np.block_1_2[mu3];
            LT_block_ex.block_1_3[mu3] += ppblock_np.block_1_3[mu3] * wm_snk_y_g5;

            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                LT_block_unex.block_0_12[mu4 * 4 + mu3] += g5_wm_y_src * ppblock_np.block_0_12[mu4 * 4 + mu3];
                LT_block_ex.block_0_12[mu4 * 4 + mu3] += ppblock_np.block_0_12[mu4 * 4 + mu3] * wm_snk_y_g5;
                LT_block_unex.block_1_1[mu4 * 4 + mu3] += g5_wm_y_src * ppblock_np.block_1_1[mu4 * 4 + mu3];
                LT_block_ex.block_1_1[mu4 * 4 + mu3] += ppblock_np.block_1_1[mu4 * 4 + mu3] * wm_snk_y_g5;

                for (int mu = 0; mu < 4; mu++)
                {
                    for (int c = 0; c < NUM_COLOR; c++)
                    {
                        for (int c1 = 0; c1 < NUM_COLOR; c1++)
                        {
                            for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                            {
                                LT_block_unex.block_0_3[mu4](mu * NUM_COLOR + c, mu3 * NUM_COLOR + c1) += g5_wm_y_src(mu * NUM_COLOR + c, mu4 * NUM_COLOR + c1p) * ppblock_np.block_0_3(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1);

                                LT_block_ex.block_0_3[mu3](mu4 * NUM_COLOR + c1p, mu * NUM_COLOR + c) += ppblock_np.block_0_3(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * wm_snk_y_g5(mu3 * NUM_COLOR + c1, mu * NUM_COLOR + c);

                                for (int mup = 0; mup < 4; mup++)
                                {
                                    LT_block_ex.block_1_2[mu4 * 4 + mu3](mup * NUM_COLOR + c1p, mu * NUM_COLOR + c) += ppblock_np.block_1_2[mu4](mup * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * wm_snk_y_g5(mu3 * NUM_COLOR + c1, mu * NUM_COLOR + c);
                                    LT_block_unex.block_1_3[mu4 * 4 + mu3](mu * NUM_COLOR + c, mup * NUM_COLOR + c1) += g5_wm_y_src(mu * NUM_COLOR + c, mu4 * NUM_COLOR + c1p) * ppblock_np.block_1_3[mu3](mu4 * NUM_COLOR + c1p, mup * NUM_COLOR + c1);
                                }
                            }
                        }
                    }
                }
            }
        }

        LT_block_unex.is_build = true;
        LT_block_ex.is_build = true;

        return;
    }

    inline void contract_swm_ppblock_np_nova(std::vector<SpinMatrix> &spin_result_list, std::vector<SpinMatrix> &spin_result_list_ex, const Leading_Twist_Block_No_Projection_Psel_Psel &LT_block_np, const Sequential_Prop_1pt_Block &block_x4snk_y, const Sequential_Prop_1pt_Block &block_x4_y_src)
    {
        TIMER_VERBOSE("contract_swm_ppblock_np_nova");

        spin_result_list.resize(5);
        spin_result_list_ex.resize(5);
        for (int topology_id = 0; topology_id < 5; topology_id++)
        {
            set_zero(spin_result_list[topology_id]);
            set_zero(spin_result_list_ex[topology_id]);
        }

        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                spin_result_list[1](mu4, mu3) = matrix_trace(LT_block_np.block_0_12[0][mu4 * 4 + mu3] * block_x4snk_y.block);
                spin_result_list_ex[1](mu4, mu3) = matrix_trace(LT_block_np.block_0_12[1][mu4 * 4 + mu3] * block_x4_y_src.block);

                spin_result_list[4](mu4, mu3) = matrix_trace(LT_block_np.block_1_1[0][mu4 * 4 + mu3] * block_x4snk_y.block);
                spin_result_list_ex[4](mu4, mu3) = matrix_trace(LT_block_np.block_1_1[1][mu4 * 4 + mu3] * block_x4_y_src.block);

                spin_result_list_ex[2](mu4, mu3) = matrix_trace(LT_block_np.block_1_2[1][mu4 * 4 + mu3] * block_x4_y_src.block);
                spin_result_list[3](mu4, mu3) = matrix_trace(LT_block_np.block_1_3[0][mu4 * 4 + mu3] * block_x4snk_y.block);

                for (int c1 = 0; c1 < NUM_COLOR; c1++)
                {
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            spin_result_list[0](mu4, mu3) += LT_block_np.block_0_3[0][mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * block_x4snk_y.block(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);
                            spin_result_list_ex[0](mu4, mu3) += LT_block_np.block_0_3[1][mu3](mu4 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * block_x4_y_src.block(mu1 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);

                            spin_result_list[2](mu4, mu3) += LT_block_np.block_1_2[0][mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * block_x4snk_y.block(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);
                            spin_result_list_ex[3](mu4, mu3) += LT_block_np.block_1_3[1][mu3](mu4 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * block_x4_y_src.block(mu1 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                        }
                    }
                }
            }
        }
        return;
    }

    

    inline void contract_swm_ppblock_np_nova_unex(
        std::vector<SpinMatrix> &spin_result_list, 
    const Leading_Twist_Block_No_Projection_Psel_Psel_Unex_Ultra &LT_block_unex, 
    const Sequential_Prop_1pt_Block &block_x4snk_y)
    {
        TIMER_VERBOSE("contract_swm_ppblock_np_nova");

        spin_result_list.resize(5);
        for (int topology_id = 0; topology_id < 5; topology_id++)
        {
            set_zero(spin_result_list[topology_id]);
        }

        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                spin_result_list[1](mu4, mu3) = matrix_trace(LT_block_unex.block_0_12[mu4 * 4 + mu3] * block_x4snk_y.block);

                spin_result_list[4](mu4, mu3) = matrix_trace(LT_block_unex.block_1_1[mu4 * 4 + mu3] * block_x4snk_y.block);

                spin_result_list[3](mu4, mu3) = matrix_trace(LT_block_unex.block_1_3[mu4 * 4 + mu3] * block_x4snk_y.block);

                for (int c1 = 0; c1 < NUM_COLOR; c1++)
                {
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            spin_result_list[0](mu4, mu3) += LT_block_unex.block_0_3[mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * block_x4snk_y.block(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);

                            spin_result_list[2](mu4, mu3) += LT_block_unex.block_1_2[mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * block_x4snk_y.block(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);
                        }
                    }
                }
            }
        }
        return;
    }

    

    inline void contract_swm_ppblock_np_nova_ex(std::vector<SpinMatrix> &spin_result_list_ex, 
    const Leading_Twist_Block_No_Projection_Psel_Psel_Ex_Ultra &LT_block_ex, 
    const Sequential_Prop_1pt_Block &block_x4_y_src)
    {
        TIMER_VERBOSE("contract_swm_ppblock_np_nova_ex");

        spin_result_list_ex.resize(5);
        for (int topology_id = 0; topology_id < 5; topology_id++)
        {
            set_zero(spin_result_list_ex[topology_id]);
        }

        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                spin_result_list_ex[1](mu4, mu3) = matrix_trace(LT_block_ex.block_0_12[mu4 * 4 + mu3] * block_x4_y_src.block);

                spin_result_list_ex[4](mu4, mu3) = matrix_trace(LT_block_ex.block_1_1[mu4 * 4 + mu3] * block_x4_y_src.block);

                spin_result_list_ex[2](mu4, mu3) = matrix_trace(LT_block_ex.block_1_2[mu4 * 4 + mu3] * block_x4_y_src.block);

                for (int c1 = 0; c1 < NUM_COLOR; c1++)
                {
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            spin_result_list_ex[0](mu4, mu3) += LT_block_ex.block_0_3[mu3](mu4 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * block_x4_y_src.block(mu1 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);

                            spin_result_list_ex[3](mu4, mu3) += LT_block_ex.block_1_3[mu3](mu4 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * block_x4_y_src.block(mu1 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                        }
                    }
                }
            }
        }
        return;
    }

    inline void contract_swm_ppblock_np(std::vector<SpinMatrix> &spin_result_list, const Baryon_Block_No_Projection_Psel_Psel &ppblock_np, const WilsonMatrix &swm)
    {
        TIMER_VERBOSE("contract_swm_ppblock_np");
        spin_result_list.resize(5);
        for (int topology_id = 0; topology_id < 5; topology_id++)
        {
            set_zero(spin_result_list[topology_id]);
        }

        // spin_result_list[0]
        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                spin_result_list[1](mu4, mu3) = matrix_trace(ppblock_np.block_0_12[mu4 * 4 + mu3] * swm);
                spin_result_list[4](mu4, mu3) = matrix_trace(ppblock_np.block_1_1[mu4 * 4 + mu3] * swm);

                for (int c1 = 0; c1 < NUM_COLOR; c1++)
                {
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        spin_result_list[0](mu4, mu3) += ppblock_np.block_0_3(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * swm(mu3 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            spin_result_list[2](mu4, mu3) += ppblock_np.block_1_2[mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * swm(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);
                            spin_result_list[3](mu4, mu3) += ppblock_np.block_1_3[mu3](mu4 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * swm(mu1 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                        }
                    }
                }
            }
        }
    }

    inline void contract_swm_spblock_np(std::vector<SpinMatrix> &spin_result_list, const Baryon_SPBlock_No_Projection_Psel_Psel &spblock_np, const WilsonMatrix &swm)
    {
        TIMER_VERBOSE("contract_swm_spblock_np");
        spin_result_list.resize(10);
        for (int topology_id = 0; topology_id < 10; topology_id++)
        {
            set_zero(spin_result_list[topology_id]);
        }

        // spin_result_list[0]
        for (int mu3 = 0; mu3 < 4; mu3++)
        {
            for (int mu4 = 0; mu4 < 4; mu4++)
            {
                spin_result_list[1](mu4, mu3) = matrix_trace(spblock_np.spblock_6[mu4 * 4 + mu3] * swm);

                spin_result_list[5](mu4, mu3) = matrix_trace(spblock_np.spblock_ex_5[mu4 * 4 + mu3] * swm);
                spin_result_list[6](mu4, mu3) = matrix_trace(spblock_np.spblock_ex_6[mu4 * 4 + mu3] * swm);
                spin_result_list[7](mu4, mu3) = matrix_trace(spblock_np.spblock_ex_7[mu4 * 4 + mu3] * swm);
                spin_result_list[8](mu4, mu3) = matrix_trace(spblock_np.spblock_ex_8[mu4 * 4 + mu3] * swm);

                for (int c1 = 0; c1 < NUM_COLOR; c1++)
                {
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        spin_result_list[0](mu4, mu3) += spblock_np.spblock_5(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * swm(mu3 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                        for (int mu1 = 0; mu1 < 4; mu1++)
                        {
                            spin_result_list[2](mu4, mu3) += spblock_np.spblock_7[mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * swm(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);
                            spin_result_list[3](mu4, mu3) += spblock_np.spblock_8[mu4](mu3 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * swm(mu1 * NUM_COLOR + c1, mu3 * NUM_COLOR + c1p);
                            spin_result_list[4](mu4, mu3) += spblock_np.spblock_9[mu4](mu1 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) * swm(mu3 * NUM_COLOR + c1, mu1 * NUM_COLOR + c1p);
                            spin_result_list[9](mu4, mu3) += spblock_np.spblock_ex_9[mu4](mu3 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) * swm(mu1 * NUM_COLOR + c1, mu3 * NUM_COLOR + c1p);
                        }
                    }
                }
            }
        }
    }
}