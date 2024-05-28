#pragma once

// #include "data-load.h"
// #include "baryon-pion-block-acc.h"
// #include "baryon-pion-coefficient.h"
#include "compute-utils.h"
// #include "baryon_pi_debug/baryon-pion-confrim-coeff.h"
// #include "baryon_pi_debug/baryon-pion-confrim-block.h"

namespace qlat
{
    inline void three_point_time_slice(const int &dtxy, const int &t_src_snk, const int &tsep_src, const int &tsep_snk, const Coordinate &total_site, int &t_meson, const int &t_baryon_src, int &t_baryon_snk, int &xt, double &anti_period)
    {
        qassert(t_src_snk >= dtxy);
        t_baryon_snk = mod(t_baryon_src + tsep_src + t_src_snk + tsep_snk, total_site[3]);
        t_meson = mod(t_baryon_src + tsep_src, total_site[3]);
        xt = mod(t_meson + dtxy, total_site[3]);

        if ((t_baryon_src <= t_baryon_snk) && (t_meson <= xt))
        {
            anti_period = 1.0;
        }
        else
        {
            anti_period = -1.0;
        }
    }

    inline void four_point_time_slice(const int &dtxy, const int &num_dtxy, const int &tsep_src, const int &tsep_snk, const Coordinate &total_site, const int &t_meson, int &t_baryon_src, int &t_baryon_snk, int &xt, double &anti_period)
    {
        // int xt = (tsep >= 0) ? mod(t_meson + dtxy, total_site[3]) : (dtxy >= 0 ? mod(t_meson - 2 * tsep + dtxy, total_site[3]) : mod(t_meson + 2 * tsep + dtxy, total_site[3]));
        if (tsep_src >= 0 && tsep_snk >= 0)
        {
            xt = mod(t_meson + dtxy, total_site[3]);
        }
        else if (tsep_src < 0 && tsep_snk < 0)
        {
            if (dtxy >= 0)
            {
                xt = mod(t_meson - tsep_snk - tsep_src + dtxy, total_site[3]);
            }
            else
            {
                xt = mod(t_meson + tsep_snk + tsep_src + dtxy, total_site[3]);
            }
        }
        else if (tsep_src >= 0 && tsep_snk < 0)
        {
            if (dtxy >= 0)
            {
                xt = mod(t_meson - tsep_snk + dtxy, total_site[3]);
            }
            else
            {
                xt = mod(t_meson + tsep_snk + dtxy, total_site[3]);
            }
        }
        else if (tsep_src < 0 && tsep_snk >= 0)
        {
            if (dtxy >= 0)
            {
                xt = mod(t_meson - tsep_snk + dtxy, total_site[3]);
            }
            else
            {
                xt = mod(t_meson + tsep_snk + dtxy, total_site[3]);
            }
        }
        else
        {
            qassert(false);
        }

        if (dtxy >= 0)
        {
            t_baryon_src = mod(t_meson - tsep_src, total_site[3]);
            t_baryon_snk = mod(xt + tsep_snk, total_site[3]);
        }
        else
        {
            t_baryon_src = mod(xt - tsep_src, total_site[3]);
            t_baryon_snk = mod(t_meson + tsep_snk, total_site[3]);
        }

        if ((t_baryon_src <= t_baryon_snk) && (t_meson <= xt))
        {
            anti_period = 1.0;
        }
        else
        {
            anti_period = -1.0;
        }
    }
}

// const std::vector<int> &tsep_src_list, const std::vector<int> &tsep_snk_list,
//     const int &tsep_src, const int &tsep_snk,

//     const int &tsep_src = tsep_src_list[ti];
// const int &tsep_snk = tsep_snk_list[ti];

// const int total_time_length = abs(dtxy) + abs(tsep_src) + abs(tsep_snk);
// const int dtxy = dt - num_dtxy / 2;
// int t_baryon_src, t_baryon_snk;
// int xt, t_baryon_src, t_baryon_snk;
// four_point_time_slice(dtxy, dt, num_dtxy, tsep_src, tsep_snk, total_site, t_meson, t_baryon_src, t_baryon_snk, xt);

//      if ((t_baryon_src <= t_baryon_snk) && (t_meson <= xt))