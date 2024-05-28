#pragma once

#include <qlat/contract-pion.h>
#include "block_contraction.h"

namespace qlat
{ //
    template <class T>
    qacc ComplexD pion_twop_block(const WilsonMatrixT<T> &m1, const WilsonMatrixT<T> &m2)
    {
        ComplexD ret = pion_2pt_block(m1, m2);
        return ret;
    };

    template <class T>
    qacc ComplexD proton_twop_block(const WilsonMatrixT<T> &m1, const WilsonMatrixT<T> &m2, const WilsonMatrixT<T> &m3, const int type)
    {
        const SpinMatrix &Cg5 = BaryonMatrixConstants::get_Cgamma5();
        const SpinMatrix &proj = BaryonMatrixConstants::get_projector();
        ComplexD ret = 0;
        if (type == 0)
        {
            const WilsonMatrix wm1 = Cg5 * m1 * Cg5;
            const WilsonMatrix wm2 = m2;
            const WilsonMatrix wm3 = m3 * proj;
            ret = baryon_2pt_block(wm1, wm2, wm3, 0);
        }
        else if (type == 1)
        {
            const WilsonMatrix wm1 = Cg5 * m1 * Cg5;
            const WilsonMatrix wm2 = m2;
            const WilsonMatrix wm3 = m3 * proj;
            ret = baryon_2pt_block(wm1, wm2, wm3, 1);
        }
        else if (type == -1)
        {
            const WilsonMatrix wm1 = Cg5 * m1 * Cg5;
            const WilsonMatrix wm2 = m2;
            const WilsonMatrix wm3 = m3 * proj;
            ret = baryon_2pt_block(wm1, wm2, wm3, 0) - baryon_2pt_block(wm1, wm2, wm3, 1);
        }
        else
        {
            qassert(false);
        }
        return ret;
    };

    template <class T>
    qacc SpinMatrix proton_twop_no_projection_block(const WilsonMatrixT<T> &m1, const WilsonMatrixT<T> &m2, const WilsonMatrixT<T> &m3, const int type)
    {
        const SpinMatrix &Cg5 = BaryonMatrixConstants::get_Cgamma5();
        const SpinMatrix &proj = BaryonMatrixConstants::get_projector();
        // SpinMatrix block;
        // set_zeros(block);

        if (type == 0)
        {
            const WilsonMatrix wm1 = Cg5 * m1 * Cg5;
            const WilsonMatrix wm2 = m2;
            const WilsonMatrix wm3 = m3;
            return baryon_2pt_no_projection_block(wm1, wm2, wm3, 0);
        }
        else if (type == 1)
        {
            const WilsonMatrix wm1 = Cg5 * m1 * Cg5;
            const WilsonMatrix wm2 = m2;
            const WilsonMatrix wm3 = m3;
            return baryon_2pt_no_projection_block(wm1, wm2, wm3, 1);
        }
        else
        {
            qassert(false);
            return SpinMatrix();
        }
        // return block;
    };

    qacc SpinMatrix baryon_pi_block_no_projection_type(const WilsonMatrix &wm_sq, const WilsonMatrix &wm, const int type)
    {
        SpinMatrix block;
        // set_zeros(block);
        switch (type)
        {
        case 0:
            block = proton_twop_no_projection_block(wm, wm, wm_sq, 0);
            break;

        case 1:
            block = proton_twop_no_projection_block(wm_sq, wm, wm, 0);
            break;

        case 2:
            block = proton_twop_no_projection_block(wm, wm_sq, wm, 1);
            break;

        case 3:
            block = proton_twop_no_projection_block(wm, wm, wm_sq, 1);
            break;

        case 4:
            block = proton_twop_no_projection_block(wm_sq, wm, wm, 1);
            break;

        default:
            qassert(false);
            break;
        }
        return block;
    }

    qacc SpinMatrix baryon_pi_sblock_no_projection_type(const WilsonMatrix &wm_x_sq, const WilsonMatrix &wm_y_sq, const WilsonMatrix &wm, const int type)
    {
        SpinMatrix block;
        // set_zeros(block);
        switch (type)
        {
        case 0:
            block = proton_twop_no_projection_block(wm, wm_y_sq, wm_x_sq, 0);
            break;

        case 1:
            block = proton_twop_no_projection_block(wm_y_sq, wm_x_sq, wm, 0);
            break;

        case 2:
            block = proton_twop_no_projection_block(wm_y_sq, wm_x_sq, wm, 1);
            break;

        case 3:
            block = proton_twop_no_projection_block(wm_y_sq, wm, wm_x_sq, 1);
            break;

        case 4:
            block = proton_twop_no_projection_block(wm, wm_x_sq, wm_y_sq, 1);
            break;

        default:
            qassert(false);
            break;
        }
        return block;
    }

    qacc WilsonMatrix baryon_pi_block_type(const WilsonMatrix &m1, const WilsonMatrix &m2, const int &type)
    {
        const SpinMatrix &Cg5 = BaryonMatrixConstants::get_Cgamma5();
        const SpinMatrix &proj = BaryonMatrixConstants::get_projector();
        WilsonMatrix block;
        WilsonMatrix wm1, wm2;
        set_zero(block);

        switch (type)
        {
        case 0:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = m2;
            block = block_contraction(wm1, wm2, proj, type);
            break;

        case 1:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = m2 * proj;
            block = block_contraction(wm1, wm2, proj, type);
            break;

        case 2:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = m2 * proj;
            block = block_contraction(wm1, wm2, proj, type);
            break;

        case 3:
            wm1 = Cg5 * m1 * Cg5;
            wm2 = proj * m2;
            block = block_contraction(wm1, wm2, proj, type);
            break;

        case 4:
            wm1 = Cg5 * m1 * proj;
            wm2 = m2 * Cg5;
            block = block_contraction(wm1, wm2, proj, type);
            break;

        default:
            qassert(false);
            break;
        }
        return block;
    };

    qacc void baryon_pi_block_pp(array<WilsonMatrix, 5> &block, const WilsonMatrix &m1, const WilsonMatrix &m2)
    {
        for (int type = 0; type < 5; type++)
        {
            block[type] = baryon_pi_block_type(m1, m2, type);
        }
    };

    qacc void baryon_pi_block_sp(array<WilsonMatrix, 5> &block, const WilsonMatrix &m1, const WilsonMatrix &m2, const WilsonMatrix &m3)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        const WilsonMatrix wm = m3 * gamma5 * m2;
        block[0] = baryon_pi_block_type(m1, wm, 0);
        block[1] = baryon_pi_block_type(wm, m1, 1);
        block[2] = baryon_pi_block_type(wm, m1, 2);
        block[3] = baryon_pi_block_type(wm, m1, 3);
        block[4] = baryon_pi_block_type(m1, wm, 2);
    };

    qacc void baryon_pi_block_sp(array<WilsonMatrix, 8 * 5> &block, const WilsonMatrix &m1, const WilsonMatrix &m2, const WilsonMatrix &m3)
    {
        const array<SpinMatrix, 8> &va_ms = get_va_matrices();
        for (int mu = 0; mu < 8; mu++)
        {
            const WilsonMatrix wm = m3 * va_ms[mu] * m2;
            block[mu + 0 * 8] = baryon_pi_block_type(m1, wm, 0);
            block[mu + 1 * 8] = baryon_pi_block_type(wm, m1, 1);
            block[mu + 2 * 8] = baryon_pi_block_type(wm, m1, 2);
            block[mu + 3 * 8] = baryon_pi_block_type(wm, m1, 3);
            block[mu + 4 * 8] = baryon_pi_block_type(m1, wm, 2);
        }
    };
}
