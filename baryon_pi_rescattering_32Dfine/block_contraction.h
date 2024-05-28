#pragma once

#include <qlat/qcd.h>
#include <qlat-utils/complex.h>

namespace qlat
{

    template <class T = Real>
    struct BaryonMatrixConstantsT
    {
        SpinMatrixT<T> Cg5;
        SpinMatrixT<T> Proj_p;

        qacc BaryonMatrixConstantsT() { init(); }
        qacc void init()
        {
            const array<SpinMatrixT<T>, 16> &gms = SpinMatrixConstants::get_cps_gms();
            const SpinMatrixT<T> &gamma5 = SpinMatrixConstants::get_gamma5();
            ComplexD ii(0.0, 1.0);
            Cg5 = gms[10] * gamma5;
            Cg5 *= (ComplexT<T>)ii;

            Proj_p = (gms[0] + gms[8]) / (ComplexT<T>)2;
            // Proj_p_x = (gms[0] + gms[8]) * ((T)ii * gamma5 * gms[1]) / (T)2; // x-direction
            // Proj_p_y = (gms[0] + gms[8]) * ((T)ii * gamma5 * gms[2]) / (T)2; // y-direction
            // Proj_p_z = (gms[0] + gms[8]) * ((T)ii * gamma5 * gms[4]) / (T)2; // z-direction
        }
        //
        static const box<BaryonMatrixConstantsT<T>> &get_instance_box()
        {
            static box<BaryonMatrixConstantsT<T>> smcs =
                box<BaryonMatrixConstantsT<T>>(BaryonMatrixConstantsT<T>());
            return smcs;
        }
        //
        static const BaryonMatrixConstantsT<T> &get_instance()
        {
            return get_instance_box()();
        }
        //
        static const SpinMatrixT<T> &get_Cgamma5()
        {
            return get_instance().Cg5;
        }
        static const SpinMatrixT<T> &get_projector()
        {
            return get_instance().Proj_p;
        }
    };

    template <class T>
    qacc ComplexD meson_2pt_block(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2)
    {
        const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
        WilsonMatrixT<T> wm2_inv_dir = gamma5 * (WilsonMatrix)matrix_adjoint(wm2) * gamma5;
        return matrix_trace(wm1 * wm2_inv_dir);
    };

    template <class T>
    qacc ComplexD pion_2pt_block(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2)
    {
        WilsonMatrixT<T> wm2_inv_dir = (WilsonMatrix)matrix_adjoint(wm2);
        return matrix_trace(wm1 * wm2_inv_dir);
    };

    template <class T>
    qacc ComplexD baryon_2pt_block(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2, const WilsonMatrixT<T> &wm3, const int type)
    {
        ComplexD ret = 0;
        ComplexD temp;
        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }
                                if (type == 0)
                                {
                                    temp = 0;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            temp += wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                        }
                                    }
                                    for (int mu3 = 0; mu3 < 4; mu3++)
                                    {
                                        ret += temp * wm3(mu3 * NUM_COLOR + c1, mu3 * NUM_COLOR + c1p);
                                    }
                                }
                                else if (type == 1)
                                {
                                    temp = 0;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                ret += wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu3 * NUM_COLOR + c1, mu2 * NUM_COLOR + c2p) * wm3(mu1 * NUM_COLOR + c2, mu3 * NUM_COLOR + c1p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    qassert(false);
                                }
                            }
                        }
                    }
                }
            }
        }
        return ret;
    };

    template <class T>
    qacc SpinMatrixT<T> baryon_2pt_no_projection_block(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2, const WilsonMatrixT<T> &wm3, const int type)
    {
        SpinMatrixT<T> block;
        ComplexD temp;
        set_zero(block);

        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }
                                if (type == 0)
                                {
                                    temp = 0;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            temp += wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                        }
                                    }
                                    for (int mu3 = 0; mu3 < 4; mu3++)
                                    {
                                        for (int mu4 = 0; mu4 < 4; mu4++)
                                        {
                                            block(mu4, mu3) += temp * wm3(mu3 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                                        }
                                    }
                                }
                                else if (type == 1)
                                {
                                    temp = 0;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                for (int mu4 = 0; mu4 < 4; mu4++)
                                                {
                                                    block(mu4, mu3) += wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu3 * NUM_COLOR + c1, mu2 * NUM_COLOR + c2p) * wm3(mu1 * NUM_COLOR + c2, mu4 * NUM_COLOR + c1p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    qassert(false);
                                }
                            }
                        }
                    }
                }
            }
        }
        return block;
    };

    template <class T>
    qacc WilsonMatrixT<T> block_contraction(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2, const SpinMatrixT<T> &proj, const int type)
    {
        WilsonMatrixT<T> block;
        ComplexD temp;
        set_zero(block);

        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }
                                if (type == 0)
                                {
                                    temp = 0;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            temp += wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                        }
                                    }
                                    for (int mu3 = 0; mu3 < 4; mu3++)
                                    {
                                        for (int mu4 = 0; mu4 < 4; mu4++)
                                        {
                                            if (proj(mu4, mu3) == 0.0)
                                            {
                                                continue;
                                            }
                                            block(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) += temp * proj(mu4, mu3);
                                        }
                                    }
                                }
                                else if (type == 1)
                                {
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * wm2(mu3 * NUM_COLOR + c1, mu3 * NUM_COLOR + c1p);

                                                block(mu2 * NUM_COLOR + c3p, mu1 * NUM_COLOR + c3) += temp * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                                else if (type == 2)
                                {
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu1 * NUM_COLOR + c2, mu3 * NUM_COLOR + c1p);

                                                block(mu2 * NUM_COLOR + c2p, mu3 * NUM_COLOR + c1) += temp * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                                else if (type == 3)
                                {
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu3 * NUM_COLOR + c1, mu2 * NUM_COLOR + c2p);

                                                block(mu3 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c2) += temp * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                                else if (type == 4)
                                {
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c2, mu3 * NUM_COLOR + c1p) * wm2(mu3 * NUM_COLOR + c1, mu2 * NUM_COLOR + c2p);

                                                block(mu2 * NUM_COLOR + c3p, mu1 * NUM_COLOR + c3) += temp * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    qassert(false);
                                }
                            }
                        }
                    }
                }
            }
        }
        return block;
    };

    qacc void baryon_pi_no_projection_block_0_12(std::vector<WilsonMatrix> &block_0_12, const WilsonMatrix &wm1, const WilsonMatrix &wm2)
    {
        qassert(block_0_12.size() == 4 * 4);

        for (int mu1 = 0; mu1 < 4; mu1++)
        {
            for (int mu2 = 0; mu2 < 4; mu2++)
            {
                set_zero(block_0_12[mu1 * 4 + mu2]);
            }
        }
        WilsonMatrix block_temp;

        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }

                                set_zero(block_temp);

                                for (int mu1 = 0; mu1 < 4; mu1++)
                                {
                                    for (int mu2 = 0; mu2 < 4; mu2++)
                                    {
                                        block_temp(mu2 * NUM_COLOR + c3p, mu1 * NUM_COLOR + c3) += wm1(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                    }
                                }

                                for (int mu3 = 0; mu3 < 4; mu3++)
                                {
                                    for (int mu4 = 0; mu4 < 4; mu4++)
                                    {
                                        block_0_12[mu4 * 4 + mu3] += block_temp * wm2(mu3 * NUM_COLOR + c1, mu4 * NUM_COLOR + c1p);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    qacc void baryon_pi_no_projection_block_0_3(WilsonMatrix &block_0_3, const WilsonMatrix &wm1, const WilsonMatrix &wm2)
    {
        ComplexD temp;
        set_zero(block_0_3);
        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }

                                temp = 0;
                                for (int mu1 = 0; mu1 < 4; mu1++)
                                {
                                    for (int mu2 = 0; mu2 < 4; mu2++)
                                    {
                                        temp += wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                    }
                                }
                                for (int mu3 = 0; mu3 < 4; mu3++)
                                {
                                    for (int mu4 = 0; mu4 < 4; mu4++)
                                    {
                                        block_0_3(mu4 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c1) += temp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    qacc void baryon_pi_no_projection_block_1_1(std::vector<WilsonMatrix> &block_1_1, const WilsonMatrix &wm1, const WilsonMatrix &wm2)
    {
        qassert(block_1_1.size() == 4 * 4);

        for (int mu1 = 0; mu1 < 4; mu1++)
        {
            for (int mu2 = 0; mu2 < 4; mu2++)
            {
                set_zero(block_1_1[mu1 * 4 + mu2]);
            }
        }

        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }

                                for (int mu1 = 0; mu1 < 4; mu1++)
                                {
                                    for (int mu2 = 0; mu2 < 4; mu2++)
                                    {
                                        for (int mu3 = 0; mu3 < 4; mu3++)
                                        {
                                            for (int mu4 = 0; mu4 < 4; mu4++)
                                            {
                                                block_1_1[mu4 * 4 + mu3](mu2 * NUM_COLOR + c3p, mu1 * NUM_COLOR + c3) += wm1(mu3 * NUM_COLOR + c1, mu2 * NUM_COLOR + c2p) * wm2(mu1 * NUM_COLOR + c2, mu4 * NUM_COLOR + c1p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    qacc void baryon_pi_no_projection_block_1_2(std::vector<WilsonMatrix> &block_1_2, const WilsonMatrix &wm1, const WilsonMatrix &wm2)
    {
        ComplexD temp;
        qassert(block_1_2.size() == 4);

        for (int mu1 = 0; mu1 < 4; mu1++)
        {
            set_zero(block_1_2[mu1]);
        }

        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }

                                for (int mu1 = 0; mu1 < 4; mu1++)
                                {
                                    for (int mu2 = 0; mu2 < 4; mu2++)
                                    {
                                        for (int mu3 = 0; mu3 < 4; mu3++)
                                        {
                                            for (int mu4 = 0; mu4 < 4; mu4++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu1 * NUM_COLOR + c2, mu4 * NUM_COLOR + c1p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);

                                                block_1_2[mu4](mu2 * NUM_COLOR + c2p, mu3 * NUM_COLOR + c1) += temp;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    qacc void baryon_pi_no_projection_block_1_3(std::vector<WilsonMatrix> &block_1_3, const WilsonMatrix &wm1, const WilsonMatrix &wm2)
    {
        ComplexD temp;

        qassert(block_1_3.size() == 4);

        for (int mu1 = 0; mu1 < 4; mu1++)
        {
            set_zero(block_1_3[mu1]);
        }
        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor_acc(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor_acc(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }

                                for (int mu1 = 0; mu1 < 4; mu1++)
                                {
                                    for (int mu2 = 0; mu2 < 4; mu2++)
                                    {
                                        for (int mu3 = 0; mu3 < 4; mu3++)
                                        {
                                            for (int mu4 = 0; mu4 < 4; mu4++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c3, mu2 * NUM_COLOR + c3p) * wm2(mu3 * NUM_COLOR + c1, mu2 * NUM_COLOR + c2p) * (ComplexD)epsilon_tensor_acc(c1, c2, c3) * (ComplexD)epsilon_tensor_acc(c1p, c2p, c3p);

                                                block_1_3[mu3](mu4 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c2) += temp;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    };
#ifndef QLAT_NO_DEFAULT_TYPE1

    typedef BaryonMatrixConstantsT<> BaryonMatrixConstants;

#endif
}
