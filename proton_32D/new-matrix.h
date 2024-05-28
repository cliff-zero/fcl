#pragma once

#include <qlat/qcd.h>

namespace qlat
{
    struct EpisilonTensorIndex3
    {
        int tensoridx[6][3];
        EpisilonTensorIndex3() { init(); }
        void init()
        {
            int i = 0;
            for (int c1 = 0; c1 < NUM_COLOR; c1++)
            {
                for (int j = c1; j < c1 + 2; j++)
                {
                    int c2 = j % 3;
                    int c3 = 3 - c1 - c2;
                    tensoridx[i][0] = c1;
                    tensoridx[i][1] = c2;
                    tensoridx[i][2] = c3;
                }
            }
        }
    };

    template <class T = ComplexT>
    struct ProtonMatrixConstantsT
    {
        SpinMatrixT<T> Cg5;
        SpinMatrixT<T> Proj_p;
        qacc ProtonMatrixConstantsT() { init(); }
        qacc void init()
        {
            const array<SpinMatrix, 16> &gms = SpinMatrixConstants::get_gms();
            const SpinMatrix &gamma5 = SpinMatrixConstants::get_gamma5();
            Complex ii(0.0, 1.0);
            Cg5 = gms[10] * gamma5;
            Cg5 *= (T)ii;
            Proj_p = (gms[0] + gms[8]) / (T)2;
        }
        //
        static const box<ProtonMatrixConstantsT<T>> &get_instance_box()
        {
            static box<ProtonMatrixConstantsT<T>> smcs =
                box<ProtonMatrixConstantsT<T>>(ProtonMatrixConstantsT<T>());
            return smcs;
        }
        //
        static const ProtonMatrixConstantsT<T> &get_instance()
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
    qacc Complex lv_wm3(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2, const WilsonMatrixT<T> &wm3, const int type)
    {
        Complex ret = 0;
        Complex temp;
        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor(c1p, c2p, c3p) == 0)
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
                                            temp += wm1(mu1 * NUM_COLOR + c1, mu2 * NUM_COLOR + c1p) * wm2(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
                                        }
                                    }
                                    for (int mu3 = 0; mu3 < 4; mu3++)
                                    {
                                        ret += temp * wm3(mu3 * NUM_COLOR + c3, mu3 * NUM_COLOR + c3p);
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
                                                ret += wm1(mu3 * NUM_COLOR + c3, mu2 * NUM_COLOR + c1p) * wm2(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * wm3(mu1 * NUM_COLOR + c1, mu3 * NUM_COLOR + c3p) * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
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
    qacc WilsonMatrixT<T> lv_wm2(const WilsonMatrixT<T> &wm1, const WilsonMatrixT<T> &wm2, const SpinMatrixT<T> &proj, const int type)
    {
        WilsonMatrixT<T> block;
        Complex temp;
        set_zero(block);
        for (int c1 = 0; c1 < NUM_COLOR; c1++)
        {
            for (int c2 = 0; c2 < NUM_COLOR; c2++)
            {
                for (int c3 = 0; c3 < NUM_COLOR; c3++)
                {
                    if (epsilon_tensor(c1, c2, c3) == 0)
                    {
                        continue;
                    }
                    for (int c1p = 0; c1p < NUM_COLOR; c1p++)
                    {
                        for (int c2p = 0; c2p < NUM_COLOR; c2p++)
                        {
                            for (int c3p = 0; c3p < NUM_COLOR; c3p++)
                            {
                                if (epsilon_tensor(c1p, c2p, c3p) == 0)
                                {
                                    continue;
                                }
                                if (type == 1)
                                {
                                    temp = 0;
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            temp += wm1(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * wm2(mu1 * NUM_COLOR + c1, mu2 * NUM_COLOR + c1p) * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
                                        }
                                    }
                                    for (int mu3 = 0; mu3 < 4; mu3++)
                                    {
                                        for (int mu4 = 0; mu4 < 4; mu4++)
                                        {
                                            if (proj(mu3, mu4) == 0.0)
                                            {
                                                continue;
                                            }
                                            block(mu4 * NUM_COLOR + c3p, mu3 * NUM_COLOR + c3) += temp * proj(mu4, mu3);
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
                                                temp = wm1(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * wm2(mu3 * NUM_COLOR + c3, mu3 * NUM_COLOR + c3p);
                                                block(mu2 * NUM_COLOR + c1p, mu1 * NUM_COLOR + c1) += temp * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
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
                                                temp = wm1(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * wm2(mu3 * NUM_COLOR + c3, mu2 * NUM_COLOR + c1p);
                                                block(mu3 * NUM_COLOR + c3p, mu1 * NUM_COLOR + c1) += temp * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
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
                                                temp = wm1(mu1 * NUM_COLOR + c2, mu2 * NUM_COLOR + c2p) * wm2(mu1 * NUM_COLOR + c1, mu3 * NUM_COLOR + c3p);
                                                block(mu2 * NUM_COLOR + c1p, mu3 * NUM_COLOR + c3) += temp * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
                                            }
                                        }
                                    }
                                }
                                else if (type == 5)
                                {
                                    for (int mu1 = 0; mu1 < 4; mu1++)
                                    {
                                        for (int mu2 = 0; mu2 < 4; mu2++)
                                        {
                                            for (int mu3 = 0; mu3 < 4; mu3++)
                                            {
                                                temp = wm1(mu1 * NUM_COLOR + c1, mu3 * NUM_COLOR + c3p) * wm2(mu3 * NUM_COLOR + c3, mu2 * NUM_COLOR + c1p);
                                                block(mu2 * NUM_COLOR + c2p, mu1 * NUM_COLOR + c2) += temp * (Complex)epsilon_tensor(c1, c2, c3) * (Complex)epsilon_tensor(c1p, c2p, c3p);
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

#ifndef QLAT_NO_DEFAULT_TYPE1

    typedef ProtonMatrixConstantsT<> ProtonMatrixConstants;

#endif
}
