import pytest

import src.pwr_tests as pwr_tests


class Test_2p:

    @staticmethod
    def test_2p_noNone() -> None:
        with pytest.raises(ValueError, match="One of h, n, sig_level or power must be None"):
            pwr_tests.pwr_2p_test(1, 1, 1, 1)

    @staticmethod
    def test_2p_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of h, n, sig_level or power may be None"):
            pwr_tests.pwr_2p_test(None, 1, None, 1)

    @staticmethod
    def test_2p_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_2p_test(1, 1, 1.5, None)
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_2p_test(1, 1, -1, None)

    @staticmethod
    def test_2p_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_2p_test(None, 1, 0.05, 1.2)
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_2p_test(None, 1, 0.05, -0.5)

    @staticmethod
    def test_2p_results() -> None:
        p_results = pwr_tests.pwr_2p_test(h=0.3, n=200, sig_level=0.05, alternative='greater')
        # pwr.2p.test(h=0.3,n=200, sig.level=0.05, alternative="greater")
        #
        #      Difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.3
        #               n = 200
        #       sig.level = 0.05
        #           power = 0.9123145
        #     alternative = greater
        #
        # NOTE: same sample sizes
        expected = 0.9123145
        assert p_results['power'] == pytest.approx(expected, 0.000001)

        h_results = pwr_tests.pwr_2p_test(n=200, sig_level=0.05, power=0.8, alternative='two-sided')
        # pwr.2p.test(n=200, sig.level=0.05, power=0.8, alternative="two.sided")
        #
        #      Difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.2801491
        #               n = 200
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = two.sided
        #
        # NOTE: same sample sizes
        expected = 0.2801491
        assert h_results['effect_size'] == pytest.approx(expected, 0.0001)

        n_results = pwr_tests.pwr_2p_test(h=-0.3, sig_level=0.05, power=0.8, alternative='less')
        # pwr.2p.test(h=-0.3, sig.level=0.05, power=0.8, alternative="less")
        #
        #      Difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = -0.3
        #               n = 137.3902
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = less
        #
        # NOTE: same sample sizes
        expected = 138
        assert n_results['n'] == expected

        sig_results = pwr_tests.pwr_2p_test(h=0.3, n=200, power=0.8, alternative='two-sided')
        # pwr.2p.test(h=0.3, n=200, sig.level=NULL, power=0.8, alternative="two.sided")
        #
        #      Difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.3
        #               n = 200
        #       sig.level = 0.03089736
        #           power = 0.8
        #     alternative = two.sided
        #
        # NOTE: same sample sizes
        expected = 0.03089736
        assert sig_results['sig_level'] == pytest.approx(expected, 0.0001)


class Test_2p2n:

    @staticmethod
    def test_2p2n_noNone() -> None:
        with pytest.raises(ValueError, match="One of h, n1, n2, sig_level or power must be None"):
            pwr_tests.pwr_2p2n_test(1, 2, 2, 1, 1)

    @staticmethod
    def test_2p2n_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of h, n1, n2, sig_level or power may be None"):
            pwr_tests.pwr_2p2n_test(None, 2, None, 1, 1)

    @staticmethod
    def test_2p2n_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_2p2n_test(1, 2, 2, 1.5, None)
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_2p2n_test(1, 2, 2, -1, None)

    @staticmethod
    def test_2p2n_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_2p2n_test(None, 2, 2, 0.05, 1.2)
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_2p2n_test(None, 2, 2, 0.05, -0.5)

    @staticmethod
    def test_2p2n_result() -> None:
        p_results = pwr_tests.pwr_2p2n_test(h=0.30, n1=80, n2=245, sig_level=0.05, alternative="greater")
        # pwr.2p2n.test(h=0.30,n1=80,n2=245,sig.level=0.05,alternative="greater")
        #
        #      difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.3
        #              n1 = 80
        #              n2 = 245
        #       sig.level = 0.05
        #           power = 0.7532924
        #     alternative = greater
        #
        # NOTE: different sample sizes
        expected = 0.7532924
        assert p_results['power'] == pytest.approx(expected, 0.0000001)

        h_results = pwr_tests.pwr_2p2n_test(n1=1600, n2=3000, power=0.8, sig_level=0.1, alternative="less")
        # pwr.2p2n.test(n1=1600,n2=3000,power=0.8,sig.level=0.1,alternative="less")
        #
        #      difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = -0.06572994
        #              n1 = 1600
        #              n2 = 3000
        #       sig.level = 0.1
        #           power = 0.8
        #     alternative = less
        #
        # NOTE: different sample sizes
        expected = -0.06572994
        assert h_results['effect_size'] == pytest.approx(expected, 0.0001)

        n1_result = pwr_tests.pwr_2p2n_test(h=0.3, n2=1000, power=0.8, sig_level=0.05, alternative='two-sided')
        # pwr.2p2n.test(h=0.3, n2=1000, power=0.8, sig.level=0.05, alternative='two.sided')
        #
        #      difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.3
        #              n1 = 95.54171
        #              n2 = 1000
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = two.sided
        #
        # NOTE: different sample sizes
        expected = 96
        assert n1_result['n1'] == expected

        n2_results = pwr_tests.pwr_2p2n_test(h=0.20, n1=1600, power=0.9, sig_level=0.01, alternative="two-sided")
        # pwr.2p2n.test(h=0.20,n1=1600,power=0.9,sig.level=0.01,alternative="two.sided")
        #
        #      difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.2
        #              n1 = 1600
        #              n2 = 484.6646
        #       sig.level = 0.01
        #           power = 0.9
        #     alternative = two.sided
        #
        # NOTE: different sample sizes
        expected = 485
        assert n2_results['n2'] == expected

        s_results = pwr_tests.pwr_2p2n_test(h=0.3, n1=100, n2=1000, power=0.8, alternative='two-sided')
        # pwr.2p2n.test(h=0.3, n1=500, n2=1000, power=0.8, sig.level=NULL, alternative='two.sided')
        #
        #      difference of proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.3
        #              n1 = 100
        #              n2 = 1000
        #       sig.level = 0.04352786
        #           power = 0.8
        #     alternative = two.sided
        #
        # NOTE: different sample sizes
        expected = 0.04352786
        assert s_results['sig_level'] == pytest.approx(expected, 0.001)


class Test_Anova:
    @staticmethod
    def test_anova_noNone() -> None:
        with pytest.raises(ValueError, match="One of k, n, f, sig_level or power must be None"):
            pwr_tests.pwr_anova_test(1, 2, 2, 0.05, 0.8)

    @staticmethod
    def test_anova_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of k, n, f, sig_level or power may be None"):
            pwr_tests.pwr_anova_test(None, 2, None, 0.05, 0.8)

    @staticmethod
    def test_anova_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_anova_test(2, 2, 2, 1.5, None)
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_anova_test(2, 2, 2, -1, None)

    @staticmethod
    def test_anova_f() -> None:
        with pytest.raises(ValueError, match="f must be positive"):
            pwr_tests.pwr_anova_test(2, 2, -0.5, 0.5, None)

    @staticmethod
    def test_anova_k() -> None:
        with pytest.raises(ValueError, match="Number of groups must be at least 2"):
            pwr_tests.pwr_anova_test(1, 2, 3, None, 0.8)

    @staticmethod
    def test_anova_n() -> None:
        with pytest.raises(ValueError, match="Number of observations must be at least 2"):
            pwr_tests.pwr_anova_test(2, 1, None, 0.05, 0.9)

    @staticmethod
    def test_anova_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_anova_test(None, 2, 2, 0.05, 1.2)
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_anova_test(None, 2, 2, 0.05, -0.5)

    @staticmethod
    def test_anova_results() -> None:
        p_results = pwr_tests.pwr_anova_test(f=0.28, k=4, n=20, sig_level=0.05)
        # pwr.anova.test(f=0.28,k=4,n=20,sig.level=0.05)
        #
        #      Balanced one-way analysis of variance power calculation
        #
        #               k = 4
        #               n = 20
        #               f = 0.28
        #       sig.level = 0.05
        #           power = 0.5149793
        #
        # NOTE: n is number in each group
        expected = 0.5149793
        assert p_results['power'] == pytest.approx(expected, 0.00001)

        k_results = pwr_tests.pwr_anova_test(f=0.1, n=50, power=0.80, sig_level=0.05)
        # pwr.anova.test(f=0.1, n=50, power=0.80, sig.level=0.05)
        #
        #      Balanced one-way analysis of variance power calculation
        #
        #               k = 70.08744
        #               n = 50
        #               f = 0.1
        #       sig.level = 0.05
        #           power = 0.8
        #
        # NOTE: n is number in each group
        expected = 71
        assert k_results['k'] == expected

        n_results = pwr_tests.pwr_anova_test(f=0.28, k=4, power=0.80, sig_level=0.05)
        # pwr.anova.test(f=0.28,k=4,power=0.80,sig.level=0.05)
        #
        #      Balanced one-way analysis of variance power calculation
        #
        #               k = 4
        #               n = 35.75789
        #               f = 0.28
        #       sig.level = 0.05
        #           power = 0.8
        #
        # NOTE: n is number in each group
        expected = 36
        assert n_results['n'] == expected

        f_results = pwr_tests.pwr_anova_test(k=5, n=10, power=0.80, sig_level=0.05)
        # pwr.anova.test(k=5, n=10, power=0.80, sig.level=0.05)
        #
        #      Balanced one-way analysis of variance power calculation
        #
        #               k = 5
        #               n = 10
        #               f = 0.5148773
        #       sig.level = 0.05
        #           power = 0.8
        #
        # NOTE: n is number in each group
        expected = 0.5148773
        assert f_results['effect_size'] == pytest.approx(expected, 0.0001)

        s_results = pwr_tests.pwr_anova_test(k=3, n=20, f=0.5, power=0.8)
        # pwr.anova.test(k=3, n=20, f=0.5, power=0.8, sig.level=NULL)
        #
        #      Balanced one-way analysis of variance power calculation
        #
        #               k = 3
        #               n = 20
        #               f = 0.5
        #       sig.level = 0.01023897
        #           power = 0.8
        #
        # NOTE: n is number in each group
        expected = 0.01023897
        assert s_results['sig_level'] == pytest.approx(expected, 0.002)


class Test_Chisq:
    @staticmethod
    def test_chisq_noNone() -> None:
        with pytest.raises(ValueError, match="One of w, n, sig_level or power must be None"):
            pwr_tests.pwr_chisq_test(1, 2, 2, 0.05, 0.8)

    @staticmethod
    def test_chisq_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of w, n, sig_level or power may be None"):
            pwr_tests.pwr_chisq_test(None, None, 2, 0.05, 0.8)

    @staticmethod
    def test_chisq_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_chisq_test(2, 2, 2, 1.5, None)
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_chisq_test(2, 2, 2, -1, None)

    @staticmethod
    def test_chisq_w() -> None:
        with pytest.raises(ValueError, match="w must be positive"):
            pwr_tests.pwr_chisq_test(-1, 2, 2, 0.5, None)

    @staticmethod
    def test_chisq_n() -> None:
        with pytest.raises(ValueError, match="Number of observations must be at least 1"):
            pwr_tests.pwr_chisq_test(0.3, 0, 3, None, 0.8)

    @staticmethod
    def test_chisq_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_chisq_test(None, 2, 2, 0.05, 1.2)
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_chisq_test(None, 2, 2, 0.05, -0.5)

    @staticmethod
    def test_chisq_results() -> None:
        p_results = pwr_tests.pwr_chisq_test(w=0.289, df=3, n=100, sig_level=0.05)
        # pwr.chisq.test(w=0.289,df=(4-1),N=100,sig.level=0.05)
        #
        #      Chi squared power calculation
        #
        #               w = 0.289
        #               N = 100
        #              df = 3
        #       sig.level = 0.05
        #           power = 0.6750777
        #
        # NOTE: N is the number of observations
        expected = 0.6750777
        assert p_results['power'] == pytest.approx(expected, 0.00001)

        w_results = pwr_tests.pwr_chisq_test(n=300, df=30, power=0.80, sig_level=0.05)
        # pwr.chisq.test(N=300, df=30, power=0.80, sig.level=0.05)
        #
        #      Chi squared power calculation
        #
        #               w = 0.2860569
        #               N = 300
        #              df = 30
        #       sig.level = 0.05
        #           power = 0.8
        #
        # NOTE: N is the number of observations
        expected = 0.2860569
        assert w_results['effect_size'] == pytest.approx(expected, 0.0002)

        n_results = pwr_tests.pwr_chisq_test(w=0.1, df=(5-1)*(6-1), power=0.80, sig_level=0.05)
        # pwr.chisq.test(w=0.1,df=(5-1)*(6-1),power=0.80,sig.level=0.05)
        #
        #      Chi squared power calculation
        #
        #               w = 0.1
        #               N = 2096.079
        #              df = 20
        #       sig.level = 0.05
        #           power = 0.8
        #
        # NOTE: N is the number of observations
        expected = 2097
        assert n_results['n'] == expected

        s_results = pwr_tests.pwr_chisq_test(w=0.25, n=300, df=30, power=0.80)
        # pwr.chisq.test(w=0.25, N=300, df=30, power=0.80, sig.level = NULL)
        #
        #      Chi squared power calculation
        #
        #               w = 0.25
        #               N = 300
        #              df = 30
        #       sig.level = 0.1304078
        #           power = 0.8
        #
        # NOTE: N is the number of observations
        expected = 0.1304078
        assert s_results['sig_level'] == pytest.approx(expected, 0.0001)


class Test_f2:
    @staticmethod
    def test_f2_noNone() -> None:
        with pytest.raises(ValueError, match="One of u, v, f2, sig_level or power must be None"):
            pwr_tests.pwr_f2_test(1, 2, 2, 0.05, 0.8)

    @staticmethod
    def test_f2_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of u, v, f2, sig_level or power may be None"):
            pwr_tests.pwr_f2_test(None, None, 2, 0.05, 0.8)

    @staticmethod
    def test_f2_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_f2_test(2, 2, 2, 1.5, None)
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_f2_test(2, 2, 2, -1, None)

    @staticmethod
    def test_f2_u() -> None:
        with pytest.raises(ValueError, match="Degrees of freedom u for numerator must be at least 1"):
            pwr_tests.pwr_f2_test(0, 2, 2, 0.5, None)

    @staticmethod
    def test_f2_v() -> None:
        with pytest.raises(ValueError, match="Degrees of freedom v for denominator must be at least 1"):
            pwr_tests.pwr_f2_test(2, 0, 3, None, 0.8)

    @staticmethod
    def test_f2_f2() -> None:
        with pytest.raises(ValueError, match="f2 must be positive"):
            pwr_tests.pwr_f2_test(2, None, -2, 0.05, 0.8)

    @staticmethod
    def test_f2_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_f2_test(None, 2, 2, 0.05, 1.2)
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_f2_test(None, 2, 2, 0.05, -0.5)

    @staticmethod
    def test_f2_results() -> None:
        p_result = pwr_tests.pwr_f2_test(u=5, v=89, f2=0.1/(1-0.1), sig_level=0.05)
        # pwr.f2.test(u=5,v=89,f2=0.1/(1-0.1),sig.level=0.05)
        #
        #      Multiple regression power calculation
        #
        #               u = 5
        #               v = 89
        #              f2 = 0.1111111
        #       sig.level = 0.05
        #           power = 0.6735858
        expected = 0.6735858
        assert p_result['power'] == pytest.approx(expected, 0.00001)

        u_result = pwr_tests.pwr_f2_test(v=90, f2=0.3, sig_level=0.05, power=0.8)
        # pwr.f2.test(v=90,f2=0.3,sig.level=0.05, power=0.8)
        #
        #      Multiple regression power calculation
        #
        #               u = 55.48244
        #               v = 90
        #              f2 = 0.3
        #       sig.level = 0.05
        #           power = 0.8
        expected = 56
        assert u_result['u'] == expected

        v_result = pwr_tests.pwr_f2_test(u=90, f2=0.01, sig_level=0.05, power=0.8)
        # pwr.f2.test(u=90,f2=0.01,sig.level=0.05, power=0.8)
        #
        #      Multiple regression power calculation
        #
        #               u = 90
        #               v = 3840.388
        #              f2 = 0.01
        #       sig.level = 0.05
        #           power = 0.8
        expected = 3841
        assert v_result['v'] == expected

        f2_result = pwr_tests.pwr_f2_test(u=100, v=1000, sig_level=0.1, power=0.8)
        # pwr.f2.test(u=100, v=1000, sig.level=0.1, power=0.8)
        #
        #      Multiple regression power calculation
        #
        #               u = 100
        #               v = 1000
        #              f2 = 0.03279811
        #       sig.level = 0.1
        #           power = 0.8
        expected = 0.03279811
        assert f2_result['effect_size'] == pytest.approx(expected, 0.001)

        s_result = pwr_tests.pwr_f2_test(f2=0.15, u=100, v=130, power=0.8)
        # pwr.f2.test(f=0.15, u=100, v=130, power=0.8, sig.level = NULL)
        #
        #      Multiple regression power calculation
        #
        #               u = 100
        #               v = 130
        #              f2 = 0.15
        #       sig.level = 0.2253544
        #           power = 0.8
        expected = 0.2253544
        assert s_result['sig_level'] == pytest.approx(expected, 0.0001)


class Test_Norm:
    @staticmethod
    def test_norm_noNone() -> None:
        with pytest.raises(ValueError, match="One of d, n, sig_level or power must be None"):
            pwr_tests.pwr_norm_test(1, 1, 1, 1, 'two-sided')

    @staticmethod
    def test_norm_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of d, n, sig_level or power may be None"):
            pwr_tests.pwr_norm_test(None, None, 2, 0.05, 'less')

    @staticmethod
    def test_norm_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_norm_test(2, 2, 1.5, None, 'two-sided')
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_norm_test(2, None, -1, 0.5, 'greater')

    @staticmethod
    def test_norm_n() -> None:
        with pytest.raises(ValueError, match="Number of observations in each group must be at least 1"):
            pwr_tests.pwr_norm_test(None, 0, 0.8, 0.5, 'less')

    @staticmethod
    def test_f2_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_norm_test(None, 2, 0.05, 1.2, 'greater')
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_norm_test(None, 2, 0.05, -0.5, 'two-sided')

    @staticmethod
    def test_norm_results() -> None:
        p_results = pwr_tests.pwr_norm_test(d=1/3, n=20, sig_level=0.05, alternative='greater')
        # pwr.norm.test(d=d,n=20,sig.level=0.05,alternative="greater")
        #
        #      Mean power calculation for normal distribution with known variance
        #
        #               d = 0.3333333
        #               n = 20
        #       sig.level = 0.05
        #           power = 0.438749
        #     alternative = greater
        expected = 0.438749
        assert p_results['power'] == pytest.approx(expected, 0.00001)

        d_results = pwr_tests.pwr_norm_test(n=30, power=0.8, sig_level=0.05, alternative='two-sided')
        # pwr.norm.test(n=30, power=0.8, sig.level=0.05, alternative='two.sided')
        #
        #      Mean power calculation for normal distribution with known variance
        #
        #               d = 0.5114965
        #               n = 30
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = two.sided
        expected = 0.5114965
        assert d_results['effect_size'] == pytest.approx(expected, 0.00001)

        n_results = pwr_tests.pwr_norm_test(d=1/3, power=0.8, sig_level=0.05, alternative='greater')
        # pwr.norm.test(d=1/3,power=0.8,sig.level=0.05,alternative="greater")
        #
        #      Mean power calculation for normal distribution with known variance
        #
        #               d = 0.3333333
        #               n = 55.64301
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = greater
        expected = 56
        assert n_results['n'] == expected

        s_results = pwr_tests.pwr_norm_test(d=0.15, n=20, power=0.8, alternative='less')
        # pwr.norm.test(d=0.15, n=20, power=0.8, alternative='less', sig.level = NULL)
        #
        #      Mean power calculation for normal distribution with known variance
        #
        #               d = 0.15
        #               n = 20
        #       sig.level = 0.934789
        #           power = 0.8
        #     alternative = less
        expected = 0.934789
        assert s_results['sig_level'] == pytest.approx(expected, 0.0001)


class Test_P:
    @staticmethod
    def test_p_noNone() -> None:
        with pytest.raises(ValueError, match="One of h, n, sig_level or power must be None"):
            pwr_tests.pwr_p_test(1, 1, 1, 1, 'two-sided')

    @staticmethod
    def test_p_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of h, n, sig_level or power may be None"):
            pwr_tests.pwr_p_test(None, None, 2, 0.05, 'less')

    @staticmethod
    def test_p_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_p_test(2, 2, 1.5, None, 'two-sided')
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_p_test(2, None, -1, 0.5, 'greater')

    @staticmethod
    def test_p_n() -> None:
        with pytest.raises(ValueError, match="Number of observations in each group must be at least 1"):
            pwr_tests.pwr_p_test(None, 0, 0.8, 0.5, 'less')

    @staticmethod
    def test_p_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_p_test(None, 2, 0.05, 1.2, 'greater')
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_p_test(None, 2, 0.05, -0.5, 'two-sided')

    @staticmethod
    def test_p_result() -> None:
        p_result = pwr_tests.pwr_p_test(h=0.2013579, n=60, sig_level=0.05, alternative='two-sided')
        # pwr.p.test(h=h,n=60,sig.level=0.05,alternative="two.sided")
        #
        #      proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.2013579
        #               n = 60
        #       sig.level = 0.05
        #           power = 0.3447014
        #     alternative = two.sided
        expected = 0.3447014
        assert p_result['power'] == pytest.approx(expected, 0.000001)

        h_result = pwr_tests.pwr_p_test(n=200, power=0.80, sig_level=0.1, alternative="less")
        # pwr.p.test(n = 200,power=0.80,sig.level=0.1,alternative="less")
        #
        #      proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = -0.150142
        #               n = 200
        #       sig.level = 0.1
        #           power = 0.8
        #     alternative = less
        expected = -0.150142
        assert h_result['effect_size'] == pytest.approx(expected, 0.0001)

        n_result = pwr_tests.pwr_p_test(h=0.2, power=0.95, sig_level=0.05, alternative="two-sided")
        # pwr.p.test(h=0.2,power=0.95,sig.level=0.05,alternative="two.sided")
        #
        #      proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.2
        #               n = 324.8677
        #       sig.level = 0.05
        #           power = 0.95
        #     alternative = two.sided
        expected = 325
        assert n_result['n'] == expected

        s_result = pwr_tests.pwr_p_test(h=0.2, n=150, power=0.95, alternative="greater")
        # pwr.p.test(h=0.2, n=150, power=0.95, alternative="greater", sig.level = NULL)
        #
        #      proportion power calculation for binomial distribution (arcsine transformation)
        #
        #               h = 0.2
        #               n = 150
        #       sig.level = 0.2105159
        #           power = 0.95
        #     alternative = greater
        expected = 0.2105159
        assert s_result['sig_level'] == pytest.approx(expected, 0.0001)


class Test_R:
    @staticmethod
    def test_r_noNone() -> None:
        with pytest.raises(ValueError, match="One of r, n, sig_level or power must be None"):
            pwr_tests.pwr_r_test(1, 1, 1, 1, 'two-sided')

    @staticmethod
    def test_r_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of r, n, sig_level or power may be None"):
            pwr_tests.pwr_r_test(None, None, 2, 0.05, 'less')

    @staticmethod
    def test_r_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_r_test(4, 0.5, 1.5, None, 'two-sided')
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_r_test(4, None, -1, 0.5, 'greater')

    @staticmethod
    def test_r_n() -> None:
        with pytest.raises(ValueError, match="Number of observations must be at least 4"):
            pwr_tests.pwr_r_test(3, None, 0.8, 0.5, 'less')

    @staticmethod
    def test_r_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_r_test(None, 0.5, 0.05, 1.2, 'greater')
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_r_test(None, 0.5, 0.05, -0.5, 'two-sided')

    @staticmethod
    def test_r_result() -> None:
        p_result = pwr_tests.pwr_r_test(r=0.3, n=50, sig_level=0.05, alternative="greater")
        # pwr.r.test(r=0.3,n=50,sig.level=0.05,alternative="greater")
        #
        #      approximate correlation power calculation (arctangh transformation)
        #
        #               n = 50
        #               r = 0.3
        #       sig.level = 0.05
        #           power = 0.6911395
        #     alternative = greater
        expected = 0.6911395
        assert p_result['power'] == pytest.approx(expected, 0.00001)

        r_result = pwr_tests.pwr_r_test(n=125, power=0.8, sig_level=0.1, alternative="less")
        # pwr.r.test(n=125, power=0.8, sig.level=0.1, alternative="less")
        #
        #      approximate correlation power calculation (arctangh transformation)
        #
        #               n = 125
        #               r = -0.1890504
        #       sig.level = 0.1
        #           power = 0.8
        #     alternative = less
        expected = -0.1890504
        assert r_result['effect_size'] == pytest.approx(expected, 0.00001)

        n_result = pwr_tests.pwr_r_test(r=0.3, power=0.80, sig_level=0.05, alternative="two-sided")
        # pwr.r.test(r=0.3,power=0.80,sig.level=0.05,alternative="two.sided")
        #
        #      approximate correlation power calculation (arctangh transformation)
        #
        #               n = 84.07364
        #               r = 0.3
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = two.sided
        expected = 85
        assert n_result['n'] == expected

        s_result = pwr_tests.pwr_r_test(r=0.3, n=125, power=0.8, alternative="two-sided")
        # pwr.r.test(r=0.3, n=125, power=0.8, sig.level=NULL, alternative="two.sided")
        #
        #      approximate correlation power calculation (arctangh transformation)
        #
        #               n = 125
        #               r = 0.3
        #       sig.level = 0.009736855
        #           power = 0.8
        #     alternative = two.sided
        expected = 0.009736855
        assert s_result['sig_level'] == pytest.approx(expected, abs=0.0001)


class Test_T:
    @staticmethod
    def test_t_noNone() -> None:
        with pytest.raises(ValueError, match="One of n, d, sig_level or power must be None"):
            pwr_tests.pwr_t_test(1, 1, 1, 1, 'paired', 'two-sided')

    @staticmethod
    def test_t_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of n, d, sig_level or power may be None"):
            pwr_tests.pwr_t_test(None, None, 2, 0.05, 'one', 'less')

    @staticmethod
    def test_t_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_t_test(4, 0.5, 1.5, None, 'two', 'two-sided')
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_t_test(4, None, -1, 0.5, 'paired', 'greater')

    @staticmethod
    def test_t_n() -> None:
        with pytest.raises(ValueError, match="Number of observations must be at least 2"):
            pwr_tests.pwr_t_test(1, None, 0.8, 0.5, 'paired', 'less')

    @staticmethod
    def test_t_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_t_test(None, 0.5, 0.05, 1.2, 'one', 'greater')
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_t_test(None, 0.5, 0.05, -0.5, 'two', 'two-sided')

    @staticmethod
    def test_t_results() -> None:
        p_results = pwr_tests.pwr_t_test(d=0.2, n=60, sig_level=0.10, test_type="one-sample", alternative="two-sided")
        # pwr.t.test(d=0.2,n=60,sig.level=0.10,type="one.sample",alternative="two.sided")
        #
        #      One-sample t test power calculation
        #
        #               n = 60
        #               d = 0.2
        #       sig.level = 0.1
        #           power = 0.4555818
        #     alternative = two.sided
        expected = 0.4555818
        assert p_results['power'] == pytest.approx(expected, 0.0001)

        d_results = pwr_tests.pwr_t_test(n=100, power=0.8, sig_level=0.05, test_type="paired", alternative="greater")
        # pwr.t.test(n=100,power=0.8,sig.level=0.05,type="paired",alternative="greater")
        #
        #      Paired t test power calculation
        #
        #               n = 100
        #               d = 0.2503641
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = greater
        #
        # NOTE: n is number of *pairs*
        expected = 0.2503641
        assert d_results['effect_size'] == pytest.approx(expected, 0.0001)

        n_results = pwr_tests.pwr_t_test(d=0.3, power=0.75, sig_level=0.05, test_type='two-sample', alternative='greater')
        # pwr.t.test(d=0.3,power=0.75,sig.level=0.05,type="two.sample",alternative="greater")
        #
        #      Two-sample t test power calculation
        #
        #               n = 120.2232
        #               d = 0.3
        #       sig.level = 0.05
        #           power = 0.75
        #     alternative = greater
        #
        # NOTE: n is number in *each* group
        expected = 121
        assert n_results['n'] == expected

        s_results = pwr_tests.pwr_t_test(d=-0.1, power=0.75, n=100, test_type='paired', alternative='less')
        # pwr.t.test(d=-0.1, power=0.75, n=100, type='paired', alternative='less', sig.level = NULL)
        #
        #      Paired t test power calculation
        #
        #               n = 100
        #               d = -0.1
        #       sig.level = 0.3725015
        #           power = 0.75
        #     alternative = less
        #
        # NOTE: n is number of *pairs*
        expected = 0.3725015
        assert s_results['sig_level'] == pytest.approx(expected, 0.001)


class Test_T2N:
    @staticmethod
    def test_t2n_noNone() -> None:
        with pytest.raises(ValueError, match="One of n1, n2, d sig_level or power must be None"):
            pwr_tests.pwr_t2n_test(10, 5, 0.5, 0.05, 0.8, 'two-sided')

    @staticmethod
    def test_t2n_multipleNone() -> None:
        with pytest.raises(ValueError, match="Only one of n1, n2, d, sig_level or power may be None"):
            pwr_tests.pwr_t2n_test(None, None, 0.3, 0.05, 0.9, 'less')

    @staticmethod
    def test_t2n_sigLevel() -> None:
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_t2n_test(4, 5, 0.5, 1.5, None, 'two-sided')
        with pytest.raises(ValueError, match="sig_level must be between 0 and 1"):
            pwr_tests.pwr_t2n_test(4, 7, None, -1, 0.5, 'greater')

    @staticmethod
    def test_t2n_n1() -> None:
        with pytest.raises(ValueError, match="Number of observations in the first group must be at least 2"):
            pwr_tests.pwr_t2n_test(1, None, 0.8, 0.5, 0.8, 'less')

    @staticmethod
    def test_t2n_n2() -> None:
        with pytest.raises(ValueError, match="Number of observations in the second group must be at least 2"):
            pwr_tests.pwr_t2n_test(None, 1, 0.8, 0.5, 0.8, 'less')

    @staticmethod
    def test_t2n_power() -> None:
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_t2n_test(None, 10, 0.5, 0.05, 1.2, 'greater')
        with pytest.raises(ValueError, match="power must be between 0 and 1"):
            pwr_tests.pwr_t2n_test(None, 15, 0.5, 0.05, -0.5, 'two-sided')

    @staticmethod
    def test_t2n_results() -> None:
        p_results = pwr_tests.pwr_t2n_test(d=0.6, n1=90, n2=60, sig_level=0.05, alternative="greater")
        # pwr.t2n.test(d=0.6,n1=90,n2=60,alternative="greater")
        #
        #      t test power calculation
        #
        #              n1 = 90
        #              n2 = 60
        #               d = 0.6
        #       sig.level = 0.05
        #           power = 0.9737262
        #     alternative = greater
        expected = 0.9737262
        assert p_results['power'] == pytest.approx(expected, 0.0001)

        d_results = pwr_tests.pwr_t2n_test(n1=85, n2=100, sig_level=0.1, power=0.9, alternative='less')
        # pwr.t2n.test(n1=85, n2=100, sig.level=0.1, alternative='less', power=0.9)
        #
        #      t test power calculation
        #
        #              n1 = 85
        #              n2 = 100
        #               d = -0.3789791
        #       sig.level = 0.1
        #           power = 0.9
        #     alternative = less
        expected = -0.3789791
        assert d_results['effect_size'] == pytest.approx(expected, 0.0001)

        n1_results = pwr_tests.pwr_t2n_test(n2=90, sig_level=0.05, d=0.5, power=0.8, alternative='two-sided')
        # pwr.t2n.test(90, sig.level=0.05, d=0.5, power=0.8, alternative='two.sided')
        #
        #      t test power calculation
        #
        #              n1 = 90
        #              n2 = 49.27233
        #               d = 0.5
        #       sig.level = 0.05
        #           power = 0.8
        #     alternative = two.sided
        expected = 50
        assert n1_results['n1'] == expected

        n2_results = pwr_tests.pwr_t2n_test(n1=1000, sig_level=0.05, d=0.5, power=0.9, alternative='two-sided')
        # pwr.t2n.test(n1=1000, sig.level=0.05, d=0.5, power=0.9, alternative='two.sided')
        #
        #      t test power calculation
        #
        #              n1 = 1000
        #              n2 = 43.95826
        #               d = 0.5
        #       sig.level = 0.05
        #           power = 0.9
        #     alternative = two.sided
        expected = 44
        assert n2_results['n2'] == expected

        s_results = pwr_tests.pwr_t2n_test(n1=100, n2=200, d=0.2, power=0.8, alternative='greater')
        # pwr.t2n.test(n1=100, n2=200, d=0.2, power=0.8, alternative='greater', sig.level = NULL)
        #
        #      t test power calculation
        #
        #              n1 = 100
        #              n2 = 200
        #               d = 0.2
        #       sig.level = 0.2146133
        #           power = 0.8
        #     alternative = greater
        expected = 0.2146133
        assert s_results['sig_level'] == pytest.approx(expected, 0.0001)
