import pytest

import src.effect_size as effect_size


class Test_EffectSize:

    @staticmethod
    def test_cohenes_result() -> None:
        actual = effect_size.cohen_es('r', 'medium')
        # cohen.ES(test="r", size="medium")
        #
        #      Conventional effect size from Cohen (1982)
        #
        #            test = r
        #            size = medium
        #     effect.size = 0.3
        expected = 0.3
        assert actual['effect_size'] == pytest.approx(expected)

    @staticmethod
    def test_esh_result() -> None:
        actual = effect_size.es_h(0.5, 0.4)
        # ES.h(0.5,0.4)
        # [1] 0.2013579
        expected = 0.2013579
        assert actual == pytest.approx(expected)

    @staticmethod
    def test_esw1_result() -> None:
        p0 = [1/4] * 4
        p1 = [3/8] + [(5/8) / 3] * 3
        actual = effect_size.es_w1(p0, p1)
        # P0<-rep(1/4,4)
        # P1<-c(0.375,rep((1-0.375)/3,3))
        # ES.w1(P0,P1)
        # [1] 0.2886751
        expected = 0.2886751
        assert actual == pytest.approx(expected)

    @staticmethod
    def test_esw2_result() -> None:
        prob = [[0.225, 0.125, 0.125, 0.125], [0.160, 0.160, 0.040, 0.040]]
        actual = effect_size.es_w2(prob)
        # > prob<-matrix(c(0.225,0.125,0.125,0.125,0.16,0.16,0.04,0.04),nrow=2,byrow=TRUE)
        # > prob
        #       [,1]  [,2]  [,3]  [,4]
        # [1,] 0.225 0.125 0.125 0.125
        # [2,] 0.160 0.160 0.040 0.040
        # > ES.w2(prob)
        # [1] 0.2558646
        expected = 0.2558646
        assert actual == pytest.approx(expected)
