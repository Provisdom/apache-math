(ns provisdom.apache-math.alternative-random-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.apache-math.alternative-random :as alt-random]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]))

;;? seconds

(set! *warn-on-reflection* true)

(ost/instrument)

;;;APACHE RANDOM NUMBER GENERATORS
(deftest quasi-rnd-vector-lazy-test
  (is (spec-check alt-random/quasi-rnd-vector-lazy))
  (is= '([0.0] [0.5] [0.75])
       (take 3 (alt-random/quasi-rnd-vector-lazy 1)))
  (is= '([0.0 0.0] [0.5 0.5] [0.75 0.25])
       (take 3 (alt-random/quasi-rnd-vector-lazy 2))))

(deftest secure-rnd-lazy-test
  (is (spec-check alt-random/secure-rnd-lazy))
  (is= 0.26673862796330083 (first (alt-random/secure-rnd-lazy 4)))
  (is= 0.9214463212165593 (first (alt-random/secure-rnd-lazy 0))))

(deftest mersenne-rnd-lazy-test
  (is (spec-check alt-random/mersenne-rnd-lazy))
  (is= 0.8335762378570932 (first (alt-random/mersenne-rnd-lazy 4)))
  (is= 0.15071724896777527 (first (alt-random/mersenne-rnd-lazy 0))))

#_(ost/unstrument)