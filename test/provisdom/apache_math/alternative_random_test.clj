(ns provisdom.apache-math.alternative-random-test
  (:require
    [clojure.spec.test.alpha :as st]
    [clojure.test :refer :all]
    [provisdom.test.core :as t]
    [provisdom.apache-math.alternative-random :as alt-random]))

;;? seconds

(set! *warn-on-reflection* true)

;;;APACHE RANDOM NUMBER GENERATORS
(deftest quasi-rnd-vector-lazy-test
  (t/with-instrument `alt-random/quasi-rnd-vector-lazy
    (is (t/spec-check alt-random/quasi-rnd-vector-lazy)))
  (t/with-instrument (st/instrumentable-syms)
    (t/is= '([0.0] [0.5] [0.75])
      (take 3 (alt-random/quasi-rnd-vector-lazy 1)))
    (t/is= '([0.0 0.0] [0.5 0.5] [0.75 0.25])
      (take 3 (alt-random/quasi-rnd-vector-lazy 2)))))

(deftest secure-rnd-lazy-test
  (t/with-instrument `alt-random/secure-rnd-lazy
    (is (t/spec-check alt-random/secure-rnd-lazy)))
  (t/with-instrument (st/instrumentable-syms)
  (t/is= 0.26673862796330083 (first (alt-random/secure-rnd-lazy 4)))
  (t/is= 0.9214463212165593 (first (alt-random/secure-rnd-lazy 0)))))

(deftest mersenne-rnd-lazy-test
  (t/with-instrument `alt-random/mersenne-rnd-lazy
    (is (t/spec-check alt-random/mersenne-rnd-lazy)))
  (t/with-instrument (st/instrumentable-syms)
  (t/is= 0.8335762378570932 (first (alt-random/mersenne-rnd-lazy 4)))
  (t/is= 0.15071724896777527 (first (alt-random/mersenne-rnd-lazy 0)))))
