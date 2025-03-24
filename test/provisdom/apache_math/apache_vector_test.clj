(ns provisdom.apache-math.apache-vector-test
  (:require
    [clojure.spec.test.alpha :as st]
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.apache-math.apache-vector :as apache-v]))

;;;2 seconds

(set! *warn-on-reflection* true)

;;;TYPES
(deftest apache-vector?-test
  (with-instrument `apache-v/apache-vector?
    (is (spec-check apache-v/apache-vector?)))
  (with-instrument (st/instrumentable-syms)
    (is (apache-v/apache-vector? (apache-v/apache-vector [])))
    (is (apache-v/apache-vector? (apache-v/apache-vector [1])))
    (is-not (apache-v/apache-vector? "A"))
    (is-not (apache-v/apache-vector? [1 2]))))

;;;CONSTRUCTORS
(deftest apache-vector-&-apache-vector->vector-test
  (with-instrument `apache-v/apache-vector
    (is (spec-check apache-v/apache-vector)))
  (with-instrument `apache-v/apache-vector->vector
    (is (spec-check apache-v/apache-vector->vector)))
  (with-instrument (st/instrumentable-syms)
    (is= [] (apache-v/apache-vector->vector (apache-v/apache-vector [])))
    (is= [1.0] (apache-v/apache-vector->vector (apache-v/apache-vector [1.0])))
    (is= [1.0 2.0]
      (apache-v/apache-vector->vector (apache-v/apache-vector [1.0 2.0])))))
