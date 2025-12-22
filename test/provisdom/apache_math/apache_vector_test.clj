(ns provisdom.apache-math.apache-vector-test
  (:require
    [clojure.test :as ct]
    [provisdom.test.core :as t]
    [provisdom.apache-math.apache-vector :as apache-v]))

;;;2 seconds

(set! *warn-on-reflection* true)

;;;TYPES
(ct/deftest apache-vector?-test
  (t/with-instrument `apache-v/apache-vector?
    (t/is-spec-check apache-v/apache-vector?))
  (t/with-instrument :all
    (t/is (apache-v/apache-vector? (apache-v/apache-vector [])))
    (t/is (apache-v/apache-vector? (apache-v/apache-vector [1])))
    (t/is-not (apache-v/apache-vector? "A"))
    (t/is-not (apache-v/apache-vector? [1 2]))))

;;;CONSTRUCTORS
(ct/deftest apache-vector-&-apache-vector->vector-test
  (t/with-instrument `apache-v/apache-vector
    (t/is-spec-check apache-v/apache-vector))
  (t/with-instrument `apache-v/apache-vector->vector
    (t/is-spec-check apache-v/apache-vector->vector))
  (t/with-instrument :all
    (t/is= [] (apache-v/apache-vector->vector (apache-v/apache-vector [])))
    (t/is= [1.0] (apache-v/apache-vector->vector (apache-v/apache-vector [1.0])))
    (t/is= [1.0 2.0]
      (apache-v/apache-vector->vector (apache-v/apache-vector [1.0 2.0])))))
