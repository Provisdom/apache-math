(ns provisdom.apache-math.apache-vector
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.vector :as vector])
  (:import
    [org.apache.commons.math3.linear
     ArrayRealVector
     RealVector]))

(declare apache-vector)

;;;APACHE VECTOR TYPES
(defn apache-vector?
  "Returns true if an Apache Commons vector."
  [x]
  (instance? ArrayRealVector x))

(s/fdef apache-vector?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::apache-vector
  (s/with-gen apache-vector?
    #(gen/fmap apache-vector (s/gen ::vector/vector))))

;;;APACHE VECTOR CONSTRUCTOR
(defn apache-vector
  "Returns a Apache Commons vector from a vector."
  [v]
  (ArrayRealVector. ^"[D" (double-array v)))

(s/fdef apache-vector
  :args (s/cat :v ::vector/vector)
  :ret ::apache-vector)

(defn apache-vector->vector
  "Converts an Apache Commons vector into a vector."
  [apache-v]
  (vec (.toArray ^RealVector apache-v)))

(s/fdef apache-vector->vector
  :args (s/cat :apache-v ::apache-vector)
  :ret ::vector/vector)
