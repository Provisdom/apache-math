(ns provisdom.apache-math.apache-vector
  "High-performance vector operations using Apache Commons Math.
  
  Provides vector operations built on Apache Commons Math ArrayRealVector
  for improved performance over standard Clojure vectors in numerical computations.
  Offers seamless conversion between Clojure vectors and Apache Commons vectors.
  
  Example usage:
    (def av (apache-vector [1 2 3]))
    (apache-vector->vector av) ;=> [1.0 2.0 3.0]"
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
  "Tests if an object is an Apache Commons vector.
  
  Parameters:
    x - object to test
  
  Returns true if x is an ArrayRealVector instance.
  
  Example:
    (apache-vector? (apache-vector [1 2 3])) ;=> true
    (apache-vector? [1 2 3]) ;=> false"
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
  "Creates an Apache Commons vector from a Clojure vector.
  
  Converts a standard Clojure vector of numbers into an Apache Commons
  ArrayRealVector for high-performance numerical operations.
  
  Parameters:
    v - Clojure vector of numbers
  
  Returns Apache Commons ArrayRealVector instance.
  
  Example:
    (apache-vector [1 2 3]) ;=> ArrayRealVector with elements [1.0 2.0 3.0]
    (apache-vector []) ;=> empty ArrayRealVector"
  [v]
  (ArrayRealVector. ^"[D" (double-array v)))

(s/fdef apache-vector
  :args (s/cat :v ::vector/vector)
  :ret ::apache-vector)

(defn apache-vector->vector
  "Converts an Apache Commons vector back to a Clojure vector.
  
  Extracts the underlying data from an Apache Commons vector and
  converts it to a standard Clojure vector of doubles.
  
  Parameters:
    apache-v - Apache Commons vector instance
  
  Returns Clojure vector of double values.
  
  Example:
    (def av (apache-vector [1 2 3]))
    (apache-vector->vector av) ;=> [1.0 2.0 3.0]"
  [apache-v]
  (vec (.toArray ^RealVector apache-v)))

(s/fdef apache-vector->vector
  :args (s/cat :apache-v ::apache-vector)
  :ret ::vector/vector)
