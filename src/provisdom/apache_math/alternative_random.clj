(ns provisdom.apache-math.alternative-random
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.math.core :as m])
  (:import
    [org.apache.commons.math3.random MersenneTwister
                                     ISAACRandom
                                     SobolSequenceGenerator]))

(def mdl 6)

(s/def ::seed ::m/long)
(s/def ::rnd ::m/prob)
(s/def ::rnd-lazy (s/every ::rnd))

(s/def ::rnd-vector
  (s/with-gen
    (s/coll-of ::rnd :kind clojure.core/vector? :into [])
    #(gen/vector (s/gen ::rnd) 0 mdl)))

(s/def ::apache-rng
  #(or (instance? MersenneTwister %)
       (instance? SobolSequenceGenerator %)
       (instance? ISAACRandom %)))

;;;APACHE RANDOM NUMBER GENERATORS
(defn quasi-rng
  "Creates an Apache RNG with better coverage but more predictable through a
  lazy sequence of vectors of size `dimensions`. Because of predictability, can
  be better for a single use simulation."
  [dimensions]
  (SobolSequenceGenerator. ^long dimensions))

(s/fdef quasi-rng
        :args (s/cat :dimensions (s/int-in 1 1000))
        :ret ::apache-rng)

(defn quasi-rnd-vector-lazy
  "Better coverage but more predictable through a lazy sequence of vectors of
  size `dimensions`. Because of predictability, can be better for a single use
  simulation."
  [dimensions]
  (let [qr (quasi-rng dimensions)]
    (repeatedly #(vec (.nextVector ^SobolSequenceGenerator qr)))))

(s/fdef quasi-rnd-vector-lazy
        :args (s/cat :dimensions (s/int-in 1 1000))
        :ret (s/every ::rnd-vector))

(defn secure-rng
  "Creates an Apache RNG that is less predictable but slower RNG than Mersenne
  Twister."
  [seed]
  (ISAACRandom. ^long seed))

(s/fdef secure-rng
        :args (s/cat :seed ::seed)
        :ret ::apache-rng)

(defn secure-rnd-lazy
  "A less predictable but slower rnd-lazy than Mersenne Twister."
  [seed]
  (repeatedly #(.nextDouble ^ISAACRandom (secure-rng seed))))

(s/fdef secure-rnd-lazy
        :args (s/cat :seed ::seed)
        :ret ::rnd-lazy)

(defn mersenne-rng
  "Creates an Apache RNG using `seed`."
  [seed]
  (MersenneTwister. ^long seed))

(s/fdef mersenne-rng
        :args (s/cat :seed ::seed)
        :ret ::apache-rng)

(defn mersenne-rnd-lazy
  "Returns rnd-lazy using `seed`."
  [seed]
  (repeatedly #(.nextDouble ^MersenneTwister (mersenne-rng seed))))

(s/fdef mersenne-rnd-lazy
        :args (s/cat :seed ::seed)
        :ret ::rnd-lazy)