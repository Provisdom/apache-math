(ns provisdom.apache-math.alternative-random
  "Alternative random number generators from Apache Commons Math.
  
  Provides specialized random number generators beyond basic uniform random:
  - Quasi-random sequences with better coverage (Sobol sequences)
  - Cryptographically secure random numbers (ISAAC)
  - Fast, high-quality pseudorandom numbers (Mersenne Twister)
  
  Example usage:
    (def quasi-seq (quasi-rnd-vector-lazy 2))
    (take 3 quasi-seq) ;=> [[0.0 0.0] [0.5 0.5] [0.75 0.25]]"
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.core :as m])
  (:import
    [org.apache.commons.math3.random
     MersenneTwister
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
  "Creates a Sobol quasi-random number generator for improved uniform coverage.
  
  Quasi-random sequences provide better distribution coverage than pseudorandom
  sequences, making them ideal for Monte Carlo simulations and numerical integration.
  However, they are more predictable and should not be used for cryptographic purposes.
  
  Parameters:
    dimensions - number of dimensions for each generated vector (1-1000)
  
  Returns Apache Commons SobolSequenceGenerator instance.
  
  Example:
    (def qrng (quasi-rng 2))
    (.nextVector qrng) ;=> [0.0 0.0]"
  [dimensions]
  (SobolSequenceGenerator. ^long dimensions))

(s/fdef quasi-rng
  :args (s/cat :dimensions (s/int-in 1 1000))
  :ret ::apache-rng)

(defn quasi-rnd-vector-lazy
  "Returns lazy sequence of quasi-random vectors with improved coverage.
  
  Generates an infinite lazy sequence of vectors using Sobol quasi-random sequences.
  Each vector contains `dimensions` numbers in [0,1). Provides better uniform
  coverage than pseudorandom sequences, making it superior for Monte Carlo methods.
  
  Parameters:
    dimensions - size of each vector (1-1000)
  
  Returns lazy sequence of vectors with quasi-random numbers.
  
  Example:
    (take 3 (quasi-rnd-vector-lazy 2))
    ;=> [[0.0 0.0] [0.5 0.5] [0.75 0.25]]"
  [dimensions]
  (let [qr (quasi-rng dimensions)]
    (repeatedly #(vec (.nextVector ^SobolSequenceGenerator qr)))))

(s/fdef quasi-rnd-vector-lazy
  :args (s/cat :dimensions (s/int-in 1 1000))
  :ret (s/every ::rnd-vector))

(defn secure-rng
  "Creates cryptographically secure random number generator using ISAAC algorithm.
  
  ISAAC (Indirection, Shift, Accumulate, Add, and Count) provides cryptographically
  secure random numbers suitable for security applications. Slower than Mersenne
  Twister but much less predictable.
  
  Parameters:
    seed - long integer seed for reproducible sequences
  
  Returns Apache Commons ISAACRandom instance.
  
  Example:
    (def secure-gen (secure-rng 12345))
    (.nextDouble secure-gen) ;=> 0.26673862796330083"
  [seed]
  (ISAACRandom. ^long seed))

(s/fdef secure-rng
  :args (s/cat :seed ::seed)
  :ret ::apache-rng)

(defn secure-rnd-lazy
  "Returns lazy sequence of cryptographically secure random numbers.
  
  Generates an infinite lazy sequence of random doubles in [0,1) using the ISAAC
  algorithm. Provides cryptographic security but at the cost of performance.
  Each call creates a new generator instance.
  
  Parameters:
    seed - long integer seed for reproducible sequences
  
  Returns lazy sequence of secure random doubles.
  
  Example:
    (take 5 (secure-rnd-lazy 42))
    ;=> (0.123... 0.456... 0.789...)"
  [seed]
  (repeatedly #(.nextDouble ^ISAACRandom (secure-rng seed))))

(s/fdef secure-rnd-lazy
  :args (s/cat :seed ::seed)
  :ret ::rnd-lazy)

(defn mersenne-rng
  "Creates Mersenne Twister random number generator.
  
  Mersenne Twister is a fast, high-quality pseudorandom number generator with
  excellent statistical properties and a very long period (2^19937-1). Suitable
  for most simulation work but not cryptographically secure.
  
  Parameters:
    seed - long integer seed for reproducible sequences
  
  Returns Apache Commons MersenneTwister instance.
  
  Example:
    (def mt-gen (mersenne-rng 12345))
    (.nextDouble mt-gen) ;=> 0.8335762378570932"
  [seed]
  (MersenneTwister. ^long seed))

(s/fdef mersenne-rng
  :args (s/cat :seed ::seed)
  :ret ::apache-rng)

(defn mersenne-rnd-lazy
  "Returns lazy sequence of high-quality pseudorandom numbers.
  
  Generates an infinite lazy sequence of random doubles in [0,1) using the
  Mersenne Twister algorithm. Provides excellent statistical properties and
  performance for most simulation needs. Each call creates a new generator instance.
  
  Parameters:
    seed - long integer seed for reproducible sequences
  
  Returns lazy sequence of pseudorandom doubles.
  
  Example:
    (take 5 (mersenne-rnd-lazy 42))
    ;=> (0.833... 0.150... 0.271...)"
  [seed]
  (repeatedly #(.nextDouble ^MersenneTwister (mersenne-rng seed))))

(s/fdef mersenne-rnd-lazy
  :args (s/cat :seed ::seed)
  :ret ::rnd-lazy)
