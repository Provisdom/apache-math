(ns provisdom.apache-math.apache-matrix
  "High-performance matrix operations using Apache Commons Math.
  
  Provides comprehensive matrix operations including:
  - Matrix creation, conversion, and manipulation
  - Linear algebra: multiplication, addition, transpose, inverse
  - Matrix decompositions: LU, QR, Cholesky, SVD, Eigenvalue
  - Specialized matrix types: positive definite, correlation, symmetric
  - Matrix properties: rank, determinant, trace, condition number
  
  Built on Apache Commons Math Array2DRowRealMatrix for performance.
  Handles edge cases like empty matrices and provides comprehensive error handling.
  
  Example usage:
    (def m (->apache-matrix [[1 2] [3 4]]))
    (def inv (inverse m))
    (mx* m inv) ;=> approximately identity matrix"
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.arrays :as ar]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.random :as random]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.vector :as vector]
    [provisdom.utility-belt.anomalies :as anomalies])
  (:import
    [org.apache.commons.math3.linear
     Array2DRowRealMatrix RealMatrix
     QRDecomposition LUDecomposition CholeskyDecomposition
     RectangularCholeskyDecomposition Array2DRowRealMatrix
     EigenDecomposition SingularValueDecomposition
     RRQRDecomposition DecompositionSolver MatrixUtils]))

(declare
  apache-matrix? apache-square-matrix? ->apache-matrix eigen-decomposition
  cholesky-decomposition lu-decomposition-with-determinant-and-inverse
  rrqr-decomposition rows columns === transpose
  pos-definite-apache-matrix-finite-by-squaring
  pos-semidefinite-apache-matrix-finite-by-squaring diagonal mx* add
  covariance-apache-matrix->correlation-apache-matrix
  correlation-apache-matrix-by-squaring some-kv get-entry assoc-entry!
  symmetric-apache-matrix-by-averaging! assoc-diagonal!)

(s/def ::rank ::m/int-non-)

;;;MATRIX TYPES
(defn apache-matrix?
  "Returns true if an Apache Commons matrix."
  [x]
  (instance? Array2DRowRealMatrix x))

(s/fdef apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::apache-matrix
  (s/with-gen
    apache-matrix?
    #(gen/fmap ->apache-matrix (s/gen ::mx/matrix))))

(defn empty-apache-matrix?
  "Returns true if an empty Apache Commons matrix."
  [x]
  (and (apache-matrix? x) (zero? (rows x))))

(s/fdef empty-apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::empty-apache-matrix
  (s/with-gen
    empty-apache-matrix?
    #(= (->apache-matrix [[]]) %)))

(defn apache-matrix-finite?
  "Returns true if an Apache Commons matrix without infinite numbers."
  [x]
  (and (apache-matrix? x)
    (nil? (some-kv (fn [_ _ number]
                     (m/inf? number))
            x))))

(s/fdef apache-matrix-finite?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::apache-matrix-finite
  (s/with-gen
    apache-matrix-finite?
    #(gen/fmap ->apache-matrix (s/gen ::mx/matrix-finite))))

(defn square-apache-matrix?
  "Returns true if a square Apache Commons matrix (i.e., same number of rows and
  columns)."
  [x]
  (and (apache-matrix? x) (.isSquare ^Array2DRowRealMatrix x)))

(s/fdef square-apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::square-apache-matrix
  (s/with-gen
    square-apache-matrix?
    #(gen/fmap ->apache-matrix (s/gen ::mx/square-matrix))))

(s/def ::square-apache-matrix-finite
  (s/with-gen
    (s/and square-apache-matrix? apache-matrix-finite?)
    #(gen/fmap ->apache-matrix (s/gen ::mx/square-matrix-finite))))

(defn diagonal-apache-matrix?
  "Returns true if a diagonal matrix (the entries outside the main diagonal are
  all zero)."
  [x]
  (and (apache-matrix? x)
    (nil? (some-kv (fn [i j e]
                     (not (or (= i j) (zero? e))))
            x))))

(s/fdef diagonal-apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::diagonal-apache-matrix
  (s/with-gen
    diagonal-apache-matrix?
    #(gen/fmap ->apache-matrix (s/gen ::mx/diagonal-matrix))))

(defn upper-triangular-apache-matrix?
  "Returns true if an upper triangular matrix (the entries below the main
  diagonal are all zero)."
  [x]
  (and (square-apache-matrix? x)
    (nil? (some-kv (fn [i j e]
                     (not (or (<= i j) (zero? e))))
            x))))

(s/fdef upper-triangular-apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::upper-triangular-apache-matrix
  (s/with-gen
    upper-triangular-apache-matrix?
    #(gen/fmap ->apache-matrix (s/gen ::mx/upper-triangular-matrix))))

(defn lower-triangular-apache-matrix?
  "Returns true if a lower triangular matrix (the entries above the main
  diagonal are all zero)."
  [x]
  (and (square-apache-matrix? x)
    (nil? (some-kv (fn [i j e]
                     (not (or (>= i j) (zero? e))))
            x))))

(s/fdef lower-triangular-apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::lower-triangular-apache-matrix
  (s/with-gen
    lower-triangular-apache-matrix?
    #(gen/fmap ->apache-matrix (s/gen ::mx/lower-triangular-matrix))))

(defn symmetric-apache-matrix?
  "Returns true is a symmetric Apache Commons matrix."
  [x]
  (and (apache-matrix? x) (MatrixUtils/isSymmetric x 0.0)))

(s/fdef symmetric-apache-matrix?
  :args (s/cat :x any?)
  :ret boolean?)

(s/def ::symmetric-apache-matrix
  (s/with-gen
    symmetric-apache-matrix?
    #(gen/fmap ->apache-matrix (s/gen ::mx/symmetric-matrix))))

(defn pos-semidefinite-apache-matrix-finite?
  "Returns true if a finite `::pos-semidefinite-apache-matrix-finite`. Larger
  `accu` creates more false positives and less false negatives."
  [x accu]
  (and (symmetric-apache-matrix? x)
    (apache-matrix-finite? x)
    (every? #(m/roughly-non-? % accu)
      (try (::eigenvalues (eigen-decomposition x))
        (catch Exception _ nil)))))

(s/fdef pos-semidefinite-apache-matrix-finite?
  :args (s/cat :x any? :accu ::m/accu)
  :ret boolean?)

(s/def ::pos-semidefinite-apache-matrix-finite
  (s/with-gen
    #(pos-semidefinite-apache-matrix-finite? % m/sgl-close)
    #(gen/fmap (fn [m]
                 (pos-semidefinite-apache-matrix-finite-by-squaring
                   (->apache-matrix m)))
       (s/gen ::mx/square-matrix-finite))))

(defn pos-definite-apache-matrix-finite?
  "Tests if a matrix is positive definite with finite values.
  
  A positive definite matrix is symmetric with all positive eigenvalues.
  Such matrices are invertible and have Cholesky decompositions. Used in
  optimization, statistics, and defining valid covariance matrices.
  
  Parameters:
    x - matrix to test
    accu - accuracy threshold for eigenvalue positivity test
  
  Returns true if matrix is symmetric, finite, has positive eigenvalues,
  admits Cholesky decomposition, and is invertible.
  
  Example:
    (pos-definite-apache-matrix-finite?
      (->apache-matrix [[2 1] [1 2]]) 1e-10) ;=> true"
  [x accu]
  (and (symmetric-apache-matrix? x)
    (apache-matrix-finite? x)
    (every? #(> % accu)
      (try (::eigenvalues (eigen-decomposition x))
        (catch Exception _ nil)))
    (not (anomalies/anomaly? (cholesky-decomposition x)))
    (some? (::inverse (lu-decomposition-with-determinant-and-inverse x)))))

(s/fdef pos-definite-apache-matrix-finite?
  :args (s/cat :x any?
          :accu ::m/accu)
  :ret boolean?)

(s/def ::pos-definite-apache-matrix-finite
  (s/with-gen
    #(pos-definite-apache-matrix-finite? % m/sgl-close)
    #(gen/fmap (fn [m]
                 (pos-definite-apache-matrix-finite-by-squaring
                   (->apache-matrix m)))
       (s/gen ::mx/square-matrix-finite))))

(defn correlation-apache-matrix?
  "Tests if a matrix is a valid correlation matrix.
  
  A correlation matrix must be positive semidefinite (valid covariance
  structure) with unit diagonal elements (variables have unit variance).
  Elements represent correlations between variables and must be in [-1, 1].
  
  Parameters:
    x - matrix to test
    accu - accuracy threshold for positive semidefinite test
  
  Returns true if matrix is positive semidefinite with 1.0 on diagonal.
  
  Example:
    (correlation-apache-matrix?
      (->apache-matrix [[1.0 0.5] [0.5 1.0]]) 1e-10) ;=> true"
  [x accu]
  (and (pos-semidefinite-apache-matrix-finite? x accu)
    (every? m/one? (diagonal x))))

(s/fdef correlation-apache-matrix?
  :args (s/cat :x any? :accu ::m/accu)
  :ret boolean?)

(s/def ::correlation-apache-matrix
  (s/with-gen
    #(correlation-apache-matrix? % m/sgl-close)
    #(gen/fmap (fn [m]
                 (correlation-apache-matrix-by-squaring (->apache-matrix m)))
       (s/gen ::mx/square-matrix-finite))))

;;;MATRIX CONSTRUCTORS
(defn- block-apache-matrix->apache-matrix
  "Converts BlockRealMatrix (or other RealMatrix) to Array2DRowRealMatrix."
  [block-apache-matrix]
  (if (zero? (.getRowDimension ^RealMatrix block-apache-matrix))
    (Array2DRowRealMatrix.)
    (Array2DRowRealMatrix. ^"[[D" (.getData ^RealMatrix block-apache-matrix))))

(defn ->apache-matrix
  "Converts a Clojure matrix to Apache Commons matrix implementation.
  
  Creates an Array2DRowRealMatrix from a nested vector matrix structure.
  Handles empty matrices gracefully.
  
  Parameters:
    m - nested vector matrix (e.g., [[1 2] [3 4]])
  
  Returns Apache Commons Array2DRowRealMatrix instance or nil for invalid input.
  
  Example:
    (->apache-matrix [[1 2] [3 4]])
    (->apache-matrix [[]]) ;=> empty matrix"
  [m]
  (if (mx/empty-matrix? m)
    (Array2DRowRealMatrix.)
    (Array2DRowRealMatrix. ^"[[D" (ar/array2D :double m))))

(s/fdef ->apache-matrix
  :args (s/cat :m ::mx/matrix)
  :ret (s/nilable ::apache-matrix))

(defn apache-matrix->matrix
  "Converts Apache Commons matrix back to Clojure nested vector matrix.
  
  Extracts the underlying data from Apache Commons matrix and converts
  to standard Clojure nested vector format.
  
  Parameters:
    apache-m - Apache Commons matrix instance
  
  Returns nested vector matrix (e.g., [[1.0 2.0] [3.0 4.0]]).
  
  Example:
    (def am (->apache-matrix [[1 2] [3 4]]))
    (apache-matrix->matrix am) ;=> [[1.0 2.0] [3.0 4.0]]"
  [apache-m]
  (if (zero? (rows apache-m))
    [[]]
    (mapv vec (.getData ^Array2DRowRealMatrix apache-m))))

(s/fdef apache-matrix->matrix
  :args (s/cat :apache-m ::apache-matrix)
  :ret ::mx/matrix)

(defn pos-semidefinite-apache-matrix-finite-by-squaring
  "Returns a `::pos-semidefinite-apache-matrix-finite` by first
  squaring `square-apache-m-finite`."
  [square-apache-m-finite]
  (let [size (rows square-apache-m-finite)
        new-m square-apache-m-finite]
    (if (zero? size)
      square-apache-m-finite
      (loop [i 0]
        (when (< i 10)
          (let [lower-m (mx* new-m
                          (->apache-matrix
                            (mx/diagonal-matrix
                              (repeat size (m/pow 10 (m/one- (m/pow 2 i)))))))
                lower-m (mx* lower-m (transpose lower-m))
                _ (symmetric-apache-matrix-by-averaging! lower-m)]
            (if (pos-semidefinite-apache-matrix-finite?
                  lower-m m/sgl-close)
              lower-m
              (recur (inc i)))))))))

(s/fdef pos-semidefinite-apache-matrix-finite-by-squaring
  :args (s/cat :square-apache-m-finite ::square-apache-matrix-finite)
  :ret (s/nilable ::pos-semidefinite-apache-matrix-finite))

(defn pos-definite-apache-matrix-finite-by-squaring
  "Returns a `::pos-definite-apache-matrix-finite` squaring it. Will
  tweak original matrix if necessary to ensure positive definite."
  [square-apache-m-finite]
  (let [size (rows square-apache-m-finite)
        new-m square-apache-m-finite]
    (if (zero? size)
      square-apache-m-finite
      (loop [i 0]
        (when (< i 10)
          (let [lower-m (mx* new-m
                          (->apache-matrix
                            (mx/diagonal-matrix
                              (repeat size (m/pow 10 (m/one- (m/pow 2 i)))))))
                lower-m (if (zero? i)
                          lower-m
                          (add lower-m
                            (->apache-matrix
                              (mx/diagonal-matrix
                                size (fn [_]
                                       (* (if (odd? i) (- 1.0) 1.0)
                                         (m/pow 10 (- i m/sgl-digits))))))))
                lower-m (mx* lower-m (transpose lower-m))
                _ (symmetric-apache-matrix-by-averaging! lower-m)]
            (if (pos-definite-apache-matrix-finite? lower-m m/sgl-close)
              lower-m
              (recur (inc i)))))))))

(s/fdef pos-definite-apache-matrix-finite-by-squaring
  :args (s/cat :square-apache-m-finite ::square-apache-matrix-finite)
  :ret (s/nilable ::pos-definite-apache-matrix-finite))

(defn correlation-apache-matrix-by-squaring
  "Returns a finite Correlation Apache Commons matrix by first squaring
  `square-apache-m-finite`."
  [square-apache-m-finite]
  (if (zero? (rows square-apache-m-finite))
    square-apache-m-finite
    (let [new-m (covariance-apache-matrix->correlation-apache-matrix
                  (pos-definite-apache-matrix-finite-by-squaring
                    square-apache-m-finite))
          size (when new-m (rows new-m))]
      (when new-m
        (loop [i 0]
          (when (< i 100)
            (let [lower-m (mx* new-m (->apache-matrix
                                       (mx/diagonal-matrix
                                         (repeat size (m/one- (* 0.01 i))))))
                  _ (symmetric-apache-matrix-by-averaging! lower-m)
                  _ (assoc-diagonal! lower-m (repeat size 1.0))]
              (if (correlation-apache-matrix? lower-m m/sgl-close)
                lower-m
                (recur (inc i))))))))))

(s/fdef correlation-apache-matrix-by-squaring
  :args (s/cat :square-apache-m-finite ::square-apache-matrix-finite)
  :ret (s/nilable ::correlation-apache-matrix))

(defn rnd-pos-definite-apache-matrix-finite!
  "Returns a `::pos-definite-apache-matrix-finite` with a random
  spectrum. The orthogonal matrices are generated by using 2 × `size` composed
  Householder reflections.
  Alternative #1: Sample from the Inverse-Wishart Distribution.
  Alternative #2: Use [[pos-definite-apache-matrix-finite-by-squaring]] with a
  random square matrix."
  [size]
  (if (zero? size)
    (->apache-matrix [[]])
    (loop [i 0]
      (when (< i 100)
        (let [m (->apache-matrix
                  (mx/rnd-spectral-matrix!
                    (vec (take size (random/rnd-lazy!)))))]
          (if (pos-definite-apache-matrix-finite? m m/sgl-close)
            m
            (recur (inc i))))))))

(s/fdef rnd-pos-definite-apache-matrix-finite!
  :args (s/cat :size ::vector/size)
  :ret ::pos-definite-apache-matrix-finite)

(defn rnd-correlation-apache-matrix!
  "Returns a correlation Apache Commons matrix from a covariance matrix
  with a random spectrum. The orthogonal matrices are generated by using
  2 * `size` composed Householder reflections.
  Alternative #1: Sample Covariance from the Inverse-Wishart Distribution.
  Alternative #2: Use [[correlation-apache-matrix-by-squaring]] with a
  random square matrix."
  [size]
  (if (zero? size)
    (->apache-matrix [[]])
    (covariance-apache-matrix->correlation-apache-matrix
      (rnd-pos-definite-apache-matrix-finite! size))))

(s/fdef rnd-correlation-apache-matrix!
  :args (s/cat :size ::vector/size)
  :ret ::correlation-apache-matrix)

;;;MATRIX INFO
(defn rows
  "Gets the number of rows in a matrix.
  
  Parameters:
    apache-m - Apache Commons matrix
  
  Returns number of rows as integer.
  
  Example:
    (rows (->apache-matrix [[1 2] [3 4] [5 6]])) ;=> 3"
  [apache-m]
  (.getRowDimension ^Array2DRowRealMatrix apache-m))

(s/fdef rows
  :args (s/cat :apache-m ::apache-matrix)
  :ret ::mx/rows)

(defn columns
  "Gets the number of columns in a matrix.
  
  Parameters:
    apache-m - Apache Commons matrix
  
  Returns number of columns as integer.
  
  Example:
    (columns (->apache-matrix [[1 2 3] [4 5 6]])) ;=> 3"
  [apache-m]
  (.getColumnDimension ^Array2DRowRealMatrix apache-m))

(s/fdef columns
  :args (s/cat :apache-m ::apache-matrix)
  :ret ::mx/columns)

(defn get-entry
  "Retrieves a single element from a matrix.
  
  Gets the element at the specified row and column position.
  Uses zero-based indexing.
  
  Parameters:
    apache-m - Apache Commons matrix
    row - row index (0-based)
    column - column index (0-based)
  
  Returns matrix element value or NaN if indices are out of bounds.
  
  Example:
    (def m (->apache-matrix [[1 2] [3 4]]))
    (get-entry m 0 1) ;=> 2.0
    (get-entry m 1 0) ;=> 3.0"
  [apache-m row column]
  (try (.getEntry ^Array2DRowRealMatrix apache-m
         row
         column)
    (catch Exception _ m/nan)))

(s/fdef get-entry
  :args (s/and (s/cat :apache-m ::apache-matrix
                 :row ::mx/row
                 :column ::mx/column)
          (fn [{:keys [apache-m row column]}]
            (and (< row (rows apache-m))
              (< column (columns apache-m)))))
  :ret ::m/number)

(defn get-row
  "Gets a `row` of an Apache Commons matrix, as a vector."
  [apache-m row]
  (try (vec (.getRow ^Array2DRowRealMatrix apache-m
              row))
    (catch Exception _ nil)))

(s/fdef get-row
  :args (s/and (s/cat :apache-m ::apache-matrix :row ::mx/row)
          (fn [{:keys [apache-m row]}]
            (< row (rows apache-m))))
  :ret (s/nilable ::vector/vector))

(defn get-column
  "Gets a `column` of an Apache Commons matrix, as a vector."
  [apache-m column]
  (try (vec (.getColumn ^Array2DRowRealMatrix apache-m
              column))
    (catch Exception _ nil)))

(s/fdef get-column
  :args (s/and (s/cat :apache-m ::apache-matrix :column ::mx/column)
          (fn [{:keys [apache-m column]}]
            (< column (columns apache-m))))
  :ret (s/nilable ::vector/vector))

(defn diagonal
  "Extracts a diagonal from a matrix as a vector.
  
  Returns elements from a diagonal line through the matrix. The main diagonal
  (k=0) goes from top-left to bottom-right. Positive k values select diagonals
  above the main diagonal, negative k values select diagonals below.
  
  Parameters:
    apache-m - Apache Commons matrix
    k - (optional) diagonal offset: 0=main, positive=above, negative=below
  
  Returns vector of diagonal elements.
  
  Example:
    (def m (->apache-matrix [[1 2 3] [4 5 6] [7 8 9]]))
    (diagonal m)   ;=> [1 5 9] (main diagonal)
    (diagonal m 1) ;=> [2 6] (first super-diagonal)
    (diagonal m -1) ;=> [4 8] (first sub-diagonal)"
  ([apache-m]
   (if (zero? (rows apache-m))
     []
     (reduce (fn [tot e]
               (conj tot (get-entry apache-m e e)))
       []
       (range (min (rows apache-m) (columns apache-m))))))
  ([apache-m k]
   (if (zero? (rows apache-m))
     []
     (let [r (if (neg? k) (- k) 0)
           c (if (pos? k) k 0)
           nc (- (columns apache-m) c)
           nr (- (rows apache-m) r)
           start (- (min r c))
           end (min nc nr)]
       (if (pos? end)
         (vec (for [i (range start end)]
                (get-entry apache-m (+ i r) (+ i c))))
         [])))))

(s/fdef diagonal
  :args (s/cat :apache-m ::apache-matrix :k (s/? ::m/int))
  :ret ::vector/vector)

(defn trace
  "Calculates the trace of a square matrix.
  
  The trace is the sum of elements on the main diagonal (top-left to bottom-right).
  Only defined for square matrices.
  
  Parameters:
    square-apache-m - square Apache Commons matrix
  
  Returns sum of diagonal elements.
  
  Example:
    (def m (->apache-matrix [[1 2] [3 4]]))
    (trace m) ;=> 5.0 (1 + 4)"
  [square-apache-m]
  (.getTrace ^Array2DRowRealMatrix square-apache-m))

(s/fdef trace
  :args (s/cat :square-apache-m ::square-apache-matrix)
  :ret ::m/number)

(defn get-slices-as-matrix
  "Performs a slice on the Apache Commons matrix given by the options.
  Options:
    `::mx/row-indices` returns all rows by default, can pass a row index or
      sequence of row indices
    `::mx/column-indices` returns all columns by default, can pass a column
      index or sequence of column indices
    `::mx/exception-row-indices` can pass a row index or sequence of row indices
      to exclude
    `::mx/exception-column-indices` can pass a column index or sequence of
      column indices to exclude.
  Exceptions override inclusions. Can be used to permute matrix through index
   sequence ordering."
  [apache-m {:keys [::mx/row-indices
                    ::mx/column-indices
                    ::mx/exception-row-indices
                    ::mx/exception-column-indices]}]
  (let [calc-fn (fn [i except-i n]
                  (cond (and (not i) (not except-i)) true
                    (not except-i) (if (number? i)
                                     (if (< i n) i [])
                                     (remove #(>= % n) i))
                    (number? i) (if (number? except-i)
                                  (if (= except-i i)
                                    []
                                    (if (< i n) i []))
                                  (if (contains? (set except-i) i)
                                    []
                                    (if (< i n) i [])))
                    :else (let [indices (or i (range n))]
                            (if (number? except-i)
                              (remove
                                (fn [index]
                                  (or (= except-i index) (>= index n)))
                                indices)
                              (reduce
                                (fn [tot e]
                                  (if (or (>= e n) (some #(= % e) except-i))
                                    tot
                                    (conj tot e)))
                                []
                                indices)))))
        n-rows (rows apache-m)
        rs (calc-fn row-indices exception-row-indices n-rows)
        cs (calc-fn column-indices exception-column-indices (columns apache-m))
        new-m (cond
                (or (and (coll? rs) (empty? rs)) (and (coll? cs) (empty? cs)))
                nil

                (and (number? rs) (number? cs))
                [[(get-entry apache-m rs cs)]]

                (and (number? rs) (coll? cs))
                (mx/row-matrix (let [row-vector (get-row apache-m rs)]
                                 (map (fn [c]
                                        (get row-vector c))
                                   cs)))

                (and (number? rs) (true? cs))
                (mx/row-matrix (get-row apache-m rs))

                (and (coll? rs) (number? cs))
                (mx/column-matrix (let [column-vector (get-column apache-m cs)]
                                    (map (fn [r]
                                           (get column-vector r))
                                      rs)))

                (and (coll? rs) (coll? cs))
                (mapv (fn [row-vector]
                        (reduce (fn [tot column]
                                  (conj tot (get row-vector column)))
                          []
                          cs))
                  (map (fn [r]
                         (get-row apache-m r))
                    rs))

                (and (coll? rs) (true? cs))
                (mapv (fn [r]
                        (get-row apache-m r))
                  rs)

                (and (true? rs) (number? cs))
                (mx/column-matrix (get-column apache-m cs))

                (and (true? rs) (coll? cs))
                (mapv (fn [row]
                        (reduce (fn [tot column]
                                  (conj tot (get-entry apache-m row column)))
                          []
                          cs))
                  (range n-rows))

                (and (true? rs) (true? cs))
                apache-m)]
    (if new-m
      (if (apache-matrix? new-m)
        new-m
        (->apache-matrix new-m))
      (->apache-matrix [[]]))))

(s/fdef get-slices-as-matrix
  :args (s/cat :apache-m ::apache-matrix
          :opts (s/keys :opt [::mx/row-indices
                              ::mx/column-indices
                              ::mx/exception-row-indices
                              ::mx/exception-column-indices]))
  :ret ::apache-matrix)

(s/def ::top-left ::apache-matrix)
(s/def ::bottom-left ::apache-matrix)
(s/def ::top-right ::apache-matrix)
(s/def ::bottom-right ::apache-matrix)

(defn matrix-partition
  "Returns a map containing the four sub-matrices labeled `::top-left`,
  `::bottom-left`, `::top-right`, and `::bottom-right`.
  `first-bottom-row` is the bottom of where the slice will occur.
  `first-right-column` is the right edge of where the slice will occur."
  [apache-m first-bottom-row first-right-column]
  {::top-left     (get-slices-as-matrix
                    apache-m {::mx/row-indices    (range first-bottom-row)
                              ::mx/column-indices (range first-right-column)})
   ::bottom-left  (get-slices-as-matrix
                    apache-m
                    {::mx/exception-row-indices (range first-bottom-row)
                     ::mx/column-indices        (range first-right-column)})
   ::top-right    (get-slices-as-matrix
                    apache-m
                    {::mx/row-indices              (range first-bottom-row)
                     ::mx/exception-column-indices (range first-right-column)})
   ::bottom-right (get-slices-as-matrix
                    apache-m
                    {::mx/exception-row-indices (range first-bottom-row)
                     ::mx/exception-column-indices
                     (range first-right-column)})})

(s/fdef matrix-partition
  :args (s/cat :apache-m ::apache-matrix
          :first-bottom-row ::mx/row
          :first-right-column ::mx/column)
  :ret (s/keys :req [::top-left ::bottom-left ::top-right ::bottom-right]))

(defn some-kv
  "Returns the first logical true value of (pred row column number) for any
  number in Apache Commons matrix, else nil.
  Options: `::mx/by-row?` (default: true)."
  ([pred apache-m] (some-kv pred apache-m {::mx/by-row? true}))
  ([pred apache-m {:keys [::mx/by-row?] :or {by-row? true}}]
   (let [mt (if by-row?
              apache-m
              (transpose apache-m))
         rows (rows mt)]
     (loop [row 0]
       (when (< row rows)
         (or (vector/some-kv (fn [column number]
                               (pred row column number))
               (get-row mt row))
           (recur (inc row))))))))

(s/fdef some-kv
  :args (s/cat :pred (s/fspec :args (s/cat :row ::mx/row
                                      :column ::mx/column
                                      :number ::m/number)
                       :ret boolean?)
          :apache-m ::apache-matrix
          :opts (s/? (s/keys :opt [::mx/by-row?])))
  :ret (s/nilable ::m/number))

;;;MATRIX MANIPULATION
(defn transpose
  "Transposes a matrix by swapping rows and columns.
  
  Creates a new matrix where element (i,j) becomes element (j,i).
  
  Parameters:
    apache-m - Apache Commons matrix to transpose
  
  Returns new transposed matrix.
  
  Example:
    (def m (->apache-matrix [[1 2 3] [4 5 6]]))
    (transpose m) ;=> [[1 4] [2 5] [3 6]]"
  [apache-m]
  (if (zero? (rows apache-m))
    apache-m
    (.transpose ^Array2DRowRealMatrix apache-m)))

(s/fdef transpose
  :args (s/cat :apache-m ::apache-matrix)
  :ret ::apache-matrix)

(defn assoc-entry!
  "Sets an entry in an Apache Commons matrix."
  [apache-m row column number]
  (when (and (< row (rows apache-m)) (< column (columns apache-m)))
    (.setEntry ^Array2DRowRealMatrix apache-m
      ^long row
      ^long column
      ^double (double number))))

(s/fdef assoc-entry!
  :args (s/cat :apache-m ::apache-matrix
          :row ::mx/row
          :column ::mx/column
          :number ::m/number)
  :ret nil)

(defn assoc-diagonal!
  "Sets a diagonal in an Apache Commons matrix using the specified numbers."
  [apache-m numbers]
  (let [v (vec numbers)
        c (count numbers)]
    (when (= c (count (diagonal apache-m)))
      (doseq [rc (range c)]
        (assoc-entry! apache-m rc rc (get v rc 0.0))))))

(s/fdef assoc-diagonal!
  :args (s/cat :apache-m ::apache-matrix :numbers ::m/numbers)
  :ret nil)

(defn symmetric-apache-matrix-by-averaging!
  "Updates symmetric Apache Commons matrix where each element above or below the
  diagonal is equal to the average of the corresponding numbers. This is useful
  to help with rounding errors."
  [square-apache-m]
  (doseq [row (range (rows square-apache-m))]
    (doseq [column (range (inc row) (columns square-apache-m))]
      (let [number (* 0.5 (+ (get-entry square-apache-m row column)
                            (get-entry square-apache-m column row)))]
        (assoc-entry! square-apache-m row column number)
        (assoc-entry! square-apache-m column row number)))))

(s/fdef symmetric-apache-matrix-by-averaging!
  :args (s/cat :square-apache-m ::square-apache-matrix)
  :ret nil)

(defn concat-rows
  "Appends rows from all the Apache Common matrices after the first to the
  first. Each matrix's column count must be the same or will return nil."
  ([] (->apache-matrix [[]]))
  ([apache-m] apache-m)
  ([apache-m & apache-ms]
   (when-let [new-m (apply mx/concat-rows
                      (apache-matrix->matrix apache-m)
                      (map apache-matrix->matrix apache-ms))]
     (->apache-matrix new-m))))

(s/fdef concat-rows
  :args (s/or :zero (s/cat)
          :one+ (s/cat :apache-m ::apache-matrix
                  :apache-ms (s/* ::apache-matrix)))
  :ret (s/nilable ::apache-matrix))

(defn concat-columns
  "Appends columns from all the Apache Common matrices after the first to the
  first. Each matrix's row count must be the same or will return nil."
  ([] (->apache-matrix [[]]))
  ([apache-m] apache-m)
  ([apache-m & apache-ms]
   (when-let [new-m (apply mx/concat-columns
                      (apache-matrix->matrix apache-m)
                      (map apache-matrix->matrix apache-ms))]
     (->apache-matrix new-m))))

(s/fdef concat-columns
  :args (s/or :zero (s/cat)
          :one+ (s/cat :apache-m ::apache-matrix
                  :apache-ms (s/* ::apache-matrix)))
  :ret (s/nilable ::apache-matrix))

(defn correlation-apache-matrix->covariance-apache-matrix
  "Returns Covariance Apache Commons matrix from a Correlation Apache Commons
  matrix."
  [correlation-apache-matrix variances]
  (when (= (count variances) (rows correlation-apache-matrix))
    (if (zero? (count variances))
      correlation-apache-matrix
      (let [sqrt-m (->apache-matrix
                     (mx/diagonal-matrix
                       (mapv m/sqrt variances)))
            cov (mx* sqrt-m correlation-apache-matrix sqrt-m)
            _ (symmetric-apache-matrix-by-averaging! cov)]
        (when (pos-definite-apache-matrix-finite? cov m/sgl-close)
          cov)))))

(s/fdef correlation-apache-matrix->covariance-apache-matrix
  :args (s/cat :correlation-apache-matrix ::correlation-apache-matrix
          :variances ::vector/vector-finite+)
  :ret (s/nilable ::pos-definite-apache-matrix-finite))

(defn covariance-apache-matrix->correlation-apache-matrix
  "Returns Correlation Apache Commons matrix from a Covariance Apache Commons
  matrix."
  [covariance-apache-matrix]
  (if (zero? (rows covariance-apache-matrix))
    covariance-apache-matrix
    (let [inv-sqrt (->apache-matrix
                     (mx/diagonal-matrix
                       (map #(m/pow % -0.5)
                         (diagonal covariance-apache-matrix))))
          corr (mx* inv-sqrt covariance-apache-matrix inv-sqrt)
          _ (symmetric-apache-matrix-by-averaging! corr)
          _ (assoc-diagonal! corr (repeat (rows inv-sqrt) 1.0))]
      (when (correlation-apache-matrix? corr m/sgl-close)
        corr))))

(s/fdef covariance-apache-matrix->correlation-apache-matrix
  :args (s/cat :covariance-apache-matrix ::pos-definite-apache-matrix-finite)
  :ret (s/nilable ::correlation-apache-matrix))

;;;MATRIX MATH
(defn ===
  "Apache Commons matrix equality that works with NaN."
  ([apache-m] true)
  ([apache-m1 apache-m2]
   (tensor/=== (apache-matrix->matrix apache-m1)
     (apache-matrix->matrix apache-m2)))
  ([apache-m1 apache-m2 & apache-ms]
   (apply tensor/=== (apache-matrix->matrix apache-m1)
     (apache-matrix->matrix apache-m2)
     (map apache-matrix->matrix apache-ms))))

(s/fdef ===
  :args (s/or :one (s/cat :apache-m ::apache-matrix)
          :two+ (s/cat :apache-m1 ::apache-matrix
                  :apache-m2 ::apache-matrix
                  :apache-ms (s/* ::apache-matrix)))
  :ret boolean?)

(defn mx*
  "Performs matrix multiplication on Apache Commons matrices.
  
  Multiplies two or more matrices using standard linear algebra rules.
  The number of columns in the first matrix must equal the number of
  rows in the second matrix for each multiplication.
  
  Parameters:
    apache-m1, apache-m2, ... - Apache Commons matrices to multiply
  
  Returns resulting matrix or nil if dimensions are incompatible.
  
  Example:
    (def a (->apache-matrix [[1 2] [3 4]]))
    (def b (->apache-matrix [[5 6] [7 8]]))
    (mx* a b) ;=> [[19 22] [43 50]]"
  ([apache-m] apache-m)
  ([apache-m1 apache-m2]
   (when (= (columns apache-m1) (rows apache-m2))
     (if (zero? (rows apache-m1))
       apache-m1
       (.multiply ^Array2DRowRealMatrix apache-m1
         ^Array2DRowRealMatrix apache-m2))))
  ([apache-m1 apache-m2 & apache-ms]
   (when-let [apache-m3 (mx* apache-m1 apache-m2)]
     (apply mx* apache-m3 apache-ms))))

(s/fdef mx*
  :args (s/or :one (s/cat :apache-m ::apache-matrix)
          :two+ (s/cat :apache-m1 ::apache-matrix
                  :apache-m2 ::apache-matrix
                  :apache-ms (s/* ::apache-matrix)))
  :ret (s/nilable ::apache-matrix))

(defn add
  "Performs element-wise addition of matrices.
  
  Adds corresponding elements of two or more matrices. All matrices
  must have the same dimensions.
  
  Parameters:
    apache-m1, apache-m2, ... - matrices to add (same dimensions required)
  
  Returns resulting sum matrix or nil if dimensions are incompatible.
  
  Example:
    (def a (->apache-matrix [[1 2] [3 4]]))
    (def b (->apache-matrix [[5 6] [7 8]]))
    (add a b) ;=> [[6 8] [10 12]]"
  ([apache-m] apache-m)
  ([apache-m1 apache-m2]
   (if (and (zero? (rows apache-m1)) (zero? (rows apache-m2)))
     apache-m1
     (try (.add ^Array2DRowRealMatrix apache-m1
            ^Array2DRowRealMatrix apache-m2)
       (catch Exception _ nil))))
  ([apache-m1 apache-m2 & apache-ms]
   (when-let [apache-m3 (add apache-m1 apache-m2)]
     (apply add apache-m3 apache-ms))))

(s/fdef add
  :args (s/or :one (s/cat :apache-m ::apache-matrix)
          :two+ (s/cat :apache-m1 ::apache-matrix
                  :apache-m2 ::apache-matrix
                  :apache-ms (s/* ::apache-matrix)))
  :ret (s/nilable ::apache-matrix))

(defn subtract
  "Apache Commons matrix subtraction."
  ([apache-m] apache-m)
  ([apache-m1 apache-m2]
   (if (and (zero? (rows apache-m1)) (zero? (rows apache-m2)))
     apache-m1
     (try (.subtract ^Array2DRowRealMatrix apache-m1
            ^Array2DRowRealMatrix apache-m2)
       (catch Exception _ nil))))
  ([apache-m1 apache-m2 & apache-ms]
   (when-let [apache-m3 (subtract apache-m1 apache-m2)]
     (apply subtract apache-m3 apache-ms))))

(s/fdef subtract
  :args (s/or :one (s/cat :apache-m ::apache-matrix)
          :two+ (s/cat :apache-m1 ::apache-matrix
                  :apache-m2 ::apache-matrix
                  :apache-ms (s/* ::apache-matrix)))
  :ret (s/nilable ::apache-matrix))

(defn scalar-add
  "Apache Commons matrix addition to a scalar."
  [apache-m number]
  (if (zero? (rows apache-m))
    apache-m
    (.scalarAdd ^Array2DRowRealMatrix apache-m
      (double number))))

(s/fdef scalar-add
  :args (s/cat :apache-m ::apache-matrix
          :number ::m/number)
  :ret ::apache-matrix)

(defn scalar-multiply
  "Apache Commons matrix multiplication to a scalar."
  [apache-m number]
  (if (zero? (rows apache-m))
    apache-m
    (.scalarMultiply ^Array2DRowRealMatrix apache-m
      (double number))))

(s/fdef scalar-multiply
  :args (s/cat :apache-m ::apache-matrix
          :number ::m/number)
  :ret ::apache-matrix)

;;;MATRIX DECOMPOSITION
(s/def ::inverse (s/nilable ::square-apache-matrix))

(defn inverse
  "Computes the inverse of a square matrix.
  
  Uses Apache Commons Math's adaptive inverse computation which selects
  the most appropriate decomposition method (QR, LU, etc.) based on
  matrix properties for optimal performance and numerical stability.
  
  Parameters:
    square-apache-m - square Apache Commons matrix
  
  Returns inverse matrix or nil if matrix is singular (non-invertible).
  
  Example:
    (def m (->apache-matrix [[2 0] [0 3]]))
    (inverse m) ;=> [[0.5 0] [0 0.333...]]"
  [square-apache-m]
  (if (zero? (rows square-apache-m))
    square-apache-m
    (try (block-apache-matrix->apache-matrix
           (MatrixUtils/inverse ^RealMatrix square-apache-m))
      (catch Exception _ nil))))

(s/fdef inverse
  :args (s/cat :square-apache-m ::square-apache-matrix)
  :ret ::inverse)

(s/def ::LU-permutation (s/nilable ::square-apache-matrix))
(s/def ::L (s/nilable ::lower-triangular-apache-matrix))
(s/def ::U (s/nilable ::upper-triangular-apache-matrix))
(s/def ::determinant ::m/number)

(defn lu-decomposition-with-determinant-and-inverse
  "Returns a map containing
      `::L` -- the lower triangular factor Apache Commons matrix
      `::U` -- the upper triangular factor Apache Commons matrix
      `::LU-permutation` -- the permutation Apache Commons matrix
      `::determinant` -- the determinant
      `::inverse` -- the inverse Apache Commons matrix."
  [square-apache-m]
  (if (zero? (rows square-apache-m))
    {::L              square-apache-m
     ::U              square-apache-m
     ::LU-permutation square-apache-m
     ::determinant    m/nan
     ::inverse        square-apache-m}
    (let [lud (LUDecomposition. ^Array2DRowRealMatrix square-apache-m)
          s (.getSolver lud)
          inverse (try (.getInverse s) (catch Exception _ nil))
          det (.getDeterminant lud)]
      {::L              (.getL lud)
       ::U              (.getU lud)
       ::LU-permutation (.getP lud)
       ::determinant    det
       ::inverse        inverse})))

(s/fdef lu-decomposition-with-determinant-and-inverse
  :args (s/cat :square-apache-m ::square-apache-matrix)
  :ret (s/keys :req [::L ::U ::LU-permutation ::determinant ::inverse]))

(defn lu-decomposition-with-determinant
  "Returns a map containing:
      `::L` -- the lower triangular factor Apache Commons matrix
      `::U` -- the upper triangular factor Apache Commons matrix
      `::LU-permutation` -- the permutation Apache Commons matrix
      `::determinant` -- the determinant."
  [square-apache-m]
  (if (zero? (rows square-apache-m))
    {::L              square-apache-m
     ::U              square-apache-m
     ::LU-permutation square-apache-m
     ::determinant    m/nan}
    (let [lud (LUDecomposition. ^Array2DRowRealMatrix square-apache-m)]
      {::L              (.getL lud)
       ::U              (.getU lud)
       ::LU-permutation (.getP lud)
       ::determinant    (.getDeterminant lud)})))

(s/fdef lu-decomposition-with-determinant
  :args (s/cat :square-apache-m ::square-apache-matrix)
  :ret (s/keys :req [::L ::U ::LU-permutation ::determinant]))

(s/def ::eigenvectorsT ::square-apache-matrix)
(s/def ::eigenvalues-matrix ::apache-matrix)
(s/def ::eigenvectors ::square-apache-matrix)
(s/def ::eigenvalues (s/nilable ::vector/vector))

(defn eigen-decomposition
  "Computes eigenvalue decomposition of a square matrix.
  
  Finds eigenvalues and eigenvectors such that A * v = λ * v, where v is an
  eigenvector and λ is the corresponding eigenvalue. The matrix is decomposed
  as A = V * D * V^(-1), where V contains eigenvectors and D contains eigenvalues.
  
  Parameters:
    square-apache-m-finite - square matrix with finite values
  
  Returns map with keys:
    ::eigenvectorsT - matrix with eigenvectors as columns
    ::eigenvalues-matrix - diagonal matrix of eigenvalues
    ::eigenvalues - vector of real eigenvalues (nil if complex eigenvalues exist)
    ::eigenvectors - matrix with eigenvectors as rows
  
  Returns anomaly map if decomposition fails.
  
  Example:
    (eigen-decomposition (->apache-matrix [[3 1] [0 2]]))
    ;=> eigenvalues [3.0 2.0], eigenvectors for 3×3 identity-like behavior"
  [square-apache-m-finite]
  (if (zero? (rows square-apache-m-finite))
    {::eigenvectorsT      square-apache-m-finite
     ::eigenvalues-matrix square-apache-m-finite
     ::eigenvalues        []
     ::eigenvectors       square-apache-m-finite}
    (try (let [r (EigenDecomposition. ^Array2DRowRealMatrix
                   square-apache-m-finite)
               eigenvalues (when-not (.hasComplexEigenvalues r)
                             (vec (.getRealEigenvalues r)))]
           {::eigenvectorsT      (.getV r)
            ::eigenvalues-matrix (.getD r)
            ::eigenvalues        eigenvalues
            ::eigenvectors       (.getVT r)})
      (catch Exception e
        {::anomalies/message  (.getMessage e)
         ::anomalies/fn       (var eigen-decomposition)
         ::anomalies/category ::anomalies/third-party}))))

(s/fdef eigen-decomposition
  :args (s/cat :square-apache-m-finite ::square-apache-matrix-finite)
  :ret (s/or :anomaly ::anomalies/anomaly
         :eigen (s/keys :req [::eigenvectorsT ::eigenvalues-matrix
                              ::eigenvalues ::eigenvectors])))

(s/def ::cholesky-L ::lower-triangular-apache-matrix)
(s/def ::cholesky-LT ::upper-triangular-apache-matrix)
(defn cholesky-decomposition
  "Computes Cholesky decomposition of a positive definite matrix.
  
  Decomposes a positive definite matrix A into A = L * L^T, where L is a
  lower triangular matrix. This is the matrix square root operation and is
  useful for solving linear systems and generating correlated random numbers.
  
  Parameters:
    pos-definite-apache-m - positive definite matrix (symmetric with all positive eigenvalues)
  
  Returns map with keys:
    ::cholesky-L - lower triangular matrix L
    ::cholesky-LT - upper triangular matrix L^T (transpose of L)
  
  Returns anomaly map if matrix is not positive definite.
  
  Example:
    (cholesky-decomposition (->apache-matrix [[4 2] [2 2]]))
    ;=> L = [[2 0] [1 1]], L^T = [[2 1] [0 1]]"
  [pos-definite-apache-m]
  (if (zero? (rows pos-definite-apache-m))
    {::cholesky-L  pos-definite-apache-m
     ::cholesky-LT pos-definite-apache-m}
    (try (let [r (CholeskyDecomposition.
                   ^Array2DRowRealMatrix pos-definite-apache-m)]
           {::cholesky-L  (.getL r)
            ::cholesky-LT (.getLT r)})
      (catch Exception e
        {::anomalies/message  (.getMessage e)
         ::anomalies/fn       (var cholesky-decomposition)
         ::anomalies/category ::anomalies/third-party}))))

(s/fdef cholesky-decomposition
  :args (s/cat :pos-semidefinite-apache-m ::pos-semidefinite-apache-matrix-finite)
  :ret (s/or :anomaly ::anomalies/anomaly
         :res (s/keys :req [::cholesky-L ::cholesky-LT])))

(s/def ::rectangular-root ::apache-matrix)

(defn rectangular-cholesky-decomposition
  "Calculates the rectangular Cholesky decomposition of a
  `::pos-semidefinite-apache-matrix-finite`. The rectangular Cholesky
  decomposition of a real ::pos-semidefinite-apache-matrix-finite consists of a
  `rectangular-root` matrix with the same number of rows such that
  ::pos-semidefinite-apache-matrix-finite is almost equal to
  `rectangular-root` × (transpose `rectangular-root`), depending on a
  user-defined tolerance. In a sense, this is the square root of
  ::pos-semidefinite-apache-matrix-finite. The difference with respect to the
  regular CholeskyDecomposition is that rows/columns may be permuted (hence the
  rectangular shape instead of the traditional triangular shape) and there is a
  threshold to ignore small diagonal elements. This is used for example to
  generate correlated random n-dimensions vectors in a p-dimension subspace
  (p < n). In other words, it allows generating random vectors from a covariance
  matrix that is only positive semidefinite, and not positive definite.

  `accu` - Diagonal elements threshold under which columns are considered to be
    dependent on previous ones and are discarded.

  Returns a map containing:
    `::rectangular-root` -- rectangular root Apache Commons matrix
    `::rank` -- rank is the number of independent rows of original matrix, and
      the number of columns of `rectangular-root` matrix."
  [pos-semidefinite-apache-m accu]
  (if (zero? (rows pos-semidefinite-apache-m))
    {::rectangular-root pos-semidefinite-apache-m
     ::rank             0}
    (try (let [r (RectangularCholeskyDecomposition.
                   ^Array2DRowRealMatrix pos-semidefinite-apache-m
                   (double accu))]
           {::rectangular-root (.getRootMatrix r)
            ::rank             (.getRank r)})
      (catch Exception e
        {::anomalies/message  (.getMessage e)
         ::anomalies/fn       (var rectangular-cholesky-decomposition)
         ::anomalies/category ::anomalies/third-party}))))

(s/fdef rectangular-cholesky-decomposition
  :args (s/cat :pos-semidefinite-apache-m ::pos-semidefinite-apache-matrix-finite
          :accu ::m/accu)
  :ret (s/or :anomaly ::anomalies/anomaly
         :res (s/keys :req [::rectangular-root ::rank])))

(s/def ::svd-left ::apache-matrix)
(s/def ::singular-values ::diagonal-apache-matrix)
(s/def ::svd-right ::apache-matrix)

(defn sv-decomposition
  "Calculates the compact Singular Value Decomposition of an Apache Commons
  matrix. The Singular Value Decomposition of `apache-m` is a set of three
   matrices: `svd-left`, `singular-values`, and `svd-right` such that:
    `apache-m` = `svd-left` × `singular-values` × `svd-right`.
   Let `apache-m` be a m × n matrix, then `svd-left` is a m × p orthogonal
   matrix of the left singular vectors, `singular-values` is a p × p diagonal
   matrix of singular values with positive or nil elements, and are ordered from
   largest to smallest. `svd-right` is a p × n orthogonal matrix of the right
   singular vectors where p = min(m,n). Note that:
   Identity Matrix = (transpose `svd-left`) × `svd-left` =
    `svd-right` × (transpose `svd-right`).
   Returns a map containing:
      `::svd-left` -- Apache Commons matrix of left singular vectors
      `::singular-values` -- diagonal Apache Commons matrix
      `::svd-right` -- transpose of Apache Commons matrix of right singular
        vectors
      `::rank` -- rank."
  [apache-m]
  (if (zero? (rows apache-m))
    {::svd-left        apache-m
     ::singular-values apache-m
     ::svd-right       apache-m
     ::rank            0}
    (let [d (SingularValueDecomposition. ^Array2DRowRealMatrix apache-m)]
      {::svd-left        (.getU d)
       ::singular-values (.getS d)
       ::svd-right       (.getVT d)
       ::rank            (.getRank d)})))

(s/fdef sv-decomposition
  :args (s/cat :apache-m ::apache-matrix)
  :ret (s/keys :req [::svd-left ::singular-values ::svd-right ::rank]))

(defn condition
  "The `singular-values-apache-matrix` is the diagonal matrix of
  `::singular-values` from [[sv-decomposition]]. Returns the norm2 condition
  number, which is the maximum element value from the
  `singular-values-apache-matrix` divided by the minimum element value."
  [singular-values-apache-matrix]
  (if (zero? (rows singular-values-apache-matrix))
    m/nan
    (let [vs (diagonal singular-values-apache-matrix)]
      (m/div (apply max vs) (apply min vs) m/nan))))

(s/fdef condition
  :args (s/cat :singular-values-apache-matrix ::singular-values)
  :ret ::m/number)

(s/def ::Q ::apache-matrix)
(s/def ::R ::apache-matrix)

(s/def ::LLS-solution
  (s/or :apache-m ::apache-matrix
    :anomaly ::anomalies/anomaly))

(s/def ::error (s/nilable ::symmetric-apache-matrix))

(defn qr-decomposition-with-LLS-and-error-matrix
  "Returns a map containing:
    `::Q` -- orthogonal factors of `apache-m1`
    `::R` -- the upper triangular factors of `apache-m2` (not necessarily
      square)
    `::LLS-solution` -- Apache Commons matrix with linear least squares solution
    `::error` -- Apache Commons matrix of errors."
  [apache-m1 apache-m2]
  (if (zero? (rows apache-m1))
    {::Q            apache-m1
     ::R            apache-m1
     ::LLS-solution apache-m1
     ::error        apache-m1}
    (let [d (QRDecomposition. ^Array2DRowRealMatrix apache-m1)
          ^DecompositionSolver s (.getSolver d)
          var-fn (var qr-decomposition-with-LLS-and-error-matrix)
          solution (try (block-apache-matrix->apache-matrix
                          (.solve s ^Array2DRowRealMatrix apache-m2))
                     (catch Exception e
                       {::anomalies/message  (.getMessage e)
                        ::anomalies/fn       var-fn
                        ::anomalies/category ::anomalies/third-party}))
          r (.getR d)
          r-rows (rows r)
          r-columns (columns r)
          error (when (and (apache-matrix? solution) (>= r-rows r-columns))
                  (let [trimmed-r (get-slices-as-matrix
                                    r
                                    {::mx/exception-row-indices
                                     (range r-columns r-rows)})
                        ri (inverse trimmed-r)]
                    (mx* ri (transpose ri))))]
      {::Q            (.getQ d)
       ::R            r
       ::LLS-solution solution
       ::error        error})))

(s/fdef qr-decomposition-with-LLS-and-error-matrix
  :args (s/cat :apache-m1 ::apache-matrix :apache-m2 ::apache-matrix)
  :ret (s/keys :req [::Q ::R ::LLS-solution ::error]))

(defn qr-decomposition-with-linear-least-squares
  "Returns a map containing:
    `::Q` -- orthogonal factors of `apache-m1`
    `::R` -- the upper triangular factors of `apache-m2`
    `::linear-least-squares` -- Apache Commons matrix with linear least squares
      solution."
  [apache-m1 apache-m2]
  (if (zero? (rows apache-m1))
    {::Q            apache-m1
     ::R            apache-m1
     ::LLS-solution apache-m1}
    (let [d (QRDecomposition. ^Array2DRowRealMatrix apache-m1)
          ^DecompositionSolver s (.getSolver d)
          var-fn (var qr-decomposition-with-linear-least-squares)
          solution (try (block-apache-matrix->apache-matrix
                          (.solve s ^Array2DRowRealMatrix apache-m2))
                     (catch Exception e
                       {::anomalies/message  (.getMessage e)
                        ::anomalies/fn       var-fn
                        ::anomalies/category ::anomalies/third-party}))]
      {::Q            (.getQ d)
       ::R            (.getR d)
       ::LLS-solution solution})))

(s/fdef qr-decomposition-with-linear-least-squares
  :args (s/cat :apache-m1 ::apache-matrix :apache-m2 ::apache-matrix)
  :ret (s/keys :req [::Q ::R ::LLS-solution]))

(defn qr-decomposition
  "Computes QR decomposition of a matrix.
  
  Decomposes a matrix A into A = Q * R, where Q is an orthogonal matrix
  (Q^T * Q = I) and R is an upper triangular matrix. Useful for solving
  linear least squares problems and computing matrix rank.
  
  Parameters:
    apache-m - matrix to decompose (any dimensions)
  
  Returns map with keys:
    ::Q - orthogonal matrix (columns are orthonormal)
    ::R - upper triangular matrix
  
  Example:
    (qr-decomposition (->apache-matrix [[1 2] [3 4] [5 6]]))
    ;=> Q has orthonormal columns, R is upper triangular"
  [apache-m]
  (if (zero? (rows apache-m))
    {::Q apache-m
     ::R apache-m}
    (let [d (QRDecomposition. ^Array2DRowRealMatrix apache-m)]
      {::Q (.getQ d)
       ::R (.getR d)})))

(s/fdef qr-decomposition
  :args (s/cat :apache-m ::apache-matrix)
  :ret (s/keys :req [::Q ::R]))

(s/def ::RRQR-permutation ::apache-matrix)
(defn rank-revealing-qr-decomposition
  "Calculates the rank-revealing QR-decomposition of a matrix, with column
  pivoting. The rank-revealing QR-decomposition of `apache-m` consists of three
  matrices `Q`, `R`, and `RRQR-permutation` such that:
  `apache-m` × `permutation` = `Q` × `R`. `Q` is orthogonal (QT × `Q` = I), and
  `R` is upper triangular. If `apache-m` is m × n, `Q` is m × m, `R` is m × n,
  and `RRQR-permutation` is n × n. QR decomposition with column pivoting
  produces a rank-revealing QR decomposition. This class computes the
  decomposition using Householder reflectors. Returns a map containing:
      `::Q` -- orthogonal factors
      `::R` -- the upper triangular factors
      `::RRQR-permutation` -- Permutation Matrix
      `::rank` -- the rank."
  [apache-m accu]
  (if (zero? (rows apache-m))
    {::Q                apache-m
     ::R                apache-m
     ::RRQR-permutation apache-m
     ::rank             0}
    (let [d (RRQRDecomposition.
              ^Array2DRowRealMatrix apache-m
              (double accu))]
      {::Q                (.getQ d)
       ::R                (.getR d)
       ::RRQR-permutation (.getP d)
       ::rank             (.getRank d accu)})))

(s/fdef rank-revealing-qr-decomposition
  :args (s/cat :apache-m ::apache-matrix :accu ::m/accu)
  :ret (s/keys :req [::Q ::R ::RRQR-permutation ::rank]))
