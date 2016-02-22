;; gorilla-repl.fileformat = 1

;; **
;;; #HYBRID WITH VARIABLE DEGREES OF FREEDOM
;;; 
;;; The mma code has all been moved into LagrangeSolverCompact.m, the Java class is LagrangianScore.java. I have followed Jony's idea of memoizing the score function. This dramatically improves speed. What this does is that it keeps a record of previous calls to the score function and the corresponding result. So if it calls with the same arguments etc then it is able to find the answer via some sort of look up table.
;;; 
;;; It also saves each zeitgeist to disk - every generation.
;; **

;; @@
(ns goliath
  (:require [gorilla-plot.core :as plot]
            [clojure.repl :as repl]
            [clojure.pprint :as pprint]
            [clojure.walk :as walk]
            [clojure.zip :as zip]
            [gorilla-repl.table :as table]
            [gorilla-repl.html :as html]
    		[darwin.core :as darwin]
            [darwin.evolution.metrics :as metrics]
            [algebolic.expression.core :as expression]
            [algebolic.expression.tree :as tree]
            [algebolic.expression.genetics :as genetics]
            [algebolic.expression.score :as score]
            [algebolic.expression.render :as spea-render]
            [algebolic.expression.interpreter :as interpreter]
            [darwin.evolution.core :as evolution]
            [darwin.evolution.metrics :as metrics]
            [darwin.evolution.reproduction :as reproduction]
            [darwin.evolution.scoring :as scoring]
            [darwin.evolution.selection :as selection]
            [darwin.evolution.transform :as transform]
            [darwin.evolution.pareto :as pareto]
            [darwin.algorithms.spea2 :as spea2]
            [criterium.core :as criterium])
  (:import [com.wolfram.jlink]
      [goliath.mathlink])
  
  
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; SET THE PATH TO THE MATHEMATICA KERNEL
;; **

;; @@
(def mathKernelSz "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mathKernelSz</span>","value":"#'goliath/mathKernelSz"}
;; <=

;; **
;;; CHOOSE THE DATA - NEEDS TO BE IN RESOURCES/
;; **

;; @@
(def experimentalDataSz "mma_sho1.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/experimentalDataSz</span>","value":"#'goliath/experimentalDataSz"}
;; <=

;; **
;;; CHOOSE THE TIME INTERVAL
;; **

;; @@
(def dt 0.1)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/dt</span>","value":"#'goliath/dt"}
;; <=

;; **
;;; 
;;; 
;;; CHOOSE THE NUMBER OF DEGREES OF FREEDOM (in this case 2 for the double prendulum)
;; **

;; @@
(def df 1)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/df</span>","value":"#'goliath/df"}
;; <=

;; @@
(def prob_inheritance 0.75) ;;to do with crossover  ;probability that terms are directly inherited from parent to child.
(def mutate_pref1 0.8)
(def mutate_pref2 0.8);;change 1 df within a gene  ;prob. that "gene-replace"/"gene-add" is called when mutate is called.
(def mutate_pref3 0.2) ;; add a whole new gene/replace a term (big change)

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mutate_pref3</span>","value":"#'goliath/mutate_pref3"}
;; <=

;; @@
(def power_sum 10)    ;the maximum sum of powers present in each term within a polynomial.
(def poly_length 10)  ;the maximum number of terms of a polynomial. 
(def power 8) 		  ;the maximum power a variable can be raised to.
(def term_length (* 2 df))   ;the number of variables within each term - 2*degree_freedom
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/term_length</span>","value":"#'goliath/term_length"}
;; <=

;; @@
(defn new-poly-term 
  
  "creates a random new poly term satisfying given constraints
   powermax - the maximum power of any given variable, 
   spmax - the max sum of the powers and length is the length of the gene to be generated. 
   spmax - checked and the function re-called if the gene is invalid"
  
  [max_pow spmax length]
  
  (let [new-gene (vec (repeatedly length #(rand-int (+ 1 max_pow)))) sp (apply + new-gene)]
    
    (if (or (< spmax sp) (= new-gene (vec (replicate length 0)) )) 
      
      (assoc (vec (repeat length 0)) (rand-int length) 1)
    
      ;(new-poly-term max_pow spmax length)
      ;[0 1];;cant have this or else it will crash the mathematica
      new-gene
      )
  ))
  
  
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/new-poly-term</span>","value":"#'goliath/new-poly-term"}
;; <=

;; @@


(defn create-genotype
  
  "creates a new individual of length between 1 and max_length decided at random.
  terms within an individual are of legnth n.
  polynomial terms are constrained to maximum power p_max"
  
  [max_length max_pow spmax term_length]
  
  (let [geno (vec (repeatedly (+ 1 (rand-int (- max_length 1))) #(new-poly-term max_pow spmax term_length)))]
    
      (vec (distinct geno))
     )    
    )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/create-genotype</span>","value":"#'goliath/create-genotype"}
;; <=

;; @@
(defn distribute-parent   
  "Distributes the alleses of a parent to children with seeded probability"
 ([parent-1 prob]
  (let [{c1 true c2 false} (group-by (fn [x] (<= prob (rand))) parent-1)]    
[c1 c2]  
   ))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/distribute-parent</span>","value":"#'goliath/distribute-parent"}
;; <=

;; **
;;; 
;; **

;; **
;;; ##Cross-over and mutate functions
;; **

;; **
;;; cross-over, gene-mutate, gene-replace and gene-tweak
;; **

;; @@
(defn cross-over
  
  "takes two individuals and crosses over the polynomial from one parent to another under a weighted selection bias" 
  
  [parent-1 parent-2] 
  
  (let [offspring (mapv vec (mapv set (mapv concat (distribute-parent parent-1 prob_inheritance) (reverse (distribute-parent parent-2 prob_inheritance)))))]
    
    (if (or (zero? (count (offspring 0))) (zero? (count (offspring 1))))
      	
      [parent-1 parent-2]
      [(offspring 0)(offspring 1)]
      )
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/cross-over</span>","value":"#'goliath/cross-over"}
;; <=

;; @@
(defn gene-replace
  
  "Replaces a gene from an individual (with a unique gene). In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  
  [indv powermax spmax n]
  
(let [length (count (first indv)) new-gene (new-poly-term powermax spmax length) new-geno (assoc indv (rand-int (count indv)) new-gene) sp (apply + new-gene)]
  
  (if (> n 10) indv  ;;
  (if (= (count (set new-geno)) (count new-geno))
    new-geno
    (gene-replace indv powermax spmax (inc n))
 
 ))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-replace</span>","value":"#'goliath/gene-replace"}
;; <=

;; @@
(defn gene-add
  
  "Adds a new gene to a GA individual. The gene will be unique addition to the individual. In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  
  [indiv powermax spmax n]
   
  (let [length (count (first indiv)) new-gene (new-poly-term powermax spmax length) new-geno (vec (conj (set indiv) new-gene)) sp (apply + new-gene)]
   
      (if (> n 10) 
        indiv (if (= (count new-geno) (count indiv))
      (gene-add indiv powermax spmax (inc n))
     new-geno
      ))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-add</span>","value":"#'goliath/gene-add"}
;; <=

;; @@
(defn gene-tweak
  [indv]
  (let [n (rand-int (count indv)) gene (nth indv n) new-gene (assoc gene (rand-int (count gene)) (rand-int power))]
    (assoc indv n new-gene)
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-tweak</span>","value":"#'goliath/gene-tweak"}
;; <=

;; @@
(defn gene-delete [indv] 
  (let [new-indv (subvec indv 1 (count indv))]
    (if (= new-indv []) indv
      new-indv))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-delete</span>","value":"#'goliath/gene-delete"}
;; <=

;; @@
(defn mutate
  "takes two genotypes and returns a mutated genotype"
  [indv]
  (let [rn (rand)]
    (if (>= rn mutate_pref1)
      (gene-delete indv)
    (if (>= rn mutate_pref2)
        (gene-tweak indv)
          (if (>= rn mutate_pref3) 
          (gene-replace indv power power_sum 1)
          (gene-add indv power power_sum 1))
      )
      )
    )
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mutate</span>","value":"#'goliath/mutate"}
;; <=

;; **
;;; #Initial populations and generation-configuration
;; **

;; **
;;; plus a little kernal opening and closing
;; **

;; @@
(defn random-initial-population 
  
   "creates an innitial population of polynimials
  size - number of individuals in the population"
   
  [size max_length max_power max_sum_powers num_terms]
 
  (repeatedly size #(create-genotype max_length max_power max_sum_powers num_terms)) 
       
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/random-initial-population</span>","value":"#'goliath/random-initial-population"}
;; <=

;; **
;;; Save the population off to disk after each generation - so that possible good L can be examined later.
;;; 
;; **

;; @@
(defn SaveZ [zg str]
  (spit str (pr-str zg))
  
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/SaveZ</span>","value":"#'goliath/SaveZ"}
;; <=

;; @@
(defn task [zg]
  ;;do something
  (print (:age zg))
  (SaveZ zg "Coupled_SHO_Run1.txt")
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/task</span>","value":"#'goliath/task"}
;; <=

;; @@
(def initial-zeitgeist (evolution/make-zeitgeist (random-initial-population 100 poly_length power power_sum term_length)))

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/initial-zeitgeist</span>","value":"#'goliath/initial-zeitgeist"}
;; <=

;; @@
(goliath.mathlink.LagrangianScore/InitFunctions
  mathKernelSz
  "resources/"
  experimentalDataSz
  dt
  df
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(defn score
  [indv] 
	(first
      (into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (vec (flatten indv))), df))
      )
  
  
  )
  

;;(into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (vec (flatten [[2 2 0 ;;0] [0 2 0 0] [1 1 0 0] [0 0 2 0] [0 0 0 2]]))), df))


;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/score</span>","value":"#'goliath/score"}
;; <=

;; @@
(def memScore (memoize score))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/memScore</span>","value":"#'goliath/memScore"}
;; <=

;; @@
(def generation-config			      
  (let [size-limit 120                   
        min-size 20
        ea-config (spea2/spea2-config
                    {:goals [:error :complexity]
                     :archive-size 100
                     :binary-ops [{:op cross-over :repeat 40}]
                     :unary-ops [{:op mutate :repeat 20}]})
        score-functions {:complexity (fn [x] (count x))
                         :error (fn [e] (memScore e))}]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/generation-config</span>","value":"#'goliath/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (task zg) (>= (:age zg) 200)))))
;; @@
;; ->
;;; 0.1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22.23.24.25.26.27.28.29.30.31.32.33.34.35.36.37.38.39.40.41.42.43.44.45.46.47.48.49.50.51.52.53.54.55.56.57.58.59.60.61.62.63.64.65.66.67.68.69.70.71.72.73.74.75.76.77.78.79.80.81.82.83.84.85.86.87.88.89.90.91.92.93.94.95.96.97.98.99.100.101.102.103.104.105.106.107.108.109.110.111.112.113.114.115.116.117.118.119.120.121.122.123.124.125.126.127.128.129.130.131.132.133.134.135.136.137.138.139.140.141.142.143.144.145.146.147.148.149.150.151.152.153.154.155.156.157.158.159.160.161.162.163.164.165.166.167.168.169.170.171.172.173.174.175.176.177.178.179.180.181.182.183.184.185.186.187.188.189.190.191.192.193.194.195.196.197.198.199.200&quot;Elapsed time: 208604.301569 msecs&quot;
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/result</span>","value":"#'goliath/result"}
;; <=

;; @@
(goliath.mathlink.LagrangianScore/Shutdown)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
;;(def result (read-string (slurp "results16LScore.txt") ))
;; @@

;; @@
(mapv #(println (:genotype %)) (sort-by :error (:elite result)))
;; @@
;; ->
;;; [[1 0] [3 0] [1 4] [6 0] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [6 0] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [6 0] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [6 0] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [1 4] [5 0] [1 2] [3 2]]
;;; [[1 0] [3 0] [0 2] [2 0] [1 2]]
;;; [[1 0] [3 0] [0 2] [2 0] [1 2]]
;;; [[1 0] [3 0] [0 2] [2 0] [1 2]]
;;; [[1 0] [3 0] [0 2] [2 0] [1 2]]
;;; [[1 0] [3 0] [1 4] [5 0]]
;;; [[1 0] [3 0] [1 4] [5 0]]
;;; [[1 0] [3 0] [1 4] [5 0]]
;;; [[1 0] [3 0] [1 4] [5 0]]
;;; [[7 2] [0 2] [2 0]]
;;; [[7 2] [0 2] [2 0]]
;;; [[7 2] [0 2] [2 0]]
;;; [[7 2] [0 2] [2 0]]
;;; [[0 2] [2 0]]
;;; [[0 2] [2 0]]
;;; [[0 2] [2 0]]
;;; [[0 2] [2 0]]
;;; [[0 2] [2 0]]
;;; [[0 2] [2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; [[2 0]]
;;; 
;; <-
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}],"value":"[nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil]"}
;; <=

;; **
;;; #Pareto Analysis
;;; 
;; **

;; **
;;; s
;; **

;; @@
(plot/list-plot (:min (:error @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"488983b2-8be6-4c34-a984-ae08455e8145","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"488983b2-8be6-4c34-a984-ae08455e8145","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"488983b2-8be6-4c34-a984-ae08455e8145"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"488983b2-8be6-4c34-a984-ae08455e8145","values":[{"x":0,"y":-12.693998876964283},{"x":1,"y":-13.217067274841154},{"x":2,"y":-13.217067274841154},{"x":3,"y":-13.217067274841154},{"x":4,"y":-13.223890324589698},{"x":5,"y":-13.223890324589698},{"x":6,"y":-14.208967812775947},{"x":7,"y":-14.208967812775947},{"x":8,"y":-14.224288720862294},{"x":9,"y":-14.224288720862294},{"x":10,"y":-14.242413452390386},{"x":11,"y":-14.242413452390386},{"x":12,"y":-14.242918728168783},{"x":13,"y":-14.248867016171468},{"x":14,"y":-14.248867016171468},{"x":15,"y":-14.248867016171468},{"x":16,"y":-14.248867016171468},{"x":17,"y":-14.248867016171468},{"x":18,"y":-14.37337508329446},{"x":19,"y":-22.438641839132654},{"x":20,"y":-22.438641839132654},{"x":21,"y":-23.04492196088355},{"x":22,"y":-23.04508164282913},{"x":23,"y":-23.04508164282913},{"x":24,"y":-23.04508164282913},{"x":25,"y":-23.04508164282913},{"x":26,"y":-23.045149139126618},{"x":27,"y":-23.045149139126618},{"x":28,"y":-23.045149139126618},{"x":29,"y":-23.045149139126618},{"x":30,"y":-23.045149139126618},{"x":31,"y":-23.045149139126618},{"x":32,"y":-23.045149139126618},{"x":33,"y":-23.045149139126618},{"x":34,"y":-23.045149139126618},{"x":35,"y":-23.045149139126618},{"x":36,"y":-23.045149139126618},{"x":37,"y":-23.045149139126618},{"x":38,"y":-23.045149139126618},{"x":39,"y":-23.161866730884505},{"x":40,"y":-23.161866730884505},{"x":41,"y":-23.161866730884505},{"x":42,"y":-23.161866730884505},{"x":43,"y":-23.161866730884505},{"x":44,"y":-23.161884246698172},{"x":45,"y":-23.161884246698172},{"x":46,"y":-23.161884246698172},{"x":47,"y":-23.161991310244147},{"x":48,"y":-23.161991310244147},{"x":49,"y":-23.161991310244147},{"x":50,"y":-23.162060882422452},{"x":51,"y":-23.162060882422452},{"x":52,"y":-23.162060882422452},{"x":53,"y":-23.162060882422452},{"x":54,"y":-23.162060882422452},{"x":55,"y":-23.162060882422452},{"x":56,"y":-23.162060882422452},{"x":57,"y":-23.162060882422452},{"x":58,"y":-23.162060882422452},{"x":59,"y":-23.162060882422452},{"x":60,"y":-23.162060882422452},{"x":61,"y":-23.162060882422452},{"x":62,"y":-23.162060882422452},{"x":63,"y":-23.162060882422452},{"x":64,"y":-23.162060882422452},{"x":65,"y":-23.16500825082244},{"x":66,"y":-23.16500825082244},{"x":67,"y":-23.165018331519192},{"x":68,"y":-23.165018331519192},{"x":69,"y":-23.16984797101876},{"x":70,"y":-23.16984797101876},{"x":71,"y":-23.16984797101876},{"x":72,"y":-23.16984797101876},{"x":73,"y":-23.16984797101876},{"x":74,"y":-23.16984797101876},{"x":75,"y":-23.16984797101876},{"x":76,"y":-23.16984797101876},{"x":77,"y":-23.16984797101876},{"x":78,"y":-23.16984797101876},{"x":79,"y":-23.16984797101876},{"x":80,"y":-23.16984797101876},{"x":81,"y":-23.16984797101876},{"x":82,"y":-23.16984797101876},{"x":83,"y":-23.16984797101876},{"x":84,"y":-23.16984797101876},{"x":85,"y":-23.16984797101876},{"x":86,"y":-23.16984797101876},{"x":87,"y":-23.16984797101876},{"x":88,"y":-23.16984797101876},{"x":89,"y":-23.16984797101876},{"x":90,"y":-23.16984797101876},{"x":91,"y":-23.16984797101876},{"x":92,"y":-23.16984797101876},{"x":93,"y":-23.16984797101876},{"x":94,"y":-23.16984797101876},{"x":95,"y":-23.16984797101876},{"x":96,"y":-23.16984797101876},{"x":97,"y":-23.16984797101876},{"x":98,"y":-23.16984797101876},{"x":99,"y":-23.16984797101876},{"x":100,"y":-23.16984797101876},{"x":101,"y":-23.16984797101876},{"x":102,"y":-23.16984797101876},{"x":103,"y":-23.16984797101876},{"x":104,"y":-23.16984797101876},{"x":105,"y":-23.16984797101876},{"x":106,"y":-23.16984797101876},{"x":107,"y":-23.16984797101876},{"x":108,"y":-23.16984797101876},{"x":109,"y":-23.16984797101876},{"x":110,"y":-23.16984797101876},{"x":111,"y":-23.16984797101876},{"x":112,"y":-24.305419824988324},{"x":113,"y":-24.305419824988324},{"x":114,"y":-24.305419824988324},{"x":115,"y":-24.305419824988324},{"x":116,"y":-24.305419824988324},{"x":117,"y":-24.305419824988324},{"x":118,"y":-24.305419824988324},{"x":119,"y":-24.305419824988324},{"x":120,"y":-24.305419824988324},{"x":121,"y":-24.305419824988324},{"x":122,"y":-24.305419824988324},{"x":123,"y":-24.305419824988324},{"x":124,"y":-24.305419824988324},{"x":125,"y":-24.305419824988324},{"x":126,"y":-24.305419824988324},{"x":127,"y":-24.305419824988324},{"x":128,"y":-24.305419824988324},{"x":129,"y":-24.305419824988324},{"x":130,"y":-24.305419824988324},{"x":131,"y":-24.305419824988324},{"x":132,"y":-24.305419824988324},{"x":133,"y":-24.305419824988324},{"x":134,"y":-24.305419824988324},{"x":135,"y":-24.305419824988324},{"x":136,"y":-24.305419824988324},{"x":137,"y":-24.305419824988324},{"x":138,"y":-24.305419824988324},{"x":139,"y":-24.305419824988324},{"x":140,"y":-24.305419824988324},{"x":141,"y":-24.305419824988324},{"x":142,"y":-24.305419824988324},{"x":143,"y":-24.305419824988324},{"x":144,"y":-24.305419824988324},{"x":145,"y":-24.305419824988324},{"x":146,"y":-24.305419824988324},{"x":147,"y":-24.305419824988324},{"x":148,"y":-24.305419824988324},{"x":149,"y":-24.305419824988324},{"x":150,"y":-24.305419824988324},{"x":151,"y":-24.305419824988324},{"x":152,"y":-24.305419824988324},{"x":153,"y":-24.305419824988324},{"x":154,"y":-24.305419824988324},{"x":155,"y":-24.305419824988324},{"x":156,"y":-24.305419824988324},{"x":157,"y":-24.305419824988324},{"x":158,"y":-24.305419824988324},{"x":159,"y":-24.305419824988324},{"x":160,"y":-24.305419824988324},{"x":161,"y":-24.305419824988324},{"x":162,"y":-24.305419824988324},{"x":163,"y":-24.305419824988324},{"x":164,"y":-24.305419824988324},{"x":165,"y":-24.305419824988324},{"x":166,"y":-24.305419824988324},{"x":167,"y":-24.305419824988324},{"x":168,"y":-24.305419824988324},{"x":169,"y":-24.305419824988324},{"x":170,"y":-24.305419824988324},{"x":171,"y":-24.305419824988324},{"x":172,"y":-24.305419824988324},{"x":173,"y":-24.305419824988324},{"x":174,"y":-24.305419824988324},{"x":175,"y":-24.305419824988324},{"x":176,"y":-24.305419824988324},{"x":177,"y":-24.305419824988324},{"x":178,"y":-24.305419824988324},{"x":179,"y":-24.305419824988324},{"x":180,"y":-24.305419824988324},{"x":181,"y":-24.305419824988324},{"x":182,"y":-24.305419824988324},{"x":183,"y":-24.305419824988324},{"x":184,"y":-24.305419824988324},{"x":185,"y":-24.305419824988324},{"x":186,"y":-24.305419824988324},{"x":187,"y":-24.305419824988324},{"x":188,"y":-24.305419824988324},{"x":189,"y":-24.305419824988324},{"x":190,"y":-24.305419824988324},{"x":191,"y":-24.305419824988324},{"x":192,"y":-24.305419824988324},{"x":193,"y":-24.305419824988324},{"x":194,"y":-24.305419824988324},{"x":195,"y":-24.305419824988324},{"x":196,"y":-24.305419824988324},{"x":197,"y":-24.305419824988324},{"x":198,"y":-24.305419824988324},{"x":199,"y":-24.305419824988324}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"488983b2-8be6-4c34-a984-ae08455e8145\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"488983b2-8be6-4c34-a984-ae08455e8145\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"488983b2-8be6-4c34-a984-ae08455e8145\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"488983b2-8be6-4c34-a984-ae08455e8145\", :values ({:x 0, :y -12.693998876964283} {:x 1, :y -13.217067274841154} {:x 2, :y -13.217067274841154} {:x 3, :y -13.217067274841154} {:x 4, :y -13.223890324589698} {:x 5, :y -13.223890324589698} {:x 6, :y -14.208967812775947} {:x 7, :y -14.208967812775947} {:x 8, :y -14.224288720862294} {:x 9, :y -14.224288720862294} {:x 10, :y -14.242413452390386} {:x 11, :y -14.242413452390386} {:x 12, :y -14.242918728168783} {:x 13, :y -14.248867016171468} {:x 14, :y -14.248867016171468} {:x 15, :y -14.248867016171468} {:x 16, :y -14.248867016171468} {:x 17, :y -14.248867016171468} {:x 18, :y -14.37337508329446} {:x 19, :y -22.438641839132654} {:x 20, :y -22.438641839132654} {:x 21, :y -23.04492196088355} {:x 22, :y -23.04508164282913} {:x 23, :y -23.04508164282913} {:x 24, :y -23.04508164282913} {:x 25, :y -23.04508164282913} {:x 26, :y -23.045149139126618} {:x 27, :y -23.045149139126618} {:x 28, :y -23.045149139126618} {:x 29, :y -23.045149139126618} {:x 30, :y -23.045149139126618} {:x 31, :y -23.045149139126618} {:x 32, :y -23.045149139126618} {:x 33, :y -23.045149139126618} {:x 34, :y -23.045149139126618} {:x 35, :y -23.045149139126618} {:x 36, :y -23.045149139126618} {:x 37, :y -23.045149139126618} {:x 38, :y -23.045149139126618} {:x 39, :y -23.161866730884505} {:x 40, :y -23.161866730884505} {:x 41, :y -23.161866730884505} {:x 42, :y -23.161866730884505} {:x 43, :y -23.161866730884505} {:x 44, :y -23.161884246698172} {:x 45, :y -23.161884246698172} {:x 46, :y -23.161884246698172} {:x 47, :y -23.161991310244147} {:x 48, :y -23.161991310244147} {:x 49, :y -23.161991310244147} {:x 50, :y -23.162060882422452} {:x 51, :y -23.162060882422452} {:x 52, :y -23.162060882422452} {:x 53, :y -23.162060882422452} {:x 54, :y -23.162060882422452} {:x 55, :y -23.162060882422452} {:x 56, :y -23.162060882422452} {:x 57, :y -23.162060882422452} {:x 58, :y -23.162060882422452} {:x 59, :y -23.162060882422452} {:x 60, :y -23.162060882422452} {:x 61, :y -23.162060882422452} {:x 62, :y -23.162060882422452} {:x 63, :y -23.162060882422452} {:x 64, :y -23.162060882422452} {:x 65, :y -23.16500825082244} {:x 66, :y -23.16500825082244} {:x 67, :y -23.165018331519192} {:x 68, :y -23.165018331519192} {:x 69, :y -23.16984797101876} {:x 70, :y -23.16984797101876} {:x 71, :y -23.16984797101876} {:x 72, :y -23.16984797101876} {:x 73, :y -23.16984797101876} {:x 74, :y -23.16984797101876} {:x 75, :y -23.16984797101876} {:x 76, :y -23.16984797101876} {:x 77, :y -23.16984797101876} {:x 78, :y -23.16984797101876} {:x 79, :y -23.16984797101876} {:x 80, :y -23.16984797101876} {:x 81, :y -23.16984797101876} {:x 82, :y -23.16984797101876} {:x 83, :y -23.16984797101876} {:x 84, :y -23.16984797101876} {:x 85, :y -23.16984797101876} {:x 86, :y -23.16984797101876} {:x 87, :y -23.16984797101876} {:x 88, :y -23.16984797101876} {:x 89, :y -23.16984797101876} {:x 90, :y -23.16984797101876} {:x 91, :y -23.16984797101876} {:x 92, :y -23.16984797101876} {:x 93, :y -23.16984797101876} {:x 94, :y -23.16984797101876} {:x 95, :y -23.16984797101876} {:x 96, :y -23.16984797101876} {:x 97, :y -23.16984797101876} {:x 98, :y -23.16984797101876} {:x 99, :y -23.16984797101876} {:x 100, :y -23.16984797101876} {:x 101, :y -23.16984797101876} {:x 102, :y -23.16984797101876} {:x 103, :y -23.16984797101876} {:x 104, :y -23.16984797101876} {:x 105, :y -23.16984797101876} {:x 106, :y -23.16984797101876} {:x 107, :y -23.16984797101876} {:x 108, :y -23.16984797101876} {:x 109, :y -23.16984797101876} {:x 110, :y -23.16984797101876} {:x 111, :y -23.16984797101876} {:x 112, :y -24.305419824988324} {:x 113, :y -24.305419824988324} {:x 114, :y -24.305419824988324} {:x 115, :y -24.305419824988324} {:x 116, :y -24.305419824988324} {:x 117, :y -24.305419824988324} {:x 118, :y -24.305419824988324} {:x 119, :y -24.305419824988324} {:x 120, :y -24.305419824988324} {:x 121, :y -24.305419824988324} {:x 122, :y -24.305419824988324} {:x 123, :y -24.305419824988324} {:x 124, :y -24.305419824988324} {:x 125, :y -24.305419824988324} {:x 126, :y -24.305419824988324} {:x 127, :y -24.305419824988324} {:x 128, :y -24.305419824988324} {:x 129, :y -24.305419824988324} {:x 130, :y -24.305419824988324} {:x 131, :y -24.305419824988324} {:x 132, :y -24.305419824988324} {:x 133, :y -24.305419824988324} {:x 134, :y -24.305419824988324} {:x 135, :y -24.305419824988324} {:x 136, :y -24.305419824988324} {:x 137, :y -24.305419824988324} {:x 138, :y -24.305419824988324} {:x 139, :y -24.305419824988324} {:x 140, :y -24.305419824988324} {:x 141, :y -24.305419824988324} {:x 142, :y -24.305419824988324} {:x 143, :y -24.305419824988324} {:x 144, :y -24.305419824988324} {:x 145, :y -24.305419824988324} {:x 146, :y -24.305419824988324} {:x 147, :y -24.305419824988324} {:x 148, :y -24.305419824988324} {:x 149, :y -24.305419824988324} {:x 150, :y -24.305419824988324} {:x 151, :y -24.305419824988324} {:x 152, :y -24.305419824988324} {:x 153, :y -24.305419824988324} {:x 154, :y -24.305419824988324} {:x 155, :y -24.305419824988324} {:x 156, :y -24.305419824988324} {:x 157, :y -24.305419824988324} {:x 158, :y -24.305419824988324} {:x 159, :y -24.305419824988324} {:x 160, :y -24.305419824988324} {:x 161, :y -24.305419824988324} {:x 162, :y -24.305419824988324} {:x 163, :y -24.305419824988324} {:x 164, :y -24.305419824988324} {:x 165, :y -24.305419824988324} {:x 166, :y -24.305419824988324} {:x 167, :y -24.305419824988324} {:x 168, :y -24.305419824988324} {:x 169, :y -24.305419824988324} {:x 170, :y -24.305419824988324} {:x 171, :y -24.305419824988324} {:x 172, :y -24.305419824988324} {:x 173, :y -24.305419824988324} {:x 174, :y -24.305419824988324} {:x 175, :y -24.305419824988324} {:x 176, :y -24.305419824988324} {:x 177, :y -24.305419824988324} {:x 178, :y -24.305419824988324} {:x 179, :y -24.305419824988324} {:x 180, :y -24.305419824988324} {:x 181, :y -24.305419824988324} {:x 182, :y -24.305419824988324} {:x 183, :y -24.305419824988324} {:x 184, :y -24.305419824988324} {:x 185, :y -24.305419824988324} {:x 186, :y -24.305419824988324} {:x 187, :y -24.305419824988324} {:x 188, :y -24.305419824988324} {:x 189, :y -24.305419824988324} {:x 190, :y -24.305419824988324} {:x 191, :y -24.305419824988324} {:x 192, :y -24.305419824988324} {:x 193, :y -24.305419824988324} {:x 194, :y -24.305419824988324} {:x 195, :y -24.305419824988324} {:x 196, :y -24.305419824988324} {:x 197, :y -24.305419824988324} {:x 198, :y -24.305419824988324} {:x 199, :y -24.305419824988324})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(defn pareto-plot-population
  [[k1 k2] result]
  (let [pareto-front (pareto/non-dominated-individuals [k1 k2] (:elite result))
        coord-extract (fn [i] [(k1 i) (k2 i)])]
    (plot/compose
      (plot/list-plot (map coord-extract (:elite result)) :colour "red")
      (plot/list-plot (map coord-extract (:rabble result)) :colour "blue")
          (plot/list-plot (map coord-extract pareto-front) :colour "#ff29d2")
  )))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/pareto-plot-population</span>","value":"#'goliath/pareto-plot-population"}
;; <=

;; @@
(pareto-plot-population [:error :complexity] result)
;; @@
;; =>
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"b5061890-573e-439c-b8bf-504b24b0e69a","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"b5061890-573e-439c-b8bf-504b24b0e69a","field":"data.y"}}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"b5061890-573e-439c-b8bf-504b24b0e69a","values":[{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-14.167399471668798,"y":4},{"x":-14.167399471668798,"y":4},{"x":-14.167399471668798,"y":4},{"x":-14.167399471668798,"y":4},{"x":-23.04492459894466,"y":6},{"x":-23.04492459894466,"y":6},{"x":-23.04492459894466,"y":6},{"x":-23.04492459894466,"y":6},{"x":-24.305419824988324,"y":7},{"x":-24.305419824988324,"y":7},{"x":-24.305419824988324,"y":7},{"x":-24.305419824988324,"y":7},{"x":-17.12456096156009,"y":5},{"x":-17.12456096156009,"y":5},{"x":-17.12456096156009,"y":5},{"x":-17.12456096156009,"y":5}]},{"name":"b7b7b894-fa2c-4fd0-a50d-c465778c8d8f","values":[{"x":-12.71896527479115,"y":5},{"x":1.063751358701247,"y":1},{"x":-13.447022880190035,"y":5},{"x":16.212236546842988,"y":2},{"x":0.4726297634946417,"y":1},{"x":-13.702330518074435,"y":3},{"x":15.809942462408959,"y":1},{"x":-0.2875411879323268,"y":2},{"x":0.4726297634946417,"y":1},{"x":10.934856417758798,"y":1},{"x":-0.4130076033079347,"y":4},{"x":-13.095138264389364,"y":7},{"x":-0.5729665423909693,"y":4},{"x":0.4726297634946417,"y":1},{"x":9.218823504075301,"y":1},{"x":0.06392783211977582,"y":2},{"x":1.5832220019188459,"y":1},{"x":-17.12456095967027,"y":6},{"x":1.5832220019188459,"y":1},{"x":9.218823504075301,"y":1},{"x":-0.10654889707454009,"y":2},{"x":-0.366962501675893,"y":3},{"x":-12.71896527479115,"y":5},{"x":-12.900038528290512,"y":4},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-0.12535219486350394,"y":3},{"x":-13.44407568640285,"y":5},{"x":-23.04492459894466,"y":6},{"x":-0.3314505329385039,"y":3},{"x":-23.04492459894466,"y":6},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":6.980877274384624,"y":2},{"x":-1.3775555968132558,"y":4},{"x":-13.716495979027698,"y":3},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.44025277883972497,"y":2},{"x":-13.702330516612514,"y":2},{"x":-13.72315454174506,"y":4},{"x":9.218823504075301,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-23.0449752372222,"y":7},{"x":1.063751358701247,"y":1},{"x":0.4726297634946417,"y":1},{"x":16.212236546842988,"y":1},{"x":-12.722428593802828,"y":5},{"x":-0.26415263948576684,"y":2},{"x":-13.702330516612514,"y":2},{"x":0.44025277883972497,"y":2},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":-14.167399471668798,"y":4},{"x":-0.6362441107207915,"y":3},{"x":-0.2875411879323268,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":6.998525253793526,"y":1},{"x":-0.6395412402714555,"y":4},{"x":-0.2875411879323268,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-12.74565134057376,"y":6},{"x":-0.2875411879323268,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-0.08391604643222912,"y":2},{"x":-0.41000188643242796,"y":3},{"x":-13.066509658824884,"y":5},{"x":0.0782210669678484,"y":2},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-0.5774665461597778,"y":3},{"x":-13.72315131939051,"y":3},{"x":0.0782210669678484,"y":2},{"x":-13.09322792846534,"y":6},{"x":0.11343839059317351,"y":3},{"x":-13.448355983255915,"y":5},{"x":0.4726297634946417,"y":1},{"x":-14.167399471668798,"y":4},{"x":-15.556564607974659,"y":6},{"x":0.4614795571597382,"y":2},{"x":0.4726297634946417,"y":1},{"x":6.980877274384624,"y":2},{"x":-13.424872457672079,"y":4},{"x":-13.226445947580999,"y":5},{"x":16.212236546842988,"y":1},{"x":0.44025277883972497,"y":2}]},{"name":"22244a15-38fc-47f5-a00b-72c9fb37a118","values":[{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":0.4726297634946417,"y":1},{"x":-13.702330516612514,"y":2},{"x":0.4726297634946417,"y":1},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-14.167399471668798,"y":4},{"x":-14.167399471668798,"y":4},{"x":-14.167399471668798,"y":4},{"x":-14.167399471668798,"y":4},{"x":-23.04492459894466,"y":6},{"x":-23.04492459894466,"y":6},{"x":-23.04492459894466,"y":6},{"x":-23.04492459894466,"y":6},{"x":-24.305419824988324,"y":7},{"x":-24.305419824988324,"y":7},{"x":-24.305419824988324,"y":7},{"x":-24.305419824988324,"y":7},{"x":-17.12456096156009,"y":5},{"x":-17.12456096156009,"y":5},{"x":-17.12456096156009,"y":5},{"x":-17.12456096156009,"y":5}]}],"marks":[{"type":"symbol","from":{"data":"b5061890-573e-439c-b8bf-504b24b0e69a"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"red"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"b7b7b894-fa2c-4fd0-a50d-c465778c8d8f"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"blue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"22244a15-38fc-47f5-a00b-72c9fb37a118"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"#ff29d2"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"b5061890-573e-439c-b8bf-504b24b0e69a\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"b5061890-573e-439c-b8bf-504b24b0e69a\", :field \"data.y\"}}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"b5061890-573e-439c-b8bf-504b24b0e69a\", :values ({:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -14.167399471668798, :y 4} {:x -14.167399471668798, :y 4} {:x -14.167399471668798, :y 4} {:x -14.167399471668798, :y 4} {:x -23.04492459894466, :y 6} {:x -23.04492459894466, :y 6} {:x -23.04492459894466, :y 6} {:x -23.04492459894466, :y 6} {:x -24.305419824988324, :y 7} {:x -24.305419824988324, :y 7} {:x -24.305419824988324, :y 7} {:x -24.305419824988324, :y 7} {:x -17.12456096156009, :y 5} {:x -17.12456096156009, :y 5} {:x -17.12456096156009, :y 5} {:x -17.12456096156009, :y 5})} {:name \"b7b7b894-fa2c-4fd0-a50d-c465778c8d8f\", :values ({:x -12.71896527479115, :y 5} {:x 1.063751358701247, :y 1} {:x -13.447022880190035, :y 5} {:x 16.212236546842988, :y 2} {:x 0.4726297634946417, :y 1} {:x -13.702330518074435, :y 3} {:x 15.809942462408959, :y 1} {:x -0.2875411879323268, :y 2} {:x 0.4726297634946417, :y 1} {:x 10.934856417758798, :y 1} {:x -0.4130076033079347, :y 4} {:x -13.095138264389364, :y 7} {:x -0.5729665423909693, :y 4} {:x 0.4726297634946417, :y 1} {:x 9.218823504075301, :y 1} {:x 0.06392783211977582, :y 2} {:x 1.5832220019188459, :y 1} {:x -17.12456095967027, :y 6} {:x 1.5832220019188459, :y 1} {:x 9.218823504075301, :y 1} {:x -0.10654889707454009, :y 2} {:x -0.366962501675893, :y 3} {:x -12.71896527479115, :y 5} {:x -12.900038528290512, :y 4} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -0.12535219486350394, :y 3} {:x -13.44407568640285, :y 5} {:x -23.04492459894466, :y 6} {:x -0.3314505329385039, :y 3} {:x -23.04492459894466, :y 6} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 6.980877274384624, :y 2} {:x -1.3775555968132558, :y 4} {:x -13.716495979027698, :y 3} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.44025277883972497, :y 2} {:x -13.702330516612514, :y 2} {:x -13.72315454174506, :y 4} {:x 9.218823504075301, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -23.0449752372222, :y 7} {:x 1.063751358701247, :y 1} {:x 0.4726297634946417, :y 1} {:x 16.212236546842988, :y 1} {:x -12.722428593802828, :y 5} {:x -0.26415263948576684, :y 2} {:x -13.702330516612514, :y 2} {:x 0.44025277883972497, :y 2} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x -14.167399471668798, :y 4} {:x -0.6362441107207915, :y 3} {:x -0.2875411879323268, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 6.998525253793526, :y 1} {:x -0.6395412402714555, :y 4} {:x -0.2875411879323268, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -12.74565134057376, :y 6} {:x -0.2875411879323268, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -0.08391604643222912, :y 2} {:x -0.41000188643242796, :y 3} {:x -13.066509658824884, :y 5} {:x 0.0782210669678484, :y 2} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -0.5774665461597778, :y 3} {:x -13.72315131939051, :y 3} {:x 0.0782210669678484, :y 2} {:x -13.09322792846534, :y 6} {:x 0.11343839059317351, :y 3} {:x -13.448355983255915, :y 5} {:x 0.4726297634946417, :y 1} {:x -14.167399471668798, :y 4} {:x -15.556564607974659, :y 6} {:x 0.4614795571597382, :y 2} {:x 0.4726297634946417, :y 1} {:x 6.980877274384624, :y 2} {:x -13.424872457672079, :y 4} {:x -13.226445947580999, :y 5} {:x 16.212236546842988, :y 1} {:x 0.44025277883972497, :y 2})} {:name \"22244a15-38fc-47f5-a00b-72c9fb37a118\", :values ({:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x 0.4726297634946417, :y 1} {:x -13.702330516612514, :y 2} {:x 0.4726297634946417, :y 1} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -14.167399471668798, :y 4} {:x -14.167399471668798, :y 4} {:x -14.167399471668798, :y 4} {:x -14.167399471668798, :y 4} {:x -23.04492459894466, :y 6} {:x -23.04492459894466, :y 6} {:x -23.04492459894466, :y 6} {:x -23.04492459894466, :y 6} {:x -24.305419824988324, :y 7} {:x -24.305419824988324, :y 7} {:x -24.305419824988324, :y 7} {:x -24.305419824988324, :y 7} {:x -17.12456096156009, :y 5} {:x -17.12456096156009, :y 5} {:x -17.12456096156009, :y 5} {:x -17.12456096156009, :y 5})}), :marks ({:type \"symbol\", :from {:data \"b5061890-573e-439c-b8bf-504b24b0e69a\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"red\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"b7b7b894-fa2c-4fd0-a50d-c465778c8d8f\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"blue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"22244a15-38fc-47f5-a00b-72c9fb37a118\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"#ff29d2\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}})}}"}
;; <=

;; **
;;; Timings of the various steps in each generation.
;;; 
;;; - "Red" - selection time
;;; - "Green" Reproduction time 
;;; - "Blue"  Scoring time
;; **

;; @@
(plot/compose
  (plot/list-plot (:time @metrics/metrics) :colour "yellow" :joined true
                  :plot-range [:all [0 (apply max (:time @metrics/metrics))]])
  (plot/list-plot (:selection-time @metrics/metrics) :colour "red" :joined true)
  (plot/list-plot (:reproduction-time @metrics/metrics) :colour "green" :joined true)
  (plot/list-plot (:scoring-time @metrics/metrics) :colour "blue" :joined true))
;; @@
;; =>
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"90342872-f717-4977-b967-3404cba3f04a","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":[0,4269]}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"90342872-f717-4977-b967-3404cba3f04a","values":[{"x":0,"y":3272},{"x":1,"y":4269},{"x":2,"y":3038},{"x":3,"y":1862},{"x":4,"y":1797},{"x":5,"y":2451},{"x":6,"y":1369},{"x":7,"y":852},{"x":8,"y":1486},{"x":9,"y":826},{"x":10,"y":1623},{"x":11,"y":1610},{"x":12,"y":561},{"x":13,"y":1383},{"x":14,"y":2232},{"x":15,"y":1130},{"x":16,"y":1648},{"x":17,"y":1186},{"x":18,"y":2439},{"x":19,"y":1663},{"x":20,"y":727},{"x":21,"y":1086},{"x":22,"y":1326},{"x":23,"y":1064},{"x":24,"y":1091},{"x":25,"y":941},{"x":26,"y":1060},{"x":27,"y":504},{"x":28,"y":437},{"x":29,"y":960},{"x":30,"y":977},{"x":31,"y":1008},{"x":32,"y":958},{"x":33,"y":544},{"x":34,"y":1108},{"x":35,"y":1103},{"x":36,"y":1089},{"x":37,"y":1137},{"x":38,"y":645},{"x":39,"y":932},{"x":40,"y":1348},{"x":41,"y":827},{"x":42,"y":1443},{"x":43,"y":1393},{"x":44,"y":1113},{"x":45,"y":589},{"x":46,"y":472},{"x":47,"y":1187},{"x":48,"y":1178},{"x":49,"y":781},{"x":50,"y":1108},{"x":51,"y":1429},{"x":52,"y":1030},{"x":53,"y":1130},{"x":54,"y":1296},{"x":55,"y":1891},{"x":56,"y":1333},{"x":57,"y":1476},{"x":58,"y":1267},{"x":59,"y":1863},{"x":60,"y":1619},{"x":61,"y":1725},{"x":62,"y":1564},{"x":63,"y":1986},{"x":64,"y":2156},{"x":65,"y":1578},{"x":66,"y":1864},{"x":67,"y":1705},{"x":68,"y":1052},{"x":69,"y":1164},{"x":70,"y":947},{"x":71,"y":858},{"x":72,"y":1035},{"x":73,"y":1063},{"x":74,"y":1006},{"x":75,"y":772},{"x":76,"y":1274},{"x":77,"y":980},{"x":78,"y":1102},{"x":79,"y":808},{"x":80,"y":1458},{"x":81,"y":1110},{"x":82,"y":1133},{"x":83,"y":878},{"x":84,"y":1365},{"x":85,"y":841},{"x":86,"y":808},{"x":87,"y":808},{"x":88,"y":583},{"x":89,"y":1030},{"x":90,"y":677},{"x":91,"y":862},{"x":92,"y":1639},{"x":93,"y":1097},{"x":94,"y":818},{"x":95,"y":1229},{"x":96,"y":657},{"x":97,"y":1136},{"x":98,"y":801},{"x":99,"y":1017},{"x":100,"y":1738},{"x":101,"y":788},{"x":102,"y":1330},{"x":103,"y":897},{"x":104,"y":538},{"x":105,"y":588},{"x":106,"y":572},{"x":107,"y":780},{"x":108,"y":824},{"x":109,"y":574},{"x":110,"y":1287},{"x":111,"y":568},{"x":112,"y":808},{"x":113,"y":656},{"x":114,"y":808},{"x":115,"y":708},{"x":116,"y":556},{"x":117,"y":907},{"x":118,"y":826},{"x":119,"y":821},{"x":120,"y":710},{"x":121,"y":780},{"x":122,"y":490},{"x":123,"y":805},{"x":124,"y":699},{"x":125,"y":542},{"x":126,"y":976},{"x":127,"y":1166},{"x":128,"y":577},{"x":129,"y":994},{"x":130,"y":623},{"x":131,"y":918},{"x":132,"y":1389},{"x":133,"y":1054},{"x":134,"y":912},{"x":135,"y":963},{"x":136,"y":849},{"x":137,"y":646},{"x":138,"y":631},{"x":139,"y":1145},{"x":140,"y":711},{"x":141,"y":666},{"x":142,"y":713},{"x":143,"y":1351},{"x":144,"y":869},{"x":145,"y":654},{"x":146,"y":760},{"x":147,"y":833},{"x":148,"y":345},{"x":149,"y":691},{"x":150,"y":547},{"x":151,"y":928},{"x":152,"y":747},{"x":153,"y":485},{"x":154,"y":942},{"x":155,"y":1112},{"x":156,"y":836},{"x":157,"y":902},{"x":158,"y":611},{"x":159,"y":761},{"x":160,"y":768},{"x":161,"y":777},{"x":162,"y":740},{"x":163,"y":1137},{"x":164,"y":588},{"x":165,"y":481},{"x":166,"y":1065},{"x":167,"y":332},{"x":168,"y":430},{"x":169,"y":619},{"x":170,"y":541},{"x":171,"y":971},{"x":172,"y":756},{"x":173,"y":310},{"x":174,"y":736},{"x":175,"y":562},{"x":176,"y":898},{"x":177,"y":633},{"x":178,"y":773},{"x":179,"y":720},{"x":180,"y":718},{"x":181,"y":367},{"x":182,"y":873},{"x":183,"y":342},{"x":184,"y":798},{"x":185,"y":465},{"x":186,"y":701},{"x":187,"y":721},{"x":188,"y":622},{"x":189,"y":624},{"x":190,"y":435},{"x":191,"y":453},{"x":192,"y":431},{"x":193,"y":605},{"x":194,"y":931},{"x":195,"y":600},{"x":196,"y":440},{"x":197,"y":682},{"x":198,"y":536},{"x":199,"y":620}]},{"name":"35e3af88-7de7-4302-a906-9ce226acf951","values":[{"x":0,"y":19},{"x":1,"y":67},{"x":2,"y":57},{"x":3,"y":53},{"x":4,"y":60},{"x":5,"y":79},{"x":6,"y":55},{"x":7,"y":54},{"x":8,"y":57},{"x":9,"y":91},{"x":10,"y":329},{"x":11,"y":250},{"x":12,"y":251},{"x":13,"y":316},{"x":14,"y":328},{"x":15,"y":217},{"x":16,"y":227},{"x":17,"y":248},{"x":18,"y":256},{"x":19,"y":193},{"x":20,"y":178},{"x":21,"y":242},{"x":22,"y":418},{"x":23,"y":227},{"x":24,"y":334},{"x":25,"y":350},{"x":26,"y":381},{"x":27,"y":218},{"x":28,"y":259},{"x":29,"y":286},{"x":30,"y":252},{"x":31,"y":228},{"x":32,"y":325},{"x":33,"y":170},{"x":34,"y":439},{"x":35,"y":300},{"x":36,"y":320},{"x":37,"y":517},{"x":38,"y":283},{"x":39,"y":358},{"x":40,"y":372},{"x":41,"y":257},{"x":42,"y":290},{"x":43,"y":258},{"x":44,"y":242},{"x":45,"y":203},{"x":46,"y":276},{"x":47,"y":370},{"x":48,"y":226},{"x":49,"y":337},{"x":50,"y":284},{"x":51,"y":343},{"x":52,"y":313},{"x":53,"y":283},{"x":54,"y":232},{"x":55,"y":291},{"x":56,"y":268},{"x":57,"y":276},{"x":58,"y":285},{"x":59,"y":186},{"x":60,"y":182},{"x":61,"y":257},{"x":62,"y":264},{"x":63,"y":269},{"x":64,"y":388},{"x":65,"y":231},{"x":66,"y":191},{"x":67,"y":253},{"x":68,"y":146},{"x":69,"y":282},{"x":70,"y":248},{"x":71,"y":308},{"x":72,"y":302},{"x":73,"y":302},{"x":74,"y":388},{"x":75,"y":297},{"x":76,"y":279},{"x":77,"y":254},{"x":78,"y":227},{"x":79,"y":215},{"x":80,"y":293},{"x":81,"y":250},{"x":82,"y":300},{"x":83,"y":360},{"x":84,"y":290},{"x":85,"y":261},{"x":86,"y":347},{"x":87,"y":290},{"x":88,"y":272},{"x":89,"y":341},{"x":90,"y":264},{"x":91,"y":295},{"x":92,"y":786},{"x":93,"y":419},{"x":94,"y":204},{"x":95,"y":266},{"x":96,"y":355},{"x":97,"y":250},{"x":98,"y":207},{"x":99,"y":239},{"x":100,"y":423},{"x":101,"y":295},{"x":102,"y":328},{"x":103,"y":269},{"x":104,"y":261},{"x":105,"y":287},{"x":106,"y":342},{"x":107,"y":283},{"x":108,"y":379},{"x":109,"y":292},{"x":110,"y":227},{"x":111,"y":350},{"x":112,"y":224},{"x":113,"y":250},{"x":114,"y":370},{"x":115,"y":464},{"x":116,"y":293},{"x":117,"y":475},{"x":118,"y":467},{"x":119,"y":450},{"x":120,"y":395},{"x":121,"y":459},{"x":122,"y":308},{"x":123,"y":535},{"x":124,"y":349},{"x":125,"y":351},{"x":126,"y":486},{"x":127,"y":259},{"x":128,"y":291},{"x":129,"y":373},{"x":130,"y":484},{"x":131,"y":479},{"x":132,"y":384},{"x":133,"y":320},{"x":134,"y":412},{"x":135,"y":310},{"x":136,"y":399},{"x":137,"y":353},{"x":138,"y":393},{"x":139,"y":341},{"x":140,"y":327},{"x":141,"y":241},{"x":142,"y":385},{"x":143,"y":353},{"x":144,"y":445},{"x":145,"y":375},{"x":146,"y":380},{"x":147,"y":365},{"x":148,"y":289},{"x":149,"y":456},{"x":150,"y":269},{"x":151,"y":325},{"x":152,"y":296},{"x":153,"y":294},{"x":154,"y":310},{"x":155,"y":373},{"x":156,"y":331},{"x":157,"y":275},{"x":158,"y":267},{"x":159,"y":264},{"x":160,"y":333},{"x":161,"y":391},{"x":162,"y":306},{"x":163,"y":206},{"x":164,"y":216},{"x":165,"y":274},{"x":166,"y":242},{"x":167,"y":234},{"x":168,"y":252},{"x":169,"y":298},{"x":170,"y":314},{"x":171,"y":389},{"x":172,"y":263},{"x":173,"y":279},{"x":174,"y":311},{"x":175,"y":421},{"x":176,"y":352},{"x":177,"y":317},{"x":178,"y":389},{"x":179,"y":559},{"x":180,"y":321},{"x":181,"y":292},{"x":182,"y":348},{"x":183,"y":290},{"x":184,"y":310},{"x":185,"y":324},{"x":186,"y":371},{"x":187,"y":314},{"x":188,"y":267},{"x":189,"y":268},{"x":190,"y":308},{"x":191,"y":266},{"x":192,"y":278},{"x":193,"y":306},{"x":194,"y":324},{"x":195,"y":307},{"x":196,"y":325},{"x":197,"y":413},{"x":198,"y":260},{"x":199,"y":327}]},{"name":"27df3a37-dc21-4846-8456-df230142d31c","values":[{"x":0,"y":4},{"x":1,"y":3},{"x":2,"y":3},{"x":3,"y":2},{"x":4,"y":3},{"x":5,"y":3},{"x":6,"y":2},{"x":7,"y":3},{"x":8,"y":2},{"x":9,"y":2},{"x":10,"y":1},{"x":11,"y":1},{"x":12,"y":1},{"x":13,"y":0},{"x":14,"y":1},{"x":15,"y":1},{"x":16,"y":1},{"x":17,"y":1},{"x":18,"y":2},{"x":19,"y":1},{"x":20,"y":2},{"x":21,"y":1},{"x":22,"y":5},{"x":23,"y":1},{"x":24,"y":2},{"x":25,"y":1},{"x":26,"y":1},{"x":27,"y":1},{"x":28,"y":1},{"x":29,"y":1},{"x":30,"y":1},{"x":31,"y":1},{"x":32,"y":1},{"x":33,"y":2},{"x":34,"y":1},{"x":35,"y":2},{"x":36,"y":1},{"x":37,"y":1},{"x":38,"y":1},{"x":39,"y":1},{"x":40,"y":1},{"x":41,"y":1},{"x":42,"y":1},{"x":43,"y":1},{"x":44,"y":1},{"x":45,"y":1},{"x":46,"y":1},{"x":47,"y":2},{"x":48,"y":2},{"x":49,"y":2},{"x":50,"y":1},{"x":51,"y":2},{"x":52,"y":2},{"x":53,"y":2},{"x":54,"y":1},{"x":55,"y":1},{"x":56,"y":1},{"x":57,"y":1},{"x":58,"y":1},{"x":59,"y":1},{"x":60,"y":1},{"x":61,"y":1},{"x":62,"y":1},{"x":63,"y":1},{"x":64,"y":1},{"x":65,"y":1},{"x":66,"y":2},{"x":67,"y":1},{"x":68,"y":1},{"x":69,"y":1},{"x":70,"y":1},{"x":71,"y":2},{"x":72,"y":1},{"x":73,"y":1},{"x":74,"y":1},{"x":75,"y":1},{"x":76,"y":1},{"x":77,"y":1},{"x":78,"y":1},{"x":79,"y":1},{"x":80,"y":2},{"x":81,"y":1},{"x":82,"y":1},{"x":83,"y":1},{"x":84,"y":1},{"x":85,"y":1},{"x":86,"y":1},{"x":87,"y":1},{"x":88,"y":1},{"x":89,"y":1},{"x":90,"y":1},{"x":91,"y":1},{"x":92,"y":3},{"x":93,"y":3},{"x":94,"y":0},{"x":95,"y":1},{"x":96,"y":1},{"x":97,"y":1},{"x":98,"y":1},{"x":99,"y":3},{"x":100,"y":1},{"x":101,"y":1},{"x":102,"y":1},{"x":103,"y":1},{"x":104,"y":1},{"x":105,"y":1},{"x":106,"y":1},{"x":107,"y":1},{"x":108,"y":1},{"x":109,"y":1},{"x":110,"y":1},{"x":111,"y":1},{"x":112,"y":2},{"x":113,"y":1},{"x":114,"y":0},{"x":115,"y":1},{"x":116,"y":1},{"x":117,"y":1},{"x":118,"y":1},{"x":119,"y":1},{"x":120,"y":1},{"x":121,"y":2},{"x":122,"y":1},{"x":123,"y":0},{"x":124,"y":1},{"x":125,"y":1},{"x":126,"y":1},{"x":127,"y":1},{"x":128,"y":1},{"x":129,"y":1},{"x":130,"y":1},{"x":131,"y":0},{"x":132,"y":2},{"x":133,"y":1},{"x":134,"y":1},{"x":135,"y":1},{"x":136,"y":2},{"x":137,"y":0},{"x":138,"y":0},{"x":139,"y":0},{"x":140,"y":2},{"x":141,"y":1},{"x":142,"y":1},{"x":143,"y":1},{"x":144,"y":2},{"x":145,"y":2},{"x":146,"y":1},{"x":147,"y":1},{"x":148,"y":1},{"x":149,"y":1},{"x":150,"y":1},{"x":151,"y":1},{"x":152,"y":1},{"x":153,"y":1},{"x":154,"y":2},{"x":155,"y":1},{"x":156,"y":1},{"x":157,"y":1},{"x":158,"y":1},{"x":159,"y":1},{"x":160,"y":1},{"x":161,"y":1},{"x":162,"y":1},{"x":163,"y":1},{"x":164,"y":1},{"x":165,"y":1},{"x":166,"y":1},{"x":167,"y":1},{"x":168,"y":1},{"x":169,"y":1},{"x":170,"y":1},{"x":171,"y":0},{"x":172,"y":1},{"x":173,"y":1},{"x":174,"y":1},{"x":175,"y":2},{"x":176,"y":2},{"x":177,"y":1},{"x":178,"y":1},{"x":179,"y":1},{"x":180,"y":1},{"x":181,"y":1},{"x":182,"y":4},{"x":183,"y":1},{"x":184,"y":1},{"x":185,"y":0},{"x":186,"y":1},{"x":187,"y":1},{"x":188,"y":1},{"x":189,"y":1},{"x":190,"y":1},{"x":191,"y":0},{"x":192,"y":1},{"x":193,"y":1},{"x":194,"y":1},{"x":195,"y":1},{"x":196,"y":1},{"x":197,"y":1},{"x":198,"y":1},{"x":199,"y":1}]},{"name":"6716777e-de4d-409e-acbb-7722a5bcb202","values":[{"x":0,"y":3249},{"x":1,"y":4199},{"x":2,"y":2978},{"x":3,"y":1807},{"x":4,"y":1734},{"x":5,"y":2369},{"x":6,"y":1312},{"x":7,"y":795},{"x":8,"y":1427},{"x":9,"y":733},{"x":10,"y":1293},{"x":11,"y":1359},{"x":12,"y":309},{"x":13,"y":1067},{"x":14,"y":1903},{"x":15,"y":912},{"x":16,"y":1420},{"x":17,"y":937},{"x":18,"y":2181},{"x":19,"y":1469},{"x":20,"y":547},{"x":21,"y":843},{"x":22,"y":903},{"x":23,"y":836},{"x":24,"y":755},{"x":25,"y":590},{"x":26,"y":678},{"x":27,"y":285},{"x":28,"y":177},{"x":29,"y":673},{"x":30,"y":724},{"x":31,"y":779},{"x":32,"y":632},{"x":33,"y":372},{"x":34,"y":668},{"x":35,"y":801},{"x":36,"y":768},{"x":37,"y":619},{"x":38,"y":361},{"x":39,"y":573},{"x":40,"y":975},{"x":41,"y":569},{"x":42,"y":1152},{"x":43,"y":1134},{"x":44,"y":870},{"x":45,"y":385},{"x":46,"y":195},{"x":47,"y":815},{"x":48,"y":950},{"x":49,"y":442},{"x":50,"y":823},{"x":51,"y":1084},{"x":52,"y":715},{"x":53,"y":845},{"x":54,"y":1063},{"x":55,"y":1599},{"x":56,"y":1064},{"x":57,"y":1199},{"x":58,"y":981},{"x":59,"y":1676},{"x":60,"y":1436},{"x":61,"y":1467},{"x":62,"y":1299},{"x":63,"y":1716},{"x":64,"y":1767},{"x":65,"y":1346},{"x":66,"y":1671},{"x":67,"y":1451},{"x":68,"y":905},{"x":69,"y":881},{"x":70,"y":698},{"x":71,"y":548},{"x":72,"y":732},{"x":73,"y":760},{"x":74,"y":617},{"x":75,"y":474},{"x":76,"y":994},{"x":77,"y":725},{"x":78,"y":874},{"x":79,"y":592},{"x":80,"y":1163},{"x":81,"y":859},{"x":82,"y":832},{"x":83,"y":517},{"x":84,"y":1074},{"x":85,"y":579},{"x":86,"y":460},{"x":87,"y":517},{"x":88,"y":310},{"x":89,"y":688},{"x":90,"y":412},{"x":91,"y":566},{"x":92,"y":850},{"x":93,"y":675},{"x":94,"y":614},{"x":95,"y":962},{"x":96,"y":301},{"x":97,"y":885},{"x":98,"y":593},{"x":99,"y":775},{"x":100,"y":1314},{"x":101,"y":492},{"x":102,"y":1001},{"x":103,"y":627},{"x":104,"y":276},{"x":105,"y":300},{"x":106,"y":229},{"x":107,"y":496},{"x":108,"y":444},{"x":109,"y":281},{"x":110,"y":1059},{"x":111,"y":217},{"x":112,"y":582},{"x":113,"y":405},{"x":114,"y":438},{"x":115,"y":243},{"x":116,"y":262},{"x":117,"y":431},{"x":118,"y":358},{"x":119,"y":370},{"x":120,"y":314},{"x":121,"y":319},{"x":122,"y":181},{"x":123,"y":270},{"x":124,"y":349},{"x":125,"y":190},{"x":126,"y":489},{"x":127,"y":906},{"x":128,"y":285},{"x":129,"y":620},{"x":130,"y":138},{"x":131,"y":439},{"x":132,"y":1003},{"x":133,"y":733},{"x":134,"y":499},{"x":135,"y":652},{"x":136,"y":448},{"x":137,"y":293},{"x":138,"y":238},{"x":139,"y":804},{"x":140,"y":382},{"x":141,"y":424},{"x":142,"y":327},{"x":143,"y":997},{"x":144,"y":422},{"x":145,"y":277},{"x":146,"y":379},{"x":147,"y":467},{"x":148,"y":55},{"x":149,"y":234},{"x":150,"y":277},{"x":151,"y":602},{"x":152,"y":450},{"x":153,"y":190},{"x":154,"y":630},{"x":155,"y":738},{"x":156,"y":504},{"x":157,"y":626},{"x":158,"y":343},{"x":159,"y":496},{"x":160,"y":434},{"x":161,"y":385},{"x":162,"y":433},{"x":163,"y":930},{"x":164,"y":371},{"x":165,"y":206},{"x":166,"y":822},{"x":167,"y":97},{"x":168,"y":177},{"x":169,"y":320},{"x":170,"y":226},{"x":171,"y":582},{"x":172,"y":492},{"x":173,"y":30},{"x":174,"y":424},{"x":175,"y":139},{"x":176,"y":544},{"x":177,"y":315},{"x":178,"y":383},{"x":179,"y":160},{"x":180,"y":396},{"x":181,"y":74},{"x":182,"y":521},{"x":183,"y":51},{"x":184,"y":487},{"x":185,"y":141},{"x":186,"y":329},{"x":187,"y":406},{"x":188,"y":354},{"x":189,"y":355},{"x":190,"y":126},{"x":191,"y":187},{"x":192,"y":152},{"x":193,"y":298},{"x":194,"y":606},{"x":195,"y":292},{"x":196,"y":114},{"x":197,"y":268},{"x":198,"y":275},{"x":199,"y":292}]}],"marks":[{"type":"line","from":{"data":"90342872-f717-4977-b967-3404cba3f04a"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"yellow"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"35e3af88-7de7-4302-a906-9ce226acf951"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"red"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"27df3a37-dc21-4846-8456-df230142d31c"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"green"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"6716777e-de4d-409e-acbb-7722a5bcb202"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"blue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"90342872-f717-4977-b967-3404cba3f04a\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain [0 4269]}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"90342872-f717-4977-b967-3404cba3f04a\", :values ({:x 0, :y 3272} {:x 1, :y 4269} {:x 2, :y 3038} {:x 3, :y 1862} {:x 4, :y 1797} {:x 5, :y 2451} {:x 6, :y 1369} {:x 7, :y 852} {:x 8, :y 1486} {:x 9, :y 826} {:x 10, :y 1623} {:x 11, :y 1610} {:x 12, :y 561} {:x 13, :y 1383} {:x 14, :y 2232} {:x 15, :y 1130} {:x 16, :y 1648} {:x 17, :y 1186} {:x 18, :y 2439} {:x 19, :y 1663} {:x 20, :y 727} {:x 21, :y 1086} {:x 22, :y 1326} {:x 23, :y 1064} {:x 24, :y 1091} {:x 25, :y 941} {:x 26, :y 1060} {:x 27, :y 504} {:x 28, :y 437} {:x 29, :y 960} {:x 30, :y 977} {:x 31, :y 1008} {:x 32, :y 958} {:x 33, :y 544} {:x 34, :y 1108} {:x 35, :y 1103} {:x 36, :y 1089} {:x 37, :y 1137} {:x 38, :y 645} {:x 39, :y 932} {:x 40, :y 1348} {:x 41, :y 827} {:x 42, :y 1443} {:x 43, :y 1393} {:x 44, :y 1113} {:x 45, :y 589} {:x 46, :y 472} {:x 47, :y 1187} {:x 48, :y 1178} {:x 49, :y 781} {:x 50, :y 1108} {:x 51, :y 1429} {:x 52, :y 1030} {:x 53, :y 1130} {:x 54, :y 1296} {:x 55, :y 1891} {:x 56, :y 1333} {:x 57, :y 1476} {:x 58, :y 1267} {:x 59, :y 1863} {:x 60, :y 1619} {:x 61, :y 1725} {:x 62, :y 1564} {:x 63, :y 1986} {:x 64, :y 2156} {:x 65, :y 1578} {:x 66, :y 1864} {:x 67, :y 1705} {:x 68, :y 1052} {:x 69, :y 1164} {:x 70, :y 947} {:x 71, :y 858} {:x 72, :y 1035} {:x 73, :y 1063} {:x 74, :y 1006} {:x 75, :y 772} {:x 76, :y 1274} {:x 77, :y 980} {:x 78, :y 1102} {:x 79, :y 808} {:x 80, :y 1458} {:x 81, :y 1110} {:x 82, :y 1133} {:x 83, :y 878} {:x 84, :y 1365} {:x 85, :y 841} {:x 86, :y 808} {:x 87, :y 808} {:x 88, :y 583} {:x 89, :y 1030} {:x 90, :y 677} {:x 91, :y 862} {:x 92, :y 1639} {:x 93, :y 1097} {:x 94, :y 818} {:x 95, :y 1229} {:x 96, :y 657} {:x 97, :y 1136} {:x 98, :y 801} {:x 99, :y 1017} {:x 100, :y 1738} {:x 101, :y 788} {:x 102, :y 1330} {:x 103, :y 897} {:x 104, :y 538} {:x 105, :y 588} {:x 106, :y 572} {:x 107, :y 780} {:x 108, :y 824} {:x 109, :y 574} {:x 110, :y 1287} {:x 111, :y 568} {:x 112, :y 808} {:x 113, :y 656} {:x 114, :y 808} {:x 115, :y 708} {:x 116, :y 556} {:x 117, :y 907} {:x 118, :y 826} {:x 119, :y 821} {:x 120, :y 710} {:x 121, :y 780} {:x 122, :y 490} {:x 123, :y 805} {:x 124, :y 699} {:x 125, :y 542} {:x 126, :y 976} {:x 127, :y 1166} {:x 128, :y 577} {:x 129, :y 994} {:x 130, :y 623} {:x 131, :y 918} {:x 132, :y 1389} {:x 133, :y 1054} {:x 134, :y 912} {:x 135, :y 963} {:x 136, :y 849} {:x 137, :y 646} {:x 138, :y 631} {:x 139, :y 1145} {:x 140, :y 711} {:x 141, :y 666} {:x 142, :y 713} {:x 143, :y 1351} {:x 144, :y 869} {:x 145, :y 654} {:x 146, :y 760} {:x 147, :y 833} {:x 148, :y 345} {:x 149, :y 691} {:x 150, :y 547} {:x 151, :y 928} {:x 152, :y 747} {:x 153, :y 485} {:x 154, :y 942} {:x 155, :y 1112} {:x 156, :y 836} {:x 157, :y 902} {:x 158, :y 611} {:x 159, :y 761} {:x 160, :y 768} {:x 161, :y 777} {:x 162, :y 740} {:x 163, :y 1137} {:x 164, :y 588} {:x 165, :y 481} {:x 166, :y 1065} {:x 167, :y 332} {:x 168, :y 430} {:x 169, :y 619} {:x 170, :y 541} {:x 171, :y 971} {:x 172, :y 756} {:x 173, :y 310} {:x 174, :y 736} {:x 175, :y 562} {:x 176, :y 898} {:x 177, :y 633} {:x 178, :y 773} {:x 179, :y 720} {:x 180, :y 718} {:x 181, :y 367} {:x 182, :y 873} {:x 183, :y 342} {:x 184, :y 798} {:x 185, :y 465} {:x 186, :y 701} {:x 187, :y 721} {:x 188, :y 622} {:x 189, :y 624} {:x 190, :y 435} {:x 191, :y 453} {:x 192, :y 431} {:x 193, :y 605} {:x 194, :y 931} {:x 195, :y 600} {:x 196, :y 440} {:x 197, :y 682} {:x 198, :y 536} {:x 199, :y 620})} {:name \"35e3af88-7de7-4302-a906-9ce226acf951\", :values ({:x 0, :y 19} {:x 1, :y 67} {:x 2, :y 57} {:x 3, :y 53} {:x 4, :y 60} {:x 5, :y 79} {:x 6, :y 55} {:x 7, :y 54} {:x 8, :y 57} {:x 9, :y 91} {:x 10, :y 329} {:x 11, :y 250} {:x 12, :y 251} {:x 13, :y 316} {:x 14, :y 328} {:x 15, :y 217} {:x 16, :y 227} {:x 17, :y 248} {:x 18, :y 256} {:x 19, :y 193} {:x 20, :y 178} {:x 21, :y 242} {:x 22, :y 418} {:x 23, :y 227} {:x 24, :y 334} {:x 25, :y 350} {:x 26, :y 381} {:x 27, :y 218} {:x 28, :y 259} {:x 29, :y 286} {:x 30, :y 252} {:x 31, :y 228} {:x 32, :y 325} {:x 33, :y 170} {:x 34, :y 439} {:x 35, :y 300} {:x 36, :y 320} {:x 37, :y 517} {:x 38, :y 283} {:x 39, :y 358} {:x 40, :y 372} {:x 41, :y 257} {:x 42, :y 290} {:x 43, :y 258} {:x 44, :y 242} {:x 45, :y 203} {:x 46, :y 276} {:x 47, :y 370} {:x 48, :y 226} {:x 49, :y 337} {:x 50, :y 284} {:x 51, :y 343} {:x 52, :y 313} {:x 53, :y 283} {:x 54, :y 232} {:x 55, :y 291} {:x 56, :y 268} {:x 57, :y 276} {:x 58, :y 285} {:x 59, :y 186} {:x 60, :y 182} {:x 61, :y 257} {:x 62, :y 264} {:x 63, :y 269} {:x 64, :y 388} {:x 65, :y 231} {:x 66, :y 191} {:x 67, :y 253} {:x 68, :y 146} {:x 69, :y 282} {:x 70, :y 248} {:x 71, :y 308} {:x 72, :y 302} {:x 73, :y 302} {:x 74, :y 388} {:x 75, :y 297} {:x 76, :y 279} {:x 77, :y 254} {:x 78, :y 227} {:x 79, :y 215} {:x 80, :y 293} {:x 81, :y 250} {:x 82, :y 300} {:x 83, :y 360} {:x 84, :y 290} {:x 85, :y 261} {:x 86, :y 347} {:x 87, :y 290} {:x 88, :y 272} {:x 89, :y 341} {:x 90, :y 264} {:x 91, :y 295} {:x 92, :y 786} {:x 93, :y 419} {:x 94, :y 204} {:x 95, :y 266} {:x 96, :y 355} {:x 97, :y 250} {:x 98, :y 207} {:x 99, :y 239} {:x 100, :y 423} {:x 101, :y 295} {:x 102, :y 328} {:x 103, :y 269} {:x 104, :y 261} {:x 105, :y 287} {:x 106, :y 342} {:x 107, :y 283} {:x 108, :y 379} {:x 109, :y 292} {:x 110, :y 227} {:x 111, :y 350} {:x 112, :y 224} {:x 113, :y 250} {:x 114, :y 370} {:x 115, :y 464} {:x 116, :y 293} {:x 117, :y 475} {:x 118, :y 467} {:x 119, :y 450} {:x 120, :y 395} {:x 121, :y 459} {:x 122, :y 308} {:x 123, :y 535} {:x 124, :y 349} {:x 125, :y 351} {:x 126, :y 486} {:x 127, :y 259} {:x 128, :y 291} {:x 129, :y 373} {:x 130, :y 484} {:x 131, :y 479} {:x 132, :y 384} {:x 133, :y 320} {:x 134, :y 412} {:x 135, :y 310} {:x 136, :y 399} {:x 137, :y 353} {:x 138, :y 393} {:x 139, :y 341} {:x 140, :y 327} {:x 141, :y 241} {:x 142, :y 385} {:x 143, :y 353} {:x 144, :y 445} {:x 145, :y 375} {:x 146, :y 380} {:x 147, :y 365} {:x 148, :y 289} {:x 149, :y 456} {:x 150, :y 269} {:x 151, :y 325} {:x 152, :y 296} {:x 153, :y 294} {:x 154, :y 310} {:x 155, :y 373} {:x 156, :y 331} {:x 157, :y 275} {:x 158, :y 267} {:x 159, :y 264} {:x 160, :y 333} {:x 161, :y 391} {:x 162, :y 306} {:x 163, :y 206} {:x 164, :y 216} {:x 165, :y 274} {:x 166, :y 242} {:x 167, :y 234} {:x 168, :y 252} {:x 169, :y 298} {:x 170, :y 314} {:x 171, :y 389} {:x 172, :y 263} {:x 173, :y 279} {:x 174, :y 311} {:x 175, :y 421} {:x 176, :y 352} {:x 177, :y 317} {:x 178, :y 389} {:x 179, :y 559} {:x 180, :y 321} {:x 181, :y 292} {:x 182, :y 348} {:x 183, :y 290} {:x 184, :y 310} {:x 185, :y 324} {:x 186, :y 371} {:x 187, :y 314} {:x 188, :y 267} {:x 189, :y 268} {:x 190, :y 308} {:x 191, :y 266} {:x 192, :y 278} {:x 193, :y 306} {:x 194, :y 324} {:x 195, :y 307} {:x 196, :y 325} {:x 197, :y 413} {:x 198, :y 260} {:x 199, :y 327})} {:name \"27df3a37-dc21-4846-8456-df230142d31c\", :values ({:x 0, :y 4} {:x 1, :y 3} {:x 2, :y 3} {:x 3, :y 2} {:x 4, :y 3} {:x 5, :y 3} {:x 6, :y 2} {:x 7, :y 3} {:x 8, :y 2} {:x 9, :y 2} {:x 10, :y 1} {:x 11, :y 1} {:x 12, :y 1} {:x 13, :y 0} {:x 14, :y 1} {:x 15, :y 1} {:x 16, :y 1} {:x 17, :y 1} {:x 18, :y 2} {:x 19, :y 1} {:x 20, :y 2} {:x 21, :y 1} {:x 22, :y 5} {:x 23, :y 1} {:x 24, :y 2} {:x 25, :y 1} {:x 26, :y 1} {:x 27, :y 1} {:x 28, :y 1} {:x 29, :y 1} {:x 30, :y 1} {:x 31, :y 1} {:x 32, :y 1} {:x 33, :y 2} {:x 34, :y 1} {:x 35, :y 2} {:x 36, :y 1} {:x 37, :y 1} {:x 38, :y 1} {:x 39, :y 1} {:x 40, :y 1} {:x 41, :y 1} {:x 42, :y 1} {:x 43, :y 1} {:x 44, :y 1} {:x 45, :y 1} {:x 46, :y 1} {:x 47, :y 2} {:x 48, :y 2} {:x 49, :y 2} {:x 50, :y 1} {:x 51, :y 2} {:x 52, :y 2} {:x 53, :y 2} {:x 54, :y 1} {:x 55, :y 1} {:x 56, :y 1} {:x 57, :y 1} {:x 58, :y 1} {:x 59, :y 1} {:x 60, :y 1} {:x 61, :y 1} {:x 62, :y 1} {:x 63, :y 1} {:x 64, :y 1} {:x 65, :y 1} {:x 66, :y 2} {:x 67, :y 1} {:x 68, :y 1} {:x 69, :y 1} {:x 70, :y 1} {:x 71, :y 2} {:x 72, :y 1} {:x 73, :y 1} {:x 74, :y 1} {:x 75, :y 1} {:x 76, :y 1} {:x 77, :y 1} {:x 78, :y 1} {:x 79, :y 1} {:x 80, :y 2} {:x 81, :y 1} {:x 82, :y 1} {:x 83, :y 1} {:x 84, :y 1} {:x 85, :y 1} {:x 86, :y 1} {:x 87, :y 1} {:x 88, :y 1} {:x 89, :y 1} {:x 90, :y 1} {:x 91, :y 1} {:x 92, :y 3} {:x 93, :y 3} {:x 94, :y 0} {:x 95, :y 1} {:x 96, :y 1} {:x 97, :y 1} {:x 98, :y 1} {:x 99, :y 3} {:x 100, :y 1} {:x 101, :y 1} {:x 102, :y 1} {:x 103, :y 1} {:x 104, :y 1} {:x 105, :y 1} {:x 106, :y 1} {:x 107, :y 1} {:x 108, :y 1} {:x 109, :y 1} {:x 110, :y 1} {:x 111, :y 1} {:x 112, :y 2} {:x 113, :y 1} {:x 114, :y 0} {:x 115, :y 1} {:x 116, :y 1} {:x 117, :y 1} {:x 118, :y 1} {:x 119, :y 1} {:x 120, :y 1} {:x 121, :y 2} {:x 122, :y 1} {:x 123, :y 0} {:x 124, :y 1} {:x 125, :y 1} {:x 126, :y 1} {:x 127, :y 1} {:x 128, :y 1} {:x 129, :y 1} {:x 130, :y 1} {:x 131, :y 0} {:x 132, :y 2} {:x 133, :y 1} {:x 134, :y 1} {:x 135, :y 1} {:x 136, :y 2} {:x 137, :y 0} {:x 138, :y 0} {:x 139, :y 0} {:x 140, :y 2} {:x 141, :y 1} {:x 142, :y 1} {:x 143, :y 1} {:x 144, :y 2} {:x 145, :y 2} {:x 146, :y 1} {:x 147, :y 1} {:x 148, :y 1} {:x 149, :y 1} {:x 150, :y 1} {:x 151, :y 1} {:x 152, :y 1} {:x 153, :y 1} {:x 154, :y 2} {:x 155, :y 1} {:x 156, :y 1} {:x 157, :y 1} {:x 158, :y 1} {:x 159, :y 1} {:x 160, :y 1} {:x 161, :y 1} {:x 162, :y 1} {:x 163, :y 1} {:x 164, :y 1} {:x 165, :y 1} {:x 166, :y 1} {:x 167, :y 1} {:x 168, :y 1} {:x 169, :y 1} {:x 170, :y 1} {:x 171, :y 0} {:x 172, :y 1} {:x 173, :y 1} {:x 174, :y 1} {:x 175, :y 2} {:x 176, :y 2} {:x 177, :y 1} {:x 178, :y 1} {:x 179, :y 1} {:x 180, :y 1} {:x 181, :y 1} {:x 182, :y 4} {:x 183, :y 1} {:x 184, :y 1} {:x 185, :y 0} {:x 186, :y 1} {:x 187, :y 1} {:x 188, :y 1} {:x 189, :y 1} {:x 190, :y 1} {:x 191, :y 0} {:x 192, :y 1} {:x 193, :y 1} {:x 194, :y 1} {:x 195, :y 1} {:x 196, :y 1} {:x 197, :y 1} {:x 198, :y 1} {:x 199, :y 1})} {:name \"6716777e-de4d-409e-acbb-7722a5bcb202\", :values ({:x 0, :y 3249} {:x 1, :y 4199} {:x 2, :y 2978} {:x 3, :y 1807} {:x 4, :y 1734} {:x 5, :y 2369} {:x 6, :y 1312} {:x 7, :y 795} {:x 8, :y 1427} {:x 9, :y 733} {:x 10, :y 1293} {:x 11, :y 1359} {:x 12, :y 309} {:x 13, :y 1067} {:x 14, :y 1903} {:x 15, :y 912} {:x 16, :y 1420} {:x 17, :y 937} {:x 18, :y 2181} {:x 19, :y 1469} {:x 20, :y 547} {:x 21, :y 843} {:x 22, :y 903} {:x 23, :y 836} {:x 24, :y 755} {:x 25, :y 590} {:x 26, :y 678} {:x 27, :y 285} {:x 28, :y 177} {:x 29, :y 673} {:x 30, :y 724} {:x 31, :y 779} {:x 32, :y 632} {:x 33, :y 372} {:x 34, :y 668} {:x 35, :y 801} {:x 36, :y 768} {:x 37, :y 619} {:x 38, :y 361} {:x 39, :y 573} {:x 40, :y 975} {:x 41, :y 569} {:x 42, :y 1152} {:x 43, :y 1134} {:x 44, :y 870} {:x 45, :y 385} {:x 46, :y 195} {:x 47, :y 815} {:x 48, :y 950} {:x 49, :y 442} {:x 50, :y 823} {:x 51, :y 1084} {:x 52, :y 715} {:x 53, :y 845} {:x 54, :y 1063} {:x 55, :y 1599} {:x 56, :y 1064} {:x 57, :y 1199} {:x 58, :y 981} {:x 59, :y 1676} {:x 60, :y 1436} {:x 61, :y 1467} {:x 62, :y 1299} {:x 63, :y 1716} {:x 64, :y 1767} {:x 65, :y 1346} {:x 66, :y 1671} {:x 67, :y 1451} {:x 68, :y 905} {:x 69, :y 881} {:x 70, :y 698} {:x 71, :y 548} {:x 72, :y 732} {:x 73, :y 760} {:x 74, :y 617} {:x 75, :y 474} {:x 76, :y 994} {:x 77, :y 725} {:x 78, :y 874} {:x 79, :y 592} {:x 80, :y 1163} {:x 81, :y 859} {:x 82, :y 832} {:x 83, :y 517} {:x 84, :y 1074} {:x 85, :y 579} {:x 86, :y 460} {:x 87, :y 517} {:x 88, :y 310} {:x 89, :y 688} {:x 90, :y 412} {:x 91, :y 566} {:x 92, :y 850} {:x 93, :y 675} {:x 94, :y 614} {:x 95, :y 962} {:x 96, :y 301} {:x 97, :y 885} {:x 98, :y 593} {:x 99, :y 775} {:x 100, :y 1314} {:x 101, :y 492} {:x 102, :y 1001} {:x 103, :y 627} {:x 104, :y 276} {:x 105, :y 300} {:x 106, :y 229} {:x 107, :y 496} {:x 108, :y 444} {:x 109, :y 281} {:x 110, :y 1059} {:x 111, :y 217} {:x 112, :y 582} {:x 113, :y 405} {:x 114, :y 438} {:x 115, :y 243} {:x 116, :y 262} {:x 117, :y 431} {:x 118, :y 358} {:x 119, :y 370} {:x 120, :y 314} {:x 121, :y 319} {:x 122, :y 181} {:x 123, :y 270} {:x 124, :y 349} {:x 125, :y 190} {:x 126, :y 489} {:x 127, :y 906} {:x 128, :y 285} {:x 129, :y 620} {:x 130, :y 138} {:x 131, :y 439} {:x 132, :y 1003} {:x 133, :y 733} {:x 134, :y 499} {:x 135, :y 652} {:x 136, :y 448} {:x 137, :y 293} {:x 138, :y 238} {:x 139, :y 804} {:x 140, :y 382} {:x 141, :y 424} {:x 142, :y 327} {:x 143, :y 997} {:x 144, :y 422} {:x 145, :y 277} {:x 146, :y 379} {:x 147, :y 467} {:x 148, :y 55} {:x 149, :y 234} {:x 150, :y 277} {:x 151, :y 602} {:x 152, :y 450} {:x 153, :y 190} {:x 154, :y 630} {:x 155, :y 738} {:x 156, :y 504} {:x 157, :y 626} {:x 158, :y 343} {:x 159, :y 496} {:x 160, :y 434} {:x 161, :y 385} {:x 162, :y 433} {:x 163, :y 930} {:x 164, :y 371} {:x 165, :y 206} {:x 166, :y 822} {:x 167, :y 97} {:x 168, :y 177} {:x 169, :y 320} {:x 170, :y 226} {:x 171, :y 582} {:x 172, :y 492} {:x 173, :y 30} {:x 174, :y 424} {:x 175, :y 139} {:x 176, :y 544} {:x 177, :y 315} {:x 178, :y 383} {:x 179, :y 160} {:x 180, :y 396} {:x 181, :y 74} {:x 182, :y 521} {:x 183, :y 51} {:x 184, :y 487} {:x 185, :y 141} {:x 186, :y 329} {:x 187, :y 406} {:x 188, :y 354} {:x 189, :y 355} {:x 190, :y 126} {:x 191, :y 187} {:x 192, :y 152} {:x 193, :y 298} {:x 194, :y 606} {:x 195, :y 292} {:x 196, :y 114} {:x 197, :y 268} {:x 198, :y 275} {:x 199, :y 292})}), :marks ({:type \"line\", :from {:data \"90342872-f717-4977-b967-3404cba3f04a\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"yellow\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"35e3af88-7de7-4302-a906-9ce226acf951\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"red\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"27df3a37-dc21-4846-8456-df230142d31c\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"green\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"6716777e-de4d-409e-acbb-7722a5bcb202\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"blue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}})}}"}
;; <=

;; **
;;; Spread of complexity in population
;; **

;; @@
(plot/histogram (map :complexity (:rabble result)))
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"c99cb365-6dc0-4884-8b15-335f8817c202","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"c99cb365-6dc0-4884-8b15-335f8817c202","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"c99cb365-6dc0-4884-8b15-335f8817c202"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"interpolate":{"value":"step-before"},"fill":{"value":"steelblue"},"fillOpacity":{"value":0.4},"stroke":{"value":"steelblue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"c99cb365-6dc0-4884-8b15-335f8817c202","values":[{"x":1.0,"y":0},{"x":1.75,"y":40.0},{"x":2.5,"y":25.0},{"x":3.25,"y":10.0},{"x":4.0,"y":0.0},{"x":4.75,"y":9.0},{"x":5.5,"y":8.0},{"x":6.25,"y":6.0},{"x":7.0,"y":0.0},{"x":7.75,"y":2.0},{"x":8.5,"y":0}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"c99cb365-6dc0-4884-8b15-335f8817c202\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"c99cb365-6dc0-4884-8b15-335f8817c202\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"c99cb365-6dc0-4884-8b15-335f8817c202\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :interpolate {:value \"step-before\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 0.4}, :stroke {:value \"steelblue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"c99cb365-6dc0-4884-8b15-335f8817c202\", :values ({:x 1.0, :y 0} {:x 1.75, :y 40.0} {:x 2.5, :y 25.0} {:x 3.25, :y 10.0} {:x 4.0, :y 0.0} {:x 4.75, :y 9.0} {:x 5.5, :y 8.0} {:x 6.25, :y 6.0} {:x 7.0, :y 0.0} {:x 7.75, :y 2.0} {:x 8.5, :y 0})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(plot/list-plot (:mean (:complexity @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"1f88aa35-bcfb-4441-8842-5fc5a21ab852","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"1f88aa35-bcfb-4441-8842-5fc5a21ab852","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"1f88aa35-bcfb-4441-8842-5fc5a21ab852"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"1f88aa35-bcfb-4441-8842-5fc5a21ab852","values":[{"x":0,"y":4.46},{"x":1,"y":3.97},{"x":2,"y":3.47},{"x":3,"y":3.16},{"x":4,"y":3.115},{"x":5,"y":3.67},{"x":6,"y":3.495},{"x":7,"y":2.58},{"x":8,"y":2.61},{"x":9,"y":1.695},{"x":10,"y":1.785},{"x":11,"y":1.79},{"x":12,"y":1.585},{"x":13,"y":1.965},{"x":14,"y":2.425},{"x":15,"y":2.22},{"x":16,"y":2.16},{"x":17,"y":2.285},{"x":18,"y":2.695},{"x":19,"y":2.665},{"x":20,"y":2.11},{"x":21,"y":2.365},{"x":22,"y":2.305},{"x":23,"y":1.89},{"x":24,"y":1.64},{"x":25,"y":1.55},{"x":26,"y":1.73},{"x":27,"y":1.71},{"x":28,"y":1.72},{"x":29,"y":1.76},{"x":30,"y":1.905},{"x":31,"y":2.03},{"x":32,"y":1.955},{"x":33,"y":1.705},{"x":34,"y":1.785},{"x":35,"y":1.865},{"x":36,"y":1.84},{"x":37,"y":1.95},{"x":38,"y":1.81},{"x":39,"y":1.87},{"x":40,"y":2.18},{"x":41,"y":2.04},{"x":42,"y":2.245},{"x":43,"y":2.395},{"x":44,"y":2.18},{"x":45,"y":1.785},{"x":46,"y":1.725},{"x":47,"y":2.23},{"x":48,"y":2.14},{"x":49,"y":2.16},{"x":50,"y":2.28},{"x":51,"y":2.56},{"x":52,"y":2.44},{"x":53,"y":2.535},{"x":54,"y":2.485},{"x":55,"y":2.79},{"x":56,"y":2.645},{"x":57,"y":2.695},{"x":58,"y":2.795},{"x":59,"y":2.785},{"x":60,"y":2.83},{"x":61,"y":2.905},{"x":62,"y":2.675},{"x":63,"y":2.845},{"x":64,"y":2.88},{"x":65,"y":2.885},{"x":66,"y":2.72},{"x":67,"y":3.065},{"x":68,"y":2.49},{"x":69,"y":2.43},{"x":70,"y":2.25},{"x":71,"y":2.195},{"x":72,"y":2.4},{"x":73,"y":2.41},{"x":74,"y":2.51},{"x":75,"y":2.46},{"x":76,"y":2.555},{"x":77,"y":2.62},{"x":78,"y":2.69},{"x":79,"y":2.58},{"x":80,"y":2.81},{"x":81,"y":2.605},{"x":82,"y":2.54},{"x":83,"y":2.5},{"x":84,"y":2.54},{"x":85,"y":2.36},{"x":86,"y":2.505},{"x":87,"y":2.435},{"x":88,"y":2.46},{"x":89,"y":2.54},{"x":90,"y":2.45},{"x":91,"y":2.58},{"x":92,"y":2.515},{"x":93,"y":2.71},{"x":94,"y":2.58},{"x":95,"y":2.605},{"x":96,"y":2.62},{"x":97,"y":2.77},{"x":98,"y":2.535},{"x":99,"y":2.45},{"x":100,"y":2.79},{"x":101,"y":2.505},{"x":102,"y":2.735},{"x":103,"y":2.575},{"x":104,"y":2.38},{"x":105,"y":2.48},{"x":106,"y":2.53},{"x":107,"y":2.47},{"x":108,"y":2.59},{"x":109,"y":2.725},{"x":110,"y":2.63},{"x":111,"y":2.51},{"x":112,"y":2.48},{"x":113,"y":1.99},{"x":114,"y":1.97},{"x":115,"y":2.07},{"x":116,"y":2.035},{"x":117,"y":1.95},{"x":118,"y":2.155},{"x":119,"y":2.115},{"x":120,"y":2.07},{"x":121,"y":2.085},{"x":122,"y":1.96},{"x":123,"y":2.04},{"x":124,"y":2.04},{"x":125,"y":2.045},{"x":126,"y":2.1},{"x":127,"y":2.15},{"x":128,"y":2.05},{"x":129,"y":2.07},{"x":130,"y":2.155},{"x":131,"y":2.28},{"x":132,"y":2.21},{"x":133,"y":2.105},{"x":134,"y":2.12},{"x":135,"y":2.095},{"x":136,"y":2.065},{"x":137,"y":2.08},{"x":138,"y":2.1},{"x":139,"y":2.12},{"x":140,"y":2.18},{"x":141,"y":2.01},{"x":142,"y":2.02},{"x":143,"y":2.21},{"x":144,"y":2.13},{"x":145,"y":2.095},{"x":146,"y":2.19},{"x":147,"y":2.305},{"x":148,"y":2.05},{"x":149,"y":2.34},{"x":150,"y":2.145},{"x":151,"y":2.36},{"x":152,"y":2.2},{"x":153,"y":2.395},{"x":154,"y":2.155},{"x":155,"y":2.28},{"x":156,"y":2.12},{"x":157,"y":2.44},{"x":158,"y":2.425},{"x":159,"y":2.155},{"x":160,"y":2.3},{"x":161,"y":2.135},{"x":162,"y":2.24},{"x":163,"y":2.245},{"x":164,"y":2.2},{"x":165,"y":2.155},{"x":166,"y":2.37},{"x":167,"y":2.24},{"x":168,"y":2.045},{"x":169,"y":2.125},{"x":170,"y":2.08},{"x":171,"y":2.37},{"x":172,"y":2.125},{"x":173,"y":2.135},{"x":174,"y":2.225},{"x":175,"y":2.135},{"x":176,"y":2.38},{"x":177,"y":2.295},{"x":178,"y":2.315},{"x":179,"y":2.13},{"x":180,"y":2.19},{"x":181,"y":2.085},{"x":182,"y":2.3},{"x":183,"y":2.115},{"x":184,"y":2.115},{"x":185,"y":2.14},{"x":186,"y":2.345},{"x":187,"y":2.295},{"x":188,"y":2.255},{"x":189,"y":2.19},{"x":190,"y":2.145},{"x":191,"y":2.245},{"x":192,"y":2.19},{"x":193,"y":2.155},{"x":194,"y":2.1},{"x":195,"y":2.165},{"x":196,"y":1.97},{"x":197,"y":2.295},{"x":198,"y":2.175},{"x":199,"y":2.16}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"1f88aa35-bcfb-4441-8842-5fc5a21ab852\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"1f88aa35-bcfb-4441-8842-5fc5a21ab852\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"1f88aa35-bcfb-4441-8842-5fc5a21ab852\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"1f88aa35-bcfb-4441-8842-5fc5a21ab852\", :values ({:x 0, :y 4.46} {:x 1, :y 3.97} {:x 2, :y 3.47} {:x 3, :y 3.16} {:x 4, :y 3.115} {:x 5, :y 3.67} {:x 6, :y 3.495} {:x 7, :y 2.58} {:x 8, :y 2.61} {:x 9, :y 1.695} {:x 10, :y 1.785} {:x 11, :y 1.79} {:x 12, :y 1.585} {:x 13, :y 1.965} {:x 14, :y 2.425} {:x 15, :y 2.22} {:x 16, :y 2.16} {:x 17, :y 2.285} {:x 18, :y 2.695} {:x 19, :y 2.665} {:x 20, :y 2.11} {:x 21, :y 2.365} {:x 22, :y 2.305} {:x 23, :y 1.89} {:x 24, :y 1.64} {:x 25, :y 1.55} {:x 26, :y 1.73} {:x 27, :y 1.71} {:x 28, :y 1.72} {:x 29, :y 1.76} {:x 30, :y 1.905} {:x 31, :y 2.03} {:x 32, :y 1.955} {:x 33, :y 1.705} {:x 34, :y 1.785} {:x 35, :y 1.865} {:x 36, :y 1.84} {:x 37, :y 1.95} {:x 38, :y 1.81} {:x 39, :y 1.87} {:x 40, :y 2.18} {:x 41, :y 2.04} {:x 42, :y 2.245} {:x 43, :y 2.395} {:x 44, :y 2.18} {:x 45, :y 1.785} {:x 46, :y 1.725} {:x 47, :y 2.23} {:x 48, :y 2.14} {:x 49, :y 2.16} {:x 50, :y 2.28} {:x 51, :y 2.56} {:x 52, :y 2.44} {:x 53, :y 2.535} {:x 54, :y 2.485} {:x 55, :y 2.79} {:x 56, :y 2.645} {:x 57, :y 2.695} {:x 58, :y 2.795} {:x 59, :y 2.785} {:x 60, :y 2.83} {:x 61, :y 2.905} {:x 62, :y 2.675} {:x 63, :y 2.845} {:x 64, :y 2.88} {:x 65, :y 2.885} {:x 66, :y 2.72} {:x 67, :y 3.065} {:x 68, :y 2.49} {:x 69, :y 2.43} {:x 70, :y 2.25} {:x 71, :y 2.195} {:x 72, :y 2.4} {:x 73, :y 2.41} {:x 74, :y 2.51} {:x 75, :y 2.46} {:x 76, :y 2.555} {:x 77, :y 2.62} {:x 78, :y 2.69} {:x 79, :y 2.58} {:x 80, :y 2.81} {:x 81, :y 2.605} {:x 82, :y 2.54} {:x 83, :y 2.5} {:x 84, :y 2.54} {:x 85, :y 2.36} {:x 86, :y 2.505} {:x 87, :y 2.435} {:x 88, :y 2.46} {:x 89, :y 2.54} {:x 90, :y 2.45} {:x 91, :y 2.58} {:x 92, :y 2.515} {:x 93, :y 2.71} {:x 94, :y 2.58} {:x 95, :y 2.605} {:x 96, :y 2.62} {:x 97, :y 2.77} {:x 98, :y 2.535} {:x 99, :y 2.45} {:x 100, :y 2.79} {:x 101, :y 2.505} {:x 102, :y 2.735} {:x 103, :y 2.575} {:x 104, :y 2.38} {:x 105, :y 2.48} {:x 106, :y 2.53} {:x 107, :y 2.47} {:x 108, :y 2.59} {:x 109, :y 2.725} {:x 110, :y 2.63} {:x 111, :y 2.51} {:x 112, :y 2.48} {:x 113, :y 1.99} {:x 114, :y 1.97} {:x 115, :y 2.07} {:x 116, :y 2.035} {:x 117, :y 1.95} {:x 118, :y 2.155} {:x 119, :y 2.115} {:x 120, :y 2.07} {:x 121, :y 2.085} {:x 122, :y 1.96} {:x 123, :y 2.04} {:x 124, :y 2.04} {:x 125, :y 2.045} {:x 126, :y 2.1} {:x 127, :y 2.15} {:x 128, :y 2.05} {:x 129, :y 2.07} {:x 130, :y 2.155} {:x 131, :y 2.28} {:x 132, :y 2.21} {:x 133, :y 2.105} {:x 134, :y 2.12} {:x 135, :y 2.095} {:x 136, :y 2.065} {:x 137, :y 2.08} {:x 138, :y 2.1} {:x 139, :y 2.12} {:x 140, :y 2.18} {:x 141, :y 2.01} {:x 142, :y 2.02} {:x 143, :y 2.21} {:x 144, :y 2.13} {:x 145, :y 2.095} {:x 146, :y 2.19} {:x 147, :y 2.305} {:x 148, :y 2.05} {:x 149, :y 2.34} {:x 150, :y 2.145} {:x 151, :y 2.36} {:x 152, :y 2.2} {:x 153, :y 2.395} {:x 154, :y 2.155} {:x 155, :y 2.28} {:x 156, :y 2.12} {:x 157, :y 2.44} {:x 158, :y 2.425} {:x 159, :y 2.155} {:x 160, :y 2.3} {:x 161, :y 2.135} {:x 162, :y 2.24} {:x 163, :y 2.245} {:x 164, :y 2.2} {:x 165, :y 2.155} {:x 166, :y 2.37} {:x 167, :y 2.24} {:x 168, :y 2.045} {:x 169, :y 2.125} {:x 170, :y 2.08} {:x 171, :y 2.37} {:x 172, :y 2.125} {:x 173, :y 2.135} {:x 174, :y 2.225} {:x 175, :y 2.135} {:x 176, :y 2.38} {:x 177, :y 2.295} {:x 178, :y 2.315} {:x 179, :y 2.13} {:x 180, :y 2.19} {:x 181, :y 2.085} {:x 182, :y 2.3} {:x 183, :y 2.115} {:x 184, :y 2.115} {:x 185, :y 2.14} {:x 186, :y 2.345} {:x 187, :y 2.295} {:x 188, :y 2.255} {:x 189, :y 2.19} {:x 190, :y 2.145} {:x 191, :y 2.245} {:x 192, :y 2.19} {:x 193, :y 2.155} {:x 194, :y 2.1} {:x 195, :y 2.165} {:x 196, :y 1.97} {:x 197, :y 2.295} {:x 198, :y 2.175} {:x 199, :y 2.16})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@

;; @@

;; @@

;; @@
