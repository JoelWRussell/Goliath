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
(def OS (System/getProperty "os.name"))
(if (= OS "Mac OS X")
(def mathKernelSz "/Applications/Mathematica.app/Contents/MacOS/MathKernel")  
(def mathKernelSz "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe"))

OS 
mathKernelSz
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;c:/program files/wolfram research/mathematica/10.0/mathkernel.exe&quot;</span>","value":"\"c:/program files/wolfram research/mathematica/10.0/mathkernel.exe\""}
;; <=

;; **
;;; CHOOSE THE DATA - NEEDS TO BE IN RESOURCES/
;; **

;; @@
(def experimentalDataSz "sho_coupled_0.1_sim.csv")
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
(def df 2)
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
;;this function uses discrete differential to find speeds 
;;(goliath.mathlink.LagrangianScore/InitFunctions
  ;;mathKernelSz
  ;;"resources/"
  ;;experimentalDataSz
  ;;dt
  ;;df
  ;;)
;; @@

;; @@
;;this function reads in the speed accel from the initial data file
;;in the form th1, th2, th3, w1, w2, w3, a1, a2, a3
(goliath.mathlink.LagrangianScore/InitFunctionsSim
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
  
;;uncomment this to score a particular lagrangian eg coupled sho
;;(into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (vec (flatten [[2 0 0 ;;0] [0 2 0 0] [1 1 0 0] [0 0 2 0] [0 0 0 2]]))), df))


;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-22.806634470647438</span>","value":"-22.806634470647438"},{"type":"html","content":"<span class='clj-double'>0.18854088335408387</span>","value":"0.18854088335408387"},{"type":"html","content":"<span class='clj-double'>0.09427039337015511</span>","value":"0.09427039337015511"},{"type":"html","content":"<span class='clj-double'>-0.188540811738729</span>","value":"-0.188540811738729"},{"type":"html","content":"<span class='clj-double'>-0.09427045853381483</span>","value":"-0.09427045853381483"},{"type":"html","content":"<span class='clj-double'>-0.09427042790594096</span>","value":"-0.09427042790594096"}],"value":"[-22.806634470647438 0.18854088335408387 0.09427039337015511 -0.188540811738729 -0.09427045853381483 -0.09427042790594096]"}
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
;;; 0.1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22.23.24.25.26.27.28.29.30.31.32.33.34.35.36.37.38.39.40.41.42.43.44.45.46.47.48.49.50.51.52.53.54.55.56.57.58.59.60.61.62.63.64.65.66.67.68.69.70.71.72.73.74.75.76.77.78.79.80.81.82.83.84.85.86.87.88.89.90.91.92.93.94.95.96.97.98.99.100.101.102.103.104.105.106.107.108.109.110.111.112.113.114.115.116.117.118.119.120.121.122.123.124.125.126.127.128.129.130.131.132.133.134.135.136.137.138.139.140.141.142.143.144.145.146.147.148.149.150.151.152.153.154.155.156.157.158.159.160.161.162.163.164.165.166.167.168.169.170.171.172.173.174.175.176.177.178.179.180.181.182.183.184.185.186.187.188.189.190.191.192.193.194.195.196.197.198.199.200&quot;Elapsed time: 332433.960161 msecs&quot;
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
;;; [[1 0] [0 6] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [0 6] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [0 6] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [0 6] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [8 1] [5 1] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [2 1] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [2 1] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [2 1] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [2 1] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [7 0] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [1 6] [5 0]]
;;; [[1 0] [3 4] [3 0] [1 4] [1 6] [5 0]]
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
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac","values":[{"x":0,"y":-13.453269848243083},{"x":1,"y":-13.730858444500203},{"x":2,"y":-13.73144385256534},{"x":3,"y":-15.584717650002066},{"x":4,"y":-15.584717650002066},{"x":5,"y":-15.584717650002066},{"x":6,"y":-17.16164893875585},{"x":7,"y":-17.16164893875585},{"x":8,"y":-17.162734363927875},{"x":9,"y":-17.16334302491667},{"x":10,"y":-17.16334302491667},{"x":11,"y":-17.16334302491667},{"x":12,"y":-17.186362013587612},{"x":13,"y":-17.186362013587612},{"x":14,"y":-17.186362039735414},{"x":15,"y":-17.189251088136636},{"x":16,"y":-17.189251088136636},{"x":17,"y":-17.19373590535055},{"x":18,"y":-17.193741084204966},{"x":19,"y":-17.193741099003947},{"x":20,"y":-17.193741099003947},{"x":21,"y":-17.193741099003947},{"x":22,"y":-17.198154753975913},{"x":23,"y":-17.198154753975913},{"x":24,"y":-17.198154753975913},{"x":25,"y":-17.198154753975913},{"x":26,"y":-17.20303580114171},{"x":27,"y":-17.20303580114171},{"x":28,"y":-17.20303580114171},{"x":29,"y":-17.20303580114171},{"x":30,"y":-17.20303580114171},{"x":31,"y":-17.20303580114171},{"x":32,"y":-21.506399651969573},{"x":33,"y":-21.506399651969573},{"x":34,"y":-21.50688878464267},{"x":35,"y":-21.50688878464267},{"x":36,"y":-21.50688878464267},{"x":37,"y":-21.52256485891552},{"x":38,"y":-21.52256485891552},{"x":39,"y":-21.52256485891552},{"x":40,"y":-21.52256485891552},{"x":41,"y":-21.52256485891552},{"x":42,"y":-21.52256485891552},{"x":43,"y":-21.52256485891552},{"x":44,"y":-21.52256485891552},{"x":45,"y":-21.52256485891552},{"x":46,"y":-21.52256485891552},{"x":47,"y":-21.52256485891552},{"x":48,"y":-21.52256485891552},{"x":49,"y":-21.52256485891552},{"x":50,"y":-21.54864509706405},{"x":51,"y":-21.54864509706405},{"x":52,"y":-21.54864509706405},{"x":53,"y":-21.54864509706405},{"x":54,"y":-21.5517089211401},{"x":55,"y":-21.551915634645198},{"x":56,"y":-22.13385893222708},{"x":57,"y":-22.13385893222708},{"x":58,"y":-22.13385893222708},{"x":59,"y":-22.13385893222708},{"x":60,"y":-22.13385893222708},{"x":61,"y":-22.13385893222708},{"x":62,"y":-22.13385893222708},{"x":63,"y":-22.13385893222708},{"x":64,"y":-22.13385893222708},{"x":65,"y":-22.13385893222708},{"x":66,"y":-22.13385893222708},{"x":67,"y":-22.13385893222708},{"x":68,"y":-22.13385893222708},{"x":69,"y":-22.13385893222708},{"x":70,"y":-22.13385893222708},{"x":71,"y":-22.13385893222708},{"x":72,"y":-22.13385893222708},{"x":73,"y":-22.13385893222708},{"x":74,"y":-22.13385893222708},{"x":75,"y":-22.13385893222708},{"x":76,"y":-22.13385893222708},{"x":77,"y":-22.13385893222708},{"x":78,"y":-22.133860092466428},{"x":79,"y":-22.133860092466428},{"x":80,"y":-22.557787962725666},{"x":81,"y":-22.557787962725666},{"x":82,"y":-22.557787962725666},{"x":83,"y":-22.557787962725666},{"x":84,"y":-22.557787962725666},{"x":85,"y":-22.55783479225662},{"x":86,"y":-22.55783479225662},{"x":87,"y":-22.55783479225662},{"x":88,"y":-22.55783479225662},{"x":89,"y":-22.55783479225662},{"x":90,"y":-22.55783479225662},{"x":91,"y":-22.55783479225662},{"x":92,"y":-22.595511390236283},{"x":93,"y":-22.595511390236283},{"x":94,"y":-22.595511390236283},{"x":95,"y":-22.595511390236283},{"x":96,"y":-22.595511390236283},{"x":97,"y":-22.595511390236283},{"x":98,"y":-22.595511390236283},{"x":99,"y":-22.982052845955696},{"x":100,"y":-22.982052845955696},{"x":101,"y":-22.982052845955696},{"x":102,"y":-23.044537969929024},{"x":103,"y":-23.044537969929024},{"x":104,"y":-23.164076180204987},{"x":105,"y":-23.164076180204987},{"x":106,"y":-23.164076180204987},{"x":107,"y":-23.164076180204987},{"x":108,"y":-23.164076180204987},{"x":109,"y":-23.16429057231151},{"x":110,"y":-23.164911642458197},{"x":111,"y":-23.164911642458197},{"x":112,"y":-23.164911642458197},{"x":113,"y":-23.164911642458197},{"x":114,"y":-23.164911642458197},{"x":115,"y":-23.164911642458197},{"x":116,"y":-23.165443835313848},{"x":117,"y":-23.175747404661305},{"x":118,"y":-23.175747404661305},{"x":119,"y":-23.175747404661305},{"x":120,"y":-23.175747404661305},{"x":121,"y":-23.175747404661305},{"x":122,"y":-23.176738636809237},{"x":123,"y":-23.176738636809237},{"x":124,"y":-23.176738636809237},{"x":125,"y":-23.176738636809237},{"x":126,"y":-23.176738636809237},{"x":127,"y":-23.176738636809237},{"x":128,"y":-23.176738636809237},{"x":129,"y":-23.176738636809237},{"x":130,"y":-23.176738636809237},{"x":131,"y":-23.176738636809237},{"x":132,"y":-23.176738636809237},{"x":133,"y":-23.176738636809237},{"x":134,"y":-23.176738636809237},{"x":135,"y":-23.176738636809237},{"x":136,"y":-23.176738636809237},{"x":137,"y":-23.176738636809237},{"x":138,"y":-23.176738636809237},{"x":139,"y":-23.176738636809237},{"x":140,"y":-23.176738636809237},{"x":141,"y":-23.218404259109793},{"x":142,"y":-23.218404259109793},{"x":143,"y":-23.218404259109793},{"x":144,"y":-23.218404259109793},{"x":145,"y":-23.218404259109793},{"x":146,"y":-23.218404259109793},{"x":147,"y":-23.218404259109793},{"x":148,"y":-23.218404259109793},{"x":149,"y":-23.218404259109793},{"x":150,"y":-23.218404259109793},{"x":151,"y":-23.221909523360445},{"x":152,"y":-23.22378476997795},{"x":153,"y":-23.22378476997795},{"x":154,"y":-23.22378476997795},{"x":155,"y":-23.22378476997795},{"x":156,"y":-23.22378476997795},{"x":157,"y":-23.22378476997795},{"x":158,"y":-23.22378476997795},{"x":159,"y":-23.22378476997795},{"x":160,"y":-23.22378476997795},{"x":161,"y":-23.22378476997795},{"x":162,"y":-23.224850190684602},{"x":163,"y":-23.224850190684602},{"x":164,"y":-23.224850190684602},{"x":165,"y":-23.224850190684602},{"x":166,"y":-23.224850190684602},{"x":167,"y":-23.224850190684602},{"x":168,"y":-23.224850190684602},{"x":169,"y":-23.224850190684602},{"x":170,"y":-23.224850190684602},{"x":171,"y":-23.225246173229735},{"x":172,"y":-23.225246173229735},{"x":173,"y":-23.225246173229735},{"x":174,"y":-23.225246173229735},{"x":175,"y":-23.225246173229735},{"x":176,"y":-23.225246173229735},{"x":177,"y":-23.225246173229735},{"x":178,"y":-23.225246173229735},{"x":179,"y":-23.225246173229735},{"x":180,"y":-23.225246173229735},{"x":181,"y":-23.225246173229735},{"x":182,"y":-23.225246173229735},{"x":183,"y":-23.225246173229735},{"x":184,"y":-23.225246173229735},{"x":185,"y":-23.225246173229735},{"x":186,"y":-23.225246173229735},{"x":187,"y":-23.225246173229735},{"x":188,"y":-23.225246173229735},{"x":189,"y":-23.225246173229735},{"x":190,"y":-23.225246173229735},{"x":191,"y":-23.225246173229735},{"x":192,"y":-23.225246173229735},{"x":193,"y":-23.225246173229735},{"x":194,"y":-23.225246173229735},{"x":195,"y":-23.225246173229735},{"x":196,"y":-23.225246173229735},{"x":197,"y":-23.225246173229735},{"x":198,"y":-23.225246173229735},{"x":199,"y":-23.225246173229735}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"8d2b5c2b-e2b7-4a28-a395-9e57d6609dac\", :values ({:x 0, :y -13.453269848243083} {:x 1, :y -13.730858444500203} {:x 2, :y -13.73144385256534} {:x 3, :y -15.584717650002066} {:x 4, :y -15.584717650002066} {:x 5, :y -15.584717650002066} {:x 6, :y -17.16164893875585} {:x 7, :y -17.16164893875585} {:x 8, :y -17.162734363927875} {:x 9, :y -17.16334302491667} {:x 10, :y -17.16334302491667} {:x 11, :y -17.16334302491667} {:x 12, :y -17.186362013587612} {:x 13, :y -17.186362013587612} {:x 14, :y -17.186362039735414} {:x 15, :y -17.189251088136636} {:x 16, :y -17.189251088136636} {:x 17, :y -17.19373590535055} {:x 18, :y -17.193741084204966} {:x 19, :y -17.193741099003947} {:x 20, :y -17.193741099003947} {:x 21, :y -17.193741099003947} {:x 22, :y -17.198154753975913} {:x 23, :y -17.198154753975913} {:x 24, :y -17.198154753975913} {:x 25, :y -17.198154753975913} {:x 26, :y -17.20303580114171} {:x 27, :y -17.20303580114171} {:x 28, :y -17.20303580114171} {:x 29, :y -17.20303580114171} {:x 30, :y -17.20303580114171} {:x 31, :y -17.20303580114171} {:x 32, :y -21.506399651969573} {:x 33, :y -21.506399651969573} {:x 34, :y -21.50688878464267} {:x 35, :y -21.50688878464267} {:x 36, :y -21.50688878464267} {:x 37, :y -21.52256485891552} {:x 38, :y -21.52256485891552} {:x 39, :y -21.52256485891552} {:x 40, :y -21.52256485891552} {:x 41, :y -21.52256485891552} {:x 42, :y -21.52256485891552} {:x 43, :y -21.52256485891552} {:x 44, :y -21.52256485891552} {:x 45, :y -21.52256485891552} {:x 46, :y -21.52256485891552} {:x 47, :y -21.52256485891552} {:x 48, :y -21.52256485891552} {:x 49, :y -21.52256485891552} {:x 50, :y -21.54864509706405} {:x 51, :y -21.54864509706405} {:x 52, :y -21.54864509706405} {:x 53, :y -21.54864509706405} {:x 54, :y -21.5517089211401} {:x 55, :y -21.551915634645198} {:x 56, :y -22.13385893222708} {:x 57, :y -22.13385893222708} {:x 58, :y -22.13385893222708} {:x 59, :y -22.13385893222708} {:x 60, :y -22.13385893222708} {:x 61, :y -22.13385893222708} {:x 62, :y -22.13385893222708} {:x 63, :y -22.13385893222708} {:x 64, :y -22.13385893222708} {:x 65, :y -22.13385893222708} {:x 66, :y -22.13385893222708} {:x 67, :y -22.13385893222708} {:x 68, :y -22.13385893222708} {:x 69, :y -22.13385893222708} {:x 70, :y -22.13385893222708} {:x 71, :y -22.13385893222708} {:x 72, :y -22.13385893222708} {:x 73, :y -22.13385893222708} {:x 74, :y -22.13385893222708} {:x 75, :y -22.13385893222708} {:x 76, :y -22.13385893222708} {:x 77, :y -22.13385893222708} {:x 78, :y -22.133860092466428} {:x 79, :y -22.133860092466428} {:x 80, :y -22.557787962725666} {:x 81, :y -22.557787962725666} {:x 82, :y -22.557787962725666} {:x 83, :y -22.557787962725666} {:x 84, :y -22.557787962725666} {:x 85, :y -22.55783479225662} {:x 86, :y -22.55783479225662} {:x 87, :y -22.55783479225662} {:x 88, :y -22.55783479225662} {:x 89, :y -22.55783479225662} {:x 90, :y -22.55783479225662} {:x 91, :y -22.55783479225662} {:x 92, :y -22.595511390236283} {:x 93, :y -22.595511390236283} {:x 94, :y -22.595511390236283} {:x 95, :y -22.595511390236283} {:x 96, :y -22.595511390236283} {:x 97, :y -22.595511390236283} {:x 98, :y -22.595511390236283} {:x 99, :y -22.982052845955696} {:x 100, :y -22.982052845955696} {:x 101, :y -22.982052845955696} {:x 102, :y -23.044537969929024} {:x 103, :y -23.044537969929024} {:x 104, :y -23.164076180204987} {:x 105, :y -23.164076180204987} {:x 106, :y -23.164076180204987} {:x 107, :y -23.164076180204987} {:x 108, :y -23.164076180204987} {:x 109, :y -23.16429057231151} {:x 110, :y -23.164911642458197} {:x 111, :y -23.164911642458197} {:x 112, :y -23.164911642458197} {:x 113, :y -23.164911642458197} {:x 114, :y -23.164911642458197} {:x 115, :y -23.164911642458197} {:x 116, :y -23.165443835313848} {:x 117, :y -23.175747404661305} {:x 118, :y -23.175747404661305} {:x 119, :y -23.175747404661305} {:x 120, :y -23.175747404661305} {:x 121, :y -23.175747404661305} {:x 122, :y -23.176738636809237} {:x 123, :y -23.176738636809237} {:x 124, :y -23.176738636809237} {:x 125, :y -23.176738636809237} {:x 126, :y -23.176738636809237} {:x 127, :y -23.176738636809237} {:x 128, :y -23.176738636809237} {:x 129, :y -23.176738636809237} {:x 130, :y -23.176738636809237} {:x 131, :y -23.176738636809237} {:x 132, :y -23.176738636809237} {:x 133, :y -23.176738636809237} {:x 134, :y -23.176738636809237} {:x 135, :y -23.176738636809237} {:x 136, :y -23.176738636809237} {:x 137, :y -23.176738636809237} {:x 138, :y -23.176738636809237} {:x 139, :y -23.176738636809237} {:x 140, :y -23.176738636809237} {:x 141, :y -23.218404259109793} {:x 142, :y -23.218404259109793} {:x 143, :y -23.218404259109793} {:x 144, :y -23.218404259109793} {:x 145, :y -23.218404259109793} {:x 146, :y -23.218404259109793} {:x 147, :y -23.218404259109793} {:x 148, :y -23.218404259109793} {:x 149, :y -23.218404259109793} {:x 150, :y -23.218404259109793} {:x 151, :y -23.221909523360445} {:x 152, :y -23.22378476997795} {:x 153, :y -23.22378476997795} {:x 154, :y -23.22378476997795} {:x 155, :y -23.22378476997795} {:x 156, :y -23.22378476997795} {:x 157, :y -23.22378476997795} {:x 158, :y -23.22378476997795} {:x 159, :y -23.22378476997795} {:x 160, :y -23.22378476997795} {:x 161, :y -23.22378476997795} {:x 162, :y -23.224850190684602} {:x 163, :y -23.224850190684602} {:x 164, :y -23.224850190684602} {:x 165, :y -23.224850190684602} {:x 166, :y -23.224850190684602} {:x 167, :y -23.224850190684602} {:x 168, :y -23.224850190684602} {:x 169, :y -23.224850190684602} {:x 170, :y -23.224850190684602} {:x 171, :y -23.225246173229735} {:x 172, :y -23.225246173229735} {:x 173, :y -23.225246173229735} {:x 174, :y -23.225246173229735} {:x 175, :y -23.225246173229735} {:x 176, :y -23.225246173229735} {:x 177, :y -23.225246173229735} {:x 178, :y -23.225246173229735} {:x 179, :y -23.225246173229735} {:x 180, :y -23.225246173229735} {:x 181, :y -23.225246173229735} {:x 182, :y -23.225246173229735} {:x 183, :y -23.225246173229735} {:x 184, :y -23.225246173229735} {:x 185, :y -23.225246173229735} {:x 186, :y -23.225246173229735} {:x 187, :y -23.225246173229735} {:x 188, :y -23.225246173229735} {:x 189, :y -23.225246173229735} {:x 190, :y -23.225246173229735} {:x 191, :y -23.225246173229735} {:x 192, :y -23.225246173229735} {:x 193, :y -23.225246173229735} {:x 194, :y -23.225246173229735} {:x 195, :y -23.225246173229735} {:x 196, :y -23.225246173229735} {:x 197, :y -23.225246173229735} {:x 198, :y -23.225246173229735} {:x 199, :y -23.225246173229735})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
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
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"5483196b-791a-4411-b3ea-019c17e07a57","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"5483196b-791a-4411-b3ea-019c17e07a57","field":"data.y"}}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"5483196b-791a-4411-b3ea-019c17e07a57","values":[{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-23.225246173229735,"y":10},{"x":-23.225246173229735,"y":10},{"x":-23.225246173229735,"y":10},{"x":-23.225246173229735,"y":10},{"x":-23.221909523360445,"y":9},{"x":-23.221909523360445,"y":9},{"x":-23.221909523360445,"y":9},{"x":-23.221909523360445,"y":9},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-23.19078196871628,"y":8},{"x":-23.19078196871628,"y":8},{"x":-23.19078196871628,"y":8},{"x":-23.19078196871628,"y":8},{"x":-23.154224200751045,"y":7},{"x":-23.154224200751045,"y":7},{"x":-23.154224200751045,"y":7},{"x":-23.154224200751045,"y":7},{"x":-14.167399480063514,"y":4},{"x":-14.167399480063514,"y":4},{"x":-14.167399480063514,"y":4},{"x":-14.167399480063514,"y":4},{"x":-19.85861158060795,"y":6},{"x":-19.85861158060795,"y":6},{"x":-19.85861158060795,"y":6},{"x":-19.85861158060795,"y":6},{"x":-17.124560966399606,"y":5},{"x":-17.124560966399606,"y":5},{"x":-17.124560966399606,"y":5},{"x":-17.124560966399606,"y":5}]},{"name":"0eadf1e0-e300-4fe2-a063-f3e8f363476f","values":[{"x":-23.153654351149367,"y":8},{"x":21.984319949170665,"y":1},{"x":-17.945000402964766,"y":9},{"x":-12.766793869191543,"y":6},{"x":1.5832220019188459,"y":1},{"x":0.47262976349464153,"y":1},{"x":-17.753572045283583,"y":8},{"x":0.47262976349464153,"y":1},{"x":-1.7395713410823415,"y":6},{"x":-13.723159580923811,"y":4},{"x":0.47262976349464153,"y":1},{"x":0.4554564099438658,"y":2},{"x":6.998525253793526,"y":1},{"x":0.4562542131205834,"y":2},{"x":0.46434696960927796,"y":2},{"x":-17.945275235674583,"y":10},{"x":0.47262976349464153,"y":1},{"x":-13.726669499420082,"y":5},{"x":-13.730178123229146,"y":4},{"x":6.7906092875466335,"y":2},{"x":0.47262976349464153,"y":1},{"x":-17.124560966399606,"y":5},{"x":0.46147955715973793,"y":2},{"x":-13.706087557928639,"y":5},{"x":-13.727771816573156,"y":4},{"x":-0.47619100041314355,"y":2},{"x":0.16182588229568615,"y":3},{"x":-23.160830581284856,"y":8},{"x":-13.723247134694176,"y":4},{"x":0.7086898265293349,"y":1},{"x":0.11343839059317351,"y":3},{"x":-13.558707859720641,"y":8},{"x":-0.740083692949504,"y":5},{"x":-2.056665513418497,"y":8},{"x":-0.8685993201272337,"y":4},{"x":-0.2780919291524708,"y":4},{"x":-1.832848857927833,"y":6},{"x":-0.47634843321060716,"y":3},{"x":0.45210327060453,"y":2},{"x":-17.945809039777764,"y":10},{"x":-0.2875411879323268,"y":2},{"x":-17.75512722767535,"y":9},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.493942291559064,"y":5},{"x":0.1421918741178639,"y":6},{"x":-13.706087557928639,"y":5},{"x":0.46147955715973793,"y":2},{"x":-14.17418524331605,"y":6},{"x":0.7254101236533876,"y":3},{"x":-23.154224200751045,"y":7},{"x":0.47262976349464153,"y":1},{"x":-19.858934041558648,"y":8},{"x":-23.154224200751045,"y":7},{"x":-1.8896083471520773,"y":8},{"x":-18.769826969991453,"y":8},{"x":-13.731358576363322,"y":4},{"x":-13.723159580923811,"y":4},{"x":-14.167399480063514,"y":4},{"x":-0.926518669700218,"y":4},{"x":-1.838307190105307,"y":7},{"x":-0.5957824517330587,"y":3},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.24055151149046236,"y":3},{"x":-13.706087557928639,"y":5},{"x":0.7436440291896712,"y":2},{"x":-0.5957824517330587,"y":3},{"x":0.22815070122309758,"y":2},{"x":-19.256328807854455,"y":6},{"x":-0.9490445270506043,"y":4},{"x":0.4520639533120918,"y":3},{"x":0.008023389163842473,"y":3},{"x":-13.60618586389685,"y":5},{"x":0.47262976349464153,"y":1},{"x":-19.85861158060795,"y":6},{"x":-0.6639467669447258,"y":4},{"x":-13.73744383394971,"y":6},{"x":-23.225246173229735,"y":10},{"x":-19.256328807854455,"y":6},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":6.980877274384624,"y":2},{"x":0.47262976349464153,"y":1},{"x":-13.441344954016571,"y":4},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":-13.44143819101242,"y":4},{"x":-13.702330517028352,"y":2}]},{"name":"67fc3182-7f4b-41ec-b577-0f112a8e9734","values":[{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-13.702330517028352,"y":2},{"x":-13.702330517028352,"y":2},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":0.47262976349464153,"y":1},{"x":-23.225246173229735,"y":10},{"x":-23.225246173229735,"y":10},{"x":-23.225246173229735,"y":10},{"x":-23.225246173229735,"y":10},{"x":-23.221909523360445,"y":9},{"x":-23.221909523360445,"y":9},{"x":-23.221909523360445,"y":9},{"x":-23.221909523360445,"y":9},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-13.727242679922018,"y":3},{"x":-23.19078196871628,"y":8},{"x":-23.19078196871628,"y":8},{"x":-23.19078196871628,"y":8},{"x":-23.19078196871628,"y":8},{"x":-23.154224200751045,"y":7},{"x":-23.154224200751045,"y":7},{"x":-23.154224200751045,"y":7},{"x":-23.154224200751045,"y":7},{"x":-14.167399480063514,"y":4},{"x":-14.167399480063514,"y":4},{"x":-14.167399480063514,"y":4},{"x":-14.167399480063514,"y":4},{"x":-19.85861158060795,"y":6},{"x":-19.85861158060795,"y":6},{"x":-19.85861158060795,"y":6},{"x":-19.85861158060795,"y":6},{"x":-17.124560966399606,"y":5},{"x":-17.124560966399606,"y":5},{"x":-17.124560966399606,"y":5},{"x":-17.124560966399606,"y":5}]}],"marks":[{"type":"symbol","from":{"data":"5483196b-791a-4411-b3ea-019c17e07a57"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"red"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"0eadf1e0-e300-4fe2-a063-f3e8f363476f"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"blue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"67fc3182-7f4b-41ec-b577-0f112a8e9734"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"#ff29d2"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"5483196b-791a-4411-b3ea-019c17e07a57\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"5483196b-791a-4411-b3ea-019c17e07a57\", :field \"data.y\"}}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"5483196b-791a-4411-b3ea-019c17e07a57\", :values ({:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -23.225246173229735, :y 10} {:x -23.225246173229735, :y 10} {:x -23.225246173229735, :y 10} {:x -23.225246173229735, :y 10} {:x -23.221909523360445, :y 9} {:x -23.221909523360445, :y 9} {:x -23.221909523360445, :y 9} {:x -23.221909523360445, :y 9} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -23.19078196871628, :y 8} {:x -23.19078196871628, :y 8} {:x -23.19078196871628, :y 8} {:x -23.19078196871628, :y 8} {:x -23.154224200751045, :y 7} {:x -23.154224200751045, :y 7} {:x -23.154224200751045, :y 7} {:x -23.154224200751045, :y 7} {:x -14.167399480063514, :y 4} {:x -14.167399480063514, :y 4} {:x -14.167399480063514, :y 4} {:x -14.167399480063514, :y 4} {:x -19.85861158060795, :y 6} {:x -19.85861158060795, :y 6} {:x -19.85861158060795, :y 6} {:x -19.85861158060795, :y 6} {:x -17.124560966399606, :y 5} {:x -17.124560966399606, :y 5} {:x -17.124560966399606, :y 5} {:x -17.124560966399606, :y 5})} {:name \"0eadf1e0-e300-4fe2-a063-f3e8f363476f\", :values ({:x -23.153654351149367, :y 8} {:x 21.984319949170665, :y 1} {:x -17.945000402964766, :y 9} {:x -12.766793869191543, :y 6} {:x 1.5832220019188459, :y 1} {:x 0.47262976349464153, :y 1} {:x -17.753572045283583, :y 8} {:x 0.47262976349464153, :y 1} {:x -1.7395713410823415, :y 6} {:x -13.723159580923811, :y 4} {:x 0.47262976349464153, :y 1} {:x 0.4554564099438658, :y 2} {:x 6.998525253793526, :y 1} {:x 0.4562542131205834, :y 2} {:x 0.46434696960927796, :y 2} {:x -17.945275235674583, :y 10} {:x 0.47262976349464153, :y 1} {:x -13.726669499420082, :y 5} {:x -13.730178123229146, :y 4} {:x 6.7906092875466335, :y 2} {:x 0.47262976349464153, :y 1} {:x -17.124560966399606, :y 5} {:x 0.46147955715973793, :y 2} {:x -13.706087557928639, :y 5} {:x -13.727771816573156, :y 4} {:x -0.47619100041314355, :y 2} {:x 0.16182588229568615, :y 3} {:x -23.160830581284856, :y 8} {:x -13.723247134694176, :y 4} {:x 0.7086898265293349, :y 1} {:x 0.11343839059317351, :y 3} {:x -13.558707859720641, :y 8} {:x -0.740083692949504, :y 5} {:x -2.056665513418497, :y 8} {:x -0.8685993201272337, :y 4} {:x -0.2780919291524708, :y 4} {:x -1.832848857927833, :y 6} {:x -0.47634843321060716, :y 3} {:x 0.45210327060453, :y 2} {:x -17.945809039777764, :y 10} {:x -0.2875411879323268, :y 2} {:x -17.75512722767535, :y 9} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.493942291559064, :y 5} {:x 0.1421918741178639, :y 6} {:x -13.706087557928639, :y 5} {:x 0.46147955715973793, :y 2} {:x -14.17418524331605, :y 6} {:x 0.7254101236533876, :y 3} {:x -23.154224200751045, :y 7} {:x 0.47262976349464153, :y 1} {:x -19.858934041558648, :y 8} {:x -23.154224200751045, :y 7} {:x -1.8896083471520773, :y 8} {:x -18.769826969991453, :y 8} {:x -13.731358576363322, :y 4} {:x -13.723159580923811, :y 4} {:x -14.167399480063514, :y 4} {:x -0.926518669700218, :y 4} {:x -1.838307190105307, :y 7} {:x -0.5957824517330587, :y 3} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.24055151149046236, :y 3} {:x -13.706087557928639, :y 5} {:x 0.7436440291896712, :y 2} {:x -0.5957824517330587, :y 3} {:x 0.22815070122309758, :y 2} {:x -19.256328807854455, :y 6} {:x -0.9490445270506043, :y 4} {:x 0.4520639533120918, :y 3} {:x 0.008023389163842473, :y 3} {:x -13.60618586389685, :y 5} {:x 0.47262976349464153, :y 1} {:x -19.85861158060795, :y 6} {:x -0.6639467669447258, :y 4} {:x -13.73744383394971, :y 6} {:x -23.225246173229735, :y 10} {:x -19.256328807854455, :y 6} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 6.980877274384624, :y 2} {:x 0.47262976349464153, :y 1} {:x -13.441344954016571, :y 4} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x -13.44143819101242, :y 4} {:x -13.702330517028352, :y 2})} {:name \"67fc3182-7f4b-41ec-b577-0f112a8e9734\", :values ({:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -13.702330517028352, :y 2} {:x -13.702330517028352, :y 2} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x 0.47262976349464153, :y 1} {:x -23.225246173229735, :y 10} {:x -23.225246173229735, :y 10} {:x -23.225246173229735, :y 10} {:x -23.225246173229735, :y 10} {:x -23.221909523360445, :y 9} {:x -23.221909523360445, :y 9} {:x -23.221909523360445, :y 9} {:x -23.221909523360445, :y 9} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -13.727242679922018, :y 3} {:x -23.19078196871628, :y 8} {:x -23.19078196871628, :y 8} {:x -23.19078196871628, :y 8} {:x -23.19078196871628, :y 8} {:x -23.154224200751045, :y 7} {:x -23.154224200751045, :y 7} {:x -23.154224200751045, :y 7} {:x -23.154224200751045, :y 7} {:x -14.167399480063514, :y 4} {:x -14.167399480063514, :y 4} {:x -14.167399480063514, :y 4} {:x -14.167399480063514, :y 4} {:x -19.85861158060795, :y 6} {:x -19.85861158060795, :y 6} {:x -19.85861158060795, :y 6} {:x -19.85861158060795, :y 6} {:x -17.124560966399606, :y 5} {:x -17.124560966399606, :y 5} {:x -17.124560966399606, :y 5} {:x -17.124560966399606, :y 5})}), :marks ({:type \"symbol\", :from {:data \"5483196b-791a-4411-b3ea-019c17e07a57\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"red\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"0eadf1e0-e300-4fe2-a063-f3e8f363476f\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"blue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"67fc3182-7f4b-41ec-b577-0f112a8e9734\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"#ff29d2\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}})}}"}
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
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":[0,4705]}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e","values":[{"x":0,"y":3031},{"x":1,"y":2637},{"x":2,"y":2838},{"x":3,"y":1869},{"x":4,"y":1351},{"x":5,"y":1734},{"x":6,"y":2059},{"x":7,"y":1612},{"x":8,"y":2100},{"x":9,"y":2165},{"x":10,"y":1197},{"x":11,"y":407},{"x":12,"y":959},{"x":13,"y":913},{"x":14,"y":848},{"x":15,"y":1009},{"x":16,"y":987},{"x":17,"y":890},{"x":18,"y":946},{"x":19,"y":966},{"x":20,"y":1400},{"x":21,"y":1718},{"x":22,"y":2031},{"x":23,"y":1650},{"x":24,"y":1348},{"x":25,"y":1258},{"x":26,"y":1490},{"x":27,"y":1899},{"x":28,"y":1674},{"x":29,"y":1100},{"x":30,"y":1399},{"x":31,"y":1214},{"x":32,"y":1460},{"x":33,"y":1084},{"x":34,"y":1215},{"x":35,"y":1594},{"x":36,"y":1194},{"x":37,"y":1239},{"x":38,"y":1771},{"x":39,"y":1526},{"x":40,"y":1051},{"x":41,"y":2470},{"x":42,"y":1664},{"x":43,"y":1225},{"x":44,"y":1762},{"x":45,"y":1188},{"x":46,"y":1796},{"x":47,"y":1497},{"x":48,"y":1914},{"x":49,"y":2222},{"x":50,"y":1302},{"x":51,"y":1136},{"x":52,"y":1460},{"x":53,"y":1611},{"x":54,"y":1331},{"x":55,"y":1318},{"x":56,"y":2424},{"x":57,"y":1350},{"x":58,"y":1283},{"x":59,"y":1921},{"x":60,"y":1702},{"x":61,"y":1019},{"x":62,"y":821},{"x":63,"y":1310},{"x":64,"y":1802},{"x":65,"y":2066},{"x":66,"y":1509},{"x":67,"y":2287},{"x":68,"y":2140},{"x":69,"y":1937},{"x":70,"y":1875},{"x":71,"y":1548},{"x":72,"y":1601},{"x":73,"y":2026},{"x":74,"y":1379},{"x":75,"y":1603},{"x":76,"y":2005},{"x":77,"y":1690},{"x":78,"y":1673},{"x":79,"y":1684},{"x":80,"y":1384},{"x":81,"y":1592},{"x":82,"y":1289},{"x":83,"y":920},{"x":84,"y":2628},{"x":85,"y":2592},{"x":86,"y":373},{"x":87,"y":545},{"x":88,"y":743},{"x":89,"y":1387},{"x":90,"y":744},{"x":91,"y":972},{"x":92,"y":1667},{"x":93,"y":2285},{"x":94,"y":1493},{"x":95,"y":1773},{"x":96,"y":1652},{"x":97,"y":821},{"x":98,"y":1127},{"x":99,"y":1863},{"x":100,"y":1325},{"x":101,"y":2153},{"x":102,"y":2383},{"x":103,"y":1900},{"x":104,"y":2430},{"x":105,"y":1021},{"x":106,"y":1964},{"x":107,"y":1894},{"x":108,"y":1609},{"x":109,"y":2735},{"x":110,"y":2225},{"x":111,"y":1434},{"x":112,"y":1407},{"x":113,"y":1866},{"x":114,"y":1911},{"x":115,"y":2302},{"x":116,"y":1728},{"x":117,"y":2156},{"x":118,"y":1484},{"x":119,"y":690},{"x":120,"y":1785},{"x":121,"y":1222},{"x":122,"y":1092},{"x":123,"y":881},{"x":124,"y":1377},{"x":125,"y":1120},{"x":126,"y":1207},{"x":127,"y":2382},{"x":128,"y":2127},{"x":129,"y":1442},{"x":130,"y":1421},{"x":131,"y":1703},{"x":132,"y":1918},{"x":133,"y":1344},{"x":134,"y":1465},{"x":135,"y":1353},{"x":136,"y":1304},{"x":137,"y":1654},{"x":138,"y":1873},{"x":139,"y":909},{"x":140,"y":1602},{"x":141,"y":1294},{"x":142,"y":1312},{"x":143,"y":1396},{"x":144,"y":1093},{"x":145,"y":1261},{"x":146,"y":1510},{"x":147,"y":1419},{"x":148,"y":1349},{"x":149,"y":1550},{"x":150,"y":1455},{"x":151,"y":1088},{"x":152,"y":759},{"x":153,"y":1376},{"x":154,"y":1554},{"x":155,"y":2036},{"x":156,"y":2141},{"x":157,"y":1327},{"x":158,"y":1599},{"x":159,"y":1744},{"x":160,"y":1822},{"x":161,"y":1973},{"x":162,"y":891},{"x":163,"y":1362},{"x":164,"y":1474},{"x":165,"y":1174},{"x":166,"y":1042},{"x":167,"y":992},{"x":168,"y":1336},{"x":169,"y":1381},{"x":170,"y":1492},{"x":171,"y":1796},{"x":172,"y":2131},{"x":173,"y":2566},{"x":174,"y":3189},{"x":175,"y":3030},{"x":176,"y":2253},{"x":177,"y":2453},{"x":178,"y":4705},{"x":179,"y":2415},{"x":180,"y":2477},{"x":181,"y":908},{"x":182,"y":1591},{"x":183,"y":3168},{"x":184,"y":1975},{"x":185,"y":1733},{"x":186,"y":1468},{"x":187,"y":2343},{"x":188,"y":4442},{"x":189,"y":2022},{"x":190,"y":1692},{"x":191,"y":2372},{"x":192,"y":1459},{"x":193,"y":1769},{"x":194,"y":1731},{"x":195,"y":1422},{"x":196,"y":2052},{"x":197,"y":1941},{"x":198,"y":2136},{"x":199,"y":1487}]},{"name":"d76170a5-22d3-412a-800a-608a1e3c2816","values":[{"x":0,"y":86},{"x":1,"y":147},{"x":2,"y":71},{"x":3,"y":78},{"x":4,"y":54},{"x":5,"y":53},{"x":6,"y":55},{"x":7,"y":49},{"x":8,"y":61},{"x":9,"y":49},{"x":10,"y":58},{"x":11,"y":56},{"x":12,"y":298},{"x":13,"y":365},{"x":14,"y":253},{"x":15,"y":305},{"x":16,"y":296},{"x":17,"y":290},{"x":18,"y":303},{"x":19,"y":243},{"x":20,"y":299},{"x":21,"y":207},{"x":22,"y":245},{"x":23,"y":228},{"x":24,"y":181},{"x":25,"y":197},{"x":26,"y":255},{"x":27,"y":193},{"x":28,"y":202},{"x":29,"y":222},{"x":30,"y":168},{"x":31,"y":188},{"x":32,"y":202},{"x":33,"y":162},{"x":34,"y":271},{"x":35,"y":260},{"x":36,"y":189},{"x":37,"y":157},{"x":38,"y":185},{"x":39,"y":212},{"x":40,"y":198},{"x":41,"y":240},{"x":42,"y":197},{"x":43,"y":164},{"x":44,"y":196},{"x":45,"y":215},{"x":46,"y":241},{"x":47,"y":192},{"x":48,"y":154},{"x":49,"y":173},{"x":50,"y":153},{"x":51,"y":240},{"x":52,"y":219},{"x":53,"y":297},{"x":54,"y":235},{"x":55,"y":180},{"x":56,"y":232},{"x":57,"y":120},{"x":58,"y":172},{"x":59,"y":233},{"x":60,"y":280},{"x":61,"y":246},{"x":62,"y":232},{"x":63,"y":193},{"x":64,"y":222},{"x":65,"y":227},{"x":66,"y":205},{"x":67,"y":171},{"x":68,"y":213},{"x":69,"y":206},{"x":70,"y":149},{"x":71,"y":205},{"x":72,"y":192},{"x":73,"y":234},{"x":74,"y":127},{"x":75,"y":228},{"x":76,"y":208},{"x":77,"y":125},{"x":78,"y":208},{"x":79,"y":183},{"x":80,"y":198},{"x":81,"y":125},{"x":82,"y":186},{"x":83,"y":229},{"x":84,"y":262},{"x":85,"y":206},{"x":86,"y":107},{"x":87,"y":267},{"x":88,"y":316},{"x":89,"y":258},{"x":90,"y":163},{"x":91,"y":218},{"x":92,"y":242},{"x":93,"y":275},{"x":94,"y":209},{"x":95,"y":217},{"x":96,"y":202},{"x":97,"y":237},{"x":98,"y":284},{"x":99,"y":193},{"x":100,"y":218},{"x":101,"y":168},{"x":102,"y":212},{"x":103,"y":128},{"x":104,"y":308},{"x":105,"y":237},{"x":106,"y":315},{"x":107,"y":258},{"x":108,"y":222},{"x":109,"y":229},{"x":110,"y":167},{"x":111,"y":159},{"x":112,"y":272},{"x":113,"y":309},{"x":114,"y":244},{"x":115,"y":207},{"x":116,"y":296},{"x":117,"y":327},{"x":118,"y":267},{"x":119,"y":259},{"x":120,"y":291},{"x":121,"y":286},{"x":122,"y":296},{"x":123,"y":380},{"x":124,"y":421},{"x":125,"y":320},{"x":126,"y":310},{"x":127,"y":271},{"x":128,"y":274},{"x":129,"y":201},{"x":130,"y":221},{"x":131,"y":208},{"x":132,"y":154},{"x":133,"y":168},{"x":134,"y":198},{"x":135,"y":154},{"x":136,"y":215},{"x":137,"y":273},{"x":138,"y":180},{"x":139,"y":119},{"x":140,"y":208},{"x":141,"y":212},{"x":142,"y":201},{"x":143,"y":206},{"x":144,"y":171},{"x":145,"y":204},{"x":146,"y":172},{"x":147,"y":206},{"x":148,"y":162},{"x":149,"y":140},{"x":150,"y":108},{"x":151,"y":167},{"x":152,"y":183},{"x":153,"y":234},{"x":154,"y":230},{"x":155,"y":179},{"x":156,"y":143},{"x":157,"y":151},{"x":158,"y":201},{"x":159,"y":166},{"x":160,"y":216},{"x":161,"y":236},{"x":162,"y":130},{"x":163,"y":215},{"x":164,"y":214},{"x":165,"y":164},{"x":166,"y":154},{"x":167,"y":206},{"x":168,"y":250},{"x":169,"y":227},{"x":170,"y":131},{"x":171,"y":163},{"x":172,"y":154},{"x":173,"y":250},{"x":174,"y":205},{"x":175,"y":217},{"x":176,"y":223},{"x":177,"y":227},{"x":178,"y":529},{"x":179,"y":206},{"x":180,"y":224},{"x":181,"y":243},{"x":182,"y":305},{"x":183,"y":214},{"x":184,"y":246},{"x":185,"y":198},{"x":186,"y":175},{"x":187,"y":179},{"x":188,"y":232},{"x":189,"y":140},{"x":190,"y":221},{"x":191,"y":190},{"x":192,"y":164},{"x":193,"y":148},{"x":194,"y":157},{"x":195,"y":203},{"x":196,"y":251},{"x":197,"y":149},{"x":198,"y":131},{"x":199,"y":183}]},{"name":"a2d08ed6-ffbe-496f-b540-c6142be30f24","values":[{"x":0,"y":25},{"x":1,"y":19},{"x":2,"y":4},{"x":3,"y":6},{"x":4,"y":6},{"x":5,"y":4},{"x":6,"y":3},{"x":7,"y":3},{"x":8,"y":3},{"x":9,"y":3},{"x":10,"y":4},{"x":11,"y":2},{"x":12,"y":3},{"x":13,"y":3},{"x":14,"y":2},{"x":15,"y":2},{"x":16,"y":2},{"x":17,"y":2},{"x":18,"y":1},{"x":19,"y":2},{"x":20,"y":3},{"x":21,"y":1},{"x":22,"y":2},{"x":23,"y":1},{"x":24,"y":1},{"x":25,"y":1},{"x":26,"y":1},{"x":27,"y":2},{"x":28,"y":2},{"x":29,"y":1},{"x":30,"y":1},{"x":31,"y":1},{"x":32,"y":2},{"x":33,"y":1},{"x":34,"y":1},{"x":35,"y":1},{"x":36,"y":0},{"x":37,"y":1},{"x":38,"y":1},{"x":39,"y":1},{"x":40,"y":1},{"x":41,"y":1},{"x":42,"y":1},{"x":43,"y":1},{"x":44,"y":1},{"x":45,"y":0},{"x":46,"y":1},{"x":47,"y":1},{"x":48,"y":2},{"x":49,"y":1},{"x":50,"y":0},{"x":51,"y":1},{"x":52,"y":1},{"x":53,"y":1},{"x":54,"y":1},{"x":55,"y":1},{"x":56,"y":1},{"x":57,"y":2},{"x":58,"y":1},{"x":59,"y":1},{"x":60,"y":1},{"x":61,"y":1},{"x":62,"y":1},{"x":63,"y":1},{"x":64,"y":2},{"x":65,"y":1},{"x":66,"y":1},{"x":67,"y":1},{"x":68,"y":1},{"x":69,"y":1},{"x":70,"y":2},{"x":71,"y":1},{"x":72,"y":1},{"x":73,"y":1},{"x":74,"y":0},{"x":75,"y":1},{"x":76,"y":0},{"x":77,"y":1},{"x":78,"y":1},{"x":79,"y":1},{"x":80,"y":1},{"x":81,"y":1},{"x":82,"y":1},{"x":83,"y":0},{"x":84,"y":1},{"x":85,"y":1},{"x":86,"y":0},{"x":87,"y":1},{"x":88,"y":0},{"x":89,"y":1},{"x":90,"y":1},{"x":91,"y":1},{"x":92,"y":2},{"x":93,"y":1},{"x":94,"y":0},{"x":95,"y":1},{"x":96,"y":1},{"x":97,"y":1},{"x":98,"y":1},{"x":99,"y":1},{"x":100,"y":0},{"x":101,"y":1},{"x":102,"y":1},{"x":103,"y":1},{"x":104,"y":1},{"x":105,"y":1},{"x":106,"y":0},{"x":107,"y":2},{"x":108,"y":1},{"x":109,"y":1},{"x":110,"y":1},{"x":111,"y":2},{"x":112,"y":1},{"x":113,"y":1},{"x":114,"y":1},{"x":115,"y":1},{"x":116,"y":1},{"x":117,"y":1},{"x":118,"y":1},{"x":119,"y":1},{"x":120,"y":1},{"x":121,"y":1},{"x":122,"y":1},{"x":123,"y":1},{"x":124,"y":1},{"x":125,"y":1},{"x":126,"y":0},{"x":127,"y":1},{"x":128,"y":1},{"x":129,"y":1},{"x":130,"y":0},{"x":131,"y":1},{"x":132,"y":1},{"x":133,"y":1},{"x":134,"y":1},{"x":135,"y":1},{"x":136,"y":0},{"x":137,"y":1},{"x":138,"y":0},{"x":139,"y":0},{"x":140,"y":0},{"x":141,"y":1},{"x":142,"y":1},{"x":143,"y":2},{"x":144,"y":1},{"x":145,"y":0},{"x":146,"y":1},{"x":147,"y":1},{"x":148,"y":1},{"x":149,"y":1},{"x":150,"y":1},{"x":151,"y":0},{"x":152,"y":1},{"x":153,"y":1},{"x":154,"y":1},{"x":155,"y":0},{"x":156,"y":1},{"x":157,"y":1},{"x":158,"y":1},{"x":159,"y":1},{"x":160,"y":1},{"x":161,"y":1},{"x":162,"y":1},{"x":163,"y":1},{"x":164,"y":0},{"x":165,"y":1},{"x":166,"y":0},{"x":167,"y":0},{"x":168,"y":0},{"x":169,"y":0},{"x":170,"y":1},{"x":171,"y":1},{"x":172,"y":1},{"x":173,"y":1},{"x":174,"y":1},{"x":175,"y":1},{"x":176,"y":2},{"x":177,"y":1},{"x":178,"y":12},{"x":179,"y":1},{"x":180,"y":5},{"x":181,"y":1},{"x":182,"y":0},{"x":183,"y":2},{"x":184,"y":4},{"x":185,"y":1},{"x":186,"y":1},{"x":187,"y":2},{"x":188,"y":1},{"x":189,"y":1},{"x":190,"y":2},{"x":191,"y":1},{"x":192,"y":1},{"x":193,"y":1},{"x":194,"y":0},{"x":195,"y":1},{"x":196,"y":1},{"x":197,"y":1},{"x":198,"y":1},{"x":199,"y":1}]},{"name":"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1","values":[{"x":0,"y":2920},{"x":1,"y":2471},{"x":2,"y":2763},{"x":3,"y":1785},{"x":4,"y":1291},{"x":5,"y":1677},{"x":6,"y":2001},{"x":7,"y":1560},{"x":8,"y":2036},{"x":9,"y":2113},{"x":10,"y":1135},{"x":11,"y":349},{"x":12,"y":658},{"x":13,"y":545},{"x":14,"y":593},{"x":15,"y":702},{"x":16,"y":689},{"x":17,"y":598},{"x":18,"y":642},{"x":19,"y":721},{"x":20,"y":1098},{"x":21,"y":1510},{"x":22,"y":1784},{"x":23,"y":1421},{"x":24,"y":1166},{"x":25,"y":1060},{"x":26,"y":1234},{"x":27,"y":1704},{"x":28,"y":1470},{"x":29,"y":877},{"x":30,"y":1230},{"x":31,"y":1025},{"x":32,"y":1256},{"x":33,"y":921},{"x":34,"y":943},{"x":35,"y":1333},{"x":36,"y":1005},{"x":37,"y":1081},{"x":38,"y":1585},{"x":39,"y":1313},{"x":40,"y":852},{"x":41,"y":2229},{"x":42,"y":1466},{"x":43,"y":1060},{"x":44,"y":1565},{"x":45,"y":973},{"x":46,"y":1554},{"x":47,"y":1304},{"x":48,"y":1758},{"x":49,"y":2048},{"x":50,"y":1149},{"x":51,"y":895},{"x":52,"y":1240},{"x":53,"y":1313},{"x":54,"y":1095},{"x":55,"y":1137},{"x":56,"y":2191},{"x":57,"y":1228},{"x":58,"y":1110},{"x":59,"y":1687},{"x":60,"y":1421},{"x":61,"y":772},{"x":62,"y":588},{"x":63,"y":1116},{"x":64,"y":1578},{"x":65,"y":1838},{"x":66,"y":1303},{"x":67,"y":2115},{"x":68,"y":1926},{"x":69,"y":1730},{"x":70,"y":1724},{"x":71,"y":1342},{"x":72,"y":1408},{"x":73,"y":1791},{"x":74,"y":1252},{"x":75,"y":1374},{"x":76,"y":1797},{"x":77,"y":1564},{"x":78,"y":1464},{"x":79,"y":1500},{"x":80,"y":1185},{"x":81,"y":1466},{"x":82,"y":1102},{"x":83,"y":691},{"x":84,"y":2365},{"x":85,"y":2385},{"x":86,"y":266},{"x":87,"y":277},{"x":88,"y":427},{"x":89,"y":1128},{"x":90,"y":580},{"x":91,"y":753},{"x":92,"y":1423},{"x":93,"y":2009},{"x":94,"y":1284},{"x":95,"y":1555},{"x":96,"y":1449},{"x":97,"y":583},{"x":98,"y":842},{"x":99,"y":1669},{"x":100,"y":1107},{"x":101,"y":1984},{"x":102,"y":2170},{"x":103,"y":1771},{"x":104,"y":2121},{"x":105,"y":783},{"x":106,"y":1649},{"x":107,"y":1634},{"x":108,"y":1386},{"x":109,"y":2505},{"x":110,"y":2057},{"x":111,"y":1273},{"x":112,"y":1134},{"x":113,"y":1556},{"x":114,"y":1666},{"x":115,"y":2094},{"x":116,"y":1431},{"x":117,"y":1828},{"x":118,"y":1216},{"x":119,"y":430},{"x":120,"y":1493},{"x":121,"y":935},{"x":122,"y":795},{"x":123,"y":500},{"x":124,"y":955},{"x":125,"y":799},{"x":126,"y":897},{"x":127,"y":2110},{"x":128,"y":1852},{"x":129,"y":1240},{"x":130,"y":1200},{"x":131,"y":1494},{"x":132,"y":1763},{"x":133,"y":1175},{"x":134,"y":1266},{"x":135,"y":1198},{"x":136,"y":1089},{"x":137,"y":1380},{"x":138,"y":1693},{"x":139,"y":790},{"x":140,"y":1394},{"x":141,"y":1081},{"x":142,"y":1110},{"x":143,"y":1188},{"x":144,"y":921},{"x":145,"y":1057},{"x":146,"y":1337},{"x":147,"y":1212},{"x":148,"y":1186},{"x":149,"y":1409},{"x":150,"y":1346},{"x":151,"y":921},{"x":152,"y":575},{"x":153,"y":1141},{"x":154,"y":1323},{"x":155,"y":1857},{"x":156,"y":1997},{"x":157,"y":1175},{"x":158,"y":1397},{"x":159,"y":1577},{"x":160,"y":1605},{"x":161,"y":1736},{"x":162,"y":760},{"x":163,"y":1146},{"x":164,"y":1260},{"x":165,"y":1009},{"x":166,"y":888},{"x":167,"y":786},{"x":168,"y":1086},{"x":169,"y":1154},{"x":170,"y":1360},{"x":171,"y":1632},{"x":172,"y":1976},{"x":173,"y":2315},{"x":174,"y":2983},{"x":175,"y":2812},{"x":176,"y":2028},{"x":177,"y":2225},{"x":178,"y":4164},{"x":179,"y":2208},{"x":180,"y":2248},{"x":181,"y":664},{"x":182,"y":1286},{"x":183,"y":2952},{"x":184,"y":1725},{"x":185,"y":1534},{"x":186,"y":1292},{"x":187,"y":2162},{"x":188,"y":4209},{"x":189,"y":1881},{"x":190,"y":1469},{"x":191,"y":2181},{"x":192,"y":1294},{"x":193,"y":1620},{"x":194,"y":1574},{"x":195,"y":1218},{"x":196,"y":1800},{"x":197,"y":1791},{"x":198,"y":2004},{"x":199,"y":1303}]}],"marks":[{"type":"line","from":{"data":"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"yellow"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"d76170a5-22d3-412a-800a-608a1e3c2816"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"red"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"a2d08ed6-ffbe-496f-b540-c6142be30f24"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"green"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"blue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain [0 4705]}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e\", :values ({:x 0, :y 3031} {:x 1, :y 2637} {:x 2, :y 2838} {:x 3, :y 1869} {:x 4, :y 1351} {:x 5, :y 1734} {:x 6, :y 2059} {:x 7, :y 1612} {:x 8, :y 2100} {:x 9, :y 2165} {:x 10, :y 1197} {:x 11, :y 407} {:x 12, :y 959} {:x 13, :y 913} {:x 14, :y 848} {:x 15, :y 1009} {:x 16, :y 987} {:x 17, :y 890} {:x 18, :y 946} {:x 19, :y 966} {:x 20, :y 1400} {:x 21, :y 1718} {:x 22, :y 2031} {:x 23, :y 1650} {:x 24, :y 1348} {:x 25, :y 1258} {:x 26, :y 1490} {:x 27, :y 1899} {:x 28, :y 1674} {:x 29, :y 1100} {:x 30, :y 1399} {:x 31, :y 1214} {:x 32, :y 1460} {:x 33, :y 1084} {:x 34, :y 1215} {:x 35, :y 1594} {:x 36, :y 1194} {:x 37, :y 1239} {:x 38, :y 1771} {:x 39, :y 1526} {:x 40, :y 1051} {:x 41, :y 2470} {:x 42, :y 1664} {:x 43, :y 1225} {:x 44, :y 1762} {:x 45, :y 1188} {:x 46, :y 1796} {:x 47, :y 1497} {:x 48, :y 1914} {:x 49, :y 2222} {:x 50, :y 1302} {:x 51, :y 1136} {:x 52, :y 1460} {:x 53, :y 1611} {:x 54, :y 1331} {:x 55, :y 1318} {:x 56, :y 2424} {:x 57, :y 1350} {:x 58, :y 1283} {:x 59, :y 1921} {:x 60, :y 1702} {:x 61, :y 1019} {:x 62, :y 821} {:x 63, :y 1310} {:x 64, :y 1802} {:x 65, :y 2066} {:x 66, :y 1509} {:x 67, :y 2287} {:x 68, :y 2140} {:x 69, :y 1937} {:x 70, :y 1875} {:x 71, :y 1548} {:x 72, :y 1601} {:x 73, :y 2026} {:x 74, :y 1379} {:x 75, :y 1603} {:x 76, :y 2005} {:x 77, :y 1690} {:x 78, :y 1673} {:x 79, :y 1684} {:x 80, :y 1384} {:x 81, :y 1592} {:x 82, :y 1289} {:x 83, :y 920} {:x 84, :y 2628} {:x 85, :y 2592} {:x 86, :y 373} {:x 87, :y 545} {:x 88, :y 743} {:x 89, :y 1387} {:x 90, :y 744} {:x 91, :y 972} {:x 92, :y 1667} {:x 93, :y 2285} {:x 94, :y 1493} {:x 95, :y 1773} {:x 96, :y 1652} {:x 97, :y 821} {:x 98, :y 1127} {:x 99, :y 1863} {:x 100, :y 1325} {:x 101, :y 2153} {:x 102, :y 2383} {:x 103, :y 1900} {:x 104, :y 2430} {:x 105, :y 1021} {:x 106, :y 1964} {:x 107, :y 1894} {:x 108, :y 1609} {:x 109, :y 2735} {:x 110, :y 2225} {:x 111, :y 1434} {:x 112, :y 1407} {:x 113, :y 1866} {:x 114, :y 1911} {:x 115, :y 2302} {:x 116, :y 1728} {:x 117, :y 2156} {:x 118, :y 1484} {:x 119, :y 690} {:x 120, :y 1785} {:x 121, :y 1222} {:x 122, :y 1092} {:x 123, :y 881} {:x 124, :y 1377} {:x 125, :y 1120} {:x 126, :y 1207} {:x 127, :y 2382} {:x 128, :y 2127} {:x 129, :y 1442} {:x 130, :y 1421} {:x 131, :y 1703} {:x 132, :y 1918} {:x 133, :y 1344} {:x 134, :y 1465} {:x 135, :y 1353} {:x 136, :y 1304} {:x 137, :y 1654} {:x 138, :y 1873} {:x 139, :y 909} {:x 140, :y 1602} {:x 141, :y 1294} {:x 142, :y 1312} {:x 143, :y 1396} {:x 144, :y 1093} {:x 145, :y 1261} {:x 146, :y 1510} {:x 147, :y 1419} {:x 148, :y 1349} {:x 149, :y 1550} {:x 150, :y 1455} {:x 151, :y 1088} {:x 152, :y 759} {:x 153, :y 1376} {:x 154, :y 1554} {:x 155, :y 2036} {:x 156, :y 2141} {:x 157, :y 1327} {:x 158, :y 1599} {:x 159, :y 1744} {:x 160, :y 1822} {:x 161, :y 1973} {:x 162, :y 891} {:x 163, :y 1362} {:x 164, :y 1474} {:x 165, :y 1174} {:x 166, :y 1042} {:x 167, :y 992} {:x 168, :y 1336} {:x 169, :y 1381} {:x 170, :y 1492} {:x 171, :y 1796} {:x 172, :y 2131} {:x 173, :y 2566} {:x 174, :y 3189} {:x 175, :y 3030} {:x 176, :y 2253} {:x 177, :y 2453} {:x 178, :y 4705} {:x 179, :y 2415} {:x 180, :y 2477} {:x 181, :y 908} {:x 182, :y 1591} {:x 183, :y 3168} {:x 184, :y 1975} {:x 185, :y 1733} {:x 186, :y 1468} {:x 187, :y 2343} {:x 188, :y 4442} {:x 189, :y 2022} {:x 190, :y 1692} {:x 191, :y 2372} {:x 192, :y 1459} {:x 193, :y 1769} {:x 194, :y 1731} {:x 195, :y 1422} {:x 196, :y 2052} {:x 197, :y 1941} {:x 198, :y 2136} {:x 199, :y 1487})} {:name \"d76170a5-22d3-412a-800a-608a1e3c2816\", :values ({:x 0, :y 86} {:x 1, :y 147} {:x 2, :y 71} {:x 3, :y 78} {:x 4, :y 54} {:x 5, :y 53} {:x 6, :y 55} {:x 7, :y 49} {:x 8, :y 61} {:x 9, :y 49} {:x 10, :y 58} {:x 11, :y 56} {:x 12, :y 298} {:x 13, :y 365} {:x 14, :y 253} {:x 15, :y 305} {:x 16, :y 296} {:x 17, :y 290} {:x 18, :y 303} {:x 19, :y 243} {:x 20, :y 299} {:x 21, :y 207} {:x 22, :y 245} {:x 23, :y 228} {:x 24, :y 181} {:x 25, :y 197} {:x 26, :y 255} {:x 27, :y 193} {:x 28, :y 202} {:x 29, :y 222} {:x 30, :y 168} {:x 31, :y 188} {:x 32, :y 202} {:x 33, :y 162} {:x 34, :y 271} {:x 35, :y 260} {:x 36, :y 189} {:x 37, :y 157} {:x 38, :y 185} {:x 39, :y 212} {:x 40, :y 198} {:x 41, :y 240} {:x 42, :y 197} {:x 43, :y 164} {:x 44, :y 196} {:x 45, :y 215} {:x 46, :y 241} {:x 47, :y 192} {:x 48, :y 154} {:x 49, :y 173} {:x 50, :y 153} {:x 51, :y 240} {:x 52, :y 219} {:x 53, :y 297} {:x 54, :y 235} {:x 55, :y 180} {:x 56, :y 232} {:x 57, :y 120} {:x 58, :y 172} {:x 59, :y 233} {:x 60, :y 280} {:x 61, :y 246} {:x 62, :y 232} {:x 63, :y 193} {:x 64, :y 222} {:x 65, :y 227} {:x 66, :y 205} {:x 67, :y 171} {:x 68, :y 213} {:x 69, :y 206} {:x 70, :y 149} {:x 71, :y 205} {:x 72, :y 192} {:x 73, :y 234} {:x 74, :y 127} {:x 75, :y 228} {:x 76, :y 208} {:x 77, :y 125} {:x 78, :y 208} {:x 79, :y 183} {:x 80, :y 198} {:x 81, :y 125} {:x 82, :y 186} {:x 83, :y 229} {:x 84, :y 262} {:x 85, :y 206} {:x 86, :y 107} {:x 87, :y 267} {:x 88, :y 316} {:x 89, :y 258} {:x 90, :y 163} {:x 91, :y 218} {:x 92, :y 242} {:x 93, :y 275} {:x 94, :y 209} {:x 95, :y 217} {:x 96, :y 202} {:x 97, :y 237} {:x 98, :y 284} {:x 99, :y 193} {:x 100, :y 218} {:x 101, :y 168} {:x 102, :y 212} {:x 103, :y 128} {:x 104, :y 308} {:x 105, :y 237} {:x 106, :y 315} {:x 107, :y 258} {:x 108, :y 222} {:x 109, :y 229} {:x 110, :y 167} {:x 111, :y 159} {:x 112, :y 272} {:x 113, :y 309} {:x 114, :y 244} {:x 115, :y 207} {:x 116, :y 296} {:x 117, :y 327} {:x 118, :y 267} {:x 119, :y 259} {:x 120, :y 291} {:x 121, :y 286} {:x 122, :y 296} {:x 123, :y 380} {:x 124, :y 421} {:x 125, :y 320} {:x 126, :y 310} {:x 127, :y 271} {:x 128, :y 274} {:x 129, :y 201} {:x 130, :y 221} {:x 131, :y 208} {:x 132, :y 154} {:x 133, :y 168} {:x 134, :y 198} {:x 135, :y 154} {:x 136, :y 215} {:x 137, :y 273} {:x 138, :y 180} {:x 139, :y 119} {:x 140, :y 208} {:x 141, :y 212} {:x 142, :y 201} {:x 143, :y 206} {:x 144, :y 171} {:x 145, :y 204} {:x 146, :y 172} {:x 147, :y 206} {:x 148, :y 162} {:x 149, :y 140} {:x 150, :y 108} {:x 151, :y 167} {:x 152, :y 183} {:x 153, :y 234} {:x 154, :y 230} {:x 155, :y 179} {:x 156, :y 143} {:x 157, :y 151} {:x 158, :y 201} {:x 159, :y 166} {:x 160, :y 216} {:x 161, :y 236} {:x 162, :y 130} {:x 163, :y 215} {:x 164, :y 214} {:x 165, :y 164} {:x 166, :y 154} {:x 167, :y 206} {:x 168, :y 250} {:x 169, :y 227} {:x 170, :y 131} {:x 171, :y 163} {:x 172, :y 154} {:x 173, :y 250} {:x 174, :y 205} {:x 175, :y 217} {:x 176, :y 223} {:x 177, :y 227} {:x 178, :y 529} {:x 179, :y 206} {:x 180, :y 224} {:x 181, :y 243} {:x 182, :y 305} {:x 183, :y 214} {:x 184, :y 246} {:x 185, :y 198} {:x 186, :y 175} {:x 187, :y 179} {:x 188, :y 232} {:x 189, :y 140} {:x 190, :y 221} {:x 191, :y 190} {:x 192, :y 164} {:x 193, :y 148} {:x 194, :y 157} {:x 195, :y 203} {:x 196, :y 251} {:x 197, :y 149} {:x 198, :y 131} {:x 199, :y 183})} {:name \"a2d08ed6-ffbe-496f-b540-c6142be30f24\", :values ({:x 0, :y 25} {:x 1, :y 19} {:x 2, :y 4} {:x 3, :y 6} {:x 4, :y 6} {:x 5, :y 4} {:x 6, :y 3} {:x 7, :y 3} {:x 8, :y 3} {:x 9, :y 3} {:x 10, :y 4} {:x 11, :y 2} {:x 12, :y 3} {:x 13, :y 3} {:x 14, :y 2} {:x 15, :y 2} {:x 16, :y 2} {:x 17, :y 2} {:x 18, :y 1} {:x 19, :y 2} {:x 20, :y 3} {:x 21, :y 1} {:x 22, :y 2} {:x 23, :y 1} {:x 24, :y 1} {:x 25, :y 1} {:x 26, :y 1} {:x 27, :y 2} {:x 28, :y 2} {:x 29, :y 1} {:x 30, :y 1} {:x 31, :y 1} {:x 32, :y 2} {:x 33, :y 1} {:x 34, :y 1} {:x 35, :y 1} {:x 36, :y 0} {:x 37, :y 1} {:x 38, :y 1} {:x 39, :y 1} {:x 40, :y 1} {:x 41, :y 1} {:x 42, :y 1} {:x 43, :y 1} {:x 44, :y 1} {:x 45, :y 0} {:x 46, :y 1} {:x 47, :y 1} {:x 48, :y 2} {:x 49, :y 1} {:x 50, :y 0} {:x 51, :y 1} {:x 52, :y 1} {:x 53, :y 1} {:x 54, :y 1} {:x 55, :y 1} {:x 56, :y 1} {:x 57, :y 2} {:x 58, :y 1} {:x 59, :y 1} {:x 60, :y 1} {:x 61, :y 1} {:x 62, :y 1} {:x 63, :y 1} {:x 64, :y 2} {:x 65, :y 1} {:x 66, :y 1} {:x 67, :y 1} {:x 68, :y 1} {:x 69, :y 1} {:x 70, :y 2} {:x 71, :y 1} {:x 72, :y 1} {:x 73, :y 1} {:x 74, :y 0} {:x 75, :y 1} {:x 76, :y 0} {:x 77, :y 1} {:x 78, :y 1} {:x 79, :y 1} {:x 80, :y 1} {:x 81, :y 1} {:x 82, :y 1} {:x 83, :y 0} {:x 84, :y 1} {:x 85, :y 1} {:x 86, :y 0} {:x 87, :y 1} {:x 88, :y 0} {:x 89, :y 1} {:x 90, :y 1} {:x 91, :y 1} {:x 92, :y 2} {:x 93, :y 1} {:x 94, :y 0} {:x 95, :y 1} {:x 96, :y 1} {:x 97, :y 1} {:x 98, :y 1} {:x 99, :y 1} {:x 100, :y 0} {:x 101, :y 1} {:x 102, :y 1} {:x 103, :y 1} {:x 104, :y 1} {:x 105, :y 1} {:x 106, :y 0} {:x 107, :y 2} {:x 108, :y 1} {:x 109, :y 1} {:x 110, :y 1} {:x 111, :y 2} {:x 112, :y 1} {:x 113, :y 1} {:x 114, :y 1} {:x 115, :y 1} {:x 116, :y 1} {:x 117, :y 1} {:x 118, :y 1} {:x 119, :y 1} {:x 120, :y 1} {:x 121, :y 1} {:x 122, :y 1} {:x 123, :y 1} {:x 124, :y 1} {:x 125, :y 1} {:x 126, :y 0} {:x 127, :y 1} {:x 128, :y 1} {:x 129, :y 1} {:x 130, :y 0} {:x 131, :y 1} {:x 132, :y 1} {:x 133, :y 1} {:x 134, :y 1} {:x 135, :y 1} {:x 136, :y 0} {:x 137, :y 1} {:x 138, :y 0} {:x 139, :y 0} {:x 140, :y 0} {:x 141, :y 1} {:x 142, :y 1} {:x 143, :y 2} {:x 144, :y 1} {:x 145, :y 0} {:x 146, :y 1} {:x 147, :y 1} {:x 148, :y 1} {:x 149, :y 1} {:x 150, :y 1} {:x 151, :y 0} {:x 152, :y 1} {:x 153, :y 1} {:x 154, :y 1} {:x 155, :y 0} {:x 156, :y 1} {:x 157, :y 1} {:x 158, :y 1} {:x 159, :y 1} {:x 160, :y 1} {:x 161, :y 1} {:x 162, :y 1} {:x 163, :y 1} {:x 164, :y 0} {:x 165, :y 1} {:x 166, :y 0} {:x 167, :y 0} {:x 168, :y 0} {:x 169, :y 0} {:x 170, :y 1} {:x 171, :y 1} {:x 172, :y 1} {:x 173, :y 1} {:x 174, :y 1} {:x 175, :y 1} {:x 176, :y 2} {:x 177, :y 1} {:x 178, :y 12} {:x 179, :y 1} {:x 180, :y 5} {:x 181, :y 1} {:x 182, :y 0} {:x 183, :y 2} {:x 184, :y 4} {:x 185, :y 1} {:x 186, :y 1} {:x 187, :y 2} {:x 188, :y 1} {:x 189, :y 1} {:x 190, :y 2} {:x 191, :y 1} {:x 192, :y 1} {:x 193, :y 1} {:x 194, :y 0} {:x 195, :y 1} {:x 196, :y 1} {:x 197, :y 1} {:x 198, :y 1} {:x 199, :y 1})} {:name \"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1\", :values ({:x 0, :y 2920} {:x 1, :y 2471} {:x 2, :y 2763} {:x 3, :y 1785} {:x 4, :y 1291} {:x 5, :y 1677} {:x 6, :y 2001} {:x 7, :y 1560} {:x 8, :y 2036} {:x 9, :y 2113} {:x 10, :y 1135} {:x 11, :y 349} {:x 12, :y 658} {:x 13, :y 545} {:x 14, :y 593} {:x 15, :y 702} {:x 16, :y 689} {:x 17, :y 598} {:x 18, :y 642} {:x 19, :y 721} {:x 20, :y 1098} {:x 21, :y 1510} {:x 22, :y 1784} {:x 23, :y 1421} {:x 24, :y 1166} {:x 25, :y 1060} {:x 26, :y 1234} {:x 27, :y 1704} {:x 28, :y 1470} {:x 29, :y 877} {:x 30, :y 1230} {:x 31, :y 1025} {:x 32, :y 1256} {:x 33, :y 921} {:x 34, :y 943} {:x 35, :y 1333} {:x 36, :y 1005} {:x 37, :y 1081} {:x 38, :y 1585} {:x 39, :y 1313} {:x 40, :y 852} {:x 41, :y 2229} {:x 42, :y 1466} {:x 43, :y 1060} {:x 44, :y 1565} {:x 45, :y 973} {:x 46, :y 1554} {:x 47, :y 1304} {:x 48, :y 1758} {:x 49, :y 2048} {:x 50, :y 1149} {:x 51, :y 895} {:x 52, :y 1240} {:x 53, :y 1313} {:x 54, :y 1095} {:x 55, :y 1137} {:x 56, :y 2191} {:x 57, :y 1228} {:x 58, :y 1110} {:x 59, :y 1687} {:x 60, :y 1421} {:x 61, :y 772} {:x 62, :y 588} {:x 63, :y 1116} {:x 64, :y 1578} {:x 65, :y 1838} {:x 66, :y 1303} {:x 67, :y 2115} {:x 68, :y 1926} {:x 69, :y 1730} {:x 70, :y 1724} {:x 71, :y 1342} {:x 72, :y 1408} {:x 73, :y 1791} {:x 74, :y 1252} {:x 75, :y 1374} {:x 76, :y 1797} {:x 77, :y 1564} {:x 78, :y 1464} {:x 79, :y 1500} {:x 80, :y 1185} {:x 81, :y 1466} {:x 82, :y 1102} {:x 83, :y 691} {:x 84, :y 2365} {:x 85, :y 2385} {:x 86, :y 266} {:x 87, :y 277} {:x 88, :y 427} {:x 89, :y 1128} {:x 90, :y 580} {:x 91, :y 753} {:x 92, :y 1423} {:x 93, :y 2009} {:x 94, :y 1284} {:x 95, :y 1555} {:x 96, :y 1449} {:x 97, :y 583} {:x 98, :y 842} {:x 99, :y 1669} {:x 100, :y 1107} {:x 101, :y 1984} {:x 102, :y 2170} {:x 103, :y 1771} {:x 104, :y 2121} {:x 105, :y 783} {:x 106, :y 1649} {:x 107, :y 1634} {:x 108, :y 1386} {:x 109, :y 2505} {:x 110, :y 2057} {:x 111, :y 1273} {:x 112, :y 1134} {:x 113, :y 1556} {:x 114, :y 1666} {:x 115, :y 2094} {:x 116, :y 1431} {:x 117, :y 1828} {:x 118, :y 1216} {:x 119, :y 430} {:x 120, :y 1493} {:x 121, :y 935} {:x 122, :y 795} {:x 123, :y 500} {:x 124, :y 955} {:x 125, :y 799} {:x 126, :y 897} {:x 127, :y 2110} {:x 128, :y 1852} {:x 129, :y 1240} {:x 130, :y 1200} {:x 131, :y 1494} {:x 132, :y 1763} {:x 133, :y 1175} {:x 134, :y 1266} {:x 135, :y 1198} {:x 136, :y 1089} {:x 137, :y 1380} {:x 138, :y 1693} {:x 139, :y 790} {:x 140, :y 1394} {:x 141, :y 1081} {:x 142, :y 1110} {:x 143, :y 1188} {:x 144, :y 921} {:x 145, :y 1057} {:x 146, :y 1337} {:x 147, :y 1212} {:x 148, :y 1186} {:x 149, :y 1409} {:x 150, :y 1346} {:x 151, :y 921} {:x 152, :y 575} {:x 153, :y 1141} {:x 154, :y 1323} {:x 155, :y 1857} {:x 156, :y 1997} {:x 157, :y 1175} {:x 158, :y 1397} {:x 159, :y 1577} {:x 160, :y 1605} {:x 161, :y 1736} {:x 162, :y 760} {:x 163, :y 1146} {:x 164, :y 1260} {:x 165, :y 1009} {:x 166, :y 888} {:x 167, :y 786} {:x 168, :y 1086} {:x 169, :y 1154} {:x 170, :y 1360} {:x 171, :y 1632} {:x 172, :y 1976} {:x 173, :y 2315} {:x 174, :y 2983} {:x 175, :y 2812} {:x 176, :y 2028} {:x 177, :y 2225} {:x 178, :y 4164} {:x 179, :y 2208} {:x 180, :y 2248} {:x 181, :y 664} {:x 182, :y 1286} {:x 183, :y 2952} {:x 184, :y 1725} {:x 185, :y 1534} {:x 186, :y 1292} {:x 187, :y 2162} {:x 188, :y 4209} {:x 189, :y 1881} {:x 190, :y 1469} {:x 191, :y 2181} {:x 192, :y 1294} {:x 193, :y 1620} {:x 194, :y 1574} {:x 195, :y 1218} {:x 196, :y 1800} {:x 197, :y 1791} {:x 198, :y 2004} {:x 199, :y 1303})}), :marks ({:type \"line\", :from {:data \"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"yellow\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"d76170a5-22d3-412a-800a-608a1e3c2816\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"red\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"a2d08ed6-ffbe-496f-b540-c6142be30f24\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"green\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"blue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}})}}"}
;; <=

;; **
;;; Spread of complexity in population
;; **

;; @@
(plot/histogram (map :complexity (:rabble result)))
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"f50f1304-d1f3-41cb-8a7f-1600db887930","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"f50f1304-d1f3-41cb-8a7f-1600db887930","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"f50f1304-d1f3-41cb-8a7f-1600db887930"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"interpolate":{"value":"step-before"},"fill":{"value":"steelblue"},"fillOpacity":{"value":0.4},"stroke":{"value":"steelblue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"f50f1304-d1f3-41cb-8a7f-1600db887930","values":[{"x":1.0,"y":0},{"x":2.125,"y":44.0},{"x":3.25,"y":9.0},{"x":4.375,"y":14.0},{"x":5.5,"y":8.0},{"x":6.625,"y":9.0},{"x":7.75,"y":3.0},{"x":8.875,"y":8.0},{"x":10.0,"y":2.0},{"x":11.125,"y":3.0},{"x":12.25,"y":0}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"f50f1304-d1f3-41cb-8a7f-1600db887930\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"f50f1304-d1f3-41cb-8a7f-1600db887930\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"f50f1304-d1f3-41cb-8a7f-1600db887930\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :interpolate {:value \"step-before\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 0.4}, :stroke {:value \"steelblue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"f50f1304-d1f3-41cb-8a7f-1600db887930\", :values ({:x 1.0, :y 0} {:x 2.125, :y 44.0} {:x 3.25, :y 9.0} {:x 4.375, :y 14.0} {:x 5.5, :y 8.0} {:x 6.625, :y 9.0} {:x 7.75, :y 3.0} {:x 8.875, :y 8.0} {:x 10.0, :y 2.0} {:x 11.125, :y 3.0} {:x 12.25, :y 0})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(plot/list-plot (:mean (:complexity @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288","values":[{"x":0,"y":4.09},{"x":1,"y":3.395},{"x":2,"y":3.715},{"x":3,"y":3.215},{"x":4,"y":3.32},{"x":5,"y":2.79},{"x":6,"y":3.565},{"x":7,"y":2.825},{"x":8,"y":3.64},{"x":9,"y":3.77},{"x":10,"y":3.19},{"x":11,"y":2.08},{"x":12,"y":1.56},{"x":13,"y":1.565},{"x":14,"y":1.59},{"x":15,"y":1.785},{"x":16,"y":1.84},{"x":17,"y":1.955},{"x":18,"y":2.0},{"x":19,"y":1.935},{"x":20,"y":2.375},{"x":21,"y":2.655},{"x":22,"y":2.865},{"x":23,"y":2.81},{"x":24,"y":2.685},{"x":25,"y":2.535},{"x":26,"y":2.79},{"x":27,"y":2.845},{"x":28,"y":2.81},{"x":29,"y":2.795},{"x":30,"y":2.685},{"x":31,"y":2.675},{"x":32,"y":2.855},{"x":33,"y":2.66},{"x":34,"y":2.52},{"x":35,"y":2.76},{"x":36,"y":2.755},{"x":37,"y":2.6},{"x":38,"y":2.73},{"x":39,"y":2.795},{"x":40,"y":2.725},{"x":41,"y":3.165},{"x":42,"y":2.945},{"x":43,"y":2.785},{"x":44,"y":3.05},{"x":45,"y":2.87},{"x":46,"y":3.18},{"x":47,"y":3.12},{"x":48,"y":3.07},{"x":49,"y":3.365},{"x":50,"y":2.935},{"x":51,"y":2.625},{"x":52,"y":2.63},{"x":53,"y":2.755},{"x":54,"y":2.74},{"x":55,"y":2.725},{"x":56,"y":3.16},{"x":57,"y":2.975},{"x":58,"y":2.575},{"x":59,"y":3.015},{"x":60,"y":2.93},{"x":61,"y":2.765},{"x":62,"y":2.715},{"x":63,"y":2.99},{"x":64,"y":2.98},{"x":65,"y":3.07},{"x":66,"y":3.19},{"x":67,"y":3.195},{"x":68,"y":3.075},{"x":69,"y":3.23},{"x":70,"y":3.085},{"x":71,"y":3.085},{"x":72,"y":3.14},{"x":73,"y":3.495},{"x":74,"y":3.15},{"x":75,"y":3.29},{"x":76,"y":3.405},{"x":77,"y":3.145},{"x":78,"y":3.185},{"x":79,"y":3.33},{"x":80,"y":3.18},{"x":81,"y":2.85},{"x":82,"y":2.975},{"x":83,"y":2.69},{"x":84,"y":3.18},{"x":85,"y":3.375},{"x":86,"y":2.475},{"x":87,"y":2.405},{"x":88,"y":2.59},{"x":89,"y":2.625},{"x":90,"y":2.615},{"x":91,"y":2.66},{"x":92,"y":2.605},{"x":93,"y":3.06},{"x":94,"y":2.82},{"x":95,"y":2.845},{"x":96,"y":2.885},{"x":97,"y":2.785},{"x":98,"y":2.865},{"x":99,"y":3.155},{"x":100,"y":2.84},{"x":101,"y":2.985},{"x":102,"y":3.085},{"x":103,"y":2.375},{"x":104,"y":2.775},{"x":105,"y":2.47},{"x":106,"y":2.77},{"x":107,"y":2.76},{"x":108,"y":2.855},{"x":109,"y":3.065},{"x":110,"y":3.1},{"x":111,"y":2.645},{"x":112,"y":2.41},{"x":113,"y":2.565},{"x":114,"y":2.695},{"x":115,"y":2.6},{"x":116,"y":2.32},{"x":117,"y":2.355},{"x":118,"y":2.085},{"x":119,"y":1.985},{"x":120,"y":2.05},{"x":121,"y":1.955},{"x":122,"y":1.91},{"x":123,"y":1.815},{"x":124,"y":2.025},{"x":125,"y":1.97},{"x":126,"y":2.11},{"x":127,"y":2.435},{"x":128,"y":2.58},{"x":129,"y":2.425},{"x":130,"y":2.4},{"x":131,"y":2.63},{"x":132,"y":2.67},{"x":133,"y":2.63},{"x":134,"y":2.685},{"x":135,"y":2.735},{"x":136,"y":2.66},{"x":137,"y":2.785},{"x":138,"y":2.87},{"x":139,"y":2.345},{"x":140,"y":2.435},{"x":141,"y":2.405},{"x":142,"y":2.645},{"x":143,"y":2.855},{"x":144,"y":2.61},{"x":145,"y":2.775},{"x":146,"y":2.81},{"x":147,"y":2.805},{"x":148,"y":3.16},{"x":149,"y":3.18},{"x":150,"y":3.075},{"x":151,"y":2.915},{"x":152,"y":2.25},{"x":153,"y":2.56},{"x":154,"y":2.795},{"x":155,"y":2.835},{"x":156,"y":2.915},{"x":157,"y":2.795},{"x":158,"y":2.785},{"x":159,"y":2.785},{"x":160,"y":2.625},{"x":161,"y":3.095},{"x":162,"y":2.565},{"x":163,"y":2.655},{"x":164,"y":2.875},{"x":165,"y":2.985},{"x":166,"y":2.53},{"x":167,"y":2.52},{"x":168,"y":2.785},{"x":169,"y":3.005},{"x":170,"y":3.06},{"x":171,"y":3.265},{"x":172,"y":3.16},{"x":173,"y":3.135},{"x":174,"y":3.225},{"x":175,"y":3.125},{"x":176,"y":3.005},{"x":177,"y":3.095},{"x":178,"y":3.41},{"x":179,"y":3.12},{"x":180,"y":3.44},{"x":181,"y":3.0},{"x":182,"y":3.17},{"x":183,"y":3.335},{"x":184,"y":3.27},{"x":185,"y":3.48},{"x":186,"y":3.41},{"x":187,"y":3.425},{"x":188,"y":3.65},{"x":189,"y":3.11},{"x":190,"y":3.42},{"x":191,"y":3.465},{"x":192,"y":3.32},{"x":193,"y":3.4},{"x":194,"y":3.38},{"x":195,"y":3.125},{"x":196,"y":3.69},{"x":197,"y":3.455},{"x":198,"y":3.25},{"x":199,"y":3.285}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\", :values ({:x 0, :y 4.09} {:x 1, :y 3.395} {:x 2, :y 3.715} {:x 3, :y 3.215} {:x 4, :y 3.32} {:x 5, :y 2.79} {:x 6, :y 3.565} {:x 7, :y 2.825} {:x 8, :y 3.64} {:x 9, :y 3.77} {:x 10, :y 3.19} {:x 11, :y 2.08} {:x 12, :y 1.56} {:x 13, :y 1.565} {:x 14, :y 1.59} {:x 15, :y 1.785} {:x 16, :y 1.84} {:x 17, :y 1.955} {:x 18, :y 2.0} {:x 19, :y 1.935} {:x 20, :y 2.375} {:x 21, :y 2.655} {:x 22, :y 2.865} {:x 23, :y 2.81} {:x 24, :y 2.685} {:x 25, :y 2.535} {:x 26, :y 2.79} {:x 27, :y 2.845} {:x 28, :y 2.81} {:x 29, :y 2.795} {:x 30, :y 2.685} {:x 31, :y 2.675} {:x 32, :y 2.855} {:x 33, :y 2.66} {:x 34, :y 2.52} {:x 35, :y 2.76} {:x 36, :y 2.755} {:x 37, :y 2.6} {:x 38, :y 2.73} {:x 39, :y 2.795} {:x 40, :y 2.725} {:x 41, :y 3.165} {:x 42, :y 2.945} {:x 43, :y 2.785} {:x 44, :y 3.05} {:x 45, :y 2.87} {:x 46, :y 3.18} {:x 47, :y 3.12} {:x 48, :y 3.07} {:x 49, :y 3.365} {:x 50, :y 2.935} {:x 51, :y 2.625} {:x 52, :y 2.63} {:x 53, :y 2.755} {:x 54, :y 2.74} {:x 55, :y 2.725} {:x 56, :y 3.16} {:x 57, :y 2.975} {:x 58, :y 2.575} {:x 59, :y 3.015} {:x 60, :y 2.93} {:x 61, :y 2.765} {:x 62, :y 2.715} {:x 63, :y 2.99} {:x 64, :y 2.98} {:x 65, :y 3.07} {:x 66, :y 3.19} {:x 67, :y 3.195} {:x 68, :y 3.075} {:x 69, :y 3.23} {:x 70, :y 3.085} {:x 71, :y 3.085} {:x 72, :y 3.14} {:x 73, :y 3.495} {:x 74, :y 3.15} {:x 75, :y 3.29} {:x 76, :y 3.405} {:x 77, :y 3.145} {:x 78, :y 3.185} {:x 79, :y 3.33} {:x 80, :y 3.18} {:x 81, :y 2.85} {:x 82, :y 2.975} {:x 83, :y 2.69} {:x 84, :y 3.18} {:x 85, :y 3.375} {:x 86, :y 2.475} {:x 87, :y 2.405} {:x 88, :y 2.59} {:x 89, :y 2.625} {:x 90, :y 2.615} {:x 91, :y 2.66} {:x 92, :y 2.605} {:x 93, :y 3.06} {:x 94, :y 2.82} {:x 95, :y 2.845} {:x 96, :y 2.885} {:x 97, :y 2.785} {:x 98, :y 2.865} {:x 99, :y 3.155} {:x 100, :y 2.84} {:x 101, :y 2.985} {:x 102, :y 3.085} {:x 103, :y 2.375} {:x 104, :y 2.775} {:x 105, :y 2.47} {:x 106, :y 2.77} {:x 107, :y 2.76} {:x 108, :y 2.855} {:x 109, :y 3.065} {:x 110, :y 3.1} {:x 111, :y 2.645} {:x 112, :y 2.41} {:x 113, :y 2.565} {:x 114, :y 2.695} {:x 115, :y 2.6} {:x 116, :y 2.32} {:x 117, :y 2.355} {:x 118, :y 2.085} {:x 119, :y 1.985} {:x 120, :y 2.05} {:x 121, :y 1.955} {:x 122, :y 1.91} {:x 123, :y 1.815} {:x 124, :y 2.025} {:x 125, :y 1.97} {:x 126, :y 2.11} {:x 127, :y 2.435} {:x 128, :y 2.58} {:x 129, :y 2.425} {:x 130, :y 2.4} {:x 131, :y 2.63} {:x 132, :y 2.67} {:x 133, :y 2.63} {:x 134, :y 2.685} {:x 135, :y 2.735} {:x 136, :y 2.66} {:x 137, :y 2.785} {:x 138, :y 2.87} {:x 139, :y 2.345} {:x 140, :y 2.435} {:x 141, :y 2.405} {:x 142, :y 2.645} {:x 143, :y 2.855} {:x 144, :y 2.61} {:x 145, :y 2.775} {:x 146, :y 2.81} {:x 147, :y 2.805} {:x 148, :y 3.16} {:x 149, :y 3.18} {:x 150, :y 3.075} {:x 151, :y 2.915} {:x 152, :y 2.25} {:x 153, :y 2.56} {:x 154, :y 2.795} {:x 155, :y 2.835} {:x 156, :y 2.915} {:x 157, :y 2.795} {:x 158, :y 2.785} {:x 159, :y 2.785} {:x 160, :y 2.625} {:x 161, :y 3.095} {:x 162, :y 2.565} {:x 163, :y 2.655} {:x 164, :y 2.875} {:x 165, :y 2.985} {:x 166, :y 2.53} {:x 167, :y 2.52} {:x 168, :y 2.785} {:x 169, :y 3.005} {:x 170, :y 3.06} {:x 171, :y 3.265} {:x 172, :y 3.16} {:x 173, :y 3.135} {:x 174, :y 3.225} {:x 175, :y 3.125} {:x 176, :y 3.005} {:x 177, :y 3.095} {:x 178, :y 3.41} {:x 179, :y 3.12} {:x 180, :y 3.44} {:x 181, :y 3.0} {:x 182, :y 3.17} {:x 183, :y 3.335} {:x 184, :y 3.27} {:x 185, :y 3.48} {:x 186, :y 3.41} {:x 187, :y 3.425} {:x 188, :y 3.65} {:x 189, :y 3.11} {:x 190, :y 3.42} {:x 191, :y 3.465} {:x 192, :y 3.32} {:x 193, :y 3.4} {:x 194, :y 3.38} {:x 195, :y 3.125} {:x 196, :y 3.69} {:x 197, :y 3.455} {:x 198, :y 3.25} {:x 199, :y 3.285})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@

;; @@

;; @@

;; @@
