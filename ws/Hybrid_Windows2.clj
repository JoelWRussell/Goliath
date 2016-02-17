;; gorilla-repl.fileformat = 1

;; **
;;; #HYBRID WITH VARIABLE DEGREES OF FREEDOM
;;; 
;;; Join up the GP search with the Mathematica search. You can also enter the number of independent df. Assuming that each column in the data set is a different coordinate.
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
(def experimentalDataSz "mma_coupled_sho.csv")
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
  

(score (create-genotype poly_length power power_sum term_length) )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-double'>0.00986814613767775</span>","value":"0.00986814613767775"}
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
                         :error (fn [e] (score e))}]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/generation-config</span>","value":"#'goliath/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (task zg) (>= (:age zg) 1000000)))))
;; @@
;; ->
;;; 0.1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22.23.24.25.26.27.28.29.30.31.32.33.34.35.36.37.38.39.40.41.42.43.44.45.46.47.48.49.50.51.52.53.54.55.56.57.58.59.60.61.62.63.64.65.66.67.68.
;; <-

;; @@
(goliath.mathlink.LagrangianScore/Shutdown)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(mapv #(println (:genotype %)) (sort-by :error (:elite result)))
;; @@
;; ->
;;; [[1 0 0 0] [0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[1 0 0 0] [0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[1 0 0 0] [0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[1 0 0 0] [0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [0 0 1 0] [1 0 3 3] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[0 1 0 0] [0 3 0 0]]
;;; [[2 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
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
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60","values":[{"x":0,"y":-0.739574232096078},{"x":1,"y":-0.739574232096078},{"x":2,"y":-0.7395742320960781},{"x":3,"y":-0.7395742320960781},{"x":4,"y":-0.7395742320960781},{"x":5,"y":-0.7395742320960781},{"x":6,"y":-0.7395742320960781},{"x":7,"y":-0.7395742320960781},{"x":8,"y":-0.7395742320960781},{"x":9,"y":-0.7395742320960781},{"x":10,"y":-0.7395742320960781},{"x":11,"y":-0.7395742320960781},{"x":12,"y":-0.7395742320960781},{"x":13,"y":-0.7395742320960781},{"x":14,"y":-0.739626463150593},{"x":15,"y":-0.739626463150593},{"x":16,"y":-0.739626463150593},{"x":17,"y":-0.739626463150593},{"x":18,"y":-0.7396264631505937},{"x":19,"y":-0.7396264631505937},{"x":20,"y":-0.7396421673609247},{"x":21,"y":-0.7396421673609248},{"x":22,"y":-0.7396421673609248},{"x":23,"y":-0.7396421673609248},{"x":24,"y":-0.7396421673609248},{"x":25,"y":-0.7396421673609248},{"x":26,"y":-0.7396421673609248},{"x":27,"y":-0.7396421673609248},{"x":28,"y":-0.7396421673609248},{"x":29,"y":-0.7396421673609248},{"x":30,"y":-0.7396421673609248},{"x":31,"y":-0.7396421673609248},{"x":32,"y":-0.7398184705955179},{"x":33,"y":-0.7398184705955181},{"x":34,"y":-0.7398184705955181},{"x":35,"y":-0.7398184705955181},{"x":36,"y":-0.7398366248226959},{"x":37,"y":-0.7398366248226959},{"x":38,"y":-0.7398366248226959},{"x":39,"y":-0.7398366248226959},{"x":40,"y":-0.756375290270588},{"x":41,"y":-0.756375290270588},{"x":42,"y":-0.756375290270588},{"x":43,"y":-0.756375290270588},{"x":44,"y":-0.7563752986087646},{"x":45,"y":-0.7563752986087646},{"x":46,"y":-0.7563752986087646},{"x":47,"y":-0.7563752986087646},{"x":48,"y":-0.7563752986087646},{"x":49,"y":-0.7563752986087646}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"cdfeb39a-363a-4eb2-8b0e-fe14d33d9c60\", :values ({:x 0, :y -0.739574232096078} {:x 1, :y -0.739574232096078} {:x 2, :y -0.7395742320960781} {:x 3, :y -0.7395742320960781} {:x 4, :y -0.7395742320960781} {:x 5, :y -0.7395742320960781} {:x 6, :y -0.7395742320960781} {:x 7, :y -0.7395742320960781} {:x 8, :y -0.7395742320960781} {:x 9, :y -0.7395742320960781} {:x 10, :y -0.7395742320960781} {:x 11, :y -0.7395742320960781} {:x 12, :y -0.7395742320960781} {:x 13, :y -0.7395742320960781} {:x 14, :y -0.739626463150593} {:x 15, :y -0.739626463150593} {:x 16, :y -0.739626463150593} {:x 17, :y -0.739626463150593} {:x 18, :y -0.7396264631505937} {:x 19, :y -0.7396264631505937} {:x 20, :y -0.7396421673609247} {:x 21, :y -0.7396421673609248} {:x 22, :y -0.7396421673609248} {:x 23, :y -0.7396421673609248} {:x 24, :y -0.7396421673609248} {:x 25, :y -0.7396421673609248} {:x 26, :y -0.7396421673609248} {:x 27, :y -0.7396421673609248} {:x 28, :y -0.7396421673609248} {:x 29, :y -0.7396421673609248} {:x 30, :y -0.7396421673609248} {:x 31, :y -0.7396421673609248} {:x 32, :y -0.7398184705955179} {:x 33, :y -0.7398184705955181} {:x 34, :y -0.7398184705955181} {:x 35, :y -0.7398184705955181} {:x 36, :y -0.7398366248226959} {:x 37, :y -0.7398366248226959} {:x 38, :y -0.7398366248226959} {:x 39, :y -0.7398366248226959} {:x 40, :y -0.756375290270588} {:x 41, :y -0.756375290270588} {:x 42, :y -0.756375290270588} {:x 43, :y -0.756375290270588} {:x 44, :y -0.7563752986087646} {:x 45, :y -0.7563752986087646} {:x 46, :y -0.7563752986087646} {:x 47, :y -0.7563752986087646} {:x 48, :y -0.7563752986087646} {:x 49, :y -0.7563752986087646})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
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
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"044cbeb3-b5be-4ed1-a2bd-d20d0579a066","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"044cbeb3-b5be-4ed1-a2bd-d20d0579a066","field":"data.y"}}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"044cbeb3-b5be-4ed1-a2bd-d20d0579a066","values":[{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.5856795786461118,"y":5},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.793474273023734,"y":8},{"x":-1.793474273023734,"y":8},{"x":-1.793474273023734,"y":8},{"x":-1.793474273023734,"y":8},{"x":-1.7125509188174473,"y":6},{"x":-1.7125509188174473,"y":6},{"x":-1.7125509188174473,"y":6},{"x":-1.7125509188174473,"y":6},{"x":0.7897878043966174,"y":1},{"x":0.7897878043966174,"y":1},{"x":0.7897878043966174,"y":1},{"x":0.7897878043966174,"y":1},{"x":-1.793474273023744,"y":9},{"x":-1.793474273023744,"y":9},{"x":-1.793474273023744,"y":9},{"x":-1.7842783704110283,"y":7},{"x":-1.7842783704110283,"y":7}]},{"name":"2afb592d-8da3-4b9b-85a8-6c7b04c3e9c7","values":[{"x":0.0,"y":2},{"x":-0.7395742320960779,"y":3},{"x":-1.4048710778595666,"y":4},{"x":1.5832220019188459,"y":1},{"x":3.462464088822519,"y":3},{"x":-1.585690375890393,"y":6},{"x":1.5832220019188459,"y":1},{"x":-1.7842783704110288,"y":8},{"x":-1.4048710778595666,"y":4},{"x":-0.27868289560261056,"y":8},{"x":-1.404871077859566,"y":5},{"x":0.7897878043966174,"y":1},{"x":4.9268724461346105,"y":1},{"x":3.6038366574673733,"y":3},{"x":-0.9714481405736893,"y":5},{"x":-1.602915639654193,"y":8},{"x":-0.18266406338947105,"y":3},{"x":-0.8298769556455966,"y":4},{"x":-1.4048710778595677,"y":5},{"x":-1.4048717304606915,"y":5},{"x":-1.3159440653984542,"y":7},{"x":-1.3027398682784601,"y":4},{"x":3.2978015814648356,"y":4},{"x":-1.7842783704110283,"y":7},{"x":-1.2940379741452406,"y":3},{"x":23.659505120148808,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-0.1826640633894712,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-0.2239984751410467,"y":4},{"x":3.644567408281659,"y":4},{"x":0.7897878043966174,"y":1},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":1.5832220019188459,"y":1},{"x":-1.4048710778595666,"y":4},{"x":-1.2942684118645715,"y":4},{"x":4.2854616766225435,"y":2},{"x":-1.2942684118645715,"y":4},{"x":3.462464088822519,"y":3},{"x":-1.2942684118645715,"y":4},{"x":-0.8298762548139279,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.793474273023734,"y":8},{"x":-1.4048710778595666,"y":4},{"x":-0.9714376087583153,"y":4},{"x":-0.1826640633894712,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960779,"y":3},{"x":-0.27868289560261084,"y":7},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-1.793474273023734,"y":8},{"x":4.9268724461346105,"y":1},{"x":-1.793474273023734,"y":8},{"x":0.7897878043966174,"y":1},{"x":-0.24961490963558822,"y":6},{"x":-1.3168230245618102,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-1.6451587311069127,"y":7},{"x":-1.11510463526607,"y":5},{"x":-1.5533118697448842,"y":6},{"x":3.207546698175711,"y":6},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.7125509188174473,"y":6},{"x":4.571826979507611,"y":5},{"x":18.69294551456987,"y":1},{"x":-0.8299303674012493,"y":4},{"x":-0.7395742320960781,"y":2},{"x":4.9268724461346105,"y":1},{"x":4.9268724461346105,"y":1},{"x":-0.203127580816704,"y":3},{"x":4.831646713702918,"y":2},{"x":-0.24297739738815866,"y":7},{"x":3.443914207095526,"y":6},{"x":-1.5856795786461162,"y":6},{"x":-1.793474273023734,"y":8},{"x":-0.8298762548139279,"y":3},{"x":-1.5856795786461118,"y":5},{"x":4.9268724461346105,"y":1},{"x":-0.014805577146587347,"y":4},{"x":3.4480588132702614,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.5856795786461118,"y":5},{"x":4.9268724461346105,"y":1},{"x":-0.7395742320960781,"y":2},{"x":0.7897878043966174,"y":1},{"x":-0.7395742320960781,"y":2},{"x":-1.6029156396541955,"y":7},{"x":0.7890696318008138,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2}]},{"name":"c796ea87-66ea-4f1c-a3c9-d973db6c0475","values":[{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.5856795786461118,"y":5},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.4048710778595666,"y":4},{"x":-1.5856795786461118,"y":5},{"x":-1.2940379741452406,"y":3},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.2940379741452406,"y":3},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-0.7395742320960781,"y":2},{"x":-1.4048710778595666,"y":4},{"x":-1.4048710778595666,"y":4},{"x":-0.7395742320960781,"y":2},{"x":-1.793474273023734,"y":8},{"x":-1.793474273023734,"y":8},{"x":-1.793474273023734,"y":8},{"x":-1.793474273023734,"y":8},{"x":-1.7125509188174473,"y":6},{"x":-1.7125509188174473,"y":6},{"x":-1.7125509188174473,"y":6},{"x":-1.7125509188174473,"y":6},{"x":0.7897878043966174,"y":1},{"x":0.7897878043966174,"y":1},{"x":0.7897878043966174,"y":1},{"x":0.7897878043966174,"y":1},{"x":-1.793474273023744,"y":9},{"x":-1.793474273023744,"y":9},{"x":-1.793474273023744,"y":9},{"x":-1.7842783704110283,"y":7},{"x":-1.7842783704110283,"y":7}]}],"marks":[{"type":"symbol","from":{"data":"044cbeb3-b5be-4ed1-a2bd-d20d0579a066"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"red"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"2afb592d-8da3-4b9b-85a8-6c7b04c3e9c7"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"blue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"c796ea87-66ea-4f1c-a3c9-d973db6c0475"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"#ff29d2"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"044cbeb3-b5be-4ed1-a2bd-d20d0579a066\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"044cbeb3-b5be-4ed1-a2bd-d20d0579a066\", :field \"data.y\"}}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"044cbeb3-b5be-4ed1-a2bd-d20d0579a066\", :values ({:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.5856795786461118, :y 5} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.793474273023734, :y 8} {:x -1.793474273023734, :y 8} {:x -1.793474273023734, :y 8} {:x -1.793474273023734, :y 8} {:x -1.7125509188174473, :y 6} {:x -1.7125509188174473, :y 6} {:x -1.7125509188174473, :y 6} {:x -1.7125509188174473, :y 6} {:x 0.7897878043966174, :y 1} {:x 0.7897878043966174, :y 1} {:x 0.7897878043966174, :y 1} {:x 0.7897878043966174, :y 1} {:x -1.793474273023744, :y 9} {:x -1.793474273023744, :y 9} {:x -1.793474273023744, :y 9} {:x -1.7842783704110283, :y 7} {:x -1.7842783704110283, :y 7})} {:name \"2afb592d-8da3-4b9b-85a8-6c7b04c3e9c7\", :values ({:x 0.0, :y 2} {:x -0.7395742320960779, :y 3} {:x -1.4048710778595666, :y 4} {:x 1.5832220019188459, :y 1} {:x 3.462464088822519, :y 3} {:x -1.585690375890393, :y 6} {:x 1.5832220019188459, :y 1} {:x -1.7842783704110288, :y 8} {:x -1.4048710778595666, :y 4} {:x -0.27868289560261056, :y 8} {:x -1.404871077859566, :y 5} {:x 0.7897878043966174, :y 1} {:x 4.9268724461346105, :y 1} {:x 3.6038366574673733, :y 3} {:x -0.9714481405736893, :y 5} {:x -1.602915639654193, :y 8} {:x -0.18266406338947105, :y 3} {:x -0.8298769556455966, :y 4} {:x -1.4048710778595677, :y 5} {:x -1.4048717304606915, :y 5} {:x -1.3159440653984542, :y 7} {:x -1.3027398682784601, :y 4} {:x 3.2978015814648356, :y 4} {:x -1.7842783704110283, :y 7} {:x -1.2940379741452406, :y 3} {:x 23.659505120148808, :y 2} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -0.1826640633894712, :y 2} {:x -1.2940379741452406, :y 3} {:x -0.2239984751410467, :y 4} {:x 3.644567408281659, :y 4} {:x 0.7897878043966174, :y 1} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x 1.5832220019188459, :y 1} {:x -1.4048710778595666, :y 4} {:x -1.2942684118645715, :y 4} {:x 4.2854616766225435, :y 2} {:x -1.2942684118645715, :y 4} {:x 3.462464088822519, :y 3} {:x -1.2942684118645715, :y 4} {:x -0.8298762548139279, :y 3} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.793474273023734, :y 8} {:x -1.4048710778595666, :y 4} {:x -0.9714376087583153, :y 4} {:x -0.1826640633894712, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960779, :y 3} {:x -0.27868289560261084, :y 7} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -1.793474273023734, :y 8} {:x 4.9268724461346105, :y 1} {:x -1.793474273023734, :y 8} {:x 0.7897878043966174, :y 1} {:x -0.24961490963558822, :y 6} {:x -1.3168230245618102, :y 4} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -1.6451587311069127, :y 7} {:x -1.11510463526607, :y 5} {:x -1.5533118697448842, :y 6} {:x 3.207546698175711, :y 6} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.7125509188174473, :y 6} {:x 4.571826979507611, :y 5} {:x 18.69294551456987, :y 1} {:x -0.8299303674012493, :y 4} {:x -0.7395742320960781, :y 2} {:x 4.9268724461346105, :y 1} {:x 4.9268724461346105, :y 1} {:x -0.203127580816704, :y 3} {:x 4.831646713702918, :y 2} {:x -0.24297739738815866, :y 7} {:x 3.443914207095526, :y 6} {:x -1.5856795786461162, :y 6} {:x -1.793474273023734, :y 8} {:x -0.8298762548139279, :y 3} {:x -1.5856795786461118, :y 5} {:x 4.9268724461346105, :y 1} {:x -0.014805577146587347, :y 4} {:x 3.4480588132702614, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.5856795786461118, :y 5} {:x 4.9268724461346105, :y 1} {:x -0.7395742320960781, :y 2} {:x 0.7897878043966174, :y 1} {:x -0.7395742320960781, :y 2} {:x -1.6029156396541955, :y 7} {:x 0.7890696318008138, :y 3} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2})} {:name \"c796ea87-66ea-4f1c-a3c9-d973db6c0475\", :values ({:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.5856795786461118, :y 5} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.4048710778595666, :y 4} {:x -1.5856795786461118, :y 5} {:x -1.2940379741452406, :y 3} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.2940379741452406, :y 3} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -0.7395742320960781, :y 2} {:x -1.4048710778595666, :y 4} {:x -1.4048710778595666, :y 4} {:x -0.7395742320960781, :y 2} {:x -1.793474273023734, :y 8} {:x -1.793474273023734, :y 8} {:x -1.793474273023734, :y 8} {:x -1.793474273023734, :y 8} {:x -1.7125509188174473, :y 6} {:x -1.7125509188174473, :y 6} {:x -1.7125509188174473, :y 6} {:x -1.7125509188174473, :y 6} {:x 0.7897878043966174, :y 1} {:x 0.7897878043966174, :y 1} {:x 0.7897878043966174, :y 1} {:x 0.7897878043966174, :y 1} {:x -1.793474273023744, :y 9} {:x -1.793474273023744, :y 9} {:x -1.793474273023744, :y 9} {:x -1.7842783704110283, :y 7} {:x -1.7842783704110283, :y 7})}), :marks ({:type \"symbol\", :from {:data \"044cbeb3-b5be-4ed1-a2bd-d20d0579a066\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"red\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"2afb592d-8da3-4b9b-85a8-6c7b04c3e9c7\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"blue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"c796ea87-66ea-4f1c-a3c9-d973db6c0475\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"#ff29d2\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}})}}"}
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
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"3305d4ee-546e-47fd-8445-9129e16ed5e1","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":[0,5308]}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"3305d4ee-546e-47fd-8445-9129e16ed5e1","values":[{"x":0,"y":2990},{"x":1,"y":2899},{"x":2,"y":2910},{"x":3,"y":2432},{"x":4,"y":2790},{"x":5,"y":2556},{"x":6,"y":3096},{"x":7,"y":2664},{"x":8,"y":2812},{"x":9,"y":3167},{"x":10,"y":3029},{"x":11,"y":2729},{"x":12,"y":2586},{"x":13,"y":3037},{"x":14,"y":2698},{"x":15,"y":2257},{"x":16,"y":2423},{"x":17,"y":2561},{"x":18,"y":2167},{"x":19,"y":2441},{"x":20,"y":1401},{"x":21,"y":1995},{"x":22,"y":1951},{"x":23,"y":2875},{"x":24,"y":2629},{"x":25,"y":1897},{"x":26,"y":2143},{"x":27,"y":2578},{"x":28,"y":3530},{"x":29,"y":3568},{"x":30,"y":2073},{"x":31,"y":1734},{"x":32,"y":2649},{"x":33,"y":2047},{"x":34,"y":3157},{"x":35,"y":1751},{"x":36,"y":2421},{"x":37,"y":2661},{"x":38,"y":3834},{"x":39,"y":3471},{"x":40,"y":3840},{"x":41,"y":5308},{"x":42,"y":4727},{"x":43,"y":4178},{"x":44,"y":3070},{"x":45,"y":2290},{"x":46,"y":1554},{"x":47,"y":1861},{"x":48,"y":1991},{"x":49,"y":2273}]},{"name":"ce3f07cd-fae4-4ef7-8b80-909be7b6388a","values":[{"x":0,"y":72},{"x":1,"y":43},{"x":2,"y":31},{"x":3,"y":26},{"x":4,"y":30},{"x":5,"y":24},{"x":6,"y":22},{"x":7,"y":19},{"x":8,"y":20},{"x":9,"y":22},{"x":10,"y":24},{"x":11,"y":19},{"x":12,"y":46},{"x":13,"y":29},{"x":14,"y":44},{"x":15,"y":49},{"x":16,"y":58},{"x":17,"y":70},{"x":18,"y":58},{"x":19,"y":52},{"x":20,"y":30},{"x":21,"y":55},{"x":22,"y":55},{"x":23,"y":54},{"x":24,"y":57},{"x":25,"y":26},{"x":26,"y":50},{"x":27,"y":59},{"x":28,"y":66},{"x":29,"y":44},{"x":30,"y":35},{"x":31,"y":49},{"x":32,"y":45},{"x":33,"y":31},{"x":34,"y":40},{"x":35,"y":36},{"x":36,"y":46},{"x":37,"y":48},{"x":38,"y":42},{"x":39,"y":36},{"x":40,"y":36},{"x":41,"y":35},{"x":42,"y":29},{"x":43,"y":25},{"x":44,"y":33},{"x":45,"y":19},{"x":46,"y":43},{"x":47,"y":50},{"x":48,"y":56},{"x":49,"y":57}]},{"name":"5a5e8321-4927-48fa-88e9-a80a959281cf","values":[{"x":0,"y":123},{"x":1,"y":3},{"x":2,"y":4},{"x":3,"y":3},{"x":4,"y":2},{"x":5,"y":4},{"x":6,"y":2},{"x":7,"y":1},{"x":8,"y":1},{"x":9,"y":2},{"x":10,"y":1},{"x":11,"y":1},{"x":12,"y":1},{"x":13,"y":1},{"x":14,"y":2},{"x":15,"y":1},{"x":16,"y":1},{"x":17,"y":1},{"x":18,"y":1},{"x":19,"y":1},{"x":20,"y":2},{"x":21,"y":1},{"x":22,"y":1},{"x":23,"y":0},{"x":24,"y":1},{"x":25,"y":1},{"x":26,"y":1},{"x":27,"y":0},{"x":28,"y":1},{"x":29,"y":1},{"x":30,"y":1},{"x":31,"y":1},{"x":32,"y":0},{"x":33,"y":0},{"x":34,"y":1},{"x":35,"y":1},{"x":36,"y":1},{"x":37,"y":1},{"x":38,"y":1},{"x":39,"y":0},{"x":40,"y":1},{"x":41,"y":1},{"x":42,"y":1},{"x":43,"y":1},{"x":44,"y":0},{"x":45,"y":1},{"x":46,"y":0},{"x":47,"y":1},{"x":48,"y":1},{"x":49,"y":1}]},{"name":"a44e7422-c37c-4fcd-9162-f5d798be50d6","values":[{"x":0,"y":2795},{"x":1,"y":2853},{"x":2,"y":2875},{"x":3,"y":2403},{"x":4,"y":2758},{"x":5,"y":2528},{"x":6,"y":3072},{"x":7,"y":2644},{"x":8,"y":2791},{"x":9,"y":3143},{"x":10,"y":3004},{"x":11,"y":2709},{"x":12,"y":2539},{"x":13,"y":3007},{"x":14,"y":2652},{"x":15,"y":2207},{"x":16,"y":2364},{"x":17,"y":2490},{"x":18,"y":2108},{"x":19,"y":2388},{"x":20,"y":1369},{"x":21,"y":1939},{"x":22,"y":1895},{"x":23,"y":2821},{"x":24,"y":2571},{"x":25,"y":1870},{"x":26,"y":2092},{"x":27,"y":2519},{"x":28,"y":3463},{"x":29,"y":3523},{"x":30,"y":2037},{"x":31,"y":1684},{"x":32,"y":2604},{"x":33,"y":2016},{"x":34,"y":3116},{"x":35,"y":1714},{"x":36,"y":2374},{"x":37,"y":2612},{"x":38,"y":3791},{"x":39,"y":3435},{"x":40,"y":3803},{"x":41,"y":5272},{"x":42,"y":4697},{"x":43,"y":4152},{"x":44,"y":3037},{"x":45,"y":2270},{"x":46,"y":1511},{"x":47,"y":1810},{"x":48,"y":1934},{"x":49,"y":2215}]}],"marks":[{"type":"line","from":{"data":"3305d4ee-546e-47fd-8445-9129e16ed5e1"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"yellow"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"ce3f07cd-fae4-4ef7-8b80-909be7b6388a"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"red"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"5a5e8321-4927-48fa-88e9-a80a959281cf"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"green"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"a44e7422-c37c-4fcd-9162-f5d798be50d6"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"blue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"3305d4ee-546e-47fd-8445-9129e16ed5e1\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain [0 5308]}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"3305d4ee-546e-47fd-8445-9129e16ed5e1\", :values ({:x 0, :y 2990} {:x 1, :y 2899} {:x 2, :y 2910} {:x 3, :y 2432} {:x 4, :y 2790} {:x 5, :y 2556} {:x 6, :y 3096} {:x 7, :y 2664} {:x 8, :y 2812} {:x 9, :y 3167} {:x 10, :y 3029} {:x 11, :y 2729} {:x 12, :y 2586} {:x 13, :y 3037} {:x 14, :y 2698} {:x 15, :y 2257} {:x 16, :y 2423} {:x 17, :y 2561} {:x 18, :y 2167} {:x 19, :y 2441} {:x 20, :y 1401} {:x 21, :y 1995} {:x 22, :y 1951} {:x 23, :y 2875} {:x 24, :y 2629} {:x 25, :y 1897} {:x 26, :y 2143} {:x 27, :y 2578} {:x 28, :y 3530} {:x 29, :y 3568} {:x 30, :y 2073} {:x 31, :y 1734} {:x 32, :y 2649} {:x 33, :y 2047} {:x 34, :y 3157} {:x 35, :y 1751} {:x 36, :y 2421} {:x 37, :y 2661} {:x 38, :y 3834} {:x 39, :y 3471} {:x 40, :y 3840} {:x 41, :y 5308} {:x 42, :y 4727} {:x 43, :y 4178} {:x 44, :y 3070} {:x 45, :y 2290} {:x 46, :y 1554} {:x 47, :y 1861} {:x 48, :y 1991} {:x 49, :y 2273})} {:name \"ce3f07cd-fae4-4ef7-8b80-909be7b6388a\", :values ({:x 0, :y 72} {:x 1, :y 43} {:x 2, :y 31} {:x 3, :y 26} {:x 4, :y 30} {:x 5, :y 24} {:x 6, :y 22} {:x 7, :y 19} {:x 8, :y 20} {:x 9, :y 22} {:x 10, :y 24} {:x 11, :y 19} {:x 12, :y 46} {:x 13, :y 29} {:x 14, :y 44} {:x 15, :y 49} {:x 16, :y 58} {:x 17, :y 70} {:x 18, :y 58} {:x 19, :y 52} {:x 20, :y 30} {:x 21, :y 55} {:x 22, :y 55} {:x 23, :y 54} {:x 24, :y 57} {:x 25, :y 26} {:x 26, :y 50} {:x 27, :y 59} {:x 28, :y 66} {:x 29, :y 44} {:x 30, :y 35} {:x 31, :y 49} {:x 32, :y 45} {:x 33, :y 31} {:x 34, :y 40} {:x 35, :y 36} {:x 36, :y 46} {:x 37, :y 48} {:x 38, :y 42} {:x 39, :y 36} {:x 40, :y 36} {:x 41, :y 35} {:x 42, :y 29} {:x 43, :y 25} {:x 44, :y 33} {:x 45, :y 19} {:x 46, :y 43} {:x 47, :y 50} {:x 48, :y 56} {:x 49, :y 57})} {:name \"5a5e8321-4927-48fa-88e9-a80a959281cf\", :values ({:x 0, :y 123} {:x 1, :y 3} {:x 2, :y 4} {:x 3, :y 3} {:x 4, :y 2} {:x 5, :y 4} {:x 6, :y 2} {:x 7, :y 1} {:x 8, :y 1} {:x 9, :y 2} {:x 10, :y 1} {:x 11, :y 1} {:x 12, :y 1} {:x 13, :y 1} {:x 14, :y 2} {:x 15, :y 1} {:x 16, :y 1} {:x 17, :y 1} {:x 18, :y 1} {:x 19, :y 1} {:x 20, :y 2} {:x 21, :y 1} {:x 22, :y 1} {:x 23, :y 0} {:x 24, :y 1} {:x 25, :y 1} {:x 26, :y 1} {:x 27, :y 0} {:x 28, :y 1} {:x 29, :y 1} {:x 30, :y 1} {:x 31, :y 1} {:x 32, :y 0} {:x 33, :y 0} {:x 34, :y 1} {:x 35, :y 1} {:x 36, :y 1} {:x 37, :y 1} {:x 38, :y 1} {:x 39, :y 0} {:x 40, :y 1} {:x 41, :y 1} {:x 42, :y 1} {:x 43, :y 1} {:x 44, :y 0} {:x 45, :y 1} {:x 46, :y 0} {:x 47, :y 1} {:x 48, :y 1} {:x 49, :y 1})} {:name \"a44e7422-c37c-4fcd-9162-f5d798be50d6\", :values ({:x 0, :y 2795} {:x 1, :y 2853} {:x 2, :y 2875} {:x 3, :y 2403} {:x 4, :y 2758} {:x 5, :y 2528} {:x 6, :y 3072} {:x 7, :y 2644} {:x 8, :y 2791} {:x 9, :y 3143} {:x 10, :y 3004} {:x 11, :y 2709} {:x 12, :y 2539} {:x 13, :y 3007} {:x 14, :y 2652} {:x 15, :y 2207} {:x 16, :y 2364} {:x 17, :y 2490} {:x 18, :y 2108} {:x 19, :y 2388} {:x 20, :y 1369} {:x 21, :y 1939} {:x 22, :y 1895} {:x 23, :y 2821} {:x 24, :y 2571} {:x 25, :y 1870} {:x 26, :y 2092} {:x 27, :y 2519} {:x 28, :y 3463} {:x 29, :y 3523} {:x 30, :y 2037} {:x 31, :y 1684} {:x 32, :y 2604} {:x 33, :y 2016} {:x 34, :y 3116} {:x 35, :y 1714} {:x 36, :y 2374} {:x 37, :y 2612} {:x 38, :y 3791} {:x 39, :y 3435} {:x 40, :y 3803} {:x 41, :y 5272} {:x 42, :y 4697} {:x 43, :y 4152} {:x 44, :y 3037} {:x 45, :y 2270} {:x 46, :y 1511} {:x 47, :y 1810} {:x 48, :y 1934} {:x 49, :y 2215})}), :marks ({:type \"line\", :from {:data \"3305d4ee-546e-47fd-8445-9129e16ed5e1\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"yellow\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"ce3f07cd-fae4-4ef7-8b80-909be7b6388a\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"red\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"5a5e8321-4927-48fa-88e9-a80a959281cf\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"green\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"a44e7422-c37c-4fcd-9162-f5d798be50d6\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"blue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}})}}"}
;; <=

;; **
;;; Spread of complexity in population
;; **

;; @@
(plot/histogram (map :complexity (:rabble result)))
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"interpolate":{"value":"step-before"},"fill":{"value":"steelblue"},"fillOpacity":{"value":0.4},"stroke":{"value":"steelblue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4","values":[{"x":1.0,"y":0},{"x":2.1428571428571432,"y":32.0},{"x":3.2857142857142865,"y":3.0},{"x":4.42857142857143,"y":6.0},{"x":5.571428571428573,"y":4.0},{"x":6.714285714285716,"y":0.0},{"x":7.857142857142859,"y":1.0},{"x":9.000000000000002,"y":4.0},{"x":10.142857142857144,"y":0}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :interpolate {:value \"step-before\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 0.4}, :stroke {:value \"steelblue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"a51a76cd-7b5e-4b29-a6d7-55795ceef0d4\", :values ({:x 1.0, :y 0} {:x 2.1428571428571432, :y 32.0} {:x 3.2857142857142865, :y 3.0} {:x 4.42857142857143, :y 6.0} {:x 5.571428571428573, :y 4.0} {:x 6.714285714285716, :y 0.0} {:x 7.857142857142859, :y 1.0} {:x 9.000000000000002, :y 4.0} {:x 10.142857142857144, :y 0})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(plot/list-plot (:mean (:complexity @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"79879f9c-758a-49f4-884e-b18b78961009","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"79879f9c-758a-49f4-884e-b18b78961009","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"79879f9c-758a-49f4-884e-b18b78961009"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"79879f9c-758a-49f4-884e-b18b78961009","values":[{"x":0,"y":4.29},{"x":1,"y":4.86},{"x":2,"y":2.71},{"x":3,"y":3.53},{"x":4,"y":4.38},{"x":5,"y":5.71},{"x":6,"y":4.69},{"x":7,"y":5.76},{"x":8,"y":5.54},{"x":9,"y":5.15},{"x":10,"y":5.89},{"x":11,"y":5.53},{"x":12,"y":5.68},{"x":13,"y":8.67},{"x":14,"y":7.66},{"x":15,"y":6.36},{"x":16,"y":7.66},{"x":17,"y":7.37},{"x":18,"y":6.11},{"x":19,"y":7.86},{"x":20,"y":6.68},{"x":21,"y":7.67},{"x":22,"y":8.01},{"x":23,"y":7.75},{"x":24,"y":6.49},{"x":25,"y":6.0},{"x":26,"y":5.67},{"x":27,"y":6.9},{"x":28,"y":5.7},{"x":29,"y":6.6},{"x":30,"y":3.93},{"x":31,"y":3.31},{"x":32,"y":3.81},{"x":33,"y":3.72},{"x":34,"y":5.6},{"x":35,"y":5.59},{"x":36,"y":5.12},{"x":37,"y":6.32},{"x":38,"y":5.58},{"x":39,"y":5.05},{"x":40,"y":5.61},{"x":41,"y":6.08},{"x":42,"y":3.88},{"x":43,"y":3.47},{"x":44,"y":3.78},{"x":45,"y":4.48},{"x":46,"y":6.35},{"x":47,"y":7.45},{"x":48,"y":6.45},{"x":49,"y":4.88},{"x":50,"y":5.14},{"x":51,"y":4.5},{"x":52,"y":3.81},{"x":53,"y":3.44},{"x":54,"y":3.6},{"x":55,"y":3.76},{"x":56,"y":3.67},{"x":57,"y":3.61},{"x":58,"y":3.75},{"x":59,"y":4.38},{"x":60,"y":4.19},{"x":61,"y":4.79},{"x":62,"y":4.67},{"x":63,"y":4.85},{"x":64,"y":5.38},{"x":65,"y":5.04},{"x":66,"y":5.27},{"x":67,"y":5.85},{"x":68,"y":4.91},{"x":69,"y":6.57},{"x":70,"y":5.46},{"x":71,"y":5.19},{"x":72,"y":4.77},{"x":73,"y":6.28},{"x":74,"y":5.38},{"x":75,"y":6.44},{"x":76,"y":6.97},{"x":77,"y":6.13},{"x":78,"y":6.55},{"x":79,"y":5.0},{"x":80,"y":4.91},{"x":81,"y":5.85},{"x":82,"y":5.99},{"x":83,"y":5.45},{"x":84,"y":5.86},{"x":85,"y":5.47},{"x":86,"y":6.35},{"x":87,"y":6.25},{"x":88,"y":5.47},{"x":89,"y":5.72},{"x":90,"y":6.24},{"x":91,"y":5.78},{"x":92,"y":6.5},{"x":93,"y":6.06},{"x":94,"y":6.79},{"x":95,"y":6.32},{"x":96,"y":6.35},{"x":97,"y":5.78},{"x":98,"y":6.27},{"x":99,"y":6.55}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"79879f9c-758a-49f4-884e-b18b78961009\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"79879f9c-758a-49f4-884e-b18b78961009\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"79879f9c-758a-49f4-884e-b18b78961009\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"79879f9c-758a-49f4-884e-b18b78961009\", :values ({:x 0, :y 4.29} {:x 1, :y 4.86} {:x 2, :y 2.71} {:x 3, :y 3.53} {:x 4, :y 4.38} {:x 5, :y 5.71} {:x 6, :y 4.69} {:x 7, :y 5.76} {:x 8, :y 5.54} {:x 9, :y 5.15} {:x 10, :y 5.89} {:x 11, :y 5.53} {:x 12, :y 5.68} {:x 13, :y 8.67} {:x 14, :y 7.66} {:x 15, :y 6.36} {:x 16, :y 7.66} {:x 17, :y 7.37} {:x 18, :y 6.11} {:x 19, :y 7.86} {:x 20, :y 6.68} {:x 21, :y 7.67} {:x 22, :y 8.01} {:x 23, :y 7.75} {:x 24, :y 6.49} {:x 25, :y 6.0} {:x 26, :y 5.67} {:x 27, :y 6.9} {:x 28, :y 5.7} {:x 29, :y 6.6} {:x 30, :y 3.93} {:x 31, :y 3.31} {:x 32, :y 3.81} {:x 33, :y 3.72} {:x 34, :y 5.6} {:x 35, :y 5.59} {:x 36, :y 5.12} {:x 37, :y 6.32} {:x 38, :y 5.58} {:x 39, :y 5.05} {:x 40, :y 5.61} {:x 41, :y 6.08} {:x 42, :y 3.88} {:x 43, :y 3.47} {:x 44, :y 3.78} {:x 45, :y 4.48} {:x 46, :y 6.35} {:x 47, :y 7.45} {:x 48, :y 6.45} {:x 49, :y 4.88} {:x 50, :y 5.14} {:x 51, :y 4.5} {:x 52, :y 3.81} {:x 53, :y 3.44} {:x 54, :y 3.6} {:x 55, :y 3.76} {:x 56, :y 3.67} {:x 57, :y 3.61} {:x 58, :y 3.75} {:x 59, :y 4.38} {:x 60, :y 4.19} {:x 61, :y 4.79} {:x 62, :y 4.67} {:x 63, :y 4.85} {:x 64, :y 5.38} {:x 65, :y 5.04} {:x 66, :y 5.27} {:x 67, :y 5.85} {:x 68, :y 4.91} {:x 69, :y 6.57} {:x 70, :y 5.46} {:x 71, :y 5.19} {:x 72, :y 4.77} {:x 73, :y 6.28} {:x 74, :y 5.38} {:x 75, :y 6.44} {:x 76, :y 6.97} {:x 77, :y 6.13} {:x 78, :y 6.55} {:x 79, :y 5.0} {:x 80, :y 4.91} {:x 81, :y 5.85} {:x 82, :y 5.99} {:x 83, :y 5.45} {:x 84, :y 5.86} {:x 85, :y 5.47} {:x 86, :y 6.35} {:x 87, :y 6.25} {:x 88, :y 5.47} {:x 89, :y 5.72} {:x 90, :y 6.24} {:x 91, :y 5.78} {:x 92, :y 6.5} {:x 93, :y 6.06} {:x 94, :y 6.79} {:x 95, :y 6.32} {:x 96, :y 6.35} {:x 97, :y 5.78} {:x 98, :y 6.27} {:x 99, :y 6.55})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@

;; @@
