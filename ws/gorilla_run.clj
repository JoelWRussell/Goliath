;; gorilla-repl.fileformat = 1

;; **
;;; #Golliath 
;; **

;; @@
(ns goliath
  (:require [gorilla-plot.core :as plot]
            [regression-test-score :as ts]
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
            [criterium.core :as criterium]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; **Genetic-algorithm variables**
;; **

;; @@
(def prob_inheritance 0.75)   ;probability that terms are directly inherited from parent to child.
(def mutate_pref 0.5)             ;probability that
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mutate_pref</span>","value":"#'goliath/mutate_pref"}
;; <=

;; **
;;; **Polynomial restrictions**
;; **

;; @@
(def power_sum 10)    ;the maximum sum of powers present in each term within a polynomial.
(def poly_length 10)  ;the maximum number of terms of a polynomial. 
(def power 8) 		  ;the maximum power a variable can be raised to.
(def term_length 2)   ;the number of variables within each term
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/term_length</span>","value":"#'goliath/term_length"}
;; <=

;; **
;;; ##Score Individual
;; **

;; @@
(defn score
  [indv]
  (ts/ind-score (:genotype indv))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/score</span>","value":"#'goliath/score"}
;; <=

;; **
;;; ##Create New Individual 
;;; 
;; **

;; **
;;; Helper function to create new polynomial term
;; **

;; @@
(defn new-poly-term 
  
  "creates a random new poly term satisfying given constraints
   powermax - the maximum power of any given variable, 
   spmax - the max sum of the powers and length is the length of the gene to be generated. 
   spmax - checked and the function re-called if the gene is invalid"
  
  [max_pow spmax length]
  
  (let [new-gene (vec (repeatedly length #(rand-int (+ 1 max_pow)))) sp (apply + new-gene)]
    
    (if (or (< spmax sp) (= new-gene (vec (replicate length 0)) )) 
      (new-poly-term max_pow spmax length)
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
(defn create-individual
  "creates a new individual in the correct format"
  
  [max_length max_power max_sum_powers term_length]
  
  (let [ind {:genotype (create-genotype max_length max_power max_sum_powers term_length)}]
    ind)
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/create-individual</span>","value":"#'goliath/create-individual"}
;; <=

;; @@
(create-individual poly_length power power_sum term_length)

(score (create-individual poly_length power power_sum term_length))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-double'>51.77966959505665</span>","value":"51.77966959505665"}
;; <=

;; **
;;; 
;;; ##Crossover function
;;; 
;; **

;; **
;;; Helper function to distribute genes of a parent according to a probability r ∈ {0,1}
;; **

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
;;; Cross-over function
;; **

;; @@
(defn cross-over
  
  "takes two individuals and crosses over the polynomial from one parent to another under a weighted selection bias" 
  
  [parent-1 parent-2] 
  
  (let [offspring (mapv vec (mapv set (mapv concat (distribute-parent parent-1 prob_inheritance) (reverse (distribute-parent parent-2 prob_inheritance)))))]
    
    [(offspring 0)(offspring 1)]
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/cross-over</span>","value":"#'goliath/cross-over"}
;; <=

;; @@
(def x (:genotype (create-individual poly_length power power_sum term_length)))
x
(def y (:genotype (create-individual poly_length power power_sum term_length)))
y


(cross-over x y)

;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>7</span>","value":"7"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[7 0]"}],"value":"[[7 0]]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"}],"value":"[2 3]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>6</span>","value":"6"}],"value":"[1 6]"}],"value":"[[2 3] [1 6]]"}],"value":"[[[7 0]] [[2 3] [1 6]]]"}
;; <=

;; **
;;; ##Mutate function
;; **

;; **
;;; C
;; **

;; @@
(defn gene-replace
  
  "Replaces a gene from an individual (with a unique gene). In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  
  [indv powermax spmax]
  
(let [length (count (first indv)) new-gene (new-poly-term powermax spmax length) new-geno (assoc indv (rand-int (count indv)) new-gene) sp (apply + new-gene)]
  
  (if (= (count (set new-geno)) (count new-geno))
    new-geno
    (gene-replace indv powermax spmax)
 
 )))    
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-replace</span>","value":"#'goliath/gene-replace"}
;; <=

;; **
;;; 
;; **

;; @@
(defn gene-add
  
  "Adds a new gene to a GA individual. The gene will be unique addition to the individual. In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  
  [indiv powermax spmax]
   
  (let [length (count (first indiv)) new-gene (new-poly-term powermax spmax length) new-geno (vec (conj (set indiv) new-gene)) sp (apply + new-gene)]
   
      (if (= (count new-geno) (count indiv))
      (gene-add indiv powermax spmax)
     new-geno
      )))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-add</span>","value":"#'goliath/gene-add"}
;; <=

;; @@
(defn mutate
  [indv]
  (let [rn (rand)]
    (if (>= rn mutate_pref) 
      (gene-replace indv power power_sum)
      (gene-add indv power power_sum))
    )
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mutate</span>","value":"#'goliath/mutate"}
;; <=

;; @@
(def x' (:genotype (create-individual poly_length power power_sum term_length)))
x'
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>7</span>","value":"7"}],"value":"[2 7]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>6</span>","value":"6"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[6 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>8</span>","value":"8"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[8 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[1 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>8</span>","value":"8"}],"value":"[1 8]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>5</span>","value":"5"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[5 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>6</span>","value":"6"}],"value":"[2 6]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>5</span>","value":"5"}],"value":"[3 5]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[4 4]"}],"value":"[[2 7] [6 4] [8 0] [1 1] [1 8] [5 0] [2 6] [3 5] [4 4]]"}
;; <=

;; @@
(mutate x')
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[1 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[1 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>8</span>","value":"8"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[8 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>8</span>","value":"8"}],"value":"[1 8]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>6</span>","value":"6"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[6 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>7</span>","value":"7"}],"value":"[2 7]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[4 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>6</span>","value":"6"}],"value":"[2 6]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>5</span>","value":"5"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[5 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>5</span>","value":"5"}],"value":"[3 5]"}],"value":"[[1 0] [1 1] [8 0] [1 8] [6 4] [2 7] [4 4] [2 6] [5 0] [3 5]]"}
;; <=

;; **
;;; g##Initial Zeitgiest / Config
;; **

;; @@
(defn random-initial-zeitgeist 
  
   "creates an innitial population of polynimials
  size - number of individuals in the population"
   
  [size max_length max_power max_sum_powers num_terms]

  {:elite [] 
   :rabble (repeatedly size #(create-individual max_length max_power max_sum_powers num_terms)) 
   :age 0}
       
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/random-initial-zeitgeist</span>","value":"#'goliath/random-initial-zeitgeist"}
;; <=

;; @@
(def initial-zeitgeist (random-initial-zeitgeist 10 poly_length power power_sum term_length))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/initial-zeitgeist</span>","value":"#'goliath/initial-zeitgeist"}
;; <=

;; @@
(def generation-config			      
  (let [size-limit 120                   
        min-size 20
        ea-config (spea2/spea2-config
                    {:goals [:score]
                     :archive-size 50
                     :binary-ops [{:op cross-over :repeat 22}]
                     :unary-ops [{:op mutate :repeat 6}]})
        score-functions {:complexity count
                         :error score }]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/generation-config</span>","value":"#'goliath/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (>= (:age zg) 500)))))
;; @@

;; @@

;; @@
