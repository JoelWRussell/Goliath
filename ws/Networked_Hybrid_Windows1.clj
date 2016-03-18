;; gorilla-repl.fileformat = 1

;; **
;;; #Networked Hybrid Windows 1
;;; This worksheet uses a class GorillaClient which talks to LagrangeServer (which can be on another computer). LagrangeServer sends the zeitgeist out to WorkerClient computers which have Mathematica installed. So LagrangeServer organises and distributes the work and collects and returns the results back to GorillaClient which returns the results back to this worksheet. 
;;; 
;;; 
;;; 
;; **

;; @@
(def experimentalDataSz "resources/mma_double.csv")
(def dt 0.1)
(def df 2)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/df</span>","value":"#'user/df"}
;; <=

;; @@
(import '[com.lagrangianmining.GorillaClient])
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; Create an instance of a GorillaClient. This knows how to find and connect with LagrangeServer. In fact it contacts the web to find the ip:port of LagrangeServer. Obviously the LagrangeServer must be running and must have at least one WorkerClient registered with it.
;; **

;; @@
(def client (com.lagrangianmining.GorillaClient.))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/client</span>","value":"#'user/client"}
;; <=

;; **
;;; A population is scored at a time via the network of computers. This function is called by the scoring functions.
;; **

;; @@
(defn networkScore [gen]
(let[
sizePoly (map #(count (flatten %)) gen)]
(.NewZeitgeist user/client (count gen) (into-array Integer/TYPE sizePoly) (into-array Integer/TYPE (flatten gen)) user/df, user/dt)
(do (into [] (.GetScores user/client)))
)
)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/networkScore</span>","value":"#'user/networkScore"}
;; <=

;; **
;;; InitMathKernel will tell LagrangeServer to launch the MathKernels of all of the WorkerClient machines.
;; **

;; @@
(.InitMathKernel client)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; **
;;; NewData sends a copy of the experimental data to all of the WorkerClients for them to use in their scoring.
;; **

;; @@
(.NewData client "resources/mma_double.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; **
;;; PrepareData tells each WorkerClient to prepare itself by loading the experimental data and the mma script
;; **

;; @@
(.PrepareData client false df dt)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>false</span>","value":"false"}
;; <=

;; @@
(require '[darwin.evolution.core :as evolution]
         '[darwin.evolution.metrics :as metrics]
         '[darwin.evolution.reproduction :as reproduction]
         '[darwin.evolution.scoring :as scoring]
         '[darwin.evolution.selection :as selection]
         '[darwin.evolution.transform :as transform]
         '[darwin.evolution.pareto :as pareto]
         '[darwin.algorithms.spea2 :as spea2]
    	 '[darwin.core :as darwin]
         '[darwin.evolution.metrics :as metrics]
         '[algebolic.expression.core :as expression]
         '[algebolic.expression.tree :as tree]
         '[algebolic.expression.genetics :as genetics]
         '[algebolic.expression.score :as score]
         '[algebolic.expression.render :as spea-render]
         '[algebolic.expression.interpreter :as interpreter]
         '[criterium.core :as criterium]
         )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; Now for the previous GA code
;; **

;; @@
(def prob_inheritance 0.75) ;;to do with crossover  ;probability that terms are directly inherited from parent to child.
(def mutate_pref1 0.95)
(def mutate_pref2 0.8);;change 1 df within a gene  ;prob. that "gene-replace"/"gene-add" is called when mutate is called.
(def mutate_pref3 0.2) ;; add a whole new gene/replace a term (big change)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/mutate_pref3</span>","value":"#'user/mutate_pref3"}
;; <=

;; @@
(def power_sum 10)    ;the maximum sum of powers present in each term within a polynomial.
(def poly_length 10)  ;the maximum number of terms of a polynomial. 
(def power 8) 		  ;the maximum power a variable can be raised to.
(def term_length (* 2 user/df))   ;the number of variables within each term - 2*degree_freedom
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/term_length</span>","value":"#'user/term_length"}
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
      
      ;(assoc (vec (repeat length 0)) (rand-int length) 1);;;;;;change here
    
      (new-poly-term max_pow spmax length)
      ;[0 1];;cant have this or else it will crash the mathematica
      new-gene
      )
  ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/new-poly-term</span>","value":"#'user/new-poly-term"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/create-genotype</span>","value":"#'user/create-genotype"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/distribute-parent</span>","value":"#'user/distribute-parent"}
;; <=

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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/cross-over</span>","value":"#'user/cross-over"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-replace</span>","value":"#'user/gene-replace"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-add</span>","value":"#'user/gene-add"}
;; <=

;; @@
(defn gene-tweak
  [indv]
  (let [n (rand-int (count indv)) gene (nth indv n) new-gene (assoc gene (rand-int (count gene)) (rand-int power))]
    (assoc indv n new-gene)
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-tweak</span>","value":"#'user/gene-tweak"}
;; <=

;; @@
(defn gene-delete [indv] 
  (let [new-indv (subvec indv 1 (count indv))]
    (if (= new-indv []) indv
      new-indv))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-delete</span>","value":"#'user/gene-delete"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/mutate</span>","value":"#'user/mutate"}
;; <=

;; @@
(defn random-initial-population 
  
   "creates an innitial population of polynimials
  size - number of individuals in the population"
   
  [size max_length max_power max_sum_powers num_terms]
 
  (repeatedly size #(create-genotype max_length max_power max_sum_powers num_terms)) 
       
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/random-initial-population</span>","value":"#'user/random-initial-population"}
;; <=

;; **
;;; Create an initial population.
;; **

;; @@
(def initial-zeitgeist (evolution/make-zeitgeist (random-initial-population 100 poly_length power power_sum term_length)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/initial-zeitgeist</span>","value":"#'user/initial-zeitgeist"}
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
        score-functions {:complexity (fn [x] (x))
                         :error (fn [e] e}]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/generation-config</span>","value":"#'user/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (>= (:age zg) 200)))))
;; @@
;; ->
;;; ........................................................................................................................................................................................................&quot;Elapsed time: 299093.669443 msecs&quot;
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/result</span>","value":"#'user/result"}
;; <=

;; **
;;; Destroy all of the MathKernels on the WorkerClients.
;; **

;; @@
(.DestroyMathKernel client)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; **
;;; Break link with LagrangeServer.
;; **

;; @@
(.DisconnectServer user/client)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; Get Results
;; **

;; @@
(mapv #(println (:genotype %)) (sort-by :error (:rabble result)))
;; @@
;; ->
;;; [[1 7 1 3]]
;;; [[2 1 6 1]]
;;; [[1 5 3 1]]
;;; [[1 4 3 1]]
;;; [[3 1 2 2] [0 1 1 5]]
;;; [[3 1 5 0]]
;;; [[2 0 0 6]]
;;; [[1 5 1 3]]
;;; [[1 0 4 1]]
;;; [[6 3 0 1]]
;;; [[4 1 1 0]]
;;; [[5 1 0 3]]
;;; [[4 1 3 0]]
;;; [[5 0 2 0]]
;;; [[3 1 0 1]]
;;; [[3 2 2 3] [6 2 1 0]]
;;; [[5 1 0 3] [0 0 0 8]]
;;; [[2 4 1 3] [0 0 0 8]]
;;; [[0 8 0 1]]
;;; [[0 6 1 3]]
;;; [[0 0 0 8]]
;;; [[3 0 4 3]]
;;; [[0 0 0 8]]
;;; [[7 6 0 0]]
;;; [[0 5 2 2]]
;;; [[0 0 0 8]]
;;; [[2 5 0 3]]
;;; [[0 1 1 4]]
;;; [[2 5 0 3]]
;;; [[2 5 2 1]]
;;; [[1 0 0 6]]
;;; [[1 0 0 0]]
;;; [[0 5 2 2]]
;;; [[0 1 4 5]]
;;; [[2 2 0 4]]
;;; [[0 0 0 8]]
;;; [[3 0 4 3]]
;;; [[2 5 2 1]]
;;; [[2 1 2 2]]
;;; [[0 0 0 8]]
;;; [[2 3 2 2]]
;;; [[0 5 2 2]]
;;; [[1 0 0 6]]
;;; [[0 5 2 2]]
;;; [[1 2 0 1]]
;;; [[2 2 0 4]]
;;; [[2 5 1 1]]
;;; [[1 5 1 0]]
;;; [[4 2 1 2]]
;;; [[0 0 0 8]]
;;; [[3 3 0 3]]
;;; [[1 6 0 2]]
;;; [[2 1 2 2]]
;;; [[1 0 0 0]]
;;; [[1 5 3 0]]
;;; [[8 0 1 0]]
;;; [[0 0 0 8]]
;;; [[1 0 3 5]]
;;; [[3 2 2 3]]
;;; [[3 0 4 3]]
;;; [[2 2 0 4]]
;;; [[3 0 4 3]]
;;; [[2 5 2 1]]
;;; [[0 7 2 0]]
;;; [[1 6 0 2]]
;;; [[2 0 1 1]]
;;; [[1 2 0 1]]
;;; [[3 0 0 7]]
;;; [[1 2 0 1]]
;;; [[2 3 2 2]]
;;; [[3 3 0 3]]
;;; [[2 5 1 1]]
;;; [[0 3 0 7]]
;;; [[1 6 0 2]]
;;; [[0 0 6 3]]
;;; [[3 0 4 3]]
;;; [[1 0 0 6]]
;;; [[1 5 3 0]]
;;; [[1 0 0 0]]
;;; [[0 4 2 0]]
;;; [[3 2 2 3]]
;;; [[2 5 2 1]]
;;; [[0 0 0 8]]
;;; [[0 0 0 8]]
;;; [[0 4 0 1]]
;;; [[3 2 2 3]]
;;; [[2 0 1 7]]
;;; [[0 5 2 2]]
;;; [[1 0 0 0]]
;;; [[0 5 2 2]]
;;; [[0 0 0 8]]
;;; [[3 0 0 7]]
;;; [[0 0 0 8]]
;;; [[3 2 2 3]]
;;; [[4 2 1 2]]
;;; [[2 3 2 2]]
;;; [[1 6 0 2]]
;;; [[1 5 3 0]]
;;; [[3 3 0 3]]
;;; [[1 2 0 1]]
;;; 
;; <-
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}],"value":"[nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil]"}
;; <=

;; @@

;; @@
